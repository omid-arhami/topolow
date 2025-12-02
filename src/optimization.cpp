// src/optimization.cpp
// Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
//
// Optimized C++ backend for Topolow Euclidean Embedding
//
// =============================================================================
// ARCHITECTURE: DECOUPLED TWO-PHASE STRATEGY
// =============================================================================
//
// This implementation approximates the original O(N²) Topolow algorithm with
// O(E + N*k) complexity per iteration, where E is the number of measured edges
// and k is the number of negative samples per node.
//
// The original algorithm processes ALL N(N-1)/2 pairs in shuffled order.
// This version achieves similar results by splitting each iteration into:
//
//   Phase 1: Measured Interactions (Springs & Conditional Repulsion)
//            - Iterates over all known edges in shuffled order
//            - Applies spring forces toward target distances
//            - Applies repulsion when threshold constraints are satisfied
//              (e.g., ">5" constraint is satisfied when distance > 5)
//
//   Phase 2: Unmeasured Interactions (Background Repulsion via Negative Sampling)
//            - Iterates over all NODES (not edges) in shuffled order
//            - For each node, samples 'k' valid non-neighbors
//            - Applies repulsive forces bilaterally (both nodes move)
//            - Ensures uniform coverage: every node gets equal sampling density
//              regardless of its degree in the measurement graph
//
// KEY DESIGN DECISIONS:
// - Bilateral force application (Newton's 3rd law): Both nodes in a pair
//   receive equal and opposite forces, matching the original algorithm
// - Constant repulsion: c_repulsion does NOT decay (only spring constant k decays)
// - Gauss-Seidel updates: Positions update immediately for faster convergence
// - Batch RNG: Random numbers pre-generated for efficiency
//
// =============================================================================

#include <RcppArmadillo.h>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/**
 * Fast O(1) lookup to check if two nodes have a measurement between them.
 * 
 * The has_measurement matrix is stored as a flattened row-major vector:
 *   index = i * n + j
 * 
 * @param i First node index (0-based)
 * @param j Second node index (0-based)
 * @param has_meas Flattened n×n measurement indicator matrix
 * @param n Total number of nodes
 * @return true if nodes i and j have a measurement, false otherwise
 */
inline bool is_connected(int i, int j, const std::vector<int>& has_meas, int n) {
  return has_meas[i * n + j] != 0;
}

/**
 * Vectorized Mean Absolute Error (MAE) calculation for convergence checking.
 * 
 * Computes the average absolute difference between target and actual distances
 * for all edges where the constraint contributes to the error:
 *   - Exact values (thresh=0): Always contribute
 *   - Greater-than (thresh=1): Contributes only if distance < target (violated)
 *   - Less-than (thresh=-1): Contributes only if distance > target (violated)
 * 
 * Uses Armadillo vectorized operations for efficiency.
 * 
 * @param pos Current position matrix (n × dim)
 * @param edge_i_vec Source node indices for all edges
 * @param edge_j_vec Target node indices for all edges
 * @param target_vec Target distances for all edges
 * @param thresh_vec Threshold type codes for all edges
 * @return Pair of (total_error, count) for computing MAE = total_error / count
 */
inline std::pair<double, int> compute_error_vectorized(
    const arma::mat& pos,
    const arma::uvec& edge_i_vec,
    const arma::uvec& edge_j_vec,
    const arma::vec& target_vec,
    const arma::ivec& thresh_vec
) {
  // Extract positions for all edge endpoints (vectorized indexing)
  arma::mat pos_i = pos.rows(edge_i_vec);  // num_edges × dim
  arma::mat pos_j = pos.rows(edge_j_vec);  // num_edges × dim
  
  // Compute displacement vectors and Euclidean distances
  arma::mat deltas = pos_j - pos_i;
  arma::vec distances = arma::sqrt(arma::sum(arma::square(deltas), 1));
  
  // Compute absolute errors for all edges
  arma::vec abs_errors = arma::abs(target_vec - distances);
  
  // Determine which edges contribute to the error based on threshold logic:
  //   thresh == 0:  Exact value → always contributes
  //   thresh == 1:  Greater-than (>) → contributes if distance < target (constraint violated)
  //   thresh == -1: Less-than (<) → contributes if distance > target (constraint violated)
  arma::uvec exact_mask = (thresh_vec == 0);
  arma::uvec gt_violated_mask = (thresh_vec == 1) % (distances < target_vec);
  arma::uvec lt_violated_mask = (thresh_vec == -1) % (distances > target_vec);
  
  // Combine masks: an edge contributes if it's exact OR its constraint is violated
  arma::uvec contributes = exact_mask + gt_violated_mask + lt_violated_mask;
  
  // Sum errors only for contributing edges
  double total_error = arma::accu(abs_errors % arma::conv_to<arma::vec>::from(contributes));
  int count = arma::accu(contributes);
  
  return {total_error, count};
}

// =============================================================================
// MAIN OPTIMIZATION FUNCTION
// =============================================================================

/**
 * Core C++ optimization routine for Topolow Euclidean embedding.
 * 
 * Implements a two-phase force-directed layout algorithm:
 *   Phase 1: Process measured edges (springs + conditional repulsion)
 *   Phase 2: Apply background repulsion via negative sampling
 * 
 * @param initial_positions Starting coordinates (n × dim matrix)
 * @param edge_i Source node indices for measured edges (0-based)
 * @param edge_j Target node indices for measured edges (0-based)
 * @param edge_dist Target distances for measured edges
 * @param edge_thresh Threshold codes: 1 for ">", -1 for "<", 0 for exact
 * @param has_measurement_flat Flattened n×n indicator matrix (row-major)
 * @param degrees Number of measurements for each node
 * @param n_points Total number of points
 * @param n_iter Maximum iterations
 * @param k0 Initial spring constant
 * @param cooling_rate Exponential decay rate for spring constant (0 < rate < 1)
 * @param c_repulsion Repulsion constant (CONSTANT throughout optimization)
 * @param n_negative_samples Number of non-neighbor samples per node per iteration
 * @param relative_epsilon Convergence threshold for relative error change
 * @param convergence_window Number of stable iterations required for convergence
 * @param convergence_check_freq How often to check convergence (every N iterations)
 * @param verbose Whether to print progress messages
 * 
 * @return List containing: positions, converged, iterations, final_mae, final_k
 */
// [[Rcpp::export]]
List optimize_layout_cpp(
    NumericMatrix initial_positions,
    const IntegerVector& edge_i,
    const IntegerVector& edge_j,
    const NumericVector& edge_dist,
    const IntegerVector& edge_thresh,
    const IntegerVector& has_measurement_flat,
    const IntegerVector& degrees,
    int n_points,
    int n_iter,
    double k0,
    double cooling_rate,
    double c_repulsion,
    int n_negative_samples,
    double relative_epsilon,
    int convergence_window,
    int convergence_check_freq,
    bool verbose
) {
  
  const int n = initial_positions.nrow();
  const int dim = initial_positions.ncol();
  
  // ===========================================================================
  // INPUT VALIDATION
  // ===========================================================================
  if (n < 2) {
    Rcpp::stop("Need at least 2 points for embedding");
  }
  if (n != n_points) {
    Rcpp::stop("Mismatch between initial_positions rows (%d) and n_points (%d)", n, n_points);
  }
  
  const int num_edges = edge_i.size();
  if (num_edges == 0) {
    Rcpp::stop("No edges provided");
  }
  
  // ===========================================================================
  // DATA STRUCTURE INITIALIZATION
  // ===========================================================================
  
  // Copy initial positions to Armadillo matrix for efficient computation
  // Armadillo uses column-major storage; pos(i, d) accesses row i, column d
  arma::mat pos(initial_positions.begin(), n, dim, true);
  
  // Convert R vectors to std::vector for faster iteration in tight loops
  // (Avoids R's SEXP overhead on each access)
  std::vector<int> e_i = as<std::vector<int>>(edge_i);
  std::vector<int> e_j = as<std::vector<int>>(edge_j);
  std::vector<double> e_target = as<std::vector<double>>(edge_dist);
  std::vector<int> e_thresh = as<std::vector<int>>(edge_thresh);
  std::vector<int> has_meas = as<std::vector<int>>(has_measurement_flat);
  std::vector<int> node_degrees = as<std::vector<int>>(degrees);
  
  // Pre-compute degree-based normalization factors: (degree + 1)
  // Used to dampen movement of highly-connected nodes, preventing oscillation
  // Matches original R code: node_degrees_1 <- node_degrees + 1
  std::vector<double> deg_plus_one(n);
  for (int i = 0; i < n; ++i) {
    deg_plus_one[i] = static_cast<double>(node_degrees[i]) + 1.0;
  }
  
  // Pre-compute number of non-neighbors for each node (used in Phase 2)
  // num_non_neighbors[i] = n - 1 - degree[i] = maximum valid negative samples
  std::vector<int> num_non_neighbors(n);
  for (int i = 0; i < n; ++i) {
    num_non_neighbors[i] = n - 1 - node_degrees[i];
  }
  
  // Prepare Armadillo vectors for vectorized error calculation
  arma::uvec edge_i_vec = arma::conv_to<arma::uvec>::from(e_i);
  arma::uvec edge_j_vec = arma::conv_to<arma::uvec>::from(e_j);
  arma::vec target_vec = arma::conv_to<arma::vec>::from(e_target);
  arma::ivec thresh_vec = arma::conv_to<arma::ivec>::from(e_thresh);
  
  // ===========================================================================
  // RANDOM NUMBER GENERATOR SETUP
  // ===========================================================================
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> node_dist(0, n - 1);
  
  // Index vectors for shuffling (reused each iteration)
  std::vector<int> edge_order(num_edges);
  std::iota(edge_order.begin(), edge_order.end(), 0);
  
  std::vector<int> node_order(n);
  std::iota(node_order.begin(), node_order.end(), 0);
  
  // Pre-allocated buffer for batch random number generation in Phase 2
  // Size: maximum possible attempts for any node
  // We regenerate this buffer for each node based on its specific needs
  const int max_possible_attempts = n_negative_samples * 5 + 20;
  std::vector<int> rand_buffer(max_possible_attempts);
  
  // ===========================================================================
  // OPTIMIZATION STATE VARIABLES
  // ===========================================================================
  
  // Spring constant k decays exponentially: k(t+1) = k(t) * (1 - cooling_rate)
  // High k initially → strong springs pull points toward target distances
  // Low k later → fine-tuning phase, repulsion shapes final layout
  double k = k0;
  
  // NOTE: c_repulsion is CONSTANT throughout optimization (no decay)
  // This matches the original algorithm where only k decays
  // Rationale: As springs weaken, constant repulsion maintains separation
  // and allows the layout to settle into a stable equilibrium
  
  // Convergence tracking
  double prev_error = std::numeric_limits<double>::max();
  int converge_count = convergence_window;  // Counts down to 0 for convergence
  bool converged = false;
  int final_iter = n_iter;
  double final_mae = 0.0;
  
  // Ensure reasonable convergence check frequency
  if (convergence_check_freq < 1) {
    convergence_check_freq = 10;
  }
  
  if (verbose) {
    Rcpp::Rcout << "=== Topolow C++ Optimization ===" << "\n";
    Rcpp::Rcout << "Points: " << n << ", Edges: " << num_edges 
                << ", Dimensions: " << dim << "\n";
    Rcpp::Rcout << "Parameters: k0=" << k0 << ", cooling=" << cooling_rate 
                << ", c_rep=" << c_repulsion << ", neg_samples=" << n_negative_samples << "\n";
    Rcpp::Rcout << "Strategy: Decoupled Two-Phase (Edge Springs + Node-Centric Repulsion)\n";
  }
  
  // ===========================================================================
  // MAIN OPTIMIZATION LOOP
  // ===========================================================================
  
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // -------------------------------------------------------------------------
    // SHUFFLE PROCESSING ORDER (Critical for stochastic optimization)
    // -------------------------------------------------------------------------
    // Shuffling prevents systematic biases and helps escape local minima
    // Different random orders each iteration → better global exploration
    
    std::shuffle(edge_order.begin(), edge_order.end(), rng);
    std::shuffle(node_order.begin(), node_order.end(), rng);
    
    // =========================================================================
    // PHASE 1: MEASURED EDGES (Springs & Conditional Repulsion)
    // =========================================================================
    // Process all edges with known measurements in shuffled order.
    // For each edge (i, j) with target distance d and threshold type t:
    //   - If constraint requires movement toward target → apply spring force
    //   - If constraint is satisfied → apply repulsion to maintain separation
    //
    // Threshold logic:
    //   t=0 (exact):  Always apply spring toward target
    //   t=1 (>d):     If current < target → spring (pull apart)
    //                 If current >= target → repulsion (constraint satisfied)
    //   t=-1 (<d):    If current > target → spring (push together)  
    //                 If current <= target → repulsion (constraint satisfied)
    
    for (int e = 0; e < num_edges; ++e) {
      const int idx = edge_order[e];
      const int i = e_i[idx];
      const int j = e_j[idx];
      const double target_dist = e_target[idx];
      const int thresh_type = e_thresh[idx];
      
      // Get mutable pointers to positions for Gauss-Seidel (immediate) updates
      // Armadillo column-major layout: element (row, col) is at colptr(col)[row]
      // For our n×dim matrix: pos(i, d) is at pos.colptr(d)[i]
      // But we access via pos.colptr(0) + i, then stride by n for each dimension
      double* pos_i = pos.colptr(0) + i;
      double* pos_j = pos.colptr(0) + j;
      
      // Compute current Euclidean distance between nodes i and j
      double dist_sq = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = pos_j[d * n] - pos_i[d * n];
        dist_sq += diff * diff;
      }
      const double dist = std::sqrt(dist_sq);
      
      // Stabilized distance to prevent division by zero
      // Adding small constant ensures numerical stability for overlapping points
      const double dist_stable = dist + 0.01;
      
      // Determine whether to apply spring or repulsion based on threshold logic
      bool apply_spring;
      if (thresh_type == 0) {
        // Exact measurement: always apply spring toward target
        apply_spring = true;
      } else if (thresh_type == 1) {
        // Greater-than constraint (>target): 
        // Apply spring if too close (need to push apart to meet constraint)
        apply_spring = (dist < target_dist);
      } else {
        // Less-than constraint (<target):
        // Apply spring if too far (need to pull together to meet constraint)
        apply_spring = (dist > target_dist);
      }
      
      // Get degree-based normalization factors
      const double deg_i = deg_plus_one[i];
      const double deg_j = deg_plus_one[j];
      
      if (apply_spring) {
        // ---------------------------------------------------------------------
        // SPRING FORCE: Pull/push toward target distance
        // ---------------------------------------------------------------------
        // Force magnitude: F = 2 * k * (target - actual) / actual
        // Positive when target > actual (pull apart), negative when target < actual (push together)
        //
        // Update formula (from original R code):
        //   adjustment_factor = 2 * k * (ideal_distance - distance) * delta / distance
        //   positions[i] -= adjustment_factor / (4 * deg_i + k)
        //   positions[j] += adjustment_factor / (4 * deg_j + k)
        //
        // The factor (4 * deg + k) provides adaptive step size:
        //   - High-degree nodes move less (more constrained by many measurements)
        //   - High k (early iterations) → smaller steps for stability
        
        const double factor = 2.0 * k * (target_dist - dist) / dist_stable;
        const double norm_i = 4.0 * deg_i + k;
        const double norm_j = 4.0 * deg_j + k;
        
        // Apply force bilaterally with Gauss-Seidel updates
        // Node i moves opposite to displacement (toward/away from j)
        // Node j moves along displacement (toward/away from i)
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_j[d * n] - pos_i[d * n];  // Vector from i to j
          double force_d = delta_d * factor;
          
          pos_i[d * n] -= force_d / norm_i;  // i moves opposite
          pos_j[d * n] += force_d / norm_j;  // j moves along
        }
        
      } else {
        // ---------------------------------------------------------------------
        // REPULSION FORCE: Constraint satisfied, maintain separation
        // ---------------------------------------------------------------------
        // When a threshold constraint is satisfied (e.g., distance > target for ">"),
        // we still apply repulsion to prevent collapse and maintain the satisfaction.
        //
        // Force magnitude: F = c_repulsion / (2 * dist^3)
        // This is a 1/r² repulsion (inverse square law), divided by r for direction
        //
        // Update formula (from original R code):
        //   force = c_repulsion / (2 * distance^2) * (delta / distance)
        //   positions[i] -= force / deg_i
        //   positions[j] += force / deg_j
        
        const double force_mag = c_repulsion / (2.0 * dist_stable * dist_stable * dist_stable);
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_j[d * n] - pos_i[d * n];  // Vector from i to j
          double force_d = delta_d * force_mag;
          
          pos_i[d * n] -= force_d / deg_i;  // i pushed away from j
          pos_j[d * n] += force_d / deg_j;  // j pushed away from i
        }
      }
    } // End Phase 1
    
    // =========================================================================
    // PHASE 2: UNMEASURED PAIRS (Background Repulsion via Negative Sampling)
    // =========================================================================
    // Apply repulsion between nodes that have no measurement.
    // This creates a "background pressure" that prevents unmeasured nodes
    // from collapsing toward each other.
    //
    // KEY DESIGN: Iterate over NODES, not edges
    // - Each node samples exactly n_negative_samples non-neighbors
    // - This ensures uniform coverage: low-degree nodes get same sampling
    //   density as high-degree nodes
    // - Prevents the variance/coverage asymmetry of edge-centric sampling
    //
    // BILATERAL FORCES (Newton's 3rd Law):
    // When node i samples node j, BOTH nodes receive force:
    //   - Node i is pushed away from j
    //   - Node j is pushed away from i
    // This matches the original algorithm's pair-wise force application
    
    for (int i : node_order) {
      double* pos_i = pos.colptr(0) + i;
      const double deg_i = deg_plus_one[i];
      
      // Determine how many valid non-neighbors exist for this node
      const int available_non_neighbors = num_non_neighbors[i];
      
      // Skip if node is connected to everyone (no non-neighbors to sample)
      if (available_non_neighbors <= 0) {
        continue;
      }
      
      // Target number of samples (can't exceed available non-neighbors)
      const int target_samples = std::min(n_negative_samples, available_non_neighbors);
      
      // Smart max_attempts: scale based on probability of finding valid sample
      // P(valid) ≈ available_non_neighbors / n
      // Expected attempts to find one valid sample ≈ n / available_non_neighbors
      // Add buffer for variance
      const int expected_attempts_per_sample = (n + available_non_neighbors - 1) / available_non_neighbors;
      const int max_attempts = target_samples * expected_attempts_per_sample * 2 + 10;
      
      // Batch random number generation for efficiency
      // Pre-generate all random node indices we might need
      const int buffer_size = std::min(max_attempts, static_cast<int>(rand_buffer.size()));
      for (int r = 0; r < buffer_size; ++r) {
        rand_buffer[r] = node_dist(rng);
      }
      
      int valid_samples = 0;
      int attempts = 0;
      
      while (valid_samples < target_samples && attempts < max_attempts) {
        // Get next random node from pre-generated buffer
        // If we exhaust buffer (rare), generate on-the-fly
        int rand_node;
        if (attempts < buffer_size) {
          rand_node = rand_buffer[attempts];
        } else {
          rand_node = node_dist(rng);
        }
        attempts++;
        
        // REJECTION CRITERIA:
        // 1. Same node as i (can't repel self)
        // 2. Has measurement with i (handled in Phase 1)
        if (rand_node == i || is_connected(i, rand_node, has_meas, n)) {
          continue;
        }
        
        // Valid non-neighbor found - apply bilateral repulsion
        double* pos_rand = pos.colptr(0) + rand_node;
        const double deg_rand = deg_plus_one[rand_node];
        
        // Compute distance between i and rand_node
        double rep_dist_sq = 0.0;
        for (int d = 0; d < dim; ++d) {
          double diff = pos_rand[d * n] - pos_i[d * n];
          rep_dist_sq += diff * diff;
        }
        const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;  // Stabilized
        
        // Repulsion force magnitude: F = c_repulsion / (2 * dist^3)
        // Same formula as Phase 1 repulsion (consistent force model)
        const double force_mag = c_repulsion / (2.0 * rep_dist * rep_dist * rep_dist);
        
        // Apply BILATERAL force (both nodes move)
        // This ensures Newton's 3rd law and matches original algorithm
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_rand[d * n] - pos_i[d * n];  // Vector from i to rand
          double force_d = delta_d * force_mag;
          
          // Node i pushed away from rand_node (opposite direction)
          pos_i[d * n] -= force_d / deg_i;
          
          // Node rand_node pushed away from i (along direction)
          pos_rand[d * n] += force_d / deg_rand;
        }
        
        valid_samples++;
      }
    } // End Phase 2
    
    // =========================================================================
    // COOLING: Decay spring constant
    // =========================================================================
    // Exponential decay: k(t+1) = k(t) * (1 - cooling_rate)
    //
    // IMPORTANT: Only k decays, NOT c_repulsion
    // Rationale from original algorithm:
    //   - Early iterations: Strong springs dominate → rapid convergence toward targets
    //   - Late iterations: Weak springs, constant repulsion → fine-tuning phase
    //   - Constant repulsion ensures stable separation even as springs weaken
    
    k *= (1.0 - cooling_rate);
    
    // =========================================================================
    // CONVERGENCE CHECK
    // =========================================================================
    // Check for convergence every convergence_check_freq iterations
    // Convergence is declared when relative error change is small for
    // 'convergence_window' consecutive or cumulative checks
    
    if ((iter + 1) % convergence_check_freq == 0 || iter == n_iter - 1) {
      
      // Compute current MAE using vectorized calculation
      std::pair<double, int> error_result = compute_error_vectorized(
        pos, edge_i_vec, edge_j_vec, target_vec, thresh_vec
      );
      
      const double current_error = (error_result.second > 0) ? 
                                    error_result.first / error_result.second : 0.0;
      
      // Verbose progress output
      if (verbose && ((iter + 1) % 100 == 0 || iter == n_iter - 1)) {
        Rcpp::Rcout << "Iteration " << (iter + 1) << "/" << n_iter 
                    << ", MAE = " << current_error 
                    << ", k = " << k << "\n";
      }
      
      // Check for convergence based on relative error change
      if (prev_error < std::numeric_limits<double>::max()) {
        const double diff = std::abs(prev_error - current_error);
        
        // Convergence criteria (matches original R logic):
        // 1. Absolute change is negligible (< 1e-12), OR
        // 2. Relative change is below threshold (< relative_epsilon * prev_error)
        bool is_stable = (diff < 1e-12) || (diff < relative_epsilon * prev_error);
        
        if (is_stable) {
          // Error is stable - decrement convergence counter
          // Note: Counter accumulates across iterations (not strictly consecutive)
          // This is intentional: allows for small fluctuations in stochastic optimization
          --converge_count;
          
          if (converge_count <= 0) {
            converged = true;
            final_iter = iter + 1;
            final_mae = current_error;
            if (verbose) {
              Rcpp::Rcout << "Convergence achieved at iteration " << final_iter 
                          << " (MAE = " << final_mae << ")\n";
            }
            break;
          }
        }
        // Note: We do NOT reset converge_count on unstable iterations
        // This allows convergence to be detected even with occasional fluctuations
        // which are expected in stochastic optimization
      }
      
      prev_error = current_error;
      final_mae = current_error;
    }
    
    // =========================================================================
    // NUMERICAL STABILITY CHECK
    // =========================================================================
    // Periodically verify that positions remain finite
    // Non-finite values indicate explosion (k0 or c_repulsion too high)
    
    if ((iter + 1) % 100 == 0) {
      if (!pos.is_finite()) {
        Rcpp::stop("Numerical instability detected at iteration %d. "
                   "Consider reducing k0 or c_repulsion.", iter + 1);
      }
    }
    
    // Allow R to process interrupts (Ctrl+C)
    if ((iter + 1) % 50 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
  } // End main iteration loop
  
  // ===========================================================================
  // RETURN RESULTS
  // ===========================================================================
  
  if (!converged) {
    final_iter = n_iter;
  }
  
  return List::create(
    Named("positions") = pos,
    Named("converged") = converged,
    Named("iterations") = final_iter,
    Named("final_mae") = final_mae,
    Named("final_k") = k
  );
}
