// src/optimization.cpp
// Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
// Optimized C++ backend for topolow euclidean embedding
// 
// Key features:
// - COO (Coordinate) format input to avoid sparse matrix zero-dropping issues
// - Negative sampling to approximate unmeasured pair repulsion
// - Immediate (Gauss-Seidel) position updates for better convergence
// - Edge shuffling for stochastic optimization
// - C++ native RNG for performance
// - Flattened has_measurement vector for fast O(1) lookup
// - Vectorized error calculation using Armadillo

#include <RcppArmadillo.h>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Fast O(1) lookup in flattened has_measurement vector
// has_measurement is stored row-major: index = i * n + j
inline bool is_connected(int i, int j, const std::vector<int>& has_meas, int n) {
  return has_meas[i * n + j] != 0;
}

// Vectorized error calculation
// Returns: {total_error, count} for computing MAE
inline std::pair<double, int> compute_error_vectorized(
    const arma::mat& pos,
    const arma::uvec& edge_i_vec,
    const arma::uvec& edge_j_vec,
    const arma::vec& target_vec,
    const arma::ivec& thresh_vec
) {
  // Extract all source and destination positions at once (vectorized indexing)
  arma::mat pos_i = pos.rows(edge_i_vec);  // num_edges x dim
  arma::mat pos_j = pos.rows(edge_j_vec);  // num_edges x dim
  
  // Compute all deltas at once
  arma::mat deltas = pos_j - pos_i;  // num_edges x dim
  
  // Compute all distances at once: sqrt(sum of squares along rows)
  arma::vec distances = arma::sqrt(arma::sum(arma::square(deltas), 1));
  
  // Compute absolute errors for all edges
  arma::vec abs_errors = arma::abs(target_vec - distances);
  
  // Determine which edges contribute to error based on threshold type
  // thresh == 0: always contributes (exact value)
  // thresh == 1 (>): contributes if distance < target (not satisfied)
  // thresh == -1 (<): contributes if distance > target (not satisfied)
  
  arma::uvec exact_mask = (thresh_vec == 0);
  arma::uvec gt_mask = (thresh_vec == 1) % (distances < target_vec);  // > not satisfied

  arma::uvec lt_mask = (thresh_vec == -1) % (distances > target_vec); // < not satisfied
  
  // Combine masks
  arma::uvec contributes = exact_mask + gt_mask + lt_mask;
  
  // Sum errors only for contributing edges
  double total_error = arma::accu(abs_errors % arma::conv_to<arma::vec>::from(contributes));
  int count = arma::accu(contributes);
  
  return {total_error, count};
}

// [[Rcpp::export]]
List optimize_layout_cpp(
    NumericMatrix initial_positions,
    const IntegerVector& edge_i,          // COO row indices (0-based, upper triangle only)
    const IntegerVector& edge_j,          // COO column indices (0-based, upper triangle only)
    const NumericVector& edge_dist,       // Target distances
    const IntegerVector& edge_thresh,     // Threshold codes: 1(>), -1(<), 0(exact)
    const IntegerVector& has_measurement_flat, // Flattened n*n vector for O(1) lookup
    const IntegerVector& degrees,         // Node degrees for normalization
    int n_points,                         // Number of points (for has_measurement indexing)
    int n_iter,
    double k0,
    double cooling_rate,
    double c_repulsion,
    int n_negative_samples,
    double relative_epsilon,
    int convergence_window,
    int convergence_check_freq,           // How often to check convergence
    bool verbose
) {
  
  const int n = initial_positions.nrow();
  const int dim = initial_positions.ncol();
  
  // Validate inputs
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
  if (edge_j.size() != num_edges || edge_dist.size() != num_edges || 
      edge_thresh.size() != num_edges) {
    Rcpp::stop("Edge vectors must have equal lengths");
  }
  if (has_measurement_flat.size() != n * n) {
    Rcpp::stop("has_measurement_flat must have n*n elements");
  }
  
  // Copy initial positions to Armadillo matrix
  arma::mat pos(initial_positions.begin(), n, dim, true);
  
  // =========================================================================
  // PREPARE EDGE DATA IN CONTIGUOUS ARRAYS FOR CACHE EFFICIENCY
  // =========================================================================
  
  // For the main loop: simple int/double arrays
  std::vector<int> e_i(num_edges);
  std::vector<int> e_j(num_edges);
  std::vector<double> e_target(num_edges);
  std::vector<int> e_thresh(num_edges);
  
  for (int e = 0; e < num_edges; ++e) {
    e_i[e] = edge_i[e];
    e_j[e] = edge_j[e];
    e_target[e] = edge_dist[e];
    e_thresh[e] = edge_thresh[e];
  }
  
  // For vectorized error calculation: Armadillo vectors
  arma::uvec edge_i_vec(num_edges);
  arma::uvec edge_j_vec(num_edges);
  arma::vec target_vec(num_edges);
  arma::ivec thresh_vec(num_edges);
  
  for (int e = 0; e < num_edges; ++e) {
    edge_i_vec(e) = edge_i[e];
    edge_j_vec(e) = edge_j[e];
    target_vec(e) = edge_dist[e];
    thresh_vec(e) = edge_thresh[e];
  }
  
  if (verbose) {
    Rcpp::Rcout << "Optimization initialized: " << n << " points, " 
                << num_edges << " edges, " << dim << " dimensions\n";
  }
  
  // =========================================================================
  // COPY has_measurement TO STD::VECTOR FOR FASTER ACCESS
  // =========================================================================
  std::vector<int> has_meas(has_measurement_flat.begin(), has_measurement_flat.end());
  
  // =========================================================================
  // PRE-COMPUTE DEGREE-BASED NORMALIZATION FACTORS
  // =========================================================================
  std::vector<double> deg_plus_one(n);
  for (int i = 0; i < n; ++i) {
    deg_plus_one[i] = static_cast<double>(degrees[i]) + 1.0;
  }
  
  // =========================================================================
  // SETUP C++ NATIVE RNG
  // =========================================================================
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> node_dist(0, n - 1);
  
  // =========================================================================
  // PRE-ALLOCATE MEMORY
  // =========================================================================
  
  // Edge order for shuffling
  std::vector<int> edge_order(num_edges);
  std::iota(edge_order.begin(), edge_order.end(), 0);
  
  // Pre-allocated vectors for negative sampling (regenerated each iteration)
  const int total_neg_samples = num_edges * n_negative_samples * 2;
  std::vector<int> neg_samples(total_neg_samples);
  
  // =========================================================================
  // CONVERGENCE TRACKING
  // =========================================================================
  double k = k0;  // Current spring constant
  // Initialize current_repulsion to decay alongside k
  double current_repulsion = c_repulsion;
  double prev_error = std::numeric_limits<double>::max();
  int converge_count = convergence_window;
  bool converged = false;
  int final_iter = n_iter;
  double final_mae = 0.0;
  
  // Ensure convergence check frequency is reasonable
  if (convergence_check_freq < 1) convergence_check_freq = 10;
  
  // =========================================================================
  // MAIN OPTIMIZATION LOOP
  // =========================================================================
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // ---------------------------------------------------------------------
    // SHUFFLE EDGE ORDER (critical for stochastic optimization)
    // ---------------------------------------------------------------------
    std::shuffle(edge_order.begin(), edge_order.end(), rng);
    
    // ---------------------------------------------------------------------
    // PRE-GENERATE RANDOM SAMPLES FOR THIS ITERATION
    // ---------------------------------------------------------------------
    for (int s = 0; s < total_neg_samples; ++s) {
      neg_samples[s] = node_dist(rng);
    }
    
    // ---------------------------------------------------------------------
    // PROCESS EDGES IN SHUFFLED ORDER
    // ---------------------------------------------------------------------
    for (int e = 0; e < num_edges; ++e) {
      const int idx = edge_order[e];
      const int i = e_i[idx];
      const int j = e_j[idx];
      const double target_dist = e_target[idx];
      const int thresh_type = e_thresh[idx];
      
      // Compute displacement vector and current distance
      // Using direct pointer access for speed
      const double* pos_i = pos.colptr(0) + i;
      const double* pos_j = pos.colptr(0) + j;
      
      double dist_sq = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = pos_j[d * n] - pos_i[d * n];
        dist_sq += diff * diff;
      }
      const double dist = std::sqrt(dist_sq);
      const double dist_stable = dist + 0.01;
      
      // -----------------------------------------------------------------
      // DETERMINE FORCE TYPE BASED ON THRESHOLD
      // -----------------------------------------------------------------
      bool apply_spring;
      if (thresh_type == 0) {
        apply_spring = true;
      } else if (thresh_type == 1) {  // Greater than (>)
        apply_spring = (dist < target_dist);
      } else {  // Less than (<)
        apply_spring = (dist > target_dist);
      }
      
      // -----------------------------------------------------------------
      // APPLY SPRING OR REPULSION FORCE
      // -----------------------------------------------------------------
      const double deg_i = deg_plus_one[i];
      const double deg_j = deg_plus_one[j];
      
      if (apply_spring) {
        // Spring force: 2 * k * (target - actual) / dist_stable
        const double factor = 2.0 * k * (target_dist - dist) / dist_stable;
        const double norm_i = 4.0 * deg_i + k;
        const double norm_j = 4.0 * deg_j + k;
        
        // IMMEDIATE UPDATE - direct memory access
        double* pos_i_mut = pos.colptr(0) + i;
        double* pos_j_mut = pos.colptr(0) + j;
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_j_mut[d * n] - pos_i_mut[d * n];
          double force_d = delta_d * factor;
          pos_i_mut[d * n] -= force_d / norm_i;
          pos_j_mut[d * n] += force_d / norm_j;
        }
        
      } else {
        // Repulsion: repulsion / (2 * dist^2) / dist_stable
        // Use current_repulsion instead of constant c_repulsion
        const double force_mag = current_repulsion / (2.0 * dist_stable * dist_stable * dist_stable);
        
        double* pos_i_mut = pos.colptr(0) + i;
        double* pos_j_mut = pos.colptr(0) + j;
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_j_mut[d * n] - pos_i_mut[d * n];
          double force_d = delta_d * force_mag;
          pos_i_mut[d * n] -= force_d / deg_i;
          pos_j_mut[d * n] += force_d / deg_j;
        }
      }
      
      // -----------------------------------------------------------------
      // NEGATIVE SAMPLING
      // -----------------------------------------------------------------
      const int sample_base_i = e * n_negative_samples * 2;
      const int sample_base_j = sample_base_i + n_negative_samples;
      
      // Negative sampling for node i
      for (int s = 0; s < n_negative_samples; ++s) {
        const int rand_node = neg_samples[sample_base_i + s];
        
        if (rand_node == i || rand_node == j || is_connected(i, rand_node, has_meas, n)) {
          continue;
        }
        
        const double* pos_rand = pos.colptr(0) + rand_node;
        double* pos_i_mut = pos.colptr(0) + i;
        
        double rep_dist_sq = 0.0;
        for (int d = 0; d < dim; ++d) {
          double diff = pos_rand[d * n] - pos_i_mut[d * n];
          rep_dist_sq += diff * diff;
        }
        const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;
        // Use current_repulsion
        const double rep_force_mag = current_repulsion / (2.0 * rep_dist * rep_dist * rep_dist * deg_i);
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_rand[d * n] - pos_i_mut[d * n];
          pos_i_mut[d * n] -= delta_d * rep_force_mag;
        }
      }
      
      // Negative sampling for node j
      for (int s = 0; s < n_negative_samples; ++s) {
        const int rand_node = neg_samples[sample_base_j + s];
        
        if (rand_node == i || rand_node == j || is_connected(j, rand_node, has_meas, n)) {
          continue;
        }
        
        const double* pos_rand = pos.colptr(0) + rand_node;
        double* pos_j_mut = pos.colptr(0) + j;
        
        double rep_dist_sq = 0.0;
        for (int d = 0; d < dim; ++d) {
          double diff = pos_rand[d * n] - pos_j_mut[d * n];
          rep_dist_sq += diff * diff;
        }
        const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;
        const double rep_force_mag = current_repulsion / (2.0 * rep_dist * rep_dist * rep_dist * deg_j);
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_rand[d * n] - pos_j_mut[d * n];
          pos_j_mut[d * n] -= delta_d * rep_force_mag;
        }
      }
    }
    
    // ---------------------------------------------------------------------
    // COOLING
    // ---------------------------------------------------------------------
    k *= (1.0 - cooling_rate);
    // Cool repulsion at the same rate to maintain force balance
    current_repulsion *= (1.0 - cooling_rate);
    // ---------------------------------------------------------------------
    // CONVERGENCE CHECK (VECTORIZED)
    // ---------------------------------------------------------------------
    if ((iter + 1) % convergence_check_freq == 0 || iter == n_iter - 1) {
      
      std::pair<double, int> error_result = compute_error_vectorized(
        pos, edge_i_vec, edge_j_vec, target_vec, thresh_vec
      );
      const double total_error = error_result.first;
      const int error_count = error_result.second;
      
      const double current_error = (error_count > 0) ? total_error / error_count : 0.0;
      
      if (verbose && ((iter + 1) % 100 == 0 || iter == n_iter - 1)) {
        Rcpp::Rcout << "Iteration " << (iter + 1) << "/" << n_iter 
                    << ", MAE = " << current_error 
                    << ", k = " << k << "\n";
      }
      
      // Check for convergence
      if (prev_error < std::numeric_limits<double>::max()) {
        
        // Calculate the absolute difference
        // This handles cases where error goes slightly UP due to stochasticity
        const double diff = std::abs(prev_error - current_error);

        // UNIFIED CONVERGENCE CRITERIA:
        // 1. (diff < 1e-12): Handles cases where error is effectively zero (perfect fit).
        // 2. (diff < relative_epsilon * prev_error): The standard relative check.
        //    Using multiplication (epsilon * prev) instead of division (diff / prev)
        //    prevents division-by-zero errors if prev_error is 0.
        bool convergence_satisfied = (diff < 1e-12) || (diff < relative_epsilon * prev_error);
        
        if (convergence_satisfied) {
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
        } else {
          // Reset counter if stability is broken to have consecutive satisfied checks
          // However, if error increased, we do not reset to full window, as the algorithm is stochastic 
          // and swings are expected, especially near the optimum.
          // converge_count = convergence_window;
        }
      }
      
      prev_error = current_error;
      final_mae = current_error;
    }
    
    // ---------------------------------------------------------------------
    // NUMERICAL STABILITY CHECK
    // ---------------------------------------------------------------------
    if ((iter + 1) % 100 == 0) {
      if (!pos.is_finite()) {
        Rcpp::stop("Numerical instability detected at iteration %d. "
                   "Consider reducing k0 or c_repulsion.", iter + 1);
      }
    }
    
    // Allow R interrupts
    if ((iter + 1) % 50 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  
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