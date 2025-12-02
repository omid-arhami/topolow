// src/optimization.cpp
// Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
//
// Optimized C++ backend for Topolow Euclidean Embedding
//
// =============================================================================
// This file contains two distinct implementations of the embedding algorithm:
//
// 1. optimize_layout_cpp (APPROXIMATE / SCALABLE)
//    - Complexity: O(E + N*k) per iteration
//    - Strategy: Decoupled Two-Phase (Phase 1: Edges, Phase 2: Node-centric sampling)
//    - Best for: Large datasets (N > 500) where O(N²) is prohibitive
//    - Trade-off: Approximates between-cluster repulsion via negative sampling
//
// 2. optimize_layout_exact_cpp (EXACT / HIGH FIDELITY)
//    - Complexity: O(N²) per iteration
//    - Strategy: Full pairwise iteration over all N(N-1)/2 pairs
//    - Best for: Small to medium datasets (N ≤ 500) or when exact topology is critical
//    - Advantage: Exactly matches the original R algorithm behavior
//
// PHYSICS NOTES (applies to both):
// - Spring constant k DECAYS exponentially: k(t+1) = k(t) * (1 - cooling_rate)
// - Repulsion constant c_repulsion is CONSTANT (does NOT decay)
// - This ensures stable separation as springs weaken in the fine-tuning phase
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
// SHARED HELPER FUNCTIONS
// =============================================================================

/**
 * Fast O(1) lookup to check if two nodes have a measurement.
 * has_measurement is stored as flattened row-major vector: index = i * n + j
 */
inline bool is_connected(int i, int j, const std::vector<int>& has_meas, int n) {
  return has_meas[i * n + j] != 0;
}

/**
 * Vectorized MAE calculation for convergence checking.
 * Computes error only on edges where constraints are violated or exact.
 */
inline std::pair<double, int> compute_error_vectorized(
    const arma::mat& pos,
    const arma::uvec& edge_i_vec,
    const arma::uvec& edge_j_vec,
    const arma::vec& target_vec,
    const arma::ivec& thresh_vec
) {
  arma::mat pos_i = pos.rows(edge_i_vec);
  arma::mat pos_j = pos.rows(edge_j_vec);
  
  arma::mat deltas = pos_j - pos_i;
  arma::vec distances = arma::sqrt(arma::sum(arma::square(deltas), 1));
  arma::vec abs_errors = arma::abs(target_vec - distances);
  
  // Contribution logic:
  // exact (0): always contributes
  // > (1): contributes if distance < target (violated)
  // < (-1): contributes if distance > target (violated)
  arma::uvec exact_mask = (thresh_vec == 0);
  arma::uvec gt_violated = (thresh_vec == 1) % (distances < target_vec);
  arma::uvec lt_violated = (thresh_vec == -1) % (distances > target_vec);
  arma::uvec contributes = exact_mask + gt_violated + lt_violated;
  
  double total_error = arma::accu(abs_errors % arma::conv_to<arma::vec>::from(contributes));
  int count = arma::accu(contributes);
  
  return {total_error, count};
}


// =============================================================================
// 2. EXACT ALGORITHM (O(N²) Full Pairwise)
// =============================================================================
//
// Processes ALL N(N-1)/2 pairs in shuffled order each iteration.
// Exactly matches the original R algorithm:
// - Measured pairs: spring force OR repulsion (based on threshold)
// - Unmeasured pairs: repulsion only
// - Gauss-Seidel (immediate) updates
// - Shuffled order for stochastic optimization
//
// Input: Dense matrices for O(1) random access
// - dissimilarity_matrix: numeric matrix with Inf for unmeasured pairs
// - threshold_matrix: integer matrix with codes (0=exact, 1=>, -1=<)

// Simple pair struct for clean code
struct PairIdx {
  int i;
  int j;
};
// dissimilarity_matrix: Dense: Inf for unmeasured
// threshold_matrix: Dense: threshold codes
// Edge vectors for fast MAE calculation only

// [[Rcpp::export]]
List optimize_layout_exact_cpp(
    NumericMatrix initial_positions,
    NumericMatrix dissimilarity_matrix,  
    IntegerMatrix threshold_matrix,      
    const IntegerVector& degrees,
    const IntegerVector& edge_i,
    const IntegerVector& edge_j,
    const NumericVector& edge_dist,
    const IntegerVector& edge_thresh,
    int n_iter,
    double k0,
    double cooling_rate,
    double c_repulsion,
    double relative_epsilon,
    int convergence_window,
    int convergence_check_freq,
    bool verbose
) {
  
  const int n = initial_positions.nrow();
  const int dim = initial_positions.ncol();
  
  if (n < 2) Rcpp::stop("Need at least 2 points for embedding");
  
  // Initialize position matrix
  arma::mat pos(initial_positions.begin(), n, dim, true);
  
  // Degree + 1 normalization
  std::vector<double> deg_plus_one(n);
  for (int i = 0; i < n; ++i) {
    deg_plus_one[i] = static_cast<double>(degrees[i]) + 1.0;
  }
  
  // Build all pairs for shuffling
  const int num_pairs = n * (n - 1) / 2;
  std::vector<PairIdx> all_pairs;
  all_pairs.reserve(num_pairs);
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      all_pairs.push_back({i, j});
    }
  }
  
  // RNG
  std::random_device rd;
  std::mt19937 rng(rd());
  
  // Dense matrix pointers for O(1) access
  // R matrices are column-major: element [i,j] is at index i + j*n
  const double* dist_ptr = dissimilarity_matrix.begin();
  const int* thresh_ptr = threshold_matrix.begin();
  
  // Armadillo vectors for MAE calculation
  arma::uvec edge_i_vec = arma::conv_to<arma::uvec>::from(as<std::vector<int>>(edge_i));
  arma::uvec edge_j_vec = arma::conv_to<arma::uvec>::from(as<std::vector<int>>(edge_j));
  arma::vec target_vec = arma::conv_to<arma::vec>::from(as<std::vector<double>>(edge_dist));
  arma::ivec thresh_vec = arma::conv_to<arma::ivec>::from(as<std::vector<int>>(edge_thresh));
  
  // State variables
  double k = k0;
  // NOTE: c_repulsion is CONSTANT throughout - no decay!
  
  double prev_error = std::numeric_limits<double>::max();
  int converge_count = convergence_window;
  bool converged = false;
  int final_iter = n_iter;
  double final_mae = 0.0;
  
  if (convergence_check_freq < 1) convergence_check_freq = 10;
  
  if (verbose) {
    Rcpp::Rcout << "=== Exact Algorithm (O(N²) Full Pairwise) ===" << "\n";
    Rcpp::Rcout << "Points: " << n << ", Pairs per iteration: " << num_pairs << "\n";
    Rcpp::Rcout << "Parameters: k0=" << k0 << ", cooling=" << cooling_rate 
                << ", c_rep=" << c_repulsion << "\n";
  }
  
  // =========================================================================
  // MAIN LOOP
  // =========================================================================
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // Shuffle ALL pairs - critical for stochastic optimization
    std::shuffle(all_pairs.begin(), all_pairs.end(), rng);
    
    // Process each pair
    for (const auto& pair : all_pairs) {
      const int i = pair.i;
      const int j = pair.j;
      
      double* pos_i = pos.colptr(0) + i;
      double* pos_j = pos.colptr(0) + j;
      
      // Compute current distance
      double dist_sq = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = pos_j[d * n] - pos_i[d * n];
        dist_sq += diff * diff;
      }
      const double dist = std::sqrt(dist_sq);
      const double dist_stable = dist + 0.01;
      
      // Look up target distance from dense matrix
      // Column-major: [i,j] is at i + j*n
      const double target_dist = dist_ptr[i + j * n];
      
      // Check if this pair has a measurement
      // In R preprocessing, unmeasured pairs are set to Inf
      const bool has_measurement = std::isfinite(target_dist);
      
      const double deg_i = deg_plus_one[i];
      const double deg_j = deg_plus_one[j];
      
      if (has_measurement) {
        // -----------------------------------------------------------------
        // MEASURED PAIR: Apply spring or repulsion based on threshold
        // -----------------------------------------------------------------
        const int thresh_type = thresh_ptr[i + j * n];
        
        // Determine if spring should be applied
        // thresh=0 (exact): always spring
        // thresh=1 (>): spring if dist < target (need to push apart)
        // thresh=-1 (<): spring if dist > target (need to pull together)
        bool apply_spring;
        if (thresh_type == 0) {
          apply_spring = true;
        } else if (thresh_type == 1) {
          apply_spring = (dist < target_dist);
        } else {
          apply_spring = (dist > target_dist);
        }
        
        if (apply_spring) {
          // SPRING FORCE toward target
          const double factor = 2.0 * k * (target_dist - dist) / dist_stable;
          const double norm_i = 4.0 * deg_i + k;
          const double norm_j = 4.0 * deg_j + k;
          
          for (int d = 0; d < dim; ++d) {
            double delta_d = pos_j[d * n] - pos_i[d * n];
            double force_d = delta_d * factor;
            pos_i[d * n] -= force_d / norm_i;
            pos_j[d * n] += force_d / norm_j;
          }
        } else {
          // REPULSION (threshold constraint satisfied)
          const double force_mag = c_repulsion / (2.0 * dist_stable * dist_stable * dist_stable);
          
          for (int d = 0; d < dim; ++d) {
            double delta_d = pos_j[d * n] - pos_i[d * n];
            double force_d = delta_d * force_mag;
            pos_i[d * n] -= force_d / deg_i;
            pos_j[d * n] += force_d / deg_j;
          }
        }
        
      } else {
        // -----------------------------------------------------------------
        // UNMEASURED PAIR: Apply repulsion only
        // -----------------------------------------------------------------
        const double force_mag = c_repulsion / (2.0 * dist_stable * dist_stable * dist_stable);
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos_j[d * n] - pos_i[d * n];
          double force_d = delta_d * force_mag;
          pos_i[d * n] -= force_d / deg_i;
          pos_j[d * n] += force_d / deg_j;
        }
      }
    } // End pair loop
    
    // -----------------------------------------------------------------------
    // COOLING: Only spring constant k decays
    // -----------------------------------------------------------------------
    // IMPORTANT: c_repulsion is CONSTANT (no decay)
    // This matches the original R algorithm behavior
    k *= (1.0 - cooling_rate);
    
    // -----------------------------------------------------------------------
    // CONVERGENCE CHECK
    // -----------------------------------------------------------------------
    if ((iter + 1) % convergence_check_freq == 0 || iter == n_iter - 1) {
      auto error_result = compute_error_vectorized(pos, edge_i_vec, edge_j_vec, target_vec, thresh_vec);
      double current_error = (error_result.second > 0) ? error_result.first / error_result.second : 0.0;
      
      if (verbose && ((iter + 1) % 100 == 0 || iter == n_iter - 1)) {
        Rcpp::Rcout << "Iter " << (iter + 1) << "/" << n_iter 
                    << ", MAE=" << current_error << ", k=" << k << "\n";
      }
      
      if (prev_error < std::numeric_limits<double>::max()) {
        double diff = std::abs(prev_error - current_error);
        if ((diff < 1e-12) || (diff < relative_epsilon * prev_error)) {
          if (--converge_count <= 0) {
            converged = true;
            final_iter = iter + 1;
            final_mae = current_error;
            if (verbose) Rcpp::Rcout << "Converged at iteration " << final_iter << "\n";
            break;
          }
        }
      }
      prev_error = current_error;
      final_mae = current_error;
    }
    
    // Numerical stability check
    if ((iter + 1) % 100 == 0 && !pos.is_finite()) {
      Rcpp::stop("Numerical instability at iteration %d. Reduce k0 or c_repulsion.", iter + 1);
    }
    
    // Allow R interrupts
    if ((iter + 1) % 50 == 0) Rcpp::checkUserInterrupt();
    
  } // End iteration loop
  
  if (!converged) final_iter = n_iter;
  
  return List::create(
    Named("positions") = pos,
    Named("converged") = converged,
    Named("iterations") = final_iter,
    Named("final_mae") = final_mae,
    Named("final_k") = k
  );
}