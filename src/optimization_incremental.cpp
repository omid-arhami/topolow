// src/optimization_incremental.cpp
// Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
// Incremental optimization backend for topolow euclidean embedding
//
// =============================================================================
// PURPOSE
// =============================================================================
// This file implements incremental embedding: adding new points to an existing
// map while keeping existing ("fixed") points stationary. The goal is to produce
// results equivalent to re-optimizing the entire map from scratch.
//
// =============================================================================
// KEY DESIGN DECISIONS AND THEIR RATIONALE
// =============================================================================
//
// DECISION 1: FIXED POINT MASK
//    Problem:  We want to add new points without disturbing existing positions.
//    Solution: Points are marked as fixed (1) or optimizable (0) via a mask.
//              Fixed points never move; their positions come from the existing map.
//              Only new points are optimized based on new measurements.
//
// DECISION 2: EDGE FILTERING (handled in R before calling C++)
//    Problem:  Processing all n² edges is wasteful when adding few new points.
//    Solution: Only edges involving at least one new point are passed to C++.
//              Fixed-Fixed edges are excluded since neither endpoint can move.
//              For 1000 fixed + 10 new points: ~20,000 edges vs ~500,000.
//
// DECISION 3: FORCE COMPENSATION FOR FIXED ENDPOINTS
//    Problem:  In standard optimization, spring forces move BOTH endpoints.
//              When one endpoint is fixed, only one point moves.
//              Result: Distance changes less per iteration, causing slow/biased
//              convergence and potentially different final positions.
//    
//    Mathematical Analysis:
//      Let F = force magnitude, norm_i and norm_j = normalization factors.
//      
//      Original (both free):
//        Δpos_i = -F/norm_i,  Δpos_j = +F/norm_j
//        Total Δdist ≈ F/norm_i + F/norm_j = F(norm_i + norm_j)/(norm_i * norm_j)
//      
//      Incremental (j fixed):
//        Δpos_i = -F/norm_i,  Δpos_j = 0
//        Total Δdist = F/norm_i
//      
//      Ratio = norm_j / (norm_i + norm_j)  [NOT 1/2 unless norm_i = norm_j]
//    
//    Solution: Scale force applied to movable point by compensation factor:
//        compensation = (norm_i + norm_j) / norm_j = 1 + norm_i/norm_j
//      
//      Then: Δdist = F * compensation / norm_i 
//                  = F(norm_i + norm_j)/(norm_i * norm_j)  ✓ (matches original)
//
// DECISION 4: EXTRA REPULSION PASS
//    Problem:  Edge filtering removes Fixed-Fixed edges, which means the edge
//              sampling loop visits fewer edges total. Since negative sampling
//              happens per-edge-visit, new points receive fewer repulsion
//              opportunities from fixed points they have no measurements with.
//    Solution: After processing edges, explicitly apply repulsion from randomly
//              sampled fixed points to each new point (n_negative_samples each).
//    Note:     No force compensation needed here - negative sampling in the
//              original algorithm also only updates one point at a time.
//    Note:     No need to skip connected pairs - additional repulsion is harmless
//              and the force is small.
//
// DECISION 5: CONTRADICTORY DATA PREVENTION (handled in R)
//    Problem:  If new_measurements contains a distance between two fixed points
//              that differs from their actual geometric distance, optimization
//              would try to satisfy an impossible constraint.
//    Solution: In R, fixed-fixed distances come only from existing positions.
//              Any fixed-fixed measurements in new_measurements are ignored.
//
// DECISION 6: IMMEDIATE (GAUSS-SEIDEL) UPDATES
//    Problem:  Batched updates (Jacobi-style) can oscillate and converge slower.
//    Solution: Positions are updated immediately when visited, matching the
//              original topolow algorithm's behavior.
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

// ---------------------------------------------------------------------------
// Helper: Fast O(1) lookup in flattened has_measurement vector
// has_measurement is stored row-major: index = i * n + j
// ---------------------------------------------------------------------------
inline bool is_connected_incr(int i, int j, const std::vector<int>& has_meas, int n) {
  return has_meas[i * n + j] != 0;
}

// ---------------------------------------------------------------------------
// Helper: Vectorized MAE calculation over all edges
// Returns: {total_absolute_error, count_of_contributing_edges}
// ---------------------------------------------------------------------------
inline std::pair<double, int> compute_error_vectorized_incr(
    const arma::mat& pos,
    const arma::uvec& edge_i_vec,
    const arma::uvec& edge_j_vec,
    const arma::vec& target_vec,
    const arma::ivec& thresh_vec
) {
  // Extract all source and destination positions at once
  arma::mat pos_i = pos.rows(edge_i_vec);
  arma::mat pos_j = pos.rows(edge_j_vec);
  
  // Compute all deltas at once
  arma::mat deltas = pos_j - pos_i;
  
  // Compute all distances: sqrt(sum of squares along rows)
  arma::vec distances = arma::sqrt(arma::sum(arma::square(deltas), 1));
  
  // Compute absolute errors for all edges
  arma::vec abs_errors = arma::abs(target_vec - distances);
  
  // Determine which edges contribute to error based on threshold type
  arma::uvec exact_mask = (thresh_vec == 0);
  arma::uvec gt_violated = (thresh_vec == 1) % (distances < target_vec);
  arma::uvec lt_violated = (thresh_vec == -1) % (distances > target_vec);
  arma::uvec contributes = exact_mask + gt_violated + lt_violated;
  // Sum errors only for contributing edges
  double total_error = arma::accu(abs_errors % arma::conv_to<arma::vec>::from(contributes));
  int count = arma::accu(contributes);
  
  return {total_error, count};
}

// ---------------------------------------------------------------------------
// Main incremental optimization function
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List optimize_layout_incremental_cpp(
    NumericMatrix initial_positions,
    const IntegerVector& edge_i,              // COO row indices (0-based)
    const IntegerVector& edge_j,              // COO column indices (0-based)
    const NumericVector& edge_dist,           // Target distances
    const IntegerVector& edge_thresh,         // Threshold: 1(>), -1(<), 0(exact)
    const IntegerVector& has_measurement_flat,// Flattened n*n connectivity matrix
    const IntegerVector& degrees,             // Node degrees for normalization
    const IntegerVector& fixed_point_mask,    // 1 = fixed, 0 = optimizable
    const IntegerVector& fixed_indices,       // Indices of fixed points (0-based)
    const IntegerVector& new_indices,         // Indices of new points (0-based)
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
  
  // =========================================================================
  // INPUT VALIDATION
  // =========================================================================
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
  if (fixed_point_mask.size() != n) {
    Rcpp::stop("fixed_point_mask must have n elements");
  }
  
  // =========================================================================
  // COPY DATA TO EFFICIENT STRUCTURES
  // =========================================================================
  arma::mat pos(initial_positions.begin(), n, dim, true);
  
  std::vector<int> fixed_mask(fixed_point_mask.begin(), fixed_point_mask.end());
  std::vector<int> fixed_idx(fixed_indices.begin(), fixed_indices.end());
  std::vector<int> new_idx(new_indices.begin(), new_indices.end());
  std::vector<int> has_meas(has_measurement_flat.begin(), has_measurement_flat.end());
  
  const int n_fixed = fixed_idx.size();
  const int n_new = new_idx.size();
  
  if (verbose) {
    Rcpp::Rcout << "Incremental optimization: " << n_fixed << " fixed, "
                << n_new << " new, " << num_edges << " edges, " 
                << dim << "D\n";
  }
  
  // Edge data in contiguous arrays
  std::vector<int> e_i(edge_i.begin(), edge_i.end());
  std::vector<int> e_j(edge_j.begin(), edge_j.end());
  std::vector<double> e_target(edge_dist.begin(), edge_dist.end());
  std::vector<int> e_thresh(edge_thresh.begin(), edge_thresh.end());
  
  // For vectorized error calculation
  arma::uvec edge_i_vec(num_edges), edge_j_vec(num_edges);
  arma::vec target_vec(num_edges);
  arma::ivec thresh_vec(num_edges);
  for (int e = 0; e < num_edges; ++e) {
    edge_i_vec(e) = edge_i[e];
    edge_j_vec(e) = edge_j[e];
    target_vec(e) = edge_dist[e];
    thresh_vec(e) = edge_thresh[e];
  }
  
  // Degree-based normalization factors
  std::vector<double> deg_plus_one(n);
  for (int i = 0; i < n; ++i) {
    deg_plus_one[i] = static_cast<double>(degrees[i]) + 1.0;
  }
  
  // =========================================================================
  // RANDOM NUMBER GENERATION SETUP
  // =========================================================================
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> node_dist(0, n - 1);
  std::uniform_int_distribution<int> fixed_dist(0, std::max(1, n_fixed) - 1);
  
  // Pre-allocate
  std::vector<int> edge_order(num_edges);
  std::iota(edge_order.begin(), edge_order.end(), 0);
  
  const int total_neg_samples = num_edges * n_negative_samples * 2;
  std::vector<int> neg_samples(total_neg_samples);
  
  // =========================================================================
  // CONVERGENCE TRACKING
  // =========================================================================
  double k = k0;
  double prev_error = std::numeric_limits<double>::max();
  int converge_count = convergence_window;
  bool converged = false;
  int final_iter = n_iter;
  double final_mae = 0.0;
  
  if (convergence_check_freq < 1) convergence_check_freq = 10;
  
  // =========================================================================
  // MAIN OPTIMIZATION LOOP
  // =========================================================================
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // Shuffle edges for stochastic optimization
    std::shuffle(edge_order.begin(), edge_order.end(), rng);
    
    // Pre-generate random samples for negative sampling
    for (int s = 0; s < total_neg_samples; ++s) {
      neg_samples[s] = node_dist(rng);
    }
    
    // -----------------------------------------------------------------------
    // PROCESS EDGES
    // -----------------------------------------------------------------------
    for (int e = 0; e < num_edges; ++e) {
      const int idx = edge_order[e];
      const int i = e_i[idx];
      const int j = e_j[idx];
      const double target_dist = e_target[idx];
      const int thresh_type = e_thresh[idx];
      
      const bool i_fixed = (fixed_mask[i] == 1);
      const bool j_fixed = (fixed_mask[j] == 1);
      
      // Skip if both fixed (shouldn't happen with R filtering, but safe)
      if (i_fixed && j_fixed) continue;
      
      // Compute current distance
      double dist_sq = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = pos(j, d) - pos(i, d);
        dist_sq += diff * diff;
      }
      const double dist = std::sqrt(dist_sq);
      const double dist_stable = dist + 0.01;  // Prevent division by zero
      
      // Determine if spring (attraction) or repulsion applies
      bool apply_spring;
      if (thresh_type == 0) {
        apply_spring = true;  // Exact: always apply spring
      } else if (thresh_type == 1) {
        apply_spring = (dist < target_dist);  // ">": spring if too close
      } else {
        apply_spring = (dist > target_dist);  // "<": spring if too far
      }
      
      // Normalization factors
      const double deg_i = deg_plus_one[i];
      const double deg_j = deg_plus_one[j];
      const double norm_i = 4.0 * deg_i + k;
      const double norm_j = 4.0 * deg_j + k;
      
      // ---------------------------------------------------------------------
      // FORCE COMPENSATION (Decision 3)
      // When one endpoint is fixed, scale force to achieve equivalent Δdist
      // compensation = (norm_movable + norm_fixed) / norm_fixed
      // ---------------------------------------------------------------------
      double compensation_i = 1.0;
      double compensation_j = 1.0;
      
      if (i_fixed && !j_fixed) {
        // Only j moves: compensate j
        compensation_j = (norm_j + norm_i) / norm_i;
      } else if (!i_fixed && j_fixed) {
        // Only i moves: compensate i
        compensation_i = (norm_i + norm_j) / norm_j;
      }
      // Both free: no compensation (shouldn't happen with filtered edges)
      
      // ---------------------------------------------------------------------
      // APPLY FORCES
      // ---------------------------------------------------------------------
      if (apply_spring) {
        // Spring force: pulls points toward target distance
        const double base_force = 2.0 * k * (target_dist - dist) / dist_stable;
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos(j, d) - pos(i, d);
          double force_d = delta_d * base_force;
          
          // Apply with compensation (immediate Gauss-Seidel update)
          if (!i_fixed) {
            pos(i, d) -= (force_d * compensation_i) / norm_i;
          }
          if (!j_fixed) {
            pos(j, d) += (force_d * compensation_j) / norm_j;
          }
        }
      } else {
        // Repulsion force: pushes points apart (threshold not satisfied)
        const double base_force = c_repulsion / (2.0 * dist_stable * dist_stable * dist_stable);
        
        for (int d = 0; d < dim; ++d) {
          double delta_d = pos(j, d) - pos(i, d);
          double force_d = delta_d * base_force;
          
          // Apply with compensation (immediate Gauss-Seidel update)
          if (!i_fixed) {
            pos(i, d) -= (force_d * compensation_i) / deg_i;
          }
          if (!j_fixed) {
            pos(j, d) += (force_d * compensation_j) / deg_j;
          }
        }
      }
      
      // ---------------------------------------------------------------------
      // NEGATIVE SAMPLING
      // No compensation needed: original algorithm also only updates one point
      // ---------------------------------------------------------------------
      const int sample_base_i = e * n_negative_samples * 2;
      const int sample_base_j = sample_base_i + n_negative_samples;
      
      // Negative sampling for node i (if not fixed)
      if (!i_fixed) {
        for (int s = 0; s < n_negative_samples; ++s) {
          const int rand_node = neg_samples[sample_base_i + s];
          if (rand_node == i || rand_node == j || 
              is_connected_incr(i, rand_node, has_meas, n)) continue;
          
          double rep_dist_sq = 0.0;
          for (int d = 0; d < dim; ++d) {
            double diff = pos(rand_node, d) - pos(i, d);
            rep_dist_sq += diff * diff;
          }
          const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;
          const double rep_force = c_repulsion / (2.0 * rep_dist * rep_dist * rep_dist * deg_i);
          
          for (int d = 0; d < dim; ++d) {
            pos(i, d) -= (pos(rand_node, d) - pos(i, d)) * rep_force;
          }
        }
      }
      
      // Negative sampling for node j (if not fixed)
      if (!j_fixed) {
        for (int s = 0; s < n_negative_samples; ++s) {
          const int rand_node = neg_samples[sample_base_j + s];
          if (rand_node == i || rand_node == j || 
              is_connected_incr(j, rand_node, has_meas, n)) continue;
          
          double rep_dist_sq = 0.0;
          for (int d = 0; d < dim; ++d) {
            double diff = pos(rand_node, d) - pos(j, d);
            rep_dist_sq += diff * diff;
          }
          const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;
          const double rep_force = c_repulsion / (2.0 * rep_dist * rep_dist * rep_dist * deg_j);
          
          for (int d = 0; d < dim; ++d) {
            pos(j, d) -= (pos(rand_node, d) - pos(j, d)) * rep_force;
          }
        }
      }
    }
    
    // -----------------------------------------------------------------------
    // EXTRA REPULSION PASS (Decision 4)
    // Compensate for reduced negative sampling due to edge filtering.
    // Apply repulsion from random fixed points to each new point.
    // No connectivity check needed - extra repulsion is harmless.
    // No compensation needed - matches original negative sampling behavior.
    // -----------------------------------------------------------------------
    if (n_new > 0 && n_fixed > 0) {
      for (int ni = 0; ni < n_new; ++ni) {
        const int new_pt = new_idx[ni];
        const double deg_new = deg_plus_one[new_pt];
        
        for (int s = 0; s < n_negative_samples; ++s) {
          const int fixed_pt = fixed_idx[fixed_dist(rng)];
          
          double rep_dist_sq = 0.0;
          for (int d = 0; d < dim; ++d) {
            double diff = pos(fixed_pt, d) - pos(new_pt, d);
            rep_dist_sq += diff * diff;
          }
          const double rep_dist = std::sqrt(rep_dist_sq) + 0.01;
          const double rep_force = c_repulsion / (2.0 * rep_dist * rep_dist * rep_dist * deg_new);
          
          for (int d = 0; d < dim; ++d) {
            pos(new_pt, d) -= (pos(fixed_pt, d) - pos(new_pt, d)) * rep_force;
          }
        }
      }
    }
    
    // -----------------------------------------------------------------------
    // COOLING
    // -----------------------------------------------------------------------
    k *= (1.0 - cooling_rate);
    
    // -----------------------------------------------------------------------
    // CONVERGENCE CHECK
    // -----------------------------------------------------------------------
    if ((iter + 1) % convergence_check_freq == 0 || iter == n_iter - 1) {
      auto [total_error, error_count] = compute_error_vectorized_incr(
        pos, edge_i_vec, edge_j_vec, target_vec, thresh_vec
      );
      const double current_error = (error_count > 0) ? total_error / error_count : 0.0;
      
      if (verbose && ((iter + 1) % 100 == 0 || iter == n_iter - 1)) {
        Rcpp::Rcout << "Iter " << (iter + 1) << "/" << n_iter 
                    << ", MAE=" << current_error << ", k=" << k << "\n";
      }
      
      if (prev_error < std::numeric_limits<double>::max() && prev_error > 1e-10) {
        double relative_change = (prev_error - current_error) / prev_error;
        if (relative_change >= 0 && relative_change < relative_epsilon) {
          if (--converge_count <= 0) {
            converged = true;
            final_iter = iter + 1;
            final_mae = current_error;
            if (verbose) {
              Rcpp::Rcout << "Converged at iteration " << final_iter << "\n";
            }
            break;
          }
        } else {
          converge_count = convergence_window;
        }
      }
      prev_error = current_error;
      final_mae = current_error;
    }
    
    // Numerical stability check
    if ((iter + 1) % 100 == 0 && !pos.is_finite()) {
      Rcpp::stop("Numerical instability at iteration %d", iter + 1);
    }
    
    // Allow R interrupts
    if ((iter + 1) % 50 == 0) Rcpp::checkUserInterrupt();
  }
  
  if (!converged) final_iter = n_iter;
  
  return List::create(
    Named("positions") = pos,
    Named("converged") = converged,
    Named("iterations") = final_iter,
    Named("final_mae") = final_mae,
    Named("final_k") = k
  );
}