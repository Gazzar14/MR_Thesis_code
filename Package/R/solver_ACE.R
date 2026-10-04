#' Solve a Single ACE Bound via Column Generation
#' @keywords internal
solve_ace_bound <- function(target_a, target_b, is_upper_bound, slack_val, 
                            tolerance, max_iter, types_Y, prep_data, verbose) {
  
  if (verbose) {
    cat("\n==================================================\n")
    cat(
      ifelse(is_upper_bound,
             paste0("UPPER BOUND for ACE(", target_a, ",", target_b, ")"),
             paste0("LOWER BOUND for ACE(", target_a, ",", target_b, ")")),
      "\n"
    )
    cat("==================================================\n")
  }
  
  b <- prep_data$b             
  num_obs <- prep_data$num_obs 
  obs_idx_array <- prep_data$obs_idx_array
  
  # Adjusted notation: m is treatment, k is instrument/group
  m <- prep_data$m 
  k <- prep_data$k 
  
  c_Y_base <- types_Y[, paste0("X", target_b)] - types_Y[, paste0("X", target_a)]
  c_Y <- if (is_upper_bound) -c_Y_base else c_Y_base
  
  M <- 1e4
  A_master <- Matrix::sparseMatrix(
    i = 1:(num_obs + 1),
    j = 1:(num_obs + 1),
    x = 1,
    dims = c(num_obs + 1, num_obs + 1)
  )
  A_master <- methods::as(A_master, "dgCMatrix")
  c_master <- rep(M, num_obs + 1)
  
  optimal_found <- FALSE
  start_time <- Sys.time()
  res <- NULL
  slack_vector <- c(rep(slack_val, num_obs), 0)
  
  for (iter in 1:max_iter) {
    res <- highs::highs_solve(
      L = c_master,
      lower = rep(0, ncol(A_master)),
      upper = rep(1, ncol(A_master)),
      A = A_master,
      lhs = b - slack_vector, 
      rhs = b + slack_vector, 
      maximum = FALSE
    )
    
    duals <- NULL
    if (!is.null(res$solver_msg) && is.list(res$solver_msg) && "row_dual" %in% names(res$solver_msg)) {
      duals <- res$solver_msg$row_dual
    } else if ("row_dual" %in% names(res)) {
      duals <- res$row_dual
    } else if ("dual_solution" %in% names(res)) {
      duals <- res$dual_solution
    } else if (!is.null(res$info) && "row_dual" %in% names(res$info)) {
      duals = res$info$row_dual
    }
    
    if (is.null(duals)) stop("highs_solve failed to return row dual variables.")
    
    pi_duals <- duals[seq_len(num_obs)]
    pi_norm  <- duals[num_obs + 1]
    
    best_rc <- Inf
    best_col_indices <- NULL
    best_c <- NULL
    
    for (i in seq_len(nrow(types_Y))) {
      current_Y <- types_Y[i, ]
      current_c <- c_Y[i]
      pi_sum <- 0
      
      # Adjusted logic: Construct optimal X^G by iterating over the k levels of G
      current_indices <- integer(k) 
      
      for (g_val in 0:(k - 1)) {
        best_pi <- -Inf
        best_idx <- NA_integer_
        
        # Iterate over the m levels of treatment X
        for (x_val in 0:(m - 1)) {
          y_val <- current_Y[x_val + 1]
          idx <- obs_idx_array[g_val + 1, x_val + 1, y_val + 1]
          
          p_val <- if (!is.na(idx)) pi_duals[idx] else 0
          
          if (p_val > best_pi) {
            best_pi <- p_val
            best_idx <- idx
          }
        }
        pi_sum <- pi_sum + best_pi
        current_indices[g_val + 1] <- best_idx
      }
      
      rc <- current_c - (pi_sum + pi_norm)
      if (rc < best_rc) {
        best_rc <- rc
        best_col_indices <- current_indices
        best_c <- current_c
      }
    }
    
    if (best_rc > -tolerance) {
      if (verbose) cat("Optimal bound found at iteration", iter, "\n")
      optimal_found <- TRUE
      break
    }
    
    i_new <- c(best_col_indices[!is.na(best_col_indices)], num_obs + 1)
    new_col_matrix <- Matrix::sparseMatrix(
      i = i_new, j = rep(1, length(i_new)), x = rep(1, length(i_new)),
      dims = c(num_obs + 1, 1)
    )
    
    A_master <- Matrix::cbind2(A_master, methods::as(new_col_matrix, "dgCMatrix"))
    c_master <- c(c_master, best_c)
    
    if (verbose && (iter %% 50 == 0 || iter == 1)) {
      cat(
        "Iter:", sprintf("%4d", iter),
        "| Obj:", sprintf("%10.5f", res$objective_value),
        "| Best RC:", sprintf("%10.5f", best_rc),
        "| Cols:", ncol(A_master), "\n"
      )
    }
  }
  
  if (!optimal_found && verbose) {
    warning("Reached max iterations without proving optimality.")
  }
  
  run_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  q_weights <- res$primal_solution
  art_vars <- q_weights[1:(num_obs + 1)]
  gen_vars <- q_weights[(num_obs + 2):length(q_weights)]
  
  sum_art <- sum(art_vars)
  if (sum_art > 1e-6) {
    if (verbose) {
      warning(sprintf(
        "Empirical distribution is outside the IV polytope by more than slack_val. Artificial variables used: %f", 
        sum_art
      ))
    }
    final_obj <- NA 
  } else {
    true_obj <- sum(gen_vars * c_master[(num_obs + 2):length(c_master)])
    final_obj <- if (is_upper_bound) -true_obj else true_obj
  }
  
  list(
    objective = final_obj,
    iterations = iter,
    time_secs = run_time,
    final_cols = ncol(A_master),
    infeasible_gap = sum_art
  )
}