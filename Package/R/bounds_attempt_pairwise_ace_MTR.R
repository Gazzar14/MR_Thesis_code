# library(dplyr)
# library(tidyr)
# library(Matrix)
# library(highs)
# 
# k <- 10
# m <- 5
# 
# data_dir <- "~/projects/Highs_bounds/data/"
# load(file.path(data_dir, "data.rda"))
# 
# # ------------------------------------------------------------
# # 0. Output folder + run metadata
# # ------------------------------------------------------------
# run_dir <- file.path(
#   data_dir,
#   paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
# )
# dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
# 
# saveRDS(
#   list(k = k, m = m, start_time = Sys.time(), pid = Sys.getpid()),
#   file.path(run_dir, "run_meta.rds")
# )
# 
# # ------------------------------------------------------------
# # 1. Data prep
# # ------------------------------------------------------------
# prob_df <- dat_10 %>%
#   mutate(X = X - 1) %>%
#   count(G, X, Y) %>%
#   complete(G = 0:(m - 1), X = 0:(k - 1), Y = 0:1, fill = list(n = 0)) %>%
#   group_by(G) %>%
#   mutate(prob = if (sum(n) == 0) NA_real_ else n / sum(n)) %>%
#   ungroup() %>%
#   mutate(param_name = paste0("p", Y, X, "_", G))
# 
# obs_grid <- expand.grid(
#   Y = c(0, 1),
#   X = 0:(k - 1),
#   G = 0:(m - 1),
#   KEEP.OUT.ATTRS = FALSE,
#   stringsAsFactors = FALSE
# )
# 
# prob_map <- obs_grid %>%
#   left_join(prob_df, by = c("Y", "X", "G")) %>%
#   mutate(prob = ifelse(is.na(prob), 0, prob))
# 
# num_obs <- nrow(obs_grid)
# b <- c(prob_map$prob, 1)  # RHS vector (observed probs + normalization)
# 
# # ------------------------------------------------------------
# # 2. Latent Y-space (with MTR)
# # ------------------------------------------------------------
# # Generate all possible response types
# types_Y_all <- as.matrix(expand.grid(rep(list(c(0, 1)), k)))
# colnames(types_Y_all) <- paste0("X", 0:(k - 1))
# 
# # Apply Monotone Treatment Response (MTR) assumption:
# # Keep only rows where the outcome is weakly increasing with treatment.
# # diff(row) >= 0 ensures that it never goes from 1 down to 0 as X increases.
# is_mtr <- apply(types_Y_all, 1, function(row) all(diff(row) >= 0))
# 
# # Filter the latent space
# types_Y <- types_Y_all[is_mtr, ]
# # ------------------------------------------------------------
# # 3. O(1) lookup array
# # ------------------------------------------------------------
# obs_idx_array <- array(NA_integer_, dim = c(m, k, 2))
# for (r in seq_len(num_obs)) {
#   g_idx <- obs_grid$G[r] + 1
#   x_idx <- obs_grid$X[r] + 1
#   y_idx <- obs_grid$Y[r] + 1
#   obs_idx_array[g_idx, x_idx, y_idx] <- r
# }
# 
# # ------------------------------------------------------------
# # 4. Column generation solver for ACE(a,b)
# # ------------------------------------------------------------
# 
# run_column_generation <- function(target_a, target_b, is_upper_bound = FALSE, slack_val = 0.001) {
# 
#   stopifnot(target_a >= 0, target_a < k)
#   stopifnot(target_b >= 0, target_b < k)
#   stopifnot(target_a != target_b)
# 
#   cat("\n==================================================\n")
#   cat(
#     ifelse(is_upper_bound,
#            paste0("UPPER BOUND for ACE(", target_a, ",", target_b, ")"),
#            paste0("LOWER BOUND for ACE(", target_a, ",", target_b, ")")),
#     "\n"
#   )
#   cat("==================================================\n")
# 
#   # CORRECTED: Objective for latent Y-types is now Y(target_b) - Y(target_a)
#   c_Y_base <- types_Y[, paste0("X", target_b)] - types_Y[, paste0("X", target_a)]
# 
#   # For upper bound, minimize the negative objective
#   c_Y <- if (is_upper_bound) -c_Y_base else c_Y_base
# 
#   # Initialize Restricted Master Problem
#   M <- 1e4 # Lowered from 1e6 for better numerical stability
#   A_master <- sparseMatrix(
#     i = 1:(num_obs + 1),
#     j = 1:(num_obs + 1),
#     x = 1,
#     dims = c(num_obs + 1, num_obs + 1)
#   )
#   A_master <- as(A_master, "dgCMatrix")
#   c_master <- rep(M, num_obs + 1)
# 
#   tolerance <- 1e-6
#   max_iter <- 2000
#   optimal_found <- FALSE
#   start_time <- Sys.time()
#   res <- NULL
#   iter <- 0
# 
#   # Create a slack vector to absorb finite sampling noise.
#   slack_vector <- c(rep(slack_val, num_obs), 0)
# 
#   for (iter in 1:max_iter) {
# 
#     # LHS and RHS bounds incorporate the slack
#     res <- highs_solve(
#       L = c_master,
#       lower = rep(0, ncol(A_master)),
#       upper = rep(1, ncol(A_master)),
#       A = A_master,
#       lhs = b - slack_vector,
#       rhs = b + slack_vector,
#       maximum = FALSE
#     )
# 
#     # Extract duals robustly
#     duals <- NULL
#     if (!is.null(res$solver_msg) && is.list(res$solver_msg) && "row_dual" %in% names(res$solver_msg)) {
#       duals <- res$solver_msg$row_dual
#     } else if ("row_dual" %in% names(res)) {
#       duals <- res$row_dual
#     } else if ("dual_solution" %in% names(res)) {
#       duals <- res$dual_solution
#     } else if (!is.null(res$info) && "row_dual" %in% names(res$info)) {
#       duals <- res$info$row_dual
#     }
# 
#     if (is.null(duals)) {
#       stop("highs_solve did not return dual variables (row_dual).")
#     }
# 
#     pi_duals <- duals[seq_len(num_obs)]
#     pi_norm  <- duals[num_obs + 1]
# 
#     # Pricing subproblem
#     best_rc <- Inf
#     best_col_indices <- NULL
#     best_c <- NULL
# 
#     for (i in seq_len(nrow(types_Y))) {
#       current_Y <- types_Y[i, ]
#       current_c <- c_Y[i]
# 
#       pi_sum <- 0
#       current_indices <- integer(m)
# 
#       for (g_val in 0:(m - 1)) {
#         best_pi <- -Inf
#         best_idx <- NA_integer_
# 
#         for (x_val in 0:(k - 1)) {
#           y_val <- current_Y[x_val + 1]
#           idx <- obs_idx_array[g_val + 1, x_val + 1, y_val + 1]
# 
#           if (!is.na(idx)) {
#             p_val <- pi_duals[idx]
#             if (p_val > best_pi) {
#               best_pi <- p_val
#               best_idx <- idx
#             }
#           }
#         }
# 
#         pi_sum <- pi_sum + best_pi
#         current_indices[g_val + 1] <- best_idx
#       }
# 
#       rc <- current_c - (pi_sum + pi_norm)
# 
#       if (rc < best_rc) {
#         best_rc <- rc
#         best_col_indices <- current_indices
#         best_c <- current_c
#       }
#     }
# 
#     # Termination
#     if (best_rc > -tolerance) {
#       cat("Optimal bound found at iteration", iter, "\n")
#       optimal_found <- TRUE
#       break
#     }
# 
#     # Add entering column
#     i_new <- c(best_col_indices[!is.na(best_col_indices)], num_obs + 1)
#     j_new <- rep(1, length(i_new))
#     x_new <- rep(1, length(i_new))
# 
#     new_col_matrix <- sparseMatrix(
#       i = i_new,
#       j = j_new,
#       x = x_new,
#       dims = c(num_obs + 1, 1)
#     )
# 
#     A_master <- cbind(A_master, as(new_col_matrix, "dgCMatrix"))
#     c_master <- c(c_master, best_c)
# 
#     if (iter %% 50 == 0 || iter == 1) {
#       cat(
#         "Iter:", sprintf("%4d", iter),
#         "| Obj:", sprintf("%10.5f", res$objective_value),
#         "| Best RC:", sprintf("%10.5f", best_rc),
#         "| Cols:", ncol(A_master), "\n"
#       )
#     }
#   }
# 
#   if (!optimal_found) {
#     warning("Reached max iterations without proving optimality.")
#   }
# 
#   run_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
#   final_obj <- res$objective_value
#   if (is_upper_bound) final_obj <- -final_obj
# 
#   list(
#     a = target_a,
#     b = target_b,
#     lower = if (is_upper_bound) NA_real_ else final_obj,
#     upper = if (is_upper_bound) final_obj else NA_real_,
#     objective = final_obj,
#     iterations = iter,
#     time_secs = run_time,
#     final_cols = ncol(A_master),
#     status = res$status,
#     status_message = res$status_message
#   )
# }
# 
# # ------------------------------------------------------------
# # 5. Run all pairwise ACE bounds
# # ------------------------------------------------------------
# pairs <- combn(0:(k - 1), 2, simplify = FALSE)
# 
# pairwise_results <- vector("list", length(pairs))
# 
# for (p in seq_along(pairs)) {
#   a <- pairs[[p]][1]
#   b0 <- pairs[[p]][2]
# 
#   cat("\n\n########## Pair", a, "vs", b0, "##########\n")
# 
#   # CORRECTED: Passed directly in order (a -> target_a, b0 -> target_b)
#   res_lower <- run_column_generation(target_a = a, target_b = b0, is_upper_bound = FALSE)
#   res_upper <- run_column_generation(target_a = a, target_b = b0, is_upper_bound = TRUE)
# 
#   pairwise_results[[p]] <- data.frame(
#     a = a,
#     b = b0,
#     lower = res_lower$objective,
#     upper = res_upper$objective,
#     lower_iter = res_lower$iterations,
#     upper_iter = res_upper$iterations,
#     lower_time = res_lower$time_secs,
#     upper_time = res_upper$time_secs,
#     lower_cols = res_lower$final_cols,
#     upper_cols = res_upper$final_cols
# 
#   )
# 
# 
# }
# 
# pairwise_df <- bind_rows(pairwise_results)
# 
# cat("\n==================================================\n")
# cat("PAIRWISE ACE BOUNDS\n")
# print(pairwise_df)
# cat("==================================================\n")
# 
# # ------------------------------------------------------------
# # 6. Save results
# # ------------------------------------------------------------
# saveRDS(
#   list(
#     k = k,
#     m = m,
#     pairwise_df = pairwise_df,
#     prob_map = prob_map,
#     obs_grid = obs_grid
#   ),
#   file.path(run_dir, "pairwise_bounds_results.rds")
# )
# 
# write.csv(
#   pairwise_df,
#   file.path(run_dir, "pairwise_bounds_summary.csv"),
#   row.names = FALSE
# )
# 
