#' Compute ACE Bounds
#' 
#' @param data A data frame containing variables for analysis.
#' @param outcome Character string naming the outcome column (multi-state or binary).
#' @param treatment Character string naming the treatment column (values 0 to m-1).
#' @param group Character string naming the grouping or instrument column (values 0 to k-1).
#' @param assumptions Character vector specifying shape assumptions. Options include "MTR". Default is NULL.
#' @param m Integer. Number of treatment levels. If NULL, inferred automatically.
#' @param k Integer. Number of group levels. If NULL, inferred automatically.
#' @param delta Numeric scalar for managing finite sample variance. Default is 0.001. 
#' @param epsilon Numeric scalar for convergence verification. Default is 1e-6. 
#' @param max_iter Integer scalar for safety limit iteration. Default is 2000.
#' @param verbose Logical. If TRUE, prints iteration progress to the console. Default is FALSE.
#' @param all_pairs Logical. If TRUE, computes bounds for all treatment combinations. If FALSE, computes only lowest (0) vs highest (m-1). Default is FALSE.
#' 
#' @return A data frame containing computed lower and upper bounds for the selected treatment pairs.
#' @export
compute_ace_bounds <- function(data, outcome, treatment, group, 
                               assumptions = NULL, m = NULL, k = NULL, 
                               delta = 0.001, epsilon = 1e-6, max_iter = 2000,
                               verbose = TRUE, all_pairs = FALSE) {
  
  if (verbose) cat("Preparing empirical probabilities...\n")
  prep_data <- prep_empirical_probs(data, outcome, treatment, group, m, k)
  
  if (verbose) {
    assump_text <- if (is.null(assumptions)) "no assumptions (full space)" else paste(assumptions, collapse=", ")
    cat("Building latent response space under", assump_text, "...\n")
  }
  
  # Adjusted notation: The latent space build now relies on m (treatment) and l (outcome)
  types_Y <- build_latent_space(
    m = prep_data$m, 
    l = prep_data$l, 
    assumptions = assumptions
  )
  
  # Adjusted logic: Iterating over treatment pairs (m)
  if (all_pairs) {
    pairs <- utils::combn(0:(prep_data$m - 1), 2, simplify = FALSE)
  } else {
    pairs <- list(c(0, prep_data$m - 1))
  }
  
  pairwise_results <- vector("list", length(pairs))
  
  for (p in seq_along(pairs)) {
    a <- pairs[[p]][1]
    b0 <- pairs[[p]][2]
    
    if (verbose) {
      cat("\n\n########## Pair", a, "vs", b0, "##########\n")
    }
    
    res_lower <- solve_ace_bound(a, b0, FALSE, delta, epsilon, max_iter, types_Y, prep_data, verbose)
    res_upper <- solve_ace_bound(a, b0, TRUE,  delta, epsilon, max_iter, types_Y, prep_data, verbose)
    
    pairwise_results[[p]] <- data.frame(
      a = a, b = b0,
      lower = res_lower$objective, upper = res_upper$objective,
      lower_iter = res_lower$iterations, upper_iter = res_upper$iterations,
      lower_time = res_lower$time_secs, upper_time = res_upper$time_secs,
      lower_cols = res_lower$final_cols, upper_cols = res_upper$final_cols
    )
  }
  
  pairwise_df <- dplyr::bind_rows(pairwise_results)
  
  if (verbose) {
    cat("\n==================================================\n")
    cat(ifelse(all_pairs, "ALL PAIRWISE ACE BOUNDS\n", "EXTREME PAIR ACE BOUNDS (0 vs m-1)\n"))
    print(pairwise_df)
    cat("==================================================\n")
  }
  
  return(pairwise_df)
}