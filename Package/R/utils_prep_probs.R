#' Prepare Empirical Probabilities
#' @importFrom dplyr mutate count group_by ungroup left_join
#' @importFrom tidyr complete
#' @importFrom stats setNames
prep_empirical_probs <- function(data, outcome, treatment, group, m = NULL, k = NULL, l = NULL) {
  recode_0_based <- function(v, expected_n = NULL, var_name = "variable") {
    lev <- sort(unique(v))
    lev <- lev[!is.na(lev)]
    
    if (!is.null(expected_n) && length(lev) != expected_n) {
      stop(
        sprintf(
          "%s has %d unique values, but expected %d.",
          var_name, length(lev), expected_n
        ),
        call. = FALSE
      )
    }
    
    out <- match(v, lev) - 1L
    if (anyNA(out)) {
      stop(sprintf("%s contains NA or unmapped values.", var_name), call. = FALSE)
    }
    out
  }
  
  # Adjusted notation: m = treatment, k = group/instrument
  if (is.null(m)) m <- length(unique(data[[treatment]]))
  if (is.null(k)) k <- length(unique(data[[group]]))
  if (is.null(l)) l <- length(unique(data[[outcome]])) 
  
  dat0 <- data %>%
    dplyr::mutate(
      X = recode_0_based(.data[[treatment]], m, treatment),
      Y = recode_0_based(.data[[outcome]], l, outcome), 
      G = recode_0_based(.data[[group]], k, group)
    )
  
  prob_df <- dat0 %>%
    dplyr::count(.data$G, .data$X, .data$Y, name = "n") %>%
    tidyr::complete(
      G = 0:(k - 1),
      X = 0:(m - 1),
      Y = 0:(l - 1), 
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(.data$G) %>%
    dplyr::mutate(prob = if (sum(.data$n) == 0) NA_real_ else .data$n / sum(.data$n)) %>%
    dplyr::ungroup()
  
  obs_grid <- expand.grid(
    G = 0:(k - 1),
    X = 0:(m - 1),
    Y = 0:(l - 1), 
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  prob_map <- obs_grid %>%
    dplyr::left_join(prob_df, by = c("G", "X", "Y")) %>%
    dplyr::mutate(prob = ifelse(is.na(.data$prob), 0, .data$prob))
  
  num_obs <- nrow(obs_grid)
  b <- c(prob_map$prob, 1)
  
  # Adjusted dimensions: k (group), m (treatment), l (outcome)
  obs_idx_array <- array(NA_integer_, dim = c(k, m, l)) 
  for (r in seq_len(num_obs)) {
    g_idx <- obs_grid$G[r] + 1
    x_idx <- obs_grid$X[r] + 1
    y_idx <- obs_grid$Y[r] + 1
    obs_idx_array[g_idx, x_idx, y_idx] <- r
  }
  
  block_sums <- prob_df %>%
    dplyr::group_by(.data$G) %>%
    dplyr::summarise(block_sum = sum(.data$prob), .groups = "drop")
  
  list(
    b = b,
    num_obs = num_obs,
    obs_grid = obs_grid,
    obs_idx_array = obs_idx_array,
    m = m, # Treatment levels
    k = k, # Group levels
    l = l, # Outcome levels
    prob_df = prob_df
  )
}