#' Build Latent Response Types Space
#' @keywords internal
#' @importFrom RcppAlgos permuteGeneral comboGeneral
build_latent_space <- function(m, l = NULL, assumptions = NULL) {
  
  if (is.null(l)) {
    l <- 2
  }
  
  # Define the state space (e.g., 0, 1 for l=2)
  states <- 0:(l - 1)
  
  # If MTR is assumed, generate only valid monotonically non-decreasing rows
  if (!is.null(assumptions) && "MTR" %in% assumptions) {
    # Mathematically equivalent to combinations with replacement
    # Note: RcppAlgos uses 'm' for the combination size, so we use m = m
    types_Y <- RcppAlgos::comboGeneral(v = states, m = m, repetition = TRUE)
    
  } else {
    # Build full latent space (cartesian product)
    # Mathematically equivalent to permutations with replacement
    types_Y <- RcppAlgos::permuteGeneral(v = states, m = m, repetition = TRUE)
  }
  
  # Ensure the matrix is stored as integers to save memory
  storage.mode(types_Y) <- "integer"
  
  # Assign column names: X0, X1, ..., X{m-1}
  colnames(types_Y) <- paste0("X", 0:(m - 1))
  
  types_Y
}