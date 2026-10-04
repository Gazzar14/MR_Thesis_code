price_chunk <- function(
    rows,
    types_Y_mat,
    c_Y,
    pi_duals,
    pi_norm,
    obs_idx_array,
    k,
    m
) {
  
  best_rc <- Inf
  best_col_indices <- NULL
  best_c <- NULL
  
  for (i in rows) {
    
    current_Y <- types_Y_mat[i, ]
    
    current_c <- c_Y[i]
    
    pi_sum <- 0
    current_indices <- integer(m)
    
    for (g_val in 0:(m - 1)) {
      
      best_pi <- -Inf
      best_idx <- NA_integer_
      
      for (x_val in 0:(k - 1)) {
        
        y_val <- current_Y[x_val + 1]
        
        idx <- obs_idx_array[
          g_val + 1,
          x_val + 1,
          y_val + 1
        ]
        
        p_val <- if (!is.na(idx))
          pi_duals[idx]
        else
          0
        
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
  
  list(
    rc = best_rc,
    col_indices = best_col_indices,
    c = best_c
  )
}