# 
# # install.packages("igraph")
# 
# library(causaloptim)
# library(igraph)
# 
# estimate_bounds_iv_5level <- function() {
# 
#   timing <- system.time({
# 
#     g <- graph_from_literal(
#       Ul -+ Z,
#       Z  -+ X,
#       X  -+ Y,
#       Ur -+ X,
#       Ur -+ Y
#     )
# 
#     # Add the attributes expected by causaloptim
#     g <- initialize_graph(g)
# 
#     V(g)$leftside <- as.integer(V(g)$name %in% c("Ul", "Z"))
#     V(g)$latent   <- as.integer(V(g)$name %in% c("Ul", "Ur"))
# 
#     # 5-level Z and X, binary Y; latent nodes are kept binary here
#     nvals_map <- c(Ul = 2L, Z = 5L, X = 5L, Y = 2L, Ur = 2L)
#     V(g)$nvals <- unname(nvals_map[as.character(V(g)$name)])
# 
#     # No right-to-left edges
#     E(g)$rlconnect <- 0L
# 
#     # Monotonicity is only available for binary edges, so leave off here
#     E(g)$edge.monotone <- 0L
# 
#     # Effect: risk difference between X = 4 and X = 0
#     effectt <- "p{Y(X = 2) = 1} - p{Y(X = 0) = 1}"
# 
#     # No extra constraints
#     obj <- analyze_graph(g, constraints = NULL, effectt = effectt)
# 
#     # Compute symbolic bounds
#     bounds <- optimize_effect_2(obj)
# 
#     # Optional: convert the symbolic bounds to a callable function
#     # bounds_fun <- interpret_bounds(bounds$bounds, obj$parameters)
# 
# 
#   })
# 
#   list(
#     graph = g,
#     objective = obj,
#     bounds = bounds,
#     timing = timing
#   )
# }
# 
# res <- estimate_bounds_iv_5level()
# 
# print(res$timing)
# print(res$bounds)
# 
# 
