test_data <- data.frame(
  Y = sample(c(0, 1), 500, replace = TRUE),
  X = sample(1:4, 500, replace = TRUE),   # k = 10 treatment options
  G = sample(0:5, 500, replace = TRUE)     # m = 5 groupings
)

# 4. Run the user-facing package wrapper
results <- compute_ace_bounds(
  data = test_data, 
  outcome = "Y", 
  treatment = "X", 
  group = "G"
)

# 5. Check out the output matrix
print(results)

# Needs some serious parallelisation 
# consider parallelising every pair of ACE 