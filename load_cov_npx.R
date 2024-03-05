# Function for loading COVID NPX data
load_cov_npx <- function(dataset = c("trimmed", "bridge_norm", "abspqn", "abspqn_excl", "abspqn_adj", "ds3_trimmed", "ds3_abspqn", "ds3_abspqn_excl", "ds3_abspqn_adj")) {
  # Load specified dataset, default is first element of input vector (given as vector just to show options when using function)
  # "bridge_norm", Bridge normalised (BN) data (datasets 1 and 2)
  # "trimmed", BN with binders with high missingness trimmed out
  # "abspqn", trimmed with AbsPQN applied per protein panel and dataset
  # "abspqn_excl", abspqn with outliers excluded
  # "abspqn_adj", abspqn_excl adjusted for age and sex
  # "ds3", dataset 3 trimmed from high-missingness binders
  # "ds3_abspqn", trimmed dataset 3 with AbsPQN applied
  # "ds3_abspqn_excl", dataset 3 with AbsPQN and outlier removal
  # "ds3_abspqn_adj", dataset 3 with AbsPQN, outliers removed, and adjusted for age and sex
  
  data_dir <- "../data_dir/"
  data_file <- switch(dataset[1],
                      "bridge_norm" = "cov_npx_bridge_norm.rds",
                      "trimmed" = "cov_npx_trimmed.rds",
                      "abspqn" = "cov_npx_abspqn.rds",
                      "abspqn_excl" = "cov_npx_abspqn_excl.rds",
                      "abspqn_adj" = "cov_npx_abspqn_adj.rds",
                      "ds3_trimmed" = "cov_npx_ds3.rds",
                      "ds3_abspqn" = "cov_npx_ds3_abspqn.rds",
                      "ds3_abspqn_excl" = "cov_npx_ds3_abspqn_excl.rds",
                      "ds3_abspqn_adj" = "cov_npx_ds3_abspqn_adj.rds")
  
  data_path <- paste0(data_dir, data_file)
  
  data_out <- readRDS(data_path)
  data_out$data_set <- dataset[1]
  data_out$data_path <- data_path
  
  return(data_out)
}