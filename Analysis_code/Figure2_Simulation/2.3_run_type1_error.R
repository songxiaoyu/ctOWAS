library(data.table);
library(foreach);
library(doParallel);
library(doRNG);
library(glmnet);
library(ctOWAS);
library(MiXcan)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")

paper_dir <- getwd()
output_dir <- file.path(paper_dir, "Results", "Simulation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("Github/Analysis_code/Figure2_Simulation/2.2_simulation_function.R")
X_pool <- data.table::fread(file.path(paper_dir, "Data", "Simulation", "chr22_20_21Mb_1million_pseudo_subjects.csv"), data.table = FALSE)

# Simulation settings ---------------------------------------------------------

n_replicates <- 2000
workers <- 15                  # number of CPU cores
n_train <- 300
n_test <- 3000
reg_scale <- 0.05

b0 <- 1
b1 <- 1
b2 <- 1
group <- "heter xy"

scenario <- data.table(scenario_id = "S0", scenario = "eta1 = 0, eta2 = 0", eta1 = 0, eta2 = 0)

# Parallel simulation across replicates ---------------------------------------

workers <- min(workers, n_replicates)
cl <- parallel::makeCluster(workers)
doParallel::registerDoParallel(cl)

results_list <- tryCatch(
  foreach(replicate_id = seq_len(n_replicates),
          .packages = c("data.table", "glmnet", "ctOWAS", "MiXcan"),
          .export = c("generate_one_dataset", "run_one_replication"),
          .options.RNG = 2026082, .errorhandling = "pass") %dorng% {

            tryCatch(
              run_one_replication(replicate_id = replicate_id, scenario_index = 1, scenario = scenario, X_pool = X_pool,
                                  n_train = n_train, n_test = n_test, b0 = b0, b1 = b1, b2 = b2, group = group,
                                  nfolds = nfolds, prediction_alpha = prediction_alpha, regularization = "fixed", reg_scale = 0.05,
                                  condition_max = 100, reg_min_scale = 0.001, reg_max_scale = 0.1, ld_reference = "test"),
              error = function(e) data.table(scenario_id = scenario$scenario_id, scenario = scenario$scenario, eta1 = scenario$eta1,
                                             eta2 = scenario$eta2, replicate = replicate_id, error = conditionMessage(e))
            )
          },
  finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() }
)

# Combine and save ------------------------------------------------------------

full_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
fwrite(full_results, file.path(output_dir, "type1_fixed_setting_full_results.csv"))
