# Power simulation: generate one dataset, train models, run associations,
# record the result, and then proceed to the next dataset.

library(data.table)
library(foreach)
library(doParallel)
library(doRNG)
library(glmnet)
library(ctOWAS)
library(MiXcan)
library(ggplot2)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")

paper_dir <- getwd()
output_dir <- file.path(paper_dir, "Results", "Simulation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source('Github/Analysis_code/Figure4_Simulation/4.2_simulation_function.R')
X_pool <- data.table::fread(file.path(paper_dir, "Data", "Simulation", "chr22_20_21Mb_1million_pseudo_subjects.csv"), data.table = FALSE)


# Simulation settings ---------------------------------------------------------

n_replicates <- as.integer(Sys.getenv("ctOWAS_POWER_N_REPLICATES", unset = "500"))
workers <- as.integer(Sys.getenv("ctOWAS_POWER_WORKERS", unset = "5"))
n_train <- as.integer(Sys.getenv("ctOWAS_POWER_N_TRAIN", unset = "300"))
n_test <- as.integer(Sys.getenv("ctOWAS_POWER_N_TEST", unset = "3000"))
reg_scale <- as.numeric(Sys.getenv("ctOWAS_POWER_REG_SCALE", unset = "0.05"))

b0 <- 1
b1 <- 1
b2 <- 1
group <- "heter xy"

scenarios <- data.table(
  scenario_id = c("S2", "S3", "S4", "S5"),
  scenario = c("eta1 = 0.2, eta2 = 0.2", "eta1 = 0.2, eta2 = 0",
               "eta1 = 0, eta2 = 0.2", "eta1 = -0.2, eta2 = 0.2"),
  eta1 = c(0.2, 0.2, 0, -0.2),
  eta2 = c(0.2, 0, 0.2, 0.2)
)


# Run one complete dataset at a time ------------------------------------------

cl <- parallel::makeCluster(min(workers, nrow(scenarios)))
doParallel::registerDoParallel(cl)
# foreach::registerDoSEQ()
all_results <- tryCatch(
  foreach(
    scenario_index = seq_len(nrow(scenarios)),
    .packages = c("data.table", "glmnet", "ctOWAS", "MiXcan"),
    .export = c("generate_one_dataset", "run_one_replication"),
    .options.RNG = 2026081,
    .errorhandling = "stop"
  ) %dorng% {
    scenario <- scenarios[scenario_index]
    scenario_results <- vector("list", n_replicates)

    for (replicate_id in seq_len(n_replicates)) {
      message(scenario$scenario_id, ": replicate ", replicate_id, "/", n_replicates)

      # This call completes data generation, prediction, and association before
      # the loop advances to the next replicate.
      scenario_results[[replicate_id]] <- tryCatch(
        run_one_replication(
          replicate_id = replicate_id, scenario_index = scenario_index, scenario = scenario,
          X_pool = X_pool, n_train = n_train, n_test = n_test,
          b0 = b0, b1 = b1, b2 = b2, group = group,
          nfolds = nfolds, prediction_alpha = prediction_alpha,
          regularization = "fixed", reg_scale = reg_scale, condition_max = 100,
          reg_min_scale = 0.001, reg_max_scale = 0.1, ld_reference = "test"
        ),
        error = function(e) data.table(
          scenario_id = scenario$scenario_id, scenario = scenario$scenario,
          eta1 = scenario$eta1, eta2 = scenario$eta2, replicate = replicate_id,
          error = conditionMessage(e)
        )
      )
    }
    rbindlist(scenario_results, use.names = TRUE, fill = TRUE)
  }, finally = {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  })

full_results <- rbindlist(all_results, use.names = TRUE, fill = TRUE)

fwrite(full_results, file.path(output_dir, "power_fixed_setting_full_results.csv"))

