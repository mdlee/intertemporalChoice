#!/usr/bin/env Rscript
# Sequential CmdStan fits of the eight entropified hierarchical execution models.
# Same Grünwald overlay as the JAGS _entrop files: hard step() rule, then
# P(LL|rule LL)=inv_logit(w), P(LL|rule SS)=inv_logit(-w).
#
# Usage:
#   Rscript run_hierarchical_stan.R
#   Rscript run_hierarchical_stan.R --scale=0.5 --tag=MuPrecHalf
#   Rscript run_hierarchical_stan.R --scale=2 --tag=MuPrecDouble
#
# Stop: create models/stan/STOP_STAN_HIERARCHICAL
# Status: models/stan/logs/run_hierarchical_stan.status

suppressPackageStartupMessages({
  library(cmdstanr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", flag, "="), "", hit[[1]])
}
mu_prec_scale <- as.numeric(parse_arg("--scale", "1"))
tag <- parse_arg("--tag", "baseline")
n_chains <- as.integer(parse_arg("--chains", "6"))
n_warmup <- as.integer(parse_arg("--warmup", "1000"))
n_sampling <- as.integer(parse_arg("--sampling", "2000"))
seed <- as.integer(parse_arg("--seed", "1"))

stan_dir <- if (sys.nframe() == 0) {
  normalizePath(getwd())
} else {
  normalizePath(dirname(sys.frame(1)$ofile))
}
# When run via Rscript, ofile is unset; locate from --file=
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg)) {
  stan_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
}
root <- normalizePath(file.path(stan_dir, "..", ".."))
data_dir <- file.path(root, "data")
out_dir <- file.path(stan_dir, "fits", tag)
log_dir <- file.path(stan_dir, "logs")
include_dir <- file.path(stan_dir, "include")
stop_file <- file.path(stan_dir, "STOP_STAN_HIERARCHICAL")
status_file <- file.path(log_dir, "run_hierarchical_stan.status")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

write_status <- function(msg) {
  line <- sprintf("%s\n%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  writeLines(line, status_file)
  cat(line)
}

source(file.path(stan_dir, "prepare_stan_data.R"))
stan_data <- prepare_stan_data(data_dir, mu_prec_scale)
jsonlite::write_json(
  stan_data,
  file.path(stan_dir, "data", sprintf("intertemporalChoice_scale%s.json", mu_prec_scale)),
  auto_unbox = TRUE, pretty = FALSE, digits = NA
)

models <- c(
  "exponentialExecutionHierarchical_entrop",
  "hyperbolicExecutionHierarchical_entrop",
  "proportionalDifferencesExecutionHierarchical_entrop",
  "directDifferencesExecutionHierarchical_entrop",
  "tradeoffExecutionHierarchical_entrop",
  "unifiedTradeoffExecutionHierarchical_entrop",
  "itchExecutionHierarchical_entrop",
  "hyperboloidExecutionHierarchical_entrop"
)

write_status(sprintf(
  "starting | tag=%s | muPrecScale=%g | %d models | %d chains | warmup=%d sampling=%d",
  tag, mu_prec_scale, length(models), n_chains, n_warmup, n_sampling
))

max_rhat <- function(fit) {
  s <- fit$summary()
  r <- s$rhat
  r <- r[is.finite(r)]
  if (!length(r)) return(NA_real_)
  max(r)
}

for (k in seq_along(models)) {
  if (file.exists(stop_file)) {
    write_status(sprintf("stop file present — exiting before %s", models[[k]]))
    quit(save = "no", status = 0)
  }
  name <- models[[k]]
  stan_file <- file.path(stan_dir, paste0(name, ".stan"))
  write_status(sprintf("[%d/%d] %s | tag=%s | compiling+sampling", k, length(models), name, tag))
  cat(sprintf("\n========== [%d/%d] %s (%s) ==========\n", k, length(models), name, tag))
  mod <- cmdstan_model(
    stan_file,
    include_paths = include_dir,
    dir = file.path(stan_dir, "output")
  )
  fit <- mod$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = seed + k,
    refresh = 100,
    adapt_delta = 0.9,
    max_treedepth = 12,
    output_dir = out_dir,
    output_basename = name
  )
  rhat <- max_rhat(fit)
  rds <- file.path(out_dir, paste0(name, ".rds"))
  fit$save_object(rds)
  summ_csv <- file.path(out_dir, paste0(name, "_summary.csv"))
  utils::write.csv(fit$summary(), summ_csv, row.names = FALSE)
  write_status(sprintf(
    "[%d/%d] %s | done | max R-hat=%.3f | n=%d | chains=%d",
    k, length(models), name, rhat, n_sampling, n_chains
  ))
  cat(sprintf("Saved %s | max R-hat=%.3f\n", rds, rhat))
}

write_status(sprintf("finished at %s | tag=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tag))
cat("\n=== run_hierarchical_stan finished ===\n")
