# Fits bayome to the HMP2 example data and saves the results used by the
# vignette (vignettes/bayome.Rmd). The fit takes several minutes, so the
# vignette reads these saved results instead of fitting the model.
#
# Run from the package root (the folder with DESCRIPTION), either
#   Rscript data-raw/hmp2_example.R
# in a terminal, or source("data-raw/hmp2_example.R") in R with the working
# directory set to the package root. It uses the development version of the
# package (devtools::load_all) when devtools is installed.
#
# Writes, in inst/extdata/:
#   HMP2_bayome_results.tsv  results table (fit$results)
#   HMP2_bayome_lambda.tsv   estimated prevalence-abundance links and run time
# and data-raw/hmp2_example_info.txt (not versioned): details of the run.

if (!file.exists("DESCRIPTION") || !file.exists("data-raw/hmp2_example.R")) {
  stop("Run this script from the bayome package root.", call. = FALSE)
}
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(".") else library(bayome)

taxa <- read.delim("inst/extdata/HMP2_taxonomy.tsv", row.names = 1)
meta <- read.delim("inst/extdata/HMP2_metadata.tsv", row.names = 1)
meta$diagnosis <- factor(meta$diagnosis, levels = c("nonIBD", "UC", "CD"))
meta$antibiotics <- factor(meta$antibiotics, levels = c("No", "Yes"))

set.seed(1)
t0 <- Sys.time()
fit <- bayome(
  input_data     = taxa,
  input_metadata = meta,
  output         = NULL,
  formula        = ~ diagnosis + antibiotics + age + reads + (1 | participant_id),
  exposure       = "diagnosis",
  normalization  = "TSS",
  transform      = "LOG",
  save_models    = TRUE,
  cores          = 4
)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
print(fit)

if (any(grepl("INLA error", fit$results$error))) {
  stop("The fit failed: ", unique(grep("INLA error", fit$results$error, value = TRUE)))
}
lambda <- fit$lambda
vb <- attr(fit$fits$joint, "bm_vb")
info <- c(
  sprintf("date: %s", format(t0)),
  sprintf("R %s, INLA %s", getRversion(), utils::packageVersion("INLA")),
  sprintf("total minutes: %.1f", minutes),
  sprintf("vb iterations: %d", vb$iterations),
  sprintf("vb changes: %s", paste(signif(vb$changes, 3), collapse = " ")),
  sprintf("seconds per INLA fit: %s",
          paste(names(vb$seconds), round(vb$seconds), sep = "=", collapse = " ")),
  sprintf("hyperparameters: %s",
          paste(rownames(fit$fits$joint$summary.hyperpar), collapse = "; ")),
  sprintf("lambda: %s", paste(names(lambda), round(lambda, 3), collapse = ", ")),
  sprintf("residual variances: d0 = %.2f, s0^2 = %.3f", vb$d0, vb$s02))
writeLines(info, "data-raw/hmp2_example_info.txt")
cat(info, sep = "\n")
write.table(fit$results, "inst/extdata/HMP2_bayome_results.tsv", sep = "\t",
            quote = TRUE, row.names = FALSE)
write.table(data.frame(name = names(lambda), lambda = unname(lambda),
                       minutes = round(minutes, 1)),
            "inst/extdata/HMP2_bayome_lambda.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)
cat(sprintf("Done in %.1f minutes\n", minutes))
