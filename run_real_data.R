# =============================================================================
# run_real_data.R
#
# Template for applying the model to a real cohort. Nothing here is specific
# to any dataset; supply a CSV with one row per participant and the columns
# named below.
#
#   Rscript examples/run_real_data.R path/to/cohort.csv
#
# Required columns
#   fa         mean fractional anisotropy per participant
#   burden     vascular burden score (any scale; it is standardised here)
#   cognition  cognitive outcome (composite or single measure)
#   age, sex, edu    covariates
#
# Extracting one FA value per participant is done upstream, for example with
# FSL: topup and eddy for distortion and motion correction, BET for brain
# extraction, DTIFIT for tensor fitting, then TBSS to obtain a skeletonised
# mean FA. See README for the full pipeline.
# =============================================================================

source("R/model.R")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript examples/run_real_data.R <cohort.csv>")

raw <- utils::read.csv(args[1], stringsAsFactors = FALSE)
cat(sprintf("Loaded %d rows from %s\n", nrow(raw), args[1]))

d   <- prepare_data(raw,
                    fa         = "fa",
                    burden     = "burden",
                    cognition  = "cognition",
                    covariates = c("age", "sex", "edu"))
fit <- fit_wm_model(d)

cat("\n--- Sample description ---\n")
cat(sprintf("N analysed   : %d\n", nrow(d)))
cat(sprintf("Age          : %.1f (SD %.1f)\n", mean(raw$age, na.rm = TRUE),
            sd(raw$age, na.rm = TRUE)))
cat(sprintf("Mean FA      : %.3f (SD %.3f)\n", mean(raw$fa, na.rm = TRUE),
            sd(raw$fa, na.rm = TRUE)))
cat(sprintf("Burden       : %.2f (SD %.2f)\n", mean(raw$burden, na.rm = TRUE),
            sd(raw$burden, na.rm = TRUE)))

cat("\n--- Fitted model ---\n")
print(format(model_table(fit), digits = 3))
cat(sprintf("\nR-squared = %.3f, adjusted = %.3f\n",
            summary(fit)$r.squared, summary(fit)$adj.r.squared))

cat("\n--- Simple slopes ---\n")
print(format(simple_slopes(fit), digits = 3))

check_assumptions(fit)

if (!dir.exists("figures")) dir.create("figures")
plot_interaction(fit, file = "figures/fig_interaction_real.png")

utils::write.csv(model_table(fit), "model_coefficients.csv", row.names = FALSE)
cat("\nCoefficients written to model_coefficients.csv\n")
cat("Figure written to figures/fig_interaction_real.png\n")
