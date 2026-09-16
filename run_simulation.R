# =============================================================================
# run_simulation.R
#
# Reproduces the full simulation study. From the repository root:
#
#   Rscript examples/run_simulation.R
#
# All results below come from synthetic data generated under assumed
# parameters. They are not observations about any participants.
# =============================================================================

source("R/model.R")
source("R/simulate.R")

set.seed(2026)
if (!dir.exists("figures")) dir.create("figures")

P <- default_params()

cat("\n==========================================================\n")
cat(" SIMULATION STUDY: FA x BURDEN MODERATION MODEL\n")
cat(" All data below are synthetic.\n")
cat("==========================================================\n")


# --- 1. A single example cohort --------------------------------------------

demo <- simulate_cohort(200, P)
d    <- prepare_data(demo, covariates = c("age", "sex", "edu"))
fit  <- fit_wm_model(d)

cat("\n--- Example synthetic cohort (N = 200) ---\n")
cat(sprintf("Age          : %.1f (SD %.1f), range %.0f-%.0f\n",
            mean(demo$age), sd(demo$age), min(demo$age), max(demo$age)))
cat(sprintf("Sex          : %d F / %d M\n",
            sum(demo$sex == 0), sum(demo$sex == 1)))
cat(sprintf("Mean FA      : %.3f (SD %.3f)\n", mean(demo$fa), sd(demo$fa)))
cat(sprintf("Burden score : %.2f (SD %.2f)\n", mean(demo$burden), sd(demo$burden)))
cat(sprintf("FA-burden r  : %.3f\n", cor(demo$fa, demo$burden)))

cat("\n--- Fitted model ---\n")
print(format(model_table(fit), digits = 3))
cat(sprintf("\nR-squared = %.3f, adjusted = %.3f\n",
            summary(fit)$r.squared, summary(fit)$adj.r.squared))

cat("\n--- Simple slopes: FA effect at each burden level ---\n")
print(format(simple_slopes(fit), digits = 3))

check_assumptions(fit)


# --- 2. Parameter recovery --------------------------------------------------

cat("\n--- Parameter recovery (1000 replications, N = 200) ---\n")
rec <- recover_parameters(n = 200, nsim = 1000, p = P)
print(format(rec, digits = 3))
cat("\nBias near zero and coverage near 0.95 indicate the estimator is\n")
cat("unbiased and the intervals are correctly calibrated.\n")


# --- 3. Power ---------------------------------------------------------------

cat("\n--- Power across sample sizes ---\n")
pwr <- power_curve(seq(50, 600, by = 50), nsim = 500, p = P)
print(format(pwr, digits = 3))

n80 <- pwr$n[which(pwr$power_interaction >= 0.80)[1]]
nfa <- pwr$n[which(pwr$power_fa >= 0.80)[1]]
if (!is.na(n80)) cat(sprintf("\nInteraction reaches 80%% power at about N = %d.\n", n80))
if (!is.na(nfa)) cat(sprintf("FA main effect reaches 80%% power at about N = %d.\n", nfa))

cat("\n--- Required N across plausible interaction sizes ---\n")
req <- required_n()
print(format(req, digits = 3))
cat("\nThe required sample depends strongly on the assumed b3, which is not\n")
cat("known in advance. Report the range rather than a single figure.\n")


# --- 4. Measurement error ---------------------------------------------------

cat("\n--- Effect of FA measurement error ---\n")
att <- attenuation()
print(format(att, digits = 3))
cat("\nLower FA reliability attenuates b3 toward zero and reduces power.\n")


# --- 5. Figures -------------------------------------------------------------

big <- prepare_data(simulate_cohort(5000, P), covariates = c("age", "sex", "edu"))
plot_interaction(fit_wm_model(big), file = "figures/fig_interaction.png")
plot_power(pwr, file = "figures/fig_power.png")

cat("\n----------------------------------------------------------\n")
cat(" Figures written to figures/\n")
cat("----------------------------------------------------------\n\n")
