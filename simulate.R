# =============================================================================
# simulate.R
#
# Data generation and power analysis for the moderated model.
#
# IMPORTANT: everything produced here is synthetic. The parameter values are
# assumptions chosen to be plausible, not estimates from any cohort. Output
# from these functions must be reported as simulation, never as empirical
# findings about participants.
# =============================================================================


#' Default generating parameters
#'
#' Effect sizes are standardised: a 1 SD change in the predictor produces this
#' many SD of change in cognition.
#'
#' @return named list of parameters
default_params <- function() {
  list(
    b0        =  0.00,   # intercept (cognition z-scored)
    b1        =  0.30,   # FA -> cognition
    b2        = -0.25,   # burden -> cognition
    b3        = -0.20,   # FA x burden interaction (the double-hit term)
    b_age     = -0.35,   # age -> cognition
    b_sex     =  0.05,   # sex (M vs F)
    b_edu     =  0.20,   # education -> cognition

    age_mean  = 65,      # target an ageing cohort
    age_sd    =  8,
    fa_mean   =  0.45,   # plausible whole-brain skeleton mean FA
    fa_sd     =  0.035,
    fa_age    = -0.0015, # FA lost per year of age
    fa_burden = -0.060,  # FA lost across the full 0-1 burden range
    resid_sd  =  0.80    # unexplained variation in cognition
  )
}


#' Simulate a synthetic cohort
#'
#' The causal ordering matters. Age and burden are generated first, FA is
#' generated partly as a consequence of them, and cognition from all three.
#' This induces realistic correlation between FA and burden, so the power
#' estimates that follow are not optimistic.
#'
#' @param n number of participants
#' @param p parameter list, see default_params()
#' @return data frame with raw and standardised variables
simulate_cohort <- function(n, p = default_params()) {

  age <- stats::rnorm(n, p$age_mean, p$age_sd)
  sex <- stats::rbinom(n, 1, 0.5)          # 0 = F, 1 = M
  edu <- stats::rnorm(n, 0, 1)             # years of education, z-scored

  # Burden: four binary risk factors scaled to 0-1, more likely at older ages
  risk_p <- stats::plogis(-1.2 + 0.04 * (age - p$age_mean))
  burden <- rowSums(matrix(stats::rbinom(n * 4, 1, rep(risk_p, 4)), nrow = n)) / 4

  fa <- p$fa_mean +
        p$fa_age * (age - p$age_mean) +
        p$fa_burden * burden +
        stats::rnorm(n, 0, p$fa_sd)

  fa_z     <- as.numeric(scale(fa))
  burden_z <- as.numeric(scale(burden))
  age_z    <- as.numeric(scale(age))

  cognition <- p$b0 +
               p$b1 * fa_z +
               p$b2 * burden_z +
               p$b3 * fa_z * burden_z +
               p$b_age * age_z +
               p$b_sex * sex +
               p$b_edu * edu +
               stats::rnorm(n, 0, p$resid_sd)

  data.frame(age, sex, edu, burden, fa, fa_z, burden_z, age_z, cognition)
}


# Internal: fit on an already-standardised simulated cohort
.fit_sim <- function(d) {
  stats::lm(cognition ~ fa_z * burden_z + age_z + sex + edu, data = d)
}


#' Parameter recovery check
#'
#' Repeatedly generates data under known coefficients and refits the model.
#' Bias near zero and coverage near the nominal level indicate the model is
#' correctly specified and the estimator behaves as intended. This validates
#' the code, not any hypothesis about brains.
#'
#' @param n    sample size per replication
#' @param nsim number of replications
#' @param p    parameter list
#' @return data frame of true values, mean estimates, bias, SE and coverage
recover_parameters <- function(n = 200, nsim = 1000, p = default_params()) {

  terms <- c("fa_z", "burden_z", "age_z", "sex", "edu", "fa_z:burden_z")
  truth <- c(p$b1, p$b2, p$b_age, p$b_sex, p$b_edu, p$b3)

  est <- se <- cov <- matrix(NA_real_, nsim, length(terms),
                             dimnames = list(NULL, terms))

  for (i in seq_len(nsim)) {
    m  <- .fit_sim(simulate_cohort(n, p))
    cf <- summary(m)$coefficients
    ci <- stats::confint(m)
    est[i, ] <- cf[terms, "Estimate"]
    se[i, ]  <- cf[terms, "Std. Error"]
    cov[i, ] <- ci[terms, 1] <= truth & ci[terms, 2] >= truth
  }

  data.frame(term      = terms,
             true      = truth,
             mean_est  = colMeans(est),
             bias      = colMeans(est) - truth,
             emp_sd    = apply(est, 2, stats::sd),
             mean_se   = colMeans(se),
             coverage  = colMeans(cov),
             row.names = NULL)
}


#' Power at a single sample size
#'
#' @param n     sample size
#' @param nsim  replications
#' @param alpha significance threshold
#' @param p     parameter list
#' @return named vector: n, power for the interaction, power for the FA main effect
power_at_n <- function(n, nsim = 500, alpha = 0.05, p = default_params()) {
  hit_int <- hit_fa <- logical(nsim)
  for (i in seq_len(nsim)) {
    cf <- summary(.fit_sim(simulate_cohort(n, p)))$coefficients
    hit_int[i] <- cf["fa_z:burden_z", "Pr(>|t|)"] < alpha
    hit_fa[i]  <- cf["fa_z", "Pr(>|t|)"] < alpha
  }
  c(n = n, power_interaction = mean(hit_int), power_fa = mean(hit_fa))
}


#' Power curve across a grid of sample sizes
#'
#' @param ns   vector of sample sizes
#' @param nsim replications per point
#' @param p    parameter list
#' @return data frame of power by n
power_curve <- function(ns = seq(50, 600, by = 50), nsim = 500,
                        p = default_params()) {
  as.data.frame(t(sapply(ns, power_at_n, nsim = nsim, p = p)))
}


#' Required sample size across plausible interaction effect sizes
#'
#' The true size of b3 is unknown before the study is run, and the required N
#' depends on it heavily. Reporting this range is more honest than quoting a
#' single number from one assumed effect size.
#'
#' @param b3_grid candidate interaction coefficients
#' @param ns      sample sizes to search
#' @param nsim    replications per point
#' @param target  power to reach
#' @return data frame of required N per assumed b3
required_n <- function(b3_grid = c(-0.30, -0.20, -0.15, -0.10, -0.05),
                       ns      = c(seq(50, 500, by = 50), seq(600, 2000, by = 200)),
                       nsim    = 200,
                       target  = 0.80) {
  out <- lapply(b3_grid, function(b3) {
    p <- default_params(); p$b3 <- b3
    pw  <- sapply(ns, function(n) power_at_n(n, nsim = nsim, p = p)[["power_interaction"]])
    hit <- ns[which(pw >= target)[1]]
    data.frame(b3 = b3, n_required = ifelse(is.na(hit), NA_integer_, hit))
  })
  do.call(rbind, out)
}


#' Effect of FA measurement error on the interaction
#'
#' Noisy FA attenuates b3 toward zero and costs power. Relevant to any
#' discussion of limitations, since FA reliability varies with acquisition
#' quality, motion, and the choice of extraction method.
#'
#' @param reliability vector of reliability values between 0 and 1
#' @param n           sample size
#' @param nsim        replications
#' @return data frame of mean b3 and power by reliability
attenuation <- function(reliability = c(1.0, 0.9, 0.8, 0.7),
                        n = 200, nsim = 300) {
  out <- lapply(reliability, function(rel) {
    est <- numeric(nsim); hit <- logical(nsim)
    for (i in seq_len(nsim)) {
      d <- simulate_cohort(n)
      d$fa_z <- as.numeric(scale(d$fa_z * sqrt(rel) +
                                 stats::rnorm(nrow(d), 0, sqrt(1 - rel))))
      cf <- summary(.fit_sim(d))$coefficients
      est[i] <- cf["fa_z:burden_z", "Estimate"]
      hit[i] <- cf["fa_z:burden_z", "Pr(>|t|)"] < 0.05
    }
    data.frame(reliability = rel, mean_b3 = mean(est), power = mean(hit))
  })
  do.call(rbind, out)
}


#' Power curve plot
#'
#' @param pwr  output of power_curve()
#' @param file optional PNG path
plot_power <- function(pwr, file = NULL) {
  if (!is.null(file)) grDevices::png(file, width = 1800, height = 1400, res = 240)
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  plot(pwr$n, pwr$power_interaction, type = "b", pch = 19, lwd = 2,
       col = "#c0392b", ylim = c(0, 1),
       xlab = "Sample size (N)", ylab = "Power",
       main = "Power to detect the FA x burden interaction")
  graphics::lines(pwr$n, pwr$power_fa, type = "b", pch = 17, lwd = 2, col = "#2c7fb8")
  graphics::abline(h = 0.80, lty = 2, col = "grey40")
  graphics::legend("bottomright",
                   legend = c("Interaction (b3)", "FA main effect (b1)", "80% power"),
                   col = c("#c0392b", "#2c7fb8", "grey40"),
                   lty = c(1, 1, 2), pch = c(19, 17, NA), lwd = 2,
                   bty = "n", cex = 0.9)
  if (!is.null(file)) grDevices::dev.off()
  invisible(NULL)
}
