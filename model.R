# =============================================================================
# model.R
#
# Core estimation code for the moderated model of white-matter integrity,
# systemic burden, and cognition:
#
#   Cognition = b0 + b1*FA + b2*Burden + b3*(FA x Burden) + covariates + e
#
# Works on simulated data or on a real cohort. Base R only.
# =============================================================================


#' Centre and scale the model variables
#'
#' Standardises FA and burden before the interaction term is formed. This is
#' not cosmetic. On raw scales, b1 and b2 describe the effect of one variable
#' when the other equals zero, and zero is not an attainable value for either
#' FA or a 0-1 burden score. Centring makes them the effect at the sample mean,
#' and greatly reduces collinearity between the product term and its parents.
#'
#' @param data       data frame holding the cohort
#' @param fa         name of the FA column
#' @param burden     name of the burden column
#' @param cognition  name of the cognitive outcome column
#' @param covariates character vector of covariate column names
#' @return data frame with columns fa_z, burden_z, cognition and covariates
prepare_data <- function(data,
                         fa         = "fa",
                         burden     = "burden",
                         cognition  = "cognition",
                         covariates = c("age", "sex", "edu")) {

  needed <- c(fa, burden, cognition, covariates)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols))
    stop("Columns not found in data: ", paste(missing_cols, collapse = ", "))

  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  dropped <- nrow(data) - nrow(d)
  if (dropped > 0)
    message(dropped, " row(s) dropped for missing values.")

  out <- data.frame(
    cognition = as.numeric(d[[cognition]]),
    fa_z      = as.numeric(scale(d[[fa]])),
    burden_z  = as.numeric(scale(d[[burden]]))
  )

  # Continuous covariates are z-scored; binary ones are left as is so their
  # coefficients stay interpretable as a group difference.
  for (cv in covariates) {
    v <- d[[cv]]
    if (is.numeric(v) && length(unique(stats::na.omit(v))) > 2) {
      out[[cv]] <- as.numeric(scale(v))
    } else {
      out[[cv]] <- v
    }
  }

  attr(out, "covariates") <- covariates
  attr(out, "fa_mean")    <- mean(d[[fa]])
  attr(out, "fa_sd")      <- stats::sd(d[[fa]])
  out
}


#' Fit the moderated model
#'
#' @param d          output of prepare_data()
#' @param covariates covariate names; defaults to those recorded by prepare_data
#' @return an lm object
fit_wm_model <- function(d, covariates = attr(d, "covariates")) {
  rhs <- "fa_z * burden_z"
  if (length(covariates)) rhs <- paste(rhs, "+", paste(covariates, collapse = " + "))
  stats::lm(stats::as.formula(paste("cognition ~", rhs)), data = d)
}


#' Tidy coefficient table with confidence intervals
#'
#' @param fit   an lm from fit_wm_model()
#' @param level confidence level
#' @return data frame ready to paste into a manuscript
model_table <- function(fit, level = 0.95) {
  cf <- summary(fit)$coefficients
  ci <- stats::confint(fit, level = level)
  data.frame(
    term      = rownames(cf),
    estimate  = cf[, "Estimate"],
    std_error = cf[, "Std. Error"],
    ci_low    = ci[, 1],
    ci_high   = ci[, 2],
    t_value   = cf[, "t value"],
    p_value   = cf[, "Pr(>|t|)"],
    row.names = NULL
  )
}


#' Simple slopes: the FA-cognition relationship at chosen burden levels
#'
#' This is how a moderation result should be reported. A significant b3 says
#' the slopes differ; simple slopes say what each slope actually is.
#'
#' @param fit    an lm from fit_wm_model()
#' @param levels burden values in SD units, default -1, 0, +1
#' @return data frame of slope, SE, and 95% CI at each burden level
simple_slopes <- function(fit, levels = c(-1, 0, 1)) {
  b  <- stats::coef(fit)
  V  <- stats::vcov(fit)
  ix <- c("fa_z", "fa_z:burden_z")
  if (!all(ix %in% names(b)))
    stop("Model does not contain the expected fa_z and fa_z:burden_z terms.")

  df <- stats::df.residual(fit)
  tc <- stats::qt(0.975, df)

  out <- lapply(levels, function(L) {
    slope <- b["fa_z"] + L * b["fa_z:burden_z"]
    # Var(a + L*b) = Var(a) + L^2 Var(b) + 2L Cov(a,b)
    v  <- V["fa_z", "fa_z"] +
          L^2 * V["fa_z:burden_z", "fa_z:burden_z"] +
          2 * L * V["fa_z", "fa_z:burden_z"]
    se <- sqrt(v)
    data.frame(burden_sd = L,
               slope     = unname(slope),
               std_error = se,
               ci_low    = unname(slope) - tc * se,
               ci_high   = unname(slope) + tc * se,
               p_value   = 2 * stats::pt(-abs(unname(slope) / se), df))
  })
  do.call(rbind, out)
}


#' Regression diagnostics
#'
#' Checks the assumptions the manuscript claims to have checked: linearity,
#' normality of residuals, homoscedasticity, and multicollinearity.
#'
#' @param fit an lm from fit_wm_model()
#' @return invisible list of diagnostic statistics
check_assumptions <- function(fit) {
  r <- stats::resid(fit)
  f <- stats::fitted(fit)

  # Shapiro-Wilk caps out at 5000 observations
  sw <- if (length(r) <= 5000) stats::shapiro.test(r) else NULL

  # Breusch-Pagan test for non-constant error variance, computed directly
  # so that no external package is required.
  u   <- r^2 / mean(r^2)
  aux <- stats::lm(u ~ f)
  bp  <- 0.5 * sum((stats::fitted(aux) - mean(u))^2)
  bp_p <- stats::pchisq(bp, df = 1, lower.tail = FALSE)

  # Variance inflation factors via the R-squared of each predictor on the rest
  X <- stats::model.matrix(fit)[, -1, drop = FALSE]
  vif <- sapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(stats::lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  })
  names(vif) <- colnames(X)

  cat("\n--- Assumption checks ---\n")
  if (!is.null(sw))
    cat(sprintf("Shapiro-Wilk (residual normality): W = %.4f, p = %.4f\n",
                sw$statistic, sw$p.value))
  cat(sprintf("Breusch-Pagan (homoscedasticity) : LM = %.4f, p = %.4f\n", bp, bp_p))
  cat("Variance inflation factors:\n")
  print(round(vif, 3))
  cat("\nVIF above about 5 indicates problematic collinearity. Note that a\n")
  cat("product term will always show some inflation; this is expected and is\n")
  cat("minimised by the centring done in prepare_data().\n")

  invisible(list(shapiro = sw, bp_stat = bp, bp_p = bp_p, vif = vif))
}


#' Interaction plot
#'
#' @param fit    an lm from fit_wm_model()
#' @param levels burden levels in SD units
#' @param file   optional path; if given, writes a PNG
plot_interaction <- function(fit, levels = c(-1, 0, 1), file = NULL) {
  if (!is.null(file)) grDevices::png(file, width = 1800, height = 1400, res = 240)
  graphics::par(mar = c(6, 4.5, 3, 1))

  grid <- seq(-2, 2, length.out = 100)
  cols <- c("#2c7fb8", "#7a7a7a", "#c0392b")
  labs <- paste0("Burden ", ifelse(levels > 0, "+", ""), levels, " SD")

  # Build prediction frames holding covariates at their reference values
  mf   <- stats::model.frame(fit)
  covs <- setdiff(names(mf), c("cognition", "fa_z", "burden_z"))

  preds <- lapply(levels, function(L) {
    nd <- data.frame(fa_z = grid, burden_z = L)
    for (cv in covs) {
      v <- mf[[cv]]
      nd[[cv]] <- if (is.numeric(v)) mean(v) else factor(levels(v)[1], levels(v))
    }
    stats::predict(fit, nd)
  })

  yr <- range(unlist(preds))
  yr <- yr + c(-0.15, 0.15) * diff(yr)

  plot(NA, xlim = c(-2, 2), ylim = yr,
       xlab = "Fractional anisotropy (SD units)",
       ylab = "Cognitive performance",
       main = "Predicted cognition by FA across burden levels")
  graphics::abline(h = 0, col = "grey85")
  for (j in seq_along(levels))
    graphics::lines(grid, preds[[j]], col = cols[j], lwd = 3)
  graphics::legend("topleft", legend = labs, col = cols, lwd = 3,
                   bty = "n", cex = 0.9)
  graphics::mtext("Diverging slopes indicate moderation of the FA-cognition relationship.",
                  side = 1, line = 4.6, cex = 0.72, col = "grey30")

  if (!is.null(file)) grDevices::dev.off()
  invisible(NULL)
}
