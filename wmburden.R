# =============================================================================
# wmburden.R
#
# Entry point. Load the whole toolkit with:
#
#   source("wmburden.R")
#
# Base R only. No packages required.
# =============================================================================

local({
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  if (is.null(here) || !nzchar(here)) here <- "."
  for (f in c("R/model.R", "R/simulate.R")) {
    path <- file.path(here, f)
    if (!file.exists(path)) path <- f
    source(path, chdir = FALSE)
  }
})

message("wm-burden-model loaded. See examples/run_simulation.R to start.")
