# wm-burden-model

A small R toolkit for a moderated model linking white-matter microstructure,
systemic vascular burden, and cognition:

```
Cognition = b0 + b1·FA + b2·Burden + b3·(FA × Burden) + covariates + e
```

The term of interest is **b3**. A negative b3 means that the relationship
between fractional anisotropy and cognition weakens as vascular burden rises —
the "double-hit" pattern, in which burden does not merely subtract from
cognition but changes how consequential white-matter damage is.

The repository provides the estimation code, a simulation study establishing
that the model is identifiable, and a power analysis giving the sample size
required to detect the interaction.

---

## ⚠️ On the data

**No participant data are included or analysed here.** Every number produced by
`examples/run_simulation.R` comes from synthetic data generated under assumed
parameter values. These validate the *code and design*; they are not empirical
findings about any cohort.

`examples/run_real_data.R` is a template for applying the same model to a real
dataset once one is available.

---

## Requirements

R (≥ 4.0). No packages. Everything uses base R, so there is nothing to install
beyond R itself.

## Usage

```bash
git clone https://github.com/<user>/wm-burden-model.git
cd wm-burden-model
Rscript examples/run_simulation.R
```

Or interactively:

```r
source("wmburden.R")

d   <- prepare_data(simulate_cohort(200))
fit <- fit_wm_model(d)

model_table(fit)       # coefficients with 95% CIs
simple_slopes(fit)     # FA effect at low / mean / high burden
check_assumptions(fit) # normality, homoscedasticity, VIF
plot_interaction(fit)
```

## What the functions do

| Function | Purpose |
|---|---|
| `prepare_data()` | Centres and scales FA and burden **before** forming the product term |
| `fit_wm_model()` | Fits the moderated regression with covariates |
| `model_table()` | Coefficient table with confidence intervals |
| `simple_slopes()` | FA–cognition slope at each burden level, with CIs |
| `check_assumptions()` | Shapiro–Wilk, Breusch–Pagan, variance inflation factors |
| `plot_interaction()` | Predicted cognition by FA across burden levels |
| `simulate_cohort()` | Generates a synthetic cohort under known parameters |
| `recover_parameters()` | Checks the estimator returns the coefficients it was given |
| `power_curve()` | Power to detect b3 across sample sizes |
| `required_n()` | Required N across a range of plausible b3 values |
| `attenuation()` | How FA measurement error shrinks b3 and costs power |

## Notes on the specification

**Centring is not optional.** Without it, b1 and b2 describe effects at FA = 0
and burden = 0 — values neither variable can take — and the product term is
strongly collinear with its components. `prepare_data()` standardises both
before the interaction is built.

**Covariates matter.** Age drives both FA decline and vascular burden, so an
unadjusted b1 partly measures age. Age, sex, and education are included by
default.

**Simple slopes, not just b3.** A significant interaction says the slopes
differ. It does not say what any slope *is*. Report both.

**Moderation, not mediation.** This model asks whether burden changes the
strength of the FA–cognition relationship. The distinct claim that burden
damages FA which then impairs cognition is a mediation hypothesis and requires
a different two-equation structure. The two should not be conflated.

## Sample size

Interaction terms need substantially larger samples than main effects of
comparable size. Under the default assumptions the FA main effect reaches 80%
power at roughly N = 100, while the interaction needs about N = 150 — and that
figure is highly sensitive to the true b3:

| assumed b3 | N for 80% power |
|---|---|
| −0.30 | ~100 |
| −0.20 | ~150 |
| −0.15 | ~250 |
| −0.10 | ~450 |
| −0.05 | ~1800 |

Since b3 is unknown before the study runs, report the range rather than a
single figure.

## Getting FA values from diffusion images

The model takes one FA value per participant. A standard FSL route:

```bash
topup   # susceptibility distortion correction (requires a field map or
        # reverse phase-encode acquisition)
eddy    # eddy-current and motion correction — note that the older
        # eddy_correct is deprecated
bet     # brain extraction
dtifit  # tensor fitting, produces voxelwise FA
tbss    # skeletonisation, gives a mean FA per participant
```

Tract-specific rather than whole-brain FA requires atlas-based extraction and
correction for multiple comparisons across tracts.

## Suitable datasets

Testing this model requires vascular phenotyping, diffusion imaging, and
cognitive measures in the same participants, ideally in an older sample.
Candidates include ADNI, UK Biobank, the Framingham Offspring imaging
substudy, and OASIS-3. All require an access application.

## License

MIT — see `LICENSE`.
