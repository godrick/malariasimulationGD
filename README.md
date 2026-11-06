# malariasimulationGD

`malariasimulationGD` is an R package for malaria transmission simulation with
individual human infection histories, mosquito genetics and releases,
metapopulation transmission, and explicit overnight human mobility. It supports
*P. falciparum* and *P. vivax*, native count-based mosquito dynamics, and legacy
mosquito backends.

The package is under development. Model equations, defaults, timing, numerical
results, and stochastic reproducibility are reviewed as part of its interface.

## Install from a checkout

Clone this repository and select the branch you want to test. Run the following
in R from the repository root. Source installation needs a compiler toolchain
appropriate for your R installation; see [R's installation manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).

```r
dev_lib <- file.path(normalizePath("."), "test-scripts", "_lib")
dir.create(dev_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(dev_lib, .libPaths()))
Sys.setenv(R_LIBS_USER = dev_lib)

repos <- c(
  CRAN = "https://cloud.r-project.org",
  MRC = "https://mrc-ide.r-universe.dev"
)
if (!requireNamespace("remotes", quietly = TRUE)) {
  utils::install.packages("remotes", lib = dev_lib, repos = repos)
}
remotes::install_local(
  ".", lib = dev_lib, dependencies = NA, upgrade = "never",
  force = TRUE, build_vignettes = FALSE, repos = repos
)
```

This installs the checkout and required dependencies into an ignored developer
library. The repository declares additional dependency sources in `DESCRIPTION`;
[`remotes` understands those declarations](https://remotes.r-lib.org/articles/dependencies.html).
These commands require network access for missing packages. They do not update
the default R library. Re-select `dev_lib` with `.libPaths()` in each new R
session, and restart R after rebuilding a package whose native library is loaded.

The recorded baseline used R 4.6.1 on macOS ARM64. The supported R/compiler
matrix is still being established: the historical R minimum in `DESCRIPTION`
is not a tested support claim, and R 4.6 ignores the current C++14 request.

## Run a small seeded example

After installing and selecting the library above:

```r
library(malariasimulationGD)
set.seed(20261105)

parameters <- get_parameters(list(
  native_mosquito_backend = TRUE,
  individual_mosquitoes = TRUE,
  human_population = 50L,
  total_M = 200,
  init_foim = 0,
  progress_bar = FALSE
))
parameters <- parameterise_total_M(parameters, parameters$total_M)

result <- run_resumable_simulation(10L, parameters = parameters)
head(result$data[c("timestep", "n_infections", "Sm_gamb_count")])
```

This uses the native stochastic count-based tau-leap mosquito backend. It is a
short execution example, not a calibrated transmission scenario. `result$data`
contains daily output; `result$state` holds the checkpoint needed to continue.
The first use can emit a warning explaining the native/legacy backend distinction.
See [the checkpoint example](DEVELOPMENT.md#continue-a-checkpoint) for endpoint
semantics and RNG restoration.

## Choose the mosquito backend deliberately

| `native_mosquito_backend` | `individual_mosquitoes` | Mosquito dynamics |
| --- | --- | --- |
| `TRUE` | `TRUE` | Native stochastic count-based tau-leap |
| `TRUE` | `FALSE` | Native deterministic count-based ODE |
| `FALSE` | `TRUE` | Legacy event-based individual mosquitoes |
| `FALSE` | `FALSE` | Legacy compartmental ODE |

The deterministic labels describe the mosquito component; human infection
processes can still be stochastic. High-level simulation functions select the
native backend when its flag is left `NULL`. Explicit flags make examples and
comparisons easier to interpret. The native tau-leap and legacy individual
backends are distinct models, not interchangeable implementations.

Explicit human mobility requires the native metapopulation backend, identity
transmission-mixing matrices, and zero border-capture matrices. Trips change
where exposure occurs while human state remains on the home record. See
[mobility constraints](DEVELOPMENT.md#movement-and-node-ordering) before combining
these features.

## Find your next step

| Goal | Start here |
| --- | --- |
| Set up development, test, and submit a change | [Contributor workflow](CONTRIBUTING.md) |
| Understand entry points, process order, and checkpoints | [Code and model guide](DEVELOPMENT.md) |
| Read parameter and output definitions | `?get_parameters`, `?run_simulation`, `?run_metapop_simulation` |
| Read the model and intervention tutorials | [User-facing vignettes](vignettes/) |
| Explore scientific scenarios or developer examples | [Developer workspace](test-scripts/README.md) |
| Review behavior changes | [NEWS](NEWS.md) |

Everything under `test-scripts/`, including examples, local libraries, data,
plots, and generated outputs, is developer-only and excluded from the final
package. Small package regressions live independently in `tests/testthat/`.
