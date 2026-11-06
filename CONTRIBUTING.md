# Development workflow

Develop on a separate branch and keep `main` stable. The maintainer requests a
merge only after reviewing the changes and completing testing.

## Set up an isolated development library

Work from the repository root. Follow [the checkout installation steps](README.md#install-from-a-checkout)
to install the current source into `test-scripts/_lib` and set `.libPaths()` and
`R_LIBS_USER` for the current R session. That directory is already ignored by
Git and excluded from package builds. Do not commit its contents or install
through `test-scripts/install_local.R` without considering its use of the default
R library.

For development, install the additional tools into that same library:

```r
utils::install.packages(
  c("pkgload", "pkgbuild", "testthat", "roxygen2", "rcmdcheck"),
  lib = dev_lib, repos = repos
)
```

For a complete package check, install the declared suggested packages as well:

```r
remotes::install_deps(
  ".", lib = dev_lib, dependencies = TRUE, upgrade = "never", repos = repos
)
```

Run these setup commands intentionally: they can download dependencies. For
focused work, install only the optional dependencies your tests need and report
any skips. Developer scenarios can have additional requirements in their own
READMEs; they do not become package dependencies merely because a developer
workflow uses them.

For terminal commands in the same checkout, select the library for child R
sessions explicitly. In a POSIX shell:

```sh
export R_LIBS_USER="$PWD/test-scripts/_lib"
```

In PowerShell, use `$env:R_LIBS_USER = Join-Path $PWD "test-scripts/_lib"`.
An existing installed dependency may still be read from another library;
installation targets the explicit developer library. `upgrade = "never"` is
not a dependency lock: the unpinned `Remotes` entries can resolve to different
upstream revisions in a fresh environment. Record versions and resolved source
revisions for scientific comparisons.

## Start a task

1. Inspect the current branch and `git status`. Preserve existing uncommitted
   work before starting anything new.
2. Create a focused branch from `main`, normally named `codex/<topic>` for
   assistant-led work. For example, from a clean `main` checkout:

   ```sh
   git switch -c codex/fix-release-validation
   ```

3. Keep unrelated work on its own branch. Do not reset, discard, or silently
   stash changes to make a branch switch possible. If related work was started
   on `main` but remains uncommitted, create the task branch while preserving
   that work before continuing.

If a task depends on an unmerged repair branch, create a separate branch from
that prerequisite and state the dependency in the review. Compare the new work
against that base; do not merge the prerequisite merely to start another task.

The assistant has standing permission to create development branches. Commits,
pushes, merges, and history changes still require an explicit request. Finishing
implementation or passing tests is not permission to merge. Routine cleanup of
fully merged local task branches is authorized after a requested merge, unless
the maintainer asks to retain them.

## Implement and test

Use [the code guide](DEVELOPMENT.md) to find the relevant entry point, backend,
state, and test file before changing behavior. Use two-space indentation,
`snake_case` for new internal names, `<-` assignment, explicit `TRUE`/`FALSE`, and
namespaced calls where ambiguity is possible. Follow the surrounding file and
avoid unrelated reformatting or renaming public functions, arguments, or outputs.

### Shared conventions and advisory linting

`.editorconfig` shares UTF-8, LF line endings, final newlines, and two-space
indentation for R/R Markdown and YAML/JSON files with editors that support it.
Markdown trailing spaces are preserved for deliberate hard line breaks. C++
indentation is not standardized by this stage; follow the surrounding file.
Generated interfaces, vendored headers, built documentation, and everything in
`test-scripts/` opt out of these editor rules. Enabling editor support does not
itself reformat existing files; review any save-time changes before committing.

`.lintr` defines a small, explicit advisory rule set rather than inheriting
every current or future default. It checks assignment, braces, commas, operator
spacing, two-space indentation, a 100-character R line limit, trailing
whitespace/newlines, `T`/`F` abbreviations, NA comparisons, and sequence patterns.
The initial rules do not enforce object names, function length, pipe conversion,
return style, or dynamic object usage. Preserve scientific symbols and public
names. Lint findings are review prompts, not proof of a model bug or permission
to change behavior.

Install the developer tools into the selected library when needed:

```r
utils::install.packages(c("lintr", "styler"), lib = dev_lib, repos = repos)
```

From the repository root, report package findings without modifying source:

```r
lints <- lintr::lint_package(".", cache = FALSE)
print(lints)
length(lints)
```

`lint_package()` covers the package's R-code locations, including tests and
vignettes; it does not lint C++ or R snippets inside ordinary `.md` files. The
configuration explicitly excludes `R/RcppExports.R`, vendored/built files, and
`test-scripts/`. Use this package command rather than recursively linting the
whole developer workspace. Direct `lintr::lint()` calls target a file explicitly;
keep the same scope when selecting one.

Advisory mode is recorded with `error_on_lint: FALSE`. Existing findings do not
block a commit or merge in this stage, and no CI or commit hook is added. Review
new findings in the area being edited, keep any `# nolint: rule_name.` exception
narrow and explained, and do not suppress whole maintained modules to reduce
the count. See [the initial lint baseline](tests/LINT-BASELINE.md) for measured
findings and tool versions. Future rule changes and enforced checks need their
own reviewed update.

### Preview formatting before applying it

Use `styler` on one maintained file at a time and keep formatting in a separate
change from scientific fixes. The [styler preview mode](https://styler.r-lib.org/reference/style_file.html)
reports whether a file would change without writing it:

```r
styler::cache_deactivate(verbose = FALSE)
styler::style_file(
  "R/output_parameters.R", dry = "on", include_roxygen_examples = FALSE
)
```

To inspect formatted text, write a candidate in a temporary directory:

```r
candidate_dir <- tempfile("malariasimulationGD-style-")
dir.create(candidate_dir)
candidate <- file.path(candidate_dir, "output_parameters.R")
stopifnot(file.copy("R/output_parameters.R", candidate))
styler::style_file(candidate, include_roxygen_examples = FALSE)
cat(paste(readLines(candidate), collapse = "\n"))
```

Compare the candidate with the original, review the diff, and run the affected
tests before accepting a later formatting change. The two tools are
complementary: `styler` does not resolve every lint, such as long lines, and its
exclusions do not automatically come from `.lintr`. Do not run a whole-tree
formatter, format generated Rcpp interfaces or vendored code, or sweep the
developer scenarios into a package-formatting change. No existing source was
formatted when introducing these configuration files.

Keep each change reviewable. Separate formatting and refactoring from changes
to equations, defaults, timing, stochastic sampling, or model behavior. Follow
the scientific and reproducibility requirements in [AGENTS.md](AGENTS.md).

Run the smallest relevant test first, then broader checks appropriate to the
risk. For example, replace the filter below with the affected test-file stem:

```sh
Rscript -e 'testthat::test_local(".", filter = "native-mosquito-backend")'
Rscript -e 'testthat::test_local(".")'
```

For interactive source loading, use `pkgload::load_all(".")`. It is a development
convenience, not evidence that installation works. `testthat::test_local()`
loads source and can compile into `src/`; package checking also tests the
installed context. Restart R after changes to a loaded native library.

| Change | Minimum relevant validation |
| --- | --- |
| README or contributor instructions | Check links and runnable snippets; verify archive exclusions when packaging rules change |
| R behavior or input validation | Focused successful/invalid-input regression, then broader affected tests |
| Backend or C++ behavior | Affected deterministic/stochastic and single-node/metapopulation cases; source and installed contexts |
| Releases or mobility | Consumed schedule/counts or identity/nontrivial/returning-trip invariants; relevant resume equivalence |
| Checkpoint or RNG changes | Continuous versus resumed comparison, explicit seeds, and documented changes to sampling |
| Performance changes | Correctness baseline first, then comparable timed scenarios with replications and core counts |

Set seeds explicitly. Reset to the same seed before each comparison arm when
configuration should be the cause of any difference. Use scientific invariants
and appropriate tolerances; do not weaken assertions to make a change pass.

Build and check a disposable source copy when compiled artifacts or documentation
outputs would otherwise affect the working tree. Use an isolated R library
when an installation is needed. Check user-facing vignette builds separately
from routine checks that disable vignette execution.

For a bounded build/check from R, with artifacts in a temporary location:

```r
check_dir <- tempfile("malariasimulationGD-check-")
dir.create(check_dir)
check_result <- rcmdcheck::rcmdcheck(
  path = ".",
  build_args = c("--no-manual", "--no-build-vignettes"),
  args = c("--no-manual", "--no-build-vignettes", "--no-vignettes"),
  check_dir = check_dir,
  error_on = "error"
)
```

Run this in a disposable copy when source build/compile side effects would be
unwelcome. [`rcmdcheck` accepts build/check arguments separately](https://rcmdcheck.r-lib.org/reference/rcmdcheck.html).
The command runs package tests and examples, but skips vignette execution and
building. A tarball without built vignettes can therefore report missing
`inst/doc` warnings; this shortcut is not a clean release check. Setting
`error_on = "error"` keeps those warnings in the returned report rather than
stopping on them. Record and inspect warnings; do not treat the absence of an
R exception as a clean package check.

For the documentation check, use `rmarkdown::pandoc_available()` to check that
rmarkdown can find Pandoc, then build the vignettes and check the resulting
tarball. From a disposable source copy:

```r
artifact_dir <- tempfile("malariasimulationGD-build-")
dir.create(artifact_dir)
tarball <- pkgbuild::build(
  ".", dest_path = artifact_dir, vignettes = TRUE, manual = FALSE
)
release_check <- rcmdcheck::rcmdcheck(
  tarball,
  args = c("--no-manual", "--no-build-vignettes", "--no-vignettes"),
  check_dir = file.path(artifact_dir, "check"),
  error_on = "warning"
)
```

This evaluates the vignette chunks allowed by their settings during the build
and checks the shipped artifacts without repeating those simulations. It does
not produce the PDF manual. The distinction between the check flags is described
in [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Checking-and-building-packages).

Record failures, warnings, notes, and skipped tests separately. `CI=true`
activates the existing `skip_on_ci()` calls; a CI run and a full local run have
different coverage. The historical R minimum and C++14 request still need a
supported-version decision; see [the recorded compatibility limits](tests/BASELINE-REPAIRS.md#remaining-compatibility-decisions-and-limits).

`test-scripts/` is the developer workspace for scientific scenarios, exploratory
examples, and outputs; it is excluded from the distributed package. Its optional
workflow tests run through [the developer test runner](test-scripts/tests/README.md).
Package regressions in `tests/testthat/` must run independently of that workspace.

## Update documentation and generated interfaces

Write public function/parameter/return-value changes in roxygen comments in
`R/`. Regenerate documentation in a disposable source copy and review only the
relevant generated diff:

```sh
Rscript -e 'roxygen2::roxygenize(".")'
```

The recorded generator version is `Config/roxygen2/version` in `DESCRIPTION`.
Keep generator-only changes separate from the behavioral change being reviewed.

Do not hand-edit `NAMESPACE` or `man/`. If an exported C++ signature changes,
regenerate and review both Rcpp interface files:

```sh
Rscript -e 'Rcpp::compileAttributes(".")'
```

Update `NEWS.md` for user-visible behavior changes and state expected numerical,
stochastic, runtime, or compatibility effects. Keep internal implementation
helpers unexported unless they intentionally become part of the public API.

## Prepare small commits

When the maintainer requests commits, group the authorized work into small,
coherent commits on the task branch. Each commit should have one understandable
purpose; keep a behavior change with its regression test and directly related
documentation. Separate independent configuration, reporting, formatting, and
scientific changes where practical, without splitting a working change into
broken intermediate states.

Write a short, plain-language title describing the result. Add a body explaining
what changed, why it was needed, and the relevant checks or limitations. Avoid
vague titles such as "updates" or "fixes". Review `git diff --cached` before each
commit, stage only the intended files or hunks, and preserve unrelated user work
and generated runtime outputs. Report the commit hashes and titles afterward.

Committing requires an explicit request. A request to commit does not authorize
pushing, merging, or rewriting existing history; those remain separate actions.

## Review and merge

Report the branch name, changed behavior, tests run, scientific or numerical
effects, and checks not performed. The maintainer then reviews and tests the
branch. Merge into `main` only when the maintainer explicitly requests it after
testing. Recheck affected behavior if changes are made after validation.

Use [the pull-request template](.github/pull_request_template.md) to describe the
triggering problem, resulting behavior, affected backends, RNG/checkpoint
effects, and validation. For documentation-only changes, say which examples or
commands were checked and which setup/download steps were not executed.

For GitHub, the recommended repository policy is pull requests into `main`
with required validation checks. Do not push development changes directly to
`main`. Documenting this policy does not configure remote repository settings.
Shared lint/editor conventions are advisory. GitHub automation and required CI
checks are a later workflow step; the template itself does not run or enforce
them.

## Clean up merged branches

Clean up local task branches after the requested merge has been verified, unless
the maintainer asks to retain them. List local candidates with:

```sh
git branch --merged main
```

Before deleting a candidate, confirm its tip is contained in `main`, it has no
unfinished work, and it is not the current branch or used by another worktree.
Preserve `main` and any branch with unmerged work. Use `git branch -d` for a
verified local candidate; do not force-delete branches automatically.

Squash merges and cherry-picks may not appear in the merged-branch list. Review
those cases separately. Deleting a remote branch also requires an explicit
request.
