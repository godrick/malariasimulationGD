# Advisory lint baseline — 2026-11-05

This is the Stage 3 starting point for shared conventions, measured before any
source-formatting cleanup. The source revision is
`115509c57b3e7e3ae5e737c089a91710f90444dc`; the development branch is
`codex/shared-conventions`.

This report records static findings, not scientific validity, runtime test
failures, or a requirement to make the whole repository lint-free immediately.
No existing R/C++ implementation, tests, vignettes, generated interfaces, or
developer scenarios were reformatted or corrected in this stage.

## Policy and scope

- `.editorconfig` provides shared text-save settings and two-space R/YAML/JSON
  indentation. Markdown hard line breaks are preserved. Generated interfaces,
  vendored headers, built documentation, and `test-scripts/` opt out. Existing
  C++ indentation remains a surrounding-file convention.
- `.lintr` names eleven rules explicitly. Naming, function length, return style,
  pipe conversion, and dynamic object usage are not checked by this initial
  policy. Future default-linter additions do not silently enter this list.
- `error_on_lint: FALSE` makes findings advisory. No GitHub workflow, commit
  hook, automatic formatter, or enforced merge gate was introduced.
- The measured command was `lintr::lint_package(".", cache = FALSE)` on a
  temporary copy of tracked package files plus the new configuration. It
  checks the package R-code locations, including tests and vignette code.
- Generated `R/RcppExports.R`, `inst/include`, `inst/doc`, `test-scripts`, `renv`,
  and `packrat` are explicitly excluded. C++ is outside lintr's language scope;
  ordinary `.md` code fences are not scanned by this command.

Developer-only tools remain outside package dependencies. Tool installation
instructions use the ignored developer library described in
[CONTRIBUTING.md](../CONTRIBUTING.md#shared-conventions-and-advisory-linting).
The configuration, this report, and the contributor/agent workflow guides are
excluded from package builds.

## Measured results

Environment: R 4.6.1, lintr 3.4.0, styler 1.11.0, macOS ARM64.
The scan completed in 6.399 seconds with **2,756 findings in 102 files**:
2,750 style findings, six warnings, and no error-type findings.

| Rule | Findings |
| --- | ---: |
| `trailing_whitespace_linter` | 785 |
| `infix_spaces_linter` | 592 |
| `line_length_linter` (100 characters) | 501 |
| `indentation_linter` (two spaces) | 321 |
| `commas_linter` | 309 |
| `brace_linter` | 185 |
| `assignment_linter` | 29 |
| `T_and_F_symbol_linter` | 20 |
| `trailing_blank_lines_linter` | 8 |
| `seq_linter` | 6 |
| `equals_na_linter` | 0 |
| **Total** | **2,756** |

| Area | Findings | Files with findings |
| --- | ---: | ---: |
| `R/` | 1,024 | 42 |
| `tests/` | 1,146 | 45 |
| `vignettes/` | 586 | 15 |

These counts describe this revision and tool version, not a permanent numeric
threshold. Changing the tool version, policy, or source can change the counts.
The elapsed time is a local workflow observation, not a performance benchmark.
Raw findings and validation logs stayed in a temporary directory.

## Warnings requiring contextual review

All six warnings concern sequence expressions and possible empty-input behavior:

| Location at the recorded revision | Rule |
| --- | --- |
| `R/compartmental.R:25` | `seq_linter` |
| `R/compatibility.R:225` | `seq_linter` |
| `R/variables.R:676` | `seq_linter` |
| `R/vector_control_parameters.R:175` | `seq_linter` |
| `tests/testthat/test-vivax-eq.R:33` | `seq_linter` |
| `vignettes/Antimalarial_Resistance.Rmd:339` | `seq_linter` |

These are not six confirmed bugs. Establish whether inputs can be empty and
what the documented behavior should be before replacing sequence constructs.
Any behavioral repair needs its own focused regression and review, including
possible effects on loop execution, random draws, or ordering.

## Configuration validation

- Temporary fixtures exercised assignment, NA-comparison, indentation,
  sequence, boolean-abbreviation, and line-length checks. A 100-character
  line was accepted; a 101-character line was reported.
- Deliberately malformed R files in generated/developer/vendor exclusions
  were skipped by directory linting, while a maintained fixture was reported.
  The configuration probes produced no warnings and returned normally despite
  advisory findings.
- An R Markdown fixture preserved a prose hard break and ignored a long prose
  line, while reporting trailing whitespace in its R chunk.
- EditorConfig structure and representative path matches were checked against
  its specification, including R/R Markdown indentation, Markdown whitespace,
  YAML settings, and generated/vendor/developer opt-outs. No editor plugin was
  installed or exercised as part of this check.
- `styler` dry mode left a maintained pilot file byte-for-byte unchanged. A
  formatted temporary copy of `R/output_parameters.R` had the same parsed R
  expression as its original. This is a limited preview check, not a guarantee
  that arbitrary formatting changes preserve every model behavior.

The source archive was built without a manual or rebuilt vignettes. Archive
contents were checked to exclude the shared configuration, workflow guides,
this report, and `test-scripts/`. The build completed with local archive
UID/GID normalization warnings.

Existing source-file hashes are checked during completion review. Numerical
simulation tests are not rerun for configuration/documentation changes that do
not alter implementation or tests. No dependencies were installed or upgraded.

## Incremental follow-up

Review findings in the area being changed. Keep a later formatting-only pilot
small, inspect its diff, and run its relevant tests before merging it. Preserve
scientific/public names and use narrow, explained exceptions where needed.
Do not apply a whole-tree formatter or treat mechanical replacements as
scientific validation.

GitHub automation and any mandatory lint checks belong to the next workflow
stage, after this policy is reviewed. The R/compiler support decision recorded
in [the repair report](BASELINE-REPAIRS.md) remains separate and unchanged.
