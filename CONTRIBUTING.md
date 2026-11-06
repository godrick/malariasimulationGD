# Development workflow

Develop on a separate branch and keep `main` stable. The maintainer requests a
merge only after reviewing the changes and completing testing.

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

The assistant has standing permission to create development branches. Commits,
pushes, merges, and history changes still require an explicit request. Finishing
implementation or passing tests is not permission to merge. Routine cleanup of
fully merged local task branches is authorized after a requested merge, unless
the maintainer asks to retain them.

## Implement and test

Keep each change reviewable. Separate formatting and refactoring from changes
to equations, defaults, timing, stochastic sampling, or model behavior. Follow
the scientific and reproducibility requirements in [AGENTS.md](AGENTS.md).

Run the smallest relevant test first, then broader checks appropriate to the
risk. For example, replace the filter below with the affected test-file stem:

```sh
Rscript -e 'testthat::test_local(".", filter = "native-mosquito-backend")'
Rscript -e 'testthat::test_local(".")'
```

Build and check a disposable source copy when compiled artifacts or documentation
outputs would otherwise affect the working tree. Use an isolated R library
when an installation is needed. Check user-facing vignette builds separately
from routine checks that disable vignette execution.

`test-scripts/` is the developer workspace for scientific scenarios, exploratory
examples, and outputs; it is excluded from the distributed package. Its optional
workflow tests run through [the developer test runner](test-scripts/tests/README.md).
Package regressions in `tests/testthat/` must run independently of that workspace.

## Review and merge

Report the branch name, changed behavior, tests run, scientific or numerical
effects, and checks not performed. The maintainer then reviews and tests the
branch. Merge into `main` only when the maintainer explicitly requests it after
testing. Recheck affected behavior if changes are made after validation.

For GitHub, the recommended repository policy is pull requests into `main`
with required validation checks. Do not push development changes directly to
`main`. Documenting this policy does not configure remote repository settings.

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
