# Runtime Benchmark Report

Date: 2026-06-29

Comparison:
- Baseline: `3a63ce6`
- Current implementation: `67acd39`

Environment:
- R version: 4.5.0 (2025-04-11)
- OS: Darwin 24.6.0 arm64
- Benchmark artifacts: `/private/tmp/msim_benchmark/run_67acd39_vs_3a63ce6`

## Method

Both commits were installed into separate temporary R libraries. Each benchmark ran
five timed replicates after one warmup run. Timings use elapsed wall-clock seconds
from `system.time()`.

Each replicate includes normal public-API setup for the scenario, including
parameter construction and equilibrium setup where used. This matches what a user
pays when running the package normally, but it means the timings are not isolated
microbenchmarks of only the changed functions.

## Direct Before/After Results

Positive percentage means the current implementation was slower than baseline.

| Scenario | Baseline median (s) | Current median (s) | Difference (s) | Change |
|---|---:|---:|---:|---:|
| `single_node_default` | 0.554 | 0.663 | +0.109 | +19.7% |
| `native_cube_rendered` | 0.191 | 0.293 | +0.102 | +53.4% |
| `mobility_rendered` | 1.958 | 2.017 | +0.059 | +3.0% |
| `mobility_state_only` | 1.837 | 1.929 | +0.092 | +5.0% |

Raw replicate times:

| Version | Scenario | Replicate times (s) |
|---|---|---|
| Baseline | `single_node_default` | 0.554, 0.540, 0.549, 0.559, 0.575 |
| Current | `single_node_default` | 0.660, 0.739, 0.658, 0.663, 0.743 |
| Baseline | `native_cube_rendered` | 0.199, 0.191, 0.189, 0.191, 0.191 |
| Current | `native_cube_rendered` | 0.295, 0.293, 0.286, 0.297, 0.286 |
| Baseline | `mobility_rendered` | 1.958, 1.976, 1.955, 1.967, 1.946 |
| Current | `mobility_rendered` | 2.004, 2.001, 2.017, 2.048, 2.047 |
| Baseline | `mobility_state_only` | 1.837, 1.837, 1.843, 1.832, 1.846 |
| Current | `mobility_state_only` | 1.914, 1.915, 1.929, 1.934, 1.968 |

## Additional Current-Only Check

The new single-node `render_output = FALSE` option did not exist in the baseline,
so it cannot be used for a direct before/after comparison. As a current-only check,
it was compared with current rendered output for the same single-node scenario.

| Current mode | Median (s) | Change vs rendered |
|---|---:|---:|
| Rendered output | 0.647 | reference |
| `render_output = FALSE` | 0.608 | -6.0% |

## Conclusion

The implementation is not faster than the previous commit in the tested direct
before/after scenarios. The current code was slower by about 3% to 53%, depending
on the scenario.

The likely reason is that the ring-buffer and native-summary cache add R-level
bookkeeping overhead that is larger than the avoided copying for these tested
problem sizes. In particular, the lag-buffer reads now reconstruct chronological
views, and the native summary cache adds environment lookup/invalidation overhead.

The only measured speedup is the new opt-in single-node state-only mode
(`render_output = FALSE`), which was about 6% faster than current rendered output
in the small single-node scenario. That is useful functionality, but it is not a
direct baseline comparison because the baseline did not expose that single-node
API option.

Recommendation: do not treat commit `67acd39` as a runtime performance win as-is.
Keep or revise only the pieces that are needed for functionality, and redesign the
lag-buffer/native-summary optimizations before claiming speed improvements.
