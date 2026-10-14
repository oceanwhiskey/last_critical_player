# Paper benchmark

This standalone repository contains the algorithms and exactly the measurement
cases shown in the four benchmark figures of *The Last Critical Player:
Computing Deegan–Packel and Johnston in Weighted Majority Games*.

Python 3.10 or later is required. From the repository directory, run:

```console
python -m pip install -r requirements.txt
python -B benchmark.py
```

The script runs new measurements and writes four CSV tables, four Markdown
tables, and the exact index vectors to `tables/<timestamp>/`. The tables are
named `uniform_n`, `fixed_q_n`, `staudacher_n`, and `staudacher_q`, in the same
order as the figures in the paper. Successful measurements report the median
runtime over the number of repetitions used in the paper. Each worker has a
120-second time limit covering all its repetitions. Timeouts are explicitly
reported.

Algorithm timing includes the complete calculation, including sorting,
diagram construction, and exact rational results. It excludes imports,
process startup, and file output. Only one worker runs at a time. A measurement
point appearing in two figures is measured once and used in both tables.
New runtimes depend on the machine and will differ from the published timings.

Only the methods plotted in the paper are measured: Staudacher DP,
Kakimura–Terai, QOBDD with memoized construction, the paper's improved weighted
residual evaluations of Bolus's formulas, and the aggregated Armijos–Toro
Johnston DP. Points omitted from the figures after earlier timeouts are not
added to the run. Points shown with timeout markers are attempted again.

## Algorithms and inputs

`power_indices.py` contains the aggregated Johnston DP, the winning-QOBDD
constructor, scalar last-critical rewards, weighted Shapley–Shubik rewards,
and the shared backward evaluations for Deegan–Packel and Johnston, together
with their helper functions.

`weighted_residual_indices.py` implements the improved `WeightedResidualShares`
evaluations discussed in the paper for one-dimensional weighted games:
Shapley–Shubik on the winning QOBDD and Deegan–Packel on the minimal-winning
QOBDD. It includes the interval-based winning-QOBDD constructor and the MinWin
transformation required by those evaluations. This is the `Bolus_improved`
curve in the figures, named `bolus_weighted` in the tables. The original
coefficient-convolution evaluations are not included or measured.

`staudacher_dp.py` and `kakimura_reference.py` contain the other comparison
algorithms. `data/cases.json` contains the complete weights, quotas, measurement
jobs, and exact reference vectors for the plotted points. Historical runtimes
are neither read nor presented as new results. Every completed calculation is
checked exactly against available reference vectors and against the other
completed methods for the same case and index.

The Python names follow the terminology and pseudocode in the paper:

| Python | Paper terminology |
| --- | --- |
| `johnston_dp` | `JohnstonDP` |
| `prefix_rows` | `PrefixRows` |
| `minimal_winner_rewards`, `last_critical_rewards`, `suffix_counts` | Preparation of alpha, h, and cnt in `LocalRewards` |
| `shared_reward_sums` | Unnormalized contributions from the backward evaluation in `SharedShares` |
| `deegan_packel`, `johnston` in `QOBDDIndices` | Normalized results of `SharedShares` |
| `banzhaf` | `Banzhaf` |
| `shapley_shubik` in `QOBDDIndices` | Weighted last-critical Shapley–Shubik evaluation |
| `shapley_shubik`, `deegan_packel` in `WeightedResidualIndices` | The two `WeightedResidualShares` evaluations |

Method names in the tables are `staudacher_dp`, `kakimura_terai`, `qobdd`,
`bolus_weighted`, and `armijos_toro_dp`. They correspond to the five figure
legends. The Shapley–Shubik index is consistently named `shapley_shubik`.

## Running a small plotted case

For a quick first run using an actual case from the figures:

```console
python -B benchmark.py --figure uniform_n --case uniform_n12
```

A full run may take several hours because of the large cases and timeouts.
The directory can be copied to another location and initialized as a new Git
repository. It requires no files from the original workspace.
