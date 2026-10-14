# Published paper measurements

This directory contains the historical measurements used in the four benchmark
figures of *The Last Critical Player: Computing Deegan–Packel and Johnston
in Weighted Majority Games*. These are fixed paper results, extracted from
the recorded measurements rather than from a new execution of the algorithms.

| Figure / table | Successful entries | Timeout entries |
| --- | ---: | ---: |
| `uniform_n` | 130 | 0 |
| `fixed_q_n` | 91 | 0 |
| `staudacher_n` | 96 | 10 |
| `staudacher_q` | 34 | 5 |

Each figure has a CSV table and an equivalent Markdown table. Entries identify
the case, power index, method, player count, quota, repetitions, status,
runtime in seconds, QOBDD node count where available, and source report.
Only the points and timeout markers actually shown in the figures are included.
Entries omitted from a plot after an earlier timeout, unsupported methods,
and methods not plotted in the paper are excluded. A point shared by two
figures appears in both tables; the counts above are table entries, not unique
algorithm calls.

`results.json` contains the same entries with the exact rational output vectors,
individual timing samples and other recorded measurement fields when available.
The `seconds` column preserves the historical value used to generate the plot;
`median_seconds` in the JSON is the corresponding numeric measurement. Timeout
entries have no completed vector or measured runtime. Their plot marker is
placed at the 120-second limit, which is not a completed runtime.

## Measurement conditions

The original measurement reports record CPython 3.10.4 on
`Windows-10-10.0.19045-SP0`. Successful base-suite measurements use the median
of three complete calls. Additional extended cases and Staudacher linear-range
fixtures use a single call. The repetition count is included in each table row.
Workers ran sequentially with a 120-second limit per worker, including startup
and all repetitions. Algorithm timings exclude imports and process startup,
and include the complete algorithm call, diagram construction where applicable,
and exact rational results.

All paper measurements were performed on an Intel Core i5-3320M processor
(nominal clock 2.60 GHz) with 8 GB RAM. CPU and RAM were not recorded in the
original reports. On 2026-10-14, the author confirmed that all measurements
used the same machine and hardware configuration; the CPU was also read from
the machine's registry. This retrospective hardware confirmation is recorded
separately in `provenance.json`, without changing the original source settings.
The scripts in this repository are a reduced reviewer edition; the historical
timings describe the recorded implementations and environments, not a timing
run of this reduced edition.

## Provenance and verification

The figure data came from the combined report identified as
`20261011_pair_free_default/combined/power_index_benchmark.json`.
That report combines measurements from multiple runs. QOBDD measurements
include the memoized construction and scalar last-critical preparation;
Shapley–Shubik was subsequently refreshed using weighted last-critical rewards.
The improved Bolus readings use the corrected unpacked weighted evaluations.
The aggregated Armijos–Toro Johnston DP was measured separately on the same
cases. Previously measured comparison values were carried forward.

Each JSON entry has a `provenance_chain` from the combined report to its original
measurement report and a `source_report` identifying that final report.
`provenance.json` records the source report settings and SHA-256 digests,
plus the digests of the CSV and LaTeX plot sources used for extraction.
Source report names are historical identifiers indexed in `provenance.json`,
not dependencies on files outside this repository. Absolute workstation paths
are not included. Original reports are identified by their digests; this
directory contains the selected measurements and provenance metadata, not full
copies of unrelated historical benchmark archives.

Method identifiers match the reviewer scripts: `staudacher_dp`,
`kakimura_terai`, `qobdd`, `bolus_weighted`, and `armijos_toro_dp`.
The Shapley–Shubik index is named `shapley_shubik`. Original method and index
identifiers remain in `recorded_method` and `recorded_index`; naming changes
do not change recorded measurements or algorithm descriptions.

The extraction was checked against every plotted coordinate and timeout marker
in the four figure sources. All included input games and completed exact vectors
were also checked against `../../data/cases.json`. Historical timing samples,
statuses and exact vectors were checked along the provenance chains.
`sha256.json` provides checksums for the files in this directory.

New benchmark runs write to `tables/<timestamp>/` and do not overwrite these
published measurements. Small development checks are not part of this dataset.
