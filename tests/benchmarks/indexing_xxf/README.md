indexing_xxf benchmark
======================

This benchmarks the indexing routines for `xxf_lo` by repeatedly calling
`ik_idx`, `il_idx`, `ie_idx`, `is_idx`, `ig_idx`, `isign_idx` and `idx`. The
number of calls is controlled by `nstep`, with default 100000.

