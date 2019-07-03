indexing_yxf benchmark
======================

This benchmarks the indexing routines for `yxf_lo` by repeatedly calling
`it_idx`, `il_idx`, `ie_idx`, `is_idx`, `ig_idx`, `isign_idx` and `idx`. The
number of calls is controlled by `nstep`, with default 100000.

