indexing benchmark
==================

This benchmarks the indexing routines for `g_lo` by repeatedly callings 
`it_idx`, `ik_idx`, `il_idx`, `ie_idx`, `is_idx` and `idx`. The number of calls
is controlled by `nstep`, with default 100000.

