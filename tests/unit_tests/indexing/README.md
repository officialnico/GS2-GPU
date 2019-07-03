Test indexing
=============

This tests the indexing routines for `g_lo`, `xxf_lo`, `yxf_lo` and `le_lo`. It
calls `it_idx`, `ik_idx`, `il_idx`, `ie_idx`, `is_idx`, `ig_idx` and `isign_idx`
to get the dimension indices from the global index, and then calls `idx` and 
checks that the global index is correctly recovered.
