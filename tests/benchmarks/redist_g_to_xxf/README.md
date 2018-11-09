redists g to xxf
----------------

This test times the redistribution from `g_lo` to `xxf_lo` which is used in the 
nonlinear term.  

It repeatedly calls gather/scatter between a g layout (lexys, xyles or yxles) 
and xxf. The number of iterations is controlled by nstep, default 10.
