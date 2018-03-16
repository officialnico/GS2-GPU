redists g to gf
---------------

This times the redistribution from g_lo to gf_lo which is used in the gf_local field solve.  It repeatedly calls gether/scatter for each of the layouts lexys, xyles and yxles.  The number of iterations is controlled by nstep, default 1000.
