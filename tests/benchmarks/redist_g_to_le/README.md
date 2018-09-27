redists g to le
---------------

This test times the redistribution from g_lo to le_lo which is used in the 
collisions module.  

It repeatedly calls gather/scatter between a g layout (lexys, xyles or yxles) 
and gf. The number of iterations is controlled by nstep, default 1000.
