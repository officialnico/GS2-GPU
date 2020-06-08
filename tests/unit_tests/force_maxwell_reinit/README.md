Test `force_maxwell_reinit`
---------------------------

This tests that the `force_maxwell_reinit` flag is correctly applied by:

1. running gs2 with `initial.in` for 100 timesteps and saving for restart;
2. restarting with `restart.in` and running for 0 timesteps;
3. comparing phi in the netcdf files `initial.out.nc` and `restart.out.nc` to
   ensure that phi is correctly reconstructed.
