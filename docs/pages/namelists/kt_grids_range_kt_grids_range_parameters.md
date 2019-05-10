title: kt_grids_range_parameters

#kt_grids_range_parameters

| Name | Type | Default | Description |
|:-----|:----:|:-------:|:------------|
| akx_max | real | 0.0 | Max kx for periodic finite kx ballooning space runs with|
| | | |\(\hat{s}=0\). |
| akx_min | real | 0.0 | Min kx for periodic finite kx ballooning space runs with|
| | | |\(\hat{s}=0\). |
| aky_max | real | 0.0 | Upper limit of \(k_y \rho\) range. Should set to something other|
| | | |than zero. |
| aky_min | real | 0.0 | Lower limit of \(k_y \rho\) range. Should typically set to|
| | | |something other than zero. |
| kyspacing_option | character(len = 20) | 'default' | Sets the type of spacing between ky grid points, available options are :|
| | | |-  'default' : Same as 'linear'|
| | | |-  'exponential' : Evenly spaced in log(ky).|
| | | |-  'linear' : Evenly spaced in ky. |
| n0_max | integer | 0 | Maximum toroidal mode number. Use instead of `aky_max`. |
| n0_min | integer | 0 | Minimum toroidal mode number. Use instead of `aky_min`. |
| naky | integer | 1 | The number of 'actual' ky modes. |
| nn0 | integer | 1 | Number of toroidal modes, only used if `n0_min`>0. Overrides `naky`|
| | | |in kt_grids_range_parameters. Note if the number of modes isn't compatible|
| | | |with the requested min and max toroidal mode numbers then we just run with|
| | | |one mode, determined by `n0_max`. |
| ntheta0 | integer | 1   | Number of \(\theta_0\) (kx) modes |
| rhostar_range | real | 1.0e-4 | Used to convert `n0_min`, `n0_max` range into `aky_min`,|
| | | |`aky_max`, if `n0_min`>0. If `n0_min` is set,|
| | | |`aky_min=n0_min*drhodpsi*rhostar_range` and|
| | | |`aky_max=n0_max*drhodpsi*rhostar_range` |
| theta0_max | real | 0.0 | Upper limit of `theta_0` range |
| theta0_min | real | 0.0 | Lower limit of `theta_0` range |