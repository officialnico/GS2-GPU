title: theta_grid_eik_knobs

#theta_grid_eik_knobs

| Name | Type | Default | Description |
|:-----|:----:|:-------:|:------------|
| alpha_input | real | 0.0 | Used in calculation of `dp_new = -alpha_input/qval**2/rmaj*drhodpsi`|
| | | |@note This gets a "smart" default. |
| beta_prime_input | real | 0.0 | The gradient of the pressure. Strictly speaking this parameter|
| | | |is not \(\frac{\partial \beta}{\partial \rho}\) but \(\beta|
| | | |\frac{1}{p}\frac{\partial p}{\partial \rho}\): in other words,|
| | | |the gradient of the magnetic field is ignored. Used only if|
| | | |`bishop` = 4 or 9.|
| | | |@note This gets a "smart" default. |
| bishop | integer | 5 | Use Bishop relations to generate metric coefficients.|
| | | |-  0: Use high-aspect ratio coefficients (only for debugging)|
| | | |-  1: Use actual equilibrium values of shat, p' recommended|
| | | |-  2: Use numerical equilibrium + s_hat_input and p_prime_input|
| | | |-  3: Use numerical equilibrium + s_hat_input and inv_Lp_input|
| | | |-  4: Use numerical equilibrium + s_hat_input and beta_prime_input|
| | | |-  5: Use numerical equilibrium + s_hat_input and alpha_input|
| | | |-  6: Use numerical equilibrium + beta_prime_input|
| | | |-  7: Use numerical equilibrium and multiply pressure gradient by dp_mult|
| | | |-  8: Use numerical equilibrium + s_hat_input and multiply pressure gradient by dp_mult|
| | | |-  9: Use numerical equilibrium + s_hat_input and beta_prime_input|
| | | |-  Otherwise: Use magnetic shear and pressure gradient as set elsewhere. |
| chs_eq | logical | .false. | Use equilbrium data from the CHEASE file ogyropsi.dat |
| delrho | real | 1e-3 | Default usually okay|
| | | |@todo Expand documentation |
| dfit_eq | logical | .false. | Vacuum magnetic dipole geometry |
| dp_mult | real | 1.0 | Used to scale the pressure gradient, only if bishop = 7 or 8. |
| efit_eq | logical | .false. | Use EFIT equilibrium (EFIT, codes with eqdsk format) |
| eqfile | character(len = 800) | "default_unset_value" | Name of file with numerical equilibrium data (if required) |
| equal_arc | logical | .true. | Change field-line coordinate. Recommended value: F|
| | | |@note We recommend `.false.` but default to `.true.`. We should consider|
| | | |changing the default. |
| gen_eq | logical | .false. | Use Toq-style NetCDF equilibrium (TOQ) |
| gs2d_eq | logical | .false. | Read Colin Roach's GS2D equilibrium file |
| idfit_eq | logical | .false. | FIXME: Add documentation |
| iflux | integer | 0 | Choose mode of operation:|
| | | |-  0: Use analytic parameterization of local toroidal MHD equilibrium (Miller).|
| | | |-  1: Read equilibrium data from output of MHD Grad Shafranov solver|
| | | |-  2: Running inside a transport code without numerical equilibrium |
| invlp_input | real | 0. | FIXME: Add documentation |
| irho | integer | 2 | Choose definition of flux surface coordinate|
| | | |-  1: rho == sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)|
| | | |-  2: rho == midplane diameter/LCFS diameter - Recommended|
| | | |-  3: rho == poloidal flux/poloidal flux of LCFS|
| | | |NB For consistency fprim, tprim, shat, uprim, g_exb etc|
| | | |*must be computed using same radial variable*,|
| | | |i.e. *depend on choice of irho*! |
| isym | integer | 0 | 0: Recommended|
| | | |1: Force up-down symmetry.|
| | | |@todo Consider changing type to logical |
| itor | integer | 1 | Do not change|
| | | |@todo Consider removing |
| local_eq | logical | .true. | If `.true.` use Miller-style local equilibrium else use other|
| | | |numerical equilibrium types |
| ppl_eq | logical | .false. | Use Menard-style NetCDF equilibrium (JSOLVER) |
| rmax | real | 1.0 | Never used |
| rmin | real | 1e-3 | Never used |
| s_hat_input | real | 0.0 | Used to overrides s_hat prescribed by the numerical|
| | | |equilibrium, but _only_ if bishop=2,3,4,5,8, or 9.|
| | | |@note This gets a "smart" default. |
| transp_eq | logical | .false. | Use PPL NetCDF equilibrium (psipqgrz equilibrium from TRANSP/TRXPL) |
| writelots | logical | .false. | Write a little extra about geometry to the screen. |