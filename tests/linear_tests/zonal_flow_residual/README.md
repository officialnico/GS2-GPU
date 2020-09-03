Zonal Flow Residual Test
========================

A basic test of the zonal flow residual: "linear collisionless
processes do not damp poloidal flows driven by ITG turbulence"
Rosenbluth-Hinton 1998, https://doi.org/10.1103/PhysRevLett.80.724

Zonal flows are linearly stable and are damped to a non-zero
residual. This tests that that residual level is correct. Actually,
follow Xiao-Catto 2007, https://doi.org/10.1063/1.2718519, who
expand to higher order in epsilon (inverse aspect ratio).

This is not a complete, thorough test, just a fast sanity-check.

Physical parameters are:
- $\epsilon$ (inverse aspect ratio) = 0.02
- $q$ (safety factor) = 1.3
- $k_\perp \rho_i$ (perpendicular wavenumber) = 0.002

The residual is the ratio of the potential at long-times ($\phi(t =
\infty)$) and the initial value. The predicted value is 0.0466.
