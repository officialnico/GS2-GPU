---
title: Initialisation Options
---

[TOC]

GS2 has a lot of possible options for the initial condition. Of these, only a
few are usually used: these are the [typical](#typical-options) and
[restarting](#restarting-options)

## Typical options

These options are the most useful for ordinary GS2 simulations:

- `"default"`: Gaussian in \(\theta\). See [[init_g:ginit_default]]. Input
  parameters:

    - Amplitude: [[init_g_knobs:phiinit]]
    - Width in \(\theta\): [[init_g_knobs:width0]]
    - Parity: [[init_g_knobs:chop_side]] and [[init_g_knobs:left]]

- `"default_odd"`: Guassian in \(\theta\) multiplied by \(\sin(\theta -
  \theta_0)\). See [[init_g:ginit_default_odd]]. Input parameters:

    - Amplitude: [[init_g_knobs:phiinit]]
    - Width in \(\theta\): [[init_g_knobs:width0]]
    - Parity: [[init_g_knobs:chop_side]] and [[init_g_knobs:left]]

- `"noise"`: Noise along field line: each mode is given a random amplitude
  between zero and one. The recommended choice for nonlinear runs. See
  [[init_g:ginit_noise]]. Input parameters:

    - Amplitude: [[init_g_knobs:phiinit]]
    - Zonal flow factor: [[init_g_knobs:zf_init]]
    - Use a fixed set of random numbers: [[init_g_knobs:constant_random_flag]]
    - Continuity in parallel direction: [[init_g_knobs:clean_init]]
    - Single \(k_x\) mode: [[init_g_knobs:ikx_init]]
    - Parity: [[init_g_knobs:chop_side]] and [[init_g_knobs:left]]
    - Seed for random numbers: [[knobs:seed]]

- `"single_parallel_mode"`: Initialise only with a single parallel mode
  \(k_\Vert\) specified by either [[init_g_knobs:ikpar_init]] for periodic
  boundary conditions or [[init_g_knobs:kpar_init]] for linked boundary
  conditions. Only makes sense for linear calculations. See
  [[init_g:ginit_single_parallel_mode]]

### Less-usual options

There are many more initialisation options for GS2, although the majority are
unmaintained and are specific to the studies they were originally written
for. Unfortunately, many are also undocumented, so it is difficult to say
precisely what they do

- `"alf"`: Ion distribution function proportional to \(v_\Vert
  \sin(\theta)\); do not use.
- `"all_modes_equal"`: Initialise with every single parallel and
  perpendicular mode given the same amplitude. Intended for linear
  calculations
- `"cgyro"`: CGYRO-like initialisation, constant amplitude for \(k_y \ne
  0\) and zero amplitude for \(k_y == 0\) (see [[ginit_cgyro]] for
  details).
- `"cont"`: Restart, using current value of `g` in memory; unsupported,
  do not use.
- `"convect"`: Unsupported; do not use. Controlled by [[knobs:k0]].
- `"gs"`: Randomly phased \(k_\Vert = 1\) sines and cosines in density,
  parallel velocity, and parallel and perpendicular
  temperatures. Controlled by:

    - Overall amplitude: [[init_g_knobs:phiinit]]
    - Amplitude of real and imaginary parts of \(\phi\):
      [[init_g_knobs:refac]], [[init_g_knobs:imfac]]
    - Density perturbation: [[init_g_knobs:den1]]
    - Parallel velocity perturbation: [[init_g_knobs:upar1]]
    - Parallel temperature perturbation: [[init_g_knobs:tpar1]]
    - Perpendicular temperature perturbation: [[init_g_knobs:tperp1]]

- `"harris"`: Controlled by [[knobs:k0]]. Unsupported; do not use.
- `"kpar"`: Initial perturbation in fields of form \(x_0 + x_1
  \cos(\theta) + x_2 \cos(2 \theta)\), except for parallel velocity
  which uses is 90 degrees out of phase and uses \(\sin\). Density,
  parallel velocity, and perpendicular and parallel temperatures can
  all be controlled independently. Controlled by:

    - Overall amplitude: [[init_g_knobs:phiinit]]
    - Amplitude of real and imaginary parts of \(\phi\):
      [[init_g_knobs:refac]], [[init_g_knobs:imfac]]
    - Gaussian envelope width in \(\theta\): [[init_g_knobs:width0]]
    - Density perturbation: [[init_g_knobs:den0]],
      [[init_g_knobs:den1]], [[init_g_knobs:den2]]
    - Parallel velocity perturbation: [[init_g_knobs:upar0]],
      [[init_g_knobs:upar1]], [[init_g_knobs:upar2]]
    - Parallel temperature perturbation: [[init_g_knobs:tpar0]],
      [[init_g_knobs:tpar1]], [[init_g_knobs:tpar2]]
    - Perpendicular temperature perturbation: [[init_g_knobs:tperp0]],
      [[init_g_knobs:tperp1]], [[init_g_knobs:tperp2]]

- `"kz0"`: Constant along a field line, that is initialise only with
  \(k_\Vert = 0\).
- `"nl"`: Two \(k_\perp\) modes set to constant amplitude. Controlled by:

    - Overall amplitude: [[init_g_knobs:phiinit]]
    - `ky` indices of modes: [[init_g_knobs:ikk]]
    - `kx` indices of modes: [[init_g_knobs:itt]]
    - Parity: even about \(\theta_0\) unless
      [[init_g_knobs:chop_side]] is true, in which case none

- `"nl2"`: Two \(k_\perp\) modes proportional to \(1 + v_\Vert
  \sin(\theta)\). Controlled by:

    - Overall amplitude: [[init_g_knobs:phiinit]]
    - `ky` indices of modes: [[init_g_knobs:ikk]]
    - `kx` indices of modes: [[init_g_knobs:itt]]
    - Parity: even about \(\theta_0\) unless
      [[init_g_knobs:chop_side]] is true, in which case none

- `"nl3"`: Two \(k_\perp\) modes with a Gaussian envelope, and
  perturbation in fields of form \(x_0 + x_1 \cos(\theta) + x_2 \cos(2
  \theta)\), except for parallel velocity which uses is 90 degrees out
  of phase and uses \(\sin\). Density, parallel velocity, and
  perpendicular and parallel temperatures can all be controlled
  independently. Controlled by:

    - Overall amplitude: [[init_g_knobs:phiinit]]
    - Amplitude of real and imaginary parts of \(\phi\):
      [[init_g_knobs:refac]], [[init_g_knobs:imfac]]
    - Gaussian envelope width in \(\theta\): [[init_g_knobs:width0]]
    - `ky` indices of modes: [[init_g_knobs:ikk]]
    - `kx` indices of modes: [[init_g_knobs:itt]]
    - Parity: even about \(\theta_0\) unless
      [[init_g_knobs:chop_side]] is true, in which case none
    - Density perturbation: [[init_g_knobs:den0]],
      [[init_g_knobs:den1]], [[init_g_knobs:den2]]
    - Parallel velocity perturbation: [[init_g_knobs:upar0]],
      [[init_g_knobs:upar1]], [[init_g_knobs:upar2]]
    - Parallel temperature perturbation: [[init_g_knobs:tpar0]],
      [[init_g_knobs:tpar1]], [[init_g_knobs:tpar2]]
    - Perpendicular temperature perturbation: [[init_g_knobs:tperp0]],
      [[init_g_knobs:tperp1]], [[init_g_knobs:tperp2]]

- `"nl3r"`: Similar to `"nl3"`, but includes `apar`.
- `"nl4"`: Unsupported and inconsistent over time; do not use.
- `"nl5"`: Unsupported; do not use.
- `"nl6"`: Unsupported; do not use.
- `"nl7"`: Unsupported; do not use.
- `"ot"`: Orszag-Tang 2D vortex problem (see [[init_g_knobs:aparamp]],
  [[init_g_knobs:phiamp]], [[init_g_knobs:ittt]], [[init_g_knobs:ikkk]])
- `"rh"`: Maxwellian perturbation in `ik == 1` mode; do not use?
- `"stationary_mode"`: Initialises the distribution function such that
  it is a time-independent solution to the gyrokinetic equation when
  collisions are neglected and only one mode (with `aky = 0`) is
  present. See Section 4 of the [Analytic Geometry
  Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
  documentation for more details.
- `"recon"`: Unsupported; do not use.
- `"recon3"`: Unsupported; do not use.
- `"test3"`: Unsupported; do not use.
- `"xi"`: Perturbation proportional to pitch angle, \(\xi\). Unlikely to
  be useful.
- `"xi2"`: Perturbation proportional to \(1 - 3\xi^2\). Unlikely to be
  useful.
- `"zero"`: Sets distribution function to zero


## Restarting options

These options are used for restarting GS2. See [Restarting a
simulation](../user_manual/index.html#restarting-a-simulation) for more details.

- `"restart"`, `"many"`, or `"file"`: Restart, using `g` from restart files. See
  [[init_g:ginit_restart_many]].
- `"eig_restart"`: Uses the restart files written by the eigensolver. File to
  read is set by [[init_g_knobs:restart_eig_id]]. See [[init_g:ginit_restart_eig]]

If GS2 has been built with parallel netCDF, then `[[init_g_knobs:read_many]]`
controls whether a single restart file or one file per processor is used.

You can scale the fields from the restart files with
[[init_g_knobs:scale]]. Note that if `scale` is negative, the fields are scaled
by `-scale / max(abs(phi))`!

### Less-usual restarting options

These options are also used for restarting, but have additional behaviours. They
are not well documented or maintained, so use at your own risk.

- `"no_zonal"`: Restart but remove the zonal flow \((k_y = 0)\) component.
- `"small"`: Add restarted fields to a background as set by the `"noise"`
  option.
- `"smallflat"`: Similar to `"small"`, but with a different random background.
- `"zonal_only"`: Restart but set all non-zonal components \((k_y \ne 0)\) of
  the potential and the distribution function to 0. Input parameters:
  
    - Noise can be added to these other components by setting
      [[init_g_knobs:phiinit]] > 0.
    - The size of the zonal flows can be adjusted using
      [[init_g_knobs:zf_init]].
