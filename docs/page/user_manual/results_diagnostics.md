---
title: Results and Diagnostics
---

[TOC]

Introduction
============

GS2 can output many different quantities. Many are written into the netCDF output file,
and unfortunately many more are written into a plethora of undocumented ascii files. Some
are always written out, and some are only written out if you set one of a number of flags
in the namelist diagnostics. This section tries to make sense of them.

There are perhaps 2 important things we want from GS2 output:

1.  Physical quantities - results of the simulation.
2.  Diagnostic quantities - quantities that determine whether the runs are numerically
    resolved, but are not of interest physically.

These are dealt with separately. Though sometimes the distinction can seem a little
arbitrary, as most diagnostic quantities are in fact physical quantities, hopefully it
makes it easier to get the full picture.

Physical Quantities
===================

Input Parameters
----------------

Suppose you don't have the input file, how do you find out what the input parameters were?
The simple answer is

``` shell
ncdump -c run_name.out.nc
```

Geometrical Factors
-------------------

In the various toroidal geometries, GS2 calculates all sorts of interesting quantities --
\(R\), \(q\), and so on. How do you find out what they were?

Grids
-----

What values of [theta](theta), [ky](ky), [ky](ky) are the physical quantities given for?

Again the simple answer is

``` shell
ncdump -c run_name.out.nc
```

Fields
------

Fluxes
------

Until this section is fleshed out try

``` shell
grep fluxes run_name.out
```

The Distribution Function
-------------------------

Diagnostic Quantities
=====================

How do I know if my run is resolved ?
-------------------------------------

There is no complete answer to this question, but here is a rough sketch. L stands for
linear runs, N stands for non-linear.

1.  Unary Checks: look at properties of that run:
    1.  Eigenfunction:
        1.  [shat](shat) non-zero
            1.  Plot the eigenfunctions for different [ky](ky) on the extended grid
                (instructions coming soon). They should go to 0 at either edge. (L&N)
    2.  \(\phi^2\) vs time (on a log-linear plot).
        1.  Linearly should be a straight line (exponential growth).
        2.  Non-linearly should initially be a straight line and then saturate at some
            level (a very crude picture - infinite possibilities!).
    3.  Velocity space
        1.  Check the velocity space diagnostics: the fraction of energy in the high k
            velocity space modes should be around 10% and should certainly not go
            to 1. (L&N)
2.  Binary Checks: run another run with different resolution parameters and compare the
    two.
    1.  Eigenfunctions
        1.  Linearly, the normalised magnitude should be approximately the same.

Velocity Space Diagnostics
--------------------------

Set [[gs2_diagnostics_config_type:write_verr]] to true.

In the file `run_name.lpc` the columns are:

1.  Time
2.  Fraction of energy in high k for pitch angle harmonics (if this is too big,
    [[le_grids_config_type:ngauss]] may be too small).
3.  Fraction of energy in high k for energy harmonics (if this is too big,
    [[le_grids_config_type:negrid]] may be too small).

In the file `run_name.verr` the columns are:

1.  Time
2.  Fraction of energy in high k for pitch angle harmonics (if this is too big,
    [[le_grids_config_type:ngauss]] may be too small).
3.  Fraction of energy in high k for energy harmonics (if this is too big,
    [[le_grids_config_type:negrid]] may be too small).

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
