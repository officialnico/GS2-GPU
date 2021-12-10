---
title: Normalisations
---

[TOC]

Introduction
============

The perennial question: *Given a normalised quantity in gs2, what factor do I need to
multiply it by to get a real quantity?*

The answer should be answered somewhat by individual variable pages, but this page gives
an overview of the gs2 normalisation philosophy. Be sure to read the section
[Definitions](#definitions), containing definitions of some key normalizing quantities,
and also the [General Notes](#general-notes), before looking at the rest of the page. Once
you have read these sections, each definition of a normalization should be self-contained.

An extensive discussion of normalisations is in Chapter 3 of [this
thesis](http://arxiv.org/pdf/1207.4419v1).

See also Chapter 26 of [the AstroGK
manual](http://www.physics.uiowa.edu/~ghowes/astrogk/source/agk_manual.pdf).  Beware:
AstroGK does not have variation of the magnetic field magnetic along the field line, but
GS2 does. Important differences exist as a result.

There is also a good paper on G. Hammett's website introducing the [GS2 geometry
module](http://w3.pppl.gov/~hammett/work/gs2/docs/geometry/g_short.pdf), from which a lot
of this page is taken. However, be aware that there is a crucial sign error in the
definition of the equilibrium magnetic field (which also propagates into other work). Use
the form given on this page.

General Notes
-------------

-   All perturbed quantities are scaled up by \(a / \rho_{ia}\), the gyrokinetic epsilon
-   Unless otherwise stated, the subscript `s` is used to denote the species of interest,
    and the subscript `1` is used to denote the reference species, species 1. Note that
    species 1 in the GS2 input file does not necessarily set the reference mass, density,
    etc. These reference quantities are determined by the definition of mass, density,
    etc. in GS2 to be the ratio of the species quantity to the reference quantity. So, if
    dens for a given species is set to 1, then the reference density is equal to the
    density for that species.  The reference mass, temperature, etc. are defined in the
    same way.
-   The flux surface label \(\rho\) should not be confused with the species gyroradii:
    \(\rho_{sa}\)

Definitions
-----------

Here we present the definitions of some key normalising quantities. Note that many of
these quantities are defined on a given flux surface: this means normalizations are in
general dependent on the flux surface label.

-   \(a\) - All lengths are normalized to the quantity \(a\), which is chosen to be half
    the diameter of the last closed flux surface, measured on the midplane (the level of
    the magnetic axis). You can relate *a* to \(R_a\) using the input parameter
    [epsl](epsl), which is equal to \(2 a / R_a\)
-   \(R_a\) - \(R_a\) is the average of the maximum and minimum of R on the flux
    surface. With circular flux surfaces this is obviously the major radius at the center
    of the flux surface. GS2 users often use \(R_a\) and \(R_0\) interchangeably.
-   \(B_a\) - the toroidal field on the given flux surface evaluated at \(R_a\):

    $$B_a = \frac{I\left( \Psi \right) }{ R_a}$$

-   \(\Omega_a\) - the gyrofrequency of the reference species on the given flux surface
    evaluated at \(R_a\):

    $$\Omega_a = \frac{ \left| e \right| B_a }{ m_1 c}$$

-   \(\rho_{1a}\) - the gyroradius of the reference species on the given flux surface
    evaluated at \(R_a\):

    $$\rho_{1a} = \frac{v_{th1} }{ \Omega_a}$$

-   \(v_{th}=v_{th1}\) - the thermal velocity of species 1 (usually ions).  Note that the
    definition of velocity changes according to the value of [norm\_option](norm_option).

    $$v_{th}= v_{th1} = \sqrt{\frac{T_1}{m_1}} or \sqrt{\frac{2T_1}{m_1}}$$

Coordinates
===========

See main article [GS2 Coordinates](coordinates.html)

Flux Surface Label
------------------

The exact definition of the flux surface label \(\rho\) depends on the value of
[irho](irho). However, [irho](irho) is almost always set to 2, in which case

$$\rho = \frac{r}{a}$$

where *r* is the half-diameter of the current flux surface.

Lengths
-------

All lengths are normalized to the quantity *a*.

Time
----

Times are normalised to \(a/v_{th1}\)

$$t_{GS2} = \frac{t}{a/v_{th1}}$$

Geometrical Factors
===================

Safety Factor
-------------

$$q = \frac{d \Psi_T}{d \Psi}$$

Magnetic Shear
--------------

$$\hat{s} = \frac{\rho}{q}\frac{dq}{d\rho}$$

Note that if [irho](irho) = 2 (the default), then this is simply

$$\hat{s} = \frac{r}{q}\frac{dq}{dr}$$

where \(r\) is the half-diameter of the flux surface.

In the high aspect ratio, zero \(\beta\) circular flux surface limit,

$$\hat{s}= \frac{q R}{L_S}$$

Miscellaneous Parameters
========================

Flow Shear
----------

The flow shear is defined as

$$\gamma_E = \frac{\rho}{q}\frac{d \omega }{d \rho}$$

where \(\rho\) is the flux surface label and \(\omega\) is the toroidal angular
velocity. For historical, numerical and ordering reasons, in GS2 flow shear is specified
rather than the more intuitive \(\omega\) (when higher flow terms such as Coriolis and
centrifugal drifts were all absent, flow shear was the only term that appeared in the
equation).

If [irho](irho) = 2, then \(\rho = r / a\) and

$$\gamma_E = \frac{r}{q}\frac{d \omega}{d r}$$

Now the flow shear has units of inverse time, so it is normalised to \(v_{th1} / a\)

So

$$\gamma_{E(GS2)} = \frac{ \gamma_E}{v_{th1}/a} =  \frac{ 1}{v_{th1}/a} \frac{\rho}{q}\frac{d \omega }{d \rho}$$

Fields
======

Magnetic Field
--------------

The equilibrium field is defined as:

$$B_0 = I \left( \Psi \right) \nabla \phi - \nabla \psi \times \nabla \phi$$

The magnetic field is normalized to \(B_a\).

$$B_{GS2} = \frac{B}{B_a}$$

Magnetic Flux
-------------

The dimensions of \(\Psi\) are length \* length \* magnetic field so

$$\Psi_{GS2} = \frac{\Psi}{B_a a^2}$$

The Perturbed Field
-------------------

The perturbed field is normalised to \(\left| e \right| / T\) and is also scaled up by \(a
/ \rho_{ia}\), i.e. gyrokinetic epsilon. The perturbed generalised electromagnetic
potential \(\chi\) is defined as:

$$\chi = \phi - \frac{\mathbf{v}}{c} . \mathbf{A}$$

and with the normalizations we find that

$$\chi_{GS2} = \chi \frac{\left| e \right|}{T_1}\frac{a}{\rho_{1a}}$$

The Electric Potential
----------------------

$$\phi_{GS2} = \phi \frac{\left| e \right|}{T_1}\frac{a}{\rho_{1a}}$$

See [\#The Perturbed Field](#The_Perturbed_Field).

Fluxes
======

Heat Flux
---------

The radial heat flux for species *s* is defined (Abel et al., forthcoming) as

$$Q_s = \frac{1}{\left\langle \left| \nabla \psi \right| \right\rangle}  \left\langle \int d^3 \mathbf{v} \frac{m_s v^2}{2} \left( \mathbf{v}_E . \nabla \psi \right) \delta f_s \left( \mathbf{R} \right) \right\rangle$$

In GS2 it is output in spectral components as well as in aggregate.

It is normalised to the GyroBohm diffusion estimate (Cowley et al, 1991):

$$D_0 \sim \omega^{\star} \rho_{i}^2 \sim \frac{p_i}{l_n}v_{thi}\rho_i\frac{T_e}{T_i}$$

Replacing \(p_i\) by \(p_1\), the pressure of the reference species, species 1, and
replacing \(l_n\) by \(a\), dropping the temperature ratio and scaling up by
\(a/\rho_{1a}\) as is normal for perturbed quantities, we arrive at the normalisation:

$$Q_{GS2} = \frac {Q}{p_1 v_{th1}}\frac{a^2}{\rho_{1a}^2}$$

Toroidal Angular Momentum Flux
------------------------------

The radial component of the toroidal angular momentum flux for species \(s\) is defined
(Abel et al., forthcoming) as

$$\Pi_s = \frac{1}{\left\langle \left| \nabla \psi \right| \right\rangle}  \left\langle \int d^3 \mathbf{v} m_s R^2 \mathbf{v}. \nabla \phi\left( \mathbf{v}_E . \nabla \psi \right) \delta f_s \left( \mathbf{R} \right) \right\rangle$$

In GS2 it is output in spectral components as well as in aggregate.

It is normalised in an equivalent way to the heat flux:

$$\Pi_{GS2} = \frac {\Pi}{m_1 a n_1 v_{th1}^2}\frac{a^2}{\rho_{1a}^2}$$

Collisions
==========

The implementation of collisions in GS2/AstroGK is described in these two papers by Abel &
Barnes:

[Abel et. al](http://link.aip.org/link/?PHPAEN/15/122509/1) - Theory

[Barnes et. al](http://link.aip.org/link/?PHPAEN/16/072107/1) -
Numerical implementation.

Collision Frequency
-------------------

All collision frequencies are usually defined in terms of the like particle collision
frequency \(\nu\).

The like particle frequency for each species is specified by the input parameter `vnewk_s`
where \(s\) is the species index.

There are many definitions of the like particle frequency floating around. Here we display
some of them with comments.

**It is very important to note that choosing a collision frequency specifies the
temperature and density you are considering**. In a flux tube calculation, temperature and
density are scaled out, except in the collision frequency.

### GS2/AstroGK

The electron and ion collision frequencies used in GS2 and AstroGK are:

vnewk\_e =
$$\frac{a}{v_{thr}} \frac{\sqrt{2}\pi n_e e^4 ln \Lambda}{T_e^{3/2} m_e^{1/2}}$$

vnewk\_i =
$$\frac{a}{v_{thr}} \frac{\sqrt{2}\pi  Z_i^2 Z_{eff} n_e e^4 \ln \Lambda}{T_i^{3/2} m_i^{1/2}}$$

Note that in the second equation where you might expect to see \(Z^4e^4n_i\) you see
\(Z^2Z_{eff}n_e\). These formulate were taken from [a
note](http://w3.pppl.gov/~hammett/work/gs2/docs/collisions_notes.pdf) of Greg Hammett's on
the collision frequency, assuming that \(v_{thr}\) here is different from the symbols
\(v_{th}=v_{th1}\) defined earlier on this page (which depend on the choice of
\"norm\_option\") and is instead equivalent to \(\sqrt{2} v_t = \sqrt{2 T_{\mathrm
ref}/m_{\mathrm ref}}\) as used in Hammett\'s notes. Note that changing the choice of the
\"norm\_option\" variable (above) changes the normalizations of most of the inputs and
outputs involving a timescale by a factor of \(\sqrt{2}\), but I think that the definition
of vnewk\_e and vnewk\_i are unchanged because they are the collision frequencies
normalized by GS2\'s internal time scale \(a/\sqrt{2 T_{\mathrm ref}/m_{\mathrm ref}}\) .

### Abel et al.

The collision frequency in Abel et al, 2009, which presents the theory for the GS2
collision operator, is also the same as that of Trebnetov, and also Rosenbluth, Macdonald
& Judd

$$\nu_{ss} = \frac{\sqrt{2} \pi  Z_s^4 n_s e^4 \ln \Lambda}{T_s^{3/2} m_s^{1/2}}$$

Thus we see that there is no difference between that of Abel and that in GS2/AstroGK,
except the normalisation and the treatment of species charge and density.

### Rosenbluth Hazeltine & Hinton

The collision frequency in Rosenbluth Hazeltine & Hinton 1972 which derive neoclassical
transport coefficients is

$$\nu_{ss} = \frac{4\sqrt{2 \pi}  Z_s^4 n_s e^4 \ln \Lambda}{3T_s^{3/2} m_s^{1/2}}$$.

There is a factor of \(4\sqrt{\pi}/3\) between this and the GS2 frequency.

### Hinton & Wong

This is taken from Braginskii.

$$\nu_{ss} = \frac{\sqrt{2 \pi}  Z_s^4 n_s e^4 \ln \Lambda}{T_s^{3/2} m_s^{1/2}}$$

There is a factor of \(\sqrt{\pi}\) between this and the GS2 frequency.

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
