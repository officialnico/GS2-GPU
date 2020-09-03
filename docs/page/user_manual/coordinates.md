---
title: Coordinates
---

There are many different ways of indicating location within a fusion device, both local
(defined with respect to a given flux surface, usually used to indicate location within a
flux tube, e.g.  \(x, y, \theta\) in GS2) and global (denoting position within the torus,
e.g. \(R, Z, \phi\) or \(r, \theta, \phi\) or \(\psi, \theta, \phi\)).

As a general rule, in GS2, local coordinates (or their spectral transforms) are used to
map perturbed quantities such as the potential, and to evolve the gyrokinetic equations,
whereas global quantities are used to specify the location of the flux tube within the
torus, and to calculate the value of the spatially varying drifts such as the curvature
drift. Thus, various sets of coordinates are used at different times, necessitating
conversions between them. Furthermore, no set of GS2 coordinates actually corresponds to a
real space set of coordinates, so if it is desired to have fields output as functions of
real space, it is necessary to work out the correspondence between GS2 coordinates and
real space coordinates.

GS2 Coordinates
---------------

The "native" coordinates of GS2 are \(x,y,\theta\), a set of field line following local
coordinates. Each of these coordinates are defined below. Their definitions are obscure,
and are therefore reduced to physically meaningful expressions for simple magnetic
geometries where this is possible.

### x

This coordinate denotes the flux surface, and is commonly known as the radial coordinate,
though this is only the case when the magnetic flux surfaces are concentric circles. Its
actual definition is

\begin{equation}
x = a \frac{q_0}{\rho_0} \left( \Psi_{0N} - \Psi_{N} \right)
\end{equation}

where a is the reference length, \(q_0\) is the magnetic safety factor on the reference
flux surface, \(\rho_0\) the normalised flux surface label and \(\Psi_N\) the normalised
flux.

#### The Concentric Circular Flux Surface Limit

In the case where the magnetic flux surfaces are concentric circles and \(\rho_n = r_0/a\)
and \(a=R_0\), we may write:

\begin{align}
x &= a \frac{q_0}{\rho_0} \frac{d \Psi_N}{d \rho} \left( \rho_0 - \rho\right) \\
&= \frac{a^2q_0}{r_0}\frac{d \Psi_N}{d r} \left( r_0 - r \right) \\
&= \frac{R_0^2q_0}{r_0}\frac{d \Psi_N}{d r} \left( r_0 - r \right)
\end{align}

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
