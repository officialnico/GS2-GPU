---
title: Timestepping Algorithm
author: Joseph Parker
date: 2017/04/21
---

[TOC]

# Introduction

The linear and collision terms are advanced using the [Kotschenreuther
algorithm][kotschenreuther1995] which is essentially a predictor-corrector method. In the
time advance steps, operator splitting is used to separate the advances of the linear
terms and the collision terms. This is a correct method, but it is not immediately obvious
why it works.

# Collisionless, linear algorithm

The linear terms are advanced using the [Kotschenreuther algorithm][kotschenreuther1995]
which is also described in [Highcock (2012)][highcock2012] and [Numata et al.
(2010)](https://arxiv.org/abs/1004.0279). This is an implicit method, allowing timesteps
comparable with ion timescales. The method also avoids inversions of matrices whose
dimensions are the problem size, requiring instead inversions of $N_{\vartheta}\times
N_{\vartheta}$ sized matrices. This is favourable for parallelization, but is the reason
that \(\vartheta\) is kept local to processor.

The algorithm is a predictor-corrector method, with two linear advances to update the
distribution function, and, between these, a field solve to update the electromagnetic
field. These steps are usually called "LFL" (see [CMR's note on layouts in
GS2](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/CMR/GS2_Layouts.pdf)).
Despite this name, the two `L` steps are *not* an operator splitting, as they do not
advance different terms in the gyrokinetic equation. Rather, "(LFL)" is a single operator
which returns the distribution function and fields at the next time step, given their
values at the current timestep.

The electrostatic, collisionless algorithm works as follows (this description is based on
[Highcock (2012)][highcock2012]. The discretized electrostatic gyrokinetic equation may be
written

\begin{equation}
A g^{n+1} + B g^n = D \phi^{n+1} + E\phi^{n},
\label{eq:gke}
\end{equation}

where \(g\) is the distribution function, \(\phi\) is the electrostatic potential, and \(A\),
\(B\), \(D\) and \(E\) are matrices. The superscript denotes the timestep. For simplicity, we
neglect gyroaveraging.

Similarly, we may write the quasineutrality condition as

\begin{equation}
F \phi^{n+1} = G g^{n+1}
\label{eq:QN}
\end{equation}

where again \(F\) and \(G\) are matrices. Here \(G\) represents the velocity space integrals in
Maxwell's equations, so \(G\) contracts with \(g^{n+1}\) to make a field-sized array.

To obtain the distribution function at the next time step, \(g^{n+1}\), we could combine
these equations

$$A g^{n+1} + B g^{n} = DF^{-1}Gg^{n+1} + E\phi^{n},$$

and solve explicitly for \(g^{n+1}\),

$$g^{n+1} = (A-DF^{-1}G)^{-1}( E\phi^{n} - B g^{n}),$$

but this involves inverting a matrix of size \(N_{\mathrm{tot}}^2\) (where
\(N_{\mathrm{tot}}=2N_{\vartheta}N_{s}N_{\lambda}N_{E}\)) at every timestep.

Instead, we use Kotschenreuther's method, which allows us to invert smaller matrices. It
is a predictor-corrector method which works as follows. We write

\begin{align}
\begin{split}
g^{n+1} &= g^P + g^C, \\
\phi^* &= \phi^{n+1} - \phi^n,
\end{split}
\label{eq:defs}
\end{align}

and substitute into the gyrokinetic equation to give,

$$A g^P + A g^C + B g^n = D \phi^* + D \phi^n + E\phi^n.$$

We then define the predictor \(g^P\) and corrector \(g^C\) to satisfy parts
of this equation independently,

\begin{align}
A g^P + B g^n &= D \phi^n + E\phi^n, \label{eq:gp}\\
A g^C  &= D {\phi}^* . \label{eq:gc}
\end{align}

Note that \(g^P\) depends only on values at the old timestep, while \(g^C\)
depends only on the change in the electrostatic potential. The latter
can be rewritten as

\begin{equation}
g^C_i = \left( \frac{\delta g}{\delta \phi} \right)_{ij} {\phi}^*_j,
\label{eq:response}
\end{equation}

that is, we can find the response matrix \((\delta g/\delta\phi)=A^{-1}D\)
explicitly.

To perform a time step:

1.  We first obtain \(g^P\) from \\eqref{eq:gp}.
2.  We then combine \\eqref{eq:defs} and \\eqref{eq:response} with the quasineutrality
    condition \\eqref{eq:QN} to obtain the updated field \(\phi^*\),
    
    \begin{align}
    F\phi^* + F\phi^n &= Gg^P + G \left(\frac{\delta g}{\delta \phi} \right) {\phi}^*, \nonumber\\
    {\phi}^* &= \left[ F - G \left( \frac{\delta g}{\delta \phi} \right)\right]^{-1} \left[ Gg^P - F\phi^n \right],
    \label{eq:FieldUpdate}
    \end{align}

    and therefore the field at the new timestep, \(\phi^{n+1} = \phi^n + \phi^*\).

3.  Finally we obtain \(g^{n+1}\) from the gyrokinetic equation \(\eqref{eq:gke}\), since we now know
    \(g^n\), \(\phi^n\) and \(\phi^{n+1}\).

## Implementation

The above three steps are referred to as "LFL". The two linear steps \(L\) are performed
by the function [[dist_fn:invert_rhs]] in [[dist_fn.fpp(file)]]. These steps are
algorithmically the same, since both can be cast as solving the gyrokinetic equation for
\(g^{n+1}\), given values for \(g^n\), \(\phi^{n+1}\) and \(\phi^n\). That is, the first
\(L\) step is the same as solving the gyrokinetic equation with \(\phi^{n+1}\) set equal
to \(\phi^n\). In practice, [[dist_fn:invert_rhs]] is always called as:

```Fortran
call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
```

with `g` as a global variable. Note that `gnew` is calculated in this function call, but
its input is not used.

Before the first linear step, the fields and distribution function are assigned to the
values at the end of the previous timestep:

```Fortran
g    = gnew
phi  = phinew
apar = aparnew
bpar = bparnew
```

so that `g` contains \(g^n\), and `phi` and `phinew` both contain \(\phi^n\), as required
for solving \\eqref{eq:gp}.

Before the second linear step, the fields are updated

```Fortran
phinew   = phinew  + phi
aparnew  = aparnew + apar
bparnew  = bparnew + bpar
```

where the field solve `getfield` had put the field updates \(\phi^*\) in `phinew`,
`aparnew` and `bparnew`. Therefore `phi` now contains \(\phi^n\), and `phinew` contains
\(\phi^{n+1}\). The distribution function `g` is not changed before the second linear
step, so `g` still contains \(g^n\). Thus calling [[dist_fn:invert_rhs]] now solves the
gyrokinetic equation.

These assignments are made in the function [[fields_implicit:advance_implicit]] in
[[fields_implicit.f90(file)]], which then calls [[dist_fn:timeadv]] in [[dist_fn.fpp]],
which in turn calls [[dist_fn:invert_rhs]].

Note that the corrector \(g^C\) is never explicitly evaluated. Note also that the
predictor \(g^P\) is only used to evaluate the field update in \(F\), but is not used
again.

## `g_lo` and `gf_lo` memory redistributes

The present default in GS2 is that the linear advances \(L\) and the field solve \(F\) are
calculated in the `g_lo` layout. This requires no memory redistributes, however the
velocity space integrals in the field solve are very expensive and appear to be the
bottleneck in typical production runs. This led to Adrian Jackson's eCSE project
\\cite{JacksoneCSE}, which introduced a new `gf_lo` layout for the field solve \(F\),
which improves the field solve at high core counts. This option is also available in the
GS2 trunk. The new approach comes at the cost of a memory redistribute between the first
\(L\) and \(F\), to transform `gnew` from `g_lo` to `gf_lo`.

Note that this only introduces one redistribute, between the first \(L\) and the
\(F\). There is no redistribute between \(F\) and the second \(L\), as the second \(L\)
uses `g` (which is in `g_lo`), not `gnew` (which is now in `gf_lo`).

# Collisions

The above timestepping method is attractive because it only involves inverting \(A\) and
\((F-G(\delta g/\delta\phi))\), matrices of size \(N_\vartheta^2\) rather than
\(N_{\textrm{tot}}^2\). Introducing collisions would introduce a velocity space dependence
into matrix \(A\), making it more expensive to invert.

Therefore collisions are implemented using operator splitting. Most of the papers which
describe the GS2 algorithm say that collisions are advanced by implicit Euler after doing
the whole collisionless time advance ([Kotschenreuther et al. 1995][kotschenreuther1995]
[Barnes et al. 2009](https://arxiv.org/pdf/0809.3945), [Highcock, 2012][highcock2012]),

\begin{align}
g^{n+1}_{\textrm{collisionless}} &= LFL( g^n, \phi^n)\\
g^{n+1} &= (I-\Delta t C)^{-1} g^{n+1}_{\textrm{collisionless}}.
\end{align}

However in the code, collisions are performed between updates in the implicit algorithm,
as per [Numata (2010)](https://arxiv.org/abs/1004.0279) and [Roach
(2016)](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/CMR/GS2_Layouts.pdf).
That is

\begin{align}
g^{n+1}_{\textrm{collisionless}} &= {L}( g^n, \phi^n, \phi^n) \label{eq:L1}\\
g^P &= (I-\Delta t C)^{-1} g^P_{\textrm{collisionless}}   \label{eq:C1}\\
\phi^{n+1} &= {F}( g^P, \phi^n)\\
g^{n+1}_{\textrm{collisionless}} &= {L}( g^n, \phi^n, \phi^{n+1})   \label{eq:L2}\\
g^{n+1} &= (I-\Delta t C)^{-1} g^{n+1}_{\textrm{collisionless}}  \label{eq:C2}
\end{align}

To see this in more detail, consider adding a collision term to the discretized
gyrokinetic equation, so that

\begin{equation}
A g^{n+1} + B g^n = D \phi^{n+1} + E\phi^n + C g^{n+1}
\label{eq:CollisionalGKE}
\end{equation}

and as before the quasineutrality condition is,

$$F\phi^{n+1} = Gg^{n+1}$$

Now inserting the predictor and corrector \\eqref{eq:defs} and decomposing the equation as
before, we have

\begin{align}
A g^P + B g^n &= D \phi^n + E\phi^n + Cg^P  ,\label{eq:Collisionalgp}\\ 
A g^C  &= D {\phi}^* + Cg^C \label{eq:Collisionalgc}
\end{align}

It is important to notice that, because \(A\) is a discretization of a time derivative, it
is legitimate to apply operator splitting to equations \\eqref{eq:CollisionalGKE},
\\eqref{eq:Collisionalgp} and \\eqref{eq:Collisionalgc} to advance the collision term
separately. Therefore \(g^P\) may be obtained as in equations \\eqref{eq:L1} and
\\eqref{eq:C1}, and \(g^{n+1}\) may be obtained as in \\eqref{eq:L2} and \\eqref{eq:C2}.

We must also ensure that collisions are included in the field solve. Notice that the
algebra in \\eqref{eq:FieldUpdate} assumes an expression for the corrector \(g^C\), but
once we apply collisions, this expression changes from \\eqref{eq:gc} to
\\eqref{eq:Collisionalgc}. However, this change is automatically accommodated in the
code. This is because the rows of the response matrix, \((\delta g/\delta\phi^*)\), are
calculated from \\eqref{eq:response} by setting \(\phi^*\) to be a delta function, and
noticing that the result of a timestep, \(g^C\), is the required row.  Therefore, all the
terms contained in the timestep \(L\) are automatically included in the response matrix.

For the field solve to be consistent, it is therefore only necessary for the timestep used
in calculating the response matrix to be the same as the timestep used in calculating the
predictor, \(g^P\). That is, in [[dist_fn:timeadv]], we must call the same routines for
`istep == 0` (the response matrix) as for `istep /= istep_last` (the predictor step).

## Counting memory redistributes

Introducing the collision term introduces four new memory redistributes,
in addition to the one redistribute in the collisionless case. In all,
we have five redistributes at the following points:

1.  Moving `gnew` from `g_lo` to `le_lo` to calculate collisions on
    the predictor \(g^P\).
2.  Moving `gnew` back from `le_lo` to `g_lo` after collisions.
3.  Moving `gnew` from `g_lo` to `gf_lo` to use the predictor in the
    field solve; there is no direct redistribute from `le_lo` to
    `gf_lo`.
4.  Moving `gnew` from `g_lo` to `le_lo` to calculate collisions on
    \(g^{n+1}\).
5.  Moving `gnew` back from `le_lo` to `g_lo` after collisions.

## Accounting for [[dist_fn_arrays:g_adjust]]

In the above discussion we have neglected the fact that collisions are applied on \\(h\\)
rather than \\(g\\). This does not affect the operator splitting in the Kotschenreuther
algorithm, but it does affect the splitting between the collisionless and collisional
advance.

## Fully-split collisions

A problem arises with [[dist_fn_arrays:g_adjust]] if collisions are fully separated from
the collisionless Kotschenreuther advance. The method may be summarized as

\begin{align}
  \begin{split}
    g^{n+1}_{\textrm{collisionless}}, \varphi^{n+1}_{\textrm{collisionless}}, \boldsymbol{A}^{n+1}_{\perp \textrm{collisionless}} &= LFL( g^n, \varphi^n) , \\
    h^{n+1}_{\textrm{collisionless}} &= g^{n+1}_{\textrm{collisionless}} + \frac{q_{0s}f_{0s}}{T_{0s}}\langle \varphi^{n+1}_{\textrm{collisionless}} - \boldsymbol{v}_{\perp}\cdot \boldsymbol{A}^{n+1}_{\perp \textrm{collisionless}} \rangle, \\
        h^{n+1} &= (I-\Delta t C)^{-1} h^{n+1}_{\textrm{collisionless}}, \\
        g^{n+1} &= h^P - \frac{q_{0s}f_{0s}}{T_{0s}}\langle \varphi^{n+1}_{\textrm{collisionless}} - \boldsymbol{v}_{\perp}\cdot \boldsymbol{A}^{n+1}_{\perp \textrm{collisionless}} \rangle, \\
  \varphi^{n+1}, \boldsymbol{A}_{\perp}^{n+1} &= \tilde{F}(g^{n+1}). 
  \end{split}
\end{align}

That is, collisionless versions of the new distribution function and fields are found from
a collisionless Kotschenreuther advance, then [[dist_fn_arrays:g_adjust]] is used to find
a collisionless \(h\), collisions are applied, [[dist_fn_arrays:g_adjust]] applied again
to find the new distribution function \(g^{n+1}\), and then finally a field solve is
called to find new fields consistent with the new distribution function.

This is the current approach in the "split collisions" version of GS2. However, it is not
quite correct, as the second [[dist_fn_arrays:g_adjust]] uses fields which are consistent
with \(h^{n+1}_\textrm{collisionless}\) not \(h^{n+1}\). The simplest correction to this
would be to find the new fields with a field solve on \(h\) before calling the second
[[dist_fn_arrays:g_adjust]],

\begin{align}
  \begin{split}
g^{n+1}_{\textrm{collisionless}}, \varphi^{n+1}_{\textrm{collisionless}}, \boldsymbol{A}^{n+1}_{\perp \textrm{collisionless}} &= LFL( g^n, \varphi^n) , \\
    h^{n+1}_{\textrm{collisionless}} &= g^{n+1}_{\textrm{collisionless}} + \frac{q_{0s}f_{0s}}{T_{0s}}\langle \varphi^{n+1}_{\textrm{collisionless}} - \boldsymbol{v}_{\perp}\cdot \boldsymbol{A}^{n+1}_{\perp \textrm{collisionless}} \rangle, \\
        h^{n+1} &= (I-\Delta t C)^{-1} h^{n+1}_{\textrm{collisionless}}, \\
\varphi^{n+1}, \boldsymbol{A}_{\perp}^{n+1} &= \bar{F}(h^{n+1}),\\
      g^{n+1} &= h^P - \frac{q_{0s}f_{0s}}{T_{0s}}\langle \varphi^{n+1} - \boldsymbol{v}_{\perp}\cdot \boldsymbol{A}^{n+1} \rangle. \\
    \end{split}
\end{align}

This method requires a new field solve \(\bar{F}\) to be written for \(h\), but like
\(\tilde{F}\) which is used in the initialization, it would not involve the response
matrix.

The new field solve would also only need to be called in electromagnetic runs, as
\(\varphi^{n+1}\) and \(\varphi^{n+1}_{\textrm{collisionless}}\) are (probably) the same
by conservation of mass.

# References

- M. Barnes, I. G. Abel, W. Dorland, D. R. Ernst, G. W. Hammett, P.  Ricci, B. N. Rogers,
  A. A. Schekochihin, and T. Tatsuno. Linearized model Fokker--Planck collision operators
  for gyrokinetic simulations. II. Numerical implementation and tests. [Physics of
  Plasmas, 16(7), 2009.](http://aip.scitation.org/doi/abs/10.1063/1.3155085) [arXiv
  link](https://arxiv.org/pdf/0809.3945)
- E. G Highcock. The Zero-Turbulence Manifold in Fusion Plasmas. PhD thesis, University of
  Oxford, 2012. [arXiv link][highcock2012]
- O. Kotschenreuther, G. Rewoldt, and W. M. Tang. Comparison of initial value and
  eigenvalue codes for kinetic toroidal plasma instabilities. [Computer Physics
  Communications, 88(2-3):128--140,
  1995.](http://www.sciencedirect.com/science/article/pii/001046559500035E) [pdf @
  PPPL](http://w3.pppl.gov/~hammett/gyrofluid/papers/1995/kotschenreuther95.pdf)
- P.  Numata, G. G. Howes, T. Tatsuno, M. Barnes, and W. Dorland. AstroGK: Astrophysical
  gyrokinetics code. [Journal of Computational Physics, 229:9347--9372,
  2010.](http://www.sciencedirect.com/science/article/pii/S0021999110005000) [arXiv
  link](https://arxiv.org/abs/1004.0279)
- C. M. Roach. On the need for and challenges associated with layouts in GS2, [Technical
  Note, 2016.](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/CMR/GS2_Layouts.pdf)


[kotschenreuther1995]: http://w3.pppl.gov/~hammett/gyrofluid/papers/1995/kotschenreuther95.pdf
[highcock2012]: https://arxiv.org/abs/1207.4419

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
