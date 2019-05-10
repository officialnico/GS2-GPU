title: theta_grid_parameters

#theta_grid_parameters

| Name | Type | Default | Description |
|:-----|:----:|:-------:|:------------|
| akappa | real | 1.0 | The flux surface elongation (only used when `geoType=0`, which selects the Generalised Miller analytic geometry specification). |
| akappri | real | 0.0 | The radial gradient flux surface elongation (only used when geoType=0, which selects the Generalised Miller analytic geometry specification).|
| | | |`akappri` = \(d\kappa/d\rho\) |
| alpmhd | real | 0.0 | Default = 0. Do not use unless you know what you are doing. |
| asurf | real | 1.0 | Minor radius of the flux surface that receives the specified|
| | | |shaping (only used when `geoType=1`, which selects the Global|
| | | |MHD analytic geometry specification). See Section 2.1.2 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| asym | real | 0.0 | No longer functional. |
| asympri | real | 0.0 | No longer functional. |
| btor_slab | real | 0.0 | In the slab limit, determines the angle between the field and|
| | | |the background flow (which flows in an imaginary toroidal|
| | | |direction). It is effectively equal to|
| | | |\(\frac{B_t}{B}=\frac{u_{\parallel}}{u}\).|
| | | |btor_slab = btor/bpol defines direction of a flow relative to B in slab|
| | | |geometry, where flow is by definition in the toroidal direction. |
| delta2 | real | 1.0 | The elongation of the flux surface labeled by aSurf (only used|
| | | |when `geoType`=1, which selects the Global MHD analytic geometry|
| | | |specifications). See Section 2.1.2 of the [Analytic Geometry|
| | | |Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| delta3 | real | 1.0 | The triangularity of the flux surface labeled by aSurf (only|
| | | |used when geoType=1, which selects the Global MHD analytic|
| | | |geometry specifications). See Section 2.1.2 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| deltam | real | 1.0 | The magnitude of the mMode shaping effect (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| deltampri | real | 0.0 | Radial derivative of the magnitude of the mMode shaping effect|
| | | |(only used when `geoType`=2 or 3, which selects the Generalised|
| | | |Elongation or Fourier analytic geometry specifications). See|
| | | |Sections 2.1.3 and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| deltan | real | 1.0 | The magnitude of the nMode shaping effect (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| deltanpri | real | 0.0 | Radial derivative of the magnitude of the nMode shaping effect|
| | | |(only used when `geoType`=2 or 3, which selects the Generalised|
| | | |Elongation or Fourier analytic geometry specifications). See|
| | | |Sections 2.1.3 and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| eps | real | 0.3 | Controls particle trapping (among other things) in simple|
| | | |geometric models. \(\epsilon = r/R\) |
| epsl | real | 0.3 | \(\epsilon_\ell = \frac{2 L_{ref}}{ R}\) Sets curvature drift|
| | | |in s-alpha model, where \(L_{ref}\) is the GS2|
| | | |equilibrium reference normalisation length and \(R\) is the major|
| | | |radius at the centre of the flux surface. |
| geotype | integer | 0  | Selects between different analytic geometry specifications|
| | | |(only used when `iflux=0` and `local_eq=T`): `geoType=0` selects the|
| | | |Generalised Miller, `geoType=1` selects the Global MHD,|
| | | |`geoType=2` selects the Generalised Elongation, and `geoType=3`|
| | | |selects the Fourier specification. See Section 2.1 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details.|
| | | |@note The default for this *MUST* be zero otherwise the Trinity interface will break. |
| kp | real | -1.0 | `kp` sets parallel wavenumber and box length in slab geometry. \(k_p = \frac{2 \pi L_{ref}}{L_z}\).|
| | | |-  always use `kp` instead of `pk` in slab geometry (if `kp > 0` then gs2 sets `pk = 2*kp`)|
| | | |- in slab limit \(shat =  \frac{L_z}{2 \pi L_s} = \frac{L_{ref}}{k_p L_s}\) : nb if `kp = 1`, \(shat = \frac{L_{ref}}{L_s}), where \(L_s\) is the magnetic shear scale length. |
| mmode | integer | 2 | First flux surface shaping mode number (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| nmode | integer | 3 | Second flux surface shaping mode number (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| nperiod | integer | 2 | Sets the number of \(2\pi\) segments along equilibrium|
| | | |magnetic field to include in simulation|
| | | |domain. \(N_{\mathrm segments} = (2n_{\mathrm period} -|
| | | |1)\).|
| | | |- Ignored in some cases|
| | | |@todo Document when this is ignored |
| ntheta | integer | 24 | Rough number of grid points along equilibrium magnetic field between \(\theta=(-\pi,\pi)\).|
| | | |Actual number of grid points determined as follows:|
| | | |-  number of points in GS2 theta grid always odd because grid must contain both bounce points of trapped particles and one grid point at \(\theta=0\)|
| | | |-  theta grid in code runs from `-ntgrid:ntgrid`, with `ntgrid=Int(ntheta/2)` for `nperiod=1`.|
| | | |-  `ntheta` is ignored for some numerical equilibria |
| pk | real | 0.3 | \(p_k = \frac{epsl}{q} = \frac{2 L_{ref}}{ q R}\) Sets q, the|
| | | |magnetic safety factor, and therefore also sets the connection|
| | | |length, i.e. the length of the box in the parallel direction,|
| | | |in terms of \(L_{ref}\). Used only in high aspect ratio|
| | | |\(s-\alpha\) equilibrium model. |
| qinp | real | 1.5 | Sets value of the safety factor when using local toroidal equilibrium model. |
| r_geo | real | 3.0 | When `iflux = 1`: Centerpoint of LCFS (normalized to|
| | | |\(L_{ref}\)) - When \(iflux \neq 1\): Major radius of|
| | | |magnetic field reference point (normalized to|
| | | |\(L_{ref}\)). Specifically, the reference magnetic field is|
| | | |defined to be the value of the toroidal magnetic field at `R=r_geo`|
| | | |on the flux surface labeled by `rhoc`. |
| raxis | real | 3.0 | Major radial location of the magnetic axis (only used when|
| | | |`geoType`=1, which selects the Global MHD analytic geometry|
| | | |specification). See Section 2.1.2 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| rhoc | real | 0.5 | `rhoc` is the flux surface label (0< `rhoc`< 1). Its exact|
| | | |meaning depends on `irho`. Usually rho = midplane|
| | | |diameter/midplane diameter of LCFS (here, "LCFS" refers to the|
| | | |flux surface at a radius of \(L_ref\). If using Miller geometry|
| | | |with \(L_{ref} \neq a\), this is \*not\* the physical last|
| | | |closed flux surface.|
| | | |- When `irho` = 1, `rhoc` = sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)|
| | | |- When `irho` = 2, `rhoc` = midplane diameter/(midplane diameter of LCFS). recommended|
| | | |- When `irho` = 3, `rhoc` = poloidal flux/(poloidal flux of LCFS) |
| rmaj | real | 3.0 | When \(iflux = 1\): Position of magnetic axis (normalized to \(L_{ref}\))|
| | | |- When \(iflux \neq 1\): Major radius of the centre of the flux surface of interest (normalized to \(L_{ref}\)) |
| shat | real | 0.75 | Sets value of magnetic shear in simple geometric models.|
| | | |- over-ridden by `s_hat_input` in [[theta_grid_eik_knobs]] for most values of `bishop`. |
| shift | real | 0.0 | shift is related to minor radial derivatives of the major|
| | | |radial location of the flux surface centers (i.e. the Shafranov|
| | | |shift), but this input variable has **different physical|
| | | |definitions** in s-alpha and other analytic equilbrium models:|
| | | ||
| | | ||
| | | |-  In s-alpha (i.e. equilibrium_option='s-alpha'), shift\|
| | | |\(\propto\) p' is a parameter for local|
| | | |\(J_{\phi}) (and *NOT* \(B_{\theta}\) which is|
| | | |constant): \(shift =|
| | | |\frac{2epsl}{pk^2}\frac{d\beta}{d\rho}=-\frac{q^2R}{L_{ref}}\frac{d\beta}{d\rho}|
| | | |> 0\)|
| | | |-  In other analytic specifications|
| | | |(i.e. equilibrium_option='eik' and `geoType`=0, 2, or 3),|
| | | |shift is a parameter for local \(B_{\theta}\) (and NOT|
| | | |for \(J_{\phi}\)): \(shift = \frac{1}{L_{ref}}|
| | | |\frac{dR}{d\rho} < 0\)|
| | | ||
| | | |NB in Miller shift contains the *1st* radial derivative of the|
| | | |Shafranov shift, BUT in s-alpha shift is related to a *2nd*|
| | | |radial derivative of the Shafranov shift.|
| | | |- in s-alpha(i.e. equilibrium_option='s-alpha'), no additional parameter|
| | | |is required as the piece of \(J_{\phi} \propto\) p' is|
| | | |specified by shift.|
| | | |- in other analytic specifications (i.e. equilibrium_option='eik'),|
| | | |an additional parameter(beta_prime) is required to specify the piece of|
| | | |\(J_{\phi} \propto\) p' |
| shiftvert | real | 0.0 | Minor radial derivative of the axial location of the flux|
| | | |surface centers (i.e. the vertical Shafranov shift). It is|
| | | |only used when equilibrium_option='eik' and `geoType`=0, 2, or 3|
| | | |(which selects the Generalised Miller, Generalised Elongation,|
| | | |or Fourier analytic geometry specifications). See Sections|
| | | |2.1.1, 2.1.3, and 2.1.4 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| theta2 | real | 0.0 | The tilt angle of the elongation of the flux surface labeled|
| | | |by `aSurf` (only used when `geoType`=1, which selects the Global|
| | | |MHD analytic geometry specification). See Section 2.1.2 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| theta3 | real | 0.0 | The tilt angle of the triangularity of the flux surface|
| | | |labeled by `aSurf` (only used when `geoType`=1, which selects the|
| | | |Global MHD analytic geometry specification). See Section 2.1.2|
| | | |of the [Analytic Geometry|
| | | |Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| thetad | real | 0.0 | The tilt angle of the triangularity (only used when `geoType`=0,|
| | | |which selects the Generalised Miller specification). See|
| | | |Section 2.1.1 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| thetak | real | 0.0 | The tilt angle of the elongation|
| | | |(only used when `geoType`=0, which selects the Generalised|
| | | |Miller specification). See Section 2.1.1 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| thetam | real | 0.0 | The tilt angle of the mMode shaping effect (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| thetan | real | 0.0 | The tilt angle of the nMode shaping effect (only used when|
| | | |`geoType`=2 or 3, which selects the Generalised Elongation or|
| | | |Fourier analytic geometry specifications). See Sections 2.1.3|
| | | |and 2.1.4 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |
| tri | real | 0.0 | The flux surface triangularity (only used when `geoType`=0, which|
| | | |selects the Generalised Miller analytic geometry|
| | | |specification).  triangularity is `tri =|
| | | |arcsin[(R(max(Z))-R_major)/r_mid]` |
| tripri | real | 0.0 | The radial gradient of the flux surface triangularity (only|
| | | |used when `geoType`=0, which selects the Generalised Miller|
| | | |analytic geometry specification). `tripri = dtri/dρ` |
| zaxis | real | 0.0 | Axial location of the magnetic axis (only used when geoType=1,|
| | | |which selects the Global MHD analytic geometry|
| | | |specification). See Section 2.1.2 of the|
| | | |[Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)|
| | | |documentation for more details. |