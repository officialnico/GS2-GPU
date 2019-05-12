!> FIXME : Add documentation
module theta_grid_params
  use abstract_config, only: abstract_config_type

  implicit none

  private

  public :: init_theta_grid_params
  public :: finish_theta_grid_params
  public :: wnml_theta_grid_params, write_trinity_parameters
  public :: rhoc, rmaj, r_geo, eps, epsl, qinp, shat, alpmhd
  public :: pk, geoType, aSurf, shift, shiftVert, mMode, nMode, deltam, deltan
  public :: deltampri, deltanpri, thetam, thetan
  public :: btor_slab, betaprim, ntheta, nperiod
  public :: set_overrides

  public :: theta_grid_parameters_config_type
  
  real :: rhoc, rmaj, r_geo, eps, epsl
  real :: qinp, shat, alpmhd, pk
  integer :: geoType, mMode, nMode
  real :: aSurf, shift, shiftVert, raxis, zaxis, deltam, deltan
  real :: akappa, tri, delta2, delta3, deltampri, deltanpri
  real :: akappri, tripri, thetam, thetan, thetak, thetad, theta2, theta3
  real :: btor_slab, betaprim
  real :: asym, asympri ! JRB - these two geometry parameters are NOT functional

  integer :: ntheta, nperiod

  logical :: initialized = .false.
  real :: kp = -1.
  logical :: exist
  
  !> Used to represent the input configuration of theta_grid
  type, extends(abstract_config_type) :: theta_grid_parameters_config_type
     ! namelist : theta_grid_parameters
     !> The flux surface elongation (only used when `geoType=0`, which selects the Generalised Miller analytic geometry specification).
     real :: akappa = 1.0
     !> The radial gradient flux surface elongation (only used when geoType=0, which selects the Generalised Miller analytic geometry specification).
     !> `akappri` = \(d\kappa/d\rho\)
     real :: akappri = 0.0
     !> Default = 0. Do not use unless you know what you are doing.
     real :: alpmhd = 0.0
     !> Minor radius of the flux surface that receives the specified
     !> shaping (only used when `geoType=1`, which selects the Global
     !> MHD analytic geometry specification). See Section 2.1.2 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: asurf = 1.0
     !> No longer functional.
     real :: asym = 0.0
     !> No longer functional.
     real :: asympri = 0.0
     !> In the slab limit, determines the angle between the field and
     !> the background flow (which flows in an imaginary toroidal
     !> direction). It is effectively equal to
     !> \(\frac{B_t}{B}=\frac{u_{\parallel}}{u}\).
     !> btor_slab = btor/bpol defines direction of a flow relative to B in slab 
     !> geometry, where flow is by definition in the toroidal direction.
     real :: btor_slab = 0.0
     !> The elongation of the flux surface labeled by aSurf (only used
     !> when `geoType`=1, which selects the Global MHD analytic geometry
     !> specifications). See Section 2.1.2 of the [Analytic Geometry
     !> Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: delta2 = 1.0
     !> The triangularity of the flux surface labeled by aSurf (only
     !> used when geoType=1, which selects the Global MHD analytic
     !> geometry specifications). See Section 2.1.2 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: delta3 = 1.0
     !> The magnitude of the mMode shaping effect (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: deltam = 1.0
     !> Radial derivative of the magnitude of the mMode shaping effect
     !> (only used when `geoType`=2 or 3, which selects the Generalised
     !> Elongation or Fourier analytic geometry specifications). See
     !> Sections 2.1.3 and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: deltampri = 0.0
     !> The magnitude of the nMode shaping effect (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: deltan = 1.0
     !> Radial derivative of the magnitude of the nMode shaping effect
     !> (only used when `geoType`=2 or 3, which selects the Generalised
     !> Elongation or Fourier analytic geometry specifications). See
     !> Sections 2.1.3 and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details. 
     real :: deltanpri = 0.0
     !> Controls particle trapping (among other things) in simple
     !> geometric models. \(\epsilon = r/R\)
     real :: eps = 0.3
     !>  \(\epsilon_\ell = \frac{2 L_{ref}}{ R}\) Sets curvature drift
     !>  in s-alpha model, where \(L_{ref}\) is the GS2
     !>  equilibrium reference normalisation length and \(R\) is the major
     !>  radius at the centre of the flux surface.
     real :: epsl = 0.3
     !> Selects between different analytic geometry specifications
     !> (only used when `iflux=0` and `local_eq=T`): `geoType=0` selects the
     !> Generalised Miller, `geoType=1` selects the Global MHD,
     !> `geoType=2` selects the Generalised Elongation, and `geoType=3`
     !> selects the Fourier specification. See Section 2.1 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     !> @note The default for this *MUST* be zero otherwise the Trinity interface will break.
     integer :: geotype = 0 
     !> `kp` sets parallel wavenumber and box length in slab geometry. \(k_p = \frac{2 \pi L_{ref}}{L_z}\).
     !> -  always use `kp` instead of `pk` in slab geometry (if `kp > 0` then gs2 sets `pk = 2*kp`)
     !> - in slab limit \(shat =  \frac{L_z}{2 \pi L_s} = \frac{L_{ref}}{k_p L_s}\) : nb if `kp = 1`, \(shat = \frac{L_{ref}}{L_s}), where \(L_s\) is the magnetic shear scale length.
     real :: kp = -1.0
     !> First flux surface shaping mode number (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     integer :: mmode = 2
     !> Second flux surface shaping mode number (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     integer :: nmode = 3
     !> Sets the number of \(2\pi\) segments along equilibrium
     !> magnetic field to include in simulation
     !> domain. \(N_{\mathrm segments} = (2n_{\mathrm period} -
     !> 1)\).
     !> - Ignored in some cases
     !> @todo Document when this is ignored
     integer :: nperiod = 2
     !> Rough number of grid points along equilibrium magnetic field between \(\theta=(-\pi,\pi)\).
     !> Actual number of grid points determined as follows:
     !> -  number of points in GS2 theta grid always odd because grid must contain both bounce points of trapped particles and one grid point at \(\theta=0\)
     !> -  theta grid in code runs from `-ntgrid:ntgrid`, with `ntgrid=Int(ntheta/2)` for `nperiod=1`.
     !> -  `ntheta` is ignored for some numerical equilibria
     integer :: ntheta = 24
     !> \(p_k = \frac{epsl}{q} = \frac{2 L_{ref}}{ q R}\) Sets q, the
     !> magnetic safety factor, and therefore also sets the connection
     !> length, i.e. the length of the box in the parallel direction,
     !> in terms of \(L_{ref}\). Used only in high aspect ratio
     !> \(s-\alpha\) equilibrium model.
     real :: pk = 0.3
     !> Sets value of the safety factor when using local toroidal equilibrium model.
     real :: qinp = 1.5
     !> When `iflux = 1`: Centerpoint of LCFS (normalized to
     !> \(L_{ref}\)) - When \(iflux \neq 1\): Major radius of
     !> magnetic field reference point (normalized to
     !> \(L_{ref}\)). Specifically, the reference magnetic field is
     !> defined to be the value of the toroidal magnetic field at `R=r_geo`
     !> on the flux surface labeled by `rhoc`.
     real :: r_geo = 3.0
     !> Major radial location of the magnetic axis (only used when
     !> `geoType`=1, which selects the Global MHD analytic geometry
     !> specification). See Section 2.1.2 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: raxis = 3.0
     !> `rhoc` is the flux surface label (0< `rhoc`< 1). Its exact
     !> meaning depends on `irho`. Usually rho = midplane
     !> diameter/midplane diameter of LCFS (here, "LCFS" refers to the
     !> flux surface at a radius of \(L_ref\). If using Miller geometry
     !> with \(L_{ref} \neq a\), this is \*not\* the physical last
     !> closed flux surface.
     !> - When `irho` = 1, `rhoc` = sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)
     !> - When `irho` = 2, `rhoc` = midplane diameter/(midplane diameter of LCFS). recommended
     !> - When `irho` = 3, `rhoc` = poloidal flux/(poloidal flux of LCFS)
     real :: rhoc = 0.5
     !> When \(iflux = 1\): Position of magnetic axis (normalized to \(L_{ref}\))
     !> - When \(iflux \neq 1\): Major radius of the centre of the flux surface of interest (normalized to \(L_{ref}\))
     real :: rmaj = 3.0
     !> Sets value of magnetic shear in simple geometric models.
     !> - over-ridden by `s_hat_input` in [[theta_grid_eik_knobs]] for most values of `bishop`.
     real :: shat = 0.75
     !> shift is related to minor radial derivatives of the major
     !>radial location of the flux surface centers (i.e. the Shafranov
     !>shift), but this input variable has **different physical
     !>definitions** in s-alpha and other analytic equilbrium models:
     !>
     !>
     !> -  In s-alpha (i.e. equilibrium_option='s-alpha'), shift\
     !>   \(\propto\) p' is a parameter for local
     !>   \(J_{\phi}) (and *NOT* \(B_{\theta}\) which is
     !>   constant): \(shift =
     !> \frac{2epsl}{pk^2}\frac{d\beta}{d\rho}=-\frac{q^2R}{L_{ref}}\frac{d\beta}{d\rho}
     !>   > 0\)
     !> -  In other analytic specifications
     !>   (i.e. equilibrium_option='eik' and `geoType`=0, 2, or 3),
     !>   shift is a parameter for local \(B_{\theta}\) (and NOT
     !>   for \(J_{\phi}\)): \(shift = \frac{1}{L_{ref}}
     !>   \frac{dR}{d\rho} < 0\)
     !>
     !> NB in Miller shift contains the *1st* radial derivative of the
     !> Shafranov shift, BUT in s-alpha shift is related to a *2nd*
     !> radial derivative of the Shafranov shift.
     !> - in s-alpha(i.e. equilibrium_option='s-alpha'), no additional parameter
     !> is required as the piece of \(J_{\phi} \propto\) p' is
     !> specified by shift.
     !> - in other analytic specifications (i.e. equilibrium_option='eik'), 
     !> an additional parameter(beta_prime) is required to specify the piece of
     !> \(J_{\phi} \propto\) p'
     real :: shift = 0.0
     !> Minor radial derivative of the axial location of the flux
     !> surface centers (i.e. the vertical Shafranov shift). It is
     !> only used when equilibrium_option='eik' and `geoType`=0, 2, or 3
     !> (which selects the Generalised Miller, Generalised Elongation,
     !> or Fourier analytic geometry specifications). See Sections
     !> 2.1.1, 2.1.3, and 2.1.4 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: shiftvert = 0.0
     !> The tilt angle of the elongation of the flux surface labeled
     !> by `aSurf` (only used when `geoType`=1, which selects the Global
     !> MHD analytic geometry specification). See Section 2.1.2 of the     
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: theta2 = 0.0
     !> The tilt angle of the triangularity of the flux surface
     !> labeled by `aSurf` (only used when `geoType`=1, which selects the
     !> Global MHD analytic geometry specification). See Section 2.1.2
     !> of the [Analytic Geometry
     !> Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: theta3 = 0.0
     !> The tilt angle of the triangularity (only used when `geoType`=0,
     !> which selects the Generalised Miller specification). See
     !> Section 2.1.1 of the     
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: thetad = 0.0
     !> The tilt angle of the elongation
     !> (only used when `geoType`=0, which selects the Generalised
     !> Miller specification). See Section 2.1.1 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: thetak = 0.0
     !> The tilt angle of the mMode shaping effect (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: thetam = 0.0
     !> The tilt angle of the nMode shaping effect (only used when
     !> `geoType`=2 or 3, which selects the Generalised Elongation or
     !> Fourier analytic geometry specifications). See Sections 2.1.3
     !> and 2.1.4 of the     
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: thetan = 0.0
     !>The flux surface triangularity (only used when `geoType`=0, which
     !>selects the Generalised Miller analytic geometry
     !>specification).  triangularity is `tri =
     !>arcsin[(R(max(Z))-R_major)/r_mid]`
     real :: tri = 0.0
     !> The radial gradient of the flux surface triangularity (only
     !> used when `geoType`=0, which selects the Generalised Miller
     !> analytic geometry specification). `tripri = dtri/dρ`
     real :: tripri = 0.0
     !> Axial location of the magnetic axis (only used when geoType=1,
     !> which selects the Global MHD analytic geometry
     !> specification). See Section 2.1.2 of the     
     !> [Analytic Geometry Specification](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
     !> documentation for more details.
     real :: zaxis = 0.0
   contains
     procedure, public :: read => read_theta_grid_parameters_config
     procedure, public :: write => write_theta_grid_parameters_config
     procedure, public :: reset => reset_theta_grid_parameters_config
     procedure, public :: broadcast => broadcast_theta_grid_parameters_config
  end type theta_grid_parameters_config_type

  type(theta_grid_parameters_config_type) :: config_internal
  
contains

  !> FIXME : Add documentation  
  subroutine init_theta_grid_params
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
!    logical, save :: initialized = .false.

    
    call debug_message(verb, "theta_grid_params::init_theta_grid_params start")
    if (initialized) return
    initialized = .true.

    call debug_message(verb, &
      "theta_grid_params::init_theta_grid_params call read_parameters")
    call read_parameters

    call debug_message(verb, "theta_grid_params::init_theta_grid_params end")
  end subroutine init_theta_grid_params

  !> FIXME : Add documentation
  subroutine finish_theta_grid_params
    implicit none
    initialized = .false.
  end subroutine finish_theta_grid_params

  !> FIXME : Add documentation
  subroutine read_parameters(config_in)
    use file_utils, only: input_unit, input_unit_exist
    use unit_tests, only: debug_message
    implicit none
    type(theta_grid_parameters_config_type), intent(in), optional :: config_in
    integer, parameter :: verb=3
    
    call debug_message(verb, "theta_grid_params::read_parameters start")
    
    
    if (present(config_in)) config_internal = config_in
    
    call config_internal%init(name = 'theta_grid_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    akappa = config_internal%akappa
    akappri = config_internal%akappri
    alpmhd = config_internal%alpmhd
    asurf = config_internal%asurf
    asym = config_internal%asym
    asympri = config_internal%asympri
    btor_slab = config_internal%btor_slab
    delta2 = config_internal%delta2
    delta3 = config_internal%delta3
    deltam = config_internal%deltam
    deltampri = config_internal%deltampri
    deltan = config_internal%deltan
    deltanpri = config_internal%deltanpri
    eps = config_internal%eps
    epsl = config_internal%epsl
    geotype = config_internal%geotype
    kp = config_internal%kp
    mmode = config_internal%mmode
    nmode = config_internal%nmode
    nperiod = config_internal%nperiod
    ntheta = config_internal%ntheta
    pk = config_internal%pk
    qinp = config_internal%qinp
    r_geo = config_internal%r_geo
    raxis = config_internal%raxis
    rhoc = config_internal%rhoc
    rmaj = config_internal%rmaj
    shat = config_internal%shat
    shift = config_internal%shift
    shiftvert = config_internal%shiftvert
    theta2 = config_internal%theta2
    theta3 = config_internal%theta3
    thetad = config_internal%thetad
    thetak = config_internal%thetak
    thetam = config_internal%thetam
    thetan = config_internal%thetan
    tri = config_internal%tri
    tripri = config_internal%tripri
    zaxis = config_internal%zaxis

    exist = config_internal%exist

    if (kp > 0.) pk = 2.*kp

    ! Print warning if the user has specified non-default values for the
    ! geometrical input parameters of other geometry types ! JRB
    select case (geoType)
       case (0)
          if (mMode/=2) write (*,*) "WARNING: ignoring value of mMode, value not needed"
          if (nMode/=3) write (*,*) "WARNING: ignoring value of nMode, value not needed"
          if (deltam/=1.0) write (*,*) "WARNING: ignoring value of deltam, did you mean akappa?"
          if (deltan/=1.0) write (*,*) "WARNING: ignoring value of deltan, did you mean tri?"
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean akappa?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean tri?"
          if (deltampri/=0.0) write (*,*) "WARNING: ignoring value of deltampri, did you mean akappri?"
          if (deltanpri/=0.0) write (*,*) "WARNING: ignoring value of deltanpri, did you mean tripri?"
          if (thetam/=0.0) write (*,*) "WARNING: ignoring value of thetam, did you mean thetak?"
          if (thetan/=0.0) write (*,*) "WARNING: ignoring value of thetan, did you mean thetad?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetak?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetad?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          deltam = akappa
          deltan = tri
          deltampri = akappri
          deltanpri = tripri
          thetam = thetak
          thetan = thetad
       case (1)
          if (mMode/=2) write (*,*) "WARNING: ignoring value of mMode, value not needed"
          if (nMode/=3) write (*,*) "WARNING: ignoring value of nMode, value not needed"
          if (deltam/=1.0) write (*,*) "WARNING: ignoring value of deltam, did you mean delta2?"
          if (deltan/=1.0) write (*,*) "WARNING: ignoring value of deltan, did you mean delta3?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean delta2?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean delta3?"
          if (deltampri/=0.0) write (*,*) "WARNING: ignoring value of deltampri, value not needed"
          if (deltanpri/=0.0) write (*,*) "WARNING: ignoring value of deltanpri, value not needed"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, value not needed"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, value not needed"
          if (thetam/=0.0) write (*,*) "WARNING: ignoring value of thetam, did you mean theta2?"
          if (thetan/=0.0) write (*,*) "WARNING: ignoring value of thetan, did you mean theta3?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean theta2?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean theta3?"
          if (shift/=0.0) write (*,*) "WARNING: ignoring value of shift, did you mean Raxis?"
          if (shiftVert/=0.0) write (*,*) "WARNING: ignoring value of shiftVert, did you mean Zaxis?"

          if (raxis>=Rmaj+aSurf) write (*,*) "WARNING: Raxis>=Rmag+aSurf, may specify non-nested flux surfaces"
          if (raxis<=Rmaj-aSurf) write (*,*) "WARNING: Raxis<=Rmag-aSurf, may specify non-nested flux surfaces"
          if (zaxis>=aSurf) write (*,*) "WARNING: Zaxis>=aSurf, may specify non-nested flux surfaces"
          if (zaxis<=-aSurf) write (*,*) "WARNING: Zaxis<=-aSurf, may specify non-nested flux surfaces"

          shift = raxis
          shiftVert = zaxis
          deltam = delta2
          deltan = delta3
          thetam = theta2
          thetan = theta3
       case (2)
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean deltam?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean deltan?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean deltam?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean deltan?"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, did you mean deltampri?"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, did you mean deltanpri?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetam?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetan?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean thetam?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean thetan?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          
       case (3)
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean deltam?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean deltan?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean deltam?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean deltan?"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, did you mean deltampri?"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, did you mean deltanpri?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetam?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetan?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean thetam?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean thetan?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select

    if (asym/=0.0) write (*,*) "WARNING: ignoring value of asym (not functional)"
    if (asympri/=0.0) write (*,*) "WARNING: ignoring value of asympri (not functional)"

  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine wnml_theta_grid_params(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_parameters"
    write (unit, fmt="(' ntheta =   ',i4)") ntheta
    write (unit, fmt="(' nperiod =  ',i4)") nperiod
    write (unit, fmt="(' rhoc =     ',f7.4)") rhoc
    write (unit, fmt="(' Rmaj =     ',f7.4)") rmaj
    write (unit, fmt="(' R_geo =    ',f7.4)") r_geo
    write (unit, fmt="(' eps =      ',f7.4)") eps
    write (unit, fmt="(' epsl =     ',f7.4)") epsl
    write (unit, fmt="(' qinp =     ',f7.4)") qinp
    write (unit, fmt="(' shat =     ',f7.4)") shat
    write (unit, fmt="(' alpmhd =   ',f7.4)") alpmhd
    write (unit, fmt="(' pk =       ',f7.4)") pk
    write (unit, fmt="(' kp =       ',f7.4)") kp
    write (unit, fmt="(' geoType =  ',i4)") geoType
    write (unit, fmt="(' aSurf =    ',f7.4)") aSurf
    write (unit, fmt="(' shift =    ',f7.4)") shift
    write (unit, fmt="(' shiftVert =',f7.4)") shiftVert
    write (unit, fmt="(' mMode =    ',i4)") mMode
    write (unit, fmt="(' nMode =    ',i4)") nMode
    write (unit, fmt="(' deltam =   ',f7.4)") deltam
    write (unit, fmt="(' deltan =   ',f7.4)") deltan
    write (unit, fmt="(' deltampri =',f7.4)") deltampri
    write (unit, fmt="(' deltanpri =',f7.4)") deltanpri
    write (unit, fmt="(' thetam =   ',f7.4)") thetam
    write (unit, fmt="(' thetan =   ',f7.4)") thetan
    write (unit, fmt="(' btor_slab =',f7.4)") btor_slab
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_params

  !> FIXME : Add documentation
  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    write (trinpars_unit, "(A22)") '&theta_grid_parameters'
    write (trinpars_unit, *) ' rhoc = ', rhoc
    write (trinpars_unit, *) ' qinp = ', qinp
    write (trinpars_unit, *) ' shat = ', shat
    write (trinpars_unit, *) ' rmaj = ', rmaj
    write (trinpars_unit, *) ' r_geo = ', r_geo
    write (trinpars_unit, *) ' geoType = ', geoType
    write (trinpars_unit, *) ' aSurf = ', aSurf
    write (trinpars_unit, *) ' shift = ', shift
    write (trinpars_unit, *) ' shiftVert = ', shiftVert
    write (trinpars_unit, *) ' mMode = ', mMode
    write (trinpars_unit, *) ' nMode = ', nMode
    write (trinpars_unit, *) ' deltam = ', deltam
    write (trinpars_unit, *) ' deltan = ', deltan
    write (trinpars_unit, *) ' deltampri = ', deltampri
    write (trinpars_unit, *) ' deltanpri = ', deltanpri
    write (trinpars_unit, *) ' thetam = ', thetam
    write (trinpars_unit, *) ' thetan = ', thetan
    write (trinpars_unit, *) ' betaprim = ', betaprim
    write (trinpars_unit, "(A1)") '/'
  end subroutine write_trinity_parameters

  !> FIXME : Add documentation
  subroutine set_overrides(mgeo_ov)
    use overrides, only: miller_geometry_overrides_type
    type(miller_geometry_overrides_type), intent(in) :: mgeo_ov
          !write (*,*) 'Calling tgpso'
    if (mgeo_ov%override_rhoc) rhoc = mgeo_ov%rhoc
    if (mgeo_ov%override_qinp) qinp = mgeo_ov%qinp
    if (mgeo_ov%override_shat) shat = mgeo_ov%shat
    if (mgeo_ov%override_rgeo_lcfs) r_geo = mgeo_ov%rgeo_lcfs
    if (mgeo_ov%override_rgeo_local) rmaj = mgeo_ov%rgeo_local
    if (mgeo_ov%override_geoType) geoType = mgeo_ov%geoType
    if (mgeo_ov%override_aSurf) aSurf = mgeo_ov%aSurf
    if (mgeo_ov%override_shift) shift = mgeo_ov%shift
    if (mgeo_ov%override_shiftVert) shiftVert = mgeo_ov%shiftVert
    if (mgeo_ov%override_mMode) mMode = mgeo_ov%mMode
    if (mgeo_ov%override_nMode) nMode = mgeo_ov%nMode
    if (mgeo_ov%override_deltam) then
      deltam = mgeo_ov%deltam
    else if (mgeo_ov%override_akappa) then
      deltam = mgeo_ov%akappa
    end if
    if (mgeo_ov%override_deltan) then
      deltan = mgeo_ov%deltan
    else if (mgeo_ov%override_tri) then
      deltan = mgeo_ov%tri
    end if
    if (mgeo_ov%override_deltampri) then 
      deltampri = mgeo_ov%deltampri
    else if (mgeo_ov%override_akappri) then 
      deltampri = mgeo_ov%akappri
    end if
    if (mgeo_ov%override_deltanpri) then 
      deltanpri = mgeo_ov%deltanpri
    else if (mgeo_ov%override_tripri) then 
      deltanpri = mgeo_ov%tripri
    end if
    if (mgeo_ov%override_thetam) thetam = mgeo_ov%thetam
    if (mgeo_ov%override_thetan) thetan = mgeo_ov%thetan
    if (mgeo_ov%override_betaprim) betaprim = mgeo_ov%betaprim
  end subroutine set_overrides


  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  
  !> Reads in the theta_grid_parameters namelist and populates the member variables
  subroutine read_theta_grid_parameters_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(theta_grid_parameters_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: geotype, mmode, nmode, nperiod, ntheta
    real :: akappa, akappri, alpmhd, asurf, asym, asympri, btor_slab, delta2, delta3, deltam, deltampri, deltan, deltanpri, eps, epsl, kp, pk, qinp, r_geo, raxis, rhoc, rmaj, shat, shift, shiftvert, theta2, theta3, thetad, thetak, thetam, thetan, tri, tripri, zaxis

    namelist /theta_grid_parameters/ akappa, akappri, alpmhd, asurf, asym, asympri, btor_slab, delta2, delta3, deltam, deltampri, deltan, deltanpri, eps, epsl, geotype, kp, mmode, nmode, nperiod, ntheta, pk, qinp, r_geo, raxis, rhoc, rmaj, shat, shift, shiftvert, theta2, theta3, thetad, thetak, thetam, thetan, tri, tripri, zaxis

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    akappa = self%akappa
    akappri = self%akappri
    alpmhd = self%alpmhd
    asurf = self%asurf
    asym = self%asym
    asympri = self%asympri
    btor_slab = self%btor_slab
    delta2 = self%delta2
    delta3 = self%delta3
    deltam = self%deltam
    deltampri = self%deltampri
    deltan = self%deltan
    deltanpri = self%deltanpri
    eps = self%eps
    epsl = self%epsl
    geotype = self%geotype
    kp = self%kp
    mmode = self%mmode
    nmode = self%nmode
    nperiod = self%nperiod
    ntheta = self%ntheta
    pk = self%pk
    qinp = self%qinp
    r_geo = self%r_geo
    raxis = self%raxis
    rhoc = self%rhoc
    rmaj = self%rmaj
    shat = self%shat
    shift = self%shift
    shiftvert = self%shiftvert
    theta2 = self%theta2
    theta3 = self%theta3
    thetad = self%thetad
    thetak = self%thetak
    thetam = self%thetam
    thetan = self%thetan
    tri = self%tri
    tripri = self%tripri
    zaxis = self%zaxis

    ! Now read in the main namelist
    if (self%requiresIndex) then
       call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
       read(in_file, nml = theta_grid_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%name, exist)
       if (exist) read(in_file, nml = theta_grid_parameters)
    endif

    ! Now copy from local variables into type members
    self%akappa = akappa
    self%akappri = akappri
    self%alpmhd = alpmhd
    self%asurf = asurf
    self%asym = asym
    self%asympri = asympri
    self%btor_slab = btor_slab
    self%delta2 = delta2
    self%delta3 = delta3
    self%deltam = deltam
    self%deltampri = deltampri
    self%deltan = deltan
    self%deltanpri = deltanpri
    self%eps = eps
    self%epsl = epsl
    self%geotype = geotype
    self%kp = kp
    self%mmode = mmode
    self%nmode = nmode
    self%nperiod = nperiod
    self%ntheta = ntheta
    self%pk = pk
    self%qinp = qinp
    self%r_geo = r_geo
    self%raxis = raxis
    self%rhoc = rhoc
    self%rmaj = rmaj
    self%shat = shat
    self%shift = shift
    self%shiftvert = shiftvert
    self%theta2 = theta2
    self%theta3 = theta3
    self%thetad = thetad
    self%thetak = thetak
    self%thetam = thetam
    self%thetan = thetan
    self%tri = tri
    self%tripri = tripri
    self%zaxis = zaxis

    self%exist = exist
  end subroutine read_theta_grid_parameters_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_theta_grid_parameters_config(self, unit)
    implicit none
    class(theta_grid_parameters_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("akappa", self%akappa, unit_internal)
    call self%write_key_val("akappri", self%akappri, unit_internal)
    call self%write_key_val("alpmhd", self%alpmhd, unit_internal)
    call self%write_key_val("asurf", self%asurf, unit_internal)
    call self%write_key_val("asym", self%asym, unit_internal)
    call self%write_key_val("asympri", self%asympri, unit_internal)
    call self%write_key_val("btor_slab", self%btor_slab, unit_internal)
    call self%write_key_val("delta2", self%delta2, unit_internal)
    call self%write_key_val("delta3", self%delta3, unit_internal)
    call self%write_key_val("deltam", self%deltam, unit_internal)
    call self%write_key_val("deltampri", self%deltampri, unit_internal)
    call self%write_key_val("deltan", self%deltan, unit_internal)
    call self%write_key_val("deltanpri", self%deltanpri, unit_internal)
    call self%write_key_val("eps", self%eps, unit_internal)
    call self%write_key_val("epsl", self%epsl, unit_internal)
    call self%write_key_val("geotype", self%geotype, unit_internal)
    call self%write_key_val("kp", self%kp, unit_internal)
    call self%write_key_val("mmode", self%mmode, unit_internal)
    call self%write_key_val("nmode", self%nmode, unit_internal)
    call self%write_key_val("nperiod", self%nperiod, unit_internal)
    call self%write_key_val("ntheta", self%ntheta, unit_internal)
    call self%write_key_val("pk", self%pk, unit_internal)
    call self%write_key_val("qinp", self%qinp, unit_internal)
    call self%write_key_val("r_geo", self%r_geo, unit_internal)
    call self%write_key_val("raxis", self%raxis, unit_internal)
    call self%write_key_val("rhoc", self%rhoc, unit_internal)
    call self%write_key_val("rmaj", self%rmaj, unit_internal)
    call self%write_key_val("shat", self%shat, unit_internal)
    call self%write_key_val("shift", self%shift, unit_internal)
    call self%write_key_val("shiftvert", self%shiftvert, unit_internal)
    call self%write_key_val("theta2", self%theta2, unit_internal)
    call self%write_key_val("theta3", self%theta3, unit_internal)
    call self%write_key_val("thetad", self%thetad, unit_internal)
    call self%write_key_val("thetak", self%thetak, unit_internal)
    call self%write_key_val("thetam", self%thetam, unit_internal)
    call self%write_key_val("thetan", self%thetan, unit_internal)
    call self%write_key_val("tri", self%tri, unit_internal)
    call self%write_key_val("tripri", self%tripri, unit_internal)
    call self%write_key_val("zaxis", self%zaxis, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_theta_grid_parameters_config

  !> Resets the config object to the initial empty state
  subroutine reset_theta_grid_parameters_config(self)
    class(theta_grid_parameters_config_type), intent(in out) :: self
    type(theta_grid_parameters_config_type) :: empty
    select type (self)
    type is (theta_grid_parameters_config_type)
       self = empty
    end select
  end subroutine reset_theta_grid_parameters_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_theta_grid_parameters_config(self)
    use mp, only: broadcast
    implicit none
    class(theta_grid_parameters_config_type), intent(in out) :: self
    call broadcast(self%akappa)
    call broadcast(self%akappri)
    call broadcast(self%alpmhd)
    call broadcast(self%asurf)
    call broadcast(self%asym)
    call broadcast(self%asympri)
    call broadcast(self%btor_slab)
    call broadcast(self%delta2)
    call broadcast(self%delta3)
    call broadcast(self%deltam)
    call broadcast(self%deltampri)
    call broadcast(self%deltan)
    call broadcast(self%deltanpri)
    call broadcast(self%eps)
    call broadcast(self%epsl)
    call broadcast(self%geotype)
    call broadcast(self%kp)
    call broadcast(self%mmode)
    call broadcast(self%nmode)
    call broadcast(self%nperiod)
    call broadcast(self%ntheta)
    call broadcast(self%pk)
    call broadcast(self%qinp)
    call broadcast(self%r_geo)
    call broadcast(self%raxis)
    call broadcast(self%rhoc)
    call broadcast(self%rmaj)
    call broadcast(self%shat)
    call broadcast(self%shift)
    call broadcast(self%shiftvert)
    call broadcast(self%theta2)
    call broadcast(self%theta3)
    call broadcast(self%thetad)
    call broadcast(self%thetak)
    call broadcast(self%thetam)
    call broadcast(self%thetan)
    call broadcast(self%tri)
    call broadcast(self%tripri)
    call broadcast(self%zaxis)

    call broadcast(self%exist)
  end subroutine broadcast_theta_grid_parameters_config
end module theta_grid_params
