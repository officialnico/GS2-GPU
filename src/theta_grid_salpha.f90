!> FIXME : Add documentation
module theta_grid_salpha
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: init_theta_grid_salpha, finish_theta_grid_salpha
  public :: check_theta_grid_salpha, wnml_theta_grid_salpha
  public :: salpha_get_sizes
  public :: salpha_get_grids

  public :: theta_grid_salpha_config_type
  public :: set_theta_grid_salpha_config
  public :: get_theta_grid_salpha_config

  ! knobs
  real :: alpmhdfac, alpha1

  ! internal variable
  integer :: model_switch
  integer, parameter :: model_salpha = 1, model_alpha1 = 2, &
       model_nocurve = 3, model_ccurv = 4, model_b2 = 5, &
       model_eps = 6, model_normal_only = 7
  real :: shift
  logical :: exist, initialized = .false.

  !> Used to represent the input configuration of theta_grid_salpha
  type, extends(abstract_config_type) :: theta_grid_salpha_config_type
     ! namelist : theta_grid_salpha_knobs
     ! indexed : false

     !> Used in conjunction with [[theta_grid_parameters:alpmhd]] to override
     !> `shift`, set as `shift=-alpmhd*alpmhdfac`.
     real :: alpmhdfac = 0.0
     !> Coefficient in model when `model_option='alpha1'` has been selected.
     real :: alpha1 = 0.0
     !> Sets the particular model for the magnetic field and related
     !> drifts. NB: All options have gbdrift = cvdrift except where
     !> noted. Can be one of
     !>
     !> - 's-alpha' - High aspect ratio toroidal equilibrium
     !> - 'default' - Same as 's-alpha'
     !> - 'alpha1' - Mainly same as 's-alpha' but with different
     !> definition of bmag and bset
     !> - 'rogers' - (aka. model_eps) From ingen output: "This model
     !> differs from the normal s-alpha model only in the curv and
     !> grad_B drifts." Indeed, cvdrift and gbdrift have an extra
     !> term, -(epsl*eps), while cvdrift0 and gbdrift0 are the same
     !> as 's-alpha'
     !> - 'b2' - From ingen output: "This model differs from the normal
     !> s-alpha model by an additional factor of 1/B(theta)**2 in the
     !> curv and grad_B drifts." Definition of bmag is also different.
     !> - 'normal_only' - Different definition of cvdrift (shat and
     !> shift terms removed) and cvdrift0 set to zero. Presumably this
     !> means that only the component of the curvature drift normal to
     !> the flux surface is retained while the component on the
     !> surface is 0. Useful for picking apart the effect of
     !> parameters on drifts and the effects of drifts on other
     !> quantities such as stability.
     !> - 'const-curv' - From ingen output: "Constant curvature is
     !> assumed. The grad-B and curvature drifts are both = epsl",
     !> i.e. shape = 'cylinder' NB: in contradiction to ingen output,
     !> gbdrift is not equal to cvdrift since cvdrift = epsl but
     !> gbdrift = cvdrift*(1.-shift). However, gbdrift0 = cvdrift0 =
     !> 0.
     !> - 'no-curvature' - From ingen output "Zero curvature is
     !> assumed", i.e. shape = 'slab'. NB: cvdrift = cvdrift0 =
     !> gbdrift0 = 0 but gbdrift is not 0 (gbdrift = epsl). NB: This
     !> does not yield the same result as cvdriftknob=0 in 's-alpha'
     !> model with non-zero epsl
     !>
     !> NB: For the final two options here ('const-curv' and 'no-curvature')
     !> (See also ingen output and contents of [[theta_grid.f90]] for further details
     character(len = 20) :: model_option = "default"
   contains
     procedure, public :: read => read_theta_grid_salpha_config
     procedure, public :: write => write_theta_grid_salpha_config
     procedure, public :: reset => reset_theta_grid_salpha_config
     procedure, public :: broadcast => broadcast_theta_grid_salpha_config
     procedure, public, nopass :: get_default_name => get_default_name_theta_grid_salpha_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_theta_grid_salpha_config
  end type theta_grid_salpha_config_type

  type(theta_grid_salpha_config_type) :: theta_grid_salpha_config

contains
  !> FIXME : Add documentation
  subroutine check_theta_grid_salpha(report_unit,alne,dbetadrho)
    use theta_grid_params, only: eps, epsl, pk, shat
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: alne, dbetadrho
!CMR input dbetadrho is computed externally (eg from species)
!    allowing consistency check
    real :: arat, qsf
!
 ! Find q, r/R, R/a
 !
    if (epsl > 0.) then
       arat = 2. / epsl

       if (epsl == 2.0) then
          write (report_unit, &
               & fmt="('Scale lengths are normalized to the major radius, R')")
       else
          write (report_unit, fmt="('The aspect ratio R/a = ',f7.4)") arat
          if (alne == 1.0) then
             write (report_unit, &
                  & fmt="('Scale lengths are normalized to the density scale length, Ln')")
          end if
       end if
       qsf = epsl/pk
       write (report_unit, fmt="('The safety factor q =      ',f7.4)") qsf
       write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") shat
       if (abs(shat) <= 1.e-5) then
          write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
       end if
       write (report_unit, fmt="('and epsilon == r/R = ',f7.4)") eps
       write (report_unit, *)
       if (eps > epsilon(0.0)) then
          write (report_unit, fmt="('Trapped particles are included.')")
       else
          write (report_unit, fmt="('Trapped particles are neglected.')")
       end if
       write (report_unit, *)

       if (shift > -epsilon(0.0)) then
          write (report_unit, fmt="('The s-alpha alpha parameter is ',f7.4)") shift
          !CMR 10/11/06: correct sign of dbeta/drho in s-alpha
          write (report_unit, fmt="('corresponding to d beta / d rho = ',f10.4)") -shift/arat/qsf**2
          !CMR 10/11/06: correct sign of dbeta/drho in s-alpha in this check
          if (abs(dbetadrho + shift/arat/qsf**2) > 1.e-2) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('This is inconsistent with beta and the pressure gradient.')")
             write (report_unit, fmt="('################# WARNING #######################')")
          end if
       else
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The s-alpha alpha parameter is less that zero.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
       end if

    else
       arat = 1.
       write (report_unit, &
            & fmt="('The radius of curvature is infinite.  This is a slab calculation.')")
    end if

    write (report_unit, *)
    select case (model_switch)

    case (model_salpha,model_b2,model_eps)
       if (epsl > 0.) then
          write (report_unit, fmt="('An s-alpha model equilibrium has been selected.')")
          write (report_unit, fmt="('The curvature and grad-B drifts are equal.')")
          write (report_unit, *)
          if (model_switch /= model_eps) then
             write (report_unit, fmt="('For theta0 = 0, each is of the form')")
             write (report_unit, *)
             write (report_unit, fmt="('  epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))')")
             write (report_unit, *)
          else
             write (report_unit, fmt="('For theta0 = 0, each is of the form')")
             write (report_unit, *)
             write (report_unit, fmt="('  epsl*(cos(theta) - eps + (shat*theta-shift*sin(theta))*sin(theta))')")
             write (report_unit, *)
          end if
          write (report_unit, fmt="('For finite theta0, there is also a term')")
          write (report_unit, *)
          write (report_unit, fmt="('  -epsl*shat*sin(theta)*theta0')")
          write (report_unit, *)
       end if
       write (report_unit, *)
       write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
       write (report_unit, *)
       write (report_unit, fmt="('  1.0 + (shat*theta-shift*sin(theta))**2')")
       write (report_unit, *)
       write (report_unit, fmt="('For finite theta0, there is also a term')")
       write (report_unit, *)
       write (report_unit, fmt="('  -shat*(shat*theta - shift*sin(theta))*theta0')")
       write (report_unit, *)
       write (report_unit, fmt="('and finally, the term')")
       write (report_unit, *)
       write (report_unit, fmt="('  shat**2 * theta0**2')")
       write (report_unit, *)
       if (model_switch == model_eps) then
          write (report_unit, *)
          write (report_unit, fmt="(' This model differs from the normal s-alpha model')")
          write (report_unit, fmt="(' only in the curv and grad_B drifts.')")
       end if
       if (model_switch == model_b2) then
          write (report_unit, *)
          write (report_unit, fmt="(' This model differs from the normal s-alpha model')")
          write (report_unit, fmt="(' by an additional factor of 1/B(theta)**2 (not shown above)')")
          write (report_unit, fmt="(' in the curv and grad_B drifts.')")
       end if
    case (model_ccurv)
       write (report_unit, fmt="('Constant curvature is assumed.')")
       write (report_unit, fmt="('The grad-B and curvature drifts are each = ',f10.4)") epsl
       write (report_unit, *)
       write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
       write (report_unit, *)
       write (report_unit, fmt="('  1.0 + (shat*theta-shift*sin(theta))**2')")
       write (report_unit, *)
       write (report_unit, fmt="('For finite theta0, there is also a term')")
       write (report_unit, *)
       write (report_unit, fmt="('  -shat*shat*theta*theta0')")
       write (report_unit, *)
       write (report_unit, fmt="('and finally, the term')")
       write (report_unit, *)
       write (report_unit, fmt="('  shat**2 * theta0**2')")
       write (report_unit, *)
    case (model_nocurve)
       write (report_unit, fmt="('Zero curvature is assumed.')")
       write (report_unit, *)
       write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
       write (report_unit, *)
       write (report_unit, fmt="('  1.0 + (shat*theta)**2')")
       write (report_unit, *)
       write (report_unit, fmt="('For finite theta0, there is also a term')")
       write (report_unit, *)
       write (report_unit, fmt="('  -shat*shat*theta*theta0')")
       write (report_unit, *)
       write (report_unit, fmt="('and finally, the term')")
       write (report_unit, *)
       write (report_unit, fmt="('  shat**2 * theta0**2')")
       write (report_unit, *)
    end select
  end subroutine check_theta_grid_salpha

  !> FIXME : Add documentation
  subroutine wnml_theta_grid_salpha(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_salpha_knobs"
    write (unit, fmt="(' alpmhdfac = ',e17.10)") alpmhdfac
    write (unit, fmt="(' alpha1 =    ',e17.10)") alpha1

    select case (model_switch)

    case (model_salpha)
       write (unit, fmt="(a)") ' model_option = "s-alpha"'

    case (model_alpha1)
       write (unit, fmt="(a)") ' model_option = "alpha1"'

    case (model_eps)
       write (unit, fmt="(a)") ' model_option = "rogers"'

    case (model_b2)
       write (unit, fmt="(a)") ' model_option = "b2"'

    case (model_normal_only)
       write (unit, fmt="(a)") ' model_option = "normal_only"'

    case (model_ccurv)
       write (unit, fmt="(a)") ' model_option = "const-curv"'

    case (model_nocurve)
       write (unit, fmt="(a)") ' model_option = "no-curvature"'

    end select
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_salpha

  !> FIXME : Add documentation
  subroutine init_theta_grid_salpha(theta_grid_salpha_config_in)
    use theta_grid_params, only: init_theta_grid_params
    implicit none
    type(theta_grid_salpha_config_type), intent(in), optional :: theta_grid_salpha_config_in

    if (initialized) return
    initialized = .false.

    call init_theta_grid_params

    ! Only required for trinity interface and not yet defined in s-alpha
    call read_parameters(theta_grid_salpha_config_in)
  end subroutine init_theta_grid_salpha

  !> FIXME : Add documentation
  subroutine finish_theta_grid_salpha
    implicit none
    initialized = .false.
    call theta_grid_salpha_config%reset()
  end subroutine finish_theta_grid_salpha

  !> FIXME : Add documentation
  subroutine read_parameters(theta_grid_salpha_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid_params, only: shift_in => shift, alpmhd
    use text_options, only: text_option, get_option_value
    implicit none
    type(theta_grid_salpha_config_type), intent(in), optional :: theta_grid_salpha_config_in

    character(20) :: model_option
    type (text_option), dimension (8), parameter :: modelopts = &
         (/ text_option('default', model_salpha), &
            text_option('s-alpha', model_salpha), &
            text_option('alpha1', model_alpha1), &
            text_option('rogers', model_eps), &
            text_option('b2', model_b2), &
            text_option('normal_only', model_normal_only), &
            text_option('const-curv', model_ccurv), &
            text_option('no-curvature', model_nocurve) /)

    integer :: ierr

    if (present(theta_grid_salpha_config_in)) theta_grid_salpha_config = theta_grid_salpha_config_in

    call theta_grid_salpha_config%init(name = 'theta_grid_salpha_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    alpha1 = theta_grid_salpha_config%alpha1
    alpmhdfac = theta_grid_salpha_config%alpmhdfac
    model_option = theta_grid_salpha_config%model_option

    exist = theta_grid_salpha_config%exist

    ierr = error_unit()
    call get_option_value &
         (model_option, modelopts, model_switch, &
         ierr, "model_option in theta_grid_salpha_knobs",.true.)

    if (alpmhdfac > epsilon(0.0)) then
       shift = - alpmhd*alpmhdfac
    else
       shift = shift_in
    end if

  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine salpha_get_sizes (nthetaout, nperiodout, nbsetout)
    use theta_grid_params, only: ntheta, nperiod
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound when alpha1 model is used
  end subroutine salpha_get_sizes

  !> FIXME : Add documentation
  subroutine salpha_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, shat, drhodpsi, kxfac, &
       qval, shape, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
    use constants, only: pi
    use theta_grid_params, only: eps, epsl, shat_param => shat, pk
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    use file_utils, only: error_unit
    use integration, only: trapezoidal_integration
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    real, intent (out) :: shat, drhodpsi, kxfac, qval, surfarea, dvdrhon, rhoc
    character (8), intent(out) :: shape
    logical, intent (in) :: gb_to_cv
    integer :: i

    theta = (/ (real(i)*2.0*pi/real(ntheta), i=-ntgrid,ntgrid) /)

! BD: dummy response for graphics in s-alpha mode until I have time to fix this:
    if (abs(epsl) > epsilon(0.)) then
       Rplot = 2./epsl*(1.+eps*cos(theta))  ; Rprime = 0.
    else
       Rplot = 1. ; Rprime = 0.
    end if
    Zplot = 1.  ; Zprime = 0.
    aplot = 1.  ; aprime = 0.

! MB : should look into changing this
    Bpol = 0.

    if (model_switch == model_alpha1) then
       bmag = 1.0-eps*cos(theta)-alpha1*cos(3.0*theta)
    else if (model_switch == model_b2) then
       bmag = 1.0 - eps*cos(theta)
    else
       bmag = 1.0/(1.0 + eps*cos(theta))
    end if

    shat = shat_param
    if (eps > epsilon(0.0)) then
       drhodpsi = 0.5*epsl**2/(pk*eps)
    else
       drhodpsi = 1.0
    end if
    kxfac = 1.0
    if (epsl > epsilon(0.0)) then
       qval = epsl/pk
    else
       qval = 1.
    end if
    select case (model_switch)
    case (model_salpha,model_alpha1,model_b2)
       cvdrift = epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))
       cvdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (model_switch == model_b2) then
          cvdrift = cvdrift/bmag**2
          cvdrift0 = cvdrift0/bmag**2
       end if
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0

    case (model_normal_only)
       cvdrift = epsl*cos(theta)
       cvdrift0 = 0.
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0

    case (model_eps)
       cvdrift = epsl*(cos(theta) -eps + (shat*theta-shift*sin(theta))*sin(theta))
       cvdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0

    case (model_ccurv,model_nocurve)
       cvdrift = epsl
       cvdrift0 = 0.0

! Some strangeness here to get straight at some point:
!    ccurv == constant curvature should be the case used for cylindrical
!             geometry, but evidently Paolo and Barrett do not like the
!             gds2 definition there, and have been using the slab
!             option (no_curvature) for their Z-pinch studies.
!
!    Simply need to look into the shift dependence of gds2
!
       if (model_switch == model_nocurve) then
!CMR, 4/6/2014:
! commented out gbdrift=0 as looked wrong, surely really want cvdrift=0
!dja fix for no curvature
!          gbdrift = 0.0
!dja end
!CMRend
          gds2 = 1.0 + (shat*theta)**2
          gds21 = -shat*shat*theta
          shape = 'slab    '
          gbdrift = cvdrift*(1.-shift)
          gbdrift0 = cvdrift0
          cvdrift=0 !CMR, 4/6/2014: surely this is what was intended?
       else
          gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
! probably should be:
!          gds2 = 1.0 + (shat*theta)**2
          gds21 = -shat*shat*theta
          shape = 'cylinder'
          gbdrift = cvdrift*(1.-shift)
          gbdrift0 = cvdrift0
       endif

       gds22 = shat*shat
       grho = 1.0

    end select
    gradpar = pk/2.0

    ! not sure about factor of epsl below...
    cdrift = 2.*epsl*(cos(theta)+shat*theta*sin(theta))
    cdrift0 = -2.*epsl*shat*sin(theta)
    ! BD: What are gds23 and gds24?  Who put this here?
    ! MB: gds23 and gds24 are geometrical factors appearing at next order in gk eqn
    ! MB: NEED TO INCLUDE SHIFT IN BELOW EXPRESSIONS
    !<DD> The following few lines will cause an issue in the (semi-)valid case where eps=0.0 so adding a guard
    !     here. These terms are used in lowflow calculations
    if(eps>epsilon(0.0))then
       gds23 = -0.5*epsl*shat*theta*(1.+2.*eps*cos(theta))/eps
       gds24_noq = 0.5*epsl*(1.+eps*cos(theta))/eps
    else
       write(error_unit(),'("Warning : Some lowflow related geometrical terms are forced to zero in cases with eps=0.")')
       gds23 = 0.
       gds24_noq = 0.
    endif
    gds24 = shat*gds24_noq

    if (model_switch /= model_alpha1) then
       bset = bmag(-ntheta/2:0)
    else
       call theta_grid_gridgen_init
       call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    end if

    ! Possibly crude approximations for s-alpha
    ! Note jacob is usually defined as jacob = 1.0/(drhodpsi*gradpar*bmag)
    ! but isn't yet stored here.
    surfarea = 2 * pi * trapezoidal_integration(theta, grho / (drhodpsi*gradpar*bmag))
    dvdrhon = 2 * pi * trapezoidal_integration(theta, 1.0 / (drhodpsi*gradpar*bmag))

    ! make rhoc consistent with eps, epsl
    if (epsl > epsilon(0.0)) then
       rhoc = 2.*eps/epsl
    else
       rhoc = 1.
    end if
  end subroutine salpha_get_grids

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_theta_grid_salpha_config(theta_grid_salpha_config_in)
    use mp, only: mp_abort
    type(theta_grid_salpha_config_type), intent(in), optional :: theta_grid_salpha_config_in
    if (initialized) then
       call mp_abort("Trying to set theta_grid_salpha_config when already initialized.", to_screen = .true.)
    end if
    if (present(theta_grid_salpha_config_in)) then
       theta_grid_salpha_config = theta_grid_salpha_config_in
    end if
  end subroutine set_theta_grid_salpha_config

  !> Get the module level config instance
  function get_theta_grid_salpha_config()
    type(theta_grid_salpha_config_type) :: get_theta_grid_salpha_config
    get_theta_grid_salpha_config = theta_grid_salpha_config
  end function get_theta_grid_salpha_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------

  !> Reads in the theta_grid_salpha_knobs namelist and populates the member variables
  subroutine read_theta_grid_salpha_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(theta_grid_salpha_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: model_option
    real :: alpha1, alpmhdfac

    namelist /theta_grid_salpha_knobs/ alpha1, alpmhdfac, model_option

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    alpha1 = self%alpha1
    alpmhdfac = self%alpmhdfac
    model_option = self%model_option

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = theta_grid_salpha_knobs)

    ! Now copy from local variables into type members
    self%alpha1 = alpha1
    self%alpmhdfac = alpmhdfac
    self%model_option = model_option

    self%exist = exist
  end subroutine read_theta_grid_salpha_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_theta_grid_salpha_config(self, unit)
    implicit none
    class(theta_grid_salpha_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("alpha1", self%alpha1, unit_internal)
    call self%write_key_val("alpmhdfac", self%alpmhdfac, unit_internal)
    call self%write_key_val("model_option", self%model_option, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_theta_grid_salpha_config

  !> Resets the config object to the initial empty state
  subroutine reset_theta_grid_salpha_config(self)
    class(theta_grid_salpha_config_type), intent(in out) :: self
    type(theta_grid_salpha_config_type) :: empty
    select type (self)
    type is (theta_grid_salpha_config_type)
       self = empty
    end select
  end subroutine reset_theta_grid_salpha_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_theta_grid_salpha_config(self)
    use mp, only: broadcast
    implicit none
    class(theta_grid_salpha_config_type), intent(in out) :: self
    call broadcast(self%alpha1)
    call broadcast(self%alpmhdfac)
    call broadcast(self%model_option)

    call broadcast(self%exist)
  end subroutine broadcast_theta_grid_salpha_config

  !> Gets the default name for this namelist
  function get_default_name_theta_grid_salpha_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_salpha_config
    get_default_name_theta_grid_salpha_config = "theta_grid_salpha_knobs"
  end function get_default_name_theta_grid_salpha_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_theta_grid_salpha_config()
    implicit none
    logical :: get_default_requires_index_theta_grid_salpha_config
    get_default_requires_index_theta_grid_salpha_config = .false.
  end function get_default_requires_index_theta_grid_salpha_config
end module theta_grid_salpha
