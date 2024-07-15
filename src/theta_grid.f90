!> FIXME : Add documentation
module theta_grid
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: init_theta_grid, finish_theta_grid, check_theta_grid, wnml_theta_grid
  public :: theta, theta2, delthet, delthet2
  public :: bset, bmag, gradpar, itor_over_B, IoB, kxfac, qval, grho
  public :: gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0
  public :: gds2, gds21, gds22, gds23, gds24, gds24_noq
  public :: bmin, bmax, eps_trapped, shat, drhodpsi, jacob
  public :: ntheta, ntgrid, nperiod, nbset, shape, gb_to_cv, surfarea, dvdrhon
  public :: Rplot, Zplot, aplot, Rprime, Zprime, aprime, Bpol, rhoc
  public :: initialized, is_effectively_zero_shear, field_line_average

  public :: theta_grid_config_type
  public :: set_theta_grid_config
  public :: get_theta_grid_config

  real, dimension (:), allocatable :: theta, theta2, delthet, delthet2
  real, dimension (:), allocatable :: bset
  real, dimension (:), allocatable :: bmag, gradpar
  real, dimension (:), allocatable :: itor_over_B, IoB
  real, dimension (:), allocatable :: gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0
  real, dimension (:), allocatable :: gds2, gds21, gds22, gds23, gds24, gds24_noq
  real, dimension (:), allocatable :: grho, jacob
  real, dimension (:), allocatable :: Rplot, Zplot, aplot, Bpol
  real, dimension (:), allocatable :: Rprime, Zprime, aprime
  real :: bmin, bmax, eps_trapped, shat, drhodpsi, kxfac, qval
  real :: cvdriftknob, gbdriftknob
  real :: surfarea, dvdrhon, rhoc
  integer :: ntheta, ntgrid, nperiod, nbset
  logical :: gb_to_cv

  real, parameter :: smallest_non_zero_shear = 1.0e-5

  ! internal variables
  integer :: eqopt_switch
  integer, parameter :: eqopt_eik = 1, eqopt_salpha = 2, eqopt_file = 3, eqopt_file_nc = 4
  character (8) :: shape
  logical :: exist, initialized = .false.
  real :: field_line_average_weight

  interface field_line_average
     module procedure field_line_average_real
     module procedure field_line_average_complex
  end interface field_line_average

  !> Used to represent the input configuration of theta_grid
  type, extends(abstract_config_type) :: theta_grid_config_type
     ! namelist : theta_grid_knobs
     ! indexed : false
     !> The equilibrium_option variable controls which geometric
     !> assumptions are used in the run. Additional namelists
     !> determine the geometric parameters according to the choice
     !> made here. Allowed values are:
     !>
     !> - 'eik' Use routines from the geometry module, controlled mainly
     !> by the subsidiary namelists theta_grid_parameters and
     !> theta_grid_eik_knob. This includes options for Miller as well
     !> as a range of different numerical equilibrium file formats.
     !> - 'default' Same as 'eik'
     !> - 's-alpha' Use high aspect-ratio toroidal equilbrium (which can
     !> be simplified to slab or cylinder), controlled by the
     !> subsidiary namelist theta_grid_parameters and
     !> theta_grid_salpha_knob
     !> - 'file' Use output from rungridgen code directly. Controlled by
     !> theta_grid_file_knobs. Note: in this case, the variables
     !> nperiod and ntheta (from the theta_grid_parameters namelist)
     !> are ignored.
     !> - 'grid.out' Same as 'file'
     !> - 'ncfile' Similar to 'file' but using a netcdf file input
     !> instead of ascii.
     !>
     character(len = 20) :: equilibrium_option = "default"
     !> Scales the grad-B drift.
     real :: gbdriftknob = 1.0
     !> Scales the curvature drift.
     real :: cvdriftknob = 1.0
     !> If true then force grad-B drift to be equal to curvature
     !> drift. This is not recommended when `fbpar` is not 0.
     logical :: gb_to_cv = .false.
   contains
     procedure, public :: read => read_theta_grid_config
     procedure, public :: write => write_theta_grid_config
     procedure, public :: reset => reset_theta_grid_config
     procedure, public :: broadcast => broadcast_theta_grid_config
     procedure, public, nopass :: get_default_name => get_default_name_theta_grid_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_theta_grid_config
  end type theta_grid_config_type

  type(theta_grid_config_type) :: theta_grid_config

contains

  !> FIXME : Add documentation
  subroutine check_theta_grid(report_unit, alne, dbetadrho)
    use theta_grid_salpha, only: check_theta_grid_salpha
    use theta_grid_eik, only: check_theta_grid_eik
    use theta_grid_file, only: check_theta_grid_file
    use theta_grid_file, only: check_theta_grid_file_nc
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: alne, dbetadrho
    select case (eqopt_switch)
    case (eqopt_salpha)
       call check_theta_grid_salpha(report_unit, alne, dbetadrho)
    case (eqopt_eik)
       call check_theta_grid_eik(report_unit,dbetadrho)
    case (eqopt_file)
       call check_theta_grid_file(report_unit)
    case (eqopt_file_nc)
       call check_theta_grid_file_nc(report_unit)
    end select
    if (gb_to_cv) then
       write (report_unit, *) 'The grad B drift coefficients have been set equal to the'
       write (report_unit, *) 'values for the curvature drift coefficients.  Do not use'
       write (report_unit, *) 'fbpar = 1.0 in this case.'
       write (report_unit, *)
       write (report_unit, *) 'You got this option by setting gb_to_cv = .true.'
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You have chosen to set the grad B drift equal to the curvature drift.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    end if
  end subroutine check_theta_grid

  !> FIXME : Add documentation
  subroutine wnml_theta_grid(unit)
    use theta_grid_salpha, only: wnml_theta_grid_salpha
    use theta_grid_eik, only: wnml_theta_grid_eik
    use theta_grid_file, only: wnml_theta_grid_file
    use theta_grid_gridgen, only: wnml_theta_grid_gridgen
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_knobs"

    select case (eqopt_switch)
    case (eqopt_eik)
       write (unit, fmt="(a)") ' equilibrium_option = "eik"'
    case (eqopt_salpha)
       write (unit, fmt="(a)") ' equilibrium_option = "s-alpha"'
    case (eqopt_file)
       write (unit, fmt="(a)") ' equilibrium_option = "file"'
    case (eqopt_file_nc)
       write (unit, fmt="(a)") ' equilibrium_option = "ncfile"'
    end select
    write (unit, fmt="(' gb_to_cv = ',L1)") gb_to_cv
    write (unit, fmt="(' /')")
    select case (eqopt_switch)
    case (eqopt_salpha)
       call wnml_theta_grid_salpha(unit)
    case (eqopt_eik)
       call wnml_theta_grid_eik(unit)
    case (eqopt_file)
       call wnml_theta_grid_file(unit)
    case (eqopt_file_nc)
       call wnml_theta_grid_file(unit)
    end select
    call wnml_theta_grid_gridgen(unit)
  end subroutine wnml_theta_grid

  !> FIXME : Add documentation
  subroutine init_theta_grid(theta_grid_config_in, theta_grid_gridgen_config_in, &
       &theta_grid_salpha_config_in, theta_grid_file_config_in, theta_grid_eik_config_in)
    use mp, only: proc0
    use unit_tests, only: debug_message
    use theta_grid_gridgen, only: theta_grid_gridgen_config_type
    use theta_grid_salpha, only: theta_grid_salpha_config_type
    use theta_grid_file, only: theta_grid_file_config_type
    use theta_grid_eik, only: theta_grid_eik_config_type
    implicit none
    type(theta_grid_config_type), intent(in), optional :: theta_grid_config_in
    type(theta_grid_gridgen_config_type), intent(in), optional :: theta_grid_gridgen_config_in
    type(theta_grid_salpha_config_type), intent(in), optional :: theta_grid_salpha_config_in
    type(theta_grid_file_config_type), intent(in), optional :: theta_grid_file_config_in
    type(theta_grid_eik_config_type), intent(in), optional :: theta_grid_eik_config_in
    integer, parameter :: verb=3
    if (initialized) return
    initialized = .true.


    call debug_message(verb, "init_theta_grid: call read_parameters")
    call read_parameters(theta_grid_config_in)
    call debug_message(verb, "init_theta_grid: call get_sizes")
    call get_sizes(theta_grid_gridgen_config_in, theta_grid_salpha_config_in, theta_grid_file_config_in, theta_grid_eik_config_in)
    if (proc0) then
       call debug_message(verb, "init_theta_grid: call allocate_arrays")
       call allocate_arrays
       call debug_message(verb, "init_theta_grid: call get_grids")
       call get_grids
       call debug_message(verb, "init_theta_grid: call finish_init")
       call finish_init
    end if
    call broadcast_results
  end subroutine init_theta_grid

  !> FIXME : Add documentation
  subroutine finish_theta_grid
    use theta_grid_gridgen, only: finish_theta_grid_gridgen
    use theta_grid_salpha, only: finish_theta_grid_salpha
    use theta_grid_eik, only: finish_theta_grid_eik
    use theta_grid_file, only: finish_theta_grid_file
    use theta_grid_params, only: finish_theta_grid_params
    implicit none

    initialized = .false.

    call finish_theta_grid_gridgen
    call finish_theta_grid_salpha
    call finish_theta_grid_eik
    call finish_theta_grid_file
    ! This is now handled separately by gs2_init
    !call finish_theta_grid_params
    call deallocate_arrays
    call theta_grid_config%reset()
  end subroutine finish_theta_grid

  !> FIXME : Add documentation
  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb = 3

    !Scalars
    call debug_message(verb, "theta_grid::broadcast_results scalars")
    call broadcast (bmin)
    call broadcast (bmax)
    call broadcast (eps_trapped)
    call broadcast (kxfac)
    call broadcast (rhoc)
    call broadcast (qval)
    call broadcast (ntheta)
    call broadcast (ntgrid)
    call broadcast (nperiod)
    call broadcast (nbset)
    call broadcast (shat)
    call broadcast (drhodpsi)
    call broadcast (gb_to_cv)

    !Arrays
    call debug_message(verb, "theta_grid::broadcast_results allocate")
    if (.not. proc0) then
       call allocate_arrays
       if (.not. allocated(theta2)) allocate (theta2(-ntgrid:ntgrid))
       if (.not. allocated(delthet)) allocate (delthet(-ntgrid:ntgrid))
       if (.not. allocated(delthet2)) allocate (delthet2(-ntgrid:ntgrid))
    end if

    call debug_message(verb, "theta_grid::broadcast_results arrays")
    call broadcast (theta)
    call broadcast (theta2)
    call broadcast (delthet)
    call broadcast (delthet2)
    call broadcast (bset)
    call broadcast (bmag)
    call broadcast (itor_over_B)
    call broadcast (IoB)
    call broadcast (gradpar)
    call broadcast (gbdrift)
    call broadcast (gbdrift0)
    call broadcast (cvdrift)
    call broadcast (cvdrift0)
    call broadcast (cdrift)
    call broadcast (cdrift0)
    call broadcast (gds2)
    call broadcast (gds21)
    call broadcast (gds22)
    call broadcast (gds23)
    call broadcast (gds24)
    call broadcast (gds24_noq)
    call broadcast (grho)
    call broadcast (jacob)
    call broadcast (field_line_average_weight)
    call broadcast (Rplot)
    call broadcast (Zplot)
    call broadcast (aplot)
    call broadcast (Rprime)
    call broadcast (Zprime)
    call broadcast (aprime)
    call broadcast (Bpol)
    call debug_message(verb, "theta_grid::broadcast_results finished")
  end subroutine broadcast_results

  !> FIXME : Add documentation
  subroutine read_parameters(theta_grid_config_in)
    use file_utils, only: input_unit, error_unit
    use text_options, only: text_option, get_option_value
    implicit none
    type(theta_grid_config_type), intent(in), optional :: theta_grid_config_in
    type (text_option), dimension (6), parameter :: eqopts = &
         (/ text_option('default', eqopt_eik), &
            text_option('eik', eqopt_eik), &
            text_option('s-alpha', eqopt_salpha), &
            text_option('grid.out', eqopt_file), &
            text_option('file', eqopt_file), &
            text_option('ncfile', eqopt_file_nc) /)
    character(20) :: equilibrium_option
    ! 'default' 'eik': call eikcoefs for parameterized equilibrium
    ! 's-alpha': s-alpha
    ! 'grid.out' 'file': read grid from grid.out file generated by rungridgen

    integer :: ierr

    if (present(theta_grid_config_in)) theta_grid_config = theta_grid_config_in

    call theta_grid_config%init(name = 'theta_grid_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    cvdriftknob = theta_grid_config%cvdriftknob
    equilibrium_option = theta_grid_config%equilibrium_option
    gb_to_cv = theta_grid_config%gb_to_cv
    gbdriftknob = theta_grid_config%gbdriftknob

    exist = theta_grid_config%exist

    ierr = error_unit()
    call get_option_value &
         (equilibrium_option, eqopts, eqopt_switch, &
         ierr, "equilibrium_option in theta_grid_knobs",.true.)
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine allocate_arrays
    implicit none
    if (.not. allocated(theta)) allocate (theta(-ntgrid:ntgrid))
    if (.not. allocated(bset)) allocate (bset(nbset))
    if (.not. allocated(bmag)) allocate (bmag(-ntgrid:ntgrid))
    if (.not. allocated(gradpar)) allocate (gradpar(-ntgrid:ntgrid))
    if (.not. allocated(itor_over_B)) allocate (itor_over_B(-ntgrid:ntgrid))
    if (.not. allocated(IoB)) allocate (IoB(-ntgrid:ntgrid))
    if (.not. allocated(gbdrift)) allocate (gbdrift(-ntgrid:ntgrid))
    if (.not. allocated(gbdrift0)) allocate (gbdrift0(-ntgrid:ntgrid))
    if (.not. allocated(cvdrift)) allocate (cvdrift(-ntgrid:ntgrid))
    if (.not. allocated(cvdrift0)) allocate (cvdrift0(-ntgrid:ntgrid))
    if (.not. allocated(cdrift)) allocate (cdrift(-ntgrid:ntgrid))
    if (.not. allocated(cdrift0)) allocate (cdrift0(-ntgrid:ntgrid))
    if (.not. allocated(gds2)) allocate (gds2(-ntgrid:ntgrid))
    if (.not. allocated(gds21)) allocate (gds21(-ntgrid:ntgrid))
    if (.not. allocated(gds22)) allocate (gds22(-ntgrid:ntgrid))
    if (.not. allocated(gds23)) allocate (gds23(-ntgrid:ntgrid))
    if (.not. allocated(gds24)) allocate (gds24(-ntgrid:ntgrid))
    if (.not. allocated(gds24_noq)) allocate (gds24_noq(-ntgrid:ntgrid))
    if (.not. allocated(grho)) allocate (grho(-ntgrid:ntgrid))
    if (.not. allocated(jacob)) allocate (jacob(-ntgrid:ntgrid))
    if (.not. allocated(Rplot)) allocate (Rplot(-ntgrid:ntgrid))
    if (.not. allocated(Rprime)) allocate (Rprime(-ntgrid:ntgrid))
    if (.not. allocated(Zplot)) allocate (Zplot(-ntgrid:ntgrid))
    if (.not. allocated(Zprime)) allocate (Zprime(-ntgrid:ntgrid))
    if (.not. allocated(aplot)) allocate (aplot(-ntgrid:ntgrid))
    if (.not. allocated(aprime)) allocate (aprime(-ntgrid:ntgrid))
    if (.not. allocated(Bpol)) allocate (Bpol(-ntgrid:ntgrid))
  end subroutine allocate_arrays

  !> FIXME : Add documentation
  subroutine deallocate_arrays
    implicit none
    if(allocated(theta)) deallocate (theta)
    if(allocated(bset)) deallocate (bset)
    if(allocated(bmag)) deallocate (bmag)
    if(allocated(gradpar)) deallocate (gradpar)
    if(allocated(itor_over_B)) deallocate (itor_over_B)
    if(allocated(IoB)) deallocate (IoB)
    if(allocated(gbdrift)) deallocate (gbdrift)
    if(allocated(gbdrift0)) deallocate (gbdrift0)
    if(allocated(cvdrift)) deallocate (cvdrift)
    if(allocated(cvdrift0)) deallocate (cvdrift0)
    if(allocated(cdrift)) deallocate (cdrift)
    if(allocated(cdrift0)) deallocate (cdrift0)
    if(allocated(gds2)) deallocate (gds2)
    if(allocated(gds21)) deallocate (gds21)
    if(allocated(gds22)) deallocate (gds22)
    if(allocated(gds23)) deallocate (gds23)
    if(allocated(gds24)) deallocate (gds24)
    if(allocated(gds24_noq)) deallocate (gds24_noq)
    if(allocated(grho)) deallocate (grho)
    if(allocated(jacob)) deallocate (jacob)
    if(allocated(Rplot)) deallocate (Rplot)
    if(allocated(Rprime)) deallocate (Rprime)
    if(allocated(Zplot)) deallocate (Zplot)
    if(allocated(Zprime)) deallocate (Zprime)
    if(allocated(aplot)) deallocate (aplot)
    if(allocated(aprime)) deallocate (aprime)
    if(allocated(Bpol)) deallocate (Bpol)

    if(allocated(theta2)) deallocate (theta2)
    if(allocated(delthet)) deallocate (delthet)
    if(allocated(delthet2)) deallocate (delthet2)

  end subroutine deallocate_arrays

  !> FIXME : Add documentation
  subroutine finish_init
    use integration, only: trapezoidal_integration
    implicit none
    real, dimension (nbset) :: bset_save
    real, dimension (-ntgrid:ntgrid) :: eik_save

    ! in case nbset changes after gridgen
    if (nbset /= size(bset)) then
       bset_save = bset(:nbset)
       deallocate (bset)
       allocate (bset(nbset))
       bset = bset_save
    end if

    ! in case ntgrid changes after gridgen
    if (ntgrid*2+1 /= size(theta)) then

       eik_save = theta(-ntgrid:ntgrid); deallocate (theta)
       allocate (theta(-ntgrid:ntgrid)); theta = eik_save

       eik_save = bmag(-ntgrid:ntgrid); deallocate (bmag)
       allocate (bmag(-ntgrid:ntgrid)); bmag = eik_save

       eik_save = gradpar(-ntgrid:ntgrid); deallocate (gradpar)
       allocate (gradpar(-ntgrid:ntgrid)); gradpar = eik_save

       eik_save = itor_over_B(-ntgrid:ntgrid); deallocate (itor_over_B)
       allocate (itor_over_B(-ntgrid:ntgrid)); itor_over_B = eik_save

       eik_save = IoB(-ntgrid:ntgrid); deallocate (IoB)
       allocate (IoB(-ntgrid:ntgrid)); IoB = eik_save

       eik_save = gbdrift(-ntgrid:ntgrid); deallocate (gbdrift)
       allocate (gbdrift(-ntgrid:ntgrid)); gbdrift = eik_save

       eik_save = gbdrift0(-ntgrid:ntgrid); deallocate (gbdrift0)
       allocate (gbdrift0(-ntgrid:ntgrid)); gbdrift0 = eik_save

       eik_save = cvdrift(-ntgrid:ntgrid); deallocate (cvdrift)
       allocate (cvdrift(-ntgrid:ntgrid)); cvdrift = eik_save

       eik_save = cvdrift0(-ntgrid:ntgrid); deallocate (cvdrift0)
       allocate (cvdrift0(-ntgrid:ntgrid)); cvdrift0 = eik_save

       eik_save = cdrift(-ntgrid:ntgrid); deallocate (cdrift)
       allocate (cdrift(-ntgrid:ntgrid)); cdrift = eik_save

       eik_save = cdrift0(-ntgrid:ntgrid); deallocate (cdrift0)
       allocate (cdrift0(-ntgrid:ntgrid)); cdrift0 = eik_save

       eik_save = gds2(-ntgrid:ntgrid); deallocate (gds2)
       allocate (gds2(-ntgrid:ntgrid)); gds2 = eik_save

       eik_save = gds21(-ntgrid:ntgrid); deallocate (gds21)
       allocate (gds21(-ntgrid:ntgrid)); gds21 = eik_save

       eik_save = gds22(-ntgrid:ntgrid); deallocate (gds22)
       allocate (gds22(-ntgrid:ntgrid)); gds22 = eik_save

       eik_save = gds23(-ntgrid:ntgrid); deallocate (gds23)
       allocate (gds23(-ntgrid:ntgrid)); gds23 = eik_save

       eik_save = gds24(-ntgrid:ntgrid); deallocate (gds24)
       allocate (gds24(-ntgrid:ntgrid)); gds24 = eik_save

       eik_save = gds24_noq(-ntgrid:ntgrid); deallocate (gds24_noq)
       allocate (gds24_noq(-ntgrid:ntgrid)); gds24_noq = eik_save

       eik_save = grho(-ntgrid:ntgrid); deallocate (grho)
       allocate (grho(-ntgrid:ntgrid)); grho = eik_save

       eik_save = jacob(-ntgrid:ntgrid); deallocate (jacob)
       allocate (jacob(-ntgrid:ntgrid)); jacob = eik_save

       eik_save = Rplot(-ntgrid:ntgrid); deallocate (Rplot)
       allocate (Rplot(-ntgrid:ntgrid)); Rplot = eik_save

       eik_save = Zplot(-ntgrid:ntgrid); deallocate (Zplot)
       allocate (Zplot(-ntgrid:ntgrid)); Zplot = eik_save

       eik_save = aplot(-ntgrid:ntgrid); deallocate (aplot)
       allocate (aplot(-ntgrid:ntgrid)); aplot = eik_save

       eik_save = Rprime(-ntgrid:ntgrid); deallocate (Rprime)
       allocate (Rprime(-ntgrid:ntgrid)); Rprime = eik_save

       eik_save = Zprime(-ntgrid:ntgrid); deallocate (Zprime)
       allocate (Zprime(-ntgrid:ntgrid)); Zprime = eik_save

       eik_save = aprime(-ntgrid:ntgrid); deallocate (aprime)
       allocate (aprime(-ntgrid:ntgrid)); aprime = eik_save

       eik_save = Bpol(-ntgrid:ntgrid); deallocate (Bpol)
       allocate (Bpol(-ntgrid:ntgrid)); Bpol = eik_save
    end if

    bmax = maxval(bmag)
    bmin = minval(bmag)
    ! ?? check Krook collision operator coding which is only place eps is used
    ! the line with bmin/bmax was the original coding.  Changed in 2002-2004 time
    ! frame to 1.0 - 1.0/bmax , now changed back (8.19.04) BD
    ! Now only used in le_grids to check if we have any trapped particles
    eps_trapped = 1.0 - sqrt(bmin/bmax)

    if (.not. allocated(theta2)) allocate (theta2(-ntgrid:ntgrid))
    if (.not. allocated(delthet)) allocate (delthet(-ntgrid:ntgrid))
    if (.not. allocated(delthet2)) allocate (delthet2(-ntgrid:ntgrid))

    theta2 = theta*theta
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)
    delthet(ntgrid) = 0.!delthet(-ntgrid)
    delthet2 = delthet*delthet

    ! Note here we redefine jacob despite the handling in the above if block
    ! suggesting this has been defined.
    jacob = 1.0/(drhodpsi*gradpar*bmag)

    ! Calculate the weight used in the field line average routine to save
    ! repeated calculation of this constant
    field_line_average_weight = 1.0/trapezoidal_integration(theta, jacob)
  end subroutine finish_init

  !> FIXME : Add documentation
  subroutine get_sizes(theta_grid_gridgen_config_in, theta_grid_salpha_config_in, theta_grid_file_config_in, theta_grid_eik_config_in)
    use theta_grid_eik, only: eik_get_sizes, init_theta_grid_eik, theta_grid_eik_config_type
    use theta_grid_salpha, only: salpha_get_sizes, init_theta_grid_salpha, theta_grid_salpha_config_type
    use theta_grid_file, only: file_get_sizes, file_nc_get_sizes, init_theta_grid_file, theta_grid_file_config_type
    use theta_grid_file, only: ntheta_file=>ntheta, nperiod_file=>nperiod
    use theta_grid_file, only: nbset_file=>nbset
    use theta_grid_gridgen, only: theta_grid_gridgen_init, theta_grid_gridgen_config_type
    implicit none
    type(theta_grid_gridgen_config_type), intent(in), optional :: theta_grid_gridgen_config_in
    type(theta_grid_salpha_config_type), intent(in), optional :: theta_grid_salpha_config_in
    type(theta_grid_file_config_type), intent(in), optional :: theta_grid_file_config_in
    type(theta_grid_eik_config_type), intent(in), optional :: theta_grid_eik_config_in

    logical, parameter :: debug=.false.
if (debug) write(6,*) 'get_sizes: eqopt_switch=',eqopt_switch
    select case (eqopt_switch)
    case (eqopt_eik)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_eik'
       call init_theta_grid_eik(theta_grid_eik_config_in)
if (debug) write(6,*) 'get_sizes: call eik_get_sizes'
       call eik_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_salpha)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_salpha'
       call init_theta_grid_salpha(theta_grid_salpha_config_in)
if (debug) write(6,*) 'get_sizes: call salpha_get_sizes'
       call salpha_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_file)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_file'
       call init_theta_grid_file(theta_grid_file_config_in)
if (debug) write(6,*) 'get_sizes: call file_get_sizes'
       call file_get_sizes
       ntheta=ntheta_file
       nperiod=nperiod_file
       nbset=nbset_file
    case (eqopt_file_nc)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_file'
       call init_theta_grid_file(theta_grid_file_config_in)
if (debug) write(6,*) 'get_sizes: call file_nc_get_sizes'
       call file_nc_get_sizes
       ntheta=ntheta_file
       nperiod=nperiod_file
       nbset=nbset_file
    end select
    ntgrid = ntheta/2 + (nperiod-1)*ntheta

    ! Make sure gridgen is setup on all procs
    call theta_grid_gridgen_init(theta_grid_gridgen_config_in)

if (debug) write(6,*) 'get_sizes: done'
  end subroutine get_sizes

  !> FIXME : Add documentation
  subroutine get_grids
    use theta_grid_eik, only: eik_get_grids
    use theta_grid_salpha, only: salpha_get_grids
    use theta_grid_file, only: file_get_grids
    use theta_grid_file, only: file_nc_get_grids
    use theta_grid_params, only: eps, btor_slab
    implicit none
    logical, parameter :: debug=.false.
    select case (eqopt_switch)
    case (eqopt_eik)
if (debug) write(6,*) 'get_grids: call eik_get_grids'
       call eik_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
       shape = 'torus   '
    case (eqopt_salpha)
if (debug) write(6,*) 'get_grids: call salpha_get_grids'
       call salpha_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, shape, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
    case (eqopt_file)
if (debug) write(6,*) 'get_grids: call file_get_grids'
       call file_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
       shape = 'torus   '
    case (eqopt_file_nc)
if (debug) write(6,*) 'get_grids: call file_nc_get_grids'
       call file_nc_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
       shape = 'torus   '
    end select
    kxfac = abs(kxfac)
    qval = abs(qval)
!CMR, 4/6/2014
!   knobs to independently control magnetic gradB and curvature drifts
    gbdrift=gbdriftknob*gbdrift
    gbdrift0=gbdriftknob*gbdrift0
    cvdrift=cvdriftknob*cvdrift
    cvdrift0=cvdriftknob*cvdrift0

    itor_over_B=0.
!CMR, 2/2/2011:
! If using slab geometry, set itor_over_B = btor_slab from "theta_grid_params":
! cleaner equivalent alternative to using btor_slab in "dist_fn_knobs", and
! sets geometric parameter itor_over_B in one place for ALL geometries.
!
    if (eqopt_switch .eq. eqopt_salpha .and. eps .eq. 0. ) then
       itor_over_B = btor_slab
    else
!CMR, 19/10/10: moved MAB's definition of geometry quantity itor_over_B from
!               dist_fn.f90 to here.
! Calculate the parallel velocity shear drive factor itor_over_B
! (which effectively depends on the angle the field lines make with the flow)
! note that the following is only valid in a torus!
! itor_over_B = (q/rho) * Rmaj*Btor/(a*B)
       IoB = sqrt(max(Rplot**2 - (grho/(bmag*drhodpsi))**2,0.))
    ! RN> 2011/1/25: fix here avoids dividing by rhoc if rhoc=0
    ! CMR, 2/2/2011: set itor_over_B=0 if rhoc=0
    !                Dropping parallel sheared flow source term in GKE
    !                using itor_over_B=0 is safer than itor_over_B=NaN!
       if (rhoc /= 0.) then
          itor_over_B = qval / rhoc * IoB
       else
          itor_over_B = btor_slab
       end if
     endif
  end subroutine get_grids

  !> Calculates the field line / theta average of a passed quantity
  real function field_line_average_real(quantity) result(integral)
    use integration, only: trapezoidal_integration
    implicit none
    real, dimension(:), intent(in) :: quantity
    integral = trapezoidal_integration(theta, quantity*jacob) * field_line_average_weight
  end function field_line_average_real

  !> Calculates the field line / theta average of a passed quantity
  complex function field_line_average_complex(quantity) result(integral)
    use integration, only: trapezoidal_integration
    implicit none
    complex, dimension(:), intent(in) :: quantity
    integral = trapezoidal_integration(theta, quantity*jacob) * field_line_average_weight
  end function field_line_average_complex

  !> Helper function to determine if the shear is small enough that we consider
  !> it to be zero for the purposes of periodicity etc.
  pure logical function is_effectively_zero_shear(shear_in) result(zero_shear)
    use optionals, only: get_option_with_default
    implicit none
    real, intent(in), optional :: shear_in
    real :: the_shear
    the_shear = get_option_with_default(shear_in, shat)
    zero_shear = abs(the_shear) <= smallest_non_zero_shear
  end function is_effectively_zero_shear

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_theta_grid_config(theta_grid_config_in)
    use mp, only: mp_abort
    type(theta_grid_config_type), intent(in), optional :: theta_grid_config_in
    if (initialized) then
       call mp_abort("Trying to set theta_grid_config when already initialized.", to_screen = .true.)
    end if
    if (present(theta_grid_config_in)) then
       theta_grid_config = theta_grid_config_in
    end if
  end subroutine set_theta_grid_config

#include "theta_grid_auto_gen.inc"
end module theta_grid
