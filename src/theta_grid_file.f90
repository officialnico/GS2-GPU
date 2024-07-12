#include "unused_dummy.inc"
!> Use output from [[rungridgen]]
module theta_grid_file
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use constants, only: run_name_size
  implicit none

  private

  public :: init_theta_grid_file, finish_theta_grid_file
  public :: check_theta_grid_file, wnml_theta_grid_file
  public :: file_get_sizes, file_get_grids
  public :: check_theta_grid_file_nc
  public :: file_nc_get_sizes, file_nc_get_grids
  public :: ntheta, nperiod, ntgrid, nbset

  public :: theta_grid_file_config_type
  public :: set_theta_grid_file_config
  public :: get_theta_grid_file_config

  character (len = run_name_size) :: gridout_file
  real :: shat_input, drhodpsi_input, kxfac_input, qval_input
  logical :: no_geo_info = .false.
  integer :: ntheta, nperiod, ntgrid, nbset
  logical :: exist, initialized = .false.

  !> Used to represent the input configuration of theta_grid
  type, extends(abstract_config_type) :: theta_grid_file_config_type
     ! namelist : theta_grid_file_knobs
     ! indexed : false
     !> Name of file with output from [[rungridgen]].
     character(len = run_name_size) :: gridout_file = "grid.out"
     !> If false, read `Rplot`, `Rprime`, `Zplot`, `Zprime`, `aplot`, `aprime`
     !> from [[theta_grid_file_knobs:gridout_file]].
     logical :: no_geo_info = .false.
   contains
     procedure, public :: read => read_theta_grid_file_config
     procedure, public :: write => write_theta_grid_file_config
     procedure, public :: reset => reset_theta_grid_file_config
     procedure, public :: broadcast => broadcast_theta_grid_file_config
     procedure, public, nopass :: get_default_name => get_default_name_theta_grid_file_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_theta_grid_file_config
  end type theta_grid_file_config_type

  type(theta_grid_file_config_type) :: theta_grid_file_config

contains

  !> FIXME : Add documentation
  subroutine wnml_theta_grid_file(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_file_knobs"
    write (unit, fmt="(' gridout_file = ',a)") '"'//trim(gridout_file)//'"'
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_file

  !> FIXME : Add documentation
  subroutine check_theta_grid_file(report_unit)
    use file_utils, only: get_unused_unit
    implicit none
    integer, intent(in) :: report_unit
    integer :: i, iunit
    real :: drhodpsi, kxfac, rmaj, shat
    character (200) :: line

    write (report_unit, *)
    write (report_unit, fmt="('Equilibrium information obtained from gridgen output file:')")
    write (report_unit, fmt="(a)") trim(gridout_file)

    call get_unused_unit (iunit)
    open (unit=iunit, file=gridout_file, status="old", err=100)
    read (unit=iunit, fmt="(a)") line
    read (unit=iunit, fmt=*) nbset
    read (unit=iunit, fmt="(a)") line
    do i = 1, nbset
       read (unit=iunit, fmt="(a)") line
    end do

    read (unit=iunit, fmt="(a)") line
    read (unit=iunit, fmt=*) ntgrid, nperiod, ntheta, &
         drhodpsi, rmaj, shat, kxfac

    close (unit=iunit)

    write (report_unit, *)
    write (report_unit, fmt="('Limited information available:')")
    write (report_unit, *)
    write (report_unit, fmt="('nbset =     ',i5)") nbset
    write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
    write (report_unit, fmt="('ntheta =    ',i5)") ntheta
    write (report_unit, fmt="('nperiod =   ',i2)") nperiod
    write (report_unit, fmt="('drhodpsi =  ',f8.4)") drhodpsi
    write (report_unit, fmt="('R =         ',f8.4)") Rmaj
    write (report_unit, fmt="('s_hat =     ',f8.4)") shat
    write (report_unit, fmt="('kxfac =     ',f8.4)") kxfac

    write (report_unit, *)
    write (report_unit, *) 'NOTE: Regardless of the values of ntheta and nperiod'
    write (report_unit, *) '      found in the theta_grid_parameters namelist,'
    write (report_unit, *) '      this calculation will use the values listed here:'
    write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
    write (report_unit, fmt="('ntheta =    ',i5)") ntheta
    write (report_unit, *) '      These were obtained from the gridgen output file.'
    write (report_unit, *)
100 continue
  end subroutine check_theta_grid_file

  !> FIXME : Add documentation
  subroutine check_theta_grid_file_nc(report_unit)
    use gs2_io_grids, only: nc_grid_file_open, nc_grid_file_close
    use gs2_io_grids, only: nc_get_grid_sizes, nc_get_grid_scalars
    use gs2_io_grids, only: nc_get_lambda_grid_size
    implicit none
    integer, intent(in) :: report_unit
    real :: drhodpsi, kxfac, rmaj, shat, qval

    write (report_unit, *)
    write (report_unit, fmt="('Equilibrium information obtained from gridgen output file:')")
    write (report_unit, fmt="(a)") trim(gridout_file)

    call nc_grid_file_open(gridout_file, "r")
    call nc_get_grid_sizes(ntheta=ntheta, ntgrid=ntgrid, nperiod=nperiod)
    call nc_get_lambda_grid_size(nbset)
    call nc_get_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj)
    call nc_grid_file_close()

    write (report_unit, *)
    write (report_unit, fmt="('Limited information available:')")
    write (report_unit, *)
    write (report_unit, fmt="('nbset =     ',i5)") nbset
    write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
    write (report_unit, fmt="('ntheta =    ',i5)") ntheta
    write (report_unit, fmt="('nperiod =   ',i2)") nperiod
    write (report_unit, fmt="('drhodpsi =  ',f8.4)") drhodpsi
    write (report_unit, fmt="('R =         ',f8.4)") Rmaj
    write (report_unit, fmt="('s_hat =     ',f8.4)") shat
    write (report_unit, fmt="('kxfac =     ',f8.4)") kxfac
    write (report_unit, fmt="('qval =     ',f8.4)") qval

    write (report_unit, *)
    write (report_unit, *) 'NOTE: Regardless of the values of ntheta and nperiod'
    write (report_unit, *) '      found in the theta_grid_parameters namelist,'
    write (report_unit, *) '      this calculation will use the values listed here:'
    write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
    write (report_unit, fmt="('ntheta =    ',i5)") ntheta
    write (report_unit, *) '      These were obtained from the gridgen output file.'
    write (report_unit, *)
  end subroutine check_theta_grid_file_nc

  !> FIXME : Add documentation
  subroutine init_theta_grid_file(theta_grid_file_config_in)
    use theta_grid_params, only: init_theta_grid_params
    implicit none
    type(theta_grid_file_config_type), intent(in), optional :: theta_grid_file_config_in

    if (initialized) return
    initialized = .true.

    call init_theta_grid_params
    call read_parameters(theta_grid_file_config_in)
  end subroutine init_theta_grid_file

  !> FIXME : Add documentation
  subroutine finish_theta_grid_file
    implicit none
    initialized = .false.
    call theta_grid_file_config%reset()
  end subroutine finish_theta_grid_file

  !> FIXME : Add documentation
  subroutine read_parameters(theta_grid_file_config_in)
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    type(theta_grid_file_config_type), intent(in), optional :: theta_grid_file_config_in

    if (present(theta_grid_file_config_in)) theta_grid_file_config = theta_grid_file_config_in

    call theta_grid_file_config%init(name = 'theta_grid_file_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    gridout_file = theta_grid_file_config%gridout_file
    no_geo_info = theta_grid_file_config%no_geo_info

    exist = theta_grid_file_config%exist
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine file_get_sizes
    use file_utils, only: get_unused_unit
    implicit none
    integer :: unit
    character(200) :: line
    integer :: i, ntgrid
    real :: rmaj

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) nbset
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt="(a)") line
    end do

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) ntgrid, nperiod, ntheta, &
         drhodpsi_input, rmaj, shat_input, kxfac_input, qval_input

    close (unit=unit)
  end subroutine file_get_sizes

  !> FIXME : Add documentation
  subroutine file_nc_get_sizes
    use gs2_io_grids, only: nc_grid_file_open
    use gs2_io_grids, only: nc_get_grid_sizes, nc_get_lambda_grid_size
    implicit none
    call nc_grid_file_open(gridout_file, "r")
    call nc_get_grid_sizes(ntheta=ntheta, ntgrid=ntgrid, nperiod=nperiod)
    call nc_get_lambda_grid_size(nbset)
  end subroutine file_nc_get_sizes

  !> FIXME : Add documentation
  subroutine file_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
    use file_utils, only: get_unused_unit
    use constants, only: pi
    use integration, only: trapezoidal_integration
    use theta_grid_params, only: rhoc_par => rhoc
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
    logical, intent (in) :: gb_to_cv
    integer :: unit
    character(200) :: line
    integer :: i
    logical, save :: first = .true.
    !Not used because we don't call gridgen_get_grids for file
    UNUSED_DUMMY(nperiod); UNUSED_DUMMY(ntheta)

    !<DD> Should jacob also be provided by this routine?

    !<DD> NOTE: Not currently settin Bpol here. This is used in flux calculations.
    !If not set then results will be funny. Add a warning message and set to zero
    !for now.
    if (first) then
       write(*,'("WARNING: When using file_get_grids, Bpol does not have a correct definition --> Currently just setting to zero, may impact some diagnostics")')
       Bpol = 0.
       first = .false.
    end if

    shat = shat_input
    drhodpsi = drhodpsi_input
    kxfac = kxfac_input
    qval = qval_input

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt=*) bset(i) ! actually alambda
    end do
    bset = 1.0/bset ! switch alambda to bset

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gbdrift(i), gradpar(i), grho(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gds21(i), gds22(i)
    end do

    ! TMP UNTIL WORK OUT HOW TO GET FROM FILE
    gds23 = 0. ; gds24 = 0. ; gds24_noq = 0.

    ! TMP UNTIL FIGURE OUT HOW TO WORK WITH FILE -- MAB
    ! set coriolis drift to zero
    cdrift = 0. ; cdrift0 = 0.

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
    end do

    if (gb_to_cv) then
       do i =-ntgrid,ntgrid
          gbdrift(i) = cvdrift(i)
          gbdrift0(i) = cvdrift0(i)
       end do
    end if

    ! Possibly crude approximations for file
    ! Note jacob is usually defined as jacob = 1.0/(drhodpsi*gradpar*bmag)
    ! but isn't yet stored here.
    surfarea = 2 * pi * trapezoidal_integration(theta, grho / (drhodpsi*gradpar*bmag))
    dvdrhon = 2 * pi * trapezoidal_integration(theta, 1.0 / (drhodpsi*gradpar*bmag))

    ! As the eqfile doesn't specify rhoc we simply set from the value
    ! from theta_grid_params here.
    rhoc = rhoc_par

    if (.not. no_geo_info) then

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) Rplot(i), Rprime(i)
       end do

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) Zplot(i), Zprime(i)
       end do

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) aplot(i), aprime(i)
       end do

    end if

    close (unit=unit)

    return

100 continue
    write(6,*) 'Error reading Rplot etc. setting to dummy values.'
! dummy values for backward compatibility
    Rplot = 1. ; Rprime = 0.
    Zplot = 1. ; Zprime = 0.
    aplot = 1. ; aprime = 0.

    close (unit=unit)

  end subroutine file_get_grids

  !> FIXME : Add documentation
  subroutine file_nc_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
    use constants, only: pi
    use integration, only: trapezoidal_integration
    use theta_grid_params, only: rhoc_par => rhoc
    use gs2_io_grids, only: nc_grid_file_close
    use gs2_io_grids, only: nc_get_grids, nc_get_grid_scalars
    use gs2_io_grids, only: nc_get_lambda_grid
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
    real :: rmaj
    logical, intent (in) :: gb_to_cv
    !Not used because we don't call gridgen_get_grids for file
    UNUSED_DUMMY(nperiod); UNUSED_DUMMY(ntheta)
    Rplot = 1.; Rprime = 0.
    Zplot = 1.; Zprime = 0.
    aplot = 1.; aprime = 0.

    call nc_get_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj)

    call nc_get_grids(ntgrid, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, theta, &
       cdrift=cdrift, cdrift0=cdrift0, &
       Rplot=Rplot, Rprime=Rprime, Zplot=Zplot, Zprime=Zprime, aplot=aplot, aprime=aprime, &
       Bpol=Bpol)

    call nc_get_lambda_grid(nbset, bset)

    call nc_grid_file_close()

    ! switch alambda to bset
    bset = 1.0/bset

    ! TMP UNTIL WORK OUT HOW TO GET FROM FILE
    gds23 = 0. ; gds24 = 0. ; gds24_noq = 0.

    ! TMP UNTIL FIGURE OUT HOW TO WORK WITH FILE -- MAB
    ! set coriolis drift to zero
    cdrift = 0. ; cdrift0 = 0.

    if (gb_to_cv) then
       gbdrift = cvdrift
       gbdrift0 = cvdrift0
    end if

    ! Possibly crude approximations for file
    ! Note jacob is usually defined as jacob = 1.0/(drhodpsi*gradpar*bmag)
    ! but isn't yet stored here.
    surfarea = 2 * pi * trapezoidal_integration(theta, grho / (drhodpsi*gradpar*bmag))
    dvdrhon = 2 * pi * trapezoidal_integration(theta, 1.0 / (drhodpsi*gradpar*bmag))

    ! As the eqfile doesn't specify rhoc we simply set from the value
    ! from theta_grid_params here.
    rhoc = rhoc_par

  end subroutine file_nc_get_grids

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_theta_grid_file_config(theta_grid_file_config_in)
    use mp, only: mp_abort
    type(theta_grid_file_config_type), intent(in), optional :: theta_grid_file_config_in
    if (initialized) then
       call mp_abort("Trying to set theta_grid_file_config when already initialized.", to_screen = .true.)
    end if
    if (present(theta_grid_file_config_in)) then
       theta_grid_file_config = theta_grid_file_config_in
    end if
  end subroutine set_theta_grid_file_config

#include "theta_grid_file_auto_gen.inc"
end module theta_grid_file
