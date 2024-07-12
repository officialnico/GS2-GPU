!> Set up values of kx and ky for linear runs that use a single k_perp mode.
module kt_grids_single
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: init_kt_grids_single, single_get_sizes, single_get_grids
  public :: check_kt_grids_single, wnml_kt_grids_single
  public :: read_parameters_single
  public :: finish_parameters_single

  public :: kt_grids_single_config_type
  public :: set_kt_grids_single_config
  public :: get_kt_grids_single_config

  real :: akx, aky, theta0, rhostar_single
  integer :: n0
  logical :: parameters_read = .false.
  logical :: initialized = .false.
  real, parameter :: default_unset_value = -12345.6789
  !> Used to represent the input configuration of kt_grids_single
  type, extends(abstract_config_type) :: kt_grids_single_config_type
     ! namelist : kt_grids_single_parameters
     ! indexed : false
     !> \(k_x \rho\) for the reference species (but recommended to set `theta0` instead).
     real :: akx = default_unset_value
     !> \(k_y \rho\) for the reference species.
     real :: aky = 0.4
     !> if `n0`>0 use toroidal mode number to override `aky` and set
     !> `aky=n0*drhodpsi*rhostar_single` where `drhodpsi` is calculated
     !> as a part of the geometry setup.
     integer :: n0 = 0
     !> Used in conjunction with `n0`:
     !> `aky=n0*drhodpsi*rhostar_single` (if `n0` is set) where
     !> `drhodpsi` is calculated as a part of the geometry setup.
     real :: rhostar_single = 1.0e-4
     !> \(\theta_0\) is the ballooning angle, sets the point in
     !> \(\theta\) where the radial wavenumber is zero
     real :: theta0 = 0.0
   contains
     procedure, public :: read => read_kt_grids_single_config
     procedure, public :: write => write_kt_grids_single_config
     procedure, public :: reset => reset_kt_grids_single_config
     procedure, public :: broadcast => broadcast_kt_grids_single_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_single_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_single_config
  end type kt_grids_single_config_type

  type(kt_grids_single_config_type) :: kt_grids_single_config

contains

  !> FIXME : Add documentation
  subroutine read_parameters_single(kt_grids_single_config_in)
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in
    logical :: exist

    if (parameters_read) return
    parameters_read = .true.

    if (present(kt_grids_single_config_in)) kt_grids_single_config = kt_grids_single_config_in

    call kt_grids_single_config%init(name = 'kt_grids_single_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    akx = kt_grids_single_config%akx
    aky = kt_grids_single_config%aky
    n0 = kt_grids_single_config%n0
    rhostar_single = kt_grids_single_config%rhostar_single
    theta0 = kt_grids_single_config%theta0

    exist = kt_grids_single_config%exist
  end subroutine read_parameters_single

  !> FIXME : Add documentation
  subroutine finish_parameters_single
    implicit none
    parameters_read = .false.
    initialized = .false.
    call kt_grids_single_config%reset()
  end subroutine finish_parameters_single

  !> FIXME : Add documentation
  subroutine init_kt_grids_single(kt_grids_single_config_in)
    !CMR, 14/10/2013:
    ! New namelist variables n0, rhostar_single to set aky using toroidal mode number.
    ! Toroidal modenumber used if n0> 0 prescribed in input file.
    use theta_grid, only: drhodpsi, shat
    implicit none
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in

    if(initialized) return
    initialized = .true.

    call read_parameters_single(kt_grids_single_config_in)

    if (n0 .gt. 0) then
       !CMR if n0>0 then override aky inputs and use n0 to determine aky
       aky=n0*drhodpsi*rhostar_single
    endif

    ! Try to ensure consistency between akx and theta0
    if (akx == default_unset_value) then
       ! If kx hasn't been set, set it from theta0
       akx = theta0 * aky * shat
    else
       ! If kx has been set, force theta0 to be calculated from this
       if (shat /= 0 .and. aky /= 0 ) then
          theta0 = akx / (aky * shat)
       end if
    end if

  end subroutine init_kt_grids_single

  !> FIXME : Add documentation
  subroutine wnml_kt_grids_single(unit)
    implicit none
    integer, intent(in) :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_single_parameters"
    write (unit, fmt="(' aky = ',e17.10)") aky
    write (unit, fmt="(' theta0 = ',e17.10)") theta0
    write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_single

  !> FIXME : Add documentation
  subroutine single_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny

    naky = 1  ;  ntheta0 = 1
    nx = 0    ;  ny = 0

  end subroutine single_get_sizes

  !> FIXME : Add documentation
  subroutine single_get_grids (aky_out, theta0_out, akx_out, ikx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out, akx_out
    real, dimension (:,:), intent (out) :: theta0_out
    !> Discrete kx wavenumber grid indices
    integer, dimension (:), intent (out) :: ikx_out

    aky_out = aky
    theta0_out = theta0
    akx_out = akx
    ikx_out = 0
  end subroutine single_get_grids

  !> FIXME : Add documentation
  subroutine check_kt_grids_single(report_unit)
    implicit none
    integer, intent(in) :: report_unit

    write (report_unit, *)
    write (report_unit, fmt="('A single k_perp will be evolved, with: ')")
    if (n0 .gt.0) write (report_unit, fmt="('ky set using toroidal mode number, n0=',i8/T24,'rhostar_single=',1pe12.4)") n0, rhostar_single
    write (report_unit, *)
    write (report_unit, fmt="('ky rho = ',f10.4)") aky
    write (report_unit, fmt="('theta_0 = ',f10.4)") theta0
    if (akx /= 0.) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The value of akx in the kt_grids_single_parameters namelist is ignored.')")
       write (report_unit, fmt="('You have set akx to a non-zero value.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    end if
  end subroutine check_kt_grids_single

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_single_config(kt_grids_single_config_in)
    use mp, only: mp_abort
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_single_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_single_config_in)) then
       kt_grids_single_config = kt_grids_single_config_in
    end if
  end subroutine set_kt_grids_single_config

#include "kt_grids_single_auto_gen.inc"  
end module kt_grids_single
