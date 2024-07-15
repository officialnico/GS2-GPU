!> A module for handling the configuration of the optimisation
!! module via the namelist optimisation_config.
module optimisation_configuration
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use overrides, only: optimisations_overrides_type
  implicit none

  private

  public :: init_optimisation_config
  public :: finish_optimisation_config
  public :: optimisation_type
  public :: optimisation_results_type

  public :: optimisation_config_type
  public :: set_optimisation_config
  public :: get_optimisation_config

  logical :: initialized = .false.
  
  type optimisation_results_type
    ! Configuration

    ! Results
    real :: time
    real :: optimal_time
    real :: cost
    real :: optimal_cost
    real :: efficiency
    integer :: nproc
    logical :: optimal = .true.

  end type optimisation_results_type

  !> A type for storing the optimisation configuration,
  !! the results
  type optimisation_type
     integer :: nproc_max
     type(optimisation_results_type) :: results
     type(optimisations_overrides_type), &
          dimension(:), allocatable :: sorted_optimisations
     type(optimisation_results_type), dimension(:), allocatable :: sorted_results
     real :: timing_rel_error
     real :: timing_max_rel_error
     integer :: outunit
     logical :: on
     logical :: auto
     logical :: measure_all
     logical :: warm_up
     logical :: estimate_timing_error
     integer :: nstep_measure
     real :: max_imbalance
     integer :: max_unused_procs
     real :: min_efficiency
  end type optimisation_type

  !> Used to represent the input configuration for GS2's optimisation
  !> procedure. When turned on, GS2 performs a scan for the given input file,
  !> varying different optimisation flags. Results of this scan are reported
  !> in <runname>.optim. The optimal parameters are stored, allowing a user to
  !> run the same input file with optimised parameters in the same execution.
  !> A user can also choose to continue with a less-than-optimal set of
  !> parameters which satisfy other constraints.
  type, extends(abstract_config_type) :: optimisation_config_type
     ! namelist : optimisation_config
     ! indexed : false     
     !> When true, automatically continues GS2 to run the input file with the
     !> optimised parameters.
     logical :: auto = .true.
     !> Estimate the absolute and relative errors in timing data
     !> FIXME: Why would we want this to be false? On small core counts it
     !> doesn't seem like a big overhead.
     logical :: estimate_timing_error = .true.
     !> The maximum fraction of unused procs to allow in the optimisation
     !> scan.
     real :: max_imbalance = -1.0
     !> The maximum number of unused procs to allow in the optimisation scan.
     integer :: max_unused_procs = 0
     !> When true, use the "advance" timer.
     !> When false, use the "timestep" timer.
     logical :: measure_all = .false.
     !> The minimum efficiency (relative to the optimal parameters)
     !> considered when looking for a constrained set of parameters.
     !> A negative value implies only the optimal parameters are considered.
     real :: min_efficiency = -1.0
     !> The number of timestep to use in the timing experiments. Must
     !> be greater than 1
     integer :: nstep_measure = 5
     !> Set true to turn on optimisation procedure
     logical :: on = .false.
     !> When true, perform a few runs before beginning the timing experiment
     logical :: warm_up = .false.
   contains
     procedure, public :: read => read_optimisation_config
     procedure, public :: write => write_optimisation_config
     procedure, public :: reset => reset_optimisation_config
     procedure, public :: broadcast => broadcast_optimisation_config
     procedure, public, nopass :: get_default_name => get_default_name_optimisation_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_optimisation_config
  end type optimisation_config_type
  
  type(optimisation_config_type) :: optimisation_config
  
contains
  subroutine init_optimisation_config(optim, optimisation_config_in)
    use file_utils, only: open_output_file
    use mp, only: nproc, proc0
    implicit none
    type(optimisation_type), intent(inout) :: optim
    type(optimisation_config_type), intent(in), optional :: optimisation_config_in
    if(initialized) return
    initialized = .true.
    call read_parameters(optim, optimisation_config_in)
    if(optim%on) then
      if (proc0) call open_output_file(optim%outunit, '.optim')
    end if
    optim%nproc_max = nproc
  end subroutine init_optimisation_config

  subroutine finish_optimisation_config(optim)
    use file_utils, only: close_output_file
    use mp, only: proc0
    implicit none
    type(optimisation_type), intent(inout) :: optim
    initialized = .false.

    if(optim%on) then
      if (proc0) call close_output_file(optim%outunit)
    end if

    call optimisation_config%reset()

  end subroutine finish_optimisation_config

  subroutine read_parameters(optim, optimisation_config_in)
    use mp, only: proc0
    use file_utils, only: error_unit
    implicit none
    type(optimisation_type), intent(inout) :: optim
    type(optimisation_config_type), intent(in), optional :: optimisation_config_in

    if (present(optimisation_config_in)) optimisation_config = optimisation_config_in

    call optimisation_config%init(name = 'optimisation_config', requires_index = .false.)

    ! Copy out internal values into module level parameters
    optim%auto = optimisation_config%auto
    optim%estimate_timing_error = optimisation_config%estimate_timing_error
    optim%max_imbalance = optimisation_config%max_imbalance
    optim%max_unused_procs = optimisation_config%max_unused_procs
    optim%measure_all = optimisation_config%measure_all
    optim%min_efficiency = optimisation_config%min_efficiency
    optim%nstep_measure = optimisation_config%nstep_measure
    optim%on = optimisation_config%on
    optim%warm_up = optimisation_config%warm_up

    if (optim%nstep_measure < 2) then
       optim%nstep_measure = 2
       if (proc0) write(error_unit(), '("Warning nstep_measure must be at least 2. Forcing nstep_measure = 2")')
    end if

  end subroutine read_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_optimisation_config(optimisation_config_in)
    use mp, only: mp_abort
    type(optimisation_config_type), intent(in), optional :: optimisation_config_in
    if (initialized) then
       call mp_abort("Trying to set optimisation_config when already initialized.", to_screen = .true.)
    end if
    if (present(optimisation_config_in)) then
       optimisation_config = optimisation_config_in
    end if
  end subroutine set_optimisation_config

#include "optimisation_auto_gen.inc"
end module optimisation_configuration

