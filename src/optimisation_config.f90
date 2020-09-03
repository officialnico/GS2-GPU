!> A module for handling the configuration of the optimisation
!! module via the namelist optimisation_config.
module optimisation_config
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use overrides, only: optimisations_overrides_type
  implicit none

  private

  public :: init_optimisation_config
  public :: finish_optimisation_config
  public :: optimisation_type
  public :: optimisation_results_type

  public :: optimisation_config_type
  public :: set_optimisation_config_config

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
          dimension(:), pointer :: sorted_optimisations
     type(optimisation_results_type), dimension(:), pointer :: sorted_results
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
  type, extends(abstract_config_type) :: optimisation_config_type
     ! namelist : optimisation_config
     ! indexed : false     
     !> When true, automatically continues GS2 to run the input file with the
     !> optimised parameters.
     logical :: auto = .true.
     !> FIXME: Add documentation
     logical :: estimate_timing_error = .true.
     !> FIXME: Add documentation
     real :: max_imbalance = -1.0
     !> FIXME: Add documentation
     integer :: max_unused_procs = 0
     !> FIXME: Add documentation
     logical :: measure_all = .false.
     !> FIXME: Add documentation
     real :: min_efficiency = -1.0
     !> FIXME: Add documentation
     integer :: nstep_measure = 5
     !> Set true to turn on optimisation procedure
     logical :: on = .false.
     !> FIXME: Add documentation
     logical :: warm_up = .false.
   contains
     procedure, public :: read => read_optimisation_config
     procedure, public :: write => write_optimisation_config
     procedure, public :: reset => reset_optimisation_config
     procedure, public :: broadcast => broadcast_optimisation_config
     procedure, public, nopass :: get_default_name => get_default_name_optimisation_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_optimisation_config
  end type optimisation_config_type
  
  type(optimisation_config_type) :: optimisation_config_config  
  
contains
  subroutine init_optimisation_config(optim, optimisation_config_config_in)
    use file_utils, only: open_output_file
    use mp, only: nproc
    implicit none
    type(optimisation_type), intent(inout) :: optim
    type(optimisation_config_type), intent(in), optional :: optimisation_config_config_in
    if(initialized) return
    initialized = .true.
    call read_parameters(optim, optimisation_config_config_in)
    call open_output_file(optim%outunit, '.optim')
    optim%nproc_max = nproc
  end subroutine init_optimisation_config

  subroutine finish_optimisation_config(optim)
    implicit none
    type(optimisation_type), intent(inout) :: optim
    initialized = .false.
    call optimisation_config_config%reset()
  end subroutine finish_optimisation_config

  subroutine read_parameters(optim, optimisation_config_config_in)
    implicit none
    type(optimisation_type), intent(inout) :: optim
    type(optimisation_config_type), intent(in), optional :: optimisation_config_config_in

    if (present(optimisation_config_config_in)) optimisation_config_config = optimisation_config_config_in

    call optimisation_config_config%init(name = 'optimisation_config', requires_index = .false.)

    ! Copy out internal values into module level parameters
    optim%auto = optimisation_config_config%auto
    optim%estimate_timing_error = optimisation_config_config%estimate_timing_error
    optim%max_imbalance = optimisation_config_config%max_imbalance
    optim%max_unused_procs = optimisation_config_config%max_unused_procs
    optim%measure_all = optimisation_config_config%measure_all
    optim%min_efficiency = optimisation_config_config%min_efficiency
    optim%nstep_measure = optimisation_config_config%nstep_measure
    optim%on = optimisation_config_config%on
    optim%warm_up = optimisation_config_config%warm_up
  end subroutine read_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_optimisation_config_config(optimisation_config_config_in)
    use mp, only: mp_abort
    type(optimisation_config_type), intent(in), optional :: optimisation_config_config_in
    if (initialized) then
       call mp_abort("Trying to set optimisation_config_config when already initialized.", to_screen = .true.)
    end if
    if (present(optimisation_config_config_in)) then
       optimisation_config_config = optimisation_config_config_in
    end if
  end subroutine set_optimisation_config_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the optimisation_config namelist and populates the member variables
  subroutine read_optimisation_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(optimisation_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: max_unused_procs, nstep_measure
    logical :: auto, estimate_timing_error, measure_all, on, warm_up
    real :: max_imbalance, min_efficiency

    namelist /optimisation_config/ auto, estimate_timing_error, max_imbalance, max_unused_procs, measure_all, min_efficiency, nstep_measure, on, warm_up

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    auto = self%auto
    estimate_timing_error = self%estimate_timing_error
    max_imbalance = self%max_imbalance
    max_unused_procs = self%max_unused_procs
    measure_all = self%measure_all
    min_efficiency = self%min_efficiency
    nstep_measure = self%nstep_measure
    on = self%on
    warm_up = self%warm_up

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = optimisation_config, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = optimisation_config)
    endif

    ! Now copy from local variables into type members
    self%auto = auto
    self%estimate_timing_error = estimate_timing_error
    self%max_imbalance = max_imbalance
    self%max_unused_procs = max_unused_procs
    self%measure_all = measure_all
    self%min_efficiency = min_efficiency
    self%nstep_measure = nstep_measure
    self%on = on
    self%warm_up = warm_up

    self%exist = exist
  end subroutine read_optimisation_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_optimisation_config(self, unit)
    implicit none
    class(optimisation_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("auto", self%auto, unit_internal)
    call self%write_key_val("estimate_timing_error", self%estimate_timing_error, unit_internal)
    call self%write_key_val("max_imbalance", self%max_imbalance, unit_internal)
    call self%write_key_val("max_unused_procs", self%max_unused_procs, unit_internal)
    call self%write_key_val("measure_all", self%measure_all, unit_internal)
    call self%write_key_val("min_efficiency", self%min_efficiency, unit_internal)
    call self%write_key_val("nstep_measure", self%nstep_measure, unit_internal)
    call self%write_key_val("on", self%on, unit_internal)
    call self%write_key_val("warm_up", self%warm_up, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_optimisation_config

  !> Resets the config object to the initial empty state
  subroutine reset_optimisation_config(self)
    class(optimisation_config_type), intent(in out) :: self
    type(optimisation_config_type) :: empty
    select type (self)
    type is (optimisation_config_type)
       self = empty
    end select
  end subroutine reset_optimisation_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_optimisation_config(self)
    use mp, only: broadcast
    implicit none
    class(optimisation_config_type), intent(in out) :: self
    call broadcast(self%auto)
    call broadcast(self%estimate_timing_error)
    call broadcast(self%max_imbalance)
    call broadcast(self%max_unused_procs)
    call broadcast(self%measure_all)
    call broadcast(self%min_efficiency)
    call broadcast(self%nstep_measure)
    call broadcast(self%on)
    call broadcast(self%warm_up)

    call broadcast(self%exist)
  end subroutine broadcast_optimisation_config
  
  !> Gets the default name for this namelist
  function get_default_name_optimisation_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_optimisation_config
    get_default_name_optimisation_config = "optimisation_config"
  end function get_default_name_optimisation_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_optimisation_config()
    implicit none
    logical :: get_default_requires_index_optimisation_config
    get_default_requires_index_optimisation_config = .false.
  end function get_default_requires_index_optimisation_config
end module optimisation_config

