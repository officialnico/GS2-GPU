
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
  integer :: in_file

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
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = optimisation_config)

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
  integer, intent(in), optional:: unit
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

!> Get the module level config instance
pure function get_optimisation_config()
  type(optimisation_config_type) :: get_optimisation_config
  get_optimisation_config = optimisation_config
end function get_optimisation_config
