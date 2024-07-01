
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the knobs namelist and populates the member variables
subroutine read_knobs_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(knobs_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: delt_option, eqzip_option
  character(len = 65000) :: user_comments
  integer :: ncheck_stop, nstep, progress_frequency, seed
  logical :: do_eigsolve, immediate_reset, neo_test, save_init_times, save_timer_statistics, use_old_diagnostics, wstar_units
  real :: avail_cpu_time, beta, delt, fapar, fbpar, fphi, k0, margin_cpu_time, max_sim_time, rhostar, tite, zeff

  namelist /knobs/ avail_cpu_time, beta, delt, delt_option, do_eigsolve, eqzip_option, fapar, fbpar, fphi, immediate_reset, k0, margin_cpu_time, max_sim_time, ncheck_stop, neo_test, nstep, progress_frequency, rhostar, &
    save_init_times, save_timer_statistics, seed, tite, use_old_diagnostics, user_comments, wstar_units, zeff

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  avail_cpu_time = self%avail_cpu_time
  beta = self%beta
  delt = self%delt
  delt_option = self%delt_option
  do_eigsolve = self%do_eigsolve
  eqzip_option = self%eqzip_option
  fapar = self%fapar
  fbpar = self%fbpar
  fphi = self%fphi
  immediate_reset = self%immediate_reset
  k0 = self%k0
  margin_cpu_time = self%margin_cpu_time
  max_sim_time = self%max_sim_time
  ncheck_stop = self%ncheck_stop
  neo_test = self%neo_test
  nstep = self%nstep
  progress_frequency = self%progress_frequency
  rhostar = self%rhostar
  save_init_times = self%save_init_times
  save_timer_statistics = self%save_timer_statistics
  seed = self%seed
  tite = self%tite
  use_old_diagnostics = self%use_old_diagnostics
  user_comments = self%user_comments
  wstar_units = self%wstar_units
  zeff = self%zeff

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = knobs)

  ! Now copy from local variables into type members
  self%avail_cpu_time = avail_cpu_time
  self%beta = beta
  self%delt = delt
  self%delt_option = delt_option
  self%do_eigsolve = do_eigsolve
  self%eqzip_option = eqzip_option
  self%fapar = fapar
  self%fbpar = fbpar
  self%fphi = fphi
  self%immediate_reset = immediate_reset
  self%k0 = k0
  self%margin_cpu_time = margin_cpu_time
  self%max_sim_time = max_sim_time
  self%ncheck_stop = ncheck_stop
  self%neo_test = neo_test
  self%nstep = nstep
  self%progress_frequency = progress_frequency
  self%rhostar = rhostar
  self%save_init_times = save_init_times
  self%save_timer_statistics = save_timer_statistics
  self%seed = seed
  self%tite = tite
  self%use_old_diagnostics = use_old_diagnostics
  self%user_comments = user_comments
  self%wstar_units = wstar_units
  self%zeff = zeff

  self%exist = exist
end subroutine read_knobs_config

!> Writes out a namelist representing the current state of the config object
subroutine write_knobs_config(self, unit)
  implicit none
  class(knobs_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("avail_cpu_time", self%avail_cpu_time, unit_internal)
  call self%write_key_val("beta", self%beta, unit_internal)
  call self%write_key_val("delt", self%delt, unit_internal)
  call self%write_key_val("delt_option", self%delt_option, unit_internal)
  call self%write_key_val("do_eigsolve", self%do_eigsolve, unit_internal)
  call self%write_key_val("eqzip_option", self%eqzip_option, unit_internal)
  call self%write_key_val("fapar", self%fapar, unit_internal)
  call self%write_key_val("fbpar", self%fbpar, unit_internal)
  call self%write_key_val("fphi", self%fphi, unit_internal)
  call self%write_key_val("immediate_reset", self%immediate_reset, unit_internal)
  call self%write_key_val("k0", self%k0, unit_internal)
  call self%write_key_val("margin_cpu_time", self%margin_cpu_time, unit_internal)
  call self%write_key_val("max_sim_time", self%max_sim_time, unit_internal)
  call self%write_key_val("ncheck_stop", self%ncheck_stop, unit_internal)
  call self%write_key_val("neo_test", self%neo_test, unit_internal)
  call self%write_key_val("nstep", self%nstep, unit_internal)
  call self%write_key_val("progress_frequency", self%progress_frequency, unit_internal)
  call self%write_key_val("rhostar", self%rhostar, unit_internal)
  call self%write_key_val("save_init_times", self%save_init_times, unit_internal)
  call self%write_key_val("save_timer_statistics", self%save_timer_statistics, unit_internal)
  call self%write_key_val("seed", self%seed, unit_internal)
  call self%write_key_val("tite", self%tite, unit_internal)
  call self%write_key_val("use_old_diagnostics", self%use_old_diagnostics, unit_internal)
  call self%write_key_val("user_comments", self%user_comments, unit_internal)
  call self%write_key_val("wstar_units", self%wstar_units, unit_internal)
  call self%write_key_val("zeff", self%zeff, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_knobs_config

!> Resets the config object to the initial empty state
subroutine reset_knobs_config(self)
  class(knobs_config_type), intent(in out) :: self
  type(knobs_config_type) :: empty
  select type (self)
  type is (knobs_config_type)
     self = empty
  end select
end subroutine reset_knobs_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_knobs_config(self)
  use mp, only: broadcast
  implicit none
  class(knobs_config_type), intent(in out) :: self
  call broadcast(self%avail_cpu_time)
  call broadcast(self%beta)
  call broadcast(self%delt)
  call broadcast(self%delt_option)
  call broadcast(self%do_eigsolve)
  call broadcast(self%eqzip_option)
  call broadcast(self%fapar)
  call broadcast(self%fbpar)
  call broadcast(self%fphi)
  call broadcast(self%immediate_reset)
  call broadcast(self%k0)
  call broadcast(self%margin_cpu_time)
  call broadcast(self%max_sim_time)
  call broadcast(self%ncheck_stop)
  call broadcast(self%neo_test)
  call broadcast(self%nstep)
  call broadcast(self%progress_frequency)
  call broadcast(self%rhostar)
  call broadcast(self%save_init_times)
  call broadcast(self%save_timer_statistics)
  call broadcast(self%seed)
  call broadcast(self%tite)
  call broadcast(self%use_old_diagnostics)
  call broadcast(self%user_comments)
  call broadcast(self%wstar_units)
  call broadcast(self%zeff)

  call broadcast(self%exist)
end subroutine broadcast_knobs_config

!> Gets the default name for this namelist
function get_default_name_knobs_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_knobs_config
  get_default_name_knobs_config = "knobs"
end function get_default_name_knobs_config

!> Gets the default requires index for this namelist
function get_default_requires_index_knobs_config()
  implicit none
  logical :: get_default_requires_index_knobs_config
  get_default_requires_index_knobs_config = .false.
end function get_default_requires_index_knobs_config

!> Get the module level config instance
pure function get_knobs_config()
  type(knobs_config_type) :: get_knobs_config
  get_knobs_config = knobs_config
end function get_knobs_config
