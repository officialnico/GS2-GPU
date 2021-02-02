!> FIXME : Add documentation
module gs2_reinit
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: reset_time_step, delt_adj
  public :: check_time_step, time_reinit
  public :: init_reinit, wnml_gs2_reinit
  public :: reduce_time_step, increase_time_step
  public :: init_gs2_reinit, finish_gs2_reinit

  !!> This subroutine reinitializes the modules which are responsible for
  !!! solving the linear, collisonal and nonlinear parts of the GK eqn,
  !!! as well as the modules which solve the field equations. This is 
  !!! typically necessary because some parameter (e.g. the timestep) which
  !!! is used to calculate the value of cached arrays (e.g. the response
  !!! matrix) in those modules has changed. Note that this will not cause
  !!! dist_fn to reread its namelist.. thus you can change g_exb (and any
  !!! other physics parameter in any other module, of course) 
  !!! before calling this function 
  !!! and expect your changes to be preserved. The function 
  !!! save_fields_and_dist_fn must be called before this, or 
  !!! otherwise the current field and dist_fn values will be lost. The
  !!! logical flag in_memory must be given the same value that was 
  !!! set in save_fields_and_dist_fn. 
  !public :: reinit_gk_and_field_equations

  !> This function overrides the in_memory flag
  !! and should only be used if you know what
  !! you are doing.
  public :: gs2_reinit_unit_test_set_in_memory

  public :: reinit_config_type
  public :: set_gs2_reinit_config
  public :: get_gs2_reinit_config
  
  real :: delt_adj, dt0
  real :: delt_cushion
  real :: delt_minimum 
  real :: time_reinit(2)=0.
  logical :: abort_rapid_time_step_change
  logical :: first=.true.
  logical :: in_memory
  logical :: initialized = .false.
  
  !> Used to represent the input configuration of reinit
  type, extends(abstract_config_type) :: reinit_config_type
     ! namelist : reinit_knobs
     ! indexed : false
     !> If `true` (default), exit if time step changes rapidly, that
     !> is, if the time step changes at four consecutive time steps.
     logical :: abort_rapid_time_step_change = .true.
     !> When the time step needs to be changed it is adjusted by this
     !> factor, i.e `dt --> dt/delt_adj` or `dt --> dt*delt_adj` when
     !> reducing/increasing the timestep. For non-linear runs
     !> good choice of `delt_adj` can make a moderate difference to
     !> efficiency. Need to balance time taken to reinitialise against
     !> frequency of time step adjustments (i.e. if your run takes a long
     !> time to initialise you probably want to set `delt_adj` to be
     !> reasonably large).
     real :: delt_adj = 2.0
     !> Used in deciding when to increase the time step to help
     !> prevent oscillations in time step around some value. We only
     !> increase the time step when it is less than the scaled cfl
     !> estimate divided by `delt_adj*delt_cushion` whilst we decrease
     !> it as soon as the time step is larger than the scaled cfl
     !> estimate.
     real :: delt_cushion = 1.5
     !> The minimum time step allowed is delt_minimum. If the code
     !> wants to drop below this value then the run will end.
     real :: delt_minimum = 1.e-5
     !> Sets the maximum value the time step can take.
     !>
     !> @note This gets a smart default of [[knobs:delt]].
     real :: dt0 = 0.0
     !> If `true` then attempts to create temporary copies of the
     !> distribution fn and fields in memory to be restored after the
     !> time step reset rather than dumping to fields.  This could be
     !> faster on machines with slow file systems. If the required
     !> memory allocation fails then we set `in_memory=.false.` and
     !> fall back to the traditional file based approach. Note must
     !> also set the `in_memory` switch in [[init_knobs]] to `true` to
     !> enable.
     logical :: in_memory = .false.
   contains
     procedure, public :: read => read_reinit_config
     procedure, public :: write => write_reinit_config
     procedure, public :: reset => reset_reinit_config
     procedure, public :: broadcast => broadcast_reinit_config
     procedure, public, nopass :: get_default_name => get_default_name_reinit_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_reinit_config
  end type reinit_config_type

  type(reinit_config_type) :: reinit_config  
contains

  !> FIXME : Add documentation  
  subroutine gs2_reinit_unit_test_set_in_memory(in_memory_in)
    use gs2_init, only: gs2_init_in_memory=>in_memory
    implicit none
    logical, intent(in) :: in_memory_in
    in_memory = in_memory_in
    gs2_init_in_memory = in_memory
  end subroutine gs2_reinit_unit_test_set_in_memory
  
  !> FIXME : Add documentation  
  subroutine wnml_gs2_reinit(unit)
    implicit none
    integer, intent(in) :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "reinit_knobs"
    write (unit, fmt="(' delt_adj = ',e17.10)") delt_adj
    write (unit, fmt="(' delt_minimum = ',e17.10)") delt_minimum
    write (unit, fmt="(' /')")       
  end subroutine wnml_gs2_reinit

  !> FIXME : Add documentation  
  subroutine reduce_time_step
    use gs2_time, only: code_dt
    implicit none
    if (first) call init_reinit
    code_dt = code_dt/delt_adj
  end subroutine reduce_time_step

  !> FIXME : Add documentation  
  subroutine increase_time_step
    use gs2_time, only: code_dt
    implicit none
    if (first) call init_reinit
    code_dt = min(code_dt*delt_adj, dt0)
  end subroutine increase_time_step

  !> FIXME : Add documentation  
  subroutine reset_time_step (current_init, istep, my_exit, job_id)
    use run_parameters, only: reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt, user_time
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt_min
    use gs2_init, only: init_type, init
    use gs2_init, only: init_level_list
    use mp, only: proc0
    use file_utils, only: error_unit
    use job_manage, only: time_message
    use nonlinear_terms, only: gryfx_zonal
    implicit none
    integer, intent(in) :: istep 
    logical, intent(inout) :: my_exit
    integer, intent (in), optional :: job_id
    logical :: reset_in
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer, save :: nconsec=0
    type(init_type), intent(inout) :: current_init
    real :: original_time_step
    real :: fac = 1.0

    if (first) call init_reinit
    first = .false.

! save fields and distribution function

! calls on consecutive time steps is probably an error
    if (istep_last + 1 == istep) then
       nconsec=nconsec+1
    else
       nconsec=0
    endif

    if (nconsec .gt. 4 .and. abort_rapid_time_step_change) then
       my_exit = .true.
       if (proc0) write(error_unit(), *) 'Time step changing rapidly.  Abort run.'
       return
    end if

    if (code_dt/delt_adj <= code_dt_min) then
       code_dt = code_dt_min  ! set it so restart is ok
       my_exit = .true.
       if (proc0) write(error_unit(), *) 'Time step wants to fall below delt_min.  Abort run.'
       return
    end if

    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    !First disable the reset flag so we can call 
    !routines needed in reinit
    reset_in=reset
    reset=.false.

    !call save_fields_and_dist_fn

    gnew = 0.

    ! Move to the correct init level
    call init(current_init, init_level_list%override_timestep)
! change timestep 

    if(gryfx_zonal%on) then
!both code_dt = dt_gs2 and code_dt_cfl = dt_cfl_gryfx are in gs2 units
!we want to check if dt_gryfx = 2*dt_gs2 is too big/small when compared to
!dt_cfl_gryfx
      fac = 1.9999
    endif

    original_time_step = code_dt
    
! If timestep is too big, make it smaller
    if (code_dt*fac > code_dt_cfl) then
       call reduce_time_step
! If timestep is too small, make it bigger
    else if (code_dt*fac < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) then
       call increase_time_step
    endif
    
    call save_dt (code_dt)

    if (proc0 .and. .not. present(job_id)) write(*,*) 'Changing time step to ', user_dt, ' from ', original_time_step, ' on step number ', istep, ' i.e. time = ', user_time

    ! Don't reset antenna here because species parameters
    ! have not changed so resetting antenna would cause
    ! an unnecessary discontinuity
    !call reinit_gk_and_field_equations(reset_antenna=.false.)
    call init(current_init, init_level_list%full)
    
    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    istep_last = istep

    !Now re-enable reset so we leave it in the same state as on entering
    reset=reset_in

  end subroutine reset_time_step

  !> FIXME : Add documentation  
  subroutine check_time_step (reset, exit)
    use gs2_time, only: code_dt_cfl, code_dt
    use nonlinear_terms, only: gryfx_zonal, nb_check_time_step_too_large
    use mp, only: broadcast
    implicit none
    logical, intent(in) :: exit
    logical, intent(out) :: reset

    real :: fac = 1.0

    if (first) call init_reinit
    first = .false.
    reset = .false.

    if(gryfx_zonal%on) then
      !code_dt_cfl is only set on proc 0 by gryfx in nlps.cu
      call broadcast(code_dt_cfl)
!both code_dt = dt_gs2 and code_dt_cfl = dt_cfl_gryfx are in gs2 units
!we want to check if dt_gryfx = 2*dt_gs2 is too big/small when compared to
!dt_cfl_gryfx
      fac = 1.9999
    endif

! nothing to do if exiting in this iteration
    if (exit) return

! If doing nonblocking CFL check, finish it here
    call nb_check_time_step_too_large

! If timestep is too big, make it smaller
    if (code_dt*fac > code_dt_cfl) reset = .true. !Note this logic is repeated in gs2_time::check_time_step_too_large
       
! If timestep is too small, make it bigger
    if (code_dt*fac < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.

  end subroutine check_time_step

  !> FIXME : Add documentation
  subroutine init_gs2_reinit
    call init_reinit
  end subroutine init_gs2_reinit

  !> FIXME : Add documentation  
  subroutine finish_gs2_reinit
    first = .true.
    initialized = .false.
    call reinit_config%reset()
  end subroutine finish_gs2_reinit

  !> FIXME : Add documentation  
  subroutine init_reinit(reinit_config_in)
    use run_parameters, only: code_delt_max
    use file_utils, only: input_unit, input_unit_exist
    use gs2_time, only: save_dt_min
    use gs2_init, only: gs2_init_in_memory=>in_memory
    implicit none
    type(reinit_config_type), intent(in), optional :: reinit_config_in    
    logical :: exist

    if(initialized) return
    initialized = .true.
    if(.not.first)return
    first=.false.

    if (present(reinit_config_in)) reinit_config = reinit_config_in

    !Smart defaults
    if (.not.reinit_config%is_initialised()) then    
       reinit_config%dt0 = code_delt_max
    endif
    
    call reinit_config%init(name = 'reinit_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    abort_rapid_time_step_change = reinit_config%abort_rapid_time_step_change
    delt_adj = reinit_config%delt_adj
    delt_cushion = reinit_config%delt_cushion
    delt_minimum = reinit_config%delt_minimum
    dt0 = reinit_config%dt0
    in_memory = reinit_config%in_memory
    
    exist = reinit_config%exist
    
    call save_dt_min (delt_minimum)
    
    ! Override with the value from gs2 init
    ! This is not transparent to the user - why do we do this?
    in_memory = gs2_init_in_memory
  end subroutine init_reinit

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_gs2_reinit_config(gs2_reinit_config_in)
    use mp, only: mp_abort
    type(reinit_config_type), intent(in), optional :: gs2_reinit_config_in
    if (initialized) then
       call mp_abort("Trying to set gs2_reinit_config when already initialized.", to_screen = .true.)
    end if
    if (present(gs2_reinit_config_in)) then
       reinit_config = gs2_reinit_config_in
    end if
  end subroutine set_gs2_reinit_config

  !> Get the module level config instance
  function get_gs2_reinit_config()
    type(reinit_config_type) :: get_gs2_reinit_config
    get_gs2_reinit_config = reinit_config
  end function get_gs2_reinit_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the reinit_knobs namelist and populates the member variables
  subroutine read_reinit_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(reinit_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: abort_rapid_time_step_change, in_memory
    real :: delt_adj, delt_cushion, delt_minimum, dt0

    namelist /reinit_knobs/ abort_rapid_time_step_change, delt_adj, delt_cushion, delt_minimum, dt0, in_memory

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    abort_rapid_time_step_change = self%abort_rapid_time_step_change
    delt_adj = self%delt_adj
    delt_cushion = self%delt_cushion
    delt_minimum = self%delt_minimum
    dt0 = self%dt0
    in_memory = self%in_memory

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = reinit_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = reinit_knobs)
    endif

    ! Now copy from local variables into type members
    self%abort_rapid_time_step_change = abort_rapid_time_step_change
    self%delt_adj = delt_adj
    self%delt_cushion = delt_cushion
    self%delt_minimum = delt_minimum
    self%dt0 = dt0
    self%in_memory = in_memory

    self%exist = exist
  end subroutine read_reinit_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_reinit_config(self, unit)
    implicit none
    class(reinit_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("abort_rapid_time_step_change", self%abort_rapid_time_step_change, unit_internal)
    call self%write_key_val("delt_adj", self%delt_adj, unit_internal)
    call self%write_key_val("delt_cushion", self%delt_cushion, unit_internal)
    call self%write_key_val("delt_minimum", self%delt_minimum, unit_internal)
    call self%write_key_val("dt0", self%dt0, unit_internal)
    call self%write_key_val("in_memory", self%in_memory, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_reinit_config

  !> Resets the config object to the initial empty state
  subroutine reset_reinit_config(self)
    class(reinit_config_type), intent(in out) :: self
    type(reinit_config_type) :: empty
    select type (self)
    type is (reinit_config_type)
       self = empty
    end select
  end subroutine reset_reinit_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_reinit_config(self)
    use mp, only: broadcast
    implicit none
    class(reinit_config_type), intent(in out) :: self
    call broadcast(self%abort_rapid_time_step_change)
    call broadcast(self%delt_adj)
    call broadcast(self%delt_cushion)
    call broadcast(self%delt_minimum)
    call broadcast(self%dt0)
    call broadcast(self%in_memory)

    call broadcast(self%exist)
  end subroutine broadcast_reinit_config

  !> Gets the default name for this namelist
  function get_default_name_reinit_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_reinit_config
    get_default_name_reinit_config = "reinit_knobs"
  end function get_default_name_reinit_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_reinit_config()
    implicit none
    logical :: get_default_requires_index_reinit_config
    get_default_requires_index_reinit_config = .false.
  end function get_default_requires_index_reinit_config
end module gs2_reinit

