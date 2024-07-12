!> FIXME : Add documentation
module gs2_reinit
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: reset_time_step, delt_adj, in_memory
  public :: check_time_step, time_reinit
  public :: init_reinit, wnml_gs2_reinit
  public :: reduce_time_step, increase_time_step
  public :: init_gs2_reinit, finish_gs2_reinit
  public :: reinit_config_type, set_gs2_reinit_config, get_reinit_config
  
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
     !> fall back to the traditional file based approach.
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
  subroutine wnml_gs2_reinit(unit)
    implicit none
    integer, intent(in) :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "reinit_knobs"
    write (unit, fmt="(' delt_adj = ',e17.10)") delt_adj
    write (unit, fmt="(' delt_minimum = ',e17.10)") delt_minimum
    write (unit, fmt="(' /')")       
  end subroutine wnml_gs2_reinit

  !> Reduces the time step by a factor `delt_adj`.
  subroutine reduce_time_step
    use gs2_time, only: code_dt
    implicit none
    if (first) call init_reinit
    code_dt = code_dt/delt_adj
  end subroutine reduce_time_step

  !> Increases the time step by a factor `delt_adj` up to a limit of
  !> `code_dt_max` set by [[reinit_config_type:dt0]] or from the
  !> restart file.
  subroutine increase_time_step
    use gs2_time, only: code_dt, code_dt_max
    implicit none
    if (first) call init_reinit
    code_dt = min(code_dt*delt_adj, code_dt_max)
  end subroutine increase_time_step

  !> FIXME : Add documentation  
  subroutine reset_time_step (current_init, istep, my_exit, job_id)
    use run_parameters, only: reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt, user_time
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt_min, code_dt_max
    use gs2_init, only: init_type, init, init_level_list
    use mp, only: proc0
    use file_utils, only: error_unit
    use job_manage, only: time_message
    use array_utils, only: zero_array
    use optionals, only: get_option_with_default
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
    integer :: timestep_change_attempts, local_job_id
    integer, parameter :: timestep_change_attempts_limit = 10
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

    local_job_id = get_option_with_default(job_id, -1)
    if (local_job_id < 0) call time_message(proc0,time_reinit,' Re-initialize')

    !First disable the reset flag so we can call 
    !routines needed in reinit
    reset_in=reset
    reset=.false.

    !call save_fields_and_dist_fn

    ! Not clear that we really need to do this?
    call zero_array(gnew)

    ! Move to the correct init level
    call init(current_init, init_level_list%override_timestep)
! change timestep 

    original_time_step = code_dt

    timestep_change_attempts = 0

    ! Keep adjusting the timestep until it satisfies the required conditions
    do while (timestep_change_attempts < timestep_change_attempts_limit)
       ! If timestep is too big, make it smaller
       if (code_dt*fac > code_dt_cfl) then
          call reduce_time_step
          ! If timestep is too small, make it bigger
       else if (code_dt*fac < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion)) then
          call increase_time_step
       else
          exit
       endif
       timestep_change_attempts = timestep_change_attempts + 1
    end do

    ! Handle the case where it took too many iterations for the time step
    ! to satisfy the required conditions.
    if(timestep_change_attempts == timestep_change_attempts_limit) then
       if (proc0) write(error_unit(), '("Attempting to change the timestep too much in one go")')
       reset = reset_in
       my_exit = .true.
       ! Note here we leave with our init level less than full. We expect to be
       ! aborting the run here, so that's probably fine but something to be aware of.
       return
    end if

    call save_dt (code_dt)

    ! Check we still aren't below the minimum step
    if (code_dt <= code_dt_min) then
       code_dt = code_dt_min  ! set it so restart is ok
       my_exit = .true.
       if (proc0) write(error_unit(), *) 'Time step wants to fall below delt_min.  Abort run.'
       ! Note here we leave with our init level less than full. We expect to be
       ! aborting the run here, so that's probably fine but something to be aware of.
       return
    end if

    if (proc0 .and. (local_job_id < 0)) write(*,*) 'Changing time step to ', user_dt, ' from ', original_time_step, ' on step number ', istep, ' i.e. time = ', user_time

    ! Don't reset antenna here because species parameters
    ! have not changed so resetting antenna would cause
    ! an unnecessary discontinuity
    !call reinit_gk_and_field_equations(reset_antenna=.false.)
    call init(current_init, init_level_list%full)
    
    if (local_job_id < 0) call time_message(proc0,time_reinit,' Re-initialize')

    istep_last = istep

    !Now re-enable reset so we leave it in the same state as on entering
    reset=reset_in

  end subroutine reset_time_step

  !> FIXME : Add documentation  
  subroutine check_time_step (reset, exit)
    use gs2_time, only: code_dt_cfl, code_dt, code_dt_max
    use nonlinear_terms, only: nb_check_time_step_too_large
    use mp, only: broadcast
    implicit none
    logical, intent(in) :: exit
    logical, intent(out) :: reset

    real :: fac = 1.0

    if (first) call init_reinit
    first = .false.
    reset = .false.

! nothing to do if exiting in this iteration
    if (exit) return

! If doing nonblocking CFL check, finish it here
    call nb_check_time_step_too_large

! If timestep is too big, make it smaller
    if (code_dt*fac > code_dt_cfl) reset = .true. !Note this logic is repeated in gs2_time::check_time_step_too_large
       
! If timestep is too small, make it bigger
    if (code_dt*fac < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.

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
    use file_utils, only: input_unit, input_unit_exist
    use gs2_time, only: save_dt_min, save_dt_max, code_dt_max
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
       reinit_config%dt0 = code_dt_max
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
    call save_dt_max (dt0)
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

#include "reinit_auto_gen.inc"
end module gs2_reinit

