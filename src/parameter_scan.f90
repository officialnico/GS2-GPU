!> A module which allows multiple values of certain parameters to be considered within a single simulation.
!!
!! In its simplest form, it starts with a given value of the parameter in question, runs for a given number of timesteps, and then changes the parameter by a given increment, until a lower limit is reached.
!!
!! In a more advanced scenario, the parameter scan continues until a different condition is satisfied: for example, zero heat flux.
!!
!! In a third scenario, the parameter is varied using a root finding algorithm.
!!
!! In addition, the condition for changing the parameter may be changed from a simple number of time steps to, for example, reaching a saturated state.
module parameter_scan
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  !These should really be moved to specific routines where required
  use parameter_scan_arrays, only: current_scan_parameter_value, write_scan_parameter, run_scan
  use parameter_scan_arrays, only: scan_parameter_switch, scan_parameter_tprim, scan_parameter_g_exb
  use parameter_scan_arrays, only: scan_spec, hflux_tot, momflux_tot, phi2_tot, nout

  implicit none

  private

  public :: init_parameter_scan
  public :: finish_parameter_scan
  public :: update_scan_parameter_value
  public :: allocate_target_arrays
  public :: deallocate_target_arrays
  public :: target_parameter_switch, scan_type_switch, scan_type_none
  public :: target_parameter_hflux_tot, target_parameter_momflux_tot, target_parameter_phi2_tot
  public :: scan_restarted

  public :: parameter_scan_config_type
  public :: set_parameter_scan_config
  public :: get_parameter_scan_config
  
  logical :: scan_restarted

  integer :: target_parameter_switch, &
             scan_type_switch, &
             increment_condition_switch


  integer, parameter :: scan_type_none = 1, &
                        scan_type_range = 2, &
                        scan_type_target =3, &
                        scan_type_root_finding = 4

  integer, parameter :: target_parameter_hflux_tot = 1, &
                        target_parameter_momflux_tot = 2, &
                        target_parameter_phi2_tot = 3 

  integer, parameter :: increment_condition_n_timesteps = 1, &
                        increment_condition_delta_t = 2, &
                        increment_condition_saturated = 3

  real :: par_start, par_end, par_inc
  integer :: nstep_init, nstep_inc
  real :: delta_t_init, delta_t_inc
  real :: target_val
  integer :: scan_output_file
  real :: current_target_value = 0.0
  logical :: initialized = .false.

  !> Used to represent the input configuration of parameter_scan
  type, extends(abstract_config_type) :: parameter_scan_config_type
     ! namelist : parameter_scan_knobs
     ! indexed : false     

     !> When the increment condition is `"delta_t"`, the parameter will be changed
     !> every time `delta_t_inc` time has elapsed.
     real :: delta_t_inc = 0.0
     !> When the increment condition is `"delta_t"`, the parameter will not be
     !> changed until `delta_t_init` time has elapsed from the beginning of the
     !> simulation. Note, that if the simulation is restarted, this parameter
     !> will measure from beginning of original simulation.
     real :: delta_t_init = 0.0
     !> Specifies the condition for incrementing the parameter. Possible values are:
     !>
     !> - "n_timesteps": change the parameter after a given number of time steps
     !> - "delta_t": change the parameter after an elapsed time
     !> - "saturated": change the parameter after the simulation has reached a
     !>   saturated state (determined using the target parameter) at the current
     !>   value of the parameter
     character(len = 20) :: inc_con = 'delta_t'
     !> When the increment condition is `'n_timesteps'`, the parameter will be
     !> changed every `nstep_inc`.
     integer :: nstep_inc = 0
     !> When the increment condition is 'n_timesteps' or 'saturated', the
     !> parameter will not be changed until nstep_init have elapsed from the
     !> beginning of the simulation. Note that if the simulation is restarted,
     !> this parameter will measure from the restart.
     integer :: nstep_init = 0
     !> If the scan is being run in 'range' mode, specifies the value of the
     !> parameter that will be reached.
     real :: par_end = 0.0
     !> If the parameter scan is being run in 'range' or 'target' modes,
     !> specifies the amount by which the parameter is varied at one go.
     real :: par_inc = 0.0
     !> Specifies the starting value for the parameter scan.
     real :: par_start = 0.0
     !> Specify the parameter to be varied. If the parameter pertains to a
     !> species, the scan_spec must be specified as well.
     character(len = 20) :: scan_par = 'tprim'
     !> Specify the parameter to be varied. If the parameter pertains to a
     !> species, the scan_spec must be specified as well.
     logical :: scan_restarted = .false.
     !> When parameter pertains to a species, specifies the index of the
     !> species.
     integer :: scan_spec = 1
     !> Specifies the way that the parameter scan is conducted. Possible values are:
     !>
     !> - 'none': do not conduct a parameter scan (default)
     !> - 'range': vary parameter in constant increments between 2 values:
     !>   par_start and par_end. The step size is given by par_inc.
     !> - 'target': start with the parameter at par_start, and then change the
     !>   parameter by par_inc until the target parameter has reached the target
     !>   value
     !> - 'root_finding': the same as target, but the increment is changed
     !>   intelligently using a Newton-like method.
     character(len = 20) :: scan_type = 'none'
     !> If the scan is being run in 'target' or 'root_finding' mode, specifies
     !> the target parameter. Possible values are:
     !>
     !> - 'hflux_tot'
     !> - 'momflux_tot'
     !> - 'phi2_tot'
     character(len = 20) :: target_par = 'hflux_tot'
     !> 
     real :: target_val = 0.0
   contains
     procedure, public :: read => read_parameter_scan_config
     procedure, public :: write => write_parameter_scan_config
     procedure, public :: reset => reset_parameter_scan_config
     procedure, public :: broadcast => broadcast_parameter_scan_config
     procedure, public, nopass :: get_default_name => get_default_name_parameter_scan_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_parameter_scan_config
  end type parameter_scan_config_type
  
  type(parameter_scan_config_type) :: parameter_scan_config
  
contains

  !> FIXME : Add documentation
  subroutine init_parameter_scan(parameter_scan_config_in)
    use gs2_save, only: restore_current_scan_parameter_value 
    use init_g, only: init_init_g
    use file_utils, only: open_output_file
    implicit none
    type(parameter_scan_config_type), intent(in), optional :: parameter_scan_config_in

    if (initialized) return
    initialized = .true.


    call read_parameters(parameter_scan_config_in)

    select case (scan_type_switch) 
    case (scan_type_none)
      write_scan_parameter = .false.
      run_scan = .false.
      return
    case (scan_type_target)
      !call allocate_target_arrays
    case (scan_type_root_finding)
      !call allocate_target_arrays
    end select


    call open_output_file(scan_output_file, ".par_scan")
      
    write_scan_parameter = .true.
    run_scan = .true.

    ! scan_restarted is set manually 
    if (scan_restarted) then
      call init_init_g ! To get restart file name
      call restore_current_scan_parameter_value(current_scan_parameter_value)
    end if
   
    !write (*,*) "initialized parameter_scan"

  end subroutine init_parameter_scan

  !> FIXME : Add documentation
  subroutine allocate_target_arrays(nwrite,write_nl_flux)
!    use gs2_diagnostics, only: nwrite, write_nl_flux
!    use gs2_diagnostics, only: write_flux_line
    use run_parameters, only: nstep
    use mp, only : mp_abort
    implicit none
    integer, intent(in) :: nwrite
    logical, intent(in) :: write_nl_flux
    !if (.not. write_flux_line) &
     !call mp_abort("write_flux_line must be set to true for target mode")

    !select case (scan_type_switch) 
    !case (scan_type_none)
      !return
    !case (scan_type_range)
      !return
    !end select

     allocate(hflux_tot(nstep/nwrite+1))
     allocate(momflux_tot(nstep/nwrite+1))
     allocate(phi2_tot(nstep/nwrite+1))
    if (scan_type_switch == scan_type_none) return

    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for hflux target mode")
     case (target_parameter_momflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for momflux target mode")
     end select 
     !write (*,*) "allocating target arrays ", nwrite, ",", nstep 
  end subroutine allocate_target_arrays

  !> FIXME : Add documentation  
  subroutine deallocate_target_arrays
    implicit none
    if (allocated(hflux_tot)) deallocate(hflux_tot)
    if (allocated(momflux_tot)) deallocate(momflux_tot)
    if (allocated(phi2_tot)) deallocate(phi2_tot)
  end subroutine deallocate_target_arrays

  !> FIXME : Add documentation
  subroutine finish_parameter_scan
    use file_utils, only: close_output_file
    implicit none
    call deallocate_target_arrays
    call close_output_file(scan_output_file)
    initialized = .false.
    call parameter_scan_config%reset()
  end subroutine finish_parameter_scan

  !> FIXME : Add documentation  
  subroutine update_scan_parameter_value(istep, reset, exit)
    use mp, only : mp_abort
    use gs2_time, only: user_time
    implicit none
    logical, intent (inout) :: exit, reset
    integer, intent (in) :: istep
    real :: increment
    logical :: increment_condition_satisfied, target_reached

    !write (*,*) "update_scan_parameter_value"
    increment_condition_satisfied = .false.
    target_reached = .false.

    select case (scan_type_switch)
    case (scan_type_none)
      return
    case (scan_type_range)
      write (scan_output_file, fmt="('time=  ',e11.4,'  parameter=  ',e11.4)") &
          user_time, current_scan_parameter_value
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = par_inc
      if ((increment .gt. 0.0 &
         .and. current_scan_parameter_value + increment .gt. par_end) &
         .or. (increment .lt. 0.0 &
         .and. current_scan_parameter_value + increment .lt. par_end)) &
         then
         exit = .true.
       else  
          call increment_scan_parameter(increment, reset)
       end if
     case (scan_type_target)
      write (scan_output_file, &
        fmt="('time=  ',e11.4,'  parameter=  ',e11.4,'  target=  ',e11.4)") &
          user_time, current_scan_parameter_value, current_target_value 
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = par_inc
      call check_target_reached(target_reached)
      if (target_reached) then
        exit = .true.
      else
        call increment_scan_parameter(increment, reset)
      end if
    case (scan_type_root_finding)
      call mp_abort("scan_type_root_finding not implemented yet!")
      !write (scan_output_file, &
        !fmt="(time=  e11.4  parameter=  e11.4  target=  e11.4)") &
          !user_time, current_scan_parameter_value, current_target_value 
      !call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      !if (.not. increment_condition_satisfied) return
      !increment = par_inc
      !call check_target_reached(target_reached)
      !if (target_reached) then
        !exit = .true.
      !else
        !call get_root_finding_increment(increment)
        !call increment_scan_parameter(increment, reset)
      !end if
    end select
    if (exit) write (scan_output_file, *) "Parameter scan is complete...."
  end subroutine update_scan_parameter_value

  !> FIXME : Add documentation  
  subroutine check_target_reached(target_reached)
    implicit none
    logical, intent (out) :: target_reached
    logical, save :: first = .true.
    real, save :: last_value = 0.0
   
    target_reached = .false.
    
    last_value = current_target_value
    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      current_target_value = hflux_tot(nout)
    case (target_parameter_momflux_tot)
      current_target_value = momflux_tot(nout)
    case (target_parameter_phi2_tot)
      current_target_value = phi2_tot(nout)
    end select

    if (first) then
      last_value = current_target_value
      first = .false.
      return
    end if

    if ((last_value .lt. target_val .and. &
          current_target_value .gt. target_val) .or. &
          (last_value .gt. target_val .and. &
          current_target_value .lt. target_val)) &
             target_reached = .true.

  end subroutine check_target_reached

  !> FIXME : Add documentation  
  subroutine increment_scan_parameter(increment, reset)
    implicit none
    real, intent (in) :: increment
    logical, intent (inout) :: reset
    current_scan_parameter_value = current_scan_parameter_value + increment
    call set_scan_parameter(reset)
  end subroutine increment_scan_parameter

  !> FIXME : Add documentation  
  subroutine set_scan_parameter(reset)
    !use parameter_scan_arrays, only: current_scan_parameter_value
    !use parameter_scan_arrays, only: scan_parameter_switch
    !use parameter_scan_arrays, only: scan_parameter_tprim
    !use parameter_scan_arrays, only: scan_parameter_g_exb
    !use parameter_scan_arrays
    use species, only: spec 
    use dist_fn, only: g_exb
    use mp, only: proc0
    implicit none
    logical, intent (inout) :: reset
     
    select case (scan_parameter_switch)
    case (scan_parameter_tprim)
       spec(scan_spec)%tprim = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter tprim_1 to ", spec(scan_spec)%tprim
       reset = .true.
    case (scan_parameter_g_exb)
       g_exb = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter g_exb to ", g_exb
       reset = .false.
    end select
  end subroutine set_scan_parameter

  !> FIXME : Add documentation  
  subroutine check_increment_condition_satisfied(istep, increment_condition_satisfied)
    use gs2_time, only: user_time
    use mp, only : mp_abort, proc0
    implicit none
    logical, intent (out) :: increment_condition_satisfied
    integer, intent (in) :: istep
    integer, save :: last_timestep = 0
    real, save :: last_time = 0.0

    select case (increment_condition_switch)
    case (increment_condition_n_timesteps)
      if (istep .lt. nstep_init) then
        last_timestep = istep
        increment_condition_satisfied = .false.
        return
      end if 
      if ((istep - last_timestep) .gt. nstep_inc) then
        last_timestep = istep
        increment_condition_satisfied = .true.
        if (proc0) write (*,*) 'Updating parameter at ', user_time 
        return
      end if
    case (increment_condition_delta_t)
      if (user_time .lt. delta_t_init) then
        last_time = user_time
        increment_condition_satisfied = .false.
        return
      end if 
      if ((user_time - last_time) .gt. delta_t_inc) then
        last_time = user_time
        increment_condition_satisfied = .true.
        if (proc0) write (*,*) 'Updating parameter at ', user_time 
        return
      end if
    case (increment_condition_saturated)
      call mp_abort("increment_condition_saturated not implemented yet!")
      !if (istep .lt. nstep_init) then
        !last_timestep = istep
        !increment_condition_satisfied = .false.
        !return
      !end if 
      !call get_target_value_array(target_value_array)
      !call check_saturated(target_value_array, &
                           !istep, &
                           !last_timestep, &
                           !increment_condition_satisfied)
      !if (increment_condition_satisfied) last_timestep = istep
    end select 
  end subroutine check_increment_condition_satisfied

  !> FIXME : Add documentation  
  subroutine read_parameters(parameter_scan_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type(parameter_scan_config_type), intent(in), optional :: parameter_scan_config_in
    
    type (text_option), dimension (2), parameter :: scan_parameter_opts = &
         (/ text_option('tprim', scan_parameter_tprim), &
            text_option('g_exb', scan_parameter_g_exb) /)
    character(20) :: scan_par
    type (text_option), dimension (4), parameter :: scan_type_opts = &
         (/ text_option('none', scan_type_none), &
            text_option('range', scan_type_range), &
            text_option('target', scan_type_target), &
            text_option('root_finding', scan_type_root_finding) /)
    character(20) :: scan_type
    type (text_option), dimension (3), parameter :: target_parameter_opts = &
         (/ text_option('hflux_tot', target_parameter_hflux_tot), &
            text_option('momflux_tot', target_parameter_momflux_tot), &
            text_option('phi2_tot', target_parameter_phi2_tot) /)
    character(20) :: target_par
    type (text_option), dimension (3), parameter :: increment_condition_opts = &
         (/ text_option('n_timesteps', increment_condition_n_timesteps), &
            text_option('delta_t', increment_condition_delta_t), &
            text_option('saturated', increment_condition_saturated) /)
    character(20) :: inc_con 
    integer :: ierr
    logical :: exist

    if (present(parameter_scan_config_in)) parameter_scan_config = parameter_scan_config_in

    call parameter_scan_config%init(name = 'parameter_scan_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    delta_t_inc = parameter_scan_config%delta_t_inc
    delta_t_init = parameter_scan_config%delta_t_init
    inc_con = parameter_scan_config%inc_con
    nstep_inc = parameter_scan_config%nstep_inc
    nstep_init = parameter_scan_config%nstep_init
    par_end = parameter_scan_config%par_end
    par_inc = parameter_scan_config%par_inc
    par_start = parameter_scan_config%par_start
    scan_par = parameter_scan_config%scan_par
    scan_restarted = parameter_scan_config%scan_restarted
    scan_spec = parameter_scan_config%scan_spec
    scan_type = parameter_scan_config%scan_type
    target_par = parameter_scan_config%target_par
    target_val = parameter_scan_config%target_val

    exist = parameter_scan_config%exist

    ierr = error_unit()
    call get_option_value &
         (scan_par, scan_parameter_opts, scan_parameter_switch, &
         ierr, "scan_par in parameter_scan_knobs",.true.)
    call get_option_value &
         (scan_type, scan_type_opts, scan_type_switch, &
         ierr, "scan_type in parameter_scan_knobs",.true.)
    call get_option_value &
         (target_par, target_parameter_opts, target_parameter_switch, &
         ierr, "target_par in parameter_scan_knobs",.true.)
    call get_option_value &
         (inc_con, increment_condition_opts, &
         increment_condition_switch, &
         ierr, "inc_con in parameter_scan_knobs",.true.)

    if (.not. scan_restarted) current_scan_parameter_value = par_start

  end subroutine read_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_parameter_scan_config(parameter_scan_config_in)
    use mp, only: mp_abort
    type(parameter_scan_config_type), intent(in), optional :: parameter_scan_config_in
    if (initialized) then
       call mp_abort("Trying to set parameter_scan_config when already initialized.", to_screen = .true.)
    end if
    if (present(parameter_scan_config_in)) then
       parameter_scan_config = parameter_scan_config_in
    end if
  end subroutine set_parameter_scan_config

  !> Get the module level config instance
  function get_parameter_scan_config()
    type(parameter_scan_config_type) :: get_parameter_scan_config
    get_parameter_scan_config = parameter_scan_config
  end function get_parameter_scan_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the parameter_scan_knobs namelist and populates the member variables
  subroutine read_parameter_scan_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(parameter_scan_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: inc_con, scan_par, scan_type, target_par
    integer :: nstep_inc, nstep_init, scan_spec
    logical :: scan_restarted
    real :: delta_t_inc, delta_t_init, par_end, par_inc, par_start, target_val

    namelist /parameter_scan_knobs/ delta_t_inc, delta_t_init, inc_con, nstep_inc, nstep_init, par_end, par_inc, par_start, scan_par, scan_restarted, scan_spec, scan_type, target_par, target_val

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    delta_t_inc = self%delta_t_inc
    delta_t_init = self%delta_t_init
    inc_con = self%inc_con
    nstep_inc = self%nstep_inc
    nstep_init = self%nstep_init
    par_end = self%par_end
    par_inc = self%par_inc
    par_start = self%par_start
    scan_par = self%scan_par
    scan_restarted = self%scan_restarted
    scan_spec = self%scan_spec
    scan_type = self%scan_type
    target_par = self%target_par
    target_val = self%target_val

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = parameter_scan_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = parameter_scan_knobs)
    endif

    ! Now copy from local variables into type members
    self%delta_t_inc = delta_t_inc
    self%delta_t_init = delta_t_init
    self%inc_con = inc_con
    self%nstep_inc = nstep_inc
    self%nstep_init = nstep_init
    self%par_end = par_end
    self%par_inc = par_inc
    self%par_start = par_start
    self%scan_par = scan_par
    self%scan_restarted = scan_restarted
    self%scan_spec = scan_spec
    self%scan_type = scan_type
    self%target_par = target_par
    self%target_val = target_val

    self%exist = exist
  end subroutine read_parameter_scan_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_parameter_scan_config(self, unit)
    implicit none
    class(parameter_scan_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("delta_t_inc", self%delta_t_inc, unit_internal)
    call self%write_key_val("delta_t_init", self%delta_t_init, unit_internal)
    call self%write_key_val("inc_con", self%inc_con, unit_internal)
    call self%write_key_val("nstep_inc", self%nstep_inc, unit_internal)
    call self%write_key_val("nstep_init", self%nstep_init, unit_internal)
    call self%write_key_val("par_end", self%par_end, unit_internal)
    call self%write_key_val("par_inc", self%par_inc, unit_internal)
    call self%write_key_val("par_start", self%par_start, unit_internal)
    call self%write_key_val("scan_par", self%scan_par, unit_internal)
    call self%write_key_val("scan_restarted", self%scan_restarted, unit_internal)
    call self%write_key_val("scan_spec", self%scan_spec, unit_internal)
    call self%write_key_val("scan_type", self%scan_type, unit_internal)
    call self%write_key_val("target_par", self%target_par, unit_internal)
    call self%write_key_val("target_val", self%target_val, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_parameter_scan_config

  !> Resets the config object to the initial empty state
  subroutine reset_parameter_scan_config(self)
    class(parameter_scan_config_type), intent(in out) :: self
    type(parameter_scan_config_type) :: empty
    select type (self)
    type is (parameter_scan_config_type)
       self = empty
    end select
  end subroutine reset_parameter_scan_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_parameter_scan_config(self)
    use mp, only: broadcast
    implicit none
    class(parameter_scan_config_type), intent(in out) :: self
    call broadcast(self%delta_t_inc)
    call broadcast(self%delta_t_init)
    call broadcast(self%inc_con)
    call broadcast(self%nstep_inc)
    call broadcast(self%nstep_init)
    call broadcast(self%par_end)
    call broadcast(self%par_inc)
    call broadcast(self%par_start)
    call broadcast(self%scan_par)
    call broadcast(self%scan_restarted)
    call broadcast(self%scan_spec)
    call broadcast(self%scan_type)
    call broadcast(self%target_par)
    call broadcast(self%target_val)

    call broadcast(self%exist)
  end subroutine broadcast_parameter_scan_config
  
  !> Gets the default name for this namelist
  function get_default_name_parameter_scan_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_parameter_scan_config
    get_default_name_parameter_scan_config = "parameter_scan_knobs"
  end function get_default_name_parameter_scan_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_parameter_scan_config()
    implicit none
    logical :: get_default_requires_index_parameter_scan_config
    get_default_requires_index_parameter_scan_config = .false.
  end function get_default_requires_index_parameter_scan_config
end module parameter_scan

