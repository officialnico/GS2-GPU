!> FIXME : Add documentation
module nonlinear_terms
  !
  ! missing factors of B_a/B(theta) in A_perp terms??
  !

  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use mp, only: mp_request_null
  implicit none

  private

  public :: init_nonlinear_terms, finish_nonlinear_terms
  public :: read_parameters, wnml_nonlinear_terms, check_nonlinear_terms
  public :: add_explicit_terms, split_nonlinear, time_add_explicit_terms_field
  public :: finish_init, reset_init, nonlin, accelerated
  public :: cfl, time_add_explicit_terms, time_add_explicit_terms_mpi
  public :: calculate_current_nl_source_and_cfl_limit
  public :: cfl_req_hand, nb_check_time_step_too_large
  public :: nonlinear_terms_testing

  public :: nonlinear_terms_config_type
  public :: set_nonlinear_terms_config
  public :: get_nonlinear_terms_config
  
  ! GS2 GPU export (gs2_export): read-only access, see gpu/spec/09_export_infra.md §4
  public :: cflx, cfly, max_vel_components, istep_last, initializing, nonlinear_mode_switch
  public :: use_cfl_limit, use_2d_cfl, error_target, istep_error_start, use_order_based_error
  public :: include_apar, include_bpar, include_phi, zip, nl_forbid_force_zero

  !> Public type allowing tests to access otherwise-private routines
  type :: nonlinear_terms_testing
    contains
      procedure, nopass :: add_nl
  end type nonlinear_terms_testing

  integer :: istep_last = 0

  ! knobs
  integer :: nonlinear_mode_switch

  integer, parameter :: nonlinear_mode_none = 1, nonlinear_mode_on = 2

#ifndef SHMEM
  real, dimension (:,:), allocatable :: ba, gb, bracket
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, dimension (:,:,:), allocatable :: aba, agb, abracket
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc

#else
  real, save, dimension (:,:), pointer :: ba => null(), gb => null(), bracket => null()
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, save, dimension (:,:,:), pointer, contiguous :: aba => null(), &
       agb => null(), abracket => null()
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc
#endif

! CFL coefficients
  real :: cfl, cflx, cfly
  !> Data on maximum velocity of the nonlinear term, the reciprocal of
  !> the maximum timestep that satisfies the CFL condition. This
  !> currently has to be module level due to potential non-blocking
  !> reduction. Stores max_vel_x, max_vel_y and max_vel_error
  !> components to allow a max reduction over each element.
  real, dimension(3) :: max_vel_components

  integer :: cfl_req_hand = mp_request_null
  logical :: use_cfl_limit, use_2d_cfl

  !> Variables related to the order based error estimate -- see config documentation
  real :: error_target
  integer :: istep_error_start
  logical :: use_order_based_error

  logical :: split_nonlinear
  logical :: include_apar, include_bpar, include_phi

  real :: time_add_explicit_terms(2) = 0., time_add_explicit_terms_mpi = 0.
  real :: time_add_explicit_terms_field(2) = 0.

  logical :: nonlin = .false.
  logical :: initialized = .false.
  logical :: initializing = .true.
  logical :: alloc = .true.
  logical :: zip = .false.
  logical :: nl_forbid_force_zero = .true.
  logical :: accelerated = .false.

  logical :: exist = .false.

  !> The "large" cfl time step limit to use when we don't
  !> have a valid cfl limit from the NL term.
  real, parameter :: dt_cfl_default_large = 1.e8

  !> Used to represent the input configuration of nonlinear_terms
  type, extends(abstract_config_type) :: nonlinear_terms_config_type
     ! namelist : nonlinear_terms_knobs
     ! indexed : false
     !> Scales the estimate CFL limit used to determine when the
     !> timestep should be changed. The maximum allowed timestep
     !> satisfies `delt < cfl * min(Delta_perp/v_perp)` where `v_perp
     !> * delt` is the maximum distance travelled in a single
     !> timestep.
     real :: cfl = 0.95
     !> Set the error threshold used to determine when the timestep should change if
     !> [[nonlinear_terms_knobs:use_order_based_error]] is true.
     real :: error_target = 0.1
     !> Flag for testing. If false do not include apar contribution to nonlinear term.
     logical :: include_apar = .true.
     !> Flag for testing. If false do not include bpar contribution to nonlinear term.
     logical :: include_bpar = .true.
     !> Flag for testing. If false do not include phi contribution to nonlinear term.
     logical :: include_phi = .true.
     !> Set the first timestep for which the order based error checks are made if
     !> [[nonlinear_terms_knobs:use_order_based_error]] is true.
     integer :: istep_error_start = 30
     !> If `true` (default) then forces the nonlinear source term to
     !> zero in the forbidden region.
     logical :: nl_forbid_force_zero = .true.
     !> Determines if the nonlinear terms should be calculated. Must
     !> be one of:
     !>
     !> - 'none' - Do not include nonlinear terms, i.e. run a linear calculation.
     !> - 'default' - The same as 'none'
     !> - 'off' - The same as 'none'
     !> - 'on' - Include nonlinear terms.
     !>
     !> Could consider changing this to a logical.
     character(len = 20) :: nonlinear_mode = 'default'
     !> Do we evolve the nonlinear term separately from the linear terms (true) or
     !> include the nonlinear term as a source in the standard algorithm (false).
     logical :: split_nonlinear = .false.
     !> If true then sum x and y cfl limiting velocities instead of taking the maximum.
     logical :: use_2d_cfl = .true.
     !> If true then use the cfl limit to set the maximum timestep allowed.
     logical :: use_cfl_limit = .true.
     !> If true then use an error estimate from comparing the nonlinear source
     !> calculated using 2nd and 3rd order Adams-Bashforth schemes to control
     !> the timestep in use. This does not disable the CFL estimate.
     logical :: use_order_based_error = .false.
     !> Not currently used, should consider removing.  Original
     !> documentation was "Experts only (for secondary/tertiary
     !> calculations)."  which suggests a close relation to the
     !> `eqzip` option of [[knobs]].
     logical :: zip = .false.
   contains
     procedure, public :: read => read_nonlinear_terms_config
     procedure, public :: write => write_nonlinear_terms_config
     procedure, public :: reset => reset_nonlinear_terms_config
     procedure, public :: broadcast => broadcast_nonlinear_terms_config
     procedure, public, nopass :: get_default_name => get_default_name_nonlinear_terms_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_nonlinear_terms_config
  end type nonlinear_terms_config_type

  type(nonlinear_terms_config_type) :: nonlinear_terms_config
  
contains

  !> FIXME : Add documentation
  subroutine check_nonlinear_terms(report_unit,delt_adj)
    use gs2_time, only: code_dt_min, code_dt_max
    use kt_grids, only: box
    use run_parameters, only: nstep, wstar_units
    use theta_grid, only: nperiod
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: delt_adj
    if (nonlin) then
       write (report_unit, *) 
       write (report_unit, fmt="('This is a nonlinear simulation.')")
       write (report_unit, *) 
       if (wstar_units) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs require wstar_units = .false. in the knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if ( .not. box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
          write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (split_nonlinear) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The nonlinear term will be evolved separately from')")
          write (report_unit, fmt="('the linear terms in a crude IMEX style.')")
          write (report_unit, fmt="('This is an experimental feature.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if

       if (.not. include_apar) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The A|| potential contributions to the nonlinear term are disabled.')")
          write (report_unit, fmt="('This is probably an error.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if

       if (.not. include_bpar) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The B|| potential contributions to the nonlinear term are disabled.')")
          write (report_unit, fmt="('This is probably an error.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if

       if (.not. include_phi) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The electrostatic potential contributions to the nonlinear term are disabled.')")
          write (report_unit, fmt="('This is probably an error.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if

       if (nperiod > 1) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs usually have nperiod = 1.')") 
          write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (use_cfl_limit) then
          write (report_unit, *)
          write (report_unit, fmt="('The timestep will be adjusted to satisfy the CFL limit.')")
          write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp). (cfl)')") cfl
          write (report_unit, *)
       else
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear usually have use_cfl_limit = T')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if

       if (use_order_based_error) then
          write (report_unit, *)
          write (report_unit, fmt="('The timestep will be adjusted to keep the estimated error below',e11.4)") error_target
          write (report_unit, fmt="('beginning on step',I0)") istep_error_start

          write (report_unit, *)
       else
          if (.not. use_cfl_limit) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('Both CFL and error methods for controlling explicit timestep are disabled.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
       end if

       write (report_unit, fmt="('The minimum delt ( code_dt_min ) = ',e11.4)") code_dt_min
       write (report_unit, fmt="('The maximum delt (code_dt_max) = ',e11.4)") code_dt_max
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
    endif
  end subroutine check_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine wnml_nonlinear_terms(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    if (nonlin) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e17.10)") cfl
       write (unit, fmt="(' nl_forbid_force_zero = ',L1)") nl_forbid_force_zero
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine init_nonlinear_terms(nonlinear_terms_config_in) 
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, nx, ny, akx, aky
    use le_grids, only: init_le_grids, nlambda, negrid
    use species, only: init_species, nspec
    use gs2_layouts, only: init_dist_fn_layouts, yxf_lo, accelx_lo
    use gs2_layouts, only: init_gs2_layouts
    use gs2_transforms, only: init_transforms
    use mp, only: nproc, iproc
    use array_utils, only: zero_array
#ifdef SHMEM
    use shm_mpi3, only : shm_alloc
#endif
    implicit none
    type(nonlinear_terms_config_type), intent(in), optional :: nonlinear_terms_config_in
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    if (debug) write(6,*) "init_nonlinear_terms: init_gs2_layouts"
    call init_gs2_layouts
    if (debug) write(6,*) "init_nonlinear_terms: init_theta_grid"
    call init_theta_grid
    if (debug) write(6,*) "init_nonlinear_terms: init_kt_grids"
    call init_kt_grids
    if (debug) write(6,*) "init_nonlinear_terms: init_le_grids"
    call init_le_grids
    if (debug) write(6,*) "init_nonlinear_terms: init_species"
    call init_species
    if (debug) write(6,*) "init_nonlinear_terms: init_dist_fn_layouts"
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nproc, iproc)
    
    call read_parameters(nonlinear_terms_config_in)

    ! Initialise maximum velocity values
    max_vel_components = 0

    if (nonlin) then
       if (debug) write(6,*) "init_nonlinear_terms: init_transforms"
       call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

       if (debug) write(6,*) "init_nonlinear_terms: allocations"
       if (alloc) then
          if (accelerated) then
#ifndef SHMEM
             allocate (     aba(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (     agb(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (abracket(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
#else
             call shm_alloc(aba, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
             call shm_alloc(agb, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
             call shm_alloc(abracket, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
#endif
             call zero_array(aba) ; call zero_array(agb) ; call zero_array(abracket)
          else
#ifndef SHMEM
             allocate (     ba(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (     gb(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (bracket(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
#else
             call shm_alloc(ba, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
             call shm_alloc(gb, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
             call shm_alloc(bracket, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
#endif
             call zero_array(ba) ; call zero_array(gb) ; call zero_array(bracket)
          end if
          alloc = .false.
       end if

       cfly = maxval(aky) * 0.5 / cfl
       cflx = maxval(akx) * 0.5 / cfl
    end if

  end subroutine init_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine read_parameters(nonlinear_terms_config_in)
    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    implicit none
    type(nonlinear_terms_config_type), intent(in), optional :: nonlinear_terms_config_in
    type (text_option), dimension (4), parameter :: nonlinearopts = &
         [ text_option('default', nonlinear_mode_none), &
            text_option('none', nonlinear_mode_none), &
            text_option('off', nonlinear_mode_none), &
            text_option('on', nonlinear_mode_on) ]
    character(len = 20) :: nonlinear_mode
    integer :: ierr

    if (present(nonlinear_terms_config_in)) nonlinear_terms_config = nonlinear_terms_config_in

    call nonlinear_terms_config%init(name = 'nonlinear_terms_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    cfl = nonlinear_terms_config%cfl
    error_target = nonlinear_terms_config%error_target
    include_apar = nonlinear_terms_config%include_apar
    include_bpar = nonlinear_terms_config%include_bpar
    include_phi = nonlinear_terms_config%include_phi
    istep_error_start = nonlinear_terms_config%istep_error_start
    nl_forbid_force_zero = nonlinear_terms_config%nl_forbid_force_zero
    nonlinear_mode = nonlinear_terms_config%nonlinear_mode
    split_nonlinear = nonlinear_terms_config%split_nonlinear
    use_2d_cfl = nonlinear_terms_config%use_2d_cfl
    use_cfl_limit = nonlinear_terms_config%use_cfl_limit
    use_order_based_error = nonlinear_terms_config%use_order_based_error
    zip = nonlinear_terms_config%zip

    exist = nonlinear_terms_config%exist

    ierr = error_unit()
    call get_option_value &
         (nonlinear_mode, nonlinearopts, nonlinear_mode_switch, &
         ierr, "nonlinear_mode in nonlinear_terms_knobs",.true.)

    nonlin = nonlinear_mode_switch == nonlinear_mode_on
  end subroutine read_parameters

  !> FIXME : Add documentation  
  subroutine add_explicit_terms (g1, g2, g3, phi, apar, bpar, istep, bd)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: save_dt_cfl
    use job_manage, only: time_message
    use mp, only: get_mp_times
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    real :: dt_cfl
    logical, parameter :: nl = .true.
    real :: mp_total_after, mp_total
    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total)

    if (nonlin) then
       if (istep /= 0) then
          call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, nl)
       end if
    else
       dt_cfl = dt_cfl_default_large
       call save_dt_cfl (dt_cfl)
#ifdef LOWFLOW
       if (istep /=0) &
            call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd)
#endif
    end if
    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total_after)
    time_add_explicit_terms_mpi = time_add_explicit_terms_mpi + (mp_total_after - mp_total)
  end subroutine add_explicit_terms

  !> FIXME : Add documentation  
  subroutine add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd,  nl)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use dist_fn_arrays, only: g
    use gs2_time, only: save_dt_cfl, code_dt, check_time_step_too_large
    use run_parameters, only: reset, immediate_reset
    use unit_tests, only: debug_message
    use mp, only: nb_max_allreduce, max_allreduce
    use array_utils, only: copy, zero_array
    use gs2_export_probe, only: probe_active, probe_g, probe_real
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    logical, intent (in), optional :: nl
    integer, parameter :: verb = 4
    integer :: iglo, ik
    real :: dt_cfl

    real :: error
    real :: target_dt

    call debug_message(verb, 'nonlinear_terms::add_explicit starting')

    if (initializing) then
       if (present(nl)) then
          dt_cfl = dt_cfl_default_large
          call save_dt_cfl (dt_cfl)

          ! Ensure the module level max_vel is consistent with
          ! the selected cfl limiting timestep.
          max_vel_components = 1.0/dt_cfl
       end if
       return
    endif

    !
    ! @todo Currently not self-starting.  Need to fix this.
    !

    call debug_message(verb, 'nonlinear_terms::add_explicit copying old terms')

    ! The following if statement ensures that the nonlinear term is calculated only once
    ! per timestep. It guards against calculating the nonlinear term in two cases:
    !   1. the second call to `timeadv` in a usual timestep
    !   2. both calls to `timeadv` after the CFL condition check has triggered a reset
    !      of this timestep (when `immediate_reset = .true.`). As the nonlinear term is
    !      independent of `delt`, it does not need to be calculated again.
    !
    ! N.B. When `immediate_reset = .false.`, the timestep is completed even though the
    ! CFL condition is broken. `delt` and the response matrix are changed at the end of
    ! the timestep; the following timestep is "normal", and the nonlinear term is
    ! calculated.
    if (istep /= istep_last) then

       istep_last = istep

       ! Shuffle time history along -- this is a place pointers might help
       ! avoid some copying.
       ! Should we only need to do this if we have a nonlinear simulation?
       ! It's possible we might have other explicit sources, in which case we
       ! probably always need to do this shuffling, however in that case this
       ! subroutine probably best belongs in a different module.
       call copy(g2, g3)
       call copy(g1, g2)

       call debug_message(verb, 'nonlinear_terms::add_explicit copied old terms')

       ! if running nonlinearly, then compute the nonlinear term at grid points
       ! and store it in g1
       if (present(nl)) then
          call debug_message(verb, 'nonlinear_terms::add_explicit calling add_nl')

          ! If we're not using various limits then we'd like to zero
          ! out the maximum velocity however, we later take 1/max_vel
          ! so let's just set it small.
          max_vel_components = epsilon(0.0)

          call add_nl (g, g1, phi, apar, bpar)
          if (probe_active) call probe_g('nl_before_center', g1)

          ! takes g1 at grid points and returns 2*g1 at cell centers
          call nl_center (g1, bd)
          if (probe_active) call probe_g('nl_after_center', g1)

          ! Currently the module level max_vel_components contains the
          ! maximum x and y velocities found from the nonlinear term
          ! on this processor. This can be used to find the cfl
          ! limting time step.

          ! Reset the x and y components of max_vel_components if we
          ! don't want to consider the cfl limit
          if (.not. use_cfl_limit) then
             max_vel_components(1:2) = epsilon(0.0)
          end if

          ! Calculate the effective maximum velocity based on the error
          ! based time step limit (really we estimate the maximum time step
          ! to keep our errors below our target and then convert to an
          ! effective velocity).
          if (use_order_based_error .and. (istep > istep_error_start)) then
             ! Get the current order based error estimate
             error = estimate_error(g1, g2, g3)

             ! Based on the current error estimate and the expected error scaling
             ! work out what time step we expect to give us the target error.
             ! Note that we don't scale by cfl here as the user has direct control
             ! over error_target already. Also note that by taking the square root
             ! we are assuming that the error looks like dt^2. We should probably
             ! look at what order schemes we're actually comparing to determine
             ! the appropriate dependency.
             target_dt =code_dt*sqrt(error_target/(error))

             ! Now we have an error limited time step lets convert this
             ! to an effective maximum velocity
             max_vel_components(3) = 1.0/target_dt
          end if

          ! Now we can reduce our local limiting velocities to get
          ! a global value for the maximum velocity. We only need to
          ! reduce if at least one time limit check is active.
          if (use_cfl_limit .or. use_order_based_error) then
             ! Estimate the global cfl limit based on max_vel calculated
             ! in `add_nl`.If immediate_reset is true, check the max
             ! allreduce immediately.  If not, overlap comms and work,
             ! and check for reset at the end of the timestep.
             if(immediate_reset)then
                call max_allreduce(max_vel_components)
                ! This call can update the value of reset.
                call set_cfl_time_step_and_check_if_too_large
                if (probe_active) call probe_real('nl_cfl', max_vel_components)
             else
                ! Do a nonblocking reduce and check later in
                ! [[gs2_reinit:check_time_step]].
                call nb_max_allreduce(max_vel_components, cfl_req_hand)
             end if
          else
             ! If we don't have any checks active then we better make sure that the
             ! limiting time step is set. We will have set max_vel earlier to
             ! epsilon(0.0) so we can just use the usual routine to set this, but
             ! can skip communications. Really we probably just need to do this once
             ! but for now just put this here to be on the safe side. The cost should
             ! be low.
             call set_cfl_time_step_and_check_if_too_large
             if (probe_active) call probe_real('nl_cfl', max_vel_components)
          end if

          if(reset) then
             return !Return if resetting
          endif

       else
          call zero_array(g1)
       end if

#ifdef LOWFLOW
       ! do something
#endif

    end if

    if (zip) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          if (ik /= 1) then
             g (:,1,iglo) = 0.
             g (:,2,iglo) = 0.
             g1(:,1,iglo) = 0.
             g1(:,2,iglo) = 0.
          end if
       end do
    end if

    istep_last = istep
  end subroutine add_explicit

  !> Get an estimate of the error on the NL source term by comparing
  !> 2nd and 3rd order calculations of the source term.
  !>
  !> Here we're calculating g at t+code_dt using g(t) = g1, g(t-dt1) = g2
  !> and g(t-dt2) = g3 using 2nd and 3rd order methods and comparing
  !> the difference to estimate the error.
  real function estimate_error(g1, g2, g3) result(error)
    use gs2_time, only: get_adams_bashforth_coefficients, time_history_available
    use mp, only: max_allreduce
    use array_utils, only: gs2_max_abs
    implicit none
    complex, dimension(:, :, :), intent(in) :: g1, g2, g3
    complex, dimension(:, :, :), allocatable :: source_order_3, source_order_2
    real, dimension(2) :: intermediate_error
    real, dimension(:), allocatable :: coeff_3rd, coeff_2nd

    allocate(coeff_3rd(min(time_history_available(), 3)))
    allocate(coeff_2nd(min(time_history_available(), 2)))

    ! Check that we actually have 3rd and 2nd order coefficients.
    ! This could fail on the first two time steps of a simulation, for
    ! example, where we don't have enough time history to allow a
    ! 2nd/3rd order calculation. For now just return an error of
    ! error_target times cfl squared. This is chosen such the target
    ! timestep remains the current timestep.
    if ( (size(coeff_3rd) /= 3) .or. (size(coeff_2nd) /= 2) ) then
       error = error_target *cfl*cfl
       return
    end if

    ! Calculate the coefficients required for variable time step schemes
    coeff_3rd = get_adams_bashforth_coefficients(maximum_order = 3)
    coeff_2nd = get_adams_bashforth_coefficients(maximum_order = 2)

    ! Calculate the source terms using 3rd and 2nd order methods
    allocate(source_order_3, mold = g1)
    allocate(source_order_2, mold = g1)
    source_order_3 = coeff_3rd(1) * g1 + coeff_3rd(2) * g2 + coeff_3rd(3) * g3
    source_order_2 = coeff_2nd(1) * g1 + coeff_2nd(2) * g2

    ! Here we find the maximum absolute error
    intermediate_error(1) = gs2_max_abs(source_order_3 - source_order_2)

    ! Now we find the maximum of source_order_2 to scale the absolute
    ! error with.  Note we use this rather than source_order_3 as the
    ! higher order method has a smaller region of stability so is
    ! perhaps more likely to blow up.
    intermediate_error(2) = gs2_max_abs(source_order_2)

    ! Reduce over all processors
    call max_allreduce(intermediate_error)

    ! Now we form a crude relative error by dividing by the maximum
    ! source value. We do this rather than forming the true relative
    ! error (i.e. dividing the entire absolute difference by the source)
    ! to avoid issues with dividing by small numbers. This will tend to
    ! underestimate the real relative error.
    error = intermediate_error(1) / intermediate_error(2)
  end function estimate_error

  !> Takes input array evaluated at theta grid points
  !! and overwrites it with array evaluated at cell centers
  !! note that there is an extra factor of 2 in output array
  subroutine nl_center (gtmp, bd)
    !
    !CMR, 13/10/2014
    ! Fixing some (cosmetic) issues from prior to R3021.
    ! (1) forbid(-ntgrid:ntgrid) refers to grid points at cell boundaries
    !     gtmp(-ntgrid:ntgrid-1) refers to cell centers
    !     => applying forbid to output gtmp did not zero the correct things!!!
    ! (2) totally trapped particles are special, so return source without upwinding is fine
    !     BUT source for other trapped particles should NOT enter forbidden region.
    !     if ig or ig+1 forbidden (i.e. ig+1 or ig is bounce point) now return gtmp(ig,:,iglo)=0
    !
    use gs2_layouts, only: g_lo, il_idx
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, grid_has_trapped_particles, jend, is_ttp
    use mp, only: mp_abort

    implicit none
    real, intent (in) :: bd    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: gtmp
    integer :: iglo, il, ig
    logical :: trapped

    trapped = grid_has_trapped_particles()

    ! factor of one-half appears elsewhere
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ig) &
    !$OMP SHARED(g_lo, nl_forbid_force_zero, forbid, gtmp, ntgrid, is_ttp, jend, trapped, bd) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)

       !
       !CMR, 7/10/2014: 
       ! Incoming gtmp SHOULD vanish in forbidden region, 
       ! New logical in nonlinear_terms_knobs namelist: nl_forbid_force_zero
       !     nl_forbid_force_zero =.t. : force zeroing    (default)
       !     nl_forbid_force_zero =.f. : NO forced zeroing
       !                                 ie assume forbidden gtmp is zero on entry 
       !

       if ( nl_forbid_force_zero ) then
          ! force spurious gtmp outside trapped boundary to be zero
          where (forbid(:,il))
             gtmp(:,1,iglo) = 0.0
             gtmp(:,2,iglo) = 0.0
          end where
       endif

       do ig = -ntgrid, ntgrid-1
          !
          !CMR, 7/10/2014: 
          ! loop sets gtmp to value at upwinded cell center RIGHT of theta(ig)
          !           except for ttp where upwinding makes no sense!
          !
          ! Note this cycle doesn't skip when ig is in the forbidden region
          ! even if the pitch angle is greater than the totally trapped one
          ! at this location. This ensures that whilst we leave the (allowed)
          ! totally trapped source alone we still zero the NL source for forbidden
          ! points. Previously we cycled simply if the pitch angle index was larger
          ! than the totally trapped one. When there are multiple totally trapped
          ! pitch angles this may not be equivalent.
          if (is_ttp(ig, il)) cycle
          ! The below might be better off using forbid instead of checking jend
          if ( trapped .and. ( il > jend(ig) .or. il > jend(ig+1)) ) then
             !
             !CMR, 7/10/2014: 
             !   if either ig or ig+1 is forbidden, no source possible in a cell RIGHT of theta(ig) 
             !   => gtmp(ig,1:2,iglo)=0
             !
             gtmp(ig,1:2,iglo) = 0.0
          else 
             !
             !CMR, 7/10/2014: 
             !    otherwise ig and ig+1 BOTH allowed, and upwinding in cell RIGHT of theta(ig) is fine
             !
             gtmp(ig,1,iglo) = (1.+bd)*gtmp(ig+1,1,iglo) + (1.-bd)*gtmp(ig,1,iglo)
             gtmp(ig,2,iglo) = (1.-bd)*gtmp(ig+1,2,iglo) + (1.+bd)*gtmp(ig,2,iglo)

          endif
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine nl_center

  !> Calculates the current nonlinear source term by calling add_nl
  !> and the associated cfl limit.
  subroutine calculate_current_nl_source_and_cfl_limit(g_in, g1, phi, apar, bpar, &
       max_vel_local, need_to_adjust, calculate_cfl_limit)
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension (:,:,:), intent (in) :: g_in
    complex, dimension (:,:,:), intent (out) :: g1
    complex, dimension (:,:,:), intent (in) :: phi, apar, bpar
    real, intent(out) :: max_vel_local
    logical, intent(in) :: need_to_adjust
    logical, intent(in), optional :: calculate_cfl_limit

    use_cfl_limit = get_option_with_default(calculate_cfl_limit, .true.)
    call add_nl(g_in, g1, phi, apar, bpar, need_to_adjust)

    if (use_cfl_limit) then
       ! This code may need changing if the handling of the cfl limit changes
       max_vel_local = get_max_vel(max_vel_components)
    else
       max_vel_local = 1.0e-8
    end if

  end subroutine calculate_current_nl_source_and_cfl_limit

  !> Calculate the nonlinear term and part of the corresponding CFL condition
  subroutine add_nl (g_in, g1, phi, apar, bpar, adjust)
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: g_adjust, from_g_gs2
    use gs2_transforms, only: transform2, inverse2
    use kt_grids, only: aky, akx
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use constants, only: zi
    use unit_tests, only: debug_message
    use optionals, only: get_option_with_default
    use array_utils, only: copy, gs2_max_abs
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1
    !> Local array used to store the adjusted input dfn if needed. Saves
    !> calling g_adjust twice.
    complex, dimension (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc) :: g2
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    logical, intent(in), optional :: adjust
    logical :: should_adjust
    integer, parameter :: verb = 4
    integer :: iglo, ik, it
    real :: max_vel_x, max_vel_y
    
    call debug_message(verb, 'nonlinear_terms::add_nl starting')
    should_adjust = get_option_with_default(adjust, .true.)

    !Form g1=i*kx*chi
    call load_kx_phi(g1, phi)    
    call load_kx_bpar(g1, bpar)    
    call load_kx_apar(g1, apar)    

    call debug_message(verb, 'nonlinear_terms::add_nl calling transform2')

    !Transform to real space
    if (accelerated) then
       call transform2 (g1, aba)
       max_vel_y = gs2_max_abs(aba)
    else
       call transform2 (g1, ba)
       max_vel_y = gs2_max_abs(ba)
    end if
    max_vel_y = max_vel_y * cfly

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx')
    
    !Form g1=i*ky*g_wesson
    call copy(g_in, g1)
    if (should_adjust) call g_adjust(g1, phi, bpar, direction = from_g_gs2)
    call copy(g1, g2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik) &
    !$OMP SHARED(g_lo, g1, aky) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*aky(ik)
    enddo
    !$OMP END PARALLEL DO

    !Transform to real space
    if (accelerated) then
       call transform2 (g1, agb)
    else
       call transform2 (g1, gb)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dg/dy')

    call calculate_bracket(.true.)
    
    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx.dg/dy')

    !Form g1=i*ky*chi
    call load_ky_phi(g1, phi)
    call load_ky_bpar(g1, bpar)
    call load_ky_apar(g1, apar)

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dy')

    !Transform to real space    
    if (accelerated) then
       call transform2 (g1, aba)
       max_vel_x = gs2_max_abs(aba)
    else
       call transform2 (g1, ba)
       max_vel_x = gs2_max_abs(ba)
    end if
    max_vel_x = max_vel_x * cflx

    !Form g1=i*kx*g_wesson noting that g2 holds our starting g_wesson

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it) &
    !$OMP SHARED(g_lo, g2, akx) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       g2(:,:,iglo) = g2(:,:,iglo) * zi * akx(it)
    enddo
    !$OMP END PARALLEL DO

    !Transform to real space
    if (accelerated) then
       call transform2 (g2, agb)
    else
       call transform2 (g2, gb)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dg/dx')

    call calculate_bracket(.false.)
    
    call debug_message(verb, &
      'nonlinear_terms::add_nl calculated dphi/dx.dg/dy.dphi/dy.dg/dx')

    !Transform NL source back to spectral space
    if (accelerated) then
       call inverse2 (abracket, g1)
    else
       call inverse2 (bracket, g1)
    end if

    !Store local max vel components for future max reduce
    max_vel_components(1) = max_vel_x
    max_vel_components(2) = max_vel_y

    ! Scale by kxfac
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo) &
    !$OMP SHARED(g_lo, g1, kxfac) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g1(:, :, iglo) = g1(:, :, iglo) * kxfac
    end do
    !$OMP END PARALLEL DO

    call debug_message(verb, &
      'nonlinear_terms::add_nl calculated nonlinear term')

  end subroutine add_nl

  !> FIXME : Add documentation    
  subroutine load_kx_phi(g1, phi)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: aj0
    use run_parameters, only: has_phi, fphi
    use kt_grids, only: akx
    use constants, only: zi
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi      
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer :: iglo, it, ik

    if (has_phi .and. include_phi) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, it, ik, fac) &
       !$OMP SHARED(g_lo, g1, akx, aj0, phi, fphi) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          fac = zi*akx(it)*aj0(:,iglo)*phi(:,it,ik)*fphi
          g1(:,1,iglo) = fac
          g1(:,2,iglo) = fac
       end do
       !$OMP END PARALLEL DO
    else
       call zero_array(g1)
    endif
  end subroutine load_kx_phi

  !> FIXME : Add documentation    
  subroutine load_ky_phi(g1, phi)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: aj0
    use run_parameters, only: has_phi, fphi
    use kt_grids, only: aky
    use constants, only: zi
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi      
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer :: iglo, it, ik

    if (has_phi .and. include_phi) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, it, ik, fac) &
       !$OMP SHARED(g_lo, g1, aky, aj0, phi, fphi) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          fac = zi*aky(ik)*aj0(:,iglo)*phi(:,it,ik)*fphi
          g1(:,1,iglo) = fac
          g1(:,2,iglo) = fac
       end do
       !$OMP END PARALLEL DO
    else
       call zero_array(g1)
    endif

  end subroutine load_ky_phi

  ! should I use vpa or vpac in next two routines??

  !> FIXME : Add documentation    
  subroutine load_kx_apar(g1, apar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vpa, aj0
    use species, only: spec
    use run_parameters, only: has_apar, fapar
    use kt_grids, only: akx
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: apar
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if(.not. (include_apar .and. has_apar) ) return

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, fac) &
    !$OMP SHARED(g_lo, g1, vpa, akx, aj0, spec, apar, fapar) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       fac = zi*akx(it)*aj0(:,iglo)*spec(is)%stm*apar(:,it,ik)*fapar
       g1(:,1,iglo) = g1(:,1,iglo) - fac*vpa(:,1,iglo)
       g1(:,2,iglo) = g1(:,2,iglo) - fac*vpa(:,2,iglo)
    end do
    !$OMP END PARALLEL DO
  end subroutine load_kx_apar

  !> FIXME : Add documentation    
  subroutine load_ky_apar(g1, apar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vpa, aj0
    use species, only: spec
    use run_parameters, only: has_apar, fapar
    use kt_grids, only: aky
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: apar            
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if(.not. (include_apar .and. has_apar)) return

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, fac) &
    !$OMP SHARED(g_lo, g1, vpa, aky, aj0, spec, apar, fapar) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       fac = zi*aky(ik)*aj0(:,iglo)*spec(is)%stm*apar(:,it,ik)*fapar
       g1(:,1,iglo) = g1(:,1,iglo) - fac*vpa(:,1,iglo)
       g1(:,2,iglo) = g1(:,2,iglo) - fac*vpa(:,2,iglo)
    end do
    !$OMP END PARALLEL DO
  end subroutine load_ky_apar

  !> FIXME : Add documentation    
  subroutine load_kx_bpar(g1, bpar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vperp2, aj1
    use species, only: spec
    use run_parameters, only: has_bpar, fbpar
    use kt_grids, only: akx
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: bpar      
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if (.not. (include_bpar .and. has_bpar)) return

    ! Is this factor of two from the old normalization?

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, fac) &
    !$OMP SHARED(g_lo, g1, akx, aj1, vperp2, spec, bpar, fbpar) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       fac = zi*akx(it)*aj1(:,iglo) &
            *2.0*vperp2(:,iglo)*spec(is)%tz*bpar(:,it,ik)*fbpar
       g1(:,1,iglo) = g1(:,1,iglo) + fac
       g1(:,2,iglo) = g1(:,2,iglo) + fac
    end do
    !$OMP END PARALLEL DO
  end subroutine load_kx_bpar

  !> FIXME : Add documentation    
  subroutine load_ky_bpar(g1, bpar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vperp2, aj1
    use species, only: spec
    use run_parameters, only: has_bpar, fbpar
    use kt_grids, only: aky
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: bpar
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer iglo, it, ik, is

    if (.not. (include_bpar .and. has_bpar)) return

    ! Is this factor of two from the old normalization?

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, fac) &
    !$OMP SHARED(g_lo, g1, aky, aj1, vperp2, spec, bpar, fbpar) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       fac =  zi*aky(ik)*aj1(:,iglo) &
            *2.0*vperp2(:,iglo)*spec(is)%tz*bpar(:,it,ik)*fbpar
       g1(:,1,iglo) = g1(:,1,iglo) + fac
       g1(:,2,iglo) = g1(:,2,iglo) + fac
    end do
    !$OMP END PARALLEL DO
  end subroutine load_ky_bpar

  !> Calculate (d Chi /dx).(d g_wesson/dy) and store in bracket if is_first_term = .true.
  !> else calculate (d Chi /dy).(d g_wesson/dx) and subtract from bracket
  !>
  !> @note The final bracket value must be scaled by kxfac to produce
  !> the correct value to be used.
  subroutine calculate_bracket(is_first_term)
    use gs2_layouts, only: accelx_lo, yxf_lo
    implicit none
    logical :: is_first_term
    integer :: iyxf
    if (is_first_term) then
       if (accelerated) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iyxf) &
          !$OMP SHARED(accelx_lo, agb, aba, abracket) &
          !$OMP SCHEDULE(static)
          do iyxf = accelx_lo%llim_proc, accelx_lo%ulim_proc
             abracket(:, :, iyxf) = aba(:, :, iyxf)*agb(:, :, iyxf)
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iyxf) &
          !$OMP SHARED(yxf_lo, bracket, ba, gb) &
          !$OMP SCHEDULE(static)
          do iyxf = yxf_lo%llim_proc, yxf_lo%ulim_proc
             bracket(:, iyxf) = ba(:, iyxf)*gb(:, iyxf)
          end do
          !$OMP END PARALLEL DO
       endif
    else
       if (accelerated) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iyxf) &
          !$OMP SHARED(accelx_lo, agb, aba, abracket) &
          !$OMP SCHEDULE(static)
          do iyxf = accelx_lo%llim_proc, accelx_lo%ulim_proc
             abracket(:, :, iyxf) = abracket(:, :, iyxf) - aba(:, :, iyxf)*agb(:, :, iyxf)
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iyxf) &
          !$OMP SHARED(yxf_lo, bracket, ba, gb) &
          !$OMP SCHEDULE(static)
          do iyxf = yxf_lo%llim_proc, yxf_lo%ulim_proc
             bracket(:, iyxf) = bracket(:, iyxf) - ba(:, iyxf)*gb(:, iyxf)
          end do
          !$OMP END PARALLEL DO
       endif
    endif
  end subroutine calculate_bracket

  !> Finish nonblocking calculation of the maximum velocity
  !> and set the cfl limit.
  subroutine nb_check_time_step_too_large
    use run_parameters, only: immediate_reset
    implicit none
    if (immediate_reset) return
    call finish_nonblocking_max_vel_reduction
    call set_cfl_time_step_and_check_if_too_large
  end subroutine nb_check_time_step_too_large

  !> Finish nonblocking calculation of the maximum velocity
  !>
  !> This receives the mpi_iallreduce of the maximum velocity.
  subroutine finish_nonblocking_max_vel_reduction
    use mp, only: wait, mp_request_null
    implicit none

    ! Wait for the nonblocking maximum operation used to find
    ! the maximum velocity.
    if (cfl_req_hand == mp_request_null) return
    call wait(cfl_req_hand)
    cfl_req_hand = mp_request_null
  end subroutine finish_nonblocking_max_vel_reduction

  !> Given the max_vel_components array return the limiting max_vel
  real pure function get_max_vel(components) result(max_vel)
    real, dimension(3), intent(in) :: components
    real :: max_vel_cfl, max_vel_error
    if (use_2d_cfl) then
       max_vel_cfl = components(1) + components(2)
    else
       max_vel_cfl = maxval(components(1:2))
    end if
    max_vel_error = components(3)
    max_vel = max(max_vel_cfl, max_vel_error)
  end function get_max_vel

  !> Calculates the cfl limit from the module level max_vel value
  !> saves it using [[gs2_time:save_dt_cfl]] and then sets the reset
  !> flag based on checking if the current time step is too large.
  subroutine set_cfl_time_step_and_check_if_too_large
    use run_parameters, only: reset
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    implicit none
    real :: dt_cfl, max_vel_cfl, max_vel_error, max_vel

    max_vel = get_max_vel(max_vel_components)

    ! Estimate the global cfl limit based on max_vel.
    if (max_vel <= 0.0) then
       ! We should only end up here in unusual cases (likely
       ! tests) as we initialise max_vel > 0 so the only way
       ! we get here is if the estimated NL velocity limit is
       ! exactly 0.
       dt_cfl = dt_cfl_default_large
    else
       dt_cfl = 1./max_vel
    end if

    call save_dt_cfl (dt_cfl)

    ! Now check to see if we've violated the cfl condition.
    call check_time_step_too_large(reset)
  end subroutine set_cfl_time_step_and_check_if_too_large

  !> FIXME : Add documentation
  subroutine reset_init
    implicit none
    initialized = .false.
    initializing = .true.
  end subroutine reset_init

  !> FIXME : Add documentation  
  subroutine finish_init
    implicit none
    initializing = .false.
  end subroutine finish_init

  !> FIXME : Add documentation  
  subroutine finish_nonlinear_terms
    use gs2_transforms, only: finish_transforms
#ifdef SHMEM
    use shm_mpi3, only : shm_free
#endif
    implicit none

#ifndef SHMEM
    if (allocated(aba)) deallocate (aba, agb, abracket)
    if (allocated(ba)) deallocate (ba, gb, bracket)
#else
    if (associated(aba)) call shm_free(aba)
    if (associated(agb)) call shm_free(agb)
    if (associated(abracket)) call shm_free(abracket)
    if (associated(ba)) then
       call shm_free(ba)
    endif
    if (associated(gb)) then
       call shm_free(gb)
    endif
    if (associated(bracket)) then
       call shm_free(bracket)
    endif
#endif
    nonlin = .false. ; alloc = .true. ; zip = .false. ; accelerated = .false.
    initialized = .false. ; initializing = .true. ; exist = .false.

    call finish_transforms
    call nonlinear_terms_config%reset()
  end subroutine finish_nonlinear_terms

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_nonlinear_terms_config(nonlinear_terms_config_in)
    use mp, only: mp_abort
    type(nonlinear_terms_config_type), intent(in), optional :: nonlinear_terms_config_in
    if (initialized) then
       call mp_abort("Trying to set nonlinear_terms_config when already initialized.", to_screen = .true.)
    end if
    if (present(nonlinear_terms_config_in)) then
       nonlinear_terms_config = nonlinear_terms_config_in
    end if
  end subroutine set_nonlinear_terms_config

#include "nonlinear_terms_auto_gen.inc"
end module nonlinear_terms
