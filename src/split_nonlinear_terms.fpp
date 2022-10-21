!> A module to deal with advancing the nonlinear term
!> separately from the linear terms.
module split_nonlinear_terms
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use rk_schemes, only: rk_scheme
  implicit none
  private

  public :: init_split_nonlinear_terms, finish_split_nonlinear_terms
  public :: advance_nonlinear_term

  logical :: initialized = .false.

  !> Interface to describe the call signatures for routines to calculate the
  !> field from gnew (invert_field_func) and to calculate the source term used
  !> in the advance_nonlinear_term methods. Used to allow us to pass in these
  !> methods, giving us a way to test the integrators on different problems.
  interface
     subroutine invert_field_func (phi, apar, bpar, gf_lo)
       use theta_grid, only: ntgrid
       complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
       logical, optional :: gf_lo
     end subroutine invert_field_func

     subroutine source_term_func(g1, phi, apar, bpar, dt_lim, calculate_cfl_limit)
       complex, dimension (:,:,:), intent (out) :: g1
       complex, dimension (:,:,:), intent (in) :: phi, apar, bpar
       real, intent(out) :: dt_lim
       logical, intent(in), optional :: calculate_cfl_limit
     end subroutine source_term_func
  end interface

  ! Related to inputs
  integer :: split_method_switch
  integer, parameter :: split_method_ab3 = 1, split_method_backwards_euler = 2, &
       split_method_rk = 3
  type(rk_scheme) :: the_rk_scheme

  real :: convergence_tolerance
  real :: error_tolerance

  !> Used to represent the input configuration of split_nonlinear_terms
  type, extends(abstract_config_type) :: split_nonlinear_terms_config_type
     ! namelist : split_nonlinear_terms_knobs
     ! indexed : false
     !> The tolerance used in convergence checks. Currently only used
     !> for halting the fixed-point iteration in the beuler method
     real :: convergence_tolerance = 1.0e-2
     !> The tolerance used in error controlled schemes. Attempt to adjust
     !> time step to keep error below this tolerance.
     real :: error_tolerance = 1.0e-4
     !> Choose which specific RK scheme to use in the RK method.
     !> Generally if error_tolerance is small a higher order scheme
     !> will be more effective than low order. At loose tolerances
     !> low order schemes are generally more efficient
     !> Valid options include:
     !>
     !>  - 'cashkarp' -- Cash-Karp scheme of order 4(5)
     !>  - 'default' -- Same as heun
     !>  - 'heun' -- Heun scheme of order 1(2)
     character(len = 20) :: rk_method = 'default'
     !> What algorithm do we use to advance the nonlinear term if
     !> split_nonlinear is true.
     !> Valid options include:
     !>
     !>  - 'AB3' -- Adams-Bashforth (up to) 3rd order. CLF controlled time step.
     !>  - 'beuler' -- Backwards Euler. Relative error controlled time step
     !>  - 'default' -- The same as 'AB3'.
     !>  - 'RK' -- RK scheme with embedded error estimate. Exact scheme can be switched out
     !>            using rk_method switch. Absolute error controlled time step
     !>
     character(len = 20) :: split_method = 'default'
   contains
     procedure, public :: read => read_split_nonlinear_terms_config
     procedure, public :: write => write_split_nonlinear_terms_config
     procedure, public :: reset => reset_split_nonlinear_terms_config
     procedure, public :: broadcast => broadcast_split_nonlinear_terms_config
     procedure, public, nopass :: get_default_name => get_default_name_split_nonlinear_terms_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_split_nonlinear_terms_config
  end type split_nonlinear_terms_config_type

  type(split_nonlinear_terms_config_type) :: split_nonlinear_terms_config

contains

  !> Initialises the split_nonlinear_terms module. Primarily just reading input
  subroutine init_split_nonlinear_terms(split_nonlinear_terms_config_in)
    implicit none
    type(split_nonlinear_terms_config_type), intent(in), optional :: split_nonlinear_terms_config_in

    if (initialized) return
    initialized = .true.

    call read_parameters(split_nonlinear_terms_config_in)

  end subroutine init_split_nonlinear_terms

  !> Read the split_nonlinear_terms namelist
  subroutine read_parameters(split_nonlinear_terms_config_in)
    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    use rk_schemes, only: get_rk_schemes_as_text_options, get_rk_scheme_by_id
    implicit none
    type(split_nonlinear_terms_config_type), intent(in), optional :: split_nonlinear_terms_config_in
    type(text_option), dimension (*), parameter :: split_method_opts = &
         [ &
         text_option('default', split_method_ab3), &
         text_option('AB3', split_method_ab3), &
         text_option('RK', split_method_rk), &
         text_option('beuler', split_method_backwards_euler)  &
         ]
    character(len = 20) :: rk_method, split_method
    logical :: exist
    integer :: ierr, rk_method_switch

    if (present(split_nonlinear_terms_config_in)) split_nonlinear_terms_config = split_nonlinear_terms_config_in

    call split_nonlinear_terms_config%init(name = 'split_nonlinear_terms_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    convergence_tolerance = split_nonlinear_terms_config%convergence_tolerance
    error_tolerance = split_nonlinear_terms_config%error_tolerance
    rk_method = split_nonlinear_terms_config%rk_method
    split_method = split_nonlinear_terms_config%split_method

    exist = split_nonlinear_terms_config%exist

    ierr = error_unit()

    call get_option_value &
         (split_method, split_method_opts, split_method_switch, &
         ierr, "split_method in split_nonlinear_terms_knobs",.true.)

    call get_option_value &
         (rk_method, get_rk_schemes_as_text_options(), rk_method_switch, &
         ierr, "rk_method in split_nonlinear_terms_knobs",.true.)
    the_rk_scheme = get_rk_scheme_by_id(rk_method_switch)
  end subroutine read_parameters

  !> Reset the module, freeing memory etc.
  subroutine finish_split_nonlinear_terms
    implicit none
    initialized = .false.
    call split_nonlinear_terms_config%reset()
  end subroutine finish_split_nonlinear_terms

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by calling
  !> specific requested method.
  !> On output g = gnew, chi = chinew.
  subroutine advance_nonlinear_term(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use dist_fn_arrays, only: gnew, g
    use mp, only: mp_abort, get_mp_times
    use job_manage, only: time_message
    use nonlinear_terms, only: time_add_explicit_terms, time_add_explicit_terms_mpi
    implicit none
    complex, dimension (:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real :: mp_total_after, mp_total
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total)

    ! Call requested method for advancing nonlinear term
    select case(split_method_switch)
    case(split_method_ab3)
       call advance_nonlinear_term_ab3(phi, apar, bpar, phinew, aparnew, bparnew, &
            fields_func, source_func)
    case(split_method_backwards_euler)
       call advance_nonlinear_term_beuler(phi, apar, bpar, phinew, aparnew, bparnew, &
            fields_func, source_func)
    case(split_method_rk)
       call advance_nonlinear_term_rk(phi, apar, bpar, phinew, aparnew, bparnew, &
            fields_func, source_func)
    case default
       call mp_abort("Invalid split_method_switch", .true.)
    end select

    ! Now that we've advanced the explicit terms by code_dt set the current state to this new result
    g = gnew
    phi = phinew
    apar = aparnew
    bpar = bparnew

    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total_after)
    time_add_explicit_terms_mpi = time_add_explicit_terms_mpi + (mp_total_after - mp_total)
  end subroutine advance_nonlinear_term

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by using
  !> AB scheme of orders up to 3rd.
  subroutine advance_nonlinear_term_ab3(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    use gs2_time, only: get_adams_bashforth_coefficients, code_dt
    use dist_fn_arrays, only: gnew, gexp_1, gexp_2, gexp_3, g
    use gs2_time, only: code_dt_min, dt_not_set, code_dt_prev1, code_dt_prev2
    use mp, only: mp_abort
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real :: internal_time, time_step
    real, dimension(:), allocatable :: ab_coeffs
    integer :: counter
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    ! Initialise internal variables
    internal_time = 0.0
    counter = 1

    ! Reset time step sizes to reflect fact we're restarting the explicit scheme
    code_dt_prev1 = dt_not_set
    code_dt_prev2 = dt_not_set

    ! Keep taking explicit steps until we have advanced by code_dt
    do while (internal_time < code_dt)
       ! Get the current nonlinear source term and cfl limit
       call source_func(gexp_1, phinew, aparnew, bparnew, &
            time_step, calculate_cfl_limit = .true.)

       ! If the cfl time step is too small then abort
       if (time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

       ! Now make sure we're don't overshoot the target time by limiting time step.
       time_step = min(time_step, code_dt - internal_time)

       ! Get the Adams-Bashforth coefficients for the current collection of time step sizes
       ab_coeffs = get_adams_bashforth_coefficients(target_dt = time_step)

       ! Advance solution by time_step
       select case(counter)
       case(1)
          gnew = gnew + time_step * 0.5 *( &
               ab_coeffs(1) * gexp_1       &
               )
       case(2)
          gnew = gnew + time_step * 0.5 *( &
               ab_coeffs(1) * gexp_1 +     &
               ab_coeffs(2) * gexp_2       &
               )
       case default
          gnew = gnew + time_step * 0.5 *( &
               ab_coeffs(1) * gexp_1 +     &
               ab_coeffs(2) * gexp_2 +     &
               ab_coeffs(3) * gexp_3       &
               )
       end select

       ! Calculate the fields consistent with the new g
       call fields_func(phinew, aparnew, bparnew)

       ! Increment internal state
       internal_time = internal_time + time_step
       counter = counter + 1

       ! Shift time history along one point
       code_dt_prev2 = code_dt_prev1
       code_dt_prev1 = time_step
       gexp_3 = gexp_2
       gexp_2 = gexp_1
    end do

  end subroutine advance_nonlinear_term_ab3

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by using
  !> backwards Euler with fixed-point iteration.
  subroutine advance_nonlinear_term_beuler(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, gexp_1, gexp_2, g, g_work
    use gs2_time, only: code_dt_min, code_dt
    use mp, only: mp_abort, max_allreduce, proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real :: internal_time, time_step, time_step_limit
    integer :: counter
    integer, parameter :: iteration_limit = 20
    integer, parameter :: iteration_lower_limit = 7
    real, parameter :: time_adjust_factor = 2
    logical, parameter :: debug = .false.
    real :: error, actual_dt
    real, dimension(2) :: errors
    real, save :: time_step_last = -1
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func
    logical :: reject_step

    ! Initialise internal variables
    if (time_step_last < 0) time_step_last = code_dt
    internal_time = 0.0
    time_step = time_step_last

    ! Keep on advancing until we have advanced by code_dt
    do while (internal_time < code_dt)
       ! Reset counter and error
       counter = 0
       error = 1.0

       ! Limit the actual time step that we take to avoid overshooting (i.e.
       ! to keep internal_time <= code_dt).
       actual_dt = min(time_step, code_dt - internal_time)

       ! Calculate the current NL source term, for use in error
       ! estimates later.
       call source_func(gexp_2, phinew, aparnew, bparnew, &
            time_step_limit, calculate_cfl_limit = .false.)

       ! Try to take a step of size actual_dt using fixed point iteration with
       ! backwards Euler.
       do while (error > convergence_tolerance .and. counter < iteration_limit)

          ! Get the current nonlinear source term and cfl limit -- note we don't use the
          ! cfl limit here
          call source_func(gexp_1, phinew, aparnew, bparnew, &
               time_step_limit, calculate_cfl_limit = .false.)

          ! Calculate the estimate of the new g
          g_work = g + actual_dt * 0.5 * gexp_1

          ! Calculate the maximum absolute error and the maximum value of the previous solution
          errors(1) = maxval(abs(gnew - g_work))
          errors(2) = maxval(abs(g_work))
          call max_allreduce(errors)

          ! Calculate an approximate relative error
          error = errors(1) / errors(2)

          ! Update internals -- increment counter and update the 'previous' solution with the new one.
          counter = counter + 1
          gnew = g_work

          ! Update the fields consistently with the current solution
          call fields_func(phinew, aparnew, bparnew)
       end do

       ! Now we need to decide if the last step was successful or not. There are two exit conditions, either the
       ! error was small enough (good step) or the iteration limit was exceeded (bad step).
       if (counter >= iteration_limit) then
          ! If it was a bad step then throw away these changes and reduce the time step
          if (proc0 .and. debug) print*,'Too many iterations, reducing step from',time_step,'to',time_step/time_adjust_factor,'with internal_time = ',internal_time

          ! Reduce the time step
          time_step = time_step / time_adjust_factor
          ! If the time step is now too small then abort
          if (time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

          ! If we've taken too many iterations then we want to drop the time step and try
          ! again. To try again we must ensure that gnew has been reset to g, and the fields
          ! are consistent with this.
          reject_step = .true.
       else
          ! If we have met the convergence condition then we should estimate the error
          ! on the solution and see if we need to reject the step or how we should adjust
          ! the step size.
          ! To estimate the error we compare the backward Euler update with an approximate
          ! Crank-Nicolson update.
          ! g_new_beuler = g + 0.5 * actual_dt * gexp_1(g_new_beuler)
          ! g_new_cn = g + 0.5 * actual_dt * (gexp_1(g_new_beuler) + gexp_1(g))/2
          ! abs_error  = |g_new_beuler - g_new_cn|
          !      = 0.5 * actual_dt * 0.5 * (gexp_1(g_new_beuler) - gexp_1(g)
          ! Where gexp_1(h) is the source func evaluated with gnew=h.

          ! Here we calculate gexp_1(g_new_beuler) exactly. We could probably
          ! rely on using the existing value of gexp_1(g_work) as, whilst this is
          ! strictly from the previous iteration of the fixed-point scheme, we know
          ! it is within the convergence tolerance of the true value.
          call source_func(gexp_1, phinew, aparnew, bparnew, &
               time_step_limit, calculate_cfl_limit = .false.)

          errors(1) = maxval(abs(gexp_1 - gexp_2)) * 0.5 * 0.5 * actual_dt
          call max_allreduce(errors)
          error = errors(1) / errors(2)

          !Calculate new time step based on error
          time_step = 0.9 * actual_dt * sqrt(error_tolerance/(error+1.0e-14))

          ! If error is too large reject the step
          if ( error > error_tolerance) then
             if (proc0 .and. debug) print*,'Attempted beuler step failed due to error, adapting time step'
             ! Reset as we're retrying the step.
             reject_step = .true.
          else
             reject_step = .false.
          end if
       end if

       if (reject_step) then
          ! If the time step is now too small then abort
          if (time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)
          gnew = g
          call fields_func(phinew, aparnew, bparnew)
       else
          ! If we have met the error condition and not taken too many steps
          ! then gnew represents our new solution so make sure we update g to reflect this
          internal_time = internal_time + actual_dt
          g = gnew
       end if

    end do

    ! Store the final time step size so we can start from this in the next callg
    time_step_last = time_step

    if(proc0.and.debug) print*,'Explicit done in ',counter,'steps with final error',error,'and step',time_step,'final time',internal_time

  end subroutine advance_nonlinear_term_beuler

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by using
  !> an Runge-Kutta scheme with embedded error estimate for error control. This method
  !> is a wrapper to the true generic RK implementation.
  subroutine advance_nonlinear_term_rk(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    call advance_nonlinear_term_rk_implementation( &
         phi, apar, bpar, phinew, aparnew, bparnew, &
         fields_func, source_func, the_rk_scheme)
  end subroutine advance_nonlinear_term_rk

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by using
  !> an Runge-Kutta scheme with embedded error estimate for error control
  subroutine advance_nonlinear_term_rk_implementation(phi, apar, bpar, phinew, &
       aparnew, bparnew, fields_func, source_func, scheme)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g, gexp_1
    use gs2_time, only: code_dt_min, code_dt
    use mp, only: mp_abort, proc0, max_allreduce
    use gs2_layouts, only: g_lo
    use rk_schemes, only: rk_scheme
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    type(rk_scheme), intent(in) :: scheme
    complex, dimension (:, :, :, :), allocatable :: stages
    real :: internal_time, time_step, time_step_limit, error, actual_time_step
    real :: new_time_step
    integer :: counter, nstages
    logical, parameter :: debug = .false.
    integer :: istage, istage_sub, validity
    real, dimension(:), allocatable :: solution_coeffs
    real, save :: time_step_last = -1
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    ! Initialise internal variables
    if (time_step_last < 0) time_step_last = code_dt
    internal_time = 0.0
    counter = 1
    time_step = time_step_last

    ! Copy the relevant coefficients for forming the solution to reduce
    ! later duplication.
    if (scheme%follow_high_order) then
       solution_coeffs = scheme%high_order_coeffs
    else
       solution_coeffs = scheme%lower_order_coeffs
    end if

    nstages = scheme%number_of_stages
    allocate(stages(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_proc, nstages))

    ! Keep taking explicit steps until we have advanced by code_dt
    do while (internal_time < code_dt)
       stages = 0.0
       if(proc0.and.debug) print*,'Taking step with size',time_step,'at time',internal_time
       do istage = 1, nstages
          gnew = g
          do istage_sub = 1, istage-1
             gnew = gnew + time_step * scheme%coeffs(istage_sub, istage) * stages(:, :, :, istage_sub)
          end do
          call fields_func(phinew, aparnew, bparnew)

          ! Get the current nonlinear source term
          call source_func(gexp_1, &
               phinew, aparnew, bparnew, &
               time_step_limit, calculate_cfl_limit = .false.)
          stages(:, :, :, istage) = 0.5 * gexp_1
       end do

       ! One can directly form the error by simply forming error =
       ! time_step * sum(stages_i * (high_order_coeffs_i -
       ! low_order_coeffs_i)). However we also need to explicitly form one of the low/high
       ! order solutions in order to be able to construct the solution and the relative error.
       ! Here we use just the absolute error, stored in gnew,
       gnew = (scheme%high_order_coeffs(1) - scheme%lower_order_coeffs(1)) * stages(:, :, :, 1)
       do istage = 2, nstages
          gnew = gnew + (scheme%high_order_coeffs(istage) - scheme%lower_order_coeffs(istage)) * stages(:, :, :, istage)
       end do

       error = maxval(abs(gnew))*time_step
       call max_allreduce(error)

       !Calculate the new_time_step
       new_time_step = 0.9 * time_step * (error_tolerance/(error+1.0e-14))**(1.0/scheme%order)
       if(proc0.and.debug) print*,'Error is ',error,'at internal_time',internal_time,'with step',time_step,'new time step is',new_time_step,'counter=',counter
       ! If the cfl time step is too small then abort
       if (new_time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

       ! If the error is too large we need to retry with a smaller step, so don't
       ! update anything
       if (error > error_tolerance) then
          if(proc0.and.debug) print*,'Time step with size',time_step,'failed at internal_time',internal_time,'and counter',counter,'retrying with step',new_time_step
       else
          !Update the solution
          g = g + solution_coeffs(1) * stages(:, :, :, 1) * time_step
          do istage = 2, nstages
             g = g + solution_coeffs(istage) * stages(:, :, :, istage) * time_step
          end do

          internal_time = internal_time + time_step
       end if

       !Now update the time step, capped to code_dt
       time_step = min(new_time_step, code_dt)

       ! Make sure we don't overshoot the target time by limiting time step.
       if(internal_time < code_dt) then
          time_step = min(time_step, code_dt - internal_time)
          counter = counter + 1
       end if

    end do
    if(proc0.and.debug) print*,'Completed rk step in',counter,'steps to',internal_time
    gnew = g
    call fields_func(phinew, aparnew, bparnew)

    time_step_last = time_step
  end subroutine advance_nonlinear_term_rk_implementation

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------

  !> Reads in the split_nonlinear_terms_knobs namelist and populates the member variables
  subroutine read_split_nonlinear_terms_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(split_nonlinear_terms_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: rk_method, split_method
    real :: convergence_tolerance, error_tolerance

    namelist /split_nonlinear_terms_knobs/ convergence_tolerance, error_tolerance, rk_method, split_method

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    convergence_tolerance = self%convergence_tolerance
    error_tolerance = self%error_tolerance
    rk_method = self%rk_method
    split_method = self%split_method

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = split_nonlinear_terms_knobs)

    ! Now copy from local variables into type members
    self%convergence_tolerance = convergence_tolerance
    self%error_tolerance = error_tolerance
    self%rk_method = rk_method
    self%split_method = split_method

    self%exist = exist
  end subroutine read_split_nonlinear_terms_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_split_nonlinear_terms_config(self, unit)
    implicit none
    class(split_nonlinear_terms_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("convergence_tolerance", self%convergence_tolerance, unit_internal)
    call self%write_key_val("error_tolerance", self%error_tolerance, unit_internal)
    call self%write_key_val("rk_method", self%rk_method, unit_internal)
    call self%write_key_val("split_method", self%split_method, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_split_nonlinear_terms_config

  !> Resets the config object to the initial empty state
  subroutine reset_split_nonlinear_terms_config(self)
    class(split_nonlinear_terms_config_type), intent(in out) :: self
    type(split_nonlinear_terms_config_type) :: empty
    select type (self)
    type is (split_nonlinear_terms_config_type)
       self = empty
    end select
  end subroutine reset_split_nonlinear_terms_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_split_nonlinear_terms_config(self)
    use mp, only: broadcast
    implicit none
    class(split_nonlinear_terms_config_type), intent(in out) :: self
    call broadcast(self%convergence_tolerance)
    call broadcast(self%error_tolerance)
    call broadcast(self%rk_method)
    call broadcast(self%split_method)

    call broadcast(self%exist)
  end subroutine broadcast_split_nonlinear_terms_config

  !> Gets the default name for this namelist
  function get_default_name_split_nonlinear_terms_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_split_nonlinear_terms_config
    get_default_name_split_nonlinear_terms_config = "split_nonlinear_terms_knobs"
  end function get_default_name_split_nonlinear_terms_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_split_nonlinear_terms_config()
    implicit none
    logical :: get_default_requires_index_split_nonlinear_terms_config
    get_default_requires_index_split_nonlinear_terms_config = .false.
  end function get_default_requires_index_split_nonlinear_terms_config

  !> Get the module level config instance
  pure function get_split_nonlinear_terms_config()
    type(split_nonlinear_terms_config_type) :: get_split_nonlinear_terms_config
    get_split_nonlinear_terms_config = split_nonlinear_terms_config
  end function get_split_nonlinear_terms_config

end module split_nonlinear_terms
