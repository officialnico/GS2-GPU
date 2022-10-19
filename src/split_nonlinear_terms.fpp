!> A module to deal with advancing the nonlinear term
!> separately from the linear terms.
module split_nonlinear_terms
  implicit none
  private

  public :: advance_nonlinear_term

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

contains

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+code_dt by calling
  !> specific requested method.
  !> On output g = gnew, chi = chinew.
  subroutine advance_nonlinear_term(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use dist_fn_arrays, only: gnew, g
    use mp, only: mp_abort, get_mp_times
    use job_manage, only: time_message
    use nonlinear_terms, only: split_method_switch, split_method_ab3, split_method_backwards_euler
    use nonlinear_terms, only: split_method_rk
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
    use nonlinear_terms, only: calculate_current_nl_source_and_cfl_limit
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
    use nonlinear_terms, only: error_target, calculate_current_nl_source_and_cfl_limit, cfl
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real :: internal_time, time_step, time_step_limit
    integer :: counter
    real :: tolerance, convergence_tolerance
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
    tolerance = error_target
    convergence_tolerance = cfl

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
          time_step = 0.9 * actual_dt * sqrt(tolerance/(error+1.0e-14))

          ! If error is too large reject the step
          if ( error > tolerance) then
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
  !> an Runge-Kutta scheme with embedded error estimate for error control
  subroutine advance_nonlinear_term_rk(phi, apar, bpar, phinew, aparnew, bparnew, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g, gexp_1
    use gs2_time, only: code_dt_min, code_dt
    use mp, only: mp_abort, proc0, max_allreduce
    use nonlinear_terms, only: calculate_current_nl_source_and_cfl_limit, error_target
    use gs2_layouts, only: g_lo
    use rk_schemes, only: rk_scheme, get_rk_scheme, rk_cashkarp
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    complex, dimension (:, :, :, :), allocatable :: stages
    real :: internal_time, time_step, time_step_limit, error, actual_time_step
    real :: tolerance, new_time_step
    integer :: counter, nstages
    logical, parameter :: debug = .false.
    ! Define the RK scheme
    type(rk_scheme) :: scheme
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
    tolerance = error_target

    scheme = get_rk_scheme(rk_cashkarp)
    validity = scheme%validate()
    if (proc0 .and. (validity /= 0)) then
       print*,'The RK scheme failed validity check with error code:',validity
       call mp_abort("RK scheme failed validity check", .true.)
    end if

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
       new_time_step = 0.9 * time_step * (tolerance/(error+1.0e-14))**(1.0/scheme%order)
       if(proc0.and.debug) print*,'Error is ',error,'at internal_time',internal_time,'with step',time_step,'new time step is',new_time_step,'counter=',counter
       ! If the cfl time step is too small then abort
       if (new_time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

       ! If the error is too large we need to retry with a smaller step, so don't
       ! update anything
       if (error > tolerance) then
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
  end subroutine advance_nonlinear_term_rk

end module split_nonlinear_terms
