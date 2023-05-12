!> A module to deal with advancing the nonlinear term
!> separately from the linear terms.
module split_nonlinear_terms
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use rk_schemes, only: rk_scheme_type
  implicit none
  private

  public :: init_split_nonlinear_terms, finish_split_nonlinear_terms
  public :: advance_nonlinear_term, advance_nonadiabatic_dfn
  public :: split_nonlinear_terms_config_type, get_split_nonlinear_terms_config

  logical :: initialized = .false.

  !> Statistics for the number of steps taken in the integrators
  integer :: nfailed_steps = 0, nrhs_calls = 0, nsuccessful_steps = 0

  !> Interface to describe the call signatures for routines to calculate the
  !> field from gnew (invert_field_func) and to calculate the source term used
  !> in the advance_nonlinear_term methods. Used to allow us to pass in these
  !> methods, giving us a way to test the integrators on different problems.
  interface
     subroutine invert_field_func (g_in, phi, apar, bpar, gf_lo, local_only)
       complex, dimension (:,:,:), intent(in) :: g_in
       complex, dimension (:,:,:), intent (out) :: phi, apar, bpar
       logical, optional :: gf_lo, local_only
     end subroutine invert_field_func

     subroutine source_term_func(g_in, g1, phi, apar, bpar, max_vel_local, &
          need_to_adjust, calculate_cfl_limit)
       complex, dimension (:,:,:), intent (in) :: g_in
       complex, dimension (:,:,:), intent (out) :: g1
       complex, dimension (:,:,:), intent (in) :: phi, apar, bpar
       real, intent(out) :: max_vel_local
       logical, intent(in) :: need_to_adjust
       logical, intent(in), optional :: calculate_cfl_limit
     end subroutine source_term_func
  end interface

  ! Related to inputs
  integer :: split_method_switch
  integer, parameter :: split_method_ab3 = 1, split_method_backwards_euler = 2, &
       split_method_rk = 3
  type(rk_scheme_type) :: the_rk_scheme
  logical :: show_statistics, advance_nonadiabatic_dfn
  real :: convergence_tolerance, time_step_safety_factor
  real :: relative_tolerance, absolute_tolerance

  !> Used to represent the input configuration of split_nonlinear_terms
  type, extends(abstract_config_type) :: split_nonlinear_terms_config_type
     ! namelist : split_nonlinear_terms_knobs
     ! indexed : false
     !> The absolute tolerance used in error controlled schemes. Attempt to adjust
     !> time step to keep error below this tolerance.
     real :: absolute_tolerance = 1.0e-3
     !> If true then nonlinear integrator advances the nonadiabatic dfn, h,
     !> rather than the modified distribution function, g.
     logical :: advance_nonadiabatic_dfn = .false.
     !> The tolerance used in convergence checks. Currently only used
     !> for halting the fixed-point iteration in the beuler method
     real :: convergence_tolerance = 1.0e-1
     !> The relative tolerance used in error controlled schemes. Attempt to adjust
     !> time step to keep error below this tolerance.
     real :: relative_tolerance = 1.0e-2
     !> Choose which specific RK scheme to use in the RK method.
     !> Generally if the tolerances are small a higher order scheme
     !> will be more effective than low order. At loose tolerances
     !> low order schemes are generally more efficient
     !> Valid options include:
     !>
     !>  - 'cashkarp' -- Cash-Karp scheme of order 4(5)
     !>  - 'default' -- Same as heun
     !>  - 'heun' -- Heun scheme of order 1(2)
     !>
     !> See [rk_schemes] for other options.
     character(len = 20) :: rk_method = 'default'
     !> If true then reports the number of right hand side (NL term) calculations,
     !> the number of failed internal steps and the number of successful steps at the
     !> end of each linear size step.
     logical :: show_statistics = .false.
     !> What algorithm do we use to advance the nonlinear term if
     !> split_nonlinear is true.
     !> Valid options include:
     !>
     !>  - 'AB3' -- Adams-Bashforth (up to) 3rd order. CLF controlled time step.
     !>  - 'beuler' -- Backwards Euler. Relative and absolute error controlled time step.
     !>            Newton iteration controlled by convergence_tolerance.
     !>  - 'default' -- The same as 'RK'.
     !>  - 'RK' -- RK scheme with embedded error estimate. Exact scheme can be switched out
     !>            using rk_method switch. Relative and absolute error controlled time step.
     !>
     character(len = 20) :: split_method = 'default'
     !> Multiplies the new time step calculated from either the error or
     !> cfl limits. Smaller values are more conservative but may help avoid
     !> repeated violation of the error/cfl limits and hence reduce the
     !> number of failed steps.
     real :: time_step_safety_factor = 0.9
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
         text_option('default', split_method_rk), &
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
    absolute_tolerance = split_nonlinear_terms_config%absolute_tolerance
    advance_nonadiabatic_dfn = split_nonlinear_terms_config%advance_nonadiabatic_dfn
    convergence_tolerance = split_nonlinear_terms_config%convergence_tolerance
    relative_tolerance = split_nonlinear_terms_config%relative_tolerance
    rk_method = split_nonlinear_terms_config%rk_method
    show_statistics = split_nonlinear_terms_config%show_statistics
    split_method = split_nonlinear_terms_config%split_method
    time_step_safety_factor = split_nonlinear_terms_config%time_step_safety_factor

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

  !> Reset the step statistics
  subroutine reset_step_statistics
    nfailed_steps = 0 ; nrhs_calls = 0 ; nsuccessful_steps = 0
  end subroutine reset_step_statistics

  subroutine report_step_statistics(istep, unit_in)
    use iso_fortran_env, only: output_unit
    use optionals, only: get_option_with_default
    implicit none
    integer, intent(in) :: istep
    integer, intent(in), optional :: unit_in
    integer :: unit
    unit = get_option_with_default(unit_in, output_unit)
    write(unit,'("Iteration : ",I0," Number of internal steps : ",I0," Number of failed steps : ",I0," Number of RHS calls : ",I0)') &
         istep, nsuccessful_steps, nfailed_steps, nrhs_calls
  end subroutine report_step_statistics

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+dt by calling
  !> specific requested method.
  !> On output g = gnew, chi = chinew.
  subroutine advance_nonlinear_term(istep, phi, apar, bpar, dt, fields_func, source_func)
    use dist_fn_arrays, only: gnew, g, g_adjust, to_g_gs2, from_g_gs2
    use mp, only: mp_abort, get_mp_times, proc0
    use job_manage, only: time_message
    use nonlinear_terms, only: time_add_explicit_terms, time_add_explicit_terms_mpi
    use run_parameters, only: has_phi, has_apar, has_bpar
    use array_utils, only: copy
    implicit none
    integer, intent(in) :: istep
    complex, dimension (:,:,:), intent (in) :: phi, apar, bpar
    real, intent(in) :: dt
    complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
    real :: mp_total_after, mp_total
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total)
    call reset_step_statistics

    ! Copy current fields to local copies
    allocate(phinew, mold = phi)
    allocate(aparnew, mold = apar)
    allocate(bparnew, mold = bpar)
    if (has_phi) call copy(phi, phinew)
    if (has_apar) call copy(apar, aparnew)
    if (has_bpar) call copy(bpar, bparnew)

    ! If we want to advance h rather than g then adjust
    ! here to get h. Note we rely on the call site also
    ! providing the correct fields_func to calculate fields
    ! from h rather than g. A different option would be
    ! to require the call site to always pass a way to calculate
    ! from g and a way from h and we then decide here which to use.
    if (advance_nonadiabatic_dfn) then
       call g_adjust(g, phinew, bparnew, direction = from_g_gs2)
       call copy(g, gnew)
    end if

    ! Call requested method for advancing nonlinear term
    select case(split_method_switch)
    case(split_method_ab3)
       call advance_nonlinear_term_ab3(phinew, aparnew, bparnew, dt, &
            fields_func, source_func)
    case(split_method_backwards_euler)
       call advance_nonlinear_term_beuler(phinew, aparnew, bparnew, dt, &
            fields_func, source_func)
    case(split_method_rk)
       call advance_nonlinear_term_rk(phinew, aparnew, bparnew, dt, &
            fields_func, source_func)
    case default
       call mp_abort("Invalid split_method_switch", .true.)
    end select

    if (proc0 .and. show_statistics) call report_step_statistics(istep)

    ! If we advanced h then we need to reconstruct g from this
    if (advance_nonadiabatic_dfn) then
       call g_adjust(gnew, phi, bpar, direction = to_g_gs2)
    end if

    ! Now that we've advanced the explicit terms by dt set the
    ! current state to this new result. Note that we do not update the
    ! fields, instead resetting them to their original values. The
    ! update to the fields is calculated as a part of the implicit
    ! field solve and including the NL contribution here leads to
    ! "double counting".
    call copy(gnew, g)

    call time_message(.false., time_add_explicit_terms, 'Explicit terms')
    call get_mp_times(total_time = mp_total_after)
    time_add_explicit_terms_mpi = time_add_explicit_terms_mpi + (mp_total_after - mp_total)
  end subroutine advance_nonlinear_term

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+dt by using
  !> AB scheme of orders up to 3rd.
  subroutine advance_nonlinear_term_ab3(phinew, aparnew, bparnew, dt, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    use gs2_time, only: get_adams_bashforth_coefficients
    use dist_fn_arrays, only: gnew, gexp_1, gexp_2, gexp_3, g
    use gs2_time, only: code_dt_min, dt_not_set, code_dt_prev1, code_dt_prev2
    use mp, only: mp_abort, max_allreduce
    use gs2_layouts, only: g_lo
    use array_utils, only: copy
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real, intent(in) :: dt
    real :: internal_time, time_step, max_vel, half_time
    real, dimension(:), allocatable :: ab_coeffs
    integer :: counter, iglo
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    ! Initialise internal variables
    internal_time = 0.0
    counter = 1

    ! Reset time step sizes to reflect fact we're restarting the explicit scheme
    code_dt_prev1 = dt_not_set
    code_dt_prev2 = dt_not_set

    ! Keep taking explicit steps until we have advanced by dt
    do while (internal_time < dt)
       ! Get the current nonlinear source term and cfl limit
       call source_func(gnew, gexp_1, phinew, aparnew, bparnew, &
            max_vel, need_to_adjust = .not. advance_nonadiabatic_dfn, &
            calculate_cfl_limit = .true.)
       nrhs_calls = nrhs_calls + 1

       ! Get the maximum velocity across all processors
       ! and convert to time_step estimate.
       call max_allreduce(max_vel)
       time_step = 1.0 / max_vel
       time_step = time_step * time_step_safety_factor

       ! If the cfl time step is too small then abort
       if (time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

       ! Now make sure we're don't overshoot the target time by limiting time step.
       time_step = min(time_step, dt - internal_time)

       ! Get the Adams-Bashforth coefficients for the current collection of time step sizes
       ab_coeffs = get_adams_bashforth_coefficients(target_dt = time_step)

       half_time = 0.5 * time_step

       ! Advance solution by time_step
       select case(counter)
       case(1)
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, gnew, half_time, ab_coeffs, gexp_1) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             gnew(:, :, iglo) = gnew(:, :, iglo) + half_time *( &
                  ab_coeffs(1) * gexp_1(:, :, iglo)       &
                  )
          end do
          !$OMP END PARALLEL DO
       case(2)
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, gnew, half_time, ab_coeffs, gexp_1, gexp_2) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             gnew(:, :, iglo) = gnew(:, :, iglo) + half_time *( &
                  ab_coeffs(1) * gexp_1(:, :, iglo) +     &
                  ab_coeffs(2) * gexp_2(:, :, iglo)       &
                  )
          end do
          !$OMP END PARALLEL DO
       case default
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, gnew, half_time, ab_coeffs, gexp_1, gexp_2, gexp_3) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             gnew(:, :, iglo) = gnew(:, :, iglo) + half_time *( &
                  ab_coeffs(1) * gexp_1(:, :, iglo) +     &
                  ab_coeffs(2) * gexp_2(:, :, iglo) +     &
                  ab_coeffs(3) * gexp_3(:, :, iglo)       &
                  )
          end do
          !$OMP END PARALLEL DO
       end select

       ! Calculate the fields consistent with the new g
       call calculate_fields(fields_func, gnew, phinew, aparnew, bparnew)

       ! Increment internal state
       internal_time = internal_time + time_step
       counter = counter + 1
       nsuccessful_steps = nsuccessful_steps + 1

       ! Shift time history along one point
       code_dt_prev2 = code_dt_prev1
       code_dt_prev1 = time_step
       call copy(gexp_2, gexp_3)
       call copy(gexp_1, gexp_2)
    end do

  end subroutine advance_nonlinear_term_ab3

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+dt by using
  !> backwards Euler with fixed-point iteration.
  subroutine advance_nonlinear_term_beuler(phinew, aparnew, bparnew, dt, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, gexp_1, gexp_2, g, g_work
    use gs2_time, only: code_dt_min
    use mp, only: mp_abort, max_allreduce, proc0
    use run_parameters, only: immediate_reset
    use gs2_layouts, only: g_lo
    use array_utils, only: copy, gs2_max_abs
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real, intent(in) :: dt
    real :: internal_time, time_step, max_vel
    integer :: counter, iglo
    integer, parameter :: iteration_limit = 20
    integer, parameter :: iteration_lower_limit = 7
    real, parameter :: time_adjust_factor = 2
    logical, parameter :: debug = .false.
    real :: inverse_error, error, actual_dt, half_time, cfl_time_step
    real, dimension(3) :: errors
    real, save :: time_step_last = -1
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func
    logical :: reject_step, cfl_limited
    complex, dimension(:, :, :), allocatable :: error_array

    ! Initialise internal variables
    if (time_step_last < 0) time_step_last = dt
    internal_time = 0.0
    time_step = time_step_last

    allocate(error_array(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))

    ! Keep on advancing until we have advanced by dt
    do while (internal_time < dt)
       ! Reset counter and error
       counter = 0
       error = 1.0
       errors = 0.0

       ! Limit the actual time step that we take to avoid overshooting (i.e.
       ! to keep internal_time <= dt).
       actual_dt = min(time_step, dt - internal_time)
       half_time = 0.5 * actual_dt

       ! Calculate the current NL source term, for use in error
       ! estimates later.
       call source_func(gnew, gexp_2, phinew, aparnew, bparnew, &
            max_vel, need_to_adjust = .not. advance_nonadiabatic_dfn, &
            calculate_cfl_limit = .false.)
       nrhs_calls = nrhs_calls + 1

       ! Try to take a step of size actual_dt using fixed point iteration with
       ! backwards Euler.
       do while (error > convergence_tolerance .and. counter < iteration_limit)

          ! Get the current nonlinear source term and cfl limit -- note we don't use the
          ! cfl limit here
          call source_func(gnew, gexp_1, phinew, aparnew, bparnew, &
               max_vel, need_to_adjust = .not. advance_nonadiabatic_dfn, &
               calculate_cfl_limit = .false.)
          nrhs_calls = nrhs_calls + 1

          ! Calculate the estimate of the new g
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, g_work, g, half_time, gexp_1) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             g_work(:, :, iglo) = g(:, :, iglo) + half_time * gexp_1(:, :, iglo)
          end do
          !$OMP END PARALLEL DO

          ! Calculate the maximum absolute error and the maximum value of the previous solution

          ! Note we _could_ merge the following loop witth the gs2_max_abs kernel
          ! so that we only loop over the error_array once (and actually don't
          ! need to explicitly store the full array)
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, error_array, gnew, g_work) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             error_array(:, :, iglo) = gnew(:, :, iglo) - g_work(:, :, iglo)
          end do
          !$OMP END PARALLEL DO
          errors(1) = gs2_max_abs(error_array)
          errors(2) = gs2_max_abs(g_work)
          call max_allreduce(errors)

          ! Calculate an approximate relative error
          error = errors(1) / errors(2)

          ! Update internals -- increment counter and update the 'previous' solution with the new one.
          counter = counter + 1
          call copy(g_work, gnew)

          ! Update the fields consistently with the current solution
          call fields_func(gnew, phinew, aparnew, bparnew, local_only = .true.)
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
          call source_func(gnew, gexp_1, phinew, aparnew, bparnew, &
               max_vel, need_to_adjust = .not. advance_nonadiabatic_dfn, &
               calculate_cfl_limit = .true.)
          nrhs_calls = nrhs_calls + 1

          ! Note we _could_ merge the following loop witth the gs2_max_abs kernel
          ! so that we only loop over the error_array once (and actually don't
          ! need to explicitly store the full array)
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, error_array, gexp_1, gexp_2) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             error_array(:, :, iglo) = gexp_1(:, :, iglo) - gexp_2(:, :, iglo)
          end do
          !$OMP END PARALLEL DO
          errors(1) = gs2_max_abs(error_array) * (0.5 * 0.5 * actual_dt)
          errors(3) = max_vel
          call max_allreduce(errors)
          inverse_error = get_inverse_error(errors)
          cfl_time_step = 1.0 / errors(3)

          !Calculate new time step based on error
          time_step = min(actual_dt * sqrt(inverse_error), cfl_time_step)
          cfl_limited = time_step < actual_dt .and. time_step == cfl_time_step
          time_step = time_step * time_step_safety_factor

          ! If error is too large reject the step
          if ( inverse_error < 1.0 .or. (cfl_limited .and. immediate_reset)) then
             if (proc0 .and. debug) print*,'Attempted beuler step failed due to error or cfl, adapting time step'
             ! Reset as we're retrying the step.
             reject_step = .true.
          else
             reject_step = .false.
          end if
       end if

       if (reject_step) then
          ! If the time step is now too small then abort
          if (time_step < code_dt_min) call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)
          call copy(g, gnew)
          ! If we store the initial field then we could perhaps just copy here instead.
          ! This could be faster as it avoids communications involved in the velocity integration.
          call calculate_fields(fields_func, gnew, phinew, aparnew, bparnew)
          nfailed_steps = nfailed_steps + 1
       else
          ! If we have met the error condition and not taken too many steps
          ! then gnew represents our new solution so make sure we update g to reflect this
          internal_time = internal_time + actual_dt
          call copy(gnew, g)
          nsuccessful_steps = nsuccessful_steps + 1
       end if

    end do

    ! Store the final time step size so we can start from this in the next callg
    time_step_last = time_step

    if(proc0.and.debug) print*,'Explicit done in ',counter,'steps with final error',inverse_error,'and step',time_step,'final time',internal_time

  end subroutine advance_nonlinear_term_beuler

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+dt by using
  !> an Runge-Kutta scheme with embedded error estimate for error control. This method
  !> is a wrapper to the true generic RK implementation.
  subroutine advance_nonlinear_term_rk(phinew, aparnew, bparnew, dt, &
       fields_func, source_func)
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real, intent(in) :: dt
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func

    call advance_nonlinear_term_rk_implementation( &
         phinew, aparnew, bparnew, dt, &
         fields_func, source_func, the_rk_scheme)
  end subroutine advance_nonlinear_term_rk

  !> Advances dg/dt = NL(g,chi) from current gnew, chi from t -> t+dt by using
  !> an Runge-Kutta scheme with embedded error estimate for error control
  subroutine advance_nonlinear_term_rk_implementation(phinew, &
       aparnew, bparnew, dt, fields_func, source_func, scheme)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g, gexp_1
    use gs2_time, only: code_dt_min
    use mp, only: mp_abort, proc0, max_allreduce
    use gs2_layouts, only: g_lo
    use rk_schemes, only: rk_scheme_type
    use run_parameters, only: immediate_reset
    use array_utils, only: copy, gs2_max_abs, zero_array
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    real, intent(in) :: dt
    type(rk_scheme_type), intent(in) :: scheme
    complex, dimension (:, :, :, :), allocatable :: stages
    real :: internal_time, time_step, max_vel, inverse_error, actual_time_step
    real :: new_time_step, cfl_time_step
    integer :: counter, nstages
    logical, parameter :: debug = .false.
    integer :: istage, istage_sub, validity, iglo
    real, dimension(:), allocatable :: solution_coeffs, error_coeffs
    real, dimension(3) :: errors
    real, save :: time_step_last = -1
    procedure(invert_field_func) :: fields_func
    procedure(source_term_func) :: source_func
    logical :: cfl_limited

    ! Initialise internal variables
    if (time_step_last < 0) time_step_last = dt
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

    ! Get the error coefficients
    error_coeffs = scheme%high_order_coeffs - scheme%lower_order_coeffs

    nstages = scheme%number_of_stages
    allocate(stages(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc, nstages))

    ! Keep taking explicit steps until we have advanced by dt
    do while (internal_time < dt)
       if(proc0.and.debug) print*,'Taking step with size',time_step,'at time',internal_time
       do istage = 1, nstages
          call copy(g, gnew)

          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(istage_sub, iglo) &
          !$OMP SHARED(g_lo, gnew, time_step, scheme, stages, istage) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do istage_sub = 1, istage-1
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                gnew(:, :, iglo) = gnew(:, :, iglo) + (0.5 * time_step * scheme%coeffs(istage_sub, istage)) * stages(:, :, iglo, istage_sub)
             end do
          end do
          !$OMP END PARALLEL DO
          call calculate_fields(fields_func, gnew, phinew, aparnew, bparnew)

          ! Get the current nonlinear source term
          call source_func(gnew, gexp_1, &
               phinew, aparnew, bparnew, &
               max_vel, need_to_adjust = .not. advance_nonadiabatic_dfn, &
               calculate_cfl_limit = (istage == 1))

          nrhs_calls = nrhs_calls + 1
          if (istage == 1) errors(3) = max_vel

          ! Could we just pass stages(:, :, :, istage) instead of gexp_1
          ! in above source_func call?
          call copy(gexp_1, stages(:, :, :, istage))
       end do

       ! One can directly form the error by simply forming error =
       ! time_step * 0.5 * sum(stages_i * (high_order_coeffs_i -
       ! low_order_coeffs_i)). Here we store this in gexp_1
       ! without the time_step * 0.5 factor (to avoid extra
       ! array operations). We multiply the final error by
       ! time_step/2 to ensure full consistency.
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo) &
       !$OMP SHARED(g_lo, gexp_1, error_coeffs, stages) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          gexp_1(:, :, iglo) = error_coeffs(1) * stages(:, :, iglo, 1)
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, istage) &
       !$OMP SHARED(g_lo, gexp_1, error_coeffs, nstages, stages) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do istage = 2, nstages
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             gexp_1(:, :, iglo) = gexp_1(:, :, iglo) + error_coeffs(istage) * stages(:, :, iglo, istage)
          end do
       end do
       !$OMP END PARALLEL DO

       ! Store the absolute error
       errors(1) = gs2_max_abs(gexp_1) * (time_step * 0.5)
       ! Estimate the absolution size of the solution using
       ! the old solution
       errors(2) = gs2_max_abs(g)
       call max_allreduce(errors)
       inverse_error = get_inverse_error(errors)
       cfl_time_step = 1.0 / errors(3)

       !Calculate the new_time_step
       new_time_step = min(time_step * (inverse_error)**(1.0/scheme%order), cfl_time_step)
       cfl_limited = new_time_step < time_step .and. new_time_step == cfl_time_step
       new_time_step = new_time_step * time_step_safety_factor

       if(proc0 .and. debug) print*,'Inverse error is ',inverse_error,'at internal_time', &
            internal_time,'with step',time_step,'new time step is',new_time_step, &
            'counter=',counter

       ! If the cfl time step is too small then abort
       if (new_time_step < code_dt_min) &
            call mp_abort("Time step has become too small in advance_nonlinear_term.", .true.)

       ! If the error is too large we need to retry with a smaller step, so don't
       ! update anything
       if (inverse_error < 1.0 .or. (cfl_limited .and. immediate_reset)) then
          if(proc0.and.debug) print*,'Time step with size',time_step, &
               'failed at internal_time',internal_time,'and counter',counter, &
               'retrying with step',new_time_step
          nfailed_steps = nfailed_steps + 1
       else
          !Update the solution. Note we could probably merge the two
          !OpenMP loops into a single OpenMP parallel region and then
          !just add OMP DO to each loop. This would save a bit of overhead.
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo) &
          !$OMP SHARED(g_lo, g, solution_coeffs, time_step, stages) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             g(:, :, iglo) = g(:, :, iglo) + (solution_coeffs(1) * time_step * 0.5) * stages(:, :, iglo, 1)
          end do
          !$OMP END PARALLEL DO

          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(istage, iglo) &
          !$OMP SHARED(g_lo, g, solution_coeffs, nstages, time_step, stages) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do istage = 2, nstages
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                g(:, :, iglo) = g(:, :, iglo) + (solution_coeffs(istage) * time_step * 0.5) * stages(:, :, iglo, istage)
             end do
          end do
          !$OMP END PARALLEL DO

          internal_time = internal_time + time_step
          nsuccessful_steps = nsuccessful_steps + 1
       end if

       !Now update the time step, capped to dt
       time_step = min(new_time_step, dt)

       ! Make sure we don't overshoot the target time by limiting time step.
       if(internal_time < dt) then
          time_step = min(time_step, dt - internal_time)
          counter = counter + 1
       end if

    end do

    if(proc0.and.debug) print*,'Completed rk step in',counter,'steps to',internal_time

    ! We may be able to avoid this copy and field solve on exit if we restructure the main
    ! loop slightly so gnew represents the new solution. This may require use of g_work
    ! or similar to hold intermediate values.
    call copy(g, gnew)

    time_step_last = time_step
  end subroutine advance_nonlinear_term_rk_implementation

  !> Returns the normalised inverse error estimate. In other words
  !> the error tolerance, rtol*errors(2) + atol, divided by the error
  !> estimate (plus a small number to avoid divide by zero).
  !>
  !> Assumes that errors contains the global magnitude of the error
  !> estimate in the first element and the global maximum of the solution
  !> in the second element.
  real pure function get_inverse_error(errors) result(inverse_error)
    implicit none
    real, dimension(:), intent(in) :: errors
    inverse_error = (relative_tolerance * errors(2) + absolute_tolerance) / &
         (errors(1) + epsilon(errors(1)))
  end function get_inverse_error

  !> Calculates the fields consistent with the passed distribution function
  subroutine calculate_fields(fields_func, g_in, phi_out, apar_out, bpar_out)
    use nonlinear_terms, only: time_add_explicit_terms_field
    use job_manage, only: time_message
    procedure(invert_field_func) :: fields_func
    complex, dimension(:, :, :), intent(in) :: g_in
    complex, dimension (:,:,:), intent (out) :: phi_out, apar_out, bpar_out
    call time_message(.false., time_add_explicit_terms_field, 'Explicit terms - field')
    call fields_func(g_in, phi_out, apar_out, bpar_out, local_only = .true.)
    call time_message(.false., time_add_explicit_terms_field, 'Explicit terms - field')
  end subroutine calculate_fields

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
    logical :: advance_nonadiabatic_dfn, show_statistics
    real :: absolute_tolerance, convergence_tolerance, relative_tolerance, time_step_safety_factor

    namelist /split_nonlinear_terms_knobs/ absolute_tolerance, advance_nonadiabatic_dfn, convergence_tolerance, relative_tolerance, rk_method, show_statistics, split_method, time_step_safety_factor

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    absolute_tolerance = self%absolute_tolerance
    advance_nonadiabatic_dfn = self%advance_nonadiabatic_dfn
    convergence_tolerance = self%convergence_tolerance
    relative_tolerance = self%relative_tolerance
    rk_method = self%rk_method
    show_statistics = self%show_statistics
    split_method = self%split_method
    time_step_safety_factor = self%time_step_safety_factor

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = split_nonlinear_terms_knobs)

    ! Now copy from local variables into type members
    self%absolute_tolerance = absolute_tolerance
    self%advance_nonadiabatic_dfn = advance_nonadiabatic_dfn
    self%convergence_tolerance = convergence_tolerance
    self%relative_tolerance = relative_tolerance
    self%rk_method = rk_method
    self%show_statistics = show_statistics
    self%split_method = split_method
    self%time_step_safety_factor = time_step_safety_factor

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
    call self%write_key_val("absolute_tolerance", self%absolute_tolerance, unit_internal)
    call self%write_key_val("advance_nonadiabatic_dfn", self%advance_nonadiabatic_dfn, unit_internal)
    call self%write_key_val("convergence_tolerance", self%convergence_tolerance, unit_internal)
    call self%write_key_val("relative_tolerance", self%relative_tolerance, unit_internal)
    call self%write_key_val("rk_method", self%rk_method, unit_internal)
    call self%write_key_val("show_statistics", self%show_statistics, unit_internal)
    call self%write_key_val("split_method", self%split_method, unit_internal)
    call self%write_key_val("time_step_safety_factor", self%time_step_safety_factor, unit_internal)
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
    call broadcast(self%absolute_tolerance)
    call broadcast(self%advance_nonadiabatic_dfn)
    call broadcast(self%convergence_tolerance)
    call broadcast(self%relative_tolerance)
    call broadcast(self%rk_method)
    call broadcast(self%show_statistics)
    call broadcast(self%split_method)
    call broadcast(self%time_step_safety_factor)

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
