program test_split_nonlinear_integrators
  use iso_fortran_env, only: output_unit
  use unit_tests
  use gs2_main, only: gs2_program_state_type
  use mp, only: init_mp, finish_mp
  use mp, only: proc0, max_allreduce
  use gs2_main, only: initialize_wall_clock_timer
  use gs2_main, only: initialize_gs2
  use gs2_main, only: initialize_equations
  use gs2_main, only: finalize_equations
  use gs2_main, only: finalize_gs2
  use split_nonlinear_terms, only: advance_nonlinear_term
  use dist_fn_arrays, only: gnew, g
  use gs2_time, only: code_dt
  use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
  use run_parameters, only: nstep
  implicit none
  type(gs2_program_state_type) :: state
  logical :: agrees
  real, parameter :: tolerance = 1.0e-6
  complex :: exponential_coefficient
  real :: error
  character(len=40) :: name
  integer :: i

  call init_mp
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call initialize_equations(state)

  call announce_module_test('split_nonlinear_integrators')

  g = 1.0
  gnew = g
  phi = 0.0
  apar = 0.0
  bpar = 0.0
  phinew = phi
  aparnew = apar
  bparnew = bpar

  exponential_coefficient = 0.1
  do i = 1, nstep
     call advance_nonlinear_term(phi, apar, bpar, &
          null_field_op, exponential_source)
  end do

  ! Form expected result, gnew = g(t=0) * exp(exponential_coefficient * dt)
  g = exp(exponential_coefficient * code_dt * nstep)
  error = maxval(abs(gnew - g))
  call max_allreduce(error)

  ! Initialise local variables
  agrees = error <= tolerance

  name = "exponential growth"

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  call process_test(agrees,'split_nonlinear_integrators -- '//trim(name))

  call close_module_test('split_nonlinear_integrators')

  call finalize_equations(state)
  call finalize_gs2(state)
  call finish_mp

contains

     subroutine null_field_op(phi, apar, bpar, gf_lo)
       use theta_grid, only: ntgrid
       complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
       logical, optional :: gf_lo
     end subroutine null_field_op

     subroutine exponential_source(g_in, g1, phi, apar, bpar, dt_lim, calculate_cfl_limit)
       use nonlinear_terms, only: cfl
       complex, dimension (:,:,:), intent (in) :: g_in
       complex, dimension (:,:,:), intent (out) :: g1
       complex, dimension (:,:,:), intent (in) :: phi, apar, bpar
       real, intent(out) :: dt_lim
       logical, intent(in), optional :: calculate_cfl_limit
       ! Factor 2 to cancel hard-coded 0.5 in regular NL term
       g1 = g_in * exponential_coefficient * 2
       dt_lim = cfl
     end subroutine exponential_source

end program test_split_nonlinear_integrators
