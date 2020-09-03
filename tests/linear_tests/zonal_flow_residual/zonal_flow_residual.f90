!> A basic test of the zonal flow residual: "linear collisionless
!> processes do not damp poloidal flows driven by ITG turbulence"
!> Rosenbluth-Hinton 1998, https://doi.org/10.1103/PhysRevLett.80.724
!>
!> Zonal flows are linearly stable and are damped to a non-zero
!> residual. This tests that that residual level is correct. Actually,
!> follow Xiao-Catto 2007, https://doi.org/10.1063/1.2718519, who
!> expand to higher order in epsilon (inverse aspect ratio).
!>
!> This is not a complete, thorough test, just a fast sanity-check
program zonal_flow_residual
  use fields_arrays, only : phinew
  use gs2_main, only: initialize_gs2, initialize_equations, initialize_diagnostics, &
       & evolve_equations, initialize_wall_clock_timer, old_iface_state, &
       & finalize_diagnostics, finalize_equations, finalize_gs2
  use mp, only: init_mp, finish_mp
  use unit_tests, only: process_test, announce_test

  implicit none

  !> The initial value of phi
  real :: initial_phi

  call init_mp

  call initialize_wall_clock_timer
  call initialize_gs2(old_iface_state)
  call initialize_equations(old_iface_state)
  call initialize_diagnostics(old_iface_state)

  initial_phi = real(phinew(0, 1, 1))

  call evolve_equations(old_iface_state, old_iface_state%nstep)

  call announce_test('results')
  call process_test(check_residual(initial_phi), 'results')

  call finalize_diagnostics(old_iface_state)
  call finalize_equations(old_iface_state)
  call finalize_gs2(old_iface_state)

  call finish_mp

contains

  !> Check the long-time value of phi has been damped to the correct
  !> value, calculated from eq. 25 of Xiao-Catto 2007,
  !> https://doi.org/10.1063/1.2718519. Using epsilon = 0.02, q = 1.3,
  !> k_perp rho_i = 0.002
  !>
  !> Has a very loose tolerance due to needing low grid resolution in
  !> order to run in a few seconds
  function check_residual(initial_phi)
    use fields_arrays, only : phinew
    use unit_tests, only : agrees_with
    logical :: check_residual
    real, intent(in) :: initial_phi
    real :: actual_residual, current_phi
    real, parameter :: expected_residual = 0.046588479325594914
    real, parameter :: tolerance = 8.e-2

    current_phi = real(phinew(0, 1, 1))
    actual_residual = current_phi / initial_phi
    check_residual = agrees_with(actual_residual, expected_residual, tolerance)
  end function check_residual

end program zonal_flow_residual
