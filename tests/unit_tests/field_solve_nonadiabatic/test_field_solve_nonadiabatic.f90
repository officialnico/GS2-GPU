program test_field_solve_nonadiabatic
  use gs2_init, only: init, init_level_list
  use unit_tests
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use mp, only: init_mp, finish_mp, proc0, mp_comm
  use iso_fortran_env, only: output_unit
  use dist_fn, only: calculate_potentials_from_nonadiabatic_dfn, get_init_field
  use dist_fn_arrays, only: g, gnew, g_adjust
  use fields_arrays, only: phinew, aparnew, bparnew, phi, apar, bpar
  use run_parameters, only: fphi, fbpar
  use init_g, only: ginit

  implicit none
  type(gs2_program_state_type) :: state
  logical :: dummy, agrees
  real, parameter :: default_tolerance = 1.0e-10
  complex, dimension(:, :, :), allocatable :: phi_init, apar_init, bpar_init
  complex, dimension(:, :, :), allocatable :: phi_recalc, apar_recalc, bpar_recalc
  logical, parameter :: debug = .false.
  ! Set up depenencies
  call init_mp

  call announce_module_test('test_field_solve_nonadiabatic')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  ! Might expect to initialize up to but not including setting the
  ! initial values but actually we stop just before fields_level_2 where
  ! the response matrices are calculated as we do not need this level and
  ! this can save a significant amount of time.
  call init(state%init, init_level_list%fields_level_2 - 1)

  agrees = .true.
  gnew = 1.0
  g = gnew
  phi_init = phinew; apar_init = aparnew; bpar_init = bparnew;
  phi_recalc = phinew; apar_recalc = aparnew; bpar_recalc = bparnew;

  call get_init_field(phi_init, apar_init, bpar_init)
  call g_adjust(gnew, phi_init, bpar_init, fphi, fbpar)
  call calculate_potentials_from_nonadiabatic_dfn(gnew, phi_recalc, apar_recalc, bpar_recalc)
  if (debug .and. proc0) then
     print*,maxval(abs(phi_init)),maxval(abs(phi_recalc)),maxval(abs(phinew))
     print*,maxval(abs(apar_init)),maxval(abs(apar_recalc)),maxval(abs(aparnew))
     print*,maxval(abs(bpar_init)),maxval(abs(bpar_recalc)),maxval(abs(bparnew))
     print*,''
  end if

  agrees = test_fields('default') .and. agrees

  gnew = 1.0
  gnew(:, 2, :) = 0.0
  g = gnew
  phi_init = phinew; apar_init = aparnew; bpar_init = bparnew;
  phi_recalc = phinew; apar_recalc = aparnew; bpar_recalc = bparnew;

  call get_init_field(phi_init, apar_init, bpar_init)
  call g_adjust(gnew, phi_init, bpar_init, fphi, fbpar)
  call calculate_potentials_from_nonadiabatic_dfn(gnew, phi_recalc, apar_recalc, bpar_recalc)
  if (debug .and. proc0) then
     print*,maxval(abs(phi_init)),maxval(abs(phi_recalc)),maxval(abs(phinew))
     print*,maxval(abs(apar_init)),maxval(abs(apar_recalc)),maxval(abs(aparnew))
     print*,maxval(abs(bpar_init)),maxval(abs(bpar_recalc)),maxval(abs(bparnew))
     print*,''
  end if

  agrees = test_fields('finite apar') .and. agrees

  call ginit(dummy)
  g = gnew
  phi_init = phinew; apar_init = aparnew; bpar_init = bparnew;
  phi_recalc = phinew; apar_recalc = aparnew; bpar_recalc = bparnew;

  call get_init_field(phi_init, apar_init, bpar_init)
  call g_adjust(gnew, phi_init, bpar_init, fphi, fbpar)
  call calculate_potentials_from_nonadiabatic_dfn(gnew, phi_recalc, apar_recalc, bpar_recalc)
  if (debug .and. proc0) then
     print*,maxval(abs(phi_init)),maxval(abs(phi_recalc)),maxval(abs(phinew))
     print*,maxval(abs(apar_init)),maxval(abs(apar_recalc)),maxval(abs(aparnew))
     print*,maxval(abs(bpar_init)),maxval(abs(bpar_recalc)),maxval(abs(bparnew))
     print*,''
  end if

  agrees = test_fields('init_g') .and. agrees

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call process_test(agrees,'test_field_solve_nonadiabatic')

  call close_module_test('test_field_solve_nonadiabatic')

  call finish_mp

contains

  logical function test_fields(case_name) result(agrees)
    implicit none
    character(len = *), intent(in) :: case_name
    logical, dimension(:), allocatable :: checks
    checks = [&
         check_values(phi_init, phi_recalc, case_name, 'phi'), &
         check_values(apar_init, apar_recalc, case_name, 'apar'), &
         check_values(bpar_init, bpar_recalc, case_name, 'bpar')  &
         ]
    agrees = all(checks)
  end function test_fields

  logical function check_values(array_a, array_b, &
       case_name, field_name, tolerance_in) result(agrees)
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension(:, :, :), intent(in) :: array_a, array_b
    character(len=*), intent(in) :: case_name, field_name
    real, intent(in), optional :: tolerance_in
    real :: tolerance
    real :: max_diff
    logical :: dummy_flag
    tolerance = get_option_with_default(tolerance_in, default_tolerance)
    max_diff = maxval(abs(array_a-array_b))/(maxval(abs(array_a))+1.0e-16)

    dummy_flag = .true.
    call announce_test('field_solve_nonadiabatic :: '//case_name//":"//field_name)

    agrees = max_diff <= tolerance
    call process_check(dummy_flag, agrees, case_name//":"//field_name)
    if (.not.agrees .and. proc0) then
       write(output_unit, *) case_name//":"//field_name,max_diff
    end if

  end function check_values

end program test_field_solve_nonadiabatic
