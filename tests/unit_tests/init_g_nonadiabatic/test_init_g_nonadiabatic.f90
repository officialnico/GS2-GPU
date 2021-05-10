program test_init_g_nonadiabatic
  use gs2_init, only: init, init_level_list
  use unit_tests
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use mp, only: init_mp, finish_mp, proc0, mp_comm
  use iso_fortran_env, only: output_unit
  use dist_fn, only: calculate_potentials_from_nonadiabatic_dfn, get_init_field
  use dist_fn_arrays, only: g, gnew, g_adjust
  use fields_arrays, only: phinew, aparnew, bparnew, phi, apar, bpar
  use fields, only: set_init_fields
  use init_g, only: ginit
  use init_g, only: ginitopt_restart_many
  use run_parameters, only: fphi, fbpar
  implicit none
  type(gs2_program_state_type) :: state
  logical :: dummy, agrees
  real, parameter :: default_tolerance = 1.0e-10
  
  ! Set up depenencies
  call init_mp

  call announce_module_test('test_init_g_nonadiabatic')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  ! Might expect to initialize up to but not including setting the
  ! initial values but actually we stop just before fields_level_2 where
  ! the response matrices are calculated as we do not need this level and
  ! this can save a significant amount of time.
  call init(state%init, init_level_list%fields_level_2 - 1)

  agrees = .true.
  agrees = test_input_condition() .and. agrees
  agrees = test_input_condition_finite_apar() .and. agrees  
  agrees = test_zero() .and. agrees
  agrees = test_zero_charge_density() .and. agrees    
  
  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call process_test(agrees,'test_init_g_nonadiabatic')
  
  call close_module_test('test_init_g_nonadiabatic')

  call finish_mp

contains
  
  !> Test initialising with the default setup from the input file
  logical function test_input_condition() result(agrees)
    implicit none
    character(len=*), parameter :: case_name = 'input_condition'
    logical, dimension(4) :: checks
    ! Use the initial condition set in the input file
    call ginit(dummy)
    call get_init_field(phinew, aparnew, bparnew)
    phi = phinew ; apar = aparnew ; bpar = bparnew
    
    ! Shift to nonadiabatic distribution function, calculate potentials
    ! directly and then shift back
    call g_adjust(gnew, phinew, bparnew, fphi, fbpar)
    call calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)
    call g_adjust(gnew, phinew, bparnew, -fphi, -fbpar)

    ! We choose to store all the results and then check in one go
    ! to avoid any possible logical short-circuiting which could
    ! lead to some checks being skipped.
    
    checks = [ &
         check_values(g, gnew, case_name, 'g'), &
         check_values(phi, phinew, case_name, 'phi'), &
         check_values(apar, aparnew, case_name, 'apar'), &
         check_values(bpar, bparnew, case_name, 'bpar') &
         ]

    agrees = all(checks)
    call reset_state
  end function test_input_condition

  !> Test initialising with the default setup from the input file
  !> but with initial apar non-zero
  logical function test_input_condition_finite_apar() result(agrees)
    implicit none
    character(len=*), parameter :: case_name = 'input_condition_finite_apar'
    logical, dimension(4) :: checks
    ! Use the initial condition set in the input file
    call ginit(dummy)
    g(:, 2, :) = -g(:, 2, :)
    gnew = g
    call get_init_field(phinew, aparnew, bparnew)
    phi = phinew ; apar = aparnew ; bpar = bparnew
    
    ! Shift to nonadiabatic distribution function, calculate potentials
    ! directly and then shift back
    call g_adjust(gnew, phinew, bparnew, fphi, fbpar)
    call calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)
    call g_adjust(gnew, phinew, bparnew, -fphi, -fbpar)

    ! We choose to store all the results and then check in one go
    ! to avoid any possible logical short-circuiting which could
    ! lead to some checks being skipped.
    
    checks = [ &
         check_values(g, gnew, case_name, 'g'), &
         check_values(phi, phinew, case_name, 'phi'), &
         check_values(apar, aparnew, case_name, 'apar'), &
         check_values(bpar, bparnew, case_name, 'bpar') &
         ]

    agrees = all(checks)
    call reset_state
  end function test_input_condition_finite_apar
  
  !> Test fields with nonadiabatic g = 0
  logical function test_zero() result(agrees)
    implicit none
    character(len=*), parameter :: case_name = 'zero'
    logical, dimension(3) :: checks
    g = 0 ; gnew = 0.
    phi = 0. ; apar = 0. ; bpar = 0.
    phinew = 0. ; aparnew = 0. ; bparnew = 0.    
    call calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)    
    
    checks = [ &
         check_values(phi, phinew, case_name, 'phi'), &
         check_values(apar, aparnew, case_name, 'apar'), &
         check_values(bpar, bparnew, case_name, 'bpar') &
         ]

    agrees = all(checks)
    call reset_state
  end function test_zero

  !> Test fields with zero charge density.
  !> Note we don't test bpar here as this should be non-zero here.
  logical function test_zero_charge_density() result(agrees)
    use gs2_layouts, only: g_lo, is_idx
    use species, only: spec
    implicit none
    character(len=*), parameter :: case_name = 'zero_charge_density'
    logical, dimension(2) :: checks
    integer :: iglo, is

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       gnew(:, :, iglo) = spec(is)%dens
    end do
    
    phi = 0. ; apar = 0. ; bpar = 0. !< This is what we expect to get
    phinew = 0. ; aparnew = 0. ; bparnew = 0.
    
    call calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)    
    
    checks = [ &
         check_values(phi, phinew, case_name, 'phi'), &
         check_values(apar, aparnew, case_name, 'apar') &
         ]

    agrees = all(checks)
    call reset_state
  end function test_zero_charge_density
  
  subroutine reset_state
    implicit none
    g = 0. ; phi =  0. ; apar = 0. ; bpar = 0.
    gnew = g ; phinew = phi ; aparnew = apar ; bparnew = bpar
  end subroutine reset_state

  logical function check_values(array_a, array_b, &
       case_name, field_name, tolerance_in) result(agrees)
    implicit none
    complex, dimension(:, :, :), intent(in) :: array_a, array_b
    character(len=*), intent(in) :: case_name, field_name
    real, intent(in), optional :: tolerance_in
    real :: tolerance
    real :: max_diff
    logical :: dummy_flag
    tolerance = default_tolerance
    if(present(tolerance_in)) tolerance = tolerance_in

    max_diff = maxval(abs(array_a-array_b))

    dummy_flag = .true.
    call announce_test('init_g_nonadiabatic :: '//case_name//":"//field_name)

    agrees = max_diff <= tolerance
    call process_check(dummy_flag, agrees, case_name//":"//field_name)
    if (.not.agrees .and. proc0) then
       write(output_unit, *) case_name//":"//field_name,max_diff
    end if
    
  end function check_values
  
end program test_init_g_nonadiabatic
