module check_mod
  implicit none
contains
  logical function checks()
    use file_utils, only: input_unit_exist
    use functional_tests, only: check_growth_rate, check_frequency
    implicit none
    real :: tolerance, gamma_expected, omega_expected
    logical :: exist
    integer :: in_file
    namelist /test_settings/ tolerance, gamma_expected, omega_expected

    tolerance = 1.0e-6
    gamma_expected = 0.0
    omega_expected = 0.0   

    in_file = input_unit_exist("test_settings", exist)
    if (exist) read(in_file, nml = test_settings)
    checks = check_growth_rate([gamma_expected], tolerance)
    checks = checks .and. check_frequency([omega_expected], tolerance)
  end function checks
end module check_mod

!> A test program to test the impact of collisions on the cyclone base case.
!> This is essentially just a golden answer test for now.
program test_cyclone_itg_collisions
  use functional_tests, only: test_gs2
  use mp, only: init_mp, finish_mp
  use check_mod, only: checks
  implicit none

  call init_mp
  call test_gs2('cyclone itg collisions cases', checks)
  call finish_mp
end program test_cyclone_itg_collisions
