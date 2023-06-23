!> A test program to test that changes in the choice of normalisation don't
!> alter results when correcting for normalisations.
program test_cyclone_normalisations
  use functional_tests, only: test_gs2
  use mp, only: init_mp, finish_mp
  implicit none

  call init_mp
  call test_gs2('cyclone itg normalisation cases', checks)
  call finish_mp
  
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
  
end program test_cyclone_normalisations

