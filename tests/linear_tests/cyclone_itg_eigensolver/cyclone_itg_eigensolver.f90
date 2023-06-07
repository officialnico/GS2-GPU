module checks_mod
  use functional_tests, only: check_growth_rate, check_frequency
  public checks
  contains
    logical function checks()
      use functional_tests, only: check_growth_rate, check_frequency
      implicit none
      real, parameter :: tolerance = 1.0e-4
      real, parameter :: omega_expected = 0.65367777
      real, parameter :: gamma_expected = 0.16566611
      checks = check_growth_rate([gamma_expected], tolerance)
      checks = checks .and. check_frequency([omega_expected], tolerance)
    end function checks
end module checks_mod

program test_cyclone_itg_eigensolver
  use functional_tests, only: test_gs2
  use checks_mod, only: checks
  use mp, only: init_mp, proc0, finish_mp
  implicit none

  call init_mp
  call test_gs2('cyclone itg eigensolver', checks)
  call finish_mp

end program test_cyclone_itg_eigensolver
