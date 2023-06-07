!> A test program to run a linear case using the grid file
!> approach to geometry. Allows us to test more complex geometry
!> such as those with local bmag minima/maxima.
!> Thanks to Jason Parisi for providing the initial inputs
!> that this test is derived from.

module checks_mod
  use functional_tests, only: check_growth_rate, check_frequency
  public checks
  contains
    logical function checks()
      use functional_tests, only: check_growth_rate, check_frequency
      implicit none
      real, parameter :: tolerance = 1.0e-2
      real, parameter :: omega_expected = 0.723
      real, parameter :: gamma_expected = 0.286
      checks = check_growth_rate([gamma_expected], tolerance)
      checks = checks .and. check_frequency([omega_expected], tolerance)
    end function checks
end module checks_mod

program test_grid_file_local_bmag_minima
  use functional_tests, only: test_gs2
  use checks_mod, only: checks
  use mp, only: init_mp, proc0, finish_mp
  implicit none

  call init_mp
  call test_gs2('grid file with local bmag minima', checks)
  call finish_mp

end program test_grid_file_local_bmag_minima
