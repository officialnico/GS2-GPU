!> A test program to run a linear cyclone itg benchmark
!> at low resolution and check the growth rate
module checks_mod
  implicit none
  real :: phi2_value
  real, parameter :: target_value = 2200.0
  real, parameter :: tolerance = 0.05
contains
  real function phi2()
    use fields_arrays, only: phinew
    use theta_grid, only: delthet, jacob
    use kt_grids, only: ntheta0, naky, aky
    implicit none
    real :: anorm, fac
    integer :: ik, it
    anorm = sum(delthet*jacob)
    phi2 = 0
    do ik = 1, naky
       fac = 0.5
       if ( aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
          phi2 = phi2 + sum(phinew(:, it, ik) * conjg(phinew(:, it, ik)) *delthet*jacob)*fac/anorm
       end do
    end do

    phi2_value = phi2
  end function phi2
  
  function checks()
    use unit_tests, only: agrees_with
    implicit none
    logical :: checks
    ! Horrible hack to keep value available after test_gs2 has finished
    phi2_value = phi2()
    checks =  agrees_with(phi2_value,target_value,tolerance)
  end function checks
end module checks_mod

program cyclone_itg_low_res
  use functional_tests, only: test_gs2
  use checks_mod, only: checks, phi2_value, target_value
  use fields_local, only: fields_local_functional
  use mp, only: init_mp, proc0, finish_mp
  use gs2_time, only: code_dt_cfl
  implicit none
  ! Note the field type is currently determined by the input file
  ! so we will need to change this if we change [[fields_knobs:field_option]]
  ! in the future.
  integer :: length, ierr

  ! If you want to see the minimal linear test look at slab_itg_low_res
  call init_mp
  call test_gs2('Nonline CBC ITG with time step changes', checks)
  if(proc0) print*,phi2_value," ",target_value
  call finish_mp
end program cyclone_itg_low_res
