!> A small test for the ballstab is_unstable method underpinning ideal_ball
!> This is really a golden answer test checking that certain cases return
!> the expected stability status.
program test_ideal_ball
  use unit_tests
  use ballstab, only: run_stability_check, init_ballstab, finish_ballstab, stability
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, run_name, input_unit_exist
  use theta_grid, only: finish_theta_grid
  use constants, only: run_name_size
  implicit none
  logical :: list, exist, passes
  integer :: in_file, n_shat, n_betaprime
  character(len=run_name_size), target :: cbuff
  logical, dimension(:, :), allocatable :: expect_unstable
  namelist/ideal_ball_test/expect_unstable

  call init_mp
  if (proc0) then
     call init_file_utils(list, name = 'gs')
     cbuff = trim(run_name)
  end if
  call broadcast(cbuff)
  run_name => cbuff
  call announce_module_test('ideal_ball')
  passes = .false.

  call init_ballstab

  ! Ballooning stability currently only done on proc0
  ! so we can only run the test on proc0
  if (proc0) then
     ! Work out how many shear and betaprime values are in use
     n_shat = size(stability, dim = 1)
     n_betaprime = size(stability, dim = 2)

     allocate(expect_unstable(n_shat, n_betaprime))

     ! Read in expected result
     expect_unstable = .false.
     in_file = input_unit_exist("ideal_ball_test", exist)
     if (exist) read(in_file, nml = ideal_ball_test)

     ! Calculate stability
     call run_stability_check
     passes = all(expect_unstable .eqv. (stability > 0))
  end if

  call broadcast(passes)
  call process_test(passes, 'ideal_ball')
  call close_module_test('ideal_ball')

  call finish_ballstab
  call finish_theta_grid
  call finish_mp
end program test_ideal_ball
