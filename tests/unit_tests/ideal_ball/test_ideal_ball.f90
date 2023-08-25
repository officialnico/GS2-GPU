!> A small test for the ballstab is_unstable method underpinning ideal_ball
!> This is really a golden answer test checking that certain cases return
!> the expected stability status.
program test_ideal_ball
  use unit_tests
  use ballstab, only: is_unstable, init_ballstab, finish_ballstab
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, run_name, input_unit_exist
  use theta_grid, only: finish_theta_grid
  use constants, only: run_name_size
  implicit none
  logical :: list, exist, expect_unstable, passes, is_actually_unstable
  integer :: in_file
  character(len=run_name_size), target :: cbuff
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

  ! Read in expected result
  expect_unstable = .false.
  if (proc0) then
     in_file = input_unit_exist("ideal_ball_test", exist)
     if (exist) read(in_file, nml = ideal_ball_test)
  end if
  call broadcast(expect_unstable)

  call init_ballstab

  if (proc0) then
     is_actually_unstable = is_unstable() > 0
     passes = (expect_unstable .eqv. is_actually_unstable)
     if (.not. passes) then
        print*,'Fail: Expected ',expect_unstable,'but found',is_actually_unstable
     end if

  end if

  call broadcast(passes)
  call process_test(passes, 'ideal_ball')
  call close_module_test('ideal_ball')

  call finish_ballstab
  call finish_theta_grid
  call finish_mp
end program test_ideal_ball
