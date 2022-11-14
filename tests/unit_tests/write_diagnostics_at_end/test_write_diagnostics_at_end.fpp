program test_write_diagnostics_at_end
  use gs2_main, only: run_gs2, finish_gs2, old_iface_state
  use gs2_main, only: finalize_diagnostics, finalize_equations, finalize_gs2
  use unit_tests, only: functional_test_flag, ilast_step, announce_test, process_test
  use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
  use gs2_diagnostics, only: finish_gs2_diagnostics
  implicit none
  integer, parameter :: nout_expected = 2
  character(len=*), parameter :: test_name = 'Write diagnostics at end'

  call init_mp

  test_driver_flag = .true.
  functional_test_flag = .true.

  call run_gs2(mp_comm)
  call finalize_diagnostics(old_iface_state)

  call announce_test('results')
  call process_test(check_number_of_time_points(), 'results')

  call finalize_equations(old_iface_state)
  call finalize_gs2(old_iface_state, quiet=.true.)

  call finish_mp
    
contains

  logical function check_number_of_time_points() result(check)
    use file_utils, only: run_name
    implicit none
    integer :: exitstat
    character(len=1) :: str_expect
    write(str_expect, '(I0)') nout_expected
    ! Check the netcdf file has the expected number of time outputs in it
    call execute_command_line( &
         "THEVAL=$(ncdump -h "//trim(run_name)//".out.nc | "// &
         "grep UNLIM | awk -F '(' '{print $2}' | awk '{print $1}') ; "// &
         "[[ ${THEVAL} == "//str_expect//" ]]", exitstat = exitstat)
    check = exitstat == 0
  end function check_number_of_time_points
  
end program test_write_diagnostics_at_end

