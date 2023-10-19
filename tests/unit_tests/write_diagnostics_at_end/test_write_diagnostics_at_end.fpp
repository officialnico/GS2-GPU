program test_write_diagnostics_at_end
  use gs2_main, only: run_gs2, gs2_program_state_type
  use unit_tests, only: announce_test, process_test
  use mp, only: init_mp, mp_comm, finish_mp
  implicit none
  integer, parameter :: nout_expected = 2
  character(len=*), parameter :: test_name = 'Write diagnostics at end'
  type(gs2_program_state_type) :: state

  call init_mp()
  state%mp_comm = mp_comm
  call run_gs2(state)

  call announce_test(test_name)
  call process_test(check_number_of_time_points(), test_name)
  call finish_mp()
    
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

