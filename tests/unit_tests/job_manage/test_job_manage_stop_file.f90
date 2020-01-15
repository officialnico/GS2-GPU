!> This unit runs a tiny grid with nstep=10000, and checks
!! that the creation of a .stop file causes GS2 to exit.
!
program test_job_manage
   use gs2_main, only: run_gs2
   use gs2_main, only: old_iface_state, finalize_gs2
   use gs2_main, only: finalize_equations
   use unit_tests
   use mp, only: init_mp, mp_comm, proc0, finish_mp
   use job_manage, only: timer_local
   use exit_codes, only: EXIT_STOP_FILE
   implicit none

   ! Initialize
   call init_mp

   ! Tell make system to run this test
   functional_test_flag = .true.

   call announce_module_test("job_manage")

   ! Run test:
   ! 1. Create stop file
   ! 2. Run gs2
   ! 3. Read gs2's exit_reason log to check that a stop file was detected
   if (proc0) call system("touch test_job_manage.stop")
   call announce_test("that gs2 doesn't run without limit when run_name.stop is present")
   call run_gs2(mp_comm)

   ! The reason for stopping should be given in "test_job_manage.exit_reason".
   if(proc0) call process_test((get_exit_code() == EXIT_STOP_FILE%code), &
     "gs2 doesn't run without limit when run_name.stop is present")

   ! Finalize gs2
   call finalize_equations(old_iface_state)
   call finalize_gs2(old_iface_state)

   call close_module_test("job_manage")

   call finish_mp

   ! Delete stop file
   if (proc0) call system("rm test_job_manage.stop")

   contains

     ! Read error code from the file "run_name.exit_reason"
     function get_exit_code()
       use file_utils, only: get_unused_unit
       implicit none
       integer get_exit_code
       integer unit

       call get_unused_unit(unit)

       open (unit=unit, file="test_job_manage.exit_reason")
       read (unit, '(I3.1)') get_exit_code
       close (unit)

     end function get_exit_code

end program test_job_manage
