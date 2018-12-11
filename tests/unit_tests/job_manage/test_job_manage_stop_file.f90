!> This unit runs a tiny grid with nstep=10000, and checks
!! that the .stop file mechanism is working
!
program test_job_manage
   use gs2_main, only: run_gs2
   use gs2_main, only: old_iface_state, finalize_gs2, finalize_diagnostics
   use gs2_main, only: finalize_equations
   use unit_tests
   use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
   use job_manage, only: timer_local
   implicit none
   real :: start_time
   real :: avail_time = 2.0

   ! Initialize
   call init_mp

   test_driver_flag = .true.
   functional_test_flag = .true.

   call announce_module_test("job_manage")

   ! Run test
   if (proc0) call system("touch test_job_manage.stop")
   call announce_test("that gs2 doesn't run without limit when run_name.stop is present")
   start_time = timer_local()
   call run_gs2(mp_comm)
   call process_test((timer_local() - start_time < avail_time), &
     "that gs2 doesn't run without limit when run_name.stop is present")

   ! Finalize
   call finalize_diagnostics(old_iface_state)
   call finalize_equations(old_iface_state)
   call finalize_gs2(old_iface_state)

   call close_module_test("job_manage")

   call finish_mp

   if (proc0) call system("rm test_job_manage.stop")

contains

end program test_job_manage
