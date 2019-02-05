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
   implicit none
   real :: start_time
   real :: avail_time = 2.0

   ! Initialize
   call init_mp

   ! Tell make system to run this test
   functional_test_flag = .true.

   call announce_module_test("job_manage")

   ! Run test:
   ! 1. Create stop file
   ! 2. Run gs2
   ! 3. Check gs2 exits before the available cpu time expires
   if (proc0) call system("touch test_job_manage.stop")
   call announce_test("that gs2 doesn't run without limit when run_name.stop is present")
   start_time = timer_local()
   call run_gs2(mp_comm)

   ! If gs2 was aborted by the stop file, the elapsed time (timer_local - start_time)
   ! will be less than avail_time.
   call process_test((timer_local() - start_time < avail_time), &
     "gs2 doesn't run without limit when run_name.stop is present")

   ! Finalize gs2
   call finalize_equations(old_iface_state)
   call finalize_gs2(old_iface_state)

   call close_module_test("job_manage")

   call finish_mp

   ! Delete stop file
   if (proc0) call system("rm test_job_manage.stop")


end program test_job_manage
