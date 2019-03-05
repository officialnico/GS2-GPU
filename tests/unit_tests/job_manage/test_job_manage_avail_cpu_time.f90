!> This unit runs a tiny grid with nstep=20000, and checks
!! that GS2 exits when run time exceeds avail_time.
!
program test_job_manage
   use gs2_main, only: run_gs2
   use gs2_main, only: old_iface_state, finalize_gs2
   use gs2_main, only: finalize_equations
   use unit_tests
   use mp, only: init_mp, mp_comm, finish_mp
   use job_manage, only: timer_local
   implicit none
   real :: start_time
   real :: avail_time = 2.0 ! expected time for this test
   real :: elapsed_time

   ! Initialize
   call init_mp

   ! Tell make system to run this test
   functional_test_flag = .true.

   call announce_module_test("job_manage")

   ! Run test:
   ! 1. Run a long job that will take more than avail_time = 2 secs to complete,
   !    and should therefore be aborted.
   ! 2. Run gs2 and at the end check that the job was aborted, i.e. that:
   !    a) more than avail_time has elapsed
   !    b) fewer than nstep timesteps have been completed
   call announce_test("that gs2 stops when time elaspsed exceeds avail_cpu_time")
   start_time = timer_local()
   call run_gs2(mp_comm)

   elapsed_time = timer_local() - start_time

   call process_test(( old_iface_state%istep_end .lt. old_iface_state%nstep ), &
     "gs2 doesn't stop when time elapsed exceeds avail_cpu_time -- too many steps")

   
   call process_test( &
     ( elapsed_time .gt. avail_time ), &
     "gs2 doesn't stop when time elapsed exceeds avail_cpu_time -- stopped too early")

   ! Finalize gs2
   call finalize_equations(old_iface_state)
   call finalize_gs2(old_iface_state)

   call close_module_test("job_manage")

   call finish_mp


end program test_job_manage
