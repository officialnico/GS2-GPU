! Timing benchmark that repeatedly calls the local field solve

! Copyright Joseph Parker 2019

! This file is part of gs2.

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
! (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify,
! merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

program time_fields_local

  use unit_tests

  use benchmarks, only: benchmark_identifier
  use fields_arrays, only: phi, apar, bpar
  use fields_local, only: fields_local_testing, fields_local_functional
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_main, only: initialize_gs2, gs2_program_state_type, finalize_gs2
  use job_manage, only: time_message
  use mp, only: init_mp, finish_mp, broadcast, mp_comm, proc0, nproc
  use run_parameters, only: nstep

  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  integer :: istep
  integer :: timing_unit
  logical :: do_gather = .true.
  logical :: do_update = .true.

  !> Instance for access to private functions in fields_local
  type(fields_local_testing) :: fields_local_tests

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  if (proc0) call announce_module_test('time_fields_local')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)
  call init(state%init, init_level_list%fields_level_2)

  ! Run benchmark
  if (.not. fields_local_functional()) then

    if (proc0) write (*,*) "WARNING: fields_local is non-functional in your build. &
      & Skipping the fields_local unit test. &
      & If you are using the PGI compilers this is to be expected. "

  else 

    if (proc0) call announce_test('advance')

    if (proc0) call time_message(.false., time_taken, "fields local time")

    do istep=1,nstep
      call fields_local_tests%getfield_local(phi, apar, bpar, do_gather, do_update)
    enddo

    if (proc0) then
      call time_message(.false., time_taken, "fields local time")
      write(*, '(" Time for ",I10," fields_local field solves on ",I6," procs: ",F5.1," s")') nstep, nproc, time_taken(1)
      call append_output_file(timing_unit, &
        benchmark_identifier())
      write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
      call close_output_file(timing_unit)
    end if

  end if 

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  if (proc0) call close_module_test('time_fields_local')

  call finish_mp

end program time_fields_local
