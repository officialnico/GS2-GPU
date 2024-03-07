!> A small program designed to enable the use of mutiple
!> input files at once and / or a queue of input files
!> in a single process. Limited features right now (e.g.
!> no load balancing etc.) but can still be useful.
program multigs2
  use mp, only: init_mp, mp_comm, finish_mp, proc0, nproc, time_message, barrier
  use mp, only: broadcast, iproc, mp_undefined, split, rank_comm, free_comm
  use gs2_main, only : run_gs2, gs2_program_state_type
  use constants, only: run_name_size
  use standard_header, only: date_iso8601
  implicit none
  ! Command line options
  integer :: nbatch !< Number of simultaneous runs to progress
  character(len=:), allocatable :: set_file !< File containing list of jobs to run
  logical :: debug
  ! End of command line options
  type(gs2_program_state_type) :: state
  integer :: sub_comm, nfiles, nproc_per_batch, local_iproc, num_my_jobs, batch_id, i
  integer :: original_comm_world
  logical :: no_work_for_this_proc, local_proc0, actual_proc0
  integer, dimension(:), allocatable :: subset_in_charge
  character(len = run_name_size), dimension(:), allocatable :: files, my_files
  real, dimension(2) :: main_timer, job_timer
  main_timer = 0. ; job_timer = 0.
  call parse_command_line()
  call init_mp
  actual_proc0 = proc0
  original_comm_world = mp_comm
  if (proc0) write(*, '("Run started at ",A)') date_iso8601()

  call time_message(.false., main_timer, '')

  if (proc0) call parse_set_file_to_files(set_file, nfiles, files)

  ! Ensure all processors know how many files there are
  call broadcast(nfiles)

  ! Limit nbatch to be no larger than the number of processors
  if (debug .and. nbatch > nproc .and. proc0) then
     write(*,'("Warning: nbatch (",I0,") > nproc (",I0,") -- reducing nbatch")') nbatch, nproc
  end if
  nbatch = min(nbatch, nproc)
  nbatch = min(nbatch, nfiles)

  ! Distribute file list to all processors
  if (.not. proc0) allocate(files(nfiles))
  call broadcast(files)

  ! Determine which processor sub-set will be responsible for each file
  allocate(subset_in_charge(nfiles))
  do i = 1, nfiles
     ! Striped / strided ownership
     subset_in_charge(i) = 1 + mod(i, nbatch)
  end do

  ! Work out how many processors there are for each subset
  nproc_per_batch = nproc / nbatch

  if (proc0) then
     write(*,'("Processing ",I0," files in ",I0," batches with ",I0," processors/batch")') &
          nfiles, nbatch, nproc_per_batch
  end if

  ! Work out which subset this processor belongs to
  batch_id = 1 + iproc / nproc_per_batch

  ! Determine the list of files this processor is involved in and count them
  my_files = pack(files, subset_in_charge == batch_id)
  num_my_jobs = size(my_files)

  ! Check if we're in a batch without any work (processors don't split
  ! perfectly into the requested number of batches)
  if (num_my_jobs < 1) then
     ! Change batch_id to mp_undefined here so that in the split call
     ! they just get a null communicator
     batch_id = mp_undefined
     no_work_for_this_proc = .true.
  else
     no_work_for_this_proc = .false.
  end if

  ! Split comm world into sub-communicators for each processor subset
  call split(batch_id, sub_comm)

  ! Copy information into the state object
  state%mp_comm = sub_comm
  state%run_name_external = .true.
  state%is_external_job = .true.
  state%print_full_timers = .false.
  state%print_times = .false.

  if (no_work_for_this_proc) then
     local_proc0 = .false.
  else
     call rank_comm(sub_comm, local_iproc)
     local_proc0 = local_iproc == 0
  end if

  if (debug .and. proc0 .and. (nproc_per_batch * nbatch /= nproc)) then
     write(*,'("Note : ",I0," processors unused.")') nproc - nproc_per_batch * nbatch
  end if

  do i = 1, num_my_jobs
     state%run_name = my_files(i)
     job_timer = 0.0
     if (debug .and. local_proc0) &
          write(*, '("Batch : ",I0," Starting job ",I0," of ",I0," : ",A)') &
          batch_id, i, num_my_jobs, trim(state%run_name)
     call time_message(.false., job_timer, '')
     call run_gs2(state, quiet = .true.)
     if (debug .and. local_proc0) &
          write(*, '("Batch : ",I0," Done job ",I0," of ",I0," : ",A," in ",0pF9.3," s")') &
          batch_id, i, num_my_jobs, trim(state%run_name), job_timer(1)
  end do
  call time_message(.false., main_timer, '')
  call free_comm(sub_comm)
  if (debug .and. local_proc0) &
       write(*, '("Batch : ",I0," finished in ",0pF9.3," s")') batch_id, main_timer(1)
  call barrier(original_comm_world)
  call finish_mp
  if (actual_proc0) write(*, '("Run finished at ",A)') date_iso8601()
contains
  !> Parse the command line to determine user options
  subroutine parse_command_line()
    integer :: arg_count, arg_n
    character(len=:), allocatable :: argument
    character(len=*), parameter :: nl = new_line('a')
    logical :: skip
    ! Set default options
    nbatch = 1 ; debug = .false.

    arg_count = command_argument_count()
    skip = .false.
    do arg_n = 1, arg_count
       if (skip) then
          skip = .false.
          cycle
       end if

       call get_arg(arg_n, argument)

       if ((argument == "--nbatch") .or. (argument == "-n")) then
          if (arg_n == arg_count) error stop "Missing nbatch value"
          call get_arg(arg_n + 1, argument)
          read(argument, *) nbatch
          skip = .true.
       else if ((argument == "--debug") .or. (argument == "-d")) then
          debug = .true.
       else
          call get_arg(arg_n, set_file)
       end if
    end do

    if (.not. allocated(set_file)) then
       error stop 'No set_file found when parsing command line'
    end if
  end subroutine parse_command_line

  subroutine get_arg(arg_n, arg)
    integer, intent(in) :: arg_n
    character(len=:), allocatable, intent(in out) :: arg
    integer :: arg_length
    call get_command_argument(arg_n, length=arg_length)
    if (allocated(arg)) deallocate(arg)
    allocate(character(len=arg_length)::arg)
    call get_command_argument(arg_n, arg)
  end subroutine get_arg

  subroutine parse_set_file_to_files(set_file, nfiles, files)
    implicit none
    character(len = *), intent(in) :: set_file
    integer, intent(out) :: nfiles
    character(len = run_name_size), allocatable, dimension(:), intent(out) :: files
    integer :: file_unit, ierr, i
    open(newunit = file_unit, file = set_file, action = 'read')
    ! First count how many jobs we have, based on how many lines in file
    nfiles = 0 ; ierr = 0
    do while (ierr == 0)
       read(unit = file_unit, fmt = *, iostat = ierr)
       nfiles = nfiles + 1
    end do
    nfiles = nfiles - 1
    ! Reset to start of file
    rewind(unit = file_unit)

    ! Read in each line as store in files array
    allocate(files(nfiles))
    do i = 1, nfiles
       read(unit = file_unit, fmt = '(a)') files(i)
    end do
    close(file_unit)
  end subroutine parse_set_file_to_files

end program multigs2
