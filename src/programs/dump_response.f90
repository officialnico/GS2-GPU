!> A program to generate response matrix data for a given input file 
!> and dump the resulting arrays to file to be used by subsequent runs.
!>
!> By default this will produce a single set of response matrix files
!> for the time step in the provided input file. The input file can
!> contain a `dump_response_knobs` namelist which specifies
!> `n_time_steps`. If this is greater than one then the program will
!> produce response matrix files for `n_time_steps`
!> timesteps. Specifically, it will produce these files for time steps
!> {delt/delt_adj^N, N=0, n_time_steps-1}. This may be useful for
!> precomputing the response matrices required for a nonlinear run
!> where the timestep is allowed to vary and take values from this sequence.
!>
!> Written by : David Dickinson (ddickinson@users.sourceforge.net)
program dump_response
  use job_manage, only: job_fork
  use mp, only: init_mp, proc0, nproc, broadcast
  use file_utils, only: init_file_utils, run_name, input_unit_exist
  use constants, only: run_name_size
  implicit none
  logical :: list, exist, jnk_exit
  logical, parameter :: quiet=.false.
  character (run_name_size), target :: cbuff
  integer :: istep, n_time_steps, in_file
  namelist /dump_response_knobs/n_time_steps

  call init_mp
  
  ! Report # of processors being used
  if (proc0) then
     if(.not.quiet)then
        if (nproc == 1) then
           write(*,'("Running on ",I0," processor")') nproc
        else
           write(*,'("Running on ",I0," processors")') nproc
        end if
     endif

     ! Call init_file_utils, ie. initialize the inputs and outputs, checking
     ! whether we are doing a run or a list of runs.
     call init_file_utils (list, name="gs")
  end if

  call broadcast (list)
  
  ! If given a list of jobs, fork
  if (list) call job_fork
  
  if (proc0) cbuff = trim(run_name)
  call broadcast (cbuff)
  if (.not. proc0) run_name => cbuff
  
  if(proc0) then
     n_time_steps=1
     in_file = input_unit_exist ("dump_response_knobs", exist)
     if(exist) read(unit=in_file,nml=dump_response_knobs)
  endif
  call broadcast(n_time_steps)

  ! This routine initialises the response matrix and dumps it to file.
  ! This call will use the time step as specified in the input file.
  call init_and_dump(quiet)
  ! Now we deal with the other time steps in the sequence if n_time_steps > 1
  do istep=2,n_time_steps
     !Reduce the time step
     call next_time_step

     !Dump response
     call init_and_dump(quiet)
  enddo

  !Clean up
  call finish_dump_response

contains

  !> Changes to the next time step in the sequence as determined by the usual
  !> rules to reduce the timestep in a nonlinear simulation.
  subroutine next_time_step
    use gs2_time, only: code_dt, save_dt
    use gs2_reinit, only: reduce_time_step
    use fields, only: f_reset => reset_init
    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    implicit none

    !First reduce the time step
    call reduce_time_step

    !Save the time step
    call save_dt(code_dt)

    !Now reset the modules
    call d_reset
    call c_reset
    call f_reset
  end subroutine next_time_step

  !> Initialise the response matrices and write to file
  subroutine init_and_dump(quiet)
    use fields, only: fields_pre_init, fields_init_response, set_dump_and_read_response, dump_response_to_file
    use job_manage, only: time_message
    use gs2_time, only: code_dt
    use mp, only: proc0
    implicit none
    logical, intent(in) :: quiet
    real :: time_measure(2)
    character (len=64) :: str_code_dt

    if(.not.quiet) write(*,'("")')
    time_measure = 0.

    !Prepare to initialise the fields
    if(proc0.and.(.not.quiet)) then 
       write(*,'("Pre-init")',advance='no') 
       call time_message(.false.,time_measure,'Pre-init')
    endif
    call fields_pre_init
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Pre-init')
       write(*,'("  --> Done in : ",ES14.6E3," seconds")') time_measure(1)
    endif
    time_measure=0.

    !Now need to disable dump_response and read_response before
    !we try to initialise the response matrix
    call set_dump_and_read_response(.false.,.false.)

    !Now call the initialisation routine
    if(proc0.and.(.not.quiet)) then
       write(*,'("Init")', advance='no') 
       call time_message(.false.,time_measure,'Init')
    endif
    call fields_init_response
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Init')
       write(*,'("  --> Done in : ",ES14.6E3," seconds")') time_measure(1)
    endif
    time_measure=0.

    !Now we make the dt string
    write(str_code_dt,'(ES14.6E3)') code_dt

    !Now write response to files
    if(proc0.and.(.not.quiet)) then
       write(*,'("Dumping to file for timestep ",A,A,A)',advance='no') "'",trim(adjustl(str_code_dt)),"'"
       call time_message(.false.,time_measure,'Dump')
    endif
    call dump_response_to_file()
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Dump')
       write(*,'("  --> Done in : ",ES14.6E3," seconds")') time_measure(1)
    endif
    time_measure=0.
  end subroutine init_and_dump

  !> Call the finialise routines for any module we have initialised.
  subroutine finish_dump_response
    use antenna, only: finish_antenna
    use collisions, only: finish_collisions
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use file_utils, only: finish_file_utils
    use hyper, only: finish_hyper
    use init_g, only: finish_init_g
    use kt_grids, only: finish_kt_grids
    use le_grids, only: finish_le_grids
    use mp, only: proc0, finish_mp
    use run_parameters, only: finish_run_parameters
    use species, only: finish_species

    implicit none

    call finish_antenna
    call finish_collisions
    call finish_dist_fn
    call finish_fields
    call finish_hyper
    call finish_init_g
    call finish_kt_grids
    call finish_le_grids
    call finish_run_parameters
    call finish_species
    if (proc0) call finish_file_utils
    call finish_mp

  end subroutine finish_dump_response
end program dump_response

