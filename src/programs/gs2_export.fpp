!> Export program for the JAX (GPU) port of GS2: dumps every array the time step reads into a
!> netCDF "export file" and, optionally, steps the Fortran code while dumping intermediate
!> stage results into a "trace file".
!>
!> Usage:
!>
!>     mpirun -np 1 gs2_export input.in [export_knobs.nml]
!>
!> Options come from the namelist &gs2_export_knobs, read from the file named by the second
!> command argument if present, otherwise from the GS2 input file (GS2 ignores the group).
!> The control flow mirrors gs2_main::run_gs2 (gpu/spec/09_export_infra.md §1.4), so the
!> initialisation path and the traced trajectory are those of bin/gs2.
!>
!> See gpu/spec/09_export_infra.md and gpu/spec/09a_exporter_notes.md.
program gs2_export
  use mp, only: init_mp, mp_comm, finish_mp, proc0, nproc, broadcast, mp_abort
  use gs2_main, only: gs2_program_state_type, initialize_gs2, initialize_equations
  use gs2_main, only: initialize_diagnostics, evolve_equations, finalize_diagnostics
  use gs2_main, only: finalize_equations, finalize_gs2, write_used_inputs_file, parse_command_line
  use standard_header, only: standard_header_type
  use file_utils, only: run_name
  use gs2_export_probe, only: probe_active
  use gs2_export_writer, only: ex_create, ex_close, ex_att, ex_group, ex_compression
  use gs2_export_static, only: export_knobs_type, write_root_dims, export_static, export_state
  use gs2_export_static, only: export_dt_dependent, export_diagnostics_group, export_gates
  use gs2_export_static, only: export_global_attributes, export_text_attributes
  use gs2_export_trace, only: tr_open, tr_close, tr_install_hooks, tr_snapshot, tr_nrecords
#ifdef NEW_DIAG
  use gs2_diagnostics_new, only: gnostics
#endif
  implicit none
  type(gs2_program_state_type) :: state
  type(standard_header_type) :: header
  type(export_knobs_type) :: knobs
  integer :: ncid
  logical :: diag_ok
  character(len=:), allocatable :: gate_status
  character(len=:), allocatable :: export_file, trace_file

  call parse_command_line()
  call init_mp
  state%mp_comm = mp_comm

  call initialize_gs2(state, header_out=header)
  call read_export_knobs(knobs)
  call check_hard_gates_before_init(state, knobs)
  ex_compression = knobs%compression

  if (len_trim(knobs%export_file) == 0) then
     export_file = trim(run_name)//'.gs2export.nc'
  else
     export_file = trim(knobs%export_file)
  end if
  if (len_trim(knobs%trace_file) == 0) then
     trace_file = trim(run_name)//'.gs2trace.nc'
  else
     trace_file = trim(knobs%trace_file)
  end if

  if (knobs%trace_response_build .and. proc0) write (*, '(A)') 'gs2_export: WARNING: '// &
       'trace_response_build (probe P-R1) is not available in this build; see gpu/spec/09a_exporter_notes.md'

  ! init(full) + init_tstart + init_gs2_reinit, exactly as bin/gs2
  call initialize_equations(state)
  call check_hard_gates_after_init(state, knobs)

  ! Static export: grids, geometry, coefficients, configuration
  ncid = ex_create(export_file, 'GS2 GPU export', header)
  call export_global_attributes(ncid)
  call write_root_dims(ncid, knobs)
  call export_static(ncid, knobs)
  call export_state(ncid, 'state_t0_prediag', -1, knobs)

  ! Diagnostics at istep 0 (g_adjust round trips on gnew and g, 09 §11.2, 10 G-05)
  if (knobs%run_diagnostics) then
     call initialize_diagnostics(state, header=header)
     diag_ok = .true.
  else
     state%skip_diagnostics = .true.
#ifdef NEW_DIAG
     gnostics%nsave = -1
#endif
     diag_ok = .false.
  end if
  call write_used_inputs_file(header=header)

  ! THE initial state of JAX step 1
  call export_state(ncid, 'state_t0', -1, knobs)
  if (diag_ok) call export_diagnostics_group(ncid)
  gate_status = export_gates(ncid, diag_ok)
  call ex_att(ncid, '', 'gate_status', gate_status)
  call ex_att(ncid, '', 'diagnostics_run_at_t0', merge(1, 0, diag_ok))
  call export_text_attributes(ncid)
  if (proc0) write (*, '(A)') 'gs2_export: gate_status = '//gate_status

  if (knobs%write_restart_t0) call write_restart_t0(header)

  if (knobs%nstep_trace == 0) then
     ! Close <run>.out.nc now (static part + istep-0 record, 10 G-05), before a dt ladder
     ! re-initialises the equations
     if (state%init%diagnostics_initialized) call finalize_diagnostics(state)
     if (knobs%dt_ladder_n > 0) call export_dt_ladder(ncid)
  end if
  call ex_close(ncid)
  if (proc0) write (*, '(A)') 'gs2_export: wrote '//export_file

  ! Fortran-stepped trace
  if (knobs%nstep_trace > 0) then
     call tr_open(trace_file, knobs, header)
     call tr_install_hooks
     probe_active = .true.
     state%print_times = .false.
     call evolve_equations(state, min(knobs%nstep_trace, state%nstep), header=header)
     probe_active = .false.
     call tr_snapshot('final', knobs)
     call tr_close
     if (proc0) write (*, '(A,I0,A)') 'gs2_export: wrote '//trace_file//' (', tr_nrecords, ' records)'
     if (state%init%diagnostics_initialized) call finalize_diagnostics(state)
  end if

  call finalize_equations(state)
  call finalize_gs2(state)
  call finish_mp

contains

  !> Read &gs2_export_knobs (09 §8.1) on proc0 and broadcast
  subroutine read_export_knobs(knobs)
    use file_utils, only: input_unit_exist
    implicit none
    type(export_knobs_type), intent(out) :: knobs
    character(len=len(knobs%export_file)) :: export_file, trace_file
    character(len=len(knobs%trace_stages)) :: trace_stages
    logical :: write_response, write_collisions, write_g_adj, write_raw_glo, write_restart_t0
    logical :: run_diagnostics, trace_response_build, allow_multi_rank, strict, exist
    integer :: nstep_trace, trace_first_step, trace_every, dt_ladder_n, compression, unit, ios
    character(len=2000) :: arg
    namelist /gs2_export_knobs/ export_file, write_response, write_collisions, write_g_adj, &
         write_raw_glo, write_restart_t0, run_diagnostics, nstep_trace, trace_file, trace_stages, &
         trace_first_step, trace_every, trace_response_build, dt_ladder_n, allow_multi_rank, &
         strict, compression

    export_file = knobs%export_file ; trace_file = knobs%trace_file ; trace_stages = knobs%trace_stages
    write_response = knobs%write_response ; write_collisions = knobs%write_collisions
    write_g_adj = knobs%write_g_adj ; write_raw_glo = knobs%write_raw_glo
    write_restart_t0 = knobs%write_restart_t0 ; run_diagnostics = knobs%run_diagnostics
    nstep_trace = knobs%nstep_trace ; trace_first_step = knobs%trace_first_step
    trace_every = knobs%trace_every ; trace_response_build = knobs%trace_response_build
    dt_ladder_n = knobs%dt_ladder_n ; allow_multi_rank = knobs%allow_multi_rank
    strict = knobs%strict ; compression = knobs%compression

    if (proc0) then
       if (command_argument_count() >= 2) then
          call get_command_argument(2, arg)
          open (newunit=unit, file=trim(arg), status='old', action='read', iostat=ios)
          if (ios /= 0) call mp_abort('gs2_export: cannot open knobs file '//trim(arg), .true.)
          read (unit, nml=gs2_export_knobs, iostat=ios)
          if (ios /= 0) call mp_abort('gs2_export: error reading &gs2_export_knobs from '//trim(arg), .true.)
          close (unit)
       else
          unit = input_unit_exist('gs2_export_knobs', exist)
          if (exist) read (unit=unit, nml=gs2_export_knobs)
       end if
    end if

    call broadcast(export_file) ; call broadcast(trace_file) ; call broadcast(trace_stages)
    call broadcast(write_response) ; call broadcast(write_collisions) ; call broadcast(write_g_adj)
    call broadcast(write_raw_glo) ; call broadcast(write_restart_t0) ; call broadcast(run_diagnostics)
    call broadcast(nstep_trace) ; call broadcast(trace_first_step) ; call broadcast(trace_every)
    call broadcast(trace_response_build) ; call broadcast(dt_ladder_n)
    call broadcast(allow_multi_rank) ; call broadcast(strict) ; call broadcast(compression)

    knobs%export_file = export_file ; knobs%trace_file = trace_file ; knobs%trace_stages = trace_stages
    knobs%write_response = write_response ; knobs%write_collisions = write_collisions
    knobs%write_g_adj = write_g_adj ; knobs%write_raw_glo = write_raw_glo
    knobs%write_restart_t0 = write_restart_t0 ; knobs%run_diagnostics = run_diagnostics
    knobs%nstep_trace = nstep_trace ; knobs%trace_first_step = trace_first_step
    knobs%trace_every = trace_every ; knobs%trace_response_build = trace_response_build
    knobs%dt_ladder_n = dt_ladder_n ; knobs%allow_multi_rank = allow_multi_rank
    knobs%strict = strict ; knobs%compression = compression
  end subroutine read_export_knobs

  !> Hard gates that can be checked before the equations are initialised (09 §8.3)
  subroutine check_hard_gates_before_init(state, knobs)
    use, intrinsic :: iso_fortran_env, only: real64
    use build_config, only: gs2_has_lowflow, gs2_is_single_precision
#ifdef OPENMP
    use omp_lib, only: omp_get_max_threads
#endif
    implicit none
    type(gs2_program_state_type), intent(in) :: state
    type(export_knobs_type), intent(in) :: knobs
    if (gs2_is_single_precision .or. kind(1.0) /= real64) &
         call mp_abort('gs2_export: single precision build (hard gate)', .true.)
#ifdef SHMEM
    call mp_abort('gs2_export: SHMEM build (hard gate)', .true.)
#endif
    if (gs2_has_lowflow) call mp_abort('gs2_export: LOWFLOW build (hard gate)', .true.)
#ifdef OPENMP
    if (omp_get_max_threads() > 1) then
       if (knobs%strict) call mp_abort('gs2_export: OpenMP with more than one thread (hard gate)', .true.)
       if (proc0) write (*, '(A)') 'gs2_export: WARNING: more than one OpenMP thread'
    end if
#endif
    if (nproc > 1) then
       ! Multi-rank export (gathers of 09 §2.5) is not implemented; reference data must come from
       ! a single rank anyway (09 §11.1).
       call mp_abort('gs2_export: nproc > 1 is not supported (run with mpirun -np 1)', .true.)
    end if
    if (state%list .or. state%nensembles > 1) &
         call mp_abort('gs2_export: list mode / ensembles are not supported (hard gate)', .true.)
    if (knobs%dt_ladder_n > 0 .and. knobs%nstep_trace > 0) &
         call mp_abort('gs2_export: dt_ladder_n > 0 and nstep_trace > 0 are mutually exclusive', .true.)
  end subroutine check_hard_gates_before_init

  !> Hard gates that need the initialised equations (09 §8.3)
  subroutine check_hard_gates_after_init(state, knobs)
    use run_parameters, only: do_eigsolve
    implicit none
    type(gs2_program_state_type), intent(in) :: state
    type(export_knobs_type), intent(inout) :: knobs
    if (do_eigsolve) call mp_abort('gs2_export: do_eigsolve = .true. (hard gate)', .true.)
    if (knobs%nstep_trace > state%nstep) then
       if (proc0) write (*, '(A,I0)') 'gs2_export: WARNING: nstep_trace clamped to nstep = ', state%nstep
       knobs%nstep_trace = state%nstep
    end if
  end subroutine check_hard_gates_after_init

  !> GS2 restart file of the t = 0 state (09 §6.4) for cross-checks of the canonical map
  subroutine write_restart_t0(header)
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: gnew
    use collisions, only: vnmult
    use run_parameters, only: has_phi, has_apar, has_bpar
    use gs2_time, only: user_time, code_dt, code_dt_prev1, code_dt_prev2, code_dt_max
    implicit none
    type(standard_header_type), intent(in) :: header
    call gs2_save_for_restart(gnew, user_time, vnmult, has_phi, has_apar, has_bpar, &
         code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, fileopt='.gs2export_t0', &
         save_glo_info_and_grids=.true., save_velocities=.true., header=header)
  end subroutine write_restart_t0

  !> dt ladder (09 §6.3, §7.8): mirrors the rebuild of reset_time_step (gs2_reinit.f90:104-225)
  !> without state preservation; the run does not step afterwards. Rungs code_dt/delt_adj**k.
  subroutine export_dt_ladder(ncid)
    use gs2_init, only: init, init_level_list
    use gs2_reinit, only: reduce_time_step, delt_adj
    use gs2_time, only: code_dt, code_dt_min, save_dt
    use dist_fn_arrays, only: gnew
    use array_utils, only: zero_array
    implicit none
    integer, intent(in) :: ncid
    integer :: k, gid, nrung
    character(len=32) :: path
    nrung = 0
    do k = 1, knobs%dt_ladder_n
       if (code_dt/delt_adj <= code_dt_min) then
          if (proc0) write (*, '(A,I0)') 'gs2_export: dt ladder stopped at code_dt_min before rung ', k
          exit
       end if
       ! As reset_time_step: zero gnew (gs2_reinit.f90:163), down to override_timestep (:166),
       ! reduce dt, save_dt, back up to full (:216)
       call zero_array(gnew)
       call init(state%init, init_level_list%override_timestep)
       call reduce_time_step
       call save_dt(code_dt)
       call init(state%init, init_level_list%full)
       write (path, '("dt_ladder/k",I0)') k
       call export_dt_dependent(ncid, trim(path), knobs)
       nrung = k
       if (proc0) write (*, '(A,I0,A,ES24.16)') 'gs2_export: dt ladder rung ', k, ': code_dt = ', code_dt
    end do
    gid = ex_group(ncid, 'dt_ladder')
    call ex_att(gid, '', 'n_rungs', nrung)
    call ex_att(gid, '', 'delt_adj', delt_adj)
    call ex_att(gid, '', 'note', 'rung k: code_dt = dt_init/delt_adj**k (downward rungs only, 10 G-19); '// &
         'arrays as rebuilt by the levels 30-32 of the init ladder; D_res and the simple3D filter '// &
         'stay at their first-init values (level 19 is not re-entered)')
  end subroutine export_dt_ladder
end program gs2_export
