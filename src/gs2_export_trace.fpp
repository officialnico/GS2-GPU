!> Probe hooks and trace file of the GS2 GPU export program `gs2_export`.
!>
!> Implements gpu/spec/09_export_infra.md §5.3 (hook semantics), §6.2 (trace file layout) and
!> §7.9. One netCDF group per recorded probe:
!>
!>     /s<istep:06d>/a<attempt>/p<seq:03d>_<tag>
!>
!> with attributes tag, istep, attempt, seq, pass, in_advance, the time block (0-D variables)
!> and the arrays listed for the tag in gpu/spec/09a_exporter_notes.md.
!>
!> The hooks only read GS2 data (dummy arrays passed by the probe, or module arrays that are
!> not associated with a dummy at the probe point, 09 §11.19).
module gs2_export_trace
  use gs2_export_writer
  use gs2_export_static, only: export_knobs_type, write_time_block, export_state
  use gs2_export_static, only: write_root_dims, export_global_attributes
  implicit none

  private

  public :: tr_open, tr_close, tr_install_hooks, tr_snapshot, tr_nrecords

  integer :: tr_ncid = -1
  integer :: last_istep = -huge(1), attempt = 0, seq = 0
  integer :: first_step = 1, every = 1
  integer :: tr_nrecords = 0
  logical :: all_stages = .true.
  character(len=2000) :: stages = ''

contains

  !> Create the trace file: same global attributes and root dimensions as the export file
  subroutine tr_open(fname, knobs, header)
    use standard_header, only: standard_header_type
    implicit none
    character(len=*), intent(in) :: fname
    type(export_knobs_type), intent(in) :: knobs
    type(standard_header_type), intent(in), optional :: header
    tr_ncid = ex_create(fname, 'GS2 GPU export trace', header)
    call export_global_attributes(tr_ncid)
    call write_root_dims(tr_ncid, knobs)
    call ex_att(tr_ncid, '', 'nstep_trace', knobs%nstep_trace)
    call ex_att(tr_ncid, '', 'trace_stages', trim(knobs%trace_stages))
    call ex_att(tr_ncid, '', 'trace_first_step', knobs%trace_first_step)
    call ex_att(tr_ncid, '', 'trace_every', knobs%trace_every)
    stages = ' '//trim(adjustl(knobs%trace_stages))//' '
    all_stages = trim(adjustl(knobs%trace_stages)) == 'all'
    if (trim(adjustl(knobs%trace_stages)) == 'step') stages = ' advance_begin advance_end step_end '
    first_step = knobs%trace_first_step
    every = max(knobs%trace_every, 1)
    last_istep = -huge(1) ; attempt = 0 ; seq = 0 ; tr_nrecords = 0
  end subroutine tr_open

  !> Close the trace file
  subroutine tr_close
    implicit none
    call ex_att(tr_ncid, '', 'n_records', tr_nrecords)
    call ex_close(tr_ncid)
    tr_ncid = -1
  end subroutine tr_close

  !> Point the probe procedure pointers at the hooks of this module
  subroutine tr_install_hooks
    use gs2_export_probe, only: probe_point_hook, probe_g_hook, probe_le_hook, probe_real_hook
    implicit none
    probe_point_hook => tr_point
    probe_g_hook => tr_g
    probe_le_hook => tr_le
    probe_real_hook => tr_real
  end subroutine tr_install_hooks

  !> State snapshot group (e.g. /final after evolve_equations returned)
  subroutine tr_snapshot(gname, knobs)
    implicit none
    character(len=*), intent(in) :: gname
    type(export_knobs_type), intent(in) :: knobs
    call export_state(tr_ncid, gname, -1, knobs)
  end subroutine tr_snapshot

  !> Tags recorded only between advance_begin and advance_end (filters the istep = 0 calls of
  !> the response-matrix build, 09 §5.3 / Pitfall 11.17)
  logical function inner_tag(tag)
    implicit none
    character(len=*), intent(in) :: tag
    select case (tag)
    case ('advance_begin', 'advance_end', 'after_reset_time_step', 'step_end_pre_diag', 'step_end')
       inner_tag = .false.
    case default
       inner_tag = .true.
    end select
  end function inner_tag

  !> Tag selection by the trace_stages knob
  logical function selected_tag(tag)
    implicit none
    character(len=*), intent(in) :: tag
    if (all_stages) then
       selected_tag = .true.
    else
       selected_tag = index(stages, ' '//trim(tag)//' ') > 0
    end if
  end function selected_tag

  !> Step selection by trace_first_step / trace_every
  logical function selected_step(istep)
    implicit none
    integer, intent(in) :: istep
    selected_step = istep >= first_step
    if (selected_step) selected_step = mod(istep - first_step, every) == 0
  end function selected_step

  !> Is the stage behind this tag active in this run? (stages that do nothing are not recorded)
  logical function stage_active(tag)
    use dist_fn_arrays, only: kx_shift, theta0_shift
    use dist_fn, only: def_parity
    use hyper, only: hyper_on
    use collisions, only: colls, split_collisions
    use fields, only: remove_zonal_flows_switch
    implicit none
    character(len=*), intent(in) :: tag
    select case (tag)
    case ('after_exb_shear')
       stage_active = allocated(kx_shift) .or. allocated(theta0_shift)
    case ('after_hyper')
       stage_active = hyper_on
    case ('after_vspace_derivatives')
       stage_active = colls
    case ('after_parity')
       stage_active = def_parity
    case ('after_remove_zonal_flows')
       stage_active = remove_zonal_flows_switch
    case ('after_collisions_advance')
       stage_active = split_collisions
    case default
       stage_active = .true.
    end select
  end function stage_active

  !> Common bookkeeping of every hook; returns the record group id, or -1 if not recorded
  integer function begin_record(tag, istep) result(gid)
    use gs2_export_probe, only: probe_in_advance, probe_pass
    use mp, only: proc0
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: istep
    character(len=128) :: path
    gid = -1
    if (tr_ncid < 0 .and. proc0) return
    if (inner_tag(tag) .and. .not. probe_in_advance) return
    if (tag == 'advance_begin') then
       if (istep == last_istep) then
          attempt = attempt + 1
       else
          attempt = 1
       end if
       last_istep = istep
       seq = 0
    end if
    if (.not. selected_tag(tag)) return
    if (.not. selected_step(istep)) return
    if (.not. stage_active(tag)) return
    seq = seq + 1
    tr_nrecords = tr_nrecords + 1
    write (path, '("s",I6.6,"/a",I0,"/p",I3.3,"_",A)') max(istep, 0), max(attempt, 1), seq, trim(tag)
    gid = ex_group(tr_ncid, trim(path))
    if (.not. proc0) gid = 0
    call ex_att(gid, '', 'tag', trim(tag))
    call ex_att(gid, '', 'istep', istep)
    call ex_att(gid, '', 'attempt', max(attempt, 1))
    call ex_att(gid, '', 'seq', seq)
    call ex_att(gid, '', 'pass', probe_pass)
    call ex_att(gid, '', 'in_advance', merge(1, 0, probe_in_advance))
    call write_time_block(gid, istep)
  end function begin_record

  !> Hook for probes without data: dump module state chosen by the tag (09 §5.3)
  subroutine tr_point(tag, istep)
    use dist_fn_arrays, only: g, gnew, g_work, gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: jump, ikx_indexed
    use hyper, only: hypervisc_filter
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: istep
    integer :: gid
    gid = begin_record(tag, istep)
    if (gid < 0) return
    select case (tag)
    case ('advance_begin', 'step_end_pre_diag')
       call state_arrays(gid, .false., .false.)
    case ('advance_end')
       call state_arrays(gid, .false., .true.)
    case ('step_end')
       call state_arrays(gid, .true., .false.)
       call ex_sync(tr_ncid)
    case ('after_exb_shear')
       call g7(gid, 'gnew', gnew)
       call fields_new(gid)
       call shifts(gid)
       call i1(gid, 'jump', jump, 'ky')
       call i1(gid, 'ikx_indexed', ikx_indexed, 'kx')
    case ('after_add_explicit_terms')
       call g7(gid, 'gexp_1', gexp_1)
       call g7(gid, 'gexp_2', gexp_2)
       call g7(gid, 'gexp_3', gexp_3)
    case ('source')
       call g7(gid, 'g_work', g_work)
    case ('after_invert_rhs_1', 'after_invert_rhs', 'after_vspace_derivatives', 'after_parity', &
         'after_collisions_advance')
       call g7(gid, 'gnew', gnew)
    case ('after_timeadv')
       call g7(gid, 'gnew', gnew)
       call fields_old(gid)
       call fields_new(gid)
    case ('after_hyper')
       call g7(gid, 'gnew', gnew)
       if (allocated(hypervisc_filter)) then
          call ex_vec(gid, 'hypervisc_filter', hypervisc_filter, 'theta kx ky', fname='hyper::hypervisc_filter')
       end if
    case ('after_getfield', 'after_field_update', 'after_remove_zonal_flows')
       call fields_new(gid)
    case ('after_reset_time_step')
       ! time block only
    case default
       call ex_att(gid, '', 'note', 'unknown tag: time block only')
    end select
  contains
    subroutine g7(gid, name, arr)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name
      complex, dimension(:,:,:), allocatable, intent(in) :: arr
      if (allocated(arr)) then
         call ex_g7c(gid, name, arr, fname='dist_fn_arrays::'//name)
      else
         call ex_absent(gid, name)
      end if
    end subroutine g7
    subroutine f3(gid, name, arr)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name
      complex, dimension(:,:,:), allocatable, intent(in) :: arr
      if (allocated(arr)) then
         call ex_vec(gid, name, arr, 'ri theta kx ky', fname='fields_arrays::'//name, tag='F3')
      else
         call ex_absent(gid, name)
      end if
    end subroutine f3
    subroutine i1(gid, name, arr, dims)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name, dims
      integer, dimension(:), allocatable, intent(in) :: arr
      if (allocated(arr)) then
         call ex_vec(gid, name, arr, dims, fname='dist_fn::'//name)
      else
         call ex_absent(gid, name)
      end if
    end subroutine i1
    subroutine fields_new(gid)
      integer, intent(in) :: gid
      call f3(gid, 'phinew', phinew)
      call f3(gid, 'aparnew', aparnew)
      call f3(gid, 'bparnew', bparnew)
    end subroutine fields_new
    subroutine fields_old(gid)
      integer, intent(in) :: gid
      call f3(gid, 'phi', phi)
      call f3(gid, 'apar', apar)
      call f3(gid, 'bpar', bpar)
    end subroutine fields_old
    subroutine shifts(gid)
      integer, intent(in) :: gid
      if (allocated(kx_shift)) then
         call ex_vec(gid, 'kx_shift', kx_shift, 'ky', fname='dist_fn_arrays::kx_shift')
      else
         call ex_absent(gid, 'kx_shift')
      end if
      if (allocated(theta0_shift)) then
         call ex_vec(gid, 'theta0_shift', theta0_shift, 'ky', fname='dist_fn_arrays::theta0_shift')
      else
         call ex_absent(gid, 'theta0_shift')
      end if
    end subroutine shifts
    subroutine state_arrays(gid, with_g, with_old_fields)
      integer, intent(in) :: gid
      logical, intent(in) :: with_g, with_old_fields
      call g7(gid, 'gnew', gnew)
      if (with_g) call g7(gid, 'g', g)
      call fields_new(gid)
      if (with_old_fields) call fields_old(gid)
      call g7(gid, 'gexp_1', gexp_1)
      call g7(gid, 'gexp_2', gexp_2)
      call g7(gid, 'gexp_3', gexp_3)
      call shifts(gid)
    end subroutine state_arrays
  end subroutine tr_point

  !> Hook for probes passing a g_lo array (dummy argument of the caller)
  subroutine tr_g(tag, istep, arr)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: istep
    complex, dimension(:,:,:), intent(in) :: arr
    integer :: gid
    character(len=8) :: name
    gid = begin_record(tag, istep)
    if (gid < 0) return
    select case (tag)
    case ('coll_h_in', 'coll_h_out')
       name = 'h'
    case ('nl_before_center', 'nl_after_center')
       name = 'g1'
    case default
       name = 'arr'
    end select
    call ex_g7c(gid, trim(name), arr, fname='probe argument ('//trim(tag)//')')
  end subroutine tr_g

  !> Hook for probes passing an le_lo array gle(nxi+1, negrid+1, le_lo)
  subroutine tr_le(tag, istep, arr)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: istep
    complex, dimension(:,:,:), intent(in) :: arr
    integer :: gid
    gid = begin_record(tag, istep)
    if (gid < 0) return
    call ex_le6c(gid, 'gle', arr, 'xi_p1', 'egrid_p1', fname='probe argument gle ('//trim(tag)//')')
  end subroutine tr_le

  !> Hook for probes passing a real vector (nl_cfl: max_vel_components)
  subroutine tr_real(tag, istep, vec)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: istep
    real, dimension(:), intent(in) :: vec
    integer :: gid
    gid = begin_record(tag, istep)
    if (gid < 0) return
    if (size(vec) == 3) then
       call ex_vec(gid, 'max_vel_components', vec, 'three', fname='nonlinear_terms::max_vel_components')
    else
       call ex_att(gid, '', 'note', 'unexpected vector size')
    end if
  end subroutine tr_real
end module gs2_export_trace
