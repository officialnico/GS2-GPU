!> Static part of the GS2 GPU export (`gs2_export`): grids, geometry, coefficients,
!> configuration, the t = 0 state and the dt-dependent coefficient set of a dt-ladder rung.
!>
!> Contract: gpu/spec/09_export_infra.md §3 (inventory), §6.1 (file layout), §7.3-7.7, with
!> the corrections of gpu/spec/10_critique.md §A. Implementation notes and deviations:
!> gpu/spec/09a_exporter_notes.md.
!>
!> Everything here only *reads* GS2 module data (09 §11.5): no physics routine is called,
!> apart from pure/elemental helpers and the public response-dump routine of fields_local.
module gs2_export_static
  use, intrinsic :: iso_fortran_env, only: real64
  use constants, only: run_name_size
  use gs2_export_writer
  implicit none

  private

  public :: export_knobs_type
  public :: write_root_dims, export_static, export_state, export_dt_dependent
  public :: write_time_block, export_diagnostics_group, export_gates, export_global_attributes
  public :: export_text_attributes

  !> Options of the &gs2_export_knobs namelist (09 §8.1)
  type :: export_knobs_type
     character(len=run_name_size) :: export_file = ''
     logical :: write_response = .true.
     logical :: write_collisions = .true.
     logical :: write_g_adj = .false.
     logical :: write_raw_glo = .false.
     logical :: write_restart_t0 = .false.
     logical :: run_diagnostics = .true.
     integer :: nstep_trace = 0
     character(len=run_name_size) :: trace_file = ''
     character(len=2000) :: trace_stages = 'all'
     integer :: trace_first_step = 1
     integer :: trace_every = 1
     logical :: trace_response_build = .false.
     integer :: dt_ladder_n = 0
     logical :: allow_multi_rank = .false.
     logical :: strict = .true.
     integer :: compression = 0
  end type export_knobs_type

  !> Version of the export format (bump the major number on incompatible changes)
  character(len=*), parameter :: export_format_version = '1.0'

contains

  !============================================================== helpers for allocatables

  subroutine wr_r1(gid, name, arr, dims, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_r1

  subroutine wr_r2(gid, name, arr, dims, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_r2

  subroutine wr_r3(gid, name, arr, dims, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_r3

  subroutine wr_r4(gid, name, arr, dims, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_r4

  subroutine wr_i1(gid, name, arr, dims, fname, src, owner, ibase)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer, intent(in), optional :: ibase
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, ibase=ibase)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_i1

  subroutine wr_i2(gid, name, arr, dims, fname, src, owner, ibase)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer, intent(in), optional :: ibase
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, ibase=ibase)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_i2

  subroutine wr_i3(gid, name, arr, dims, fname, src, owner, ibase)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer, intent(in), optional :: ibase
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner, ibase=ibase)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_i3

  subroutine wr_l1(gid, name, arr, dims, fname, src, owner)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    logical, dimension(:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_l1

  subroutine wr_l2(gid, name, arr, dims, fname, src, owner)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    logical, dimension(:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    if (allocated(arr)) then
       call ex_vec(gid, name, arr, dims, fname=fname, src=src, owner=owner)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_l2

  subroutine wr_g7c(gid, name, arr, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    complex, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_g7c(gid, name, arr, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_g7c

  subroutine wr_g7r(gid, name, arr, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    real, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_g7r(gid, name, arr, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_g7r

  subroutine wr_g6r(gid, name, arr, fname, src, owner)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    real, dimension(:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    if (allocated(arr)) then
       call ex_g6r(gid, name, arr, fname=fname, src=src, owner=owner)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_g6r

  subroutine wr_g7m1c(gid, name, arr, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    complex, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_g7m1c(gid, name, arr, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_g7m1c

  subroutine wr_le6r(gid, name, arr, xdim, edim, fname, src, owner, dtdep)
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, xdim, edim
    real, dimension(:,:,:), allocatable, intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    if (allocated(arr)) then
       call ex_le6r(gid, name, arr, xdim, edim, fname=fname, src=src, owner=owner, dtdep=dtdep)
    else
       call ex_absent(gid, name)
    end if
  end subroutine wr_le6r

  !> Integer as text
  function itoa(i) result(s)
    implicit none
    integer, intent(in) :: i
    character(len=16) :: s
    write (s, '(I0)') i
  end function itoa

  !> Are collisions active (09 §6.1.3: /collisions only if collision_model_switch /= none and colls)
  logical function collisions_active()
    use collisions, only: colls, collision_model_switch, collision_model_none
    implicit none
    collisions_active = colls .and. (collision_model_switch /= collision_model_none)
  end function collisions_active

  !============================================================== root dimensions

  !> Root dimensions and coordinate variables (09 §6.1.1)
  subroutine write_root_dims(ncid, knobs)
    use theta_grid, only: ntgrid, theta, nbset, bset
    use kt_grids, only: naky, ntheta0, akx, aky
    use le_grids, only: nlambda, negrid, nxi, al
    use species, only: nspec
    use collisions, only: nxi_lim
    use gs2_layouts, only: g_lo
    use dist_fn, only: i_class, M_class, N_class, g_adj
    use fields, only: fieldopt_switch, fieldopt_implicit
    use fields_implicit, only: nidx
    implicit none
    integer, intent(in) :: ncid
    type(export_knobs_type), intent(in) :: knobs
    integer :: i, c
    call ex_add_dim(ncid, 'ri', 2, ivalues=[1, 2])
    call ex_add_dim(ncid, 'theta', 2*ntgrid+1, rvalues=theta)
    call ex_add_dim(ncid, 'theta_m1', 2*ntgrid, ivalues=[(i, i = 1, 2*ntgrid)])
    call ex_add_dim(ncid, 'sign', 2, ivalues=[1, 2])
    call ex_add_dim(ncid, 'kx', ntheta0, rvalues=akx)
    call ex_add_dim(ncid, 'ky', naky, rvalues=aky)
    call ex_add_dim(ncid, 'lambda', nlambda, rvalues=al)
    call ex_add_dim(ncid, 'egrid', negrid, ivalues=[(i, i = 1, negrid)])
    call ex_add_dim(ncid, 'egrid_p1', negrid+1, ivalues=[(i, i = 1, negrid+1)])
    call ex_add_dim(ncid, 'egrid_m1', max(negrid-1, 1), ivalues=[(i, i = 1, max(negrid-1, 1))])
    call ex_add_dim(ncid, 'species', nspec, ivalues=[(i, i = 1, nspec)])
    call ex_add_dim(ncid, 'xi2', 2*nlambda, ivalues=[(i, i = 1, 2*nlambda)])
    call ex_add_dim(ncid, 'xi_p1', nxi+1, ivalues=[(i, i = 1, nxi+1)])
    if (nxi_lim > 0) call ex_add_dim(ncid, 'xi_lim', nxi_lim, ivalues=[(i, i = 1, nxi_lim)])
    if (nbset > 0 .and. allocated(bset)) call ex_add_dim(ncid, 'bset', nbset, rvalues=bset)
    call ex_add_dim(ncid, 'two', 2)
    call ex_add_dim(ncid, 'three', 3)
    if (knobs%write_raw_glo) call ex_add_dim(ncid, 'glo', g_lo%ulim_world+1, &
         ivalues=[(i, i = 0, g_lo%ulim_world)])
    if (knobs%write_g_adj .and. allocated(g_adj)) call ex_add_dim(ncid, 'nlinks_max', size(g_adj, 1))
    if (allocated(M_class) .and. allocated(N_class) .and. i_class > 0) then
       call ex_add_dim(ncid, 'n_class', i_class, ivalues=[(i, i = 1, i_class)])
       call ex_add_dim(ncid, 'sc_Mmax', maxval(M_class))
       call ex_add_dim(ncid, 'sc_Nmax', maxval(N_class))
       if (fieldopt_switch == fieldopt_implicit .and. nidx > 0) then
          do c = 1, i_class
             call ex_add_dim(ncid, 'sc_R_c'//trim(itoa(c)), nidx*N_class(c))
             call ex_add_dim(ncid, 'sc_M_c'//trim(itoa(c)), M_class(c))
             call ex_add_dim(ncid, 'sc_RM_c'//trim(itoa(c)), nidx*N_class(c)*M_class(c))
          end do
       end if
    end if
  end subroutine write_root_dims

  !============================================================== global attributes

  !> Global attributes of the export and trace files (09 §6.1.1); gate_status, used_inputs and
  !> input_file_text are added later
  subroutine export_global_attributes(ncid)
    use git_version_mod, only: get_git_hash, get_git_state
    use build_config, only: formatted_build_config
    use file_utils, only: run_name
    use mp, only: nproc
    use gs2_layouts, only: g_lo, le_lo, layout
    implicit none
    integer, intent(in) :: ncid
    call ex_att(ncid, '', 'export_format_version', export_format_version)
    call ex_att(ncid, '', 'git_hash', get_git_hash())
    call ex_att(ncid, '', 'git_state', get_git_state())
    call ex_att(ncid, '', 'build_config', formatted_build_config())
    call ex_att(ncid, '', 'run_name', trim(run_name))
    call ex_att(ncid, '', 'input_file', trim(run_name)//'.in')
    call ex_att(ncid, '', 'nproc', nproc)
    call ex_att(ncid, '', 'nproc_g_lo', g_lo%nproc)
    call ex_att(ncid, '', 'nproc_le_lo', le_lo%nproc)
    call ex_att(ncid, '', 'layout', trim(layout))
    call ex_att(ncid, '', 'index_convention', 'Fortran dimension order; NumPy shows reversed dims. '// &
         'Canonical g: (theta,sign,kx,ky,lambda,egrid,species). Complex: leading ri dim (1=re,2=im). '// &
         'Logical: int8. Integers keep Fortran values (1-based indices, -1 sentinels).')
    call ex_att(ncid, '', 'real_kind_bytes', 8)
    call ex_att(ncid, '', 'group_renames', 'species:species_')
  end subroutine export_global_attributes

  !> Read a whole text file (proc0)
  subroutine read_text_file(fname, text, ok)
    implicit none
    character(len=*), intent(in) :: fname
    character(len=:), allocatable, intent(out) :: text
    logical, intent(out) :: ok
    integer :: unit, ios, n
    logical :: exists
    ok = .false.
    inquire (file=trim(fname), exist=exists, size=n)
    if (.not. exists .or. n < 0) then
       text = ''
       return
    end if
    open (newunit=unit, file=trim(fname), access='stream', form='unformatted', &
         action='read', status='old', iostat=ios)
    if (ios /= 0) then
       text = ''
       return
    end if
    allocate (character(len=n) :: text)
    if (n > 0) read (unit, iostat=ios) text
    close (unit)
    ok = (ios == 0)
  end subroutine read_text_file

  !> used_inputs and input_file_text global attributes (09 §6.1.1)
  subroutine export_text_attributes(ncid)
    use file_utils, only: run_name
    use mp, only: proc0
    implicit none
    integer, intent(in) :: ncid
    character(len=:), allocatable :: text
    logical :: ok
    if (.not. proc0) return
    call read_text_file(trim(run_name)//'.used_inputs.in', text, ok)
    if (ok) then
       call ex_att(ncid, '', 'used_inputs', text)
    else
       call ex_att(ncid, '', 'used_inputs', 'unavailable')
    end if
    call read_text_file(trim(run_name)//'.in', text, ok)
    if (ok) then
       call ex_att(ncid, '', 'input_file_text', text)
    else
       call ex_att(ncid, '', 'input_file_text', 'unavailable')
    end if
  end subroutine export_text_attributes

  !============================================================== static export

  !> Write every static group (09 §7.3); the export file must be open with root dims defined
  subroutine export_static(ncid, knobs)
    implicit none
    integer, intent(in) :: ncid
    type(export_knobs_type), intent(in) :: knobs
    call export_theta_grid(ncid)
    call export_kt_grids(ncid)
    call export_le_grids(ncid)
    call export_species(ncid)
    call export_run_parameters(ncid)
    call export_gs2_time(ncid, '/gs2_time')
    call export_gs2_reinit(ncid)
    call export_init_g(ncid)
    call export_fields(ncid)
    call export_gs2_layouts(ncid)
    call export_dist_fn_arrays(ncid, knobs)
    call export_dist_fn(ncid, knobs)
    call export_field_matrices(ncid, '', knobs)
    call export_collisions(ncid, '', knobs, .true.)
    call export_hyper(ncid, '', .true.)
    call export_nonlinear_terms(ncid, '', .true.)
  end subroutine export_static

  !> /theta_grid (09 §3.1, 10 §A.2)
  subroutine export_theta_grid(ncid)
    use theta_grid, only: theta, theta2, delthet, delthet2, bset, bmag, gradpar, itor_over_B, IoB
    use theta_grid, only: gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0
    use theta_grid, only: gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime, Bpol
    use theta_grid, only: bmin, bmax, eps_trapped, shat, drhodpsi, kxfac, qval, surfarea, dvdrhon, rhoc
    use theta_grid, only: ntheta, ntgrid, nperiod, nbset, shape, gb_to_cv, is_effectively_zero_shear
    use theta_grid, only: field_line_average_weight, cvdriftknob, gbdriftknob, eqopt_switch
    use integration, only: trapezoidal_integration
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    character(len=*), parameter :: s = 'theta_grid.f90 finish_init/geometry'
    gid = ex_group(ncid, 'theta_grid')
    call wr_r1(gid, 'theta', theta, 'theta', 'theta_grid::theta', s, '01')
    call wr_r1(gid, 'theta2', theta2, 'theta', 'theta_grid::theta2', s)
    call wr_r1(gid, 'delthet', delthet, 'theta', 'theta_grid::delthet', 'theta_grid.f90:510-511', '02')
    call wr_r1(gid, 'delthet2', delthet2, 'theta', 'theta_grid::delthet2', s)
    if (nbset > 0) then
       call wr_r1(gid, 'bset', bset, 'bset', 'theta_grid::bset', s)
    else
       call ex_absent(gid, 'bset')
    end if
    call wr_r1(gid, 'bmag', bmag, 'theta', 'theta_grid::bmag', s)
    call wr_r1(gid, 'gradpar', gradpar, 'theta', 'theta_grid::gradpar', s)
    call wr_r1(gid, 'itor_over_B', itor_over_B, 'theta', 'theta_grid::itor_over_B', &
         'geometry; btor_slab override dist_fn.fpp read_parameters')
    call wr_r1(gid, 'IoB', IoB, 'theta', 'theta_grid::IoB', s)
    call wr_r1(gid, 'gbdrift', gbdrift, 'theta', 'theta_grid::gbdrift', 'theta_grid.f90:623-626 (knobs applied)')
    call wr_r1(gid, 'gbdrift0', gbdrift0, 'theta', 'theta_grid::gbdrift0', 'theta_grid.f90:623-626 (knobs applied)')
    call wr_r1(gid, 'cvdrift', cvdrift, 'theta', 'theta_grid::cvdrift', 'theta_grid.f90:623-626 (knobs applied)')
    call wr_r1(gid, 'cvdrift0', cvdrift0, 'theta', 'theta_grid::cvdrift0', 'theta_grid.f90:623-626 (knobs applied)')
    call wr_r1(gid, 'cdrift', cdrift, 'theta', 'theta_grid::cdrift', s)
    call wr_r1(gid, 'cdrift0', cdrift0, 'theta', 'theta_grid::cdrift0', s)
    call wr_r1(gid, 'gds2', gds2, 'theta', 'theta_grid::gds2', s)
    call wr_r1(gid, 'gds21', gds21, 'theta', 'theta_grid::gds21', s)
    call wr_r1(gid, 'gds22', gds22, 'theta', 'theta_grid::gds22', s)
    call wr_r1(gid, 'gds23', gds23, 'theta', 'theta_grid::gds23', s)
    call wr_r1(gid, 'gds24', gds24, 'theta', 'theta_grid::gds24', s)
    call wr_r1(gid, 'gds24_noq', gds24_noq, 'theta', 'theta_grid::gds24_noq', s)
    call wr_r1(gid, 'grho', grho, 'theta', 'theta_grid::grho', s)
    call wr_r1(gid, 'jacob', jacob, 'theta', 'theta_grid::jacob', 'theta_grid.f90:516')
    call wr_r1(gid, 'Rplot', Rplot, 'theta', 'theta_grid::Rplot', s)
    call wr_r1(gid, 'Zplot', Zplot, 'theta', 'theta_grid::Zplot', s)
    call wr_r1(gid, 'aplot', aplot, 'theta', 'theta_grid::aplot', s)
    call wr_r1(gid, 'Rprime', Rprime, 'theta', 'theta_grid::Rprime', s)
    call wr_r1(gid, 'Zprime', Zprime, 'theta', 'theta_grid::Zprime', s)
    call wr_r1(gid, 'aprime', aprime, 'theta', 'theta_grid::aprime', s)
    call wr_r1(gid, 'Bpol', Bpol, 'theta', 'theta_grid::Bpol', s)
    call ex_scalar(gid, 'bmin', bmin)
    call ex_scalar(gid, 'bmax', bmax)
    call ex_scalar(gid, 'eps_trapped', eps_trapped)
    call ex_scalar(gid, 'shat', shat)
    call ex_scalar(gid, 'drhodpsi', drhodpsi)
    call ex_scalar(gid, 'kxfac', kxfac)
    call ex_scalar(gid, 'qval', qval)
    call ex_scalar(gid, 'surfarea', surfarea)
    call ex_scalar(gid, 'dvdrhon', dvdrhon)
    call ex_scalar(gid, 'rhoc', rhoc)
    call ex_scalar(gid, 'field_line_average_weight', field_line_average_weight, &
         fname='theta_grid::field_line_average_weight', src='theta_grid.f90:520')
    ! Exporter recompute with the same expression; must be bitwise equal (09 §7.7)
    if (allocated(theta) .and. allocated(jacob)) &
         call ex_scalar(gid, 'field_line_average_weight_check', 1.0/trapezoidal_integration(theta, jacob), &
         src='gs2_export: 1.0/trapezoidal_integration(theta, jacob)')
    call ex_scalar(gid, 'cvdriftknob', cvdriftknob)
    call ex_scalar(gid, 'gbdriftknob', gbdriftknob)
    call ex_scalar(gid, 'ntheta', ntheta)
    call ex_scalar(gid, 'ntgrid', ntgrid)
    call ex_scalar(gid, 'nperiod', nperiod)
    call ex_scalar(gid, 'nbset', nbset)
    call ex_scalar(gid, 'eqopt_switch', eqopt_switch)
    call ex_scalar(gid, 'gb_to_cv', gb_to_cv)
    call ex_scalar(gid, 'is_effectively_zero_shear', is_effectively_zero_shear())
    call ex_att(gid, '', 'shape', trim(shape))
  end subroutine export_theta_grid

  !> /kt_grids (09 §3.2, 10 §A.3); kwork_filter is the value after init_fields
  subroutine export_kt_grids(ncid)
    use kt_grids, only: aky, akx, ikx, theta0, kperp2, inv_kperp2, kwork_filter
    use kt_grids, only: naky, ntheta0, nx, ny, jtwist_out, gridopt_switch, reality, box, grid_option
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    gid = ex_group(ncid, 'kt_grids')
    call wr_r1(gid, 'aky', aky, 'ky', 'kt_grids::aky', 'kt_grids_box.f90:378-380 / range')
    call wr_r1(gid, 'akx', akx, 'kx', 'kt_grids::akx', 'kt_grids_box.f90:382-389')
    call wr_i1(gid, 'ikx', ikx, 'kx', 'kt_grids::ikx', 'kt_grids_box.f90:382-389')
    call wr_r2(gid, 'theta0', theta0, 'kx ky', 'kt_grids::theta0', 'kt_grids_box.f90:391-401')
    call wr_r3(gid, 'kperp2', kperp2, 'theta kx ky', 'kt_grids::kperp2', 'kt_grids.f90:313-349')
    call wr_r3(gid, 'inv_kperp2', inv_kperp2, 'theta kx ky', 'kt_grids::inv_kperp2', 'kt_grids.f90:313-349')
    call wr_l2(gid, 'kwork_filter', kwork_filter, 'kx ky', 'kt_grids::kwork_filter', &
         'kt_grids.f90:188-189; fields.f90:376 (box: (1,1)=.true.)')
    call ex_scalar(gid, 'naky', naky)
    call ex_scalar(gid, 'ntheta0', ntheta0)
    call ex_scalar(gid, 'nx', nx)
    call ex_scalar(gid, 'ny', ny)
    call ex_scalar(gid, 'jtwist_out', jtwist_out)
    call ex_scalar(gid, 'gridopt_switch', gridopt_switch)
    call ex_scalar(gid, 'reality', reality)
    call ex_scalar(gid, 'box', box)
    call ex_att(gid, '', 'grid_option', trim(grid_option))
  end subroutine export_kt_grids

  !> /le_grids (09 §3.3, 10 §A.4) incl. derived masks (09 §7.7)
  subroutine export_le_grids(ncid)
    use le_grids, only: energy, w, speed, energy_maxw, w_maxw, speed_maxw, al, wl, jend
    use le_grids, only: forbid, is_ttp, is_bounce_point, can_be_ttp, xi, wxi, ixi_to_il, ixi_to_isgn, sgn
    use le_grids, only: negrid, nlambda, ng2, lmax, nxi, nesub, wdim, vcut
    use le_grids, only: new_trap_int, passing_wfb, trapped_wfb, mixed_wfb
    use le_grids, only: bouncefuzz, trapped_particles, wfbbc_option_switch, ngauss, npassing, nesuper
    use le_grids, only: genquad, radau_gauss_grid, split_passing_region, new_trap_int_split
    use le_grids, only: is_lower_bounce_point, is_upper_bounce_point, il_is_wfb, il_is_passing
    use le_grids, only: il_is_trapped, is_passing_hybrid_electron, grid_has_trapped_particles
    use theta_grid, only: ntgrid
    use kt_grids, only: naky
    use species, only: nspec
    implicit none
    integer, intent(in) :: ncid
    integer :: gid, ig, il, is, ik
    logical, dimension(:,:), allocatable :: lmask
    logical, dimension(:), allocatable :: lvec
    logical, dimension(:,:,:), allocatable :: lhyb
    gid = ex_group(ncid, 'le_grids')
    call wr_r2(gid, 'energy', energy, 'egrid species', 'le_grids::energy', 'le_grids.f90:285-296')
    call wr_r2(gid, 'w', w, 'egrid species', 'le_grids::w', 'le_grids.f90:285-296')
    call wr_r2(gid, 'speed', speed, 'egrid species', 'le_grids::speed', 'le_grids.f90:285-296')
    call wr_r1(gid, 'energy_maxw', energy_maxw, 'egrid', 'le_grids::energy_maxw', 'le_grids.f90:311-313')
    call wr_r1(gid, 'w_maxw', w_maxw, 'egrid', 'le_grids::w_maxw', 'le_grids.f90:311-313')
    call wr_r1(gid, 'speed_maxw', speed_maxw, 'egrid', 'le_grids::speed_maxw', 'le_grids.f90:311-313')
    call wr_r1(gid, 'al', al, 'lambda', 'le_grids::al', 'lgridset')
    call wr_r2(gid, 'wl', wl, 'theta lambda', 'le_grids::wl', 'lgridset (0 where forbidden)')
    call wr_i1(gid, 'jend', jend, 'theta', 'le_grids::jend', 'le_grids.f90:2871-2903')
    call wr_l2(gid, 'forbid', forbid, 'theta lambda', 'le_grids::forbid', 'le_grids.f90:2715-2738')
    call wr_l2(gid, 'is_ttp', is_ttp, 'theta lambda', 'le_grids::is_ttp', 'le_grids.f90:2907-2946')
    call wr_l2(gid, 'is_bounce_point', is_bounce_point, 'theta lambda', 'le_grids::is_bounce_point', &
         'le_grids.f90:2746-2771')
    call wr_l1(gid, 'can_be_ttp', can_be_ttp, 'lambda', 'le_grids::can_be_ttp', 'le_grids.f90:2907-2946')
    call wr_r2(gid, 'xi', xi, 'xi2 theta', 'le_grids::xi', 'le_grids.f90:2955-3020')
    call wr_r2(gid, 'wxi', wxi, 'xi2 theta', 'le_grids::wxi', 'le_grids.f90:2955-3020')
    call wr_i2(gid, 'ixi_to_il', ixi_to_il, 'xi2 theta', 'le_grids::ixi_to_il', 'le_grids.f90:2972-3011', ibase=1)
    call wr_i2(gid, 'ixi_to_isgn', ixi_to_isgn, 'xi2 theta', 'le_grids::ixi_to_isgn', &
         'le_grids.f90:2972-3011', ibase=1)
    call ex_vec(gid, 'sgn', sgn, 'two', fname='le_grids::sgn')
    ! Derived masks (public elemental functions, 09 §7.7)
    allocate (lmask(-ntgrid:ntgrid, nlambda))
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          lmask(ig, il) = is_lower_bounce_point(ig, il)
       end do
    end do
    call ex_vec(gid, 'is_lower_bounce_point', lmask, 'theta lambda', &
         fname='le_grids::is_lower_bounce_point(ig,il)', src='le_grids.f90:2775-2819')
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          lmask(ig, il) = is_upper_bounce_point(ig, il)
       end do
    end do
    call ex_vec(gid, 'is_upper_bounce_point', lmask, 'theta lambda', &
         fname='le_grids::is_upper_bounce_point(ig,il)', src='le_grids.f90:2823-2867')
    allocate (lvec(nlambda))
    do il = 1, nlambda
       lvec(il) = il_is_wfb(il)
    end do
    call ex_vec(gid, 'il_is_wfb', lvec, 'lambda', fname='le_grids::il_is_wfb(il)', src='le_grids.f90:484-491')
    do il = 1, nlambda
       lvec(il) = il_is_passing(il)
    end do
    call ex_vec(gid, 'il_is_passing', lvec, 'lambda', fname='le_grids::il_is_passing(il)', src='le_grids.f90:467-480')
    do il = 1, nlambda
       lvec(il) = il_is_trapped(il)
    end do
    call ex_vec(gid, 'il_is_trapped', lvec, 'lambda', fname='le_grids::il_is_trapped(il)', src='le_grids.f90:450-463')
    allocate (lhyb(nlambda, naky, nspec))
    do is = 1, nspec
       do ik = 1, naky
          do il = 1, nlambda
             lhyb(il, ik, is) = is_passing_hybrid_electron(is, ik, il)
          end do
       end do
    end do
    call ex_vec(gid, 'is_passing_hybrid_electron', lhyb, 'lambda ky species', &
         fname='le_grids::is_passing_hybrid_electron(is,ik,il)', src='le_grids.f90:4490-4499')
    call ex_scalar(gid, 'grid_has_trapped_particles', grid_has_trapped_particles())
    call ex_scalar(gid, 'negrid', negrid)
    call ex_scalar(gid, 'nlambda', nlambda)
    call ex_scalar(gid, 'ng2', ng2)
    call ex_scalar(gid, 'lmax', lmax)
    call ex_scalar(gid, 'nxi', nxi)
    call ex_scalar(gid, 'nesub', nesub)
    call ex_scalar(gid, 'nesuper', nesuper)
    call ex_scalar(gid, 'ngauss', ngauss)
    call ex_scalar(gid, 'npassing', npassing)
    call ex_scalar(gid, 'wdim', wdim)
    call ex_scalar(gid, 'vcut', vcut)
    call ex_scalar(gid, 'bouncefuzz', bouncefuzz)
    call ex_scalar(gid, 'wfbbc_option_switch', wfbbc_option_switch)
    call ex_scalar(gid, 'new_trap_int', new_trap_int)
    call ex_scalar(gid, 'new_trap_int_split', new_trap_int_split)
    call ex_scalar(gid, 'passing_wfb', passing_wfb)
    call ex_scalar(gid, 'trapped_wfb', trapped_wfb)
    call ex_scalar(gid, 'mixed_wfb', mixed_wfb)
    call ex_scalar(gid, 'trapped_particles', trapped_particles)
    call ex_scalar(gid, 'genquad', genquad)
    call ex_scalar(gid, 'radau_gauss_grid', radau_gauss_grid)
    call ex_scalar(gid, 'split_passing_region', split_passing_region)
  end subroutine export_le_grids

  !> /species (09 §3.4, 10 §A.4): one variable per specie component
  subroutine export_species(ncid)
    use species, only: nspec, spec, nonmaxw_corr, f0_values, dlnf0drho
    use species, only: has_electron_species, has_ion_species, has_hybrid_electron_species
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    character(len=*), parameter :: s = 'species.f90:448-480 (overrides 991-1013)'
    ! netCDF-4 cannot hold a group and a dimension both named 'species' in the root group (HDF5
    ! name clash), so this group is /species_ (root attribute group_renames)
    gid = ex_group(ncid, 'species_')
    call ex_scalar(gid, 'nspec', nspec)
    call ex_vec(gid, 'spec_z', spec(:)%z, 'species', fname='species::spec%z', src=s)
    call ex_vec(gid, 'spec_mass', spec(:)%mass, 'species', fname='species::spec%mass', src=s)
    call ex_vec(gid, 'spec_dens', spec(:)%dens, 'species', fname='species::spec%dens', src=s)
    call ex_vec(gid, 'spec_temp', spec(:)%temp, 'species', fname='species::spec%temp', src=s)
    call ex_vec(gid, 'spec_tprim', spec(:)%tprim, 'species', fname='species::spec%tprim', src=s)
    call ex_vec(gid, 'spec_fprim', spec(:)%fprim, 'species', fname='species::spec%fprim', src=s)
    call ex_vec(gid, 'spec_vnewk', spec(:)%vnewk, 'species', fname='species::spec%vnewk', src=s)
    call ex_vec(gid, 'spec_stm', spec(:)%stm, 'species', fname='species::spec%stm', src=s)
    call ex_vec(gid, 'spec_zstm', spec(:)%zstm, 'species', fname='species::spec%zstm', src=s)
    call ex_vec(gid, 'spec_tz', spec(:)%tz, 'species', fname='species::spec%tz', src=s)
    call ex_vec(gid, 'spec_smz', spec(:)%smz, 'species', fname='species::spec%smz', src=s)
    call ex_vec(gid, 'spec_zt', spec(:)%zt, 'species', fname='species::spec%zt', src=s)
    call ex_vec(gid, 'spec_uprim', spec(:)%uprim, 'species', fname='species::spec%uprim', src=s)
    call ex_vec(gid, 'spec_uprim2', spec(:)%uprim2, 'species', fname='species::spec%uprim2', src=s)
    call ex_vec(gid, 'spec_dens0', spec(:)%dens0, 'species', fname='species::spec%dens0', src=s)
    call ex_vec(gid, 'spec_u0', spec(:)%u0, 'species', fname='species::spec%u0', src=s)
    call ex_vec(gid, 'spec_tpar0', spec(:)%tpar0, 'species', fname='species::spec%tpar0', src=s)
    call ex_vec(gid, 'spec_tperp0', spec(:)%tperp0, 'species', fname='species::spec%tperp0', src=s)
    call ex_vec(gid, 'spec_nu_h', spec(:)%nu_h, 'species', fname='species::spec%nu_h', src=s)
    call ex_vec(gid, 'spec_bess_fac', spec(:)%bess_fac, 'species', fname='species::spec%bess_fac', src=s)
    call ex_vec(gid, 'spec_vcrit', spec(:)%vcrit, 'species', fname='species::spec%vcrit', src=s)
    call ex_vec(gid, 'spec_vcprim', spec(:)%vcprim, 'species', fname='species::spec%vcprim', src=s)
    call ex_vec(gid, 'spec_type', spec(:)%type, 'species', fname='species::spec%type', src=s)
    call ex_vec(gid, 'spec_f0type', spec(:)%f0type, 'species', fname='species::spec%f0type', src=s)
    call wr_r2(gid, 'nonmaxw_corr', nonmaxw_corr, 'egrid species', 'species::nonmaxw_corr', 'species.f90:563-597')
    call wr_r2(gid, 'f0_values', f0_values, 'egrid species', 'species::f0_values', 'species.f90:563-597')
    call wr_r2(gid, 'dlnf0drho', dlnf0drho, 'egrid species', 'species::dlnf0drho', 'species.f90:563-597')
    call ex_scalar(gid, 'has_electron_species', has_electron_species(spec))
    call ex_scalar(gid, 'has_ion_species', has_ion_species(spec))
    call ex_scalar(gid, 'has_hybrid_electron_species', has_hybrid_electron_species(spec))
  end subroutine export_species

  !> /run_parameters (09 §3.5, 10 §A.5)
  subroutine export_run_parameters(ncid)
    use run_parameters, only: beta, zeff, tite, fphi, fapar, fbpar, k0, rhostar, delt, max_sim_time
    use run_parameters, only: avail_cpu_time, margin_cpu_time, has_phi, has_apar, has_bpar
    use run_parameters, only: nstep, ncheck_stop, progress_frequency, delt_option_switch
    use run_parameters, only: immediate_reset, wstar_units, neo_test, do_eigsolve, use_old_diagnostics
    use run_parameters, only: reset, ieqzip, eqzip_option_switch, seed
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    gid = ex_group(ncid, 'run_parameters')
    call ex_scalar(gid, 'beta', beta)
    call ex_scalar(gid, 'zeff', zeff)
    call ex_scalar(gid, 'tite', tite)
    call ex_scalar(gid, 'fphi', fphi)
    call ex_scalar(gid, 'fapar', fapar, src='run_parameters.f90:340-349 (beta==0 forces fapar=fbpar=0)')
    call ex_scalar(gid, 'fbpar', fbpar, src='run_parameters.f90:340-349 (beta==0 forces fapar=fbpar=0)')
    call ex_scalar(gid, 'k0', k0)
    call ex_scalar(gid, 'rhostar', rhostar)
    call ex_scalar(gid, 'delt', delt)
    call ex_scalar(gid, 'max_sim_time', max_sim_time)
    call ex_scalar(gid, 'avail_cpu_time', avail_cpu_time)
    call ex_scalar(gid, 'margin_cpu_time', margin_cpu_time)
    call ex_scalar(gid, 'has_phi', has_phi, src='run_parameters.f90:351-353')
    call ex_scalar(gid, 'has_apar', has_apar, src='run_parameters.f90:351-353')
    call ex_scalar(gid, 'has_bpar', has_bpar, src='run_parameters.f90:351-353')
    call ex_scalar(gid, 'nstep', nstep)
    call ex_scalar(gid, 'ncheck_stop', ncheck_stop)
    call ex_scalar(gid, 'progress_frequency', progress_frequency)
    call ex_scalar(gid, 'delt_option_switch', delt_option_switch)
    call ex_scalar(gid, 'immediate_reset', immediate_reset)
    call ex_scalar(gid, 'wstar_units', wstar_units)
    call ex_scalar(gid, 'neo_test', neo_test)
    call ex_scalar(gid, 'do_eigsolve', do_eigsolve)
    call ex_scalar(gid, 'use_old_diagnostics', use_old_diagnostics)
    call ex_scalar(gid, 'reset', reset)
    call wr_l2(gid, 'ieqzip', ieqzip, 'kx ky', 'run_parameters::ieqzip', 'run_parameters.f90:259-272')
    call ex_scalar(gid, 'eqzip_option_switch', eqzip_option_switch)
    call ex_scalar(gid, 'seed', seed)
  end subroutine export_run_parameters

  !> The time block (09 §7.6) as 0-D variables in group gid
  subroutine write_time_block(gid, istep)
    use gs2_time, only: code_dt, user_dt, code_dt_prev1, code_dt_prev2, code_dt_old, code_dt_cfl
    use gs2_time, only: user_dt_cfl, code_dt_min, code_dt_max, code_time, user_time
    use run_parameters, only: reset
    use collisions, only: vnmult
    use nonlinear_terms, only: nl_istep_last => istep_last, nl_initializing => initializing
    use dist_fn, only: exb_first
    implicit none
    integer, intent(in) :: gid, istep
    call ex_scalar(gid, 'code_dt', code_dt, fname='gs2_time::code_dt')
    call ex_scalar(gid, 'user_dt', user_dt, fname='gs2_time::user_dt')
    call ex_scalar(gid, 'code_dt_prev1', code_dt_prev1, fname='gs2_time::code_dt_prev1')
    call ex_scalar(gid, 'code_dt_prev2', code_dt_prev2, fname='gs2_time::code_dt_prev2')
    call ex_scalar(gid, 'code_dt_old', code_dt_old, fname='gs2_time::code_dt_old')
    call ex_scalar(gid, 'code_dt_cfl', code_dt_cfl, fname='gs2_time::code_dt_cfl')
    call ex_scalar(gid, 'user_dt_cfl', user_dt_cfl, fname='gs2_time::user_dt_cfl')
    call ex_scalar(gid, 'code_dt_min', code_dt_min, fname='gs2_time::code_dt_min')
    call ex_scalar(gid, 'code_dt_max', code_dt_max, fname='gs2_time::code_dt_max')
    call ex_scalar(gid, 'code_time', code_time, fname='gs2_time::code_time')
    call ex_scalar(gid, 'user_time', user_time, fname='gs2_time::user_time')
    call ex_scalar(gid, 'istep', istep)
    call ex_scalar(gid, 'reset', reset, fname='run_parameters::reset')
    call ex_vec(gid, 'vnmult', vnmult, 'two', fname='collisions::vnmult')
    call ex_scalar(gid, 'nl_istep_last', nl_istep_last, fname='nonlinear_terms::istep_last')
    call ex_scalar(gid, 'nl_initializing', nl_initializing, fname='nonlinear_terms::initializing')
    call ex_scalar(gid, 'exb_first', exb_first, fname='dist_fn::exb_first')
  end subroutine write_time_block

  !> /gs2_time (09 §3.5): time normalisation and the time block at export time
  subroutine export_gs2_time(ncid, path)
    use gs2_time, only: tunits, wunits, woutunits, dt_not_set
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: path
    integer :: gid
    gid = ex_group(ncid, path)
    call wr_r1(gid, 'tunits', tunits, 'ky', 'gs2_time::tunits', 'gs2_time.f90:83-126')
    call wr_r1(gid, 'wunits', wunits, 'ky', 'gs2_time::wunits', 'gs2_time.f90:83-126')
    call wr_r1(gid, 'woutunits', woutunits, 'ky', 'gs2_time::woutunits', 'gs2_time.f90:83-126')
    call ex_scalar(gid, 'dt_not_set', dt_not_set, fname='gs2_time::dt_not_set')
    call ex_att(gid, '', 'dt_not_set', dt_not_set)
    call write_time_block(gid, -1)
  end subroutine export_gs2_time

  !> /gs2_reinit (09 §3.5): config read before any ladder rung (G-20)
  subroutine export_gs2_reinit(ncid)
    use gs2_reinit, only: delt_adj, in_memory, get_reinit_config, reinit_config_type
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    type(reinit_config_type) :: cfg
    gid = ex_group(ncid, 'gs2_reinit')
    cfg = get_reinit_config()
    call ex_scalar(gid, 'delt_adj', delt_adj)
    call ex_scalar(gid, 'in_memory', in_memory)
    call ex_scalar(gid, 'delt_cushion', cfg%delt_cushion, src='get_reinit_config()')
    call ex_scalar(gid, 'delt_minimum', cfg%delt_minimum, src='get_reinit_config()')
    call ex_scalar(gid, 'dt0', cfg%dt0, src='get_reinit_config() (smart default code_dt_max resolved)')
    call ex_scalar(gid, 'abort_rapid_time_step_change', cfg%abort_rapid_time_step_change, &
         src='get_reinit_config()')
  end subroutine export_gs2_reinit

  !> /init_g (09 §3.5)
  subroutine export_init_g(ncid)
    use init_g, only: tstart, new_field_init, initial_condition_is_nonadiabatic_dfn, restart_file
    use init_g, only: ginitopt_switch, get_init_g_config, init_g_config_type
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    type(init_g_config_type) :: cfg
    gid = ex_group(ncid, 'init_g')
    cfg = get_init_g_config()
    call ex_scalar(gid, 'tstart', tstart)
    call ex_scalar(gid, 'new_field_init', new_field_init)
    call ex_scalar(gid, 'initial_condition_is_nonadiabatic_dfn', initial_condition_is_nonadiabatic_dfn)
    call ex_scalar(gid, 'ginitopt_switch', ginitopt_switch)
    call ex_scalar(gid, 'scale', cfg%scale, src='get_init_g_config() (input value)')
    call ex_scalar(gid, 'read_many', cfg%read_many, src='get_init_g_config()')
    call ex_scalar(gid, 'include_explicit_source_in_restart', cfg%include_explicit_source_in_restart, &
         src='get_init_g_config()')
    call ex_att(gid, '', 'ginit_option', trim(cfg%ginit_option))
    call ex_att(gid, '', 'restart_file', trim(restart_file))
    call ex_att(gid, '', 'restart_dir', trim(cfg%restart_dir))
  end subroutine export_init_g

  !> /fields (09 §3.5)
  subroutine export_fields(ncid)
    use fields, only: fieldopt_switch, remove_zonal_flows_switch, force_maxwell_reinit
    use fields, only: get_fields_config, fields_config_type
    use fields_arrays, only: response_file
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    type(fields_config_type) :: cfg
    gid = ex_group(ncid, 'fields')
    cfg = get_fields_config()
    call ex_scalar(gid, 'fieldopt_switch', fieldopt_switch)
    call ex_scalar(gid, 'remove_zonal_flows_switch', remove_zonal_flows_switch)
    call ex_scalar(gid, 'force_maxwell_reinit', force_maxwell_reinit)
    call ex_scalar(gid, 'dump_response', cfg%dump_response)
    call ex_scalar(gid, 'read_response', cfg%read_response)
    call ex_scalar(gid, 'field_subgath', cfg%field_subgath)
    call ex_scalar(gid, 'do_smart_update', cfg%do_smart_update)
    call ex_scalar(gid, 'minnrow', cfg%minnrow)
    call ex_att(gid, '', 'field_option', trim(cfg%field_option))
    call ex_att(gid, '', 'response_dir', trim(cfg%response_dir))
    call ex_att(gid, '', 'response_file', trim(response_file))
  end subroutine export_fields

  !> /gs2_layouts (09 §2.1-2.2, §3.6): orders and the summation-order maps iglo_of, ile_of
  subroutine export_gs2_layouts(ncid)
    use gs2_layouts, only: g_lo, le_lo, layout, idx, local_field_solve, intmom_sub, intspec_sub
    use gs2_layouts, only: gf_local_fields, opt_local_copy, fft_use_wisdom, fft_measure_plan
    use gs2_layouts, only: fft_wisdom_file, simple_gf_decomposition
    use gs2_layouts, only: get_layouts_config, layouts_config_type
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    implicit none
    integer, intent(in) :: ncid
    integer :: gid, it, ik, il, ie, is, ig
    integer, dimension(:,:,:,:,:), allocatable :: iglo_of
    integer, dimension(:,:,:,:), allocatable :: ile_of
    type(layouts_config_type) :: cfg
    gid = ex_group(ncid, 'gs2_layouts')
    cfg = get_layouts_config()
    call ex_att(gid, '', 'layout', trim(layout))
    call ex_scalar(gid, 'g_lo_ik_ord', g_lo%ik_ord)
    call ex_scalar(gid, 'g_lo_it_ord', g_lo%it_ord)
    call ex_scalar(gid, 'g_lo_il_ord', g_lo%il_ord)
    call ex_scalar(gid, 'g_lo_ie_ord', g_lo%ie_ord)
    call ex_scalar(gid, 'g_lo_is_ord', g_lo%is_ord)
    call ex_scalar(gid, 'g_lo_ulim_world', g_lo%ulim_world)
    call ex_scalar(gid, 'g_lo_blocksize', g_lo%blocksize)
    call ex_scalar(gid, 'g_lo_nproc', g_lo%nproc)
    ! le_lo is only initialised when collisions are set up (init_le_layouts)
    if (le_lo%ig_ord > 0) then
       call ex_scalar(gid, 'le_lo_ig_ord', le_lo%ig_ord)
       call ex_scalar(gid, 'le_lo_it_ord', le_lo%it_ord)
       call ex_scalar(gid, 'le_lo_ik_ord', le_lo%ik_ord)
       call ex_scalar(gid, 'le_lo_is_ord', le_lo%is_ord)
       call ex_scalar(gid, 'le_lo_ulim_world', le_lo%ulim_world)
    else
       call ex_absent(gid, 'le_lo_ig_ord le_lo_it_ord le_lo_ik_ord le_lo_is_ord le_lo_ulim_world ile_of')
    end if
    allocate (iglo_of(ntheta0, naky, nlambda, negrid, nspec))
    do is = 1, nspec ; do ie = 1, negrid ; do il = 1, nlambda ; do ik = 1, naky ; do it = 1, ntheta0
       iglo_of(it, ik, il, ie, is) = idx(g_lo, ik, it, il, ie, is)
    end do ; end do ; end do ; end do ; end do
    call put_i5(gid, 'iglo_of', iglo_of)
    if (le_lo%ig_ord > 0) then
       allocate (ile_of(-ntgrid:ntgrid, ntheta0, naky, nspec))
       do is = 1, nspec ; do ik = 1, naky ; do it = 1, ntheta0 ; do ig = -ntgrid, ntgrid
          ile_of(ig, it, ik, is) = idx(le_lo, ig, ik, it, is)
       end do ; end do ; end do ; end do
       call put_i4(gid, 'ile_of', ile_of)
    end if
    call ex_scalar(gid, 'local_field_solve', local_field_solve)
    call ex_scalar(gid, 'intmom_sub', intmom_sub)
    call ex_scalar(gid, 'intspec_sub', intspec_sub)
    call ex_scalar(gid, 'gf_local_fields', gf_local_fields)
    call ex_scalar(gid, 'opt_local_copy', opt_local_copy)
    call ex_scalar(gid, 'fft_use_wisdom', fft_use_wisdom)
    call ex_scalar(gid, 'fft_measure_plan', fft_measure_plan)
    call ex_scalar(gid, 'simple_gf_decomposition', simple_gf_decomposition)
    call ex_att(gid, '', 'fft_wisdom_file', trim(fft_wisdom_file))
    call ex_att(gid, '', 'nproc_g_lo', cfg%nproc_g_lo)
    call ex_att(gid, '', 'nproc_le_lo', cfg%nproc_le_lo)
    call ex_att(gid, '', 'nproc_lz_lo', cfg%nproc_lz_lo)
    call ex_att(gid, '', 'nproc_e_lo', cfg%nproc_e_lo)
    call ex_att(gid, '', 'nproc_xxf_lo', cfg%nproc_xxf_lo)
    call ex_att(gid, '', 'nproc_yxf_lo', cfg%nproc_yxf_lo)
  contains
    subroutine put_i5(gid, name, arr)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name
      integer, dimension(:,:,:,:,:), intent(in) :: arr
      call ex_vec(gid, name, reshape(arr, [size(arr)]), 'kx ky lambda egrid species', &
           fname='idx(g_lo, ik, it, il, ie, is)', src='gs2_layouts.fpp:1606-1617', ibase=0, tag='G5')
    end subroutine put_i5
    subroutine put_i4(gid, name, arr)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name
      integer, dimension(:,:,:,:), intent(in) :: arr
      call ex_vec(gid, name, reshape(arr, [size(arr)]), 'theta kx ky species', &
           fname='idx(le_lo, ig, ik, it, is)', src='gs2_layouts.fpp:2160-2170', ibase=0, tag='V')
    end subroutine put_i4
  end subroutine export_gs2_layouts

  !> /dist_fn_arrays static part (09 §3.7): vpa, vpac, vperp2, aj0, aj1, wstar
  subroutine export_dist_fn_arrays(ncid, knobs)
    use dist_fn_arrays, only: vpa, vpac, vperp2, aj0, aj1
    use gs2_layouts, only: g_lo
    implicit none
    integer, intent(in) :: ncid
    type(export_knobs_type), intent(in) :: knobs
    integer :: gid
    gid = ex_group(ncid, 'dist_fn_arrays')
    call wr_g7r(gid, 'vpa', vpa, 'dist_fn_arrays::vpa', 'dist_fn.fpp:1592-1678 init_vpa_vpac', '02')
    call wr_g7r(gid, 'vpac', vpac, 'dist_fn_arrays::vpac', 'dist_fn.fpp:1592-1678 init_vpa_vpac', '02')
    call wr_g6r(gid, 'vperp2', vperp2, 'dist_fn_arrays::vperp2', 'dist_fn.fpp:1560-1589 init_vperp2', '02')
    call wr_g6r(gid, 'aj0', aj0, 'dist_fn_arrays::aj0', 'dist_fn.fpp:1756-1823 init_bessel', '02')
    call wr_g6r(gid, 'aj1', aj1, 'dist_fn_arrays::aj1', 'dist_fn.fpp:1756-1823 init_bessel (J1(x)/x)', '02')
    call export_wstar(gid)
    if (knobs%write_raw_glo) call ex_att(gid, '', 'raw_glo_note', 'see /raw for glo-ordered copies')
    call ex_att(gid, '', 'g_lo_ulim_world', g_lo%ulim_world)
  end subroutine export_dist_fn_arrays

  subroutine export_wstar(gid)
    use dist_fn_arrays, only: wstar
    implicit none
    integer, intent(in) :: gid
    call wr_r3(gid, 'wstar', wstar, 'ky egrid species', 'dist_fn_arrays::wstar', &
         'dist_fn.fpp:1732-1750 init_wstar', '02', dtdep=.true.)
  end subroutine export_wstar

  !> dt-dependent dist_fn coefficients (10 §A.8) into group gid
  subroutine export_dist_fn_coeffs(gid)
    use dist_fn, only: vpar, wdrift, wdriftttp, dfa => a, dfb => b, dfr => r, ainv, g_h
    use dist_fn, only: source_coeffs_phim, source_coeffs_phip, source_coeffs_aparm, source_coeffs_aparp
    implicit none
    integer, intent(in) :: gid
    call wr_g7r(gid, 'vpar', vpar, 'dist_fn::vpar', 'dist_fn.fpp:1681-1729 init_vpar', '02', .true.)
    call wr_g7r(gid, 'wdrift', wdrift, 'dist_fn::wdrift', 'dist_fn.fpp:1405-1498 init_wdrift', '02', .true.)
    call wr_g7r(gid, 'wdriftttp', wdriftttp, 'dist_fn::wdriftttp', 'dist_fn.fpp:1451-1480', '02', .true.)
    call wr_g7c(gid, 'a', dfa, 'dist_fn::a', 'dist_fn.fpp:1892-1897 init_invert_rhs', '03', .true.)
    call wr_g7c(gid, 'b', dfb, 'dist_fn::b', 'dist_fn.fpp:1892-1897 init_invert_rhs', '03', .true.)
    call wr_g7c(gid, 'r', dfr, 'dist_fn::r', 'dist_fn.fpp:1885-1891 init_invert_rhs', '03', .true.)
    call wr_g7c(gid, 'ainv', ainv, 'dist_fn::ainv', 'dist_fn.fpp:1899-1921 init_invert_rhs', '03', .true.)
    call wr_g7c(gid, 'g_h', g_h, 'dist_fn::g_h', 'dist_fn.fpp:1937-2125 init_homogeneous_g', '03', .true.)
    call wr_g7m1c(gid, 'source_coeffs_phim', source_coeffs_phim, 'dist_fn::source_coeffs_phim', &
         'dist_fn.fpp:2194-2264 init_source_term (opt_source)', '02', .true.)
    call wr_g7m1c(gid, 'source_coeffs_phip', source_coeffs_phip, 'dist_fn::source_coeffs_phip', &
         'dist_fn.fpp:2194-2264 init_source_term (opt_source)', '02', .true.)
    call wr_g7m1c(gid, 'source_coeffs_aparm', source_coeffs_aparm, 'dist_fn::source_coeffs_aparm', &
         'dist_fn.fpp:2194-2264 init_source_term (opt_source)', '02', .true.)
    call wr_g7m1c(gid, 'source_coeffs_aparp', source_coeffs_aparp, 'dist_fn::source_coeffs_aparp', &
         'dist_fn.fpp:2194-2264 init_source_term (opt_source)', '02', .true.)
  end subroutine export_dist_fn_coeffs

  !> /dist_fn (09 §3.8, 10 §A.7-A.9)
  subroutine export_dist_fn(ncid, knobs)
    use dist_fn, only: g_adj, ufac, gamtot, gamtot1, gamtot2, gamtot3, gridfac1
    use dist_fn, only: inv_phi_denominator_g, inv_bpar_denominator_g, awgt
    use dist_fn, only: itleft, itright, l_links, r_links, n_links, M_class, N_class, i_class
    use dist_fn, only: supercell_labels, n_cells, n_supercells, connections, save_h, no_connections
    use dist_fn, only: jump, ikx_indexed, exb_first, get_leftmost_it, get_rightmost_it
    use dist_fn, only: boundary_option_switch, boundary_option_linked, nonad_zero
    use dist_fn, only: fexp, bkdiff, gridfac, apfac, driftknob, tpdriftknob, poisfac, vparknob
    use dist_fn, only: t0, omega0, gamma0, source0, phi_ext, afilter, exponential_boundary_factor
    use dist_fn, only: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
    use dist_fn, only: g_exb_start_time, g_exb_error_limit, g_exb_start_timestep, hyper_in_initialisation
    use dist_fn, only: start_from_previous_solution, boundary_off_grid, exponential_boundary
    use dist_fn, only: lf_default, lf_decompose, esv, opt_source, source_option_switch
    use dist_fn, only: adiabatic_option_switch, def_parity, even, zero_forbid, gf_lo_integrate, mult_imp
    use gs2_layouts, only: g_lo, it_idx
    use kt_grids, only: naky, ntheta0
    use le_grids, only: mixed_wfb
    implicit none
    integer, intent(in) :: ncid
    type(export_knobs_type), intent(in) :: knobs
    integer :: gid, iglo, it, ik, nloc
    integer, dimension(:), allocatable :: ibuf
    logical, dimension(:), allocatable :: lbuf
    integer, dimension(:,:), allocatable :: ik2
    gid = ex_group(ncid, 'dist_fn')
    call export_dist_fn_coeffs(gid)
    if (knobs%write_g_adj .and. allocated(g_adj)) then
       call ex_glc(gid, 'g_adj', g_adj, fname='dist_fn::g_adj', src='dist_fn.fpp:2120-2123 (h-slots at init)', &
            owner='03', dtdep=.true.)
    else
       call ex_absent(gid, 'g_adj')
    end if
    call wr_r2(gid, 'ufac', ufac, 'egrid species', 'dist_fn::ufac', 'dist_fn.fpp:2180-2187 init_source_term', '02')
    call wr_r3(gid, 'gamtot', gamtot, 'theta kx ky', 'dist_fn::gamtot', 'dist_fn.fpp:6470-6572 init_fieldeq', '04')
    call wr_r3(gid, 'gamtot1', gamtot1, 'theta kx ky', 'dist_fn::gamtot1', 'dist_fn.fpp:6470-6572 init_fieldeq', '04')
    call wr_r3(gid, 'gamtot2', gamtot2, 'theta kx ky', 'dist_fn::gamtot2', 'dist_fn.fpp:6470-6572 init_fieldeq', '04')
    call wr_r3(gid, 'gamtot3', gamtot3, 'theta kx ky', 'dist_fn::gamtot3', 'dist_fn.fpp:6564-6565 (fieldlineavg)', '04')
    call wr_r3(gid, 'gridfac1', gridfac1, 'theta kx ky', 'dist_fn::gridfac1', 'dist_fn.fpp:6486-6506', '04')
    call wr_r3(gid, 'inv_phi_denominator_g', inv_phi_denominator_g, 'theta kx ky', &
         'dist_fn::inv_phi_denominator_g', 'dist_fn.fpp:7457-7585 (lazy)', '04')
    call wr_r3(gid, 'inv_bpar_denominator_g', inv_bpar_denominator_g, 'theta kx ky', &
         'dist_fn::inv_bpar_denominator_g', 'dist_fn.fpp:7457-7585 (lazy)', '04')
    call wr_r2(gid, 'awgt', awgt, 'kx ky', 'dist_fn::awgt', 'dist_fn.fpp:8006-8020 (lazy, frozen)', '04')
    call wr_i2(gid, 'itleft', itleft, 'kx ky', 'dist_fn::itleft', 'dist_fn.fpp:2441-2513', ibase=1)
    call wr_i2(gid, 'itright', itright, 'kx ky', 'dist_fn::itright', 'dist_fn.fpp:2441-2513', ibase=1)
    call wr_i2(gid, 'l_links', l_links, 'kx ky', 'dist_fn::l_links', 'dist_fn.fpp:2608-2672 count_links')
    call wr_i2(gid, 'r_links', r_links, 'kx ky', 'dist_fn::r_links', 'dist_fn.fpp:2608-2672 count_links')
    call wr_i3(gid, 'n_links', n_links, 'two kx ky', 'dist_fn::n_links', 'dist_fn.fpp:2608-2672 count_links')
    call wr_i1(gid, 'M_class', M_class, 'n_class', 'dist_fn::M_class', 'dist_fn.fpp:2676-2755')
    call wr_i1(gid, 'N_class', N_class, 'n_class', 'dist_fn::N_class', 'dist_fn.fpp:2676-2755')
    call ex_scalar(gid, 'i_class', i_class)
    call wr_i2(gid, 'supercell_labels', supercell_labels, 'kx ky', 'dist_fn::supercell_labels', &
         'dist_fn.fpp:2355-2397')
    call wr_i2(gid, 'n_cells', n_cells, 'kx ky', 'dist_fn::n_cells', 'dist_fn.fpp:2355-2397')
    call wr_i1(gid, 'n_supercells', n_supercells, 'ky', 'dist_fn::n_supercells', 'dist_fn.fpp:2355-2397')
    ! Topology helpers (09 §7.7)
    allocate (ik2(ntheta0, naky))
    do ik = 1, naky
       do it = 1, ntheta0
          ik2(it, ik) = get_leftmost_it(it, ik)
       end do
    end do
    call ex_vec(gid, 'itleftmost', ik2, 'kx ky', fname='dist_fn::get_leftmost_it(it,ik)', &
         src='dist_fn.fpp:3889-3903', ibase=1)
    do ik = 1, naky
       do it = 1, ntheta0
          ik2(it, ik) = get_rightmost_it(it, ik)
       end do
    end do
    call ex_vec(gid, 'itrightmost', ik2, 'kx ky', fname='dist_fn::get_rightmost_it(it,ik)', &
         src='dist_fn.fpp:3906-3920', ibase=1)
    if (allocated(connections)) then
       nloc = g_lo%ulim_alloc - g_lo%llim_proc + 1
       allocate (lbuf(nloc), ibuf(nloc))
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          lbuf(iglo - g_lo%llim_proc + 1) = connections(iglo)%neighbor
       end do
       call ex_g5l(gid, 'connections_neighbor', lbuf, fname='dist_fn::connections%neighbor', &
            src='dist_fn.fpp:2518-2565')
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          if (connections(iglo)%iglo_left >= 0) then
             ibuf(iglo - g_lo%llim_proc + 1) = it_idx(g_lo, connections(iglo)%iglo_left)
          else
             ibuf(iglo - g_lo%llim_proc + 1) = -1
          end if
       end do
       call ex_g5i(gid, 'connections_it_left', ibuf, fname='it_idx(g_lo, connections%iglo_left) or -1', &
            src='dist_fn.fpp:2518-2565', ibase=1)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          if (connections(iglo)%iglo_right >= 0) then
             ibuf(iglo - g_lo%llim_proc + 1) = it_idx(g_lo, connections(iglo)%iglo_right)
          else
             ibuf(iglo - g_lo%llim_proc + 1) = -1
          end if
       end do
       call ex_g5i(gid, 'connections_it_right', ibuf, fname='it_idx(g_lo, connections%iglo_right) or -1', &
            src='dist_fn.fpp:2518-2565', ibase=1)
    else
       call ex_absent(gid, 'connections_neighbor')
       call ex_absent(gid, 'connections_it_left')
       call ex_absent(gid, 'connections_it_right')
    end if
    if (allocated(save_h)) then
       call ex_g5sl(gid, 'save_h', save_h, fname='dist_fn::save_h', src='dist_fn.fpp:2570-2604')
    else
       call ex_absent(gid, 'save_h')
    end if
    call ex_scalar(gid, 'no_connections', no_connections)
    call wr_i1(gid, 'jump', jump, 'ky', 'dist_fn::jump', 'dist_fn.fpp:4237-4249 (first shear call)')
    call wr_i1(gid, 'ikx_indexed', ikx_indexed, 'kx', 'dist_fn::ikx_indexed', 'dist_fn.fpp:4237-4249')
    call ex_scalar(gid, 'exb_first', exb_first)
    ! Mixed-wfb end Bessel factors (09 §7.7, 10 G-04)
    if (boundary_option_switch == boundary_option_linked .and. mixed_wfb .and. nonad_zero) then
       call export_wfb_end_factors(gid)
    else
       call ex_absent(gid, 'wfb_aj0l')
       call ex_absent(gid, 'wfb_aj1l')
       call ex_absent(gid, 'wfb_aj0r')
       call ex_absent(gid, 'wfb_aj1r')
       call ex_absent(gid, 'wfb_vperp2left')
       call ex_absent(gid, 'wfb_vperp2right')
    end if
    ! Effective knobs (read_parameters; forced values already applied)
    call wr_r1(gid, 'fexp', fexp, 'species', 'dist_fn::fexp', 'dist_fn.fpp read_species_knobs')
    call wr_r1(gid, 'bkdiff', bkdiff, 'species', 'dist_fn::bkdiff', 'dist_fn.fpp:1319-1344 (forced equal unless mult_imp)')
    call ex_scalar(gid, 'gridfac', gridfac)
    call ex_scalar(gid, 'apfac', apfac)
    call ex_scalar(gid, 'driftknob', driftknob)
    call ex_scalar(gid, 'tpdriftknob', tpdriftknob)
    call ex_scalar(gid, 'poisfac', poisfac)
    call ex_scalar(gid, 'vparknob', vparknob)
    call ex_scalar(gid, 't0', t0)
    call ex_scalar(gid, 'omega0', omega0)
    call ex_scalar(gid, 'gamma0', gamma0)
    call ex_scalar(gid, 'source0', source0)
    call ex_scalar(gid, 'phi_ext', phi_ext)
    call ex_scalar(gid, 'afilter', afilter)
    call ex_scalar(gid, 'exponential_boundary_factor', exponential_boundary_factor)
    call ex_scalar(gid, 'wfb', wfb)
    call ex_scalar(gid, 'g_exb', g_exb)
    call ex_scalar(gid, 'g_exbfac', g_exbfac)
    call ex_scalar(gid, 'omprimfac', omprimfac)
    call ex_scalar(gid, 'btor_slab', btor_slab)
    call ex_scalar(gid, 'mach', mach)
    call ex_scalar(gid, 'g_exb_start_time', g_exb_start_time)
    call ex_scalar(gid, 'g_exb_error_limit', g_exb_error_limit)
    call ex_scalar(gid, 'g_exb_start_timestep', g_exb_start_timestep)
    call ex_scalar(gid, 'hyper_in_initialisation', hyper_in_initialisation)
    call ex_scalar(gid, 'start_from_previous_solution', start_from_previous_solution)
    call ex_scalar(gid, 'boundary_off_grid', boundary_off_grid)
    call ex_scalar(gid, 'exponential_boundary', exponential_boundary)
    call ex_scalar(gid, 'nonad_zero', nonad_zero)
    call ex_scalar(gid, 'lf_default', lf_default)
    call ex_scalar(gid, 'lf_decompose', lf_decompose)
    call ex_scalar(gid, 'esv', esv, src='dist_fn.fpp:1327 (esv only with linked boundaries)')
    call ex_scalar(gid, 'opt_source', opt_source)
    call ex_scalar(gid, 'source_option_switch', source_option_switch)
    call ex_scalar(gid, 'boundary_option_switch', boundary_option_switch)
    call ex_scalar(gid, 'adiabatic_option_switch', adiabatic_option_switch)
    call ex_scalar(gid, 'def_parity', def_parity)
    call ex_scalar(gid, 'even', even)
    call ex_scalar(gid, 'zero_forbid', zero_forbid)
    call ex_scalar(gid, 'gf_lo_integrate', gf_lo_integrate)
    call ex_scalar(gid, 'mult_imp', mult_imp)
  end subroutine export_dist_fn

  !> Verbatim recompute of the mixed-wfb end Bessel factors of
  !> apply_linked_boundary_conditions (dist_fn.fpp:6016-6028), pre-gathered at the chain ends
  !> itl = get_leftmost_it(it,ik), itr = get_rightmost_it(it,ik) (10 G-04)
  subroutine export_wfb_end_factors(gid)
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy, al, il_is_wfb, nlambda, negrid
    use kt_grids, only: kperp2, naky, ntheta0
    use species, only: spec, nspec
    use spfunc, only: j0, j1
    use dist_fn, only: get_leftmost_it, get_rightmost_it
    implicit none
    integer, intent(in) :: gid
    real, dimension(:,:,:,:), allocatable :: aj0l, aj1l, aj0r, aj1r
    real, dimension(:,:), allocatable :: vperp2left, vperp2right
    real :: argl, argr
    integer :: il, ie, is, ik, it, itl, itr, iwfb
    iwfb = -1
    do il = 1, nlambda
       if (il_is_wfb(il)) iwfb = il
    end do
    if (iwfb < 0) then
       call ex_absent(gid, 'wfb_aj0l')
       call ex_absent(gid, 'wfb_aj1l')
       call ex_absent(gid, 'wfb_aj0r')
       call ex_absent(gid, 'wfb_aj1r')
       call ex_absent(gid, 'wfb_vperp2left')
       call ex_absent(gid, 'wfb_vperp2right')
       return
    end if
    il = iwfb
    allocate (aj0l(negrid, ntheta0, naky, nspec), aj1l(negrid, ntheta0, naky, nspec))
    allocate (aj0r(negrid, ntheta0, naky, nspec), aj1r(negrid, ntheta0, naky, nspec))
    allocate (vperp2left(negrid, nspec), vperp2right(negrid, nspec))
    do is = 1, nspec
       do ie = 1, negrid
          vperp2left(ie, is) = bmag(-ntgrid)*al(il)*energy(ie,is)
          vperp2right(ie, is) = bmag(ntgrid)*al(il)*energy(ie,is)
          do ik = 1, naky
             do it = 1, ntheta0
                itl=get_leftmost_it(it,ik)
                itr=get_rightmost_it(it,ik)
                argl = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(-ntgrid)*kperp2(-ntgrid,itl,ik))
                argr = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(ntgrid)*kperp2(ntgrid,itr,ik))
                aj0l(ie, it, ik, is) = j0(argl)
                aj0r(ie, it, ik, is) = j0(argr)
                aj1l(ie, it, ik, is) = j1(argl)
                aj1r(ie, it, ik, is) = j1(argr)
             end do
          end do
       end do
    end do
    call ex_vec(gid, 'wfb_aj0l', aj0l, 'egrid kx ky species', fname='aj0l (local of apply_linked_boundary_conditions)', &
         src='dist_fn.fpp:6016-6028 (verbatim, at itl=get_leftmost_it)', owner='03')
    call ex_vec(gid, 'wfb_aj1l', aj1l, 'egrid kx ky species', fname='aj1l (local of apply_linked_boundary_conditions)', &
         src='dist_fn.fpp:6016-6028 (verbatim, at itl=get_leftmost_it)', owner='03')
    call ex_vec(gid, 'wfb_aj0r', aj0r, 'egrid kx ky species', fname='aj0r (local of apply_linked_boundary_conditions)', &
         src='dist_fn.fpp:6016-6028 (verbatim, at itr=get_rightmost_it)', owner='03')
    call ex_vec(gid, 'wfb_aj1r', aj1r, 'egrid kx ky species', fname='aj1r (local of apply_linked_boundary_conditions)', &
         src='dist_fn.fpp:6016-6028 (verbatim, at itr=get_rightmost_it)', owner='03')
    call ex_vec(gid, 'wfb_vperp2left', vperp2left, 'egrid species', fname='vperp2left', &
         src='dist_fn.fpp:6019', owner='03')
    call ex_vec(gid, 'wfb_vperp2right', vperp2right, 'egrid species', fname='vperp2right', &
         src='dist_fn.fpp:6020', owner='03')
    call ex_att(gid, 'wfb_aj0l', 'il_wfb', il)
  end subroutine export_wfb_end_factors

  !> Field-solve matrices: /fields_implicit (09 §7.4) or /fields_local (09 §7.5)
  subroutine export_field_matrices(ncid, prefix, knobs)
    use fields, only: fieldopt_switch, fieldopt_implicit, fieldopt_local
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prefix
    type(export_knobs_type), intent(in) :: knobs
    integer :: gid
    if (fieldopt_switch == fieldopt_implicit) then
       gid = ex_group(ncid, trim(prefix)//'fields_implicit')
       call export_fields_implicit(gid, knobs)
    else if (fieldopt_switch == fieldopt_local) then
       gid = ex_group(ncid, trim(prefix)//'fields_local')
       call export_fields_local(ncid, gid, knobs)
    end if
  end subroutine export_field_matrices

  !> /fields_implicit: raw runtime inverse response rows aminv_c<c> (09 §7.4) + class tables
  subroutine export_fields_implicit(gid, knobs)
    use fields_implicit, only: aminv, nfield, nidx
    use gs2_layouts, only: f_lo, jf_lo, ij_idx
    use dist_fn, only: M_class, N_class, i_class
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    integer, intent(in) :: gid
    type(export_knobs_type), intent(in) :: knobs
    real(real64), dimension(:), allocatable :: buf
    integer, dimension(:,:,:), allocatable :: cls_ik, cls_it
    integer, dimension(:,:), allocatable :: sc_class, sc_member, sc_cell
    integer :: c, m, n, f, ig, ik, it, jflo, row, rr, nth, nn, mm, dc
    integer(kind=8) :: k
    character(len=16) :: cs
    call ex_scalar(gid, 'nfield', nfield, fname='fields_implicit::nfield')
    call ex_scalar(gid, 'nidx', nidx, fname='fields_implicit::nidx')
    if (.not. allocated(M_class)) return
    nth = 2*ntgrid + 1
    ! Class member tables (padded with -1) and their inverse
    allocate (cls_ik(maxval(M_class), maxval(N_class), i_class), cls_it(maxval(M_class), maxval(N_class), i_class))
    allocate (sc_class(ntheta0, naky), sc_member(ntheta0, naky), sc_cell(ntheta0, naky))
    cls_ik = -1 ; cls_it = -1 ; sc_class = -1 ; sc_member = -1 ; sc_cell = -1
    if (allocated(f_lo)) then
       do c = 1, i_class
          do n = 1, N_class(c)
             do m = 1, M_class(c)
                cls_ik(m, n, c) = f_lo(c)%ik(m, n)
                cls_it(m, n, c) = f_lo(c)%it(m, n)
                sc_class(f_lo(c)%it(m, n), f_lo(c)%ik(m, n)) = c
                sc_member(f_lo(c)%it(m, n), f_lo(c)%ik(m, n)) = m
                sc_cell(f_lo(c)%it(m, n), f_lo(c)%ik(m, n)) = n
             end do
          end do
       end do
    end if
    call ex_vec(gid, 'cls_ik', cls_ik, 'sc_Mmax sc_Nmax n_class', fname='gs2_layouts::f_lo(c)%ik(m,n)', &
         src='fields_implicit.f90:1074-1133 finish_fields_layouts', ibase=1)
    call ex_vec(gid, 'cls_it', cls_it, 'sc_Mmax sc_Nmax n_class', fname='gs2_layouts::f_lo(c)%it(m,n)', &
         src='fields_implicit.f90:1074-1133 finish_fields_layouts', ibase=1)
    call ex_vec(gid, 'sc_class', sc_class, 'kx ky', fname='class c of (it,ik)', ibase=1)
    call ex_vec(gid, 'sc_member', sc_member, 'kx ky', fname='member m of (it,ik)', ibase=1)
    call ex_vec(gid, 'sc_cell', sc_cell, 'kx ky', fname='cell n of (it,ik)', ibase=1)
    if (knobs%trace_response_build) call ex_att(gid, '', 'am_note', &
         'P-R1 (pre-inversion am) is not available: any probe inside init_response_matrix changes the '// &
         'code generation of the inlined Gauss-Jordan inversion (see gpu/spec/09a_exporter_notes.md)')
    call ex_absent(gid, 'am_c*')
    if (.not. knobs%write_response) then
       call ex_absent(gid, 'aminv_c*')
       return
    end if
    if (.not. allocated(aminv) .or. .not. allocated(f_lo)) then
       call ex_absent(gid, 'aminv_c*')
       return
    end if
    do c = 1, i_class
       write (cs, '(I0)') c
       rr = nidx*N_class(c)
       mm = M_class(c)
       if (.not. allocated(aminv(c)%dcell)) then
          call ex_absent(gid, 'aminv_c'//trim(cs))
          cycle
       end if
       allocate (buf(2_8*rr*rr*mm))
       buf = 0.0
       do m = 1, mm
          ik = f_lo(c)%ik(m, 1)
          do n = 1, N_class(c)
             it = f_lo(c)%it(m, n)
             do f = 1, nfield
                do ig = -ntgrid, ntgrid
                   jflo = ij_idx(jf_lo, ig, f, ik, it)
                   row = 1 + (ig+ntgrid) + nth*(f-1) + nidx*(n-1)
                   dc = jf_lo%dj(c, jflo)
                   if (dc < 1) call ex_fail('missing dcell in aminv export')
                   ! Fortran (ri, col, row, m)
                   do nn = 1, rr
                      k = 2_8*((nn-1) + int(rr, 8)*((row-1) + int(rr, 8)*(m-1)))
                      buf(k+1) = real(aminv(c)%dcell(dc)%supercell(nn), real64)
                      buf(k+2) = aimag(aminv(c)%dcell(dc)%supercell(nn))
                   end do
                end do
             end do
          end do
       end do
       call ex_vec(gid, 'aminv_c'//trim(cs), buf, 'ri sc_R_c'//trim(cs)//' sc_R_c'//trim(cs)//' sc_M_c'//trim(cs), &
            fname='fields_implicit::aminv(c)%dcell(jf_lo%dj(c,jflo))%supercell(:)', &
            src='fields_implicit.f90:722-1071 init_inverse_matrix', owner='04', dtdep=.true., tag='aminv')
       call ex_att(gid, 'aminv_c'//trim(cs), 'fortran_type', 'complex')
       call ex_att(gid, 'aminv_c'//trim(cs), 'layout', &
            'Fortran (ri, col, row, m): row = output (ig+ntgrid)+nth*(ifield-1)+nidx*(n-1)+1 (B order); '// &
            'col = input in the same B order; m = class member')
       deallocate (buf)
    end do
  end subroutine export_fields_implicit

  !> /fields_local: per-supercell response matrices via fields::dump_response_to_file (09 §7.5)
  subroutine export_fields_local(ncid, gid, knobs)
    use fields, only: dump_response_to_file
    use fields_local, only: fieldmat
    use fields_arrays, only: response_file, get_specific_response_file_name
    use gs2_save, only: gs2_restore_response
    use gs2_time, only: code_dt
    use file_utils, only: run_name
    use mp, only: proc0, barrier
    implicit none
    integer, intent(in) :: ncid, gid
    type(export_knobs_type), intent(in) :: knobs
    character(len=len(response_file)) :: saved
    character(len=256) :: fname
    character(len=64) :: vname, dname
    complex, dimension(:,:), allocatable :: resp
    real :: dt_file, condition
    integer :: ik, is, n, unit
    logical :: exists
    if (.not. knobs%write_response) then
       call ex_absent(gid, 'resp_*')
       return
    end if
    saved = response_file
    response_file = trim(run_name)//'.gs2export_tmp'
    call dump_response_to_file(suffix='.resp')
    call barrier
    do ik = 1, fieldmat%naky
       do is = 1, fieldmat%kyb(ik)%nsupercell
          write (vname, '("resp_ik",I0,"_is",I0)') ik, is
          fname = get_specific_response_file_name(ik, is, code_dt, '.resp')
          if (.not. proc0) cycle
          inquire (file=trim(fname), exist=exists)
          if (.not. exists) then
             call ex_absent(gid, trim(vname))
             cycle
          end if
          n = fieldmat%kyb(ik)%supercells(is)%nrow
          allocate (resp(n, n))
          call gs2_restore_response(resp, trim(fname), dt_file, condition)
          write (dname, '("lr_ax_ik",I0,"_is",I0)') ik, is
          if (.not. ex_has_dim(ncid, dname)) call ex_add_dim(ncid, trim(dname), n)
          call ex_vec(gid, trim(vname), resp, 'ri '//trim(dname)//' '//trim(dname), &
               fname='fields_local::fieldmat%kyb(ik)%supercells(is) (pull_rows_to_arr)', &
               src='fields_local.f90:987-1017 sc_dump_to_file round trip', owner='04', dtdep=.true.)
          call ex_att(gid, trim(vname), 'ik', ik)
          call ex_att(gid, trim(vname), 'is', is)
          call ex_att(gid, trim(vname), 'it_leftmost', fieldmat%kyb(ik)%supercells(is)%it_leftmost)
          call ex_att(gid, trim(vname), 'condition', condition)
          call ex_att(gid, trim(vname), 'dt', dt_file)
          call ex_att(gid, trim(vname), 'layout', 'Fortran (ri, ax1, ax2) = response(row, col) as written by '// &
               'gs2_save_response; NumPy [ax2, ax1, ri]')
          write (vname, '("it_leftmost_ik",I0,"_is",I0)') ik, is
          call ex_scalar(gid, trim(vname), fieldmat%kyb(ik)%supercells(is)%it_leftmost)
          deallocate (resp)
          open (newunit=unit, file=trim(fname), status='old')
          close (unit, status='delete')
       end do
    end do
    response_file = saved
  end subroutine export_fields_local

  !> /collisions (09 §3.10, 10 §A.11); dt-dependent part only if .not. full
  subroutine export_collisions(ncid, prefix, knobs, full)
    use collisions, only: vnew, vnew_s, vnew_D, vnew_E, delvnew, vnewh, pitch_weights, vpdiff, vpdiffle
    use collisions, only: c1le, betaale, qle, ec1le, ebetaale, eqle
    use collisions, only: z0le, s0le, w0le, bz0le, bs0le, bw0le, aj0le, vperp_aj1le, vpa_aj0_le
    use collisions, only: s0, w0, z0, bs0, bw0, bz0, vnmult
    use collisions, only: nxi_lim, cfac, drag, has_lorentz, has_diffuse, const_v, conserve_moments
    use collisions, only: conservative, resistivity, lorentz_switch, ediff_switch, ei_coll_only
    use collisions, only: special_wfb_lorentz, vpar_zero_mean, conserve_forbid_zero, force_collisions
    use collisions, only: hypermult, vnfac, vnslow, vnm_init, timesteps_between_collisions
    use collisions, only: collision_model_switch, colls, use_le_layout, adjust, heating, hyper_colls
    use collisions, only: split_collisions, vary_vnew, ncheck
    use le_grids, only: energy_maxw, speed_maxw, negrid
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prefix
    type(export_knobs_type), intent(in) :: knobs
    logical, intent(in) :: full
    integer :: gid, ie, j
    real, dimension(:), allocatable :: e15, capgh
    gid = ex_group(ncid, trim(prefix)//'collisions')
    call ex_scalar(gid, 'collision_model_switch', collision_model_switch)
    call ex_scalar(gid, 'colls', colls)
    ! always written: collisions_advance re-solves the fields from h whenever split_collisions,
    ! even if the collision operator itself is inactive (dist_fn.fpp collisions_advance)
    call ex_scalar(gid, 'split_collisions', split_collisions)
    call ex_vec(gid, 'vnmult', vnmult, 'two', fname='collisions::vnmult', src='collisions.fpp:1723-1733 (post-init)')
    if (.not. collisions_active() .or. .not. knobs%write_collisions) then
       call ex_att(gid, '', 'note', 'collisions inactive or write_collisions=.false.: arrays not exported')
       return
    end if
    ! dt-dependent matrices and conservation vectors (LE6, 06 §2.4)
    call wr_le6r(gid, 'c1le', c1le, 'xi_lim', 'egrid', 'collisions::c1le', 'collisions.fpp:2240-2285 init_lorentz', '06', .true.)
    call wr_le6r(gid, 'betaale', betaale, 'xi_lim', 'egrid', 'collisions::betaale', 'collisions.fpp:2240-2285', '06', .true.)
    call wr_le6r(gid, 'qle', qle, 'xi_lim', 'egrid', 'collisions::qle', 'collisions.fpp:2240-2285', '06', .true.)
    call wr_le6r(gid, 'ec1le', ec1le, 'xi_lim', 'egrid', 'collisions::ec1le', 'collisions.fpp:1981-1999 init_ediffuse', '06', .true.)
    call wr_le6r(gid, 'ebetaale', ebetaale, 'xi_lim', 'egrid', 'collisions::ebetaale', 'collisions.fpp:1981-1999', '06', .true.)
    call wr_le6r(gid, 'eqle', eqle, 'xi_lim', 'egrid', 'collisions::eqle', 'collisions.fpp:1981-1999', '06', .true.)
    call wr_le6r(gid, 'z0le', z0le, 'xi_lim', 'egrid_p1', 'collisions::z0le', 'collisions.fpp:691-1194', '06', .true.)
    call wr_le6r(gid, 's0le', s0le, 'xi_lim', 'egrid_p1', 'collisions::s0le', 'collisions.fpp:691-1194', '06', .true.)
    call wr_le6r(gid, 'w0le', w0le, 'xi_lim', 'egrid_p1', 'collisions::w0le', 'collisions.fpp:691-1194', '06', .true.)
    call wr_le6r(gid, 'bz0le', bz0le, 'xi_lim', 'egrid_p1', 'collisions::bz0le', 'collisions.fpp:1198-1679', '06', .true.)
    call wr_le6r(gid, 'bs0le', bs0le, 'xi_lim', 'egrid_p1', 'collisions::bs0le', 'collisions.fpp:1198-1679', '06', .true.)
    call wr_le6r(gid, 'bw0le', bw0le, 'xi_lim', 'egrid_p1', 'collisions::bw0le', 'collisions.fpp:1198-1679', '06', .true.)
    ! standard layout (use_le_layout = .false.; gated, G-24): g-like vectors only
    call wr_g7r(gid, 's0', s0, 'collisions::s0', 'collisions.fpp:1172-1192', '06', .true.)
    call wr_g7r(gid, 'w0', w0, 'collisions::w0', 'collisions.fpp:1172-1192', '06', .true.)
    call wr_g7r(gid, 'z0', z0, 'collisions::z0', 'collisions.fpp:1172-1192', '06', .true.)
    call wr_g7r(gid, 'bs0', bs0, 'collisions::bs0', 'collisions.fpp:1657-1677', '06', .true.)
    call wr_g7r(gid, 'bw0', bw0, 'collisions::bw0', 'collisions.fpp:1657-1677', '06', .true.)
    call wr_g7r(gid, 'bz0', bz0, 'collisions::bz0', 'collisions.fpp:1657-1677', '06', .true.)
    if (.not. use_le_layout) call ex_absent(gid, 'c1 betaa ql ec1 ebetaa eql')
    if (.not. full) return
    call wr_r3(gid, 'vnew', vnew, 'ky egrid species', 'collisions::vnew', 'collisions.fpp:1707-1876 init_vnew')
    call wr_r3(gid, 'vnew_s', vnew_s, 'ky egrid species', 'collisions::vnew_s', 'collisions.fpp:1707-1876')
    call wr_r3(gid, 'vnew_D', vnew_D, 'ky egrid species', 'collisions::vnew_D', 'collisions.fpp:1707-1876')
    call wr_r3(gid, 'vnew_E', vnew_E, 'ky egrid species', 'collisions::vnew_E', 'collisions.fpp:1707-1876')
    call wr_r3(gid, 'delvnew', delvnew, 'ky egrid species', 'collisions::delvnew', 'collisions.fpp:1707-1876')
    call wr_r4(gid, 'vnewh', vnewh, 'theta kx ky species', 'collisions::vnewh', 'collisions.fpp:1853-1863 (hyper_colls)')
    call wr_r2(gid, 'pitch_weights', pitch_weights, 'lambda theta', 'collisions::pitch_weights', &
         'collisions.fpp:2199-2210')
    call wr_r3(gid, 'vpdiff', vpdiff, 'theta sign lambda', 'collisions::vpdiff', 'collisions.fpp:4088-4140 init_vpdiff')
    call wr_r2(gid, 'vpdiffle', vpdiffle, 'xi_lim theta', 'collisions::vpdiffle', 'collisions.fpp:4142-4152')
    call wr_le6r(gid, 'aj0le', aj0le, 'xi_lim', 'egrid_p1', 'collisions::aj0le', 'collisions.fpp:634-687 init_le_bessel', '06')
    call wr_le6r(gid, 'vperp_aj1le', vperp_aj1le, 'xi_lim', 'egrid_p1', 'collisions::vperp_aj1le', &
         'collisions.fpp:634-687 init_le_bessel', '06')
    call wr_le6r(gid, 'vpa_aj0_le', vpa_aj0_le, 'xi_lim', 'egrid_p1', 'collisions::vpa_aj0_le', &
         'collisions.fpp:634-687 init_le_bessel', '06')
    ! Helper tables (06 §3.4): same expressions as GS2
    if (allocated(energy_maxw)) then
       allocate (e15(negrid))
       do ie = 1, negrid
          e15(ie) = energy_maxw(ie)**1.5
       end do
       call ex_vec(gid, 'e15', e15, 'egrid', fname='energy_maxw(ie)**1.5', src='collisions.fpp:3075 (exporter)', owner='06')
    end if
    if (allocated(speed_maxw) .and. negrid > 1) then
       allocate (capgh(negrid-1))
       do j = 1, negrid-1
          capgh(j) = capg((speed_maxw(j)+speed_maxw(j+1))*0.5)
       end do
       call ex_vec(gid, 'capg_half', capgh, 'egrid_m1', fname='capg((speed_maxw(j)+speed_maxw(j+1))*0.5)', &
            src='collisions.fpp:2165-2176 (verbatim copy)', owner='06')
    end if
    call ex_scalar(gid, 'nxi_lim', nxi_lim)
    call ex_scalar(gid, 'cfac', cfac)
    call ex_scalar(gid, 'drag', drag)
    call ex_scalar(gid, 'has_lorentz', has_lorentz)
    call ex_scalar(gid, 'has_diffuse', has_diffuse)
    call ex_scalar(gid, 'const_v', const_v)
    call ex_scalar(gid, 'conserve_moments', conserve_moments)
    call ex_scalar(gid, 'conservative', conservative)
    call ex_scalar(gid, 'resistivity', resistivity)
    call ex_scalar(gid, 'lorentz_switch', lorentz_switch)
    call ex_scalar(gid, 'ediff_switch', ediff_switch)
    call ex_scalar(gid, 'ei_coll_only', ei_coll_only)
    call ex_scalar(gid, 'special_wfb_lorentz', special_wfb_lorentz)
    call ex_scalar(gid, 'vpar_zero_mean', vpar_zero_mean)
    call ex_scalar(gid, 'conserve_forbid_zero', conserve_forbid_zero)
    call ex_scalar(gid, 'force_collisions', force_collisions)
    call ex_scalar(gid, 'hypermult', hypermult)
    call ex_scalar(gid, 'vnfac', vnfac)
    call ex_scalar(gid, 'vnslow', vnslow)
    call ex_vec(gid, 'vnm_init', vnm_init, 'two', fname='collisions::vnm_init')
    call ex_scalar(gid, 'timesteps_between_collisions', timesteps_between_collisions)
    call ex_scalar(gid, 'use_le_layout', use_le_layout)
    call ex_scalar(gid, 'adjust', adjust)
    call ex_scalar(gid, 'heating', heating)
    call ex_scalar(gid, 'hyper_colls', hyper_colls)
    call ex_scalar(gid, 'vary_vnew', vary_vnew)
    call ex_scalar(gid, 'ncheck', ncheck)
  contains
    ! Verbatim copies of the internal procedures of get_ediffuse_matrix (collisions.fpp:2165-2176)
    elemental real function capg_kernel(xe)
      use constants, only: sqrt_pi
      real, intent(in) :: xe
      ! This is the same as hsg_func(xe) * 2 * xe**2
      capg_kernel =  erf(xe) - 2 * xe * exp(-xe**2) / sqrt_pi
    end function capg_kernel

    elemental real function capg(xe)
      use constants, only: sqrt_pi
      real, intent(in) :: xe
      capg = 2 * exp(-xe**2) * capg_kernel(xe) / (xe * sqrt_pi)
    end function capg
  end subroutine export_collisions

  !> /hyper (09 §3.11, 10 §A.12); D_res and the simple3D filter are frozen at the first-init dt
  subroutine export_hyper(ncid, prefix, full)
    use hyper, only: D_res, hypervisc_filter, D_v, D_eta, nexp
    use hyper, only: D_hypervisc, D_hyperres, omega_osc, D_hyper
    use hyper, only: akx4_max, aky4_max, aky_max, akperp4_max, kperp2_max, kperp2_max_zonal
    use hyper, only: hyper_option_switch, const_amp, include_kpar, isotropic_shear, damp_zonal_only
    use hyper, only: hyper_on, gridnorm, D_hyper3D, P_hyper3D, ky_cut, kx_cut, isotropic_model
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prefix
    logical, intent(in) :: full
    integer, parameter :: hyper_option_simple3D = 5
    integer :: gid
    gid = ex_group(ncid, trim(prefix)//'hyper')
    call wr_r2(gid, 'D_res', D_res, 'kx ky', 'hyper::D_res', 'hyper.f90:400-411 init_hyper (frozen at first-init dt)', &
         '07', .true.)
    if (hyper_option_switch == hyper_option_simple3D) then
       call wr_r3(gid, 'hypervisc_filter', hypervisc_filter, 'theta kx ky', 'hyper::hypervisc_filter', &
            'hyper.f90:419-445 (simple3D, frozen at first-init dt)', '07', .true.)
    else
       call ex_absent(gid, 'hypervisc_filter')
    end if
    if (.not. full) return
    call ex_scalar(gid, 'D_v', D_v)
    call ex_scalar(gid, 'D_eta', D_eta)
    call ex_scalar(gid, 'nexp', nexp)
    call ex_scalar(gid, 'D_hypervisc', D_hypervisc)
    call ex_scalar(gid, 'D_hyperres', D_hyperres)
    call ex_scalar(gid, 'omega_osc', omega_osc)
    call ex_scalar(gid, 'D_hyper', D_hyper)
    call ex_scalar(gid, 'akx4_max', akx4_max)
    call ex_scalar(gid, 'aky4_max', aky4_max)
    call ex_scalar(gid, 'aky_max', aky_max)
    call ex_scalar(gid, 'akperp4_max', akperp4_max)
    call ex_scalar(gid, 'kperp2_max', kperp2_max)
    call ex_scalar(gid, 'kperp2_max_zonal', kperp2_max_zonal)
    call ex_scalar(gid, 'hyper_option_switch', hyper_option_switch)
    call ex_scalar(gid, 'const_amp', const_amp)
    call ex_scalar(gid, 'include_kpar', include_kpar)
    call ex_scalar(gid, 'isotropic_shear', isotropic_shear)
    call ex_scalar(gid, 'damp_zonal_only', damp_zonal_only)
    call ex_scalar(gid, 'hyper_on', hyper_on)
    call ex_scalar(gid, 'gridnorm', gridnorm)
    call ex_scalar(gid, 'D_hyper3D', D_hyper3D)
    call ex_scalar(gid, 'P_hyper3D', P_hyper3D)
    call ex_scalar(gid, 'ky_cut', ky_cut)
    call ex_scalar(gid, 'kx_cut', kx_cut)
    call ex_scalar(gid, 'isotropic_model', isotropic_model)
  end subroutine export_hyper

  !> /nonlinear_terms (09 §3.12)
  subroutine export_nonlinear_terms(ncid, prefix, full)
    use nonlinear_terms, only: nonlin, split_nonlinear, accelerated, cfl, cflx, cfly, max_vel_components
    use nonlinear_terms, only: nonlinear_mode_switch, use_cfl_limit, use_2d_cfl, error_target
    use nonlinear_terms, only: istep_error_start, use_order_based_error, include_apar, include_bpar
    use nonlinear_terms, only: include_phi, zip, nl_forbid_force_zero
    use nonlinear_terms, only: nl_istep_last => istep_last, nl_initializing => initializing
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prefix
    logical, intent(in) :: full
    integer :: gid
    gid = ex_group(ncid, trim(prefix)//'nonlinear_terms')
    call ex_scalar(gid, 'cflx', cflx, src='nonlinear_terms.fpp:357-358')
    call ex_scalar(gid, 'cfly', cfly, src='nonlinear_terms.fpp:357-358')
    if (.not. full) return
    call ex_scalar(gid, 'nonlin', nonlin)
    call ex_scalar(gid, 'split_nonlinear', split_nonlinear)
    call ex_scalar(gid, 'accelerated', accelerated)
    call ex_scalar(gid, 'cfl', cfl)
    call ex_vec(gid, 'max_vel_components', max_vel_components, 'three', fname='nonlinear_terms::max_vel_components')
    call ex_scalar(gid, 'nonlinear_mode_switch', nonlinear_mode_switch)
    call ex_scalar(gid, 'use_cfl_limit', use_cfl_limit)
    call ex_scalar(gid, 'use_2d_cfl', use_2d_cfl)
    call ex_scalar(gid, 'error_target', error_target)
    call ex_scalar(gid, 'istep_error_start', istep_error_start)
    call ex_scalar(gid, 'use_order_based_error', use_order_based_error)
    call ex_scalar(gid, 'include_phi', include_phi)
    call ex_scalar(gid, 'include_apar', include_apar)
    call ex_scalar(gid, 'include_bpar', include_bpar)
    call ex_scalar(gid, 'zip', zip)
    call ex_scalar(gid, 'nl_forbid_force_zero', nl_forbid_force_zero)
    call ex_scalar(gid, 'nl_istep_last', nl_istep_last, fname='nonlinear_terms::istep_last')
    call ex_scalar(gid, 'nl_initializing', nl_initializing, fname='nonlinear_terms::initializing')
  end subroutine export_nonlinear_terms

  !============================================================== state

  !> State group (09 §7.6): time block + g, gnew, fields, gexp_*, shifts
  subroutine export_state(ncid, gname, istep, knobs)
    use dist_fn_arrays, only: g, gnew, gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use run_parameters, only: has_phi, has_apar, has_bpar
    use gs2_layouts, only: g_lo
    implicit none
    integer, intent(in) :: ncid, istep
    character(len=*), intent(in) :: gname
    type(export_knobs_type), intent(in) :: knobs
    integer :: gid
    gid = ex_group(ncid, gname)
    call write_time_block(gid, istep)
    call wr_g7c(gid, 'g', g, 'dist_fn_arrays::g')
    call wr_g7c(gid, 'gnew', gnew, 'dist_fn_arrays::gnew')
    call put_field(gid, 'phi', phi, has_phi)
    call put_field(gid, 'apar', apar, has_apar)
    call put_field(gid, 'bpar', bpar, has_bpar)
    call put_field(gid, 'phinew', phinew, has_phi)
    call put_field(gid, 'aparnew', aparnew, has_apar)
    call put_field(gid, 'bparnew', bparnew, has_bpar)
    call wr_g7c(gid, 'gexp_1', gexp_1, 'dist_fn_arrays::gexp_1')
    call wr_g7c(gid, 'gexp_2', gexp_2, 'dist_fn_arrays::gexp_2')
    call wr_g7c(gid, 'gexp_3', gexp_3, 'dist_fn_arrays::gexp_3')
    call wr_r1(gid, 'kx_shift', kx_shift, 'ky', 'dist_fn_arrays::kx_shift')
    call wr_r1(gid, 'theta0_shift', theta0_shift, 'ky', 'dist_fn_arrays::theta0_shift')
    if (knobs%write_raw_glo .and. allocated(gnew)) then
       call ex_vec(gid, 'gnew_glo', gnew(:, :, g_lo%llim_proc:g_lo%ulim_proc), 'ri theta sign glo', &
            fname='dist_fn_arrays::gnew (glo order)', tag='raw')
    end if
  contains
    subroutine put_field(gid, name, arr, present_flag)
      integer, intent(in) :: gid
      character(len=*), intent(in) :: name
      complex, dimension(:,:,:), allocatable, intent(in) :: arr
      logical, intent(in) :: present_flag
      if (allocated(arr)) then
         call ex_vec(gid, name, arr, 'ri theta kx ky', fname='fields_arrays::'//name, tag='F3')
         call ex_att(gid, name, 'present', merge(1, 0, present_flag))
      else
         call ex_absent(gid, name)
      end if
    end subroutine put_field
  end subroutine export_state

  !============================================================== dt-dependent set

  !> The dt-dependent arrays of one dt-ladder rung (09 §6.3) under group `prefix`
  subroutine export_dt_dependent(ncid, prefix, knobs)
    use gs2_time, only: code_dt
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prefix
    type(export_knobs_type), intent(in) :: knobs
    integer :: gid
    gid = ex_group(ncid, trim(prefix))
    call ex_att(gid, '', 'code_dt', code_dt)
    call write_time_block(gid, -1)
    gid = ex_group(ncid, trim(prefix)//'/dist_fn')
    call export_dist_fn_coeffs(gid)
    gid = ex_group(ncid, trim(prefix)//'/dist_fn_arrays')
    call export_wstar(gid)
    call export_collisions(ncid, trim(prefix)//'/', knobs, .false.)
    call export_field_matrices(ncid, trim(prefix)//'/', knobs)
    call export_hyper(ncid, trim(prefix)//'/', .false.)
    call export_nonlinear_terms(ncid, trim(prefix)//'/', .false.)
    call export_gs2_time(ncid, trim(prefix)//'/gs2_time')
  end subroutine export_dt_dependent

  !============================================================== diagnostics

  !> /diagnostics: parsed gnostics values after init_gs2_diagnostics_new (08 §2.3)
  subroutine export_diagnostics_group(ncid)
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif
    implicit none
    integer, intent(in) :: ncid
    integer :: gid
    gid = ex_group(ncid, 'diagnostics')
#ifdef NEW_DIAG
    call ex_scalar(gid, 'nout', gnostics%nout)
    call ex_scalar(gid, 'istep', gnostics%istep)
    call ex_scalar(gid, 'nwrite', gnostics%nwrite)
    call ex_scalar(gid, 'nwrite_mult', gnostics%nwrite_mult)
    call ex_scalar(gid, 'navg', gnostics%navg)
    call ex_scalar(gid, 'igomega', gnostics%igomega)
    call ex_scalar(gid, 'omegatinst', gnostics%omegatinst)
    call ex_scalar(gid, 'omegatol', gnostics%omegatol)
    call ex_scalar(gid, 'exit_when_converged', gnostics%exit_when_converged)
    call ex_scalar(gid, 'calculate_fluxes', gnostics%calculate_fluxes)
    call ex_scalar(gid, 'user_time_old', gnostics%user_time_old)
    call ex_scalar(gid, 'start_time', gnostics%start_time)
    call ex_scalar(gid, 'ncheck', gnostics%ncheck)
    call ex_scalar(gid, 'nsave', gnostics%nsave)
    call ex_scalar(gid, 'nc_sync_freq', gnostics%nc_sync_freq)
    call ex_scalar(gid, 'nmovie', gnostics%nmovie)
    call ex_scalar(gid, 'conv_nstep_av', gnostics%conv_nstep_av)
    call ex_scalar(gid, 'conv_test_multiplier', gnostics%conv_test_multiplier)
    call ex_scalar(gid, 'conv_min_step', gnostics%conv_min_step)
    call ex_scalar(gid, 'conv_max_step', gnostics%conv_max_step)
    call ex_scalar(gid, 'conv_nsteps_converged', gnostics%conv_nsteps_converged)
    call ex_scalar(gid, 'write_any', gnostics%write_any)
    call ex_scalar(gid, 'append_old', gnostics%append_old)
    call ex_scalar(gid, 'enable_parallel', gnostics%enable_parallel)
    call ex_scalar(gid, 'serial_netcdf4', gnostics%serial_netcdf4)
    call ex_scalar(gid, 'print_line', gnostics%print_line)
    call ex_scalar(gid, 'print_flux_line', gnostics%print_flux_line)
    call ex_scalar(gid, 'write_line', gnostics%write_line)
    call ex_scalar(gid, 'write_flux_line', gnostics%write_flux_line)
    call ex_scalar(gid, 'write_fields', gnostics%write_fields)
    call ex_scalar(gid, 'write_phi_over_time', gnostics%write_phi_over_time)
    call ex_scalar(gid, 'write_apar_over_time', gnostics%write_apar_over_time)
    call ex_scalar(gid, 'write_bpar_over_time', gnostics%write_bpar_over_time)
    call ex_scalar(gid, 'make_movie', gnostics%make_movie)
    call ex_scalar(gid, 'dump_fields_periodically', gnostics%dump_fields_periodically)
    call ex_scalar(gid, 'write_moments', gnostics%write_moments)
    call ex_scalar(gid, 'write_full_moments_notgc', gnostics%write_full_moments_notgc)
    call ex_scalar(gid, 'write_ntot_over_time', gnostics%write_ntot_over_time)
    call ex_scalar(gid, 'write_density_over_time', gnostics%write_density_over_time)
    call ex_scalar(gid, 'write_upar_over_time', gnostics%write_upar_over_time)
    call ex_scalar(gid, 'write_tperp_over_time', gnostics%write_tperp_over_time)
    call ex_scalar(gid, 'write_fluxes', gnostics%write_fluxes)
    call ex_scalar(gid, 'write_fluxes_by_mode', gnostics%write_fluxes_by_mode)
    call ex_scalar(gid, 'write_symmetry', gnostics%write_symmetry)
    call ex_scalar(gid, 'write_nl_flux_dist', gnostics%write_nl_flux_dist)
    call ex_scalar(gid, 'write_parity', gnostics%write_parity)
    call ex_scalar(gid, 'write_omega', gnostics%write_omega)
    call ex_scalar(gid, 'write_ql_metric', gnostics%write_ql_metric)
    call ex_scalar(gid, 'write_verr', gnostics%write_verr)
    call ex_scalar(gid, 'write_cerr', gnostics%write_cerr)
    call ex_scalar(gid, 'write_max_verr', gnostics%write_max_verr)
    call ex_scalar(gid, 'write_heating', gnostics%write_heating)
    call ex_scalar(gid, 'write_ascii', gnostics%write_ascii)
    call ex_scalar(gid, 'write_gyx', gnostics%write_gyx)
    call ex_scalar(gid, 'write_g', gnostics%write_g)
    call ex_scalar(gid, 'use_nonlin_convergence', gnostics%use_nonlin_convergence)
    call ex_scalar(gid, 'write_cross_phase', gnostics%write_cross_phase)
    call ex_scalar(gid, 'write_correlation', gnostics%write_correlation)
    call ex_scalar(gid, 'write_correlation_extend', gnostics%write_correlation_extend)
    call ex_scalar(gid, 'write_jext', gnostics%write_jext)
    call ex_scalar(gid, 'write_lorentzian', gnostics%write_lorentzian)
    call ex_scalar(gid, 'write_eigenfunc', gnostics%write_eigenfunc)
    call ex_scalar(gid, 'write_final_fields', gnostics%write_final_fields)
    call ex_scalar(gid, 'write_kpar', gnostics%write_kpar)
    call ex_scalar(gid, 'write_final_epar', gnostics%write_final_epar)
    call ex_scalar(gid, 'write_final_db', gnostics%write_final_db)
    call ex_scalar(gid, 'write_final_moments', gnostics%write_final_moments)
    call ex_scalar(gid, 'write_final_antot', gnostics%write_final_antot)
    call ex_scalar(gid, 'write_gs', gnostics%write_gs)
    call ex_scalar(gid, 'save_for_restart', gnostics%save_for_restart)
    call ex_scalar(gid, 'save_many', gnostics%save_many)
    call ex_scalar(gid, 'file_safety_check', gnostics%file_safety_check)
    call ex_scalar(gid, 'save_distfn', gnostics%save_distfn)
    call ex_scalar(gid, 'save_glo_info_and_grids', gnostics%save_glo_info_and_grids)
    call ex_scalar(gid, 'save_velocities', gnostics%save_velocities)
    call ex_scalar(gid, 'write_zonal_transfer', gnostics%write_zonal_transfer)
    call ex_scalar(gid, 'write_kinetic_energy_transfer', gnostics%write_kinetic_energy_transfer)
    call ex_scalar(gid, 'write_collisional', gnostics%write_collisional)
    call ex_scalar(gid, 'write_omavg', gnostics%write_omavg)
    call ex_scalar(gid, 'ob_midplane', gnostics%ob_midplane)
    call ex_scalar(gid, 'write_avg_moments', gnostics%write_avg_moments)
#endif
  end subroutine export_diagnostics_group

  !============================================================== gates

  !> /gates (09 §8.3, 10 §B, 00_decisions D7): one int8 0-D variable per gate (1 = violated);
  !> returns the gate_status string ("ok" or a ';'-separated list of violated soft gates).
  !> diag_ok: diagnostics were initialised (diagnostics gates evaluated from gnostics).
  function export_gates(ncid, diag_ok) result(status)
    use mp, only: nproc
    use antenna, only: no_driver
    use run_parameters, only: ieqzip, eqzip_option_switch, wstar_units, use_old_diagnostics, do_eigsolve
    use dist_fn, only: afilter, gf_lo_integrate, source_option_switch, phi_ext, opt_source
    use dist_fn, only: boundary_option_switch, boundary_option_linked, g_exb, g_exbfac
    use dist_fn, only: boundary_off_grid, exponential_boundary, start_from_previous_solution
    use nonlinear_terms, only: split_nonlinear, nonlin, zip
    use fields, only: fieldopt_switch, fieldopt_test, fieldopt_gf_local
    use gs2_layouts, only: fft_measure_plan, gf_local_fields
    use init_g, only: new_field_init, get_init_g_config, init_g_config_type
    use gs2_reinit, only: in_memory
    use kt_grids, only: aky, ntheta0, gridopt_switch, gridopt_box
    use theta_grid, only: shat, bmag, bmax, ntgrid
    use le_grids, only: al, nlambda, bouncefuzz, grid_has_trapped_particles
    use hyper, only: hyper_on, hyper_option_switch, gridnorm, akperp4_max, isotropic_shear, akx4_max, aky_max
    use hyper, only: kperp2_max, kperp2_max_zonal
    use species, only: spec, has_hybrid_electron_species, f0_maxwellian
    use collisions, only: use_le_layout, has_lorentz, has_diffuse, conserve_moments, conservative
    use collisions, only: lorentz_switch, ediff_switch
    use build_config, only: gs2_has_lowflow, gs2_is_single_precision
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif
    implicit none
    integer, intent(in) :: ncid
    logical, intent(in) :: diag_ok
    character(len=:), allocatable :: status
    integer :: gid, ig, il
    type(init_g_config_type) :: icfg
    logical :: colls_on, flag, restart_type
    integer, parameter :: hyper_visc = 2, hyper_res = 3, hyper_both = 4, hyper_s3d = 5
    integer, parameter :: eqzip_none = 1, source_phiext_full = 5
    status = ''
    gid = ex_group(ncid, 'gates')
    colls_on = collisions_active()
    icfg = get_init_g_config()

    ! hard gates (the program aborts before export if violated; recorded for completeness)
    call gate('hard_single_precision', gs2_is_single_precision .or. kind(1.0) /= real64, 'single precision build')
    call gate('hard_lowflow', gs2_has_lowflow, 'LOWFLOW build (wdrift complex)')
    call gate('hard_multi_rank', nproc > 1, 'nproc > 1 (reference data must be single rank)')
    call gate('hard_do_eigsolve', do_eigsolve, 'do_eigsolve (not a time stepper)')
    ! soft gates (10 §B; 09 §8.3)
    call gate('antenna_on', .not. no_driver, 'antenna active (RNG stream, 10 B.2)')
    call gate('eqzip', eqzip_option_switch /= eqzip_none .or. anyl2(ieqzip), 'eqzip_option / ieqzip (G-07)')
    call gate('afilter', abs(afilter) > epsilon(0.0), 'afilter /= 0 (debug knob)')
    call gate('split_nonlinear', split_nonlinear, 'split_nonlinear (G-08/05 §4.13)')
    call gate('zip', zip, 'zip (G-08)')
    call gate('field_option_test_or_gf_local', fieldopt_switch == fieldopt_test .or. &
         fieldopt_switch == fieldopt_gf_local, "field_option 'test' or 'gf_local'")
    call gate('gf_lo_integrate', gf_lo_integrate, 'gf_lo_integrate')
    call gate('gf_local_fields', gf_local_fields, 'gf_local_fields')
    call gate('fft_measure_plan', fft_measure_plan, 'fft_measure_plan (non-reproducible FFTW plans)')
    call gate('new_field_init_false', .not. new_field_init, 'new_field_init = .false.')
    call gate('wstar_units', wstar_units, 'wstar_units (G-16, gated in v1)')
    call gate('use_old_diagnostics', use_old_diagnostics, 'use_old_diagnostics (G-09)')
    call gate('phiext_full_with_phi_ext', source_option_switch == source_phiext_full .and. phi_ext /= 0.0, &
         "source_option='phiext_full' with phi_ext /= 0 (G-14)")
    call gate('opt_source', opt_source, 'opt_source')
    call gate('linked_non_box', boundary_option_switch == boundary_option_linked .and. &
         gridopt_switch /= gridopt_box, "boundary_option='linked' without grid_option='box'")
    call gate('experimental_boundary', boundary_off_grid .or. exponential_boundary .or. &
         start_from_previous_solution, 'boundary_off_grid / exponential_boundary / start_from_previous_solution')
    call gate('nonlinear_non_box', nonlin .and. gridopt_switch /= gridopt_box, 'nonlinear_mode on without box grid')
    call gate('exb_single_kx', abs(g_exb*g_exbfac) > epsilon(0.) .and. ntheta0 == 1, 'g_exb /= 0 with ntheta0 == 1')
    ! hyperviscosity zero normalisers (07 App. B)
    flag = .false.
    if (hyper_on) then
       select case (hyper_option_switch)
       case (hyper_visc, hyper_res, hyper_both)
          if (gridnorm .and. akperp4_max == 0.0) flag = .true.
          if (.not. isotropic_shear .and. (akx4_max == 0.0 .or. aky_max == 0.0)) flag = .true.
       case (hyper_s3d)
          if (shat == 0.0 .or. kperp2_max == 0.0 .or. kperp2_max_zonal == 0.0) flag = .true.
       end select
    end if
    call gate('hyper_zero_normaliser', flag, 'hyperviscosity with a zero normaliser (07 App. B)')
    ! scale quirk (01 §6.12, 10 B.1)
    select case (trim(icfg%ginit_option))
    case ('restart', 'many', 'small', 'file', 'cont', 'smallflat', 'zonal_only', 'no_zonal', 'eig_restart')
       restart_type = .true.
    case default
       restart_type = .false.
    end select
    call gate('scale_not_one', icfg%scale /= 1.0 .and. .not. restart_type .and. .not. in_memory, &
         'init_g scale /= 1 with a non-restart ginit and file-based resets (01 §6.12)')
    ! collisions
    call gate('hybrid_with_collisions', colls_on .and. has_hybrid_electron_species(spec), &
         'hybrid electrons with collisions (G-25, Q-12)')
    call gate('nonmaxwellian_with_collisions', colls_on .and. any(spec(:)%f0type /= f0_maxwellian), &
         'non-Maxwellian f0 with collisions (Q-12)')
    call gate('use_le_layout_false', colls_on .and. .not. use_le_layout, 'use_le_layout = .false. (G-24)')
    call gate('no_trapped_lorentz_conserve', colls_on .and. has_lorentz .and. .not. grid_has_trapped_particles() &
         .and. (conserve_moments .or. .not. use_le_layout), 'no trapped particles with Lorentz conservation (06 §8.1)')
    call gate('conservative_false', colls_on .and. .not. conservative, 'conservative = .false.')
    call gate('old_collision_scheme', colls_on .and. ((has_lorentz .and. lorentz_switch /= 1) .or. &
         (has_diffuse .and. ediff_switch /= 1)), "lorentz_scheme/ediff_scheme = 'old'")
    ! 00_decisions D7
    flag = .false.
    if (grid_has_trapped_particles()) flag = (bmag(-ntgrid) /= bmax) .or. (bmag(ntgrid) /= bmax)
    call gate('Q06_bmax_not_at_ends', flag, 'bmax not at ig = +-ntgrid with trapped particles (D7/Q-06)')
    call gate('Q07_zero_shear_with_zonal', shat == 0.0 .and. any(aky == 0.0), 'shat == 0 with an aky == 0 mode (D7/Q-07)')
    flag = .false.
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          if (1.0 - al(il)*bmag(ig) == -bouncefuzz) flag = .true.
       end do
    end do
    call gate('Q09_bouncefuzz_edge', flag, '1 - al*bmag == -bouncefuzz at a grid point (D7/Q-09)')
    ! diagnostics (10 B.7)
#ifdef NEW_DIAG
    if (diag_ok) then
       call gate('diag_use_nonlin_convergence', gnostics%use_nonlin_convergence, 'use_nonlin_convergence (G-22)')
       call gate('diag_save_distfn', gnostics%save_distfn, 'save_distfn')
       call gate('diag_append_old', gnostics%append_old, 'append_old')
       call gate('diag_make_movie', gnostics%make_movie, 'make_movie')
       call gate('diag_dump_fields_periodically', gnostics%dump_fields_periodically, 'dump_fields_periodically')
       call gate('diag_write_heating', gnostics%write_heating, 'write_heating')
       call gate('diag_write_collisional', gnostics%write_collisional, 'write_collisional')
       call gate('diag_write_verr', gnostics%write_verr, 'write_verr (feeds vnmult via vary_vnew)')
       call gate('diag_write_max_verr', gnostics%write_max_verr, 'write_max_verr')
       call gate('diag_write_cerr', gnostics%write_cerr, 'write_cerr')
       call gate('diag_write_symmetry', gnostics%write_symmetry, 'write_symmetry')
       call gate('diag_write_nl_flux_dist', gnostics%write_nl_flux_dist, 'write_nl_flux_dist')
       call gate('diag_write_parity', gnostics%write_parity, 'write_parity')
       call gate('diag_write_cross_phase', gnostics%write_cross_phase, 'write_cross_phase')
       call gate('diag_write_correlation_extend', gnostics%write_correlation_extend, 'write_correlation_extend')
       call gate('diag_write_jext', gnostics%write_jext, 'write_jext')
       call gate('diag_write_lorentzian', gnostics%write_lorentzian, 'write_lorentzian')
       call gate('diag_write_zonal_transfer', gnostics%write_zonal_transfer, 'write_zonal_transfer')
       call gate('diag_write_kinetic_energy_transfer', gnostics%write_kinetic_energy_transfer, &
            'write_kinetic_energy_transfer')
       call gate('diag_write_full_moments_notgc', gnostics%write_full_moments_notgc, 'write_full_moments_notgc')
       call gate('diag_write_g', gnostics%write_g, 'write_g')
       call gate('diag_write_gyx', gnostics%write_gyx, 'write_gyx')
       call gate('diag_write_gs', gnostics%write_gs, 'write_gs')
       call gate('diag_write_kpar', gnostics%write_kpar, 'write_kpar')
       call gate('diag_write_final_db', gnostics%write_final_db, 'write_final_db')
    else
       call ex_att(gid, '', 'diagnostics_gates', 'not evaluated (diagnostics not initialised)')
    end if
#endif
    if (len(status) == 0) status = 'ok'
  contains
    subroutine gate(name, violated, description)
      character(len=*), intent(in) :: name, description
      logical, intent(in) :: violated
      call ex_scalar(gid, name, violated)
      call ex_att(gid, name, 'description', description)
      if (violated .and. name(1:min(5, len(name))) /= 'hard_') then
         if (len(status) == 0) then
            status = name
         else
            status = status//';'//name
         end if
      end if
    end subroutine gate
    logical function anyl2(x)
      logical, dimension(:,:), allocatable, intent(in) :: x
      anyl2 = .false.
      if (allocated(x)) anyl2 = any(x)
    end function anyl2
  end function export_gates
end module gs2_export_static
