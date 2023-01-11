!> FIXME : Add documentation
module gs2_save

  use mp, only: mp_comm, mp_info
  use constants, only: run_name_size

# ifdef NETCDF
# ifdef NETCDF_PARALLEL
  use netcdf, only: NF90_MPIIO
  use netcdf, only: nf90_var_par_access, NF90_COLLECTIVE
  use netcdf, only: nf90_put_att, NF90_GLOBAL, nf90_get_att
  use netcdf, only: nf90_create_par
# endif
  use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR, NF90_UNLIMITED, NF90_CHAR, NF90_NETCDF4
  use netcdf, only: nf90_create, nf90_open, nf90_sync, nf90_close, nf90_enddef
  use netcdf, only: nf90_put_var, nf90_get_var, nf90_strerror
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid, nf90_inquire_variable
  
  use netcdf_utils, only: get_netcdf_code_precision
  use netcdf_utils, only: check_netcdf_file_precision
  use netcdf_utils, only: netcdf_error
  use netcdf_utils, only: netcdf_real, netcdf_int, kind_nf
# endif

  implicit none

  private

  public :: gs2_restore, gs2_save_for_restart, finish_gs2_save

  public :: read_many, save_many, gs2_save_response, gs2_restore_response
  public :: restore_current_scan_parameter_value, proc_to_save_fields
  public :: init_gs2_save, init_dt, read_t0_from_restart_file, init_ant_amp
  public :: set_restart_file, include_explicit_source_in_restart
  public :: init_vnm, restart_writable, EigNetcdfID
  public :: init_eigenfunc_file, finish_eigenfunc_file, add_eigenpair_to_file
  public :: gs2_save_slice

  interface gs2_restore
     module procedure gs2_restore_many
  end interface

  !> A custom type to look after the netcdf ids for the eigenvalue file
  type EigNetcdfID
     integer :: ncid, conv_dim_id, theta_dim_id, ri_dim_id
     integer :: omega_var_id, theta_var_id, phi_var_id
     integer :: apar_var_id, bpar_var_id, conv_var_id
     integer :: nconv_count
  end type EigNetcdfID

  !> Read and write single or multiple restart files.
  !> Effectively ignored if not built with parallel netcdf
  logical :: read_many = .false., save_many = .false.

  !> Do we want to include the explicit source terms and
  !> related timesteps when saving/restoring from the restart file.
  logical :: include_explicit_source_in_restart = .true.

  !> Which processor should save the potentials. If <0 (default) then
  !> all processors save their potentials. If >= nproc then will be
  !> set to mod(proc_to_save_fields, nproc)
  integer :: proc_to_save_fields

  !> Used to hold the base of the restart file name.
  character (run_name_size) :: restart_file

#ifdef NETCDF
  real, allocatable, dimension(:,:,:) :: tmpr, tmpi, ftmpr, ftmpi
  real, allocatable, dimension(:) :: stmp, atmp
  ! Netcdf ids for quantities referenced in multiple routines.
  integer (kind_nf) :: ncid, thetaid, signid, gloid, kyid, kxid, nk_stir_dim
  integer (kind_nf) :: phir_id, phii_id, aparr_id, apari_id, bparr_id, bpari_id
  integer (kind_nf) :: t0id, gr_id, gi_id, vnm1id, vnm2id
  integer (kind_nf) :: delt0id, delt1id, delt2id, delt_max_id
  integer (kind_nf) :: gexp_1_r_id, gexp_1_i_id
  integer (kind_nf) :: gexp_2_r_id, gexp_2_i_id
  integer (kind_nf) :: gexp_3_r_id, gexp_3_i_id
  integer (kind_nf) :: layout_id, kx_shift_id, theta0_shift_id
  integer (kind_nf) :: current_scan_parameter_value_id
  integer (kind_nf) :: a_antr_id, b_antr_id, a_anti_id, b_anti_id
  !> Do we save the scan parameter or not?
  logical, parameter :: include_parameter_scan = .true.
  !> Open new netCDF files with netCDF-4 model, and overwrite existing files
  integer, parameter :: create_file_mode = ior(NF90_NETCDF4, NF90_CLOBBER)

  !> The NETCDF_PARALLEL directives include code for parallel
  !> netcdf using HDF5 to write the output to a single restart file
  !> The save_many/read_many flag allows the old style multiple file output
  !> Here we set a run time flag to enable us to handle some of the different
  !> paths through the code at run time rather than with ifdef mazes
#ifdef NETCDF_PARALLEL
  logical, parameter :: has_netcdf_parallel = .true.
#else
  logical, parameter :: has_netcdf_parallel = .false.
#endif
#endif
contains

  !> Create and fill the restart files
  subroutine gs2_save_for_restart &
       (g, t0, delt0, vnm, istatus, fphi, fapar, fbpar, fileopt, &
        save_glo_info_and_grids, save_velocities, header)
# ifdef NETCDF
    use fields_arrays, only: phinew, aparnew, bparnew
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use kt_grids, only: naky, ntheta0, jtwist_out, akx, aky
    use antenna_data, only: nk_stir, a_ant, b_ant, ant_on
    use netcdf_utils, only: ensure_netcdf_var_exists, ensure_netcdf_dim_exists
    use gs2_metadata, only: create_metadata
# else
    use mp, only: proc0
# endif    
    use mp, only: iproc, barrier, proc0
    use theta_grid, only: ntgrid, theta
    use gs2_layouts, only: g_lo, layout
    use gs2_layouts, only: proc_id,ik_idx,it_idx,is_idx,ie_idx,il_idx
    use layouts_type, only: g_layout_type
    use file_utils, only: error_unit
    use le_grids, only: energy, al, negrid, nlambda
    use species, only: nspec
    use dist_fn_arrays, only: vpa, vperp2
    use parameter_scan_arrays, only: current_scan_parameter_value
    use unit_tests, only: debug_message
    use gs2_time, only: code_dt, code_dt_prev1, code_dt_prev2, code_dt_max
    use optionals, only: get_option_with_default
    use standard_header, only : standard_header_type
    implicit none
    !> The distribution function \(g\)
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    !> Current simulation time
    real, intent (in) :: t0
    !> Current timestep size
    real, intent (in) :: delt0
    !> FIXME: Collisionality multiplier?
    real, dimension (2), intent (in) :: vnm
    !> Electrostatic potential, \(\phi\), scaling factor
    real, intent (in) :: fphi
    !> Vector potential, \(A_\parallel\), scaling factor
    real, intent (in) :: fapar
    !> Parallel magnetic field, \(B_\parallel\), scaling factor
    real, intent (in) :: fbpar
    !> Last netCDF status: zero on success, non-zero on failure
    integer, intent (out) :: istatus
    !> Optional extra filename infix
    character (len=*), intent (in), optional :: fileopt
    !> Control if layout and diagnostic information is saved
    logical, intent (in), optional :: save_glo_info_and_grids
    !> Control if \(v_\parallel, v_\perp\) information is saved
    logical, intent (in), optional :: save_velocities
    !> Header for files with build and run information
    type(standard_header_type), optional, intent(in) :: header
# ifdef NETCDF
    !> Arrays holding information about the layout for plotting
    integer, dimension (:,:,:,:,:), allocatable :: procid_kykxsel, iglo_kykxsel
    character (run_name_size) :: file_proc
    integer :: i, n_elements, ierr
    integer :: total_elements
    integer, dimension(3) :: start_pos, counts
    integer, parameter :: tmpunit = 348
    logical :: is_one_file_per_processor, this_proc_saves_fields
    logical :: has_gexp_1
    logical :: save_glo_info_and_grids_local, save_velocities_info
    logical :: should_define_single_entry_variables
    integer, parameter :: verb = 3
    integer :: iglo
    ! Netcdf ids for variables only handled in this routine
    integer (kind_nf) :: iglo_kykxsel_id, procid_kykxsel_id, iglo_llim_proc_id
    integer (kind_nf) :: theta_grid_id, kx_grid_id, ky_grid_id, jtwist_id
    integer (kind_nf) :: egridid,lgridid, vpa_id, vperp2_id
    integer (kind_nf) :: energy_id, lambda_id
    integer (kind_nf) :: nspecid, spec_id, char5_dim
    integer (kind_nf) :: ntheta0_id, naky_id, nlambda_id, negrid_id, nspec_id

    type(standard_header_type) :: local_header

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart initialising')

    istatus = 0
    ierr = error_unit()

    ! Handle some optional inputs
    save_glo_info_and_grids_local = get_option_with_default(save_glo_info_and_grids, .false.)
    save_velocities_info = get_option_with_default(save_velocities, .false.)

    ! Work out the problem size
    n_elements = max(g_lo%ulim_proc-g_lo%llim_proc+1, 0)
    total_elements = g_lo%ulim_world+1

    ! Decide if we want one file per processor or not
    is_one_file_per_processor = save_many .or. (.not. has_netcdf_parallel)

    ! Decide if we are saving the fields
    this_proc_saves_fields = (iproc == proc_to_save_fields) .or. (proc_to_save_fields < 0)

    ! Open the file and define dimensions and variables
    file_proc = get_file_proc(is_one_file_per_processor, fileopt)

    call debug_message(verb, 'gs2_save::gs2_save_for_restart opening file')
    call debug_message(verb, 'file proc is')
    call debug_message(verb, trim(adjustl(file_proc)))

    if (is_one_file_per_processor) then
       istatus = nf90_create (file_proc, create_file_mode, ncid)
    else
       call debug_message(verb, &
            'gs2_save::gs2_save_for_restart calling barrier before delete file')
       call barrier
       call debug_message(verb, &
            'gs2_save::gs2_save_for_restart called barrier before delete file')

       if (proc0) then
          open(unit=tmpunit, file=file_proc)
          close(unit=tmpunit, status='delete')
       end if

       call barrier
       call debug_message(verb, &
            'gs2_save::gs2_save_for_restart called barrier before opening')

# ifdef NETCDF_PARALLEL          
       istatus = nf90_create_par (file_proc, ior(create_file_mode, NF90_MPIIO), ncid=ncid, comm=mp_comm, info=mp_info)
# endif          
    end if

    if (istatus /= NF90_NOERR) then
       write(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
       goto 1
    end if

    call create_metadata(ncid, "GS2 restart file", local_header)

    if (.not.is_one_file_per_processor) then
# ifdef NETCDF_PARALLEL
       istatus = nf90_put_att(ncid, NF90_GLOBAL, 'layout', layout)
       if (istatus /= NF90_NOERR) then
          write(ierr,*) "nf90_put_attr error: ", nf90_strerror(istatus)
          goto 1
       end if
# endif
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart defining dimensions')

    if(is_one_file_per_processor) then
       call ensure_netcdf_dim_exists(ncid, "glo", n_elements, gloid)
    else
       call ensure_netcdf_dim_exists(ncid, "glo", total_elements, gloid)
    endif

    call ensure_netcdf_dim_exists(ncid, "theta", 2*ntgrid+1, thetaid)
    call ensure_netcdf_dim_exists(ncid, "sign", 2, signid)
    call ensure_netcdf_dim_exists(ncid, "aky", naky, kyid)
    call ensure_netcdf_dim_exists(ncid, "akx", ntheta0, kxid)

    ! Decide if this processor should take part in defining items
    ! which should only appear once in outputs.
    should_define_single_entry_variables = proc0 .or. (.not. is_one_file_per_processor)

    !For saving distribution function
    IF (save_glo_info_and_grids_local .and. should_define_single_entry_variables) THEN
       call ensure_netcdf_dim_exists(ncid, "negrid", negrid, egridid)
       call ensure_netcdf_dim_exists(ncid, "nlambda", nlambda, lgridid)
       call ensure_netcdf_dim_exists(ncid, "nspec", nspec, nspecid)
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart defining variables')

    call ensure_netcdf_dim_exists(ncid, "char5", 5, char5_dim)

    call ensure_netcdf_var_exists(ncid, "layout", nf90_char, char5_dim, layout_id)
    call ensure_netcdf_var_exists(ncid, "ntheta0", netcdf_int, ntheta0_id)
    call ensure_netcdf_var_exists(ncid, "naky", netcdf_int, naky_id)
    call ensure_netcdf_var_exists(ncid, "nlambda", netcdf_int, nlambda_id)
    call ensure_netcdf_var_exists(ncid, "negrid", netcdf_int, negrid_id)
    call ensure_netcdf_var_exists(ncid, "nspec", netcdf_int, nspec_id)

    call ensure_netcdf_var_exists(ncid, "t0", netcdf_real, t0id)
    call ensure_netcdf_var_exists(ncid, "delt0", netcdf_real, delt0id)
    call ensure_netcdf_var_exists(ncid, "delt1", netcdf_real, delt1id)
    call ensure_netcdf_var_exists(ncid, "delt2", netcdf_real, delt2id)
    call ensure_netcdf_var_exists(ncid, "delt_max", netcdf_real, delt_max_id)

    if (include_parameter_scan) then
       call ensure_netcdf_var_exists(ncid, "current_scan_parameter_value", &
            netcdf_real, current_scan_parameter_value_id)
    end if

    call ensure_netcdf_var_exists(ncid, "vnm1", netcdf_real, vnm1id)
    call ensure_netcdf_var_exists(ncid, "vnm2", netcdf_real, vnm2id)

    if (ant_on) then
       call ensure_netcdf_dim_exists(ncid, "nk_stir", nk_stir, nk_stir_dim)
       call ensure_netcdf_var_exists(ncid, "a_ant_r", netcdf_real, nk_stir_dim, a_antr_id)
       call ensure_netcdf_var_exists(ncid, "a_ant_i", netcdf_real, nk_stir_dim, a_anti_id)
       call ensure_netcdf_var_exists(ncid, "b_ant_r", netcdf_real, nk_stir_dim, b_antr_id)
       call ensure_netcdf_var_exists(ncid, "b_ant_i", netcdf_real, nk_stir_dim, b_anti_id)
    end if

    call ensure_netcdf_var_exists(ncid, "gr", netcdf_real, &
         [thetaid, signid, gloid], gr_id)
    call ensure_netcdf_var_exists(ncid, "gi", netcdf_real, &
         [thetaid, signid, gloid], gi_id)

    if(include_explicit_source_in_restart) then

#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          call ensure_netcdf_var_exists(ncid, "gexp_1_r", netcdf_real, &
               [thetaid, signid, gloid], gexp_1_r_id)
          call ensure_netcdf_var_exists(ncid, "gexp_1_i", netcdf_real, &
               [thetaid, signid, gloid], gexp_1_i_id)
       endif

       if(allocated(gexp_2)) then
          call ensure_netcdf_var_exists(ncid, "gexp_2_r", netcdf_real, &
               [thetaid, signid, gloid], gexp_2_r_id)
          call ensure_netcdf_var_exists(ncid, "gexp_2_i", netcdf_real, &
               [thetaid, signid, gloid], gexp_2_i_id)
       endif

       if(allocated(gexp_3)) then
          call ensure_netcdf_var_exists(ncid, "gexp_3_r", netcdf_real, &
               [thetaid, signid, gloid], gexp_3_r_id)
          call ensure_netcdf_var_exists(ncid, "gexp_3_i", netcdf_real, &
               [thetaid, signid, gloid], gexp_3_i_id)
       endif
    end if

    if (this_proc_saves_fields) then
       if (fphi > epsilon(0.)) then
          call ensure_netcdf_var_exists(ncid, "phi_r", netcdf_real, &
               [thetaid, kxid, kyid], phir_id)
          call ensure_netcdf_var_exists(ncid, "phi_i", netcdf_real, &
               [thetaid, kxid, kyid], phii_id)
       end if

       if (fapar > epsilon(0.)) then
          call ensure_netcdf_var_exists(ncid, "apar_r", netcdf_real, &
               [thetaid, kxid, kyid], aparr_id)
          call ensure_netcdf_var_exists(ncid, "apar_i", netcdf_real, &
               [thetaid, kxid, kyid], apari_id)
       end if

       if (fbpar > epsilon(0.)) then
          call ensure_netcdf_var_exists(ncid, "bpar_r", netcdf_real, &
               [thetaid, kxid, kyid], bparr_id)
          call ensure_netcdf_var_exists(ncid, "bpar_i", netcdf_real, &
               [thetaid, kxid, kyid], bpari_id)
       end if
    end if

    !For saving distribution function
    IF (save_glo_info_and_grids_local) THEN
       ! only store layout information on processor 0
       if (should_define_single_entry_variables) then
          if(.not. allocated(iglo_kykxsel)) allocate(iglo_kykxsel(naky,ntheta0,nspec,negrid,nlambda))
          if(.not. allocated(procid_kykxsel)) allocate(procid_kykxsel(naky,ntheta0,nspec,negrid,nlambda))

          do iglo=g_lo%llim_world,g_lo%ulim_world
             iglo_kykxsel(ik_idx(g_lo,iglo),it_idx(g_lo,iglo),is_idx(g_lo,iglo), &
                  ie_idx(g_lo,iglo),il_idx(g_lo,iglo)) = iglo
          end do

          do iglo=g_lo%llim_world,g_lo%ulim_world
             procid_kykxsel(ik_idx(g_lo,iglo),it_idx(g_lo,iglo),is_idx(g_lo,iglo), &
                  ie_idx(g_lo,iglo),il_idx(g_lo,iglo)) = proc_id(g_lo,iglo)
          end do

          !define layout information for plotting
          call ensure_netcdf_var_exists(ncid, "iglo_kykxsel", netcdf_int, &
               [kyid, kxid, nspecid, egridid, lgridid], iglo_kykxsel_id)
          call ensure_netcdf_var_exists(ncid, "procid_kykxsel", netcdf_int, &
               [kyid, kxid, nspecid, egridid, lgridid], procid_kykxsel_id)
       endif

       call ensure_netcdf_var_exists(ncid, "iglo_llim_proc", netcdf_int, iglo_llim_proc_id)

       ! only store layout information on processor 0
       if(should_define_single_entry_variables) then
          call ensure_netcdf_var_exists(ncid, "jtwist", netcdf_int, jtwist_id)
          call ensure_netcdf_var_exists(ncid, "theta_grid", netcdf_real, [thetaid], theta_grid_id)
          call ensure_netcdf_var_exists(ncid, "ky_grid", netcdf_real, [kyid], ky_grid_id)
          call ensure_netcdf_var_exists(ncid, "kx_grid", netcdf_real, [kxid], kx_grid_id)
          call ensure_netcdf_var_exists(ncid, "energy", netcdf_real, [egridid], energy_id)
          call ensure_netcdf_var_exists(ncid, "lambda", netcdf_real, [lgridid], lambda_id)
          call ensure_netcdf_var_exists(ncid, "species", netcdf_real, [nspecid], spec_id)
       endif

    END IF

    if(save_velocities_info) then
       call ensure_netcdf_var_exists(ncid, "vpa", netcdf_real, &
            [thetaid, signid, gloid], vpa_id)

       call ensure_netcdf_var_exists(ncid, "vperp2", netcdf_real, &
            [thetaid, gloid], vperp2_id)
    endif

    ! Because we want to be able to restart using exb shear from a
    ! case which does not have exb shear we always add kx_shift and
    ! theta0_shift even if no exb shear present in simulation)
    call ensure_netcdf_var_exists(ncid, "kx_shift", netcdf_real, &
         [kyid], kx_shift_id)
    call ensure_netcdf_var_exists(ncid, "theta0_shift", netcdf_real, &
         [kyid], theta0_shift_id)

    ! Now we've finished defining things switch to put mode
    call debug_message(verb, 'gs2_save::gs2_save_for_restart nf90_enddef')
    istatus = nf90_enddef (ncid)
    if (istatus /= NF90_NOERR) then
       write (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
       goto 1
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart writing scalars')
    
    if (is_one_file_per_processor .or. proc0) then

       istatus = nf90_put_var (ncid, delt0id, code_dt)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, delt1id, code_dt_prev1)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var delt1 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, delt2id, code_dt_prev2)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var delt2 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, delt_max_id, code_dt_max)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, t0id, t0)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, layout_id, layout)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var layout error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, ntheta0_id, ntheta0)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var ntheta0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, naky_id, naky)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var naky error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, nlambda_id, nlambda)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var nlambda error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, negrid_id, negrid)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var negrid error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, nspec_id, nspec)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var nspec error: ", nf90_strerror(istatus)
          goto 1
       end if

       if (include_parameter_scan) then
          ! <EGH see parameter_scan.f90
          istatus = nf90_put_var (ncid, &
               current_scan_parameter_value_id, current_scan_parameter_value)
          if (istatus /= NF90_NOERR) then
             write (ierr,*) "nf90_put_var current_scan_parameter_value error: ", nf90_strerror(istatus)
             goto 1
          end if
          ! EGH>
       end if

       istatus = nf90_put_var (ncid, vnm1id, vnm(1))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var vnm(1) error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, vnm2id, vnm(2))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var vnm(2) error: ", nf90_strerror(istatus)
          goto 1
       end if

    endif

    ! Here is the label which we head to if we detect an error at any point before here.
    ! Note this appears mid-way through our variable put phase, which seems a little odd.
    ! We therefore don't exit early if we encounter an error in later put statements.
1   continue

    if (istatus /= NF90_NOERR) then
       call debug_message(verb, 'gs2_save::gs2_save_for_restart closing file at the continue point')
       i = nf90_close (ncid)

       ! Deallocate module arrays
       call deallocate_arrays

       return
    end if

    if (ant_on) then

       if (.not. allocated(atmp)) allocate (atmp(nk_stir))
       atmp = real(a_ant)
       istatus = nf90_put_var (ncid, a_antr_id, atmp)

       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var a_antr error: ", &
               nf90_strerror(istatus), ' ', iproc
       end if

       atmp = aimag(a_ant)
       istatus = nf90_put_var (ncid, a_anti_id, atmp)

       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var a_anti error: ", &
               nf90_strerror(istatus), ' ', iproc
       end if

       atmp = real(b_ant)
       istatus = nf90_put_var (ncid, b_antr_id, atmp)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var b_antr error: ", &
               nf90_strerror(istatus), ' ', iproc
       end if

       atmp = aimag(b_ant)
       istatus = nf90_put_var (ncid, b_anti_id, atmp)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var b_anti error: ", &
               nf90_strerror(istatus), ' ', iproc
       end if
       deallocate (atmp)
    end if

    if (.not. allocated(tmpr)) &
         allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    tmpr = real(g)
    call debug_message(verb, 'gs2_save::gs2_save_for_restart writing dist fn')


    if (is_one_file_per_processor) then
       istatus = nf90_put_var (ncid, gr_id, tmpr)
    else
#ifdef NETCDF_PARALLEL       
       istatus = nf90_var_par_access(ncid, gr_id, NF90_COLLECTIVE)
       istatus = nf90_var_par_access(ncid, gi_id, NF90_COLLECTIVE)
# endif
       start_pos = (/1,1,g_lo%llim_proc+1/)
       counts = (/2*ntgrid+1, 2, n_elements/)

       istatus = nf90_put_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
    endif

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id, &
         var='gr put')

    tmpr = aimag(g)

    if (is_one_file_per_processor) then
       istatus = nf90_put_var (ncid, gi_id, tmpr)
    else
       istatus = nf90_put_var (ncid, gi_id, tmpr, start=start_pos, count=counts)
    endif

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id, &
         var = 'gi put')


    if (include_explicit_source_in_restart) then
       !Explicit source term gexp_1
#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          tmpr = real(gexp_1)

          if (is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_1_r_id, tmpr)
          else
#ifdef NETCDF_PARALLEL
             istatus = nf90_var_par_access(ncid, gexp_1_r_id, NF90_COLLECTIVE)
             istatus = nf90_var_par_access(ncid, gexp_1_i_id, NF90_COLLECTIVE)
#endif             
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             
             istatus = nf90_put_var (ncid, gexp_1_r_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_1_r_id, &
               var='gexp_1_r put')
          
          tmpr = aimag(gexp_1)

          if (is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_1_i_id, tmpr)
          else
             istatus = nf90_put_var (ncid, gexp_1_i_id, tmpr, start=start_pos, count=counts)
          endif

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_1_i_id, &
               var = 'gexp_1_i put')
       end if

       !Explicit source term gexp_2
       if(allocated(gexp_2)) then
          tmpr = real(gexp_2)
          
          if(is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_2_r_id, tmpr)
          else
#ifdef NETCDF_PARALLEL             
             istatus = nf90_var_par_access(ncid, gexp_2_r_id, NF90_COLLECTIVE)
             istatus = nf90_var_par_access(ncid, gexp_2_i_id, NF90_COLLECTIVE)
#endif             
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             
             istatus = nf90_put_var (ncid, gexp_2_r_id, tmpr, start=start_pos, count=counts)
          endif

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_2_r_id, &
               var='gexp_2_r put')

          tmpr = aimag(gexp_2)

          if (is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_2_i_id, tmpr)
          else
             istatus = nf90_put_var (ncid, gexp_2_i_id, tmpr, start=start_pos, count=counts)
          endif

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_2_i_id, &
               var = 'gexp_2_i put')
       end if

       !Explicit source term gexp_3
       if(allocated(gexp_3)) then
          tmpr = real(gexp_3)

          if (is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_3_r_id, tmpr)
          else
#ifdef NETCDF_PARALLEL            
             istatus = nf90_var_par_access(ncid, gexp_3_r_id, NF90_COLLECTIVE)
             istatus = nf90_var_par_access(ncid, gexp_3_i_id, NF90_COLLECTIVE)
# endif
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             
             istatus = nf90_put_var (ncid, gexp_3_r_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_3_r_id, &
               var='gexp_3_r put')
          
          tmpr = aimag(gexp_3)

          if (is_one_file_per_processor) then
             istatus = nf90_put_var (ncid, gexp_3_i_id, tmpr)
          else
             istatus = nf90_put_var (ncid, gexp_3_i_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_3_i_id, &
               var = 'gexp_3_i put')
       end if
    end if

    ! For saving distribution function
    if (save_glo_info_and_grids_local) then
       ! Fill processor and layout information for plotting purposes
       if (proc0) then
          istatus = nf90_put_var (ncid, iglo_kykxsel_id, iglo_kykxsel)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, iglo_kykxsel_id)


          istatus = nf90_put_var (ncid, procid_kykxsel_id, procid_kykxsel)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, procid_kykxsel_id)
       end if

       !Store variable iglo_llim_proc
       istatus = nf90_put_var (ncid, iglo_llim_proc_id, g_lo%llim_proc)

       !Check store was successful
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, iglo_llim_proc_id)

       if (proc0) then ! only save grids on proc0
          !Store variable jtwist
          istatus = nf90_put_var (ncid, jtwist_id, jtwist_out)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, jtwist_id)

          !Store variable theta_grid
          istatus = nf90_put_var (ncid, theta_grid_id, theta)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, theta_grid_id)

          !Store variable ky_grid
          istatus = nf90_put_var (ncid, ky_grid_id, aky)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, ky_grid_id)

          !Store variable kx_grid
          istatus = nf90_put_var (ncid, kx_grid_id, akx)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_grid_id)
          ! Fill energy and lambda information

          !Store variable energy
          istatus = nf90_put_var (ncid, energy_id, energy)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, energy_id)

          !Store variable lambda
          istatus = nf90_put_var (ncid, lambda_id, al)

          !Check store was successful
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, lambda_id)
       end if
    end if

    if(save_velocities_info) then
       !Fill velocity variables

       !Store variable vpa
       istatus = nf90_put_var (ncid, vpa_id, vpa)

       !Check store was successful
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, vpa_id)

       !Store variable vperp2
       istatus = nf90_put_var (ncid, vperp2_id, vperp2)

       !Check store was successful
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, vperp2_id)
    end if

    if (is_one_file_per_processor .or. proc0) then
       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))

       if (this_proc_saves_fields) then
          if (fphi > epsilon(0.)) then
             ftmpr = real(phinew)
             istatus = nf90_put_var (ncid, phir_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id, &
                  var = 'phir put')

             ftmpr = aimag(phinew)
             istatus = nf90_put_var (ncid, phii_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id, &
                  var = 'phir put')
          end if

          if (fapar > epsilon(0.)) then
             ftmpr = real(aparnew)
             istatus = nf90_put_var (ncid, aparr_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)

             ftmpr = aimag(aparnew)
             istatus = nf90_put_var (ncid, apari_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)
          end if

          if (fbpar > epsilon(0.)) then
             ftmpr = real(bparnew)
             istatus = nf90_put_var (ncid, bparr_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)

             ftmpr = aimag(bparnew)
             istatus = nf90_put_var (ncid, bpari_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)
          end if
       end if

       if (.not. allocated(stmp)) allocate (stmp(naky))
       if (allocated(kx_shift)) then
          stmp = kx_shift
       else
          stmp = 0.
       end if
       istatus = nf90_put_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error &
            (istatus, ncid, kx_shift_id, var = 'kx_shift put')

       if (allocated(theta0_shift)) then
          stmp = theta0_shift
       else
          stmp = 0.
       end if
       istatus = nf90_put_var (ncid, theta0_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error &
            (istatus, ncid, theta0_shift_id, var = 'theta0_shift put')
    end if

    ! Close the file
    i = nf90_close (ncid)
    if (i /= NF90_NOERR) &
         call netcdf_error (istatus, message='nf90_close error')

    ! Deallocate local arrays
    if (allocated(iglo_kykxsel)) deallocate(iglo_kykxsel)
    if (allocated(procid_kykxsel)) deallocate(procid_kykxsel)

    ! Deallocate module arrays
    call deallocate_arrays

# else
    if (proc0) write (error_unit(),*) &
         'WARNING: gs2_save_for_restart is called without netcdf library'
# endif
  end subroutine gs2_save_for_restart

  !> Save g slices
  subroutine gs2_save_slice(g_in, it, ik, il, ie, is, time,istatus, nout, fileopt, header)
# ifdef NETCDF
    use netcdf_utils, only: ensure_netcdf_var_exists, ensure_netcdf_dim_exists, netcdf_real
    use gs2_metadata, only: create_metadata
    use netcdf, only: NF90_CLOBBER, nf90_create, nf90_open, NF90_WRITE, NF90_NETCDF4, nf90_put_var
    use file_utils, only: error_unit
    use unit_tests, only: debug_message
    use mp, only: iproc, proc0
    use optionals, only: get_option_with_default
    use convert, only: c2r
    use kt_grids, only: aky, akx, naky, ntheta0
    use le_grids, only: energy, al, negrid, nlambda
    use species, only: nspec
# endif
    use theta_grid, only: ntgrid, theta
    use gs2_layouts, only: g_lo, proc_id, idx
    use standard_header, only : standard_header_type
    implicit none
    !> The distribution function \(g\)
    complex, dimension (-ntgrid:, :, g_lo%llim_proc:), intent (in) :: g_in
    !> The indices of the slice to save
    integer, intent(in) :: it, ik, il, ie, is
    !> The current time
    real, intent(in) :: time
    !> Last netCDF status: zero on success, non-zero on failure
    integer, intent (out) :: istatus
    !> The current time output index
    integer, intent(in), optional :: nout
    !> Optional extra filename infix
    character (len=*), intent (in), optional :: fileopt
    !> Header for files with build and run information
    type(standard_header_type), optional, intent(in) :: header
# ifdef NETCDF
    character (run_name_size) :: file_proc
    integer :: ierr, the_iglo, the_proc, r_handle, s_handle
    integer, parameter :: verb = 3
    ! Netcdf ids for variables only handled in this routine
    integer (kind_nf) :: theta_id, time_id, gslice_id
    integer (kind_nf) :: it_id, ik_id, il_id, ie_id, is_id, the_iglo_id, the_proc_id
    integer (kind_nf) :: kx_id, ky_id, lambda_id, energy_id
    integer (kind_nf) :: ri_dim, sign_dim, theta_dim, time_dim
    !> Open new netCDF files with netCDF-4 model, and overwrite existing files
    integer :: create_file_mode
    integer, dimension(4) :: gslice_dims
    type(standard_header_type) :: local_header
    character (len=*), parameter :: file_opt_base = '.slice'
    logical :: appending, previous_file_exists
    real, dimension(:, :, :), allocatable :: ri_gslice
    integer, save :: nout_local = 1
    logical :: error_exit
    istatus = 0
    call debug_message(verb, 'gs2_save::gs2_save_slice initialising')

    ! Sanity check the inputs
    error_exit = .false.
    if (it < 1 .or. it > ntheta0) error_exit = .true.
    if (ik < 1 .or. ik > naky) error_exit = .true.
    if (il < 1 .or. il > nlambda) error_exit = .true.
    if (ie < 1 .or. ie > negrid) error_exit = .true.
    if (is < 1 .or. is > nspec) error_exit = .true.
    if (proc0 .and. error_exit) write(6,'("Error: Invalid index based to gs2_save_slice.")')
    if (error_exit) return

    the_iglo = idx(g_lo, ik, it, il, ie, is)
    the_proc = proc_id(g_lo, the_iglo)

    ! The proc which owns this slice is responsible for writing and
    ! does the rest of the work so let us return if that's not this
    ! proc
    if (iproc /= the_proc) return

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    call debug_message(verb, 'gs2_save::gs2_save_slice initialising')

    ierr = error_unit()

    ! Open the file and define dimensions and variables
    file_proc = get_file_proc(.false., file_opt_base // fileopt)

    call debug_message(verb, 'gs2_save::gs2_save_slice opening file')
    call debug_message(verb, 'file proc is')
    call debug_message(verb, trim(adjustl(file_proc)))

    nout_local = get_option_with_default(nout, nout_local)

    ! See if the output file already exists, so we know how to open
    ! the file.  This would be simpler if netCDF had a "open if it
    ! exists, create if it doesn't" mode, but it doesn't so we have
    ! to do this instead
    inquire(file=trim(file_proc), exist=previous_file_exists)
    appending = (previous_file_exists .and. nout_local > 1)

    if (.not. appending) then
       create_file_mode = ior(NF90_NETCDF4, NF90_CLOBBER)
       istatus = nf90_create (file_proc, create_file_mode, ncid)
    else
       create_file_mode = NF90_WRITE
       istatus = nf90_open (file_proc, create_file_mode, ncid)
    end if

    if (istatus /= NF90_NOERR) then
       write(ierr,*) "nf90_create/open error: ", nf90_strerror(istatus)
       goto 1
    end if

    call create_metadata(ncid, "GS2 slice file", local_header)

    call debug_message(verb, 'gs2_save::gs2_save_slice defining dimensions')

    call ensure_netcdf_dim_exists(ncid, 'ri', 2, ri_dim)
    call ensure_netcdf_dim_exists(ncid, "theta", 2*ntgrid+1, theta_dim)
    call ensure_netcdf_dim_exists(ncid, "sign", 2, sign_dim)
    call ensure_netcdf_dim_exists(ncid, 't', NF90_UNLIMITED, time_dim)

    gslice_dims = [ri_dim, theta_dim, sign_dim, time_dim]

    call ensure_netcdf_var_exists(ncid, "t", netcdf_real, [time_dim], time_id)
    call ensure_netcdf_var_exists(ncid, "g_slice", netcdf_real, gslice_dims, gslice_id)
    call ensure_netcdf_var_exists(ncid, "theta", netcdf_real, [theta_dim], theta_id)

    call ensure_netcdf_var_exists(ncid, 'it', netcdf_int, it_id)
    call ensure_netcdf_var_exists(ncid, 'ik', netcdf_int, ik_id)
    call ensure_netcdf_var_exists(ncid, 'il', netcdf_int, il_id)
    call ensure_netcdf_var_exists(ncid, 'ie', netcdf_int, ie_id)
    call ensure_netcdf_var_exists(ncid, 'is', netcdf_int, is_id)

    call ensure_netcdf_var_exists(ncid, 'iglo', netcdf_int, the_iglo_id)
    call ensure_netcdf_var_exists(ncid, 'iproc', netcdf_int, the_proc_id)

    call ensure_netcdf_var_exists(ncid, 'kx', netcdf_real, kx_id)
    call ensure_netcdf_var_exists(ncid, 'ky', netcdf_real, ky_id)
    call ensure_netcdf_var_exists(ncid, 'lambda', netcdf_real, lambda_id)
    call ensure_netcdf_var_exists(ncid, 'energy', netcdf_real, energy_id)

    call debug_message(verb, 'gs2_save::gs2_save_slice writing variables')

    ! Write constant data
    if (nout_local == 1) then
       istatus = nf90_put_var (ncid, theta_id, theta)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var theta error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, it_id, it)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var it error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, ik_id, ik)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var ik error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, il_id, il)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var il error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, ie_id, ie)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var ie error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, is_id, is)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var is error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, the_iglo_id, the_iglo)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var iglo error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, the_proc_id, the_proc)
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var iproc error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, kx_id, akx(it))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var kx error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, ky_id, aky(ik))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var ky error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, lambda_id, al(il))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var lambda error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, energy_id, energy(ie, is))
       if (istatus /= NF90_NOERR) then
          write (ierr,*) "nf90_put_var energy error: ", nf90_strerror(istatus)
          goto 1
       end if

    end if

    istatus = nf90_put_var (ncid, time_id, time, start = [nout_local])
    if (istatus /= NF90_NOERR) then
       write (ierr,*) "nf90_put_var time error: ", nf90_strerror(istatus)
       goto 1
    end if

    ! Convert from complex to real
    allocate(ri_gslice(2, size(theta), 2))
    call c2r(g_in(:, :, the_iglo), ri_gslice)
    istatus = nf90_put_var (ncid, gslice_id, ri_gslice, &
         start = [1, 1, 1, nout_local], &
         count = [2, size(theta), 2, 1] &
         )
    if (istatus /= NF90_NOERR) then
       write (ierr,*) "nf90_put_var gslice error: ", nf90_strerror(istatus)
       goto 1
    end if

    ! Here is the label which we head to if we detect an error at any point before here.
    ! Note this appears mid-way through our variable put phase, which seems a little odd.
    ! We therefore don't exit early if we encounter an error in later put statements.
1   continue

    if (istatus /= NF90_NOERR) then
       call debug_message(verb, 'gs2_save::gs2_save_slice closing file at the continue point')
       istatus = nf90_close (ncid)
       return
    end if

    ! Close the file
    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) &
         call netcdf_error (istatus, message='nf90_close error')

    nout_local = nout_local + 1
# endif
  end subroutine gs2_save_slice

  !> FIXME : Add documentation
  subroutine gs2_restore_many (g, scale, istatus, fphi, fapar, fbpar, fileopt)
!MR, 2007: restore kx_shift array if already allocated
# ifdef NETCDF
    use mp, only: iproc, mp_abort, proc0, broadcast
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: layout
# endif
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use file_utils, only: error_unit
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    logical :: this_proc_reads_fields, need_to_broadcast_fields
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar, fbpar
    character (len=*), intent(in), optional :: fileopt
# ifdef NETCDF
    !> This parameter controls whether or not gs2
    !> aborts if it cannot read the restart file (and
    !> it has been told to load g from  a restart).
    !> It is set to .true..
    !> I can think of no conceivable reason why this
    !> would ever need to be .false., but comments welcome. EGH
    logical, parameter :: abort_on_restart_fail = .true.
    integer, dimension(3) :: counts, start_pos
    character (run_name_size) :: file_proc
    character (len = 5) :: restart_layout
    integer :: i, n_elements
    real :: fac
    logical :: has_explicit_source_terms
    logical :: is_one_file_per_processor
    logical :: has_gexp_1
    logical, save :: explicit_warning_given = .false.
    n_elements = max(g_lo%ulim_proc-g_lo%llim_proc+1, 0)

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    ! Decide if we are reading the fields
    need_to_broadcast_fields = proc_to_save_fields >= 0
    this_proc_reads_fields = (iproc == proc_to_save_fields) .or. (.not. need_to_broadcast_fields)
    
    file_proc = get_file_proc(is_one_file_per_processor, fileopt)

    if (is_one_file_per_processor) then
       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
    else
# ifdef NETCDF_PARALLEL          
       ! If using netcdf version 4.1.2 deleted NF90_MPIIO and the associated IOR
       istatus = nf90_open (file_proc, IOR(NF90_NOWRITE, NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
# endif          
    end if

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)

    call check_netcdf_file_precision (ncid)

    istatus = nf90_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='theta',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_dimid (ncid, "sign", signid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='sign',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_dimid (ncid, "glo", gloid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='glo',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_dimid (ncid, "aky", kyid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='aky',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_dimid (ncid, "akx", kxid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='akx',&
         abort=abort_on_restart_fail)

    istatus = nf90_inquire_dimension (ncid, thetaid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=thetaid,&
         abort=abort_on_restart_fail, dim='theta')
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, signid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=signid,&
         abort=abort_on_restart_fail)
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, gloid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gloid,&
         abort=abort_on_restart_fail)

    if (is_one_file_per_processor) then
       if (i /= g_lo%ulim_proc-g_lo%llim_proc+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
    else
       if (i /= g_lo%ulim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
    endif

    istatus = nf90_inquire_dimension (ncid, kyid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kyid,&
         abort=abort_on_restart_fail)
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, kxid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kxid,&
         abort=abort_on_restart_fail)
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    if (this_proc_reads_fields) then
       if (fphi > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "phi_r", phir_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_r')

          istatus = nf90_inq_varid (ncid, "phi_i", phii_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_i')
       end if

       if (fapar > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_r')

          istatus = nf90_inq_varid (ncid, "apar_i", apari_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_i')
       end if

       if (fbpar > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_r')

          istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_i')
       end if
    end if

    if (allocated(kx_shift)) then   ! MR begin
       istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='kx_shift')
    endif   ! MR end

    if (allocated(theta0_shift)) then
       istatus = nf90_inq_varid (ncid, "theta0_shift", theta0_shift_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='theta0_shift')
    endif

    if (include_explicit_source_in_restart) then
       has_explicit_source_terms = .true.

#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          istatus = nf90_inq_varid (ncid, "gexp_1_r", gexp_1_r_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.

          istatus = nf90_inq_varid (ncid, "gexp_1_i", gexp_1_i_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.
       endif

       if (allocated(gexp_2)) then
          istatus = nf90_inq_varid (ncid, "gexp_2_r", gexp_2_r_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.

          istatus = nf90_inq_varid (ncid, "gexp_2_i", gexp_2_i_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.
       endif

       if (allocated(gexp_3)) then
          istatus = nf90_inq_varid (ncid, "gexp_3_r", gexp_3_r_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.

          istatus = nf90_inq_varid (ncid, "gexp_3_i", gexp_3_i_id)
          if (istatus /= NF90_NOERR) has_explicit_source_terms = .false.
       endif

       if ((.not. explicit_warning_given).and. proc0 .and. .not.has_explicit_source_terms) then
          write(*, '("Warning: At least one explicit source term absent in restart file.")')
          write(*, '("  Will not load these terms. This is probably OK unless you expect the restart file to contain these terms.")')
          write(*, '("  This warning will not be repeated.")')

          ! Update this flag to ensure we don't give this warning again
          explicit_warning_given = .true.
       end if

    else
       ! Always skip trying to get the explicit source terms if we're not
       ! including them
       has_explicit_source_terms = .false.
    end if

    istatus = nf90_inq_varid (ncid, "gr", gr_id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gr',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_varid (ncid, "gi", gi_id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gi',&
         abort=abort_on_restart_fail)

    istatus = nf90_inq_varid (ncid, "layout", layout_id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='layout_r')

    istatus = nf90_get_var (ncid, layout_id, restart_layout)

    if (istatus /= NF90_NOERR) then
       call netcdf_error (istatus, ncid, layout_id, message=' in gs2_restore_many')
       restart_layout = layout
    endif

    ! Abort if the layouts don't match
    if(restart_layout /= layout) then
       call mp_abort("Incompatible layouts in restart file ("//restart_layout//") and simulation ("//layout//")",.true.)
    end if

    if (.not. allocated(tmpr)) &
         allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))
    if (.not. allocated(tmpi)) &
         allocate (tmpi(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    tmpr = 0.; tmpi = 0.

    if (is_one_file_per_processor) then
       istatus = nf90_get_var (ncid, gr_id, tmpr)
    else
       start_pos = (/1,1,g_lo%llim_proc+1/)
       counts = (/2*ntgrid+1, 2, n_elements/)
       istatus = nf90_get_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
    end if

   if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id,&
         abort=abort_on_restart_fail, var='gr')

    if(is_one_file_per_processor) then
       istatus = nf90_get_var (ncid, gi_id, tmpi)
    else
       istatus = nf90_get_var (ncid, gi_id, tmpi, start=start_pos, count=counts)
    end if

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id,&
         abort=abort_on_restart_fail, var='gi')

    g = cmplx(tmpr, tmpi)

    if (include_explicit_source_in_restart .and. has_explicit_source_terms) then
       ! Explicit source term
#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          if (is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_1_r_id, tmpr)
          else
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             istatus = nf90_get_var (ncid, gexp_1_r_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_1_r_id,&
               var='gexp_1_r')

          if(is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_1_i_id, tmpi)
          else
             istatus = nf90_get_var (ncid, gexp_1_i_id, tmpi, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_1_i_id,&
               var='gexp_1_i')

          gexp_1 = cmplx(tmpr, tmpi)

       endif

       ! Explicit source term
       if(allocated(gexp_2)) then
          if (is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_2_r_id, tmpr)
          else
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             istatus = nf90_get_var (ncid, gexp_2_r_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_2_r_id,&
               var='gexp_2_r')

          if(is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_2_i_id, tmpi)
          else
             istatus = nf90_get_var (ncid, gexp_2_i_id, tmpi, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_2_i_id,&
               var='gexp_2_i')

          gexp_2 = cmplx(tmpr, tmpi)

       endif

       ! Explicit source term
       if(allocated(gexp_3)) then
          if (is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_3_r_id, tmpr)
          else
             start_pos = (/1,1,g_lo%llim_proc+1/)
             counts = (/2*ntgrid+1, 2, n_elements/)
             istatus = nf90_get_var (ncid, gexp_3_r_id, tmpr, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_3_r_id,&
               var='gexp_3_r')

          if (is_one_file_per_processor) then
             istatus = nf90_get_var (ncid, gexp_3_i_id, tmpi)
          else
             istatus = nf90_get_var (ncid, gexp_3_i_id, tmpi, start=start_pos, count=counts)
          end if

          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gexp_3_i_id,&
               var='gexp_3_i')

          gexp_3 = cmplx(tmpr, tmpi)

       endif
    end if

    if (allocated(kx_shift)) then   ! MR begin
       if (.not. allocated(stmp)) allocate (stmp(naky))   ! MR 
       istatus = nf90_get_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) &
         call netcdf_error (istatus, ncid, kx_shift_id, var='kx_shift')
       kx_shift = stmp
    endif   ! MR end

    if (allocated(theta0_shift)) then
       if (.not. allocated(stmp)) allocate (stmp(naky))
       istatus = nf90_get_var (ncid, theta0_shift_id, stmp)
       if (istatus /= NF90_NOERR) &
         call netcdf_error (istatus, ncid, theta0_shift_id, var='theta0_shift')
       theta0_shift = stmp
    endif

    if (this_proc_reads_fields) then
       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
       if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

       if (fphi > epsilon(0.)) then
          istatus = nf90_get_var (ncid, phir_id, ftmpr)
          if (istatus /= NF90_NOERR) &
               call netcdf_error (istatus, ncid, phir_id, var= 'phir')

          istatus = nf90_get_var (ncid, phii_id, ftmpi)
          if (istatus /= NF90_NOERR) &
               call netcdf_error (istatus, ncid, phii_id, var='phii')

          phinew = cmplx(ftmpr, ftmpi)
       end if

       if (fapar > epsilon(0.)) then
          istatus = nf90_get_var (ncid, aparr_id, ftmpr)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)

          istatus = nf90_get_var (ncid, apari_id, ftmpi)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)

          aparnew = cmplx(ftmpr, ftmpi)
       end if

       if (fbpar > epsilon(0.)) then
          istatus = nf90_get_var (ncid, bparr_id, ftmpr)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)

          istatus = nf90_get_var (ncid, bpari_id, ftmpi)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)

          bparnew = cmplx(ftmpr, ftmpi)
       end if
    end if

    if (need_to_broadcast_fields) then
       if (fphi > epsilon(0.)) call broadcast(phinew, proc_to_save_fields)
       if (fapar > epsilon(0.)) call broadcast(aparnew, proc_to_save_fields)
       if (fbpar > epsilon(0.)) call broadcast(bparnew, proc_to_save_fields)
    end if

    if (fphi > epsilon(0.)) phi = 0.
    if (fapar > epsilon(0.)) apar = 0.
    if (fbpar > epsilon(0.)) bpar = 0.

    if (scale > 0.) then
       g = g*scale
       phinew = phinew*scale
       aparnew = aparnew*scale
       bparnew = bparnew*scale
    else
       fac = - scale/(maxval(abs(phinew)))
       g = g*fac
       phinew = phinew*fac
       aparnew = aparnew*fac
       bparnew = bparnew*fac
    end if

    ! RN 2008/05/23: this was commented out. why? HJL 2013/05/15 Because it stops future writing to the file, it's now back in after setting initialized to false
    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) then
       write(error_unit(),*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    end if

    ! Deallocate modul level arrays
    call deallocate_arrays
    
# else
    write (error_unit(),*) &
         'ERROR: gs2_restore_many is called without netcdf'
# endif
  end subroutine gs2_restore_many

  !> This routine writes a passed square complex array to a file
  !> with passed name
  subroutine gs2_save_response(resp, fname, condition, header)
    use file_utils, only: error_unit
    use standard_header, only: standard_header_type
#ifdef NETCDF
    use gs2_time, only: code_dt
    use convert, only: c2r
    use netcdf_utils, only: ensure_netcdf_dim_exists, ensure_netcdf_var_exists
    use gs2_metadata, only: create_metadata
#else
    use file_utils, only: get_unused_unit
#endif
    implicit none
    complex,dimension(:,:), intent(in) :: resp
    character(len=*), intent(in) :: fname
    real, intent(in), optional :: condition
    !> Header for files with build and run information
    type(standard_header_type), optional, intent(in) :: header

    integer :: sz
#ifdef NETCDF
    real :: condition_local
    integer :: ierr, ax1id,ax2id,riid, respid,dtid, condition_id, ncid
    real, dimension(:,:,:), allocatable :: ri_resp
    ! Actual value for optional `header` input
    type(standard_header_type) :: local_header
#else
    integer :: unit
#endif
    !Currently only support serial writing, but could be by any proc
    !so we have to make sure only one proc calls this routine

    !Verify we have a square array
    sz=size(resp(:,1))
    if(sz.ne.size(resp(1,:))) then
       write(error_unit(),'("Error: gs2_save_response expects a square array input.")')
       return
    endif

#ifdef NETCDF
    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    !/Make file
    ierr=nf90_create(fname, create_file_mode, ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

    call create_metadata(ncid, "GS2 response matrix", local_header)

    !/Define dimensions
    call ensure_netcdf_dim_exists(ncid, "ri", 2, riid)
    call ensure_netcdf_dim_exists(ncid, "ax1", sz, ax1id)
    call ensure_netcdf_dim_exists(ncid, "ax2", sz, ax2id)

    !/Define variables
    call ensure_netcdf_var_exists(ncid, "response", netcdf_real, [riid,ax1id,ax2id], respid)
    call ensure_netcdf_var_exists(ncid, "dt", netcdf_real, dtid)
    call ensure_netcdf_var_exists(ncid, "condition", netcdf_real, condition_id)

    !End definitions
    ierr=nf90_enddef(ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

    !Now we can place our data in the file
    ierr=nf90_put_var(ncid,dtid,code_dt)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="dt")

    condition_local = -1.0
    if(present(condition)) condition_local = condition
    ierr=nf90_put_var(ncid,condition_id,condition_local)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="condtion")

    !/Convert complex to ri and write
    allocate(ri_resp(2,sz,sz))
    call c2r(resp,ri_resp)
    ierr=nf90_put_var(ncid,respid,ri_resp)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="response")
    deallocate(ri_resp)

    !/Now close the file
    ierr=nf90_close(ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

#else
!Fall back on binary output if no NETCDF
    !Get a free unit
    call get_unused_unit(unit)

    !Open file and write
    open(unit=unit,file=fname,form="unformatted")
    write(unit) resp
    close(unit)
#endif
  end subroutine gs2_save_response

  !> This routine reads a square complex array from a file
  !! with passed name
  subroutine gs2_restore_response(resp, fname, condition)
    use file_utils, only: error_unit
#ifdef NETCDF
    use convert, only: r2c
#else
    use file_utils, only: get_unused_unit
#endif
    implicit none
    complex, dimension(:,:), intent(out) :: resp
    character(len=*), intent(in) :: fname
    real, intent(out), optional :: condition
    integer :: sz
#ifdef NETCDF
    integer :: ierr, respid, ncid, condition_id
    real, dimension(:,:,:), allocatable :: ri_resp
#else
    integer :: unit
#endif
    !Currently only support serial reading, but could be by any proc
    !so we have to make sure only one proc calls this routine

    !Verify we have a square array
    sz=size(resp(:,1))
    if(sz.ne.size(resp(1,:))) then
       write(error_unit(),'("Error: gs2_restore_response expects a square array output.")')
       return
    endif

#ifdef NETCDF

    !/Open file
    ierr=nf90_open(fname,NF90_NOWRITE,ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

    !/Get variable id
    ierr=nf90_inq_varid(ncid,"response",respid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="response")

    !/Read and convert ri to complex
    allocate(ri_resp(2,sz,sz))
    ierr=nf90_get_var(ncid,respid,ri_resp)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="response")
    call r2c(resp,ri_resp)
    deallocate(ri_resp)

    !/Get the condition number
    if (present(condition)) then
       ierr = nf90_inq_varid(ncid, "condition", condition_id)
       if(ierr/=NF90_NOERR) call netcdf_error(ierr, var="condition")
       ierr=nf90_get_var(ncid, condition_id, condition)
       if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="condition")
    end if

    !/Now close the file
    ierr=nf90_close(ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

#else
!Fall back on binary output if no NETCDF
    !Get a free unit
    call get_unused_unit(unit)
    
    !Open file and write
    open(unit=unit,file=fname,form="unformatted")
    read(unit) resp
    close(unit)
#endif
  end subroutine gs2_restore_response

  !> Initialises a file for saving output of eigensolver to netcdf
  subroutine init_eigenfunc_file(fname,fphi,fapar,fbpar,IDs)
    use file_utils, only: error_unit
    use theta_grid, only: ntgrid, theta
#ifdef NETCDF
    use netcdf_utils, only: ensure_netcdf_dim_exists, ensure_netcdf_var_exists
#endif
    implicit none
    character(len=*), intent(in) :: fname
    type(EigNetcdfID), intent(inout) :: IDs
    real, intent(in) :: fphi, fapar, fbpar
#ifdef NETCDF
    integer :: ierr
#endif

    !Set nconv counter to 0
    IDs%nconv_count=0

#ifdef NETCDF
    !/Make file
    ierr=nf90_create(fname, create_file_mode, IDs%ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

    !/Define dimensions
    call ensure_netcdf_dim_exists(IDs%ncid, "ri", 2, IDs%ri_dim_id)
    call ensure_netcdf_dim_exists(IDs%ncid, "theta", 2*ntgrid+1, IDs%theta_dim_id)
    call ensure_netcdf_dim_exists(IDs%ncid, "nconv", NF90_UNLIMITED, IDs%conv_dim_id)

    !/Define variables
    !Dimensions
    call ensure_netcdf_var_exists(IDs%ncid, "theta", netcdf_real, &
         IDs%theta_dim_id, IDs%theta_var_id)
    call ensure_netcdf_var_exists(IDs%ncid, "conv", netcdf_real, &
         IDs%conv_dim_id, IDs%conv_var_id)
    !Fields
    if(fphi.gt.epsilon(0.0))then
       call ensure_netcdf_var_exists(IDs%ncid, "phi", netcdf_real, &
            [IDs%ri_dim_id, IDs%theta_dim_id, IDs%conv_dim_id], IDs%phi_var_id)
    endif
    if(fapar.gt.epsilon(0.0))then
       call ensure_netcdf_var_exists(IDs%ncid, "apar", netcdf_real, &
            [IDs%ri_dim_id, IDs%theta_dim_id, IDs%conv_dim_id], IDs%apar_var_id)
    endif
    if(fbpar.gt.epsilon(0.0))then
       call ensure_netcdf_var_exists(IDs%ncid, "bpar", netcdf_real, &
            [IDs%ri_dim_id, IDs%theta_dim_id, IDs%conv_dim_id], IDs%bpar_var_id)
    endif
    !Omega
    call ensure_netcdf_var_exists(IDs%ncid, "omega", netcdf_real, &
         [IDs%ri_dim_id, IDs%conv_dim_id], IDs%omega_var_id)
    
    !End definitions
    ierr=nf90_enddef(IDs%ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,file=fname)

    !Now we can place some data in the file
    ierr=nf90_put_var(IDs%ncid,IDs%theta_var_id,theta)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr,var="theta")

#endif
  end subroutine init_eigenfunc_file

  !> Add an eigenpairs data to file
  subroutine add_eigenpair_to_file(eval,fphi,fapar,fbpar,IDs,my_conv)
    use fields_arrays, only: phinew, aparnew, bparnew
    use convert, only: c2r
    use theta_grid, only: ntgrid
    complex, intent(in) :: eval !Note just use fields to get eigenvectors
    real, intent(in), optional :: my_conv
    type(EigNetcdfID), intent(inout) :: IDs
    real, intent(in) :: fphi, fapar, fbpar
#ifdef NETCDF
    real, dimension(2) :: ri_omega
    real, dimension(:,:), allocatable :: ri_field
    integer, dimension(3) :: start_field, count_field
    integer, dimension(2) :: start_omega
    integer :: ierr
    real :: local_conv
#endif

    !First increment counter
    IDs%nconv_count=IDs%nconv_count+1

#ifdef NETCDF
    !Now we make start values
    start_field(1)=1 ; start_field(2)=1 ; start_field(3)=IDs%nconv_count
    count_field(1)=2 ; count_field(2)=2*ntgrid+1 ; count_field(3)=1
    start_omega(1)=1 ; start_omega(2)=IDs%nconv_count

    !Set the conv value
    if(present(my_conv))then
       local_conv=my_conv
    else
       local_conv=IDs%nconv_count*1.0
    endif

    !Now we can write data
    !/Conv
    ierr=nf90_put_var(IDs%ncid,IDs%conv_var_id,local_conv, start=(/IDs%nconv_count/))!,count=(/1/))
    if(ierr/=NF90_NOERR) call netcdf_error (ierr, IDs%ncid, IDs%conv_var_id)

    !/Omega
    ri_omega(1)=real(eval)
    ri_omega(2)=aimag(eval)
    ierr=nf90_put_var(IDs%ncid,IDs%omega_var_id,ri_omega, start=start_omega,count=(/2,1/))
    if(ierr/=NF90_NOERR) call netcdf_error (ierr, IDs%ncid, IDs%omega_var_id)

    !/Fields
    allocate(ri_field(2,2*ntgrid+1))
    if(fphi.gt.epsilon(0.0))then
       call c2r(phinew(:,1,1),ri_field)
       ierr=nf90_put_var(IDs%ncid,IDs%phi_var_id,ri_field, start=start_field,count=count_field)
       if(ierr/=NF90_NOERR) call netcdf_error (ierr, IDs%ncid, IDs%phi_var_id)
    endif
    if(fapar.gt.epsilon(0.0))then
       call c2r(aparnew(:,1,1),ri_field)
       ierr=nf90_put_var(IDs%ncid,IDs%apar_var_id,ri_field, start=start_field,count=count_field)
       if(ierr/=NF90_NOERR) call netcdf_error (ierr, IDs%ncid, IDs%apar_var_id)
    endif
    if(fbpar.gt.epsilon(0.0))then
       call c2r(bparnew(:,1,1),ri_field)
       ierr=nf90_put_var(IDs%ncid,IDs%bpar_var_id,ri_field, start=start_field,count=count_field)
       if(ierr/=NF90_NOERR) call netcdf_error (ierr, IDs%ncid, IDs%bpar_var_id)
    endif
    deallocate(ri_field)       
#endif
  end subroutine add_eigenpair_to_file

  !> Close the eigenfunction file
  subroutine finish_eigenfunc_file(IDs)
    implicit none
    type(EigNetcdfID), intent(inout) :: IDs
#ifdef NETCDF
    integer :: ierr

    !/Now close the file
    ierr=nf90_close(IDs%ncid)
    if(ierr/=NF90_NOERR) call netcdf_error(ierr)
#endif
  end subroutine finish_eigenfunc_file

  !> This function checks to see if we can create a file with name
  !! <restart_file>//<SomeSuffix> if not then our restarts are not
  !! going to be possible and we return false. Can also be used to check
  !! that we can read from the restart file (which assumes it exists).
  function restart_writable(read_only, my_file, error_message)
    use mp, only: proc0, broadcast
    use file_utils, only: get_unused_unit
    use constants, only: run_name_size
    implicit none
    !> If present and true, only check that files can be read
    logical, intent(in), optional :: read_only
    !> An optional specific filename to check
    character(len=*),intent(in),optional::my_file
    !> Error message returned from `open` if there was a problem
    character(len=:), allocatable, optional :: error_message

    character(len=*), parameter :: SuffixTmp = '.ThisIsATestFile'
    character(9) :: open_mode
    character(6) :: close_mode
    character(run_name_size) :: local_file, filename
    logical :: restart_writable
    integer :: unit, ierr
    character(len=200) :: io_err_msg

    ierr=-200
    local_file = trim(restart_file)
    if (present(my_file)) local_file = trim(my_file)

    ! On proc0 try to open tmp file for writing
    if (proc0)then
       call get_unused_unit(unit)

       ! Default open and close modes
       open_mode = 'readwrite'
       close_mode = 'delete'

       ! If we want to test write capability then do it with an unusual
       ! file name to prevent clobber
       filename = trim(local_file) // SuffixTmp

       ! Set filemode to READ if read_only=T
       if (present(read_only)) then
         if (read_only) then
           ! If we're only checking a file can be read, then don't
           ! delete the file after we're done
           open_mode = 'read'
           close_mode = 'keep'
           ! If checking readonly then we need to make sure we try to
           ! read from an existing file
           filename = local_file
         end if
       endif

       open(unit=unit, file=trim(filename), &
            iostat=ierr, action=open_mode, &
            iomsg=io_err_msg)

       ! If open was successful then we can close the file and delete it
       if (ierr == 0) close(unit=unit, status=close_mode)
    endif

    ! Now make sure everyone knows the answer
    call broadcast(ierr)
    restart_writable = (ierr == 0)

    if (.not. present(error_message)) return

    if (restart_writable) then
      error_message = ""
    else
      error_message = trim(io_err_msg)
    end if
  end function restart_writable

  !> FIXME : Add documentation
  subroutine init_gs2_save
  end subroutine init_gs2_save

  !> Sets the base of the restart file to be used when
  !> writing/reading files.
  subroutine set_restart_file (file)
    use constants, only: run_name_size
    use mp, only: proc0
    use file_utils, only: error_unit
    implicit none
    character(len=*), intent (in) :: file
    if (proc0 .and. len_trim(file) > run_name_size) then
       write( error_unit(), '("Argument to set_restart_file exceeds restart_file size so truncating")')
    end if
    restart_file = file(1:min(len_trim(file), run_name_size))
  end subroutine set_restart_file

  !> FIXME : Add documentation  
  subroutine finish_gs2_save
    implicit none
    call deallocate_arrays
  end subroutine finish_gs2_save

  !> Deallocate all module level arrays
  subroutine deallocate_arrays
#ifdef NETCDF    
    if (allocated(tmpr)) deallocate(tmpr)
    if (allocated(tmpi)) deallocate(tmpi)
    if (allocated(ftmpr)) deallocate(ftmpr)
    if (allocated(ftmpi)) deallocate(ftmpi)
    if (allocated(stmp)) deallocate(stmp)
    if (allocated(atmp)) deallocate(atmp)
#endif
  end subroutine deallocate_arrays

  !> Returns the file corresponding to restart file in current setup
  function get_file_proc(is_one_file_per_processor, fileopt) result(file_proc)
    use mp, only: iproc
    implicit none
    logical, intent(in) :: is_one_file_per_processor
    character(len=*), intent(in), optional :: fileopt
    character(run_name_size) :: file_proc
    file_proc = trim(restart_file)
    if (present(fileopt)) file_proc = trim(file_proc) // trim(adjustl(fileopt))
    if (is_one_file_per_processor) write(file_proc, '(A,".",I0)') trim(file_proc), iproc
  end function get_file_proc

  !> FIXME : Add documentation
  subroutine restore_current_scan_parameter_value(current_scan_parameter_value)
# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    integer :: istatus, current_scan_parameter_value_id_local
    integer :: ncid_local
    real, intent (out) :: current_scan_parameter_value
# ifdef NETCDF
    character (run_name_size) :: file_proc
    logical :: is_one_file_per_processor

    if (.not. include_parameter_scan) return

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       file_proc = get_file_proc(is_one_file_per_processor)

       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid_local)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus,file=file_proc)
       istatus = nf90_inq_varid (ncid_local, &
            "current_scan_parameter_value", &
            current_scan_parameter_value_id_local)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus,&
            var='current_scan_parameter_value_id')

       istatus = nf90_get_var (ncid_local, &
                                 current_scan_parameter_value_id_local, &
                                 current_scan_parameter_value)
       
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid_local,&
               current_scan_parameter_value_id_local)          
       endif
        
       istatus = nf90_close (ncid_local)
    endif

    call broadcast (current_scan_parameter_value)

# endif
  end subroutine restore_current_scan_parameter_value

  !> FIXME : Add documentation  
  subroutine init_dt (delt0, delt1, delt2, delt_max, istatus, not_set_value)
# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
    use optionals, only: get_option_with_default
# endif
    implicit none
    real, intent (in out) :: delt0, delt1, delt2, delt_max
    integer, intent (out) :: istatus
    real, intent(in), optional :: not_set_value
# ifdef NETCDF
    character (run_name_size) :: file_proc
    real :: not_set_value_to_use
    logical :: is_one_file_per_processor

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       file_proc = get_file_proc(is_one_file_per_processor)

       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus,file=file_proc)

       not_set_value_to_use = get_option_with_default(not_set_value, -1.0)

       ! Note unlike the explicit source terms all three time steps should always
       ! be available in the restart file so we don't silence the error messages here.
       ! The only situation we are likely to come across where the delt1 and delt2
       ! values aren't available is where we are trying to read an old restart file.
       ! This will then lead to the error/warning being displayed but the code should
       ! carry on as intended and the missing steps will be set to a special value to
       ! indicate they have not been set yet.
       istatus = nf90_inq_varid (ncid, "delt0", delt0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt0')
       
       istatus = nf90_get_var (ncid, delt0id, delt0)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt0id, message=' in init_dt')
          delt0 = not_set_value_to_use
       endif           

       istatus = nf90_inq_varid (ncid, "delt1", delt1id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt1')
       
       istatus = nf90_get_var (ncid, delt1id, delt1)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt1id, message=' in init_dt')
          delt1 = not_set_value_to_use
       endif           

       istatus = nf90_inq_varid (ncid, "delt2", delt2id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt2')
       
       istatus = nf90_get_var (ncid, delt2id, delt2)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt2id, message=' in init_dt')
          delt2 = not_set_value_to_use
       endif           

       istatus = nf90_inq_varid (ncid, "delt_max", delt_max_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt_max')

       istatus = nf90_get_var (ncid, delt_max_id, delt_max)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt_max_id, message=' in init_dt')
          delt_max = not_set_value_to_use
       endif

       istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (delt0)
    call broadcast (delt1)
    call broadcast (delt2)
    call broadcast (delt_max)

# endif
  end subroutine init_dt

!> FIXME : Add documentation  
  subroutine init_vnm (vnm, istatus)
# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, dimension(2), intent (in out) :: vnm
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (run_name_size) :: file_proc
    logical :: is_one_file_per_processor

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       file_proc = get_file_proc(is_one_file_per_processor)

       istatus = nf90_open (file_proc, 0, ncid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)

       istatus = nf90_inq_varid (ncid, "vnm1", vnm1id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='vnm1')

       istatus = nf90_inq_varid (ncid, "vnm2", vnm2id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='vnm2')

       istatus = nf90_get_var (ncid, vnm1id, vnm(1))

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, vnm1id, message=' in init_vnm')
          vnm(1) = 0.
       endif           

       istatus = nf90_get_var (ncid, vnm2id, vnm(2))

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, vnm2id, message=' in init_vnm')
          vnm(2) = 0.
       endif           

       istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (vnm)

# endif

  end subroutine init_vnm

  !> FIXME : Add documentation
  !!
  !! @note This routine gets a_ant and b_ant for proc 0 only!
  subroutine init_ant_amp (a_ant, b_ant, nk_stir, istatus)
# ifdef NETCDF
    use file_utils, only: error_unit
    use constants, only: zi
# endif
    use mp, only: proc0
    implicit none
    complex, dimension(:), intent (in out) :: a_ant, b_ant
    integer, intent (in) :: nk_stir
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (run_name_size) :: file_proc
    integer :: ierr, i
    logical :: is_one_file_per_processor

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       a_ant = 0. ; b_ant = 0.

       file_proc = get_file_proc(is_one_file_per_processor)

       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_open in init_ant_amp error: ", nf90_strerror(istatus) 
          write(ierr,*) "If you did not intend for this to be a restarted run with an external antenna,"
          write(ierr,*) "you may ignore the error message above."
          return
       endif

       istatus = nf90_inq_dimid (ncid, "nk_stir", nk_stir_dim)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='nk_stir')
       
       istatus = nf90_inquire_dimension (ncid, nk_stir_dim, len=i)
       if (istatus /= NF90_NOERR) &
            call netcdf_error (istatus, ncid, dimid=nk_stir_dim)
       if (i /= nk_stir) write(*,*) 'Restart error: nk_stir=? ',i,' : ',nk_stir

       istatus = nf90_inq_varid (ncid, "a_ant_r", a_antr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a_ant_r')

       istatus = nf90_inq_varid (ncid, "a_ant_i", a_anti_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a_ant_i')

       istatus = nf90_inq_varid (ncid, "b_ant_r", b_antr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='b_ant_r')

       istatus = nf90_inq_varid (ncid, "b_ant_i", b_anti_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='b_ant_i')

       if (.not. allocated(atmp)) allocate (atmp(nk_stir))
       atmp = 0.

       istatus = nf90_get_var (ncid, a_antr_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, a_antr_id)
       a_ant = atmp

       istatus = nf90_get_var (ncid, a_anti_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, a_anti_id)
       a_ant = a_ant + zi * atmp

       istatus = nf90_get_var (ncid, b_antr_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, b_antr_id)
       b_ant = atmp

       istatus = nf90_get_var (ncid, b_anti_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, b_anti_id)
       b_ant = b_ant + zi * atmp

       deallocate (atmp)
       istatus = nf90_close (ncid)
    endif

# else

    if (proc0) istatus = 2

# endif

  end subroutine init_ant_amp

  !> FIXME : Add documentation  
  subroutine read_t0_from_restart_file (tstart, istatus)
# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: tstart
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (run_name_size) :: file_proc
    logical :: is_one_file_per_processor

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       file_proc = get_file_proc(is_one_file_per_processor)
          
       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)
          
       istatus = nf90_inq_varid (ncid, "t0", t0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='t0')

       istatus = nf90_get_var (ncid, t0id, tstart)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, t0id, message=' in init_tstart')
          tstart = -1.
       end if

       istatus = nf90_close (ncid)

    end if

    call broadcast (istatus)
    call broadcast (tstart)

# endif

  end subroutine read_t0_from_restart_file
end module gs2_save
