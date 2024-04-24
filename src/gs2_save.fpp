!> FIXME : Add documentation
module gs2_save

  use mp, only: mp_comm, mp_info
  use constants, only: run_name_size

# ifdef NETCDF
# ifdef NETCDF_PARALLEL
  use netcdf, only: nf90_put_att, NF90_GLOBAL
# endif
  use netcdf, only: NF90_NOWRITE, NF90_NOERR
  use netcdf, only: nf90_open, nf90_close
  use netcdf, only: nf90_get_var, nf90_strerror
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  
  use netcdf_utils, only: check_netcdf_file_precision
  use netcdf_utils, only: netcdf_error
  use netcdf_utils, only: kind_nf
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

  !> Do we save the scan parameter or not?
  logical, parameter :: include_parameter_scan = .true.

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

  !> Names of dimensions for distribution-sized variables in netCDF files
  character(len=5), dimension(*), parameter :: g_dim_names = ["theta", "sign ", "glo  "]
  character(len=5), dimension(*), parameter :: field_dim_names = ["theta", "akx  ", "aky  "]
contains

  !> Create and fill the restart files
  subroutine gs2_save_for_restart &
       (g, t0, vnm, istatus, has_phi, has_apar, has_bpar, &
       code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, &
       fileopt, &
       save_glo_info_and_grids, save_velocities, header)
# ifdef NETCDF
    use fields_arrays, only: phinew, aparnew, bparnew
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use kt_grids, only: naky, ntheta0, jtwist_out, akx, aky
    use antenna_data, only: nk_stir, a_ant, b_ant, ant_on
    use netcdf_utils, only: ensure_netcdf_var_exists, ensure_netcdf_dim_exists
    use gs2_metadata, only: create_metadata
    use netcdf, only: nf90_collective, nf90_independent
    use neasyf, only: neasyf_open, neasyf_dim, neasyf_write, neasyf_close
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
    use optionals, only: get_option_with_default
    use standard_header, only : standard_header_type
    implicit none
    !> The distribution function \(g\)
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    !> Current simulation time
    real, intent (in) :: t0
    !> FIXME: Collisionality multiplier?
    real, dimension (2), intent (in) :: vnm
    !> The current time step, the two previous time steps and the
    !> maximum allowed time step.
    real, intent(in) :: code_dt, code_dt_prev1, code_dt_prev2, code_dt_max
    !> Do we include electrostatic potential, \(\phi\)
    logical, intent (in) :: has_phi
    !> Do we include vector potential, \(A_\parallel\)
    logical, intent (in) :: has_apar
    !> Do we include parallel magnetic field, \(B_\parallel\)
    logical, intent (in) :: has_bpar
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
    integer :: n_elements, ierr
    integer, dimension(3) :: start_pos, counts
    integer, parameter :: tmpunit = 348
    logical :: is_one_file_per_processor, this_proc_saves_fields
    logical :: has_gexp_1
    logical :: save_glo_info_and_grids_local, save_velocities_info
    logical :: should_define_single_entry_variables
    integer, parameter :: verb = 3
    integer :: iglo
    ! Whether to use collective (parallel) or independent access to netCDF variables
    integer :: par_access
    ! File handle
    integer :: ncid

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

    ! Decide if we want one file per processor or not
    is_one_file_per_processor = save_many .or. (.not. has_netcdf_parallel)

    ! Decide if we are saving the fields
    this_proc_saves_fields = (iproc == proc_to_save_fields) .or. (proc_to_save_fields < 0) &
         .or. .not. is_one_file_per_processor

    ! Open the file and define dimensions and variables
    file_proc = get_file_proc(is_one_file_per_processor, fileopt)

    call debug_message(verb, 'gs2_save::gs2_save_for_restart opening file')
    call debug_message(verb, 'file proc is')
    call debug_message(verb, trim(adjustl(file_proc)))

    if (is_one_file_per_processor) then
       ncid = neasyf_open(file_proc, "w")
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
       ncid = neasyf_open(file_proc, "w", comm=mp_comm, info=mp_info)
# endif          
    end if

    call create_metadata(ncid, "GS2 restart file", local_header)

    if (.not.is_one_file_per_processor) then
# ifdef NETCDF_PARALLEL
       istatus = nf90_put_att(ncid, NF90_GLOBAL, 'layout', layout)
# endif
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart defining dimensions')

    ! Work out the problem size
    if (is_one_file_per_processor) then
       n_elements = max(g_lo%ulim_proc-g_lo%llim_proc+1, 0)
    else
       n_elements = g_lo%ulim_world+1
    end if
    call neasyf_dim(ncid, "glo", dim_size=n_elements)

    call neasyf_dim(ncid, "theta", theta)
    call neasyf_dim(ncid, "sign", dim_size=2)
    call neasyf_dim(ncid, "aky", aky)
    call neasyf_dim(ncid, "akx", akx)

    ! Decide if this processor should take part in defining items
    ! which should only appear once in outputs.
    should_define_single_entry_variables = proc0 .or. (.not. is_one_file_per_processor)

    !For saving distribution function
    if (save_glo_info_and_grids_local .and. should_define_single_entry_variables) then
       call neasyf_dim(ncid, "lambda", values=al, long_name="Energy/magnetic moment", units="1/(2 B_a)")
       call neasyf_dim(ncid, "species", dim_size=nspec)

       ! Energy coordinate is a function of species, but the dimension must be 1D
       call neasyf_dim(ncid, "egrid", dim_size=negrid, long_name="See 'energy' for energy grid values")
       call neasyf_write(ncid, "energy", energy, dim_names=["egrid  ", "species"], &
            long_name="Energy grid values")
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart defining variables')

    if (ant_on) then
       call neasyf_dim(ncid, "nk_stir", dim_size=nk_stir)
    end if

    !For saving distribution function
    if (save_glo_info_and_grids_local) then
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
       endif
    end if

    call debug_message(verb, 'gs2_save::gs2_save_for_restart writing scalars')
    
    if (is_one_file_per_processor .or. proc0) then
       call neasyf_write(ncid, "delt0", code_dt)
       call neasyf_write(ncid, "delt1", code_dt_prev1)
       call neasyf_write(ncid, "delt2", code_dt_prev2)
       call neasyf_write(ncid, "delt_max", code_dt_max)
       call neasyf_write(ncid, "t0", t0)
       call neasyf_write(ncid, "layout", layout)
       call neasyf_write(ncid, "ntheta0", ntheta0)
       call neasyf_write(ncid, "naky", naky)
       call neasyf_write(ncid, "nlambda", nlambda)
       call neasyf_write(ncid, "negrid", negrid)
       call neasyf_write(ncid, "nspec", nspec)

       if (include_parameter_scan) then
          call neasyf_write(ncid, "current_scan_parameter_value", current_scan_parameter_value)
       end if

       call neasyf_write(ncid, "vnm1", vnm(1))
       call neasyf_write(ncid, "vnm2", vnm(2))
    endif

    if (ant_on) then

       if (.not. allocated(atmp)) allocate (atmp(nk_stir))
       atmp = real(a_ant)
       call neasyf_write(ncid, "a_ant_r", atmp, dim_names=["nk_stir"])
       atmp = aimag(a_ant)
       call neasyf_write(ncid, "a_ant_i", atmp, dim_names=["nk_stir"])

       atmp = real(b_ant)
       call neasyf_write(ncid, "b_ant_r", atmp, dim_names=["nk_stir"])
       atmp = aimag(b_ant)
       call neasyf_write(ncid, "b_ant_i", atmp, dim_names=["nk_stir"])
       deallocate (atmp)
    end if

    if (.not. allocated(tmpr)) &
         allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    tmpr = real(g)
    call debug_message(verb, 'gs2_save::gs2_save_for_restart writing dist fn')

    if (is_one_file_per_processor) then
       start_pos = [1, 1, 1]
       par_access = nf90_independent
    else
       start_pos = [1, 1, g_lo%llim_proc+1]
       par_access = nf90_collective
    end if
    counts = shape(tmpr)

    call neasyf_write(ncid, "gr", tmpr, dim_names=g_dim_names, &
         start=start_pos, count=counts, par_access=par_access)

    tmpr = aimag(g)

    call neasyf_write(ncid, "gi", tmpr, dim_names=g_dim_names, &
         start=start_pos, count=counts, par_access=par_access)

    if (include_explicit_source_in_restart) then
       !Explicit source term gexp_1
#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          tmpr = real(gexp_1)
          call neasyf_write(ncid, "gexp_1_r", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)
          
          tmpr = aimag(gexp_1)
          call neasyf_write(ncid, "gexp_1_i", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)
       end if

       !Explicit source term gexp_2
       if(allocated(gexp_2)) then
          tmpr = real(gexp_2)
          call neasyf_write(ncid, "gexp_2_r", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)
          
          tmpr = aimag(gexp_2)
          call neasyf_write(ncid, "gexp_2_i", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)
       end if

       !Explicit source term gexp_3
       if(allocated(gexp_3)) then
          tmpr = real(gexp_3)
          call neasyf_write(ncid, "gexp_3_r", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)

          tmpr = aimag(gexp_3)
          call neasyf_write(ncid, "gexp_3_i", tmpr, dim_names=g_dim_names, &
                  start=start_pos, count=counts, par_access=par_access)
       end if
    end if

    ! For saving distribution function
    if (save_glo_info_and_grids_local) then
       ! Fill processor and layout information for plotting purposes
       if (proc0) then
          call neasyf_write(ncid, "iglo_kykxsel", iglo_kykxsel, &
               dim_names=[character(len=7)::"aky", "akx", "species", "egrid", "lambda"])
          call neasyf_write(ncid, "procid_kykxsel", procid_kykxsel, &
               dim_names=[character(len=7)::"aky", "akx", "species", "egrid", "lambda"])
       end if

       call neasyf_write(ncid, "iglo_llim_proc", g_lo%llim_proc)

       if (proc0) then ! only save grids on proc0
          call neasyf_write(ncid, "jtwist", jtwist_out)
       end if
    end if

    if (save_velocities_info) then
       call neasyf_write(ncid, "vpa", vpa, dim_names=["theta", "sign ", "glo  "])
       call neasyf_write(ncid, "vperp2", vperp2, dim_names=["theta", "glo  "])
    end if

    if (is_one_file_per_processor .or. proc0) then
       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))

       if (this_proc_saves_fields) then
          if (has_phi) then
             ftmpr = real(phinew)
             call neasyf_write(ncid, "phi_r", ftmpr, dim_names=field_dim_names)
             
             ftmpr = aimag(phinew)
             call neasyf_write(ncid, "phi_i", ftmpr, dim_names=field_dim_names)
          end if

          if (has_apar) then
             ftmpr = real(aparnew)
             call neasyf_write(ncid, "apar_r", ftmpr, dim_names=field_dim_names)
             
             ftmpr = aimag(aparnew)
             call neasyf_write(ncid, "apar_i", ftmpr, dim_names=field_dim_names)
          end if

          if (has_bpar) then
             ftmpr = real(bparnew)
             call neasyf_write(ncid, "bpar_r", ftmpr, dim_names=field_dim_names)
             
             ftmpr = aimag(bparnew)
             call neasyf_write(ncid, "bpar_i", ftmpr, dim_names=field_dim_names)
          end if
       end if

       if (.not. allocated(stmp)) allocate (stmp(naky))

       ! Because we want to be able to restart using exb shear from a
       ! case which does not have exb shear we always add kx_shift and
       ! theta0_shift even if no exb shear present in simulation)
       if (allocated(kx_shift)) then
          stmp = kx_shift
       else
          stmp = 0.
       end if
       call neasyf_write(ncid, "kx_shift", stmp, dim_names=["aky"])
       
       if (allocated(theta0_shift)) then
          stmp = theta0_shift
       else
          stmp = 0.
       end if
       call neasyf_write(ncid, "theta0_shift", stmp, dim_names=["aky"])
    end if

    ! Close the file
    call neasyf_close(ncid)

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
  subroutine gs2_save_slice(g_in, it, ik, il, ie, is, time, nout, fileopt, header)
# ifdef NETCDF
    use neasyf, only: neasyf_open, neasyf_dim, neasyf_write, neasyf_close
    use gs2_metadata, only: create_metadata
    use unit_tests, only: debug_message
    use mp, only: iproc, proc0
    use optionals, only: get_option_with_default
    use convert, only: c2r
    use kt_grids, only: aky, akx, naky, ntheta0
    use le_grids, only: energy, al, negrid, nlambda
    use species, only: nspec
    use theta_grid, only: theta
    use gs2_layouts, only: proc_id, idx
# endif
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use standard_header, only : standard_header_type
    implicit none
    !> The distribution function \(g\)
    complex, dimension (-ntgrid:, :, g_lo%llim_proc:), intent (in) :: g_in
    !> The indices of the slice to save
    integer, intent(in) :: it, ik, il, ie, is
    !> The current time
    real, intent(in) :: time
    !> The current time output index
    integer, intent(in), optional :: nout
    !> Optional extra filename infix
    character (len=*), intent (in), optional :: fileopt
    !> Header for files with build and run information
    type(standard_header_type), optional, intent(in) :: header
# ifdef NETCDF
    character (run_name_size) :: file_proc
    integer :: the_iglo, the_proc
    integer, parameter :: verb = 3
    type(standard_header_type) :: local_header
    character (len=*), parameter :: file_opt_base = '.slice'
    logical :: appending, previous_file_exists
    real, dimension(:, :, :), allocatable :: ri_gslice
    integer, save :: nout_local = 1
    logical :: error_exit
    ! File handle
    integer :: ncid

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

    if (appending) then
       ncid = neasyf_open(trim(file_proc), 'rw')
    else
       ncid = neasyf_open(trim(file_proc), 'w')
       call create_metadata(ncid, "GS2 slice file", local_header)

       call debug_message(verb, 'gs2_save::gs2_save_slice defining dimensions')
       call neasyf_dim(ncid, "ri", dim_size = 2, long_name = "real/imaginary components")
       call neasyf_dim(ncid, "theta", values = theta, long_name="Poloidal angle", units="rad")
       call neasyf_dim(ncid, "sign",  dim_size = 2, long_name = "sign of v||")
       call neasyf_dim(ncid, "t", unlimited = .true., long_name = "Time", units = "L/vt")

       call debug_message(verb, 'gs2_save::gs2_save_slice defining constant variables')
       call neasyf_write(ncid, "it", it)
       call neasyf_write(ncid, "ik", ik)
       call neasyf_write(ncid, "il", il)
       call neasyf_write(ncid, "ie", ie)
       call neasyf_write(ncid, "is", is)

       call neasyf_write(ncid, "iglo", the_iglo)
       call neasyf_write(ncid, "iproc", the_proc)

       call neasyf_write(ncid, "kx", akx(it))
       call neasyf_write(ncid, "ky", aky(ik))
       call neasyf_write(ncid, "lambda", al(il))
       call neasyf_write(ncid, "energy", energy(ie, is))
    end if

    call neasyf_write(ncid, "t", time, start = [nout_local])

    ! Convert from complex to real
    allocate(ri_gslice(2, size(theta), 2))
    call c2r(g_in(:, :, the_iglo), ri_gslice)
    call neasyf_write(ncid, "g_slice", ri_gslice, start = [1, 1, 1, nout_local])

    ! Close the file
    call neasyf_close(ncid)

    nout_local = nout_local + 1
# endif
  end subroutine gs2_save_slice

  !> FIXME : Add documentation
  subroutine gs2_restore_many (g, scale, istatus, has_phi, has_apar, has_bpar, fileopt)
!MR, 2007: restore kx_shift array if already allocated
# ifdef NETCDF
    use mp, only: iproc, mp_abort, proc0, broadcast
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3, kx_shift, theta0_shift
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: layout
    use neasyf, only: neasyf_open, neasyf_read, neasyf_close
# endif
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use file_utils, only: error_unit
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    logical :: this_proc_reads_fields, need_to_broadcast_fields
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    logical, intent (in) :: has_phi, has_apar, has_bpar
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
    integer :: i
    real :: fac
    logical :: has_explicit_source_terms
    logical :: is_one_file_per_processor
    logical :: has_gexp_1
    logical, save :: explicit_warning_given = .false.
    ! File handle
    integer :: ncid
    ! Variable handles
    integer (kind_nf) :: thetaid, signid, gloid, kyid, kxid
    integer (kind_nf) :: gexp_1_r_id, gexp_1_i_id
    integer (kind_nf) :: gexp_2_r_id, gexp_2_i_id
    integer (kind_nf) :: gexp_3_r_id, gexp_3_i_id

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    ! Decide if we are reading the fields
    need_to_broadcast_fields = proc_to_save_fields >= 0
    this_proc_reads_fields = (iproc == proc_to_save_fields) .or. (.not. need_to_broadcast_fields)
    
    file_proc = get_file_proc(is_one_file_per_processor, fileopt)

    if (is_one_file_per_processor) then
       ncid = neasyf_open(file_proc, "r")
    else
       ncid = neasyf_open(file_proc, "r", comm=mp_comm, info=mp_info)
    end if

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

    call neasyf_read(ncid, "layout", restart_layout)

    ! Abort if the layouts don't match
    if(restart_layout /= layout) then
       call mp_abort("Incompatible layouts in restart file ("//restart_layout//") and simulation ("//layout//")",.true.)
    end if

    if (.not. allocated(tmpr)) &
         allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))
    if (.not. allocated(tmpi)) &
         allocate (tmpi(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    call zero_array(tmpr); call zero_array(tmpi)

    if (is_one_file_per_processor) then
       start_pos = [1, 1, 1]
    else
       start_pos = [1, 1, g_lo%llim_proc+1]
    end if
    counts = shape(tmpr)

    call neasyf_read(ncid, "gr", tmpr, start=start_pos, count=counts)
    call neasyf_read(ncid, "gi", tmpi, start=start_pos, count=counts)
    g = cmplx(tmpr, tmpi)

    if (include_explicit_source_in_restart .and. has_explicit_source_terms) then
       ! Explicit source term
#ifndef SHMEM
       has_gexp_1 = allocated(gexp_1)
#else
       has_gexp_1 = associated(gexp_1)
#endif

       if (has_gexp_1) then
          call neasyf_read(ncid, "gexp_1_r", tmpr, start=start_pos, count=counts)
          call neasyf_read(ncid, "gexp_1_i", tmpi, start=start_pos, count=counts)
          gexp_1 = cmplx(tmpr, tmpi)
       endif

       ! Explicit source term
       if(allocated(gexp_2)) then
          call neasyf_read(ncid, "gexp_2_r", tmpr, start=start_pos, count=counts)
          call neasyf_read(ncid, "gexp_2_i", tmpi, start=start_pos, count=counts)
          gexp_2 = cmplx(tmpr, tmpi)
       endif

       ! Explicit source term
       if(allocated(gexp_3)) then
          call neasyf_read(ncid, "gexp_3_r", tmpr, start=start_pos, count=counts)
          call neasyf_read(ncid, "gexp_3_i", tmpi, start=start_pos, count=counts)
          gexp_3 = cmplx(tmpr, tmpi)
       endif
    end if

    if (allocated(kx_shift)) then   ! MR begin
       if (.not. allocated(stmp)) allocate (stmp(naky))   ! MR 
       call neasyf_read(ncid, "kx_shift", stmp)
       kx_shift = stmp
    endif   ! MR end

    if (allocated(theta0_shift)) then
       if (.not. allocated(stmp)) allocate (stmp(naky))
       call neasyf_read(ncid, "theta0_shift", stmp)
       theta0_shift = stmp
    endif

    if (this_proc_reads_fields) then
       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
       if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

       if (has_phi) then
          call neasyf_read(ncid, "phi_r", ftmpr)
          call neasyf_read(ncid, "phi_i", ftmpi)
          phinew = cmplx(ftmpr, ftmpi)
       end if

       if (has_apar) then
          call neasyf_read(ncid, "apar_r", ftmpr)
          call neasyf_read(ncid, "apar_i", ftmpi)
          aparnew = cmplx(ftmpr, ftmpi)
       end if

       if (has_bpar) then
          call neasyf_read(ncid, "bpar_r", ftmpr)
          call neasyf_read(ncid, "bpar_i", ftmpi)
          bparnew = cmplx(ftmpr, ftmpi)
       end if
    end if

    if (need_to_broadcast_fields) then
       if (has_phi) call broadcast(phinew, proc_to_save_fields)
       if (has_apar) call broadcast(aparnew, proc_to_save_fields)
       if (has_bpar) call broadcast(bparnew, proc_to_save_fields)
    end if

    if (has_phi) call zero_array(phi)
    if (has_apar) call zero_array(apar)
    if (has_bpar) call zero_array(bpar)

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

    call neasyf_close(ncid)

    ! Deallocate modul level arrays
    call deallocate_arrays
    
# else
    write (error_unit(),*) &
         'ERROR: gs2_restore_many is called without netcdf'
# endif
  end subroutine gs2_restore_many

  !> This routine writes a passed square complex array to a file
  !> with passed name
  subroutine gs2_save_response(resp, fname, code_dt, condition, header)
    use file_utils, only: error_unit
    use standard_header, only: standard_header_type
    use optionals, only: get_option_with_default
#ifdef NETCDF
    use convert, only: c2r
    use netcdf_utils, only: ensure_netcdf_dim_exists, ensure_netcdf_var_exists
    use gs2_metadata, only: create_metadata
    use neasyf, only: neasyf_open, neasyf_dim, neasyf_write, neasyf_close
#else
    use file_utils, only: get_unused_unit
#endif
    implicit none
    complex,dimension(:,:), intent(in) :: resp
    character(len=*), intent(in) :: fname
    real, intent(in) :: code_dt
    real, intent(in), optional :: condition
    !> Header for files with build and run information
    type(standard_header_type), optional, intent(in) :: header

    integer :: sz
#ifdef NETCDF
    integer :: ncid
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

    ncid = neasyf_open(fname, "w")
    call create_metadata(ncid, "GS2 response matrix", local_header)

    call neasyf_dim(ncid, "ri", dim_size=2)
    call neasyf_dim(ncid, "ax1", dim_size=sz)
    call neasyf_dim(ncid, "ax2", dim_size=sz)

    call neasyf_write(ncid, "dt", code_dt)
    call neasyf_write(ncid, "condition", get_option_with_default(condition, -1.0))

    !/Convert complex to ri and write
    allocate(ri_resp(2,sz,sz))
    call c2r(resp,ri_resp)
    call neasyf_write(ncid, "response", ri_resp, dim_names=["ri ", "ax1", "ax2"])
    deallocate(ri_resp)

    !/Now close the file
    call neasyf_close(ncid)

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
  subroutine gs2_restore_response(resp, fname, code_dt, condition)
    use file_utils, only: error_unit
#ifdef NETCDF
    use neasyf, only: neasyf_open, neasyf_close, neasyf_read
    use convert, only: r2c
#else
    use file_utils, only: get_unused_unit
#endif
    implicit none
    complex, dimension(:,:), intent(out) :: resp
    character(len=*), intent(in) :: fname
    real, intent(out) :: code_dt
    real, intent(out), optional :: condition
    integer :: sz
#ifdef NETCDF
    integer :: ncid
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
    ncid = neasyf_open(fname, "r")

    !/Read and convert ri to complex
    allocate(ri_resp(2,sz,sz))
    call neasyf_read(ncid, "response", ri_resp)
    call r2c(resp,ri_resp)
    deallocate(ri_resp)

    !/Get the time step
    call neasyf_read(ncid, "dt", code_dt)

    !/Get the condition number
    if (present(condition)) then
       call neasyf_read(ncid, "condition", condition)
    end if

    !/Now close the file
    call neasyf_close(ncid)
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
  subroutine init_eigenfunc_file(fname, IDs)
#ifdef NETCDF
    use theta_grid, only: theta
    use netcdf_utils, only: ensure_netcdf_dim_exists, ensure_netcdf_var_exists
    use neasyf, only: neasyf_open, neasyf_dim
#endif
    implicit none
    character(len=*), intent(in) :: fname
    type(EigNetcdfID), intent(inout) :: IDs

    !Set nconv counter to 0
    IDs%nconv_count=0

#ifdef NETCDF
    IDs%ncid = neasyf_open(fname, "w")

    call neasyf_dim(IDs%ncid, "ri", dim_size=2)
    call neasyf_dim(IDs%ncid, "theta", theta)
    call neasyf_dim(IDs%ncid, "conv", unlimited=.true.)
#endif
  end subroutine init_eigenfunc_file

  !> Add an eigenpairs data to file
  subroutine add_eigenpair_to_file(eval, has_phi, has_apar, has_bpar,IDs,my_conv)
    use fields_arrays, only: phinew, aparnew, bparnew
    use convert, only: c2r
    use theta_grid, only: ntgrid
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    complex, intent(in) :: eval !Note just use fields to get eigenvectors
    real, intent(in), optional :: my_conv
    type(EigNetcdfID), intent(inout) :: IDs
    logical, intent(in) :: has_phi, has_apar, has_bpar
#ifdef NETCDF
    real, dimension(2) :: ri_omega
    real, dimension(:,:), allocatable :: ri_field
    integer, dimension(3) :: start_field, count_field
    real :: local_conv
    character(*), dimension(*), parameter :: eigen_field_dim_names = [character(len=5)::"ri", "theta", "conv"]
#endif

    !First increment counter
    IDs%nconv_count=IDs%nconv_count+1

#ifdef NETCDF
    start_field = [1, 1, IDs%nconv_count]
    count_field = [2, 2*ntgrid+1, 1]

    !Set the conv value
    if(present(my_conv))then
       local_conv=my_conv
    else
       local_conv=IDs%nconv_count*1.0
    endif

    !Now we can write data
    !/Conv
    call neasyf_write(IDs%ncid, "conv", local_conv, dim_names=["conv"], &
         start=[IDs%nconv_count])

    !/Omega
    ri_omega(1)=real(eval)
    ri_omega(2)=aimag(eval)
    call neasyf_write(IDs%ncid, "omega", ri_omega, dim_names=["ri  ", "conv"], &
         start=[1, IDs%nconv_count], count=[2,1])

    !/Fields
    allocate(ri_field(2,2*ntgrid+1))
    if(has_phi)then
       call c2r(phinew(:,1,1),ri_field)
       call neasyf_write(IDs%ncid, "phi", ri_field, dim_names=eigen_field_dim_names, &
            start=start_field, count=count_field)
    endif
    if(has_apar)then
       call c2r(aparnew(:,1,1),ri_field)
       call neasyf_write(IDs%ncid, "apar", ri_field, dim_names=eigen_field_dim_names, &
            start=start_field, count=count_field)
    endif
    if(has_bpar)then
       call c2r(bparnew(:,1,1),ri_field)
       call neasyf_write(IDs%ncid, "bpar", ri_field, dim_names=eigen_field_dim_names, &
            start=start_field, count=count_field)
    endif
    deallocate(ri_field)       
#endif
  end subroutine add_eigenpair_to_file

  !> Close the eigenfunction file
  subroutine finish_eigenfunc_file(IDs)
#ifdef NETCDF
    use neasyf, only: neasyf_close
#endif
    implicit none
    type(EigNetcdfID), intent(inout) :: IDs
#ifdef NETCDF
    call neasyf_close(IDs%ncid)
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
    use neasyf, only: neasyf_open, neasyf_read, neasyf_close
# endif
    implicit none
    integer :: ncid_local
    real, intent (out) :: current_scan_parameter_value
# ifdef NETCDF
    character (run_name_size) :: file_proc
    logical :: is_one_file_per_processor

    if (.not. include_parameter_scan) return

    is_one_file_per_processor = read_many .or. (.not. has_netcdf_parallel)

    if (proc0) then
       file_proc = get_file_proc(is_one_file_per_processor)

       ncid_local = neasyf_open(file_proc, "r")
       call neasyf_read(ncid_local, "current_scan_parameter_value", current_scan_parameter_value)
       call neasyf_close(ncid_local)
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
    use neasyf, only: neasyf_open, neasyf_read, neasyf_close
# endif
    implicit none
    real, intent (in out) :: delt0, delt1, delt2, delt_max
    integer, intent (out) :: istatus
    real, intent(in), optional :: not_set_value
# ifdef NETCDF
    character (run_name_size) :: file_proc
    real :: not_set_value_to_use
    logical :: is_one_file_per_processor
    ! NetCDF handles
    integer(kind_nf) :: ncid, delt0id, delt1id, delt2id, delt_max_id

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
    ! NetCDF handles
    integer(kind_nf) :: ncid, vnm1id, vnm2id

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
    use neasyf, only: neasyf_read, neasyf_close
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
    ! NetCDF handles
    integer(kind_nf) :: ncid, nk_stir_dim

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

       if (.not. allocated(atmp)) allocate (atmp(nk_stir))
       atmp = 0.

       call neasyf_read(ncid, "a_ant_r", atmp)
       a_ant = atmp

       call neasyf_read(ncid, "a_ant_i", atmp)
       a_ant = a_ant + zi * atmp

       call neasyf_read(ncid, "b_ant_r", atmp)
       b_ant = atmp

       call neasyf_read(ncid, "b_ant_i", atmp)
       b_ant = b_ant + zi * atmp

       deallocate (atmp)
       call neasyf_close(ncid)
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
    ! NetCDF handles
    integer(kind_nf) :: ncid, t0id

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
