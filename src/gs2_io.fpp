!> A module for writing to netcdf.
module gs2_io

# ifdef NETCDF
  use netcdf, only: NF90_NOERR
  use netcdf_utils, only: netcdf_error, kind_nf
# endif

  implicit none

  private

  public :: init_gs2_io, nc_eigenfunc, nc_final_fields, nc_final_epar
  public :: nc_final_moments, nc_final_an, nc_finish
  public :: nc_qflux, nc_vflux, nc_pflux, nc_pflux_tormom, nc_loop, nc_loop_moments
  public :: nc_loop_movie, nc_write_moments, nc_exchange
  public :: nc_loop_fullmom, nc_loop_dist, nc_loop_sym, nc_loop_corr, nc_loop_corr_extend
  public :: nc_loop_partsym_tormom, nc_sync, get_netcdf_file_id

# ifdef NETCDF
  integer (kind_nf) :: ncid
  integer (kind_nf) :: naky_dim, nakx_dim, nttot_dim, negrid_dim, nlambda_dim, nspec_dim, ncoord_dim, ncoordt_dim
  integer (kind_nf) :: nsign_dim, time_dim, char5_dim, char10_dim, ri_dim, nheat_dim, vpa_dim
  integer (kind_nf) :: nttotext_dim, time_big_dim

  character(*), dimension(*), parameter :: final_field_dims = [character(len=5)::"ri", "theta", "kx", "ky"]

  integer, dimension (6) :: lp_dim
  integer, dimension (5) :: field_dim, final_mom_dim, heatk_dim, err_dim
  integer, dimension (5) :: phi_corr_dim
  integer, dimension (4) :: omega_dim, fluxk_dim, sym_dim
  integer, dimension (4) :: phi2_extend_dim, phi_corr_2pi_dim
  integer, dimension (3) :: fluxx_dim
  integer, dimension (3) :: mode_dim, heat_dim, flux_dist_dim, phi_dist_dim
  integer, dimension (2) :: kx_dim, ky_dim, om_dim, flux_dim, fmode_dim

  ! added by EAB 03/05/04 for movies
  logical :: my_make_movie, io_write_corr_extend
  integer :: ncid_movie
  integer :: nx_dim, ny_dim, nth_dim, time_movie_dim
  integer, dimension (4) :: xmode_dim
  integer :: nx_id, ny_id, nth_id, x_id, y_id, th_id, time_movie_id
  integer :: phi_by_xmode_id, apar_by_xmode_id, bpar_by_xmode_id

  integer :: phi2_id, apar2_id, bpar2_id
  integer :: current_scan_parameter_value_id
  integer :: phi2_by_mode_id, apar2_by_mode_id, bpar2_by_mode_id, errest_by_mode_id, lpcoef_by_mode_id
  integer :: phtot_id
  integer :: phi2_by_kx_id, apar2_by_kx_id, bpar2_by_kx_id
  integer :: phi2_by_ky_id, apar2_by_ky_id, bpar2_by_ky_id
  integer :: phi0_id, apar0_id, bpar0_id
  integer :: omega_id, omegaavg_id
  integer :: es_part_tormom_flux_id
  integer :: es_part_flux_dist_id, es_mom_flux_par_dist_id, es_mom_flux_perp_dist_id, es_heat_flux_dist_id, es_phi_dist_id
  integer :: es_part_tormom_by_k_id
  integer :: es_mom_sym_id, phi_corr_id, phi2_extend_id, phi_corr_2pi_id
  integer :: es_part_sym_id 
  integer :: phi_t_id, apar_t_id, bpar_t_id
  integer :: code_id
  integer :: h_energy_id, h_energy_dot_id, h_antenna_id
  integer :: h_eapar_id, h_ebpar_id
  integer :: h_delfs2_id, h_hs2_id, h_phis2_id
  integer :: h_hypervisc_id, h_hyperres_id, h_collisions_id
  integer :: h_gradients_id, h_hypercoll_id, h_heating_id
  integer :: h_imp_colls_id
  integer :: hk_energy_id, hk_energy_dot_id, hk_antenna_id
  integer :: hk_eapar_id, hk_ebpar_id
  integer :: hk_delfs2_id, hk_hs2_id, hk_phis2_id
  integer :: hk_hypervisc_id, hk_hyperres_id, hk_collisions_id
  integer :: hk_gradients_id, hk_hypercoll_id, hk_heating_id, hk_imp_colls_id
  logical :: write_apar_t, write_phi_t, write_bpar_t ! Should the fields be written out every nwrite?

  !> Write a `complex` array to netcdf
  !>
  !> Converts the `complex` array to a `real` array with an extra dimension
  interface netcdf_write_complex
    module procedure write_complex_rank1, write_complex_rank2
    module procedure write_complex_rank3, write_complex_rank4
  end interface netcdf_write_complex

#else
  integer, parameter :: ncid = 0
#endif
  real, parameter :: zero = epsilon(0.0)

contains
  
  !> Provides read-only access to the netCDF file ID
  integer function get_netcdf_file_id()
    get_netcdf_file_id = ncid
  end function get_netcdf_file_id

  !> Create and initialise the output netCDF files
  !>
  !> Create the main and movie (real-space fields) netCDF files,
  !> creating the dimensions and main variables
  subroutine init_gs2_io (append_old, write_nl_flux, write_omega, &
      write_heating, write_final_antot, write_eigenfunc, &
      make_movie, nmovie_tot, write_verr, &
      write_moments,&
      write_full_moments_notgc, write_sym, write_pflux_sym, write_nl_flux_dist, &
      write_pflux_tormom, write_correlation, nwrite_big_tot, write_correlation_extend, & 
      write_phi_over_time, write_apar_over_time, write_bpar_over_time, &
      write_avg_moments, write_final_moments, write_final_epar, nout, &
      ob_midplane, header)
    use mp, only: proc0
    use standard_header, only : standard_header_type
# ifdef NETCDF
    use mp, only: barrier
    use file_utils, only: run_name, error_unit
    use gs2_transforms, only: init_transforms
    use kt_grids, only: naky, ntheta0,nx,ny
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use netcdf, only: NF90_CLOBBER, nf90_create, nf90_open, NF90_WRITE, NF90_NETCDF4
    use netcdf_utils, only: get_netcdf_code_precision, netcdf_real
    use netcdf_utils, only: get_netcdf_dim_length
    use constants, only: run_name_size
    use gs2_metadata, only: create_metadata
    use run_parameters, only: user_comments
    use neasyf, only: neasyf_open
# endif
    implicit none
    !> See [[gs2_diagnostics_knobs]] for details of these inputs
    logical, intent(in) :: append_old
    logical, intent(in) :: write_nl_flux, write_nl_flux_dist, write_omega
    logical, intent(in) :: write_heating, make_movie, write_moments
    logical, intent(in) :: write_final_antot, write_eigenfunc, write_verr
    logical, intent(in) :: write_full_moments_notgc, write_sym, write_correlation
    logical, intent(in) :: write_pflux_sym, write_pflux_tormom 
    logical, intent(inout) :: write_correlation_extend
    integer, intent (in) :: nmovie_tot, nwrite_big_tot
    logical, intent(in) :: write_phi_over_time
    logical, intent(in) :: write_apar_over_time,  write_bpar_over_time
    logical, intent(in) :: write_avg_moments, write_final_moments, write_final_epar
    !> If true, write out moments at the outboard midplane only,
    !> otherwise write them out along the entire flux surface
    logical, intent(in), optional :: ob_midplane
    integer, intent(out) :: nout
    !> Header for files with build and run information
    type(standard_header_type), intent(in) :: header

# ifdef NETCDF
    logical :: accelerated
    character(run_name_size) :: filename
    integer :: status
    logical :: previous_file_exists, appending
    integer, parameter :: create_file_mode = ior(NF90_CLOBBER, NF90_NETCDF4)

    !################
    !# Update internal flags

    write_phi_t = write_phi_over_time
    write_apar_t = write_apar_over_time
    write_bpar_t = write_bpar_over_time

    my_make_movie = make_movie
    io_write_corr_extend = write_correlation_extend

    !----------------

    status = NF90_NOERR

    !Create file
    filename = trim(trim(run_name)//'.out.nc')

    ! only proc0 opens the file:
    if (proc0) then
       ! See if the output file already exists, so we know how to open
       ! the file.  This would be simpler if netCDF had a "open if it
       ! exists, create if it doesn't" mode, but it doesn't so we have
       ! to do this instead
       inquire(file=trim(filename), exist=previous_file_exists)
       if (append_old .and. .not. previous_file_exists) then
         write(error_unit(), '(a, a, a)') "WARNING: 'append_old' requested, but '", trim(filename), &
              "' does not exist. Creating new file"
       end if
       appending = (previous_file_exists .and. append_old)

       if (appending) then
         ncid = neasyf_open(trim(filename), "rw")
       else
         ncid = neasyf_open(trim(filename), "w")
       end if

       call create_metadata(ncid, "GS2 Simulation Data", header, user_comments)
    end if

    !Create the movie file
    if (my_make_movie) then
       !perhaps put a call to init_y_transform
       !    call init_y_transform_layouts  
       call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

       ! only proc0 opens the file:
       if (proc0) then
         ncid_movie = neasyf_open(trim(run_name)//'.movie.nc', "w")
         call create_metadata(ncid_movie, 'GS2 Simulation x,y,theta Data for Movies', header, user_comments)
       end if
    endif

    if (proc0) then
        !Setup main file by creating netCDF file
        call define_dims (ncid, nwrite_big_tot, io_write_corr_extend)
        call define_vars (write_nl_flux, write_omega, &
                          write_heating, write_final_antot, write_eigenfunc, write_verr, &
                          write_moments, write_full_moments_notgc, write_sym, write_pflux_sym, &
                          write_nl_flux_dist, write_pflux_tormom, write_correlation, &
                          write_correlation_extend, write_avg_moments, write_final_moments, write_final_epar, &
                          ob_midplane=ob_midplane)

        if (appending) then
            nout = get_dim_length(ncid, "t")
        else
            !Write constant data
            call nc_norms(ncid)
            call nc_species(ncid)
            call nc_geo(ncid)
            call save_input(ncid)
            nout = 1
        end if

        if (my_make_movie) then
          !Setup movie file
          call define_dims_mymovie (nmovie_tot)
          call define_vars_mymovie

          !Write constant movie data
          call nc_grids_mymovie
        end if
    end if
# endif
  end subroutine init_gs2_io

  !> Ensure the netCDF file contains all the dimensions and grids,
  !> creating them if necessary
  !>
  !> Dimension variable IDs are module level variables
  subroutine define_dims (file_id, nwrite_big_tot, create_extended_correlation)
# ifdef NETCDF
    use file_utils, only: num_input_lines
    use kt_grids, only: aky, akx, theta0, jtwist_out, ntheta0, naky, nx, ny
    use theta_grid, only: theta, ntgrid
    use le_grids, only: al, energy, negrid, nlambda
    use species, only: nspec
    use mp, only: nproc
    use nonlinear_terms, only: nonlin
    use gs2_layouts, only: layout
    use neasyf, only: neasyf_dim, neasyf_write
    use netcdf_utils, only: ensure_netcdf_dim_exists
# endif
    implicit none
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Total number of output timesteps for "large" diagnostics
    integer, intent(in) :: nwrite_big_tot
    !> Create the dimensions for the extended correlation diagnostics
    logical, intent(in) :: create_extended_correlation
# ifdef NETCDF
    real :: nmesh

    call neasyf_dim(file_id, "kx", values=akx, dimid=nakx_dim, long_name="kx rho")
    call neasyf_dim(file_id, "ky", values=aky, dimid=naky_dim, long_name="ky rho")
    call neasyf_dim(file_id, "theta", values=theta, dimid=nttot_dim)
    call neasyf_dim(file_id, "lambda", values=al, dimid=nlambda_dim)
    call neasyf_dim(file_id, "species", dim_size=nspec, dimid=nspec_dim)

    ! Energy is 2D, and neasyf only handles 1D dimensions
    call ensure_netcdf_dim_exists(ncid, 'energy', negrid, negrid_dim)
    call neasyf_write(file_id, "energy", energy, dim_names=["energy ", "species"])

    call neasyf_dim(file_id, "sign", dim_size=2, dimid=nsign_dim)
    call neasyf_dim(file_id, "t", unlimited=.true., dimid=time_dim, &
         long_name="Time", units="L/vt")
    call neasyf_dim(file_id, "char5", dim_size=5, dimid=char5_dim)
    call neasyf_dim(file_id, "char10", dim_size=10, dimid=char10_dim)
    call neasyf_dim(file_id, "char200", dim_size=200)
    call neasyf_dim(file_id, "nlines", dim_size=num_input_lines, long_name="Input file line number")
    call neasyf_dim(file_id, "ri", dim_size=2, dimid=ri_dim, &
         long_name="Complex components", units="(real, imaginary)")
    call neasyf_dim(file_id, "nheat", dim_size=7, dimid=nheat_dim)
    call neasyf_dim(file_id, "vpa", dim_size=negrid*nlambda, dimid=vpa_dim)

    if (create_extended_correlation) then
       call neasyf_dim(file_id, "theta_ext", dim_size=(2*ntgrid+1)*((ntheta0-1)/jtwist_out+1), dimid=nttotext_dim)
       call neasyf_dim(file_id, "tb", dim_size=nwrite_big_tot, dimid=time_big_dim)
    end if

    call neasyf_write(file_id, "layout", layout)
    call neasyf_write(file_id, "nproc", nproc, long_name="Number of processors")
    call neasyf_write(file_id, "ntheta_tot", 2*ntgrid+1, long_name="Total number of theta points")
    call neasyf_write(file_id, "nkx", ntheta0, long_name="Number of kx points")
    call neasyf_write(file_id, "nky", naky, long_name="Number of ky points")
    call neasyf_write(file_id, "nspecies", nspec, long_name="Number of species")
    call neasyf_write(file_id, "nenergy", negrid, long_name="Number of energy points")
    call neasyf_write(file_id, "nlambda", nlambda, long_name="Number of lambda points")

    call neasyf_write(file_id, "theta0", theta0, dim_names=["kx", "ky"], long_name="Theta_0")

    if (nonlin) then
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
    else
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
    end if

    call neasyf_write(file_id, "nmesh", nmesh, long_name="Total number of mesh points")

# endif
  end subroutine define_dims

  !> Define the dimensions used in the movie file. Some duplication with
  !> main file -- may want to merge movie data into main file
  !> Main reason for using a new file is we can only have one unlimited dimension
  !> per file (with new versions of NETCDF this is no longer true) and we
  !> probably don't want to write movie data as often as other data.
  !> However, at the moment we precalculate the number of outputs, so it
  !> isn't actually unlimited
  subroutine define_dims_mymovie(nmovie_tot)
# ifdef NETCDF
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use netcdf, only: nf90_def_dim
    use netcdf_utils, only: ensure_netcdf_dim_exists
# endif
    implicit none
    integer, intent(in) :: nmovie_tot
# ifdef NETCDF
    integer :: nx, ny

    ! added by EAB 03/05/04 for GS2 movies
    ! <doc> If make_movie = .true., define some extra dimensions for movie output; x with dimension yxf_lo%nx; y with dimension yxf_lo%ny, and tm, (i.e. time) with dimension nmovie_tot = nstep/nmovie </doc>
    if (my_make_movie) then
       nx=yxf_lo%nx
       ny=yxf_lo%ny
       call ensure_netcdf_dim_exists(ncid_movie, 'x', nx, nx_dim)
       call ensure_netcdf_dim_exists(ncid_movie, 'y', ny, ny_dim)
       call ensure_netcdf_dim_exists(ncid_movie, 'theta', 2*ntgrid+1, nth_dim)
       call ensure_netcdf_dim_exists(ncid_movie, 'tm', nmovie_tot, time_movie_dim)
    endif
# endif
  end subroutine define_dims_mymovie

  !> Write the physical normalisations to the output netCDF file
  !> Write the physical normalisations to the output netCDF file
  subroutine nc_norms(file_id)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use normalisations, only: norms
#endif
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id

#ifdef NETCDF
    call neasyf_write(file_id, "mref", norms%get_value("mref"), &
         long_name="Reference mass in atomic mass units ", units="a.m.u.")
    call neasyf_write(file_id, "zref", norms%get_value("zref"), &
         long_name="Reference charge", units="proton charge")
    call neasyf_write(file_id, "nref", norms%get_value("nref"), &
         long_name="The reference density ", units="m^-3")
    call neasyf_write(file_id, "tref", norms%get_value("tref"), &
         long_name="The reference temperature ", units="eV")
    call neasyf_write(file_id, "aref", norms%get_value("aref"), &
         long_name="The reference length ", units="m")
    call neasyf_write(file_id, "vref", norms%get_value("vref"), &
         long_name="The reference (thermal) velocity ", units="m/s")
    call neasyf_write(file_id, "bref", norms%get_value("bref"), &
         long_name="The reference magnetic field ", units="Tesla")
    call neasyf_write(file_id, "rhoref", norms%get_value("rhoref"), &
         long_name="The reference larmour radius ", units="m")
#endif
  end subroutine nc_norms

  !> Puts the grid data into the movie file
  subroutine nc_grids_mymovie
#ifdef NETCDF
    use constants, only: pi
    use kt_grids, only: akx, aky
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid, theta
    use netcdf, only: nf90_put_var
#endif
    implicit none
# ifdef NETCDF
    real, dimension(:), allocatable :: x,y
    integer :: it, ik, status
    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then

       allocate(x(yxf_lo%nx),y(yxf_lo%ny))

       do it = 1, yxf_lo%nx 
          x(it) = 2.0*pi/akx(2)*(-0.5+ real(it-1)/real(yxf_lo%nx))
       end do
       do ik = 1, yxf_lo%ny
          y(ik) = 2.0*pi/aky(2)*(-0.5+ real(ik-1)/real(yxf_lo%ny))
       end do
       
       status = nf90_put_var (ncid_movie, nx_id, yxf_lo%nx)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, nx_id)
       status = nf90_put_var (ncid_movie, ny_id, yxf_lo%ny)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, ny_id)
       status = nf90_put_var (ncid_movie, nth_id, 2*ntgrid+1)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, nth_id)
       status = nf90_put_var (ncid_movie, x_id, x)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, x_id)
       status = nf90_put_var (ncid_movie, y_id, y)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, y_id)
       status = nf90_put_var (ncid_movie, th_id, theta)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, th_id)

       deallocate(x,y)
    endif
# endif
  end subroutine nc_grids_mymovie

  !> Close the output netCDF files
  subroutine nc_finish
# ifdef NETCDF
    use mp, only: proc0
    use netcdf, only: nf90_close
    use netcdf_utils, only: netcdf_error
    implicit none
    integer :: status

    if (proc0) then
       status = nf90_close (ncid)
       if (status /= NF90_NOERR) call netcdf_error (status)
       if(my_make_movie) then
          status = nf90_close (ncid_movie)
          if (status /= NF90_NOERR) call netcdf_error (status)
       endif
    end if
# endif
  end subroutine nc_finish

  !> Save the input file in the NetCDF file
  subroutine save_input(file_id)
#ifdef NETCDF
    use file_utils, only: num_input_lines, get_input_unit
    use neasyf, only: neasyf_write
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
#ifdef NETCDF
    integer, parameter :: line_length = 200
    character(line_length), dimension(:), allocatable ::  input_file_array
    integer :: n, unit

    allocate(input_file_array(num_input_lines))

    call get_input_unit(unit)
    rewind (unit=unit)
    do n = 1, num_input_lines
      read (unit=unit, fmt="(a)") input_file_array(n)
    end do

    call neasyf_write(file_id, "input_file", input_file_array, &
         long_name="Input file", dim_names=["char200", "nlines "])
#endif
  end subroutine save_input

  !> Create (but don't fill) the various diagnostics variables in the
  !> output netCDF file
  !>
  !> This also calls `nf90_enddef` which leaves "define" mode, so no
  !> more variables can be added to the file after this call
  subroutine define_vars (write_nl_flux, write_omega, &
       write_heating, write_final_antot, write_eigenfunc, write_verr, &
       write_moments, write_full_moments_notgc, write_sym, write_pflux_sym, &
       write_nl_flux_dist, write_pflux_tormom, write_correlation, &
       write_correlation_extend, write_avg_moments, write_final_moments, write_final_epar, &
       ob_midplane)
# ifdef NETCDF
    use file_utils, only: error_unit
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use parameter_scan_arrays, only: write_scan_parameter
    use netcdf, only: NF90_CHAR, NF90_INT, NF90_GLOBAL
    use netcdf, only: nf90_enddef, nf90_inq_libvers, nf90_put_var
    use netcdf_utils, only: netcdf_real, get_netcdf_dim_length
    use netcdf_utils, only: ensure_netcdf_att_exists, ensure_netcdf_var_exists
# endif
    implicit none
    logical, intent(in) :: write_nl_flux, write_nl_flux_dist, write_omega
    logical, intent(in) :: write_heating, write_final_antot
    logical, intent(in) :: write_eigenfunc, write_verr, write_moments
    logical, intent(in) :: write_full_moments_notgc, write_sym, write_correlation
    logical, intent(in) :: write_pflux_sym, write_pflux_tormom
    logical, intent(inout) :: write_correlation_extend
    logical, intent(in) :: write_avg_moments, write_final_moments, write_final_epar
    logical, optional, intent(in) :: ob_midplane
# ifdef NETCDF
    character (20) :: datestamp, timestamp, timezone
    integer :: status, ierr, dimsize(5)
    real    :: total
    logical :: obmid=.false.

    !Handle optional input
    if (present(ob_midplane)) then
       obmid=ob_midplane
    endif

    !Define dimensions/shapes for different data types
    fmode_dim       = [nakx_dim, naky_dim]
    mode_dim        = [nakx_dim, naky_dim, time_dim]
    err_dim         = [nttot_dim, nakx_dim, naky_dim, ncoordt_dim, time_dim]
    lp_dim          = [nttot_dim, nakx_dim, naky_dim, nspec_dim, ncoord_dim, time_dim]
    kx_dim          = [nakx_dim, time_dim]
    ky_dim          = [naky_dim, time_dim]
    om_dim          = [ri_dim, time_dim]
    omega_dim       = [ri_dim, nakx_dim, naky_dim, time_dim]
    flux_dim        = [nspec_dim, time_dim]
    fluxk_dim       = [nakx_dim, naky_dim, nspec_dim, time_dim]
    fluxx_dim       = [nakx_dim, nspec_dim, time_dim]
    heat_dim        = [nspec_dim, nheat_dim, time_dim]
    heatk_dim       = [nakx_dim, naky_dim, nspec_dim, nheat_dim, time_dim]
    field_dim       = [ri_dim, nttot_dim, nakx_dim, naky_dim, time_dim]
    final_mom_dim   = [ri_dim, nttot_dim, nakx_dim, naky_dim, nspec_dim]
    flux_dist_dim   = [nttot_dim, nspec_dim, time_dim]
    phi_dist_dim    = [ri_dim, nttot_dim, time_dim]
    sym_dim         = [nttot_dim, vpa_dim, nspec_dim, time_dim]


    ! HJL This is a hack to stop the netcdf diagnostic output failing when the variable sizes are too big
    ! The netcdf limit is 2^31-4 bytes, however this limit still causes a problem so 178180992 is used
    ! This limit was found by running the code with different sizes, it is not exact but this is a hack
    ! and is not worth pursuing at the moment
    !<DD>Added protective if statement as nttotext_dim is only defined if io_write_corr_extend
    if(io_write_corr_extend)then
       dimsize(1) = get_netcdf_dim_length(ncid, ri_dim)
       dimsize(2) = get_netcdf_dim_length(ncid, nttotext_dim)
       dimsize(3) = get_netcdf_dim_length(ncid, nakx_dim)
       dimsize(4) = get_netcdf_dim_length(ncid, naky_dim)
       dimsize(5) = get_netcdf_dim_length(ncid, time_big_dim)
       ierr = error_unit()
       total = product(real(dimsize(1:5)))
       if(total .gt. 178180992) then ! number of reals allowed in netcdf
          write(ierr, *) 'WARNING: phi_corr is larger than maximum netcdf size, removing.'
          write(ierr, *) 'Try increasing nwrite_mult to avoid this.'
          io_write_corr_extend = .false.
          write_correlation_extend = .false.
       endif
    endif
    ! HJL End of hack
    if (io_write_corr_extend) then
       phi2_extend_dim = [nttotext_dim, nakx_dim, naky_dim, time_big_dim]
       phi_corr_dim = [ri_dim, nttotext_dim, nakx_dim, naky_dim, time_big_dim]
    end if

    phi_corr_2pi_dim = [ri_dim, nttot_dim, naky_dim, time_dim]

    datestamp(:) = ' '
    timestamp(:) = ' '
    timezone(:) = ' '
    call date_and_time (datestamp, timestamp, timezone)

    !################################################
    !# START WITH BASIC ATTRIBUTES

    call ensure_netcdf_var_exists(ncid, "code_info", NF90_CHAR, char10_dim, code_id)
    call ensure_netcdf_att_exists(ncid, code_id, "long_name", "GS2")
    call ensure_netcdf_att_exists(ncid, code_id, "c1", 'Date: '//trim(datestamp))
    call ensure_netcdf_att_exists(ncid, code_id, "c2", &
         'Time: '//trim(timestamp)//' '//trim(timezone))
    call ensure_netcdf_att_exists(ncid, code_id, "c3", &
         'netCDF version '//trim(nf90_inq_libvers()))
    call ensure_netcdf_att_exists(ncid, code_id, "c4", &
         'Units are determined with respect to reference temperature (T_ref),')
    call ensure_netcdf_att_exists(ncid, code_id, "c5", &
         'reference charge (q_ref), reference mass (mass_ref),')
    call ensure_netcdf_att_exists(ncid, code_id, "c6", &
         'reference field (B_ref), and reference length (a_ref)')
    call ensure_netcdf_att_exists(ncid, code_id, "c7", &
         'from which one may construct rho_ref and vt_ref/a,')
    call ensure_netcdf_att_exists(ncid, code_id, "c8", &
         'which are the basic units of perpendicular length and time.')
    call ensure_netcdf_att_exists(ncid, code_id, "c9", &
         'Macroscopic lengths are normalized to the minor radius.')
    call ensure_netcdf_att_exists(ncid, code_id, "c10", &
         'The difference between rho (normalized minor radius) and rho (gyroradius)')
    call ensure_netcdf_att_exists(ncid, code_id, "c11", &
         'should be clear from the context in which they appear below.')

    !################################################
    !# NOW WE HANDLE THE DIFFERENT OPTIONAL OUTPUT

    !Phi dependent data
    if (fphi > zero) then
       !Volume averaged amplitude data
       call ensure_netcdf_var_exists(ncid, 'phi2', netcdf_real, time_dim, phi2_id)
       call ensure_netcdf_att_exists(ncid, phi2_id, 'long_name', '|Potential**2|')
       call ensure_netcdf_att_exists(ncid, phi2_id, 'units', '(T/q rho/L)**2')

       !Amplitude data per Fourier mode
       call ensure_netcdf_var_exists(ncid, 'phi2_by_mode', netcdf_real, mode_dim, phi2_by_mode_id)

       !Amplitude data summed over ky
       if (ntheta0 > 1) then
          call ensure_netcdf_var_exists(ncid, 'phi2_by_kx', netcdf_real, kx_dim, phi2_by_kx_id)
       end if

       !Amplitude data summed over theta0
       if (naky > 1) then
          call ensure_netcdf_var_exists(ncid, 'phi2_by_ky', netcdf_real, ky_dim, phi2_by_ky_id)
       end if

       !Field at theta=igomega
       call ensure_netcdf_var_exists(ncid, 'phi0', netcdf_real, omega_dim, phi0_id)

       if (write_nl_flux_dist) then
          call ensure_netcdf_var_exists(ncid, 'es_part_flux_dist',  netcdf_real, flux_dist_dim, es_part_flux_dist_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom_flux_par_dist',  netcdf_real, flux_dist_dim, es_mom_flux_par_dist_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom_flux_perp_dist',  netcdf_real, flux_dist_dim, es_mom_flux_perp_dist_id)
          call ensure_netcdf_var_exists(ncid, 'es_heat_flux_dist',  netcdf_real, flux_dist_dim, es_heat_flux_dist_id)
          call ensure_netcdf_var_exists(ncid, 'es_phi_dist',  netcdf_real, phi_dist_dim, es_phi_dist_id)
       end if

       !????
       if (write_sym) then
          call ensure_netcdf_var_exists(ncid, 'es_mom_sym',  netcdf_real, sym_dim, es_mom_sym_id)
       end if
       !????
       if (write_pflux_sym) then !JPL
          call ensure_netcdf_var_exists(ncid, 'es_part_tormom_sym',  netcdf_real, sym_dim, es_part_sym_id)
       end if
       !????
       if (write_pflux_tormom) then !JPL
          call ensure_netcdf_var_exists(ncid, 'es_part_tormom_flux', netcdf_real, flux_dim, es_part_tormom_flux_id)
          call ensure_netcdf_var_exists(ncid, 'es_part_tormom_by_k', netcdf_real, fluxk_dim, es_part_tormom_by_k_id)
       end if

       !Correlation over 2pi range
       if (write_correlation) then
          call ensure_netcdf_var_exists(ncid, 'phi_corr_2pi',  netcdf_real, phi_corr_2pi_dim, phi_corr_2pi_id)
       end if

       !Field as a function of time
       if (write_phi_t) then
          call ensure_netcdf_var_exists(ncid, 'phi_t', netcdf_real, field_dim, phi_t_id)
          call ensure_netcdf_att_exists(ncid, phi_t_id, 'long_name', 'Electrostatic Potential over time')
       end if
    end if

    !Apar dependent output
    if (fapar > zero) then
       !Volume averaged amplitude data
       call ensure_netcdf_var_exists(ncid, 'apar2', netcdf_real, time_dim, apar2_id)
       call ensure_netcdf_att_exists(ncid, apar2_id, 'long_name', 'Total A_par squared')

       !K-resolved amplitude data
       call ensure_netcdf_var_exists(ncid, 'apar2_by_mode', netcdf_real, mode_dim, apar2_by_mode_id)
       call ensure_netcdf_att_exists(ncid, apar2_by_mode_id, 'long_name', 'Apar squared')

       !Amplitude data summed over ky
       if (ntheta0 > 1) then
          call ensure_netcdf_var_exists(ncid, 'apar2_by_kx', netcdf_real, kx_dim, apar2_by_kx_id)
       end if

       !Amplitude data summed over theta0
       if (naky > 1) then
          call ensure_netcdf_var_exists(ncid, 'apar2_by_ky', netcdf_real, ky_dim, apar2_by_ky_id)
       end if

       !Field at theta=igomega
       call ensure_netcdf_var_exists(ncid, 'apar0', netcdf_real, omega_dim, apar0_id)

       !Field as a function of time
       if (write_apar_t) then
          call ensure_netcdf_var_exists(ncid, 'apar_t', netcdf_real, field_dim, apar_t_id)
          call ensure_netcdf_att_exists(ncid, apar_t_id, 'long_name', 'Parallel Magnetic Potential over Time')
       end if
    end if

    !Bpar related data
    if (fbpar > zero) then
       !Volume averaged amplitude data
       call ensure_netcdf_var_exists(ncid, 'bpar2', netcdf_real, time_dim, bpar2_id)
       call ensure_netcdf_att_exists(ncid, bpar2_id, 'long_name', 'Total A_perp squared')

       !K-resolved amplitude data
       call ensure_netcdf_var_exists(ncid, 'bpar2_by_mode', netcdf_real, mode_dim, bpar2_by_mode_id)
       call ensure_netcdf_att_exists(ncid, bpar2_by_mode_id, 'long_name', 'A_perp squared')

       !Amplitude data summed over ky
       if (ntheta0 > 1) then
          call ensure_netcdf_var_exists(ncid, 'bpar2_by_kx', netcdf_real, kx_dim, bpar2_by_kx_id)
       end if

       !Amplitude data summed over theta0
       if (naky > 1) then
          call ensure_netcdf_var_exists(ncid, 'bpar2_by_ky', netcdf_real, ky_dim, bpar2_by_ky_id)
       end if

       !Field at theta=igomega
       call ensure_netcdf_var_exists(ncid, 'bpar0', netcdf_real, omega_dim, bpar0_id)

       !Field as a function of time
       if (write_bpar_t) then
          call ensure_netcdf_var_exists(ncid, 'bpar_t', netcdf_real, field_dim, bpar_t_id)
          call ensure_netcdf_att_exists(ncid, bpar_t_id, 'long_name', 'delta B Parallel over time')
       end if
    end if

    if (write_verr) then
!       call ensure_netcdf_var_exists(ncid, 'errest_by_mode', netcdf_real, err_dim, errest_by_mode_id)
!       call ensure_netcdf_var_exists(ncid, 'lpcoef_by_mode', netcdf_real, lp_dim, lpcoef_by_mode_id)
    end if

    !The current value of the scanned parameter
    if (write_scan_parameter) then
       call ensure_netcdf_var_exists(ncid,'scan_parameter_value',netcdf_real,time_dim,current_scan_parameter_value_id)
       call ensure_netcdf_att_exists(ncid,current_scan_parameter_value_id,'long_name','The current value of the scan parameter')
    end if

    !Roughly the field line averaged summed field amplitudes
    call ensure_netcdf_var_exists(ncid, 'phtot', netcdf_real, mode_dim, phtot_id)
    call ensure_netcdf_att_exists(ncid, phtot_id, 'long_name', 'Field amplitude')

    !Frequency and growth rate estimates
    if (write_omega) then
       call ensure_netcdf_var_exists(ncid, 'omega', netcdf_real, omega_dim, omega_id)
       call ensure_netcdf_var_exists(ncid, 'omegaavg', netcdf_real, omega_dim, omegaavg_id)
    end if

    !Heating rate diagnostics
    if (write_heating) then
       call ensure_netcdf_var_exists(ncid, 'h_energy', netcdf_real, time_dim, h_energy_id)
       call ensure_netcdf_var_exists(ncid, 'h_energy_dot', netcdf_real, time_dim, h_energy_dot_id)
       call ensure_netcdf_var_exists(ncid, 'h_antenna', netcdf_real, time_dim, h_antenna_id)
       call ensure_netcdf_var_exists(ncid, 'h_eapar', netcdf_real, time_dim, h_eapar_id)
       call ensure_netcdf_var_exists(ncid, 'h_ebpar', netcdf_real, time_dim, h_ebpar_id)

       call ensure_netcdf_var_exists(ncid, 'h_delfs2', netcdf_real, flux_dim, h_delfs2_id)
       call ensure_netcdf_var_exists(ncid, 'h_hs2', netcdf_real, flux_dim, h_hs2_id)
       call ensure_netcdf_var_exists(ncid, 'h_phis2', netcdf_real, flux_dim, h_phis2_id)
       call ensure_netcdf_var_exists(ncid, 'h_hypervisc', netcdf_real, flux_dim, h_hypervisc_id)
       call ensure_netcdf_var_exists(ncid, 'h_hyperres', netcdf_real, flux_dim, h_hyperres_id)
       call ensure_netcdf_var_exists(ncid, 'h_collisions', netcdf_real, flux_dim, h_collisions_id)
       call ensure_netcdf_var_exists(ncid, 'h_gradients', netcdf_real, flux_dim, h_gradients_id)
       call ensure_netcdf_var_exists(ncid, 'h_hypercoll', netcdf_real, flux_dim, h_hypercoll_id)
       call ensure_netcdf_var_exists(ncid, 'h_imp_colls', netcdf_real, flux_dim, h_imp_colls_id)
       call ensure_netcdf_var_exists(ncid, 'h_heating', netcdf_real, flux_dim, h_heating_id)
       call ensure_netcdf_att_exists(ncid, h_heating_id, 'long_name', 'Total heating by species')

       call ensure_netcdf_var_exists(ncid, 'hk_energy', netcdf_real, mode_dim, hk_energy_id)
       call ensure_netcdf_var_exists(ncid, 'hk_energy_dot', netcdf_real, mode_dim, hk_energy_dot_id)
       call ensure_netcdf_var_exists(ncid, 'hk_antenna', netcdf_real, mode_dim, hk_antenna_id)
       call ensure_netcdf_var_exists(ncid, 'hk_eapar', netcdf_real, mode_dim, hk_eapar_id)
       call ensure_netcdf_var_exists(ncid, 'hk_ebpar', netcdf_real, mode_dim, hk_ebpar_id)

       call ensure_netcdf_var_exists(ncid, 'hk_delfs2', netcdf_real, fluxk_dim, hk_delfs2_id)
       call ensure_netcdf_var_exists(ncid, 'hk_hs2', netcdf_real, fluxk_dim, hk_hs2_id)
       call ensure_netcdf_var_exists(ncid, 'hk_phis2', netcdf_real, fluxk_dim, hk_phis2_id)
       call ensure_netcdf_var_exists(ncid, 'hk_hypervisc', netcdf_real, fluxk_dim, hk_hypervisc_id)
       call ensure_netcdf_var_exists(ncid, 'hk_hyperres', netcdf_real, fluxk_dim, hk_hyperres_id)
       call ensure_netcdf_var_exists(ncid, 'hk_collisions', netcdf_real, fluxk_dim, hk_collisions_id)
       call ensure_netcdf_var_exists(ncid, 'hk_gradients', netcdf_real, fluxk_dim, hk_gradients_id)
       call ensure_netcdf_var_exists(ncid, 'hk_imp_colls', netcdf_real, fluxk_dim, hk_imp_colls_id)
       call ensure_netcdf_var_exists(ncid, 'hk_hypercoll', netcdf_real, fluxk_dim, hk_hypercoll_id)
       call ensure_netcdf_var_exists(ncid, 'hk_heating', netcdf_real, fluxk_dim, hk_heating_id)
       call ensure_netcdf_att_exists(ncid, hk_heating_id, 'long_name', 'Total heating by species and mode')
    end if

    !Correlation function (phi only so should probably include if fphi>0)
    if (io_write_corr_extend) then
       call ensure_netcdf_var_exists(ncid, 'phi2_extend',  netcdf_real, phi2_extend_dim, phi2_extend_id)
       call ensure_netcdf_var_exists(ncid, 'phi_corr',  netcdf_real, phi_corr_dim, phi_corr_id)
    end if
    !------------------------------------------------

    !################################################
    !# NOW EXIT DEFINITION MODE
    status = nf90_enddef (ncid)  ! out of definition mode
    if (status /= NF90_NOERR) call netcdf_error (status)

# endif
  end subroutine define_vars

  !> Define the vars in the movie file. A fair bit of duplication wrt main
  !> file -- may want to think about merging movie data into main file.
  !> Main reason for using a new file is we can only have one unlimited dimension
  !> per file (with new versions of NETCDF this is no longer true) and we
  !> probably don't want to write movie data as often as other data.
  !> However, at the moment we precalculate the number of outputs, so it
  !> isn't actually unlimited
  subroutine define_vars_mymovie
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: NF90_INT, NF90_GLOBAL
    use netcdf, only: nf90_def_var, nf90_put_att, nf90_enddef
    use netcdf_utils, only: netcdf_real, ensure_netcdf_att_exists, ensure_netcdf_var_exists
    implicit none
    integer :: status

    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then
       !Setup field shape
       xmode_dim = [nx_dim, ny_dim, nth_dim, time_movie_dim]

       !Grid sizes
       call ensure_netcdf_var_exists(ncid_movie, 'nx', NF90_INT, nx_id)
       call ensure_netcdf_var_exists(ncid_movie, 'ny', NF90_INT, ny_id)
       call ensure_netcdf_var_exists(ncid_movie, 'ntheta', NF90_INT, nth_id)

       !Grids
       call ensure_netcdf_var_exists(ncid_movie, 'tm', netcdf_real, time_movie_dim, time_movie_id)
       call ensure_netcdf_att_exists(ncid_movie, time_movie_id, 'long_name', 'Time')
       call ensure_netcdf_att_exists(ncid_movie, time_movie_id, 'units', 'L/vt')
       call ensure_netcdf_var_exists(ncid_movie, 'x', netcdf_real, nx_dim, x_id)
       call ensure_netcdf_att_exists(ncid_movie, x_id, 'long_name', 'x / rho')
       call ensure_netcdf_var_exists(ncid_movie, 'y', netcdf_real, ny_dim, y_id)
       call ensure_netcdf_att_exists(ncid_movie, y_id, 'long_name', 'y / rho')
       call ensure_netcdf_var_exists(ncid_movie, 'theta', netcdf_real, nth_dim, th_id)

       !Field data for movie
       if(fphi > zero) then
          call ensure_netcdf_var_exists(ncid_movie, 'phi_by_xmode', netcdf_real, xmode_dim, phi_by_xmode_id)
       end if
       ! TEMPORARY BD | Disable DD to match gs2_diagnostics
       !call ensure_netcdf_var_exists(ncid_movie, 'density_by_xmode', netcdf_real, xmode_dim, apar_by_xmode_id)

       if(fapar > zero) then
          call ensure_netcdf_var_exists(ncid_movie, 'apar_by_xmode', netcdf_real, xmode_dim, apar_by_xmode_id)
       end if

       if(fbpar > zero) then
          call ensure_netcdf_var_exists(ncid_movie, 'bpar_by_xmode', netcdf_real, xmode_dim, bpar_by_xmode_id)
       end if

       !Finish defining
       status = nf90_enddef (ncid_movie)  ! out of definition mode
       if (status /= NF90_NOERR) call netcdf_error (status)
    endif
#endif
  end subroutine define_vars_mymovie

  !> Write the normalised fields to netCDF
  subroutine nc_eigenfunc(file_id, phase)
#ifdef NETCDF
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Normalising phase factor (nominally \(\phi\) at the outboard midplane)
    complex, dimension(:,:), intent (in) :: phase
# ifdef NETCDF

    complex, dimension(-ntgrid:ntgrid, ntheta0, naky) :: tmp
    integer :: ig

    call netcdf_write_complex(file_id, "phase", phase, ["ri", "kx", "ky"], long_name="Normalising phase")

    if (fphi > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = phi(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "phi_norm", tmp, final_field_dims)
    end if

    if (fapar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = apar(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "apar_norm", tmp, final_field_dims)
    end if

    if (fbpar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = bpar(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "bpar_norm", tmp, final_field_dims)
    end if
# endif
  end subroutine nc_eigenfunc

  !> Write various moments to netcdf
  subroutine nc_write_moments (nout, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1, ob_midplane) 
    use theta_grid, only: ntgrid
    implicit none
    !> Total perturbed density
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: ntot
    !> Non-adiabatic part of the perturbed density
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: density
    !> Parallel velocity
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: upar
    !> Parallel temperature
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: tpar
    !> Perpendicular temperature
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: tperp
    !> Parallel heat flux
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: qparflux
    !> Pressure part of the particle flux
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: pperpj1
    !> Pressure part of the heat flux
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: qpperpj1
    !> Current output timestep
    integer, intent (in) :: nout
    !> If true, write out moments at the outboard midplane only,
    !> otherwise write them out along the entire flux surface
    logical, optional, intent(in) :: ob_midplane
#ifdef NETCDF
    integer, dimension (6) :: start6
    integer, dimension (5) :: start5

    character(*), dimension(*), parameter :: mom_dim = [character(len=7)::"ri", "kx", "ky", "species", "t"]
    character(*), dimension(*), parameter :: mom_t_dim = [character(len=7)::"ri", "theta", "kx", "ky", "species", "t"]

    if (present(ob_midplane)) then
       if (ob_midplane ) then
          start5 = [1, 1, 1, 1, nout]

          call netcdf_write_complex(ncid, "ntot0", ntot(0, :, :, :), mom_dim, start=start5, &
               long_name="Total perturbed density (theta=0) over time")
          call netcdf_write_complex(ncid, "density0", density(0, :, :, :), mom_dim, start=start5, &
               long_name="Nonadiabatic piece of perturbed density (theta=0) over time")
          call netcdf_write_complex(ncid, "upar0", upar(0, :, :, :), mom_dim, start=start5, &
               long_name="Upar (theta=0) over time")
          call netcdf_write_complex(ncid, "tpar0", tpar(0, :, :, :), mom_dim, start=start5, &
               long_name="Tpar (theta=0) over time")
          call netcdf_write_complex(ncid, "tperp0", tperp(0, :, :, :), mom_dim, start=start5, &
               long_name="Tperp (theta=0) over time")
          call netcdf_write_complex(ncid, "qparflux0", qparflux(0, :, :, :), mom_dim, start=start5, &
               long_name="Qparflux (theta=0) over time")
          call netcdf_write_complex(ncid, "pperpj10", pperpj1(0, :, :, :), mom_dim, start=start5, &
               long_name="Pperpj1 (theta=0) over time")
          call netcdf_write_complex(ncid, "qpperpj10", qpperpj1(0, :, :, :), mom_dim, start=start5, &
               long_name="Qpperpj1 (theta=0) over time")
          return
       endif
    endif

    start6 = [1, 1, 1, 1, 1, nout]

    call netcdf_write_complex(ncid, 'ntot_t', ntot, mom_t_dim, start=start6, &
          long_name="Total perturbed density over time")
    call netcdf_write_complex(ncid, 'density_t', density, mom_t_dim, start=start6, &
          long_name="Nonadiabatic piece of perturbed density over time")
    call netcdf_write_complex(ncid, 'upar_t', upar, mom_t_dim, start=start6, &
          long_name="Upar over time")
    call netcdf_write_complex(ncid, 'tpar_t', tpar, mom_t_dim, start=start6, &
          long_name="Tpar over time")
    call netcdf_write_complex(ncid, 'tperp_t', tperp, mom_t_dim, start=start6, &
          long_name="Tperp over time")
    call netcdf_write_complex(ncid, 'qparflux_t', qparflux, mom_t_dim, start=start6, &
          long_name="Qparflux over time")
    call netcdf_write_complex(ncid, 'pperpj1_t', pperpj1, mom_t_dim, start=start6, &
          long_name="Pperpj1 over time")
    call netcdf_write_complex(ncid, 'qpperpj1_t', qpperpj1, mom_t_dim, start=start6, &
          long_name="Qpperpj1 over time")

# endif
  end subroutine nc_write_moments

  !> Write the fields to netCDF, overwriting existing values
  subroutine nc_final_fields
#ifdef NETCDF
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    implicit none

    if (fphi > zero) then
       call netcdf_write_complex(ncid, "phi", phinew, final_field_dims, &
            long_name="Electrostatic Potential", units="T/q rho/L")
    end if

    if (fapar > zero) then
       call netcdf_write_complex(ncid, "apar", aparnew, final_field_dims, &
            long_name="Parallel Magnetic Potential")
    end if

    if (fbpar > zero) then
       call netcdf_write_complex(ncid, "bpar", bparnew, final_field_dims, &
            long_name="delta B Parallel")
    end if
# endif
  end subroutine nc_final_fields

  !> Write \(E_\parallel\) to netCDF, overwriting existing values
  subroutine nc_final_epar (epar)
    implicit none
    complex, dimension (:,:,:), intent (in) :: epar
# ifdef NETCDF
    character(*), dimension(*), parameter :: dim_names = [character(5)::"ri", "theta", "kx", "ky"]

    call netcdf_write_complex(ncid, "epar", epar, dim_names)
# endif
  end subroutine nc_final_epar

  !> Write the various low-order moments to the output netCDF file,
  !> overwriting existing values
  subroutine nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    implicit none
    complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
    complex, dimension (:,:,:,:), intent (in) :: qparflux, pperpj1, qpperpj1
    character(*), dimension(*), parameter :: dim_names = [character(7)::"ri", "theta", "kx", "ky", "species"]

#ifdef NETCDF
    call netcdf_write_complex(ncid, 'ntot', ntot, dim_names)
    call netcdf_write_complex(ncid, 'density', density, dim_names)
    call netcdf_write_complex(ncid, 'upar', upar, dim_names)
    call netcdf_write_complex(ncid, 'tpar', tpar, dim_names)
    call netcdf_write_complex(ncid, 'tperp', tperp, dim_names)
    call netcdf_write_complex(ncid, 'qparflux', qparflux, dim_names)
    call netcdf_write_complex(ncid, 'pperpj1', pperpj1, dim_names)
    call netcdf_write_complex(ncid, 'qpperpj1', qpperpj1, dim_names)
#endif
  end subroutine nc_final_moments

  !> Write moments as: volume averaged, field line averaged, and evaulated at theta=0
  subroutine nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, ntot20_by_mode, &
       phi00, ntot00, density00, upar00, tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: ntot2, ntot20
    real, dimension (:,:,:), intent (in) :: ntot2_by_mode, ntot20_by_mode
    real, dimension (:,:,:), intent (in) :: tpar2_by_mode, tperp2_by_mode
    complex, dimension (:), intent (in) :: phi00
    complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
#ifdef NETCDF
    integer, dimension (2) :: start2
    integer, dimension (3) :: start3
    integer, dimension (4) :: start4

    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]
    character(*), dimension(*), parameter :: loop_phi_dim = [character(7)::"ri", "kx", "t"]
    character(*), dimension(*), parameter :: loop_mom_dim = [character(7)::"ri", "kx", "species", "t"]

    start2 = [1, nout]
    start3 = [1, 1, nout]
    start4 = [1, 1, 1, nout]
    
    ! Volume averaged function of theta
    ! <DD>Should we only do this if .not.obmid?
    call neasyf_write (ncid, "ntot2", ntot2, dim_names=flux_dim, start=start2)
    call neasyf_write (ncid, "ntot2_by_mode", ntot2_by_mode, dim_names=fluxk_dim, start=start4)
    call neasyf_write (ncid, "tpar2_by_mode", tpar2_by_mode, dim_names=fluxk_dim, start=start4)
    call neasyf_write (ncid, "tperp2_by_mode", tperp2_by_mode, dim_names=fluxk_dim, start=start4)

    ! "Volume averaged" data from outboard midplane (really just evaluated at theta=0
    ! and then weighted with Fourier coeff factor (typically 0.5).
    ! <DD>Should we only do this if obmid?
    call neasyf_write (ncid, "ntot20", ntot20, dim_names=flux_dim, start=start2, &
         long_name='Density**2 at theta=0')
    call neasyf_write (ncid, "ntot20_by_mode", ntot20_by_mode, dim_names=fluxk_dim, start=start4)

    ! Fieldline averaged data, note this includes phi which is not a moment
    call netcdf_write_complex(ncid, "phi00", phi00, loop_phi_dim, start=start3)
    call netcdf_write_complex(ncid, "ntot00", ntot00, loop_mom_dim, start=start4)
    call netcdf_write_complex(ncid, "density00", density00, loop_mom_dim, start=start4)
    call netcdf_write_complex(ncid, "upar00", upar00, loop_mom_dim, start=start4)
    call netcdf_write_complex(ncid, "tpar00", tpar00, loop_mom_dim, start=start4)
    call netcdf_write_complex(ncid, "tperp00", tperp00, loop_mom_dim, start=start4)
#endif
  end subroutine nc_loop_moments

  !> Write the right-hand sides of the field equations to netCDF,
  !> overwriting existing values
  subroutine nc_final_an (antot, antota, antotp)
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
# endif
    implicit none
    complex, dimension (:,:,:), intent(in) :: antot, antota, antotp
# ifdef NETCDF
    if (fphi > zero) then
       call netcdf_write_complex(ncid, "antot", antot, final_field_dims, &
            long_name="Velocity moment of g used in calculating phi")
    end if

    if (fapar > zero) then
       call netcdf_write_complex(ncid, "antota", antota, final_field_dims, &
            long_name="Velocity space moment of g used in calculating apar")
    end if

    if (fbpar > zero) then
       call netcdf_write_complex(ncid, "antotp", antotp, final_field_dims, &
            long_name="Velocity space moment of g used in calculating bpar")
    end if
# endif
  end subroutine nc_final_an

  !> Write the energy exchange to netCDF
  subroutine nc_exchange (nout, exchange, exchange_avg)
# ifdef NETCDF
    use run_parameters, only: fphi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    !> Energy exchange
    real, dimension (:,:,:), intent (in) :: exchange
    !> Volume-averaged energy exchange
    real, dimension (:), intent (in) :: exchange_avg
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]

    if (fphi > zero) then
       call neasyf_write (ncid, "es_exchange", exchange_avg, &
            dim_names=flux_dim, start=[1, nout])
       call neasyf_write (ncid, "es_exchange_by_k", exchange, &
            dim_names=fluxk_dim, start=[1, 1, 1, nout])
    end if
# endif
  end subroutine nc_exchange

  !> Write various heat fluxes to netCDF
  subroutine nc_qflux (nout, qheat, qmheat, qbheat, &
       heat_par,  mheat_par,  bheat_par, &
       heat_perp, mheat_perp, bheat_perp, &
       heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot)
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: qheat, qmheat, qbheat
    real, dimension (:), intent (in) :: heat_par, mheat_par, bheat_par
    real, dimension (:), intent (in) :: heat_perp, mheat_perp, bheat_perp
    real, dimension (:), intent (in) :: heat_fluxes, mheat_fluxes, bheat_fluxes
    real, dimension (:,:), intent (in) :: x_qmflux
    real, intent (in) :: hflux_tot
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]
    character(*), dimension(*), parameter :: fluxx_dim = [character(7)::"kx", "species", "t"]

    integer, dimension (2) :: start2
    integer, dimension (3) :: start3
    integer, dimension (4) :: start4

    start2 = [1, nout]
    start3 = [1, 1, nout]
    start4 = [1, 1, 1, nout]

    if (fphi > zero) then
       call neasyf_write(ncid, "es_heat_flux", heat_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_heat_par", heat_par, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_heat_perp", heat_perp, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_heat_by_k", qheat, dim_names=fluxk_dim, start=start4)
    end if

    if (fapar > zero) then
       call neasyf_write(ncid, "apar_heat_flux", mheat_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "apar_heat_par", mheat_par, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "apar_heat_perp", mheat_perp, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "apar_heat_by_k", qmheat, dim_names=fluxk_dim, start=start4)
       call neasyf_write(ncid, "apar_heat_by_x", x_qmflux, dim_names=fluxx_dim, start=start3)
    end if

    if (fbpar > zero) then
       call neasyf_write(ncid, "bpar_heat_flux", bheat_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "bpar_heat_par", bheat_par, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "bpar_heat_perp", bheat_perp, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "bpar_heat_by_k", qbheat, dim_names=fluxk_dim, start=start4)
    end if
    
    call neasyf_write(ncid, "hflux_tot", hflux_tot, dim_names=["t"], start=[nout])
# endif
  end subroutine nc_qflux

  !> FIXME : Add documentation  
  subroutine nc_pflux (nout, pflux, pmflux, pbflux, &
       part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
    real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
    real, intent (in) :: zflux_tot
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]

    integer, dimension (2) :: start2
    integer, dimension (4) :: start4

    start2 = [1, nout]
    start4 = [1, 1, 1, nout]
    
    if (fphi > zero) then
       call neasyf_write(ncid, "es_part_flux", part_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_part_by_k", pflux, dim_names=fluxk_dim, start=start4)
    end if

    if (fapar > zero) then
       call neasyf_write(ncid, "apar_part_flux", mpart_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "apar_part_by_k", pmflux, dim_names=fluxk_dim, start=start4)
    end if

    if (fbpar > zero) then
       call neasyf_write(ncid, "bpar_part_flux", bpart_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "bpar_part_by_k", pbflux, dim_names=fluxk_dim, start=start4)
    end if
    
    call neasyf_write(ncid, "zflux_tot", zflux_tot, dim_names=["t"], start=[nout])
# endif
  end subroutine nc_pflux

  !> Write toroidal angular momentum flux to netCDF
  subroutine nc_pflux_tormom (nout, pflux_tormom, part_tormom_fluxes)
# ifdef NETCDF
    use run_parameters, only: fphi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux_tormom
    real, dimension(:), intent (in) :: part_tormom_fluxes
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]

    integer, dimension (2) :: start2
    integer, dimension (4) :: start4

    start2 = [1, nout]
    start4 = [1, 1, 1, nout]

    if (fphi > zero) then
       call neasyf_write(ncid, "es_part_tormom_flux", part_tormom_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_part_tormom_by_k", pflux_tormom, dim_names=fluxk_dim, start=start4)
    end if

# endif
  end subroutine nc_pflux_tormom

  !> FIXME : Add documentation  
  subroutine nc_vflux (nout, vflux, vmflux, vbflux, &
       mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
       vflux_par, vflux_perp, vflux0, vflux1)
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
    real, dimension (:,:,:), intent (in) :: vflux_par, vflux_perp
    real, dimension (:,:,:), intent (in) :: vflux0, vflux1
    real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
    real, intent (in) :: vflux_tot
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dim = [character(7)::"species", "t"]
    character(*), dimension(*), parameter :: fluxk_dim = [character(7)::"kx", "ky", "species", "t"]

    integer, dimension (2) :: start2
    integer, dimension (4) :: start4

    start2 = [1, nout]
    start4 = [1, 1, 1, nout]

    if (fphi > zero) then
       call neasyf_write(ncid, "es_mom_flux", mom_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "es_mom_by_k", vflux, dim_names=fluxk_dim, start=start4)
       call neasyf_write(ncid, "es_parmom_by_k", vflux_par, dim_names=fluxk_dim, start=start4)
       call neasyf_write(ncid, "es_perpmom_by_k", vflux_perp, dim_names=fluxk_dim, start=start4)
       call neasyf_write(ncid, "es_mom0_by_k", vflux0, dim_names=fluxk_dim, start=start4)
       call neasyf_write(ncid, "es_mom1_by_k", vflux1, dim_names=fluxk_dim, start=start4)
    end if

    if (fapar > zero) then
       call neasyf_write(ncid, "apar_mom_flux", mmom_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "apar_mom_by_k", vmflux, dim_names=fluxk_dim, start=start4)
    end if

    if (fbpar > zero) then
       call neasyf_write(ncid, "bpar_mom_flux", bmom_fluxes, dim_names=flux_dim, start=start2)
       call neasyf_write(ncid, "bpar_mom_by_k", vbflux, dim_names=fluxk_dim, start=start4)
    end if
    
    call neasyf_write(ncid, "vflux_tot", vflux_tot, dim_names=["t"], start=[nout])
# endif
  end subroutine nc_vflux

  !> Output the poloidal distributions of the fluxes of particles, perpendicular
  !> momentum, parallel momentum, and energy ! JRB
  subroutine nc_loop_dist (nout, part_fluxes_dist, mom_fluxes_par_dist, mom_fluxes_perp_dist, heat_fluxes_dist, phi_dist)

    use species, only: nspec
    use run_parameters, only: fphi
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:), intent (in) :: part_fluxes_dist, mom_fluxes_par_dist, mom_fluxes_perp_dist, heat_fluxes_dist
    real, dimension (2,-ntgrid:ntgrid), intent (in) :: phi_dist
# ifdef NETCDF
    integer, dimension (3) :: start1, count1
    integer, dimension (3) :: start2, count2
    integer :: status

    start1(1) = 1
    start1(2) = 1
    start1(3) = nout

    count1(1) = 2*ntgrid+1
    count1(2) = nspec
    count1(3) = 1

    start2(1) = 1
    start2(2) = 1
    start2(3) = nout

    count2(1) = 2
    count2(2) = 2*ntgrid+1
    count2(3) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_part_flux_dist_id, part_fluxes_dist, start=start1, count=count1)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_flux_dist_id)
       status = nf90_put_var (ncid, es_mom_flux_par_dist_id, mom_fluxes_par_dist, start=start1, count=count1)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_flux_par_dist_id)
       status = nf90_put_var (ncid, es_mom_flux_perp_dist_id, mom_fluxes_perp_dist, start=start1, count=count1)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_flux_perp_dist_id)
       status = nf90_put_var (ncid, es_heat_flux_dist_id, heat_fluxes_dist, start=start1, count=count1)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_flux_dist_id)
       status = nf90_put_var (ncid, es_phi_dist_id, phi_dist, start=start2, count=count2)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_phi_dist_id)
    end if

# endif
  end subroutine nc_loop_dist
  
  !> Write the momentum flux as a function of \((v_\parallel, \theta, t)\) to netCDF
  subroutine nc_loop_sym (nout, vflx_sym)
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use species, only: nspec
    use run_parameters, only: fphi
    use le_grids, only: negrid, nlambda
    use netcdf, only: nf90_put_var
    implicit none
# endif
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:,:), intent (in) :: vflx_sym
# ifdef NETCDF
    integer, dimension (4) :: start, count
    integer :: status

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout
    
    count(1) = 2*ntgrid+1
    count(2) = nlambda*negrid
    count(3) = nspec
    count(4) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_mom_sym_id, vflx_sym, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_sym_id)
    end if

# endif
  end subroutine nc_loop_sym
  
  !> Write the particle flux contribution to toroidal momentum flux to netCDF
  subroutine nc_loop_partsym_tormom (nout, partflx_sym) 
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use species, only: nspec
    use run_parameters, only: fphi
    use le_grids, only: negrid, nlambda
    use netcdf, only: nf90_put_var
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:,:), intent (in) :: partflx_sym
# ifdef NETCDF
    integer, dimension (4) :: start, count
    integer :: status

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout
    
    count(1) = 2*ntgrid+1
    count(2) = nlambda*negrid
    count(3) = nspec
    count(4) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_part_sym_id, partflx_sym, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_sym_id)
    end if

# endif
  end subroutine nc_loop_partsym_tormom

  !> Write the correlation function over the physical domain to netCDF
  subroutine nc_loop_corr (nout,phi_2pi_corr)
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use convert, only: c2r
    use run_parameters, only: fphi
    use kt_grids, only: naky
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    !> The correlation function over the physical domain
    complex, dimension (-ntgrid:,:), intent (in) :: phi_2pi_corr
# ifdef NETCDF
    integer, dimension (4) :: start4, count4
    integer :: status

    real, dimension (2, size(phi_2pi_corr,1), size(phi_2pi_corr,2)) :: ri2

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout

    count4(1) = 2
    count4(2) = 2*ntgrid+1
    count4(3) = naky
    count4(4) = 1

    if (fphi > zero) then
       call c2r (phi_2pi_corr, ri2)
       status = nf90_put_var (ncid, phi_corr_2pi_id, ri2, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_corr_2pi_id)
    end if

# endif
  end subroutine nc_loop_corr

  !> Write the correlation function over the extended domain to netCDF
  subroutine nc_loop_corr_extend (nout_big,phi_corr,phi2_extend)
# ifdef NETCDF
    use file_utils, only: error_unit
    use convert, only: c2r
    use run_parameters, only: fphi
    use kt_grids, only: ntheta0, naky, jtwist_out
    use theta_grid, only: ntgrid
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> The current output timestep for "large" diagnostics
    integer, intent (in) :: nout_big
    real, dimension (:,:,:), intent (in) :: phi2_extend
    complex, dimension (:,:,:), intent (in) :: phi_corr
# ifdef NETCDF
    integer, dimension (4) :: start4, count4
    integer, dimension (5) :: start5, count5
    integer :: status, ierr
    real, dimension (2, size(phi_corr,1), size(phi_corr,2), size(phi_corr,3)) :: ri3

    ierr = error_unit()

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout_big
    
    count4(1) = (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1)
    count4(2) = ntheta0
    count4(3) = naky
    count4(4) = 1

    start5(1) = 1
    start5(2) = 1
    start5(3) = 1
    start5(4) = 1
    start5(5) = nout_big

    count5(1) = 2
    count5(2) = (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1)
    count5(3) = ntheta0
    count5(4) = naky
    count5(5) = 1

    if (fphi > zero) then
       call c2r (phi_corr, ri3)
       status = nf90_put_var (ncid, phi_corr_id, ri3, start=start5, count=count5)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_corr_id)
       status = nf90_put_var (ncid, phi2_extend_id, phi2_extend, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_extend_id)
    end if

# endif
  end subroutine nc_loop_corr_extend

  !> Write various fields, heating, and frequency diagnostics to netCDF
  subroutine nc_loop (nout, time, fluxfac, &
       phi0,   phi2,   phi2_by_mode, &
       apar0,  apar2,  apar2_by_mode, &
       bpar0, bpar2, bpar2_by_mode, &
       h, hk, omega, omegaavg, woutunits, phitot, write_omega, write_heating)
    use gs2_heating, only: heating_diagnostics
# ifdef NETCDF
    use gs2_heating, only: hk_repack
    use run_parameters, only: fphi, fapar, fbpar
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use species, only: nspec
    use convert, only: c2r
    use fields_arrays, only: phinew, aparnew, bparnew
    use parameter_scan_arrays, only: write_scan_parameter,current_scan_parameter_value
    use netcdf, only: nf90_put_var
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    !> Current simulation time, fields squared
    real, intent (in) :: time, phi2, apar2, bpar2
    !> Scaling factor for modes (FIXME?), normalisation factor for frequency
    real, dimension (:), intent (in) :: fluxfac, woutunits
    !> Fields, frequency and time-averaged frequency
    complex, dimension(:,:), intent (in) :: phi0, apar0, bpar0, omega, omegaavg
    !> FIXME: add documentation
    real, dimension(:,:), intent (in) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode, phitot
    !> Heating diagnostics
    type(heating_diagnostics), intent (in) :: h
    !> Heating per-mode
    type(heating_diagnostics), dimension(:,:), intent (in) :: hk
    logical, intent(in) :: write_omega, write_heating
# ifdef NETCDF
    real, dimension (ntheta0) :: field2_by_kx
    real, dimension (naky) :: field2_by_ky
    real, dimension (2, ntheta0, naky) :: ri2
    real, dimension (ntheta0, naky, nspec) :: tmps
    complex, dimension (ntheta0, naky) :: tmp
    integer, dimension (4) :: start0, count0, start4, count4
    integer, dimension (3) :: start, count, starth, counth
    integer, dimension (2) :: startx, countx, starty, county, starts, counts
    integer :: status, it, ik
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer, dimension (5) :: start5, count5

    start5(1) = 1
    start5(2) = 1
    start5(3) = 1
    start5(4) = 1
    start5(5) = nout
  
    count5(1) = 2
    count5(2) = 2*ntgrid+1
    count5(3) = ntheta0
    count5(4) = naky
    count5(5) = 1

    call neasyf_write(ncid, "t", time, dim_names=["t"], start=[nout])

    start(1) = 1
    start(2) = 1
    start(3) = nout

    count(1) = ntheta0
    count(2) = naky
    count(3) = 1

    start0(1) = 1
    start0(2) = 1
    start0(3) = 1
    start0(4) = nout

    count0(1) = 2
    count0(2) = ntheta0
    count0(3) = naky
    count0(4) = 1

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout

    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    starty(1) = 1
    starty(2) = nout

    county(1) = naky
    county(2) = 1

    startx(1) = 1
    startx(2) = nout

    countx(1) = ntheta0
    countx(2) = 1

    starts(1) = 1
    starts(2) = nout

    counts(1) = nspec
    counts(2) = 1

    starth(1) = 1 
    starth(2) = 1
    starth(3) = nout

    counth(1) = nspec
    counth(2) = 7
    counth(3) = 1

    if (fphi > zero) then

       if(write_phi_t) then
          call c2r (phinew, ri3)
          !ri_phi_t(:,:,:,:,1) = ri3(:,:,:,:)
          status = nf90_put_var (ncid, phi_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_t_id)
       end if

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(phi2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, phi2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(phi2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, phi2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_ky_id)
       end if

       call c2r (phi0, ri2) 
       status = nf90_put_var (ncid, phi0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi0_id)
       status = nf90_put_var (ncid, phi2_by_mode_id, phi2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_mode_id)
       status = nf90_put_var (ncid, phi2_id, phi2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_id)
    end if

    if (fapar > zero) then

       if(write_apar_t) then
          call c2r (aparnew, ri3)
          !ri_apar_t(:,:,:,:,1) = ri3(:,:,:,:)
          status = nf90_put_var (ncid, apar_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_t_id)
       end if

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(apar2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, apar2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(apar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, apar2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_ky_id)
       end if

       call c2r (apar0, ri2) 
       status = nf90_put_var (ncid, apar0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar0_id)
       status = nf90_put_var (ncid, apar2_by_mode_id, apar2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_mode_id)
       status = nf90_put_var (ncid, apar2_id, apar2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_id)
    end if

    if (fbpar > zero) then

       if(write_bpar_t) then
          call c2r (bparnew, ri3)
          !ri_bpar_t(:,:,:,:,1) = ri3(:,:,:,:)
          status = nf90_put_var (ncid, bpar_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_t_id)
       end if

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(bpar2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, bpar2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(bpar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, bpar2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_ky_id)
       end if

       call c2r (bpar0, ri2) 
       status = nf90_put_var (ncid, bpar0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar0_id)
       status = nf90_put_var (ncid, bpar2_by_mode_id, bpar2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_mode_id)
       status = nf90_put_var (ncid, bpar2_id, bpar2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_id)
    end if

    if (write_heating) then
       status = nf90_put_var (ncid, h_energy_id, h%energy, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_energy_id)
       status = nf90_put_var (ncid, h_energy_dot_id, h%energy_dot, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_energy_dot_id)
       status = nf90_put_var (ncid, h_antenna_id, h%antenna, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_antenna_id)
       status = nf90_put_var (ncid, h_eapar_id, h%eapar, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_eapar_id)
       status = nf90_put_var (ncid, h_ebpar_id, h%ebpar, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_ebpar_id)

       status = nf90_put_var (ncid, h_delfs2_id, h%delfs2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_delfs2_id)
       status = nf90_put_var (ncid, h_hs2_id, h%hs2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hs2_id)
       status = nf90_put_var (ncid, h_phis2_id, h%phis2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_phis2_id)
       status = nf90_put_var (ncid, h_hypervisc_id, h%hypervisc, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hypervisc_id)
       status = nf90_put_var (ncid, h_hyperres_id, h%hyperres, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hyperres_id)
       status = nf90_put_var (ncid, h_collisions_id, h%collisions, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_collisions_id)
       status = nf90_put_var (ncid, h_gradients_id, h%gradients, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_gradients_id)
       status = nf90_put_var (ncid, h_imp_colls_id, h%imp_colls, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_imp_colls_id)
       status = nf90_put_var (ncid, h_hypercoll_id, h%hypercoll, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hypercoll_id)
       status = nf90_put_var (ncid, h_heating_id, h%heating, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_heating_id)

       status = nf90_put_var (ncid, hk_energy_id, hk%energy, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_energy_id)
       status = nf90_put_var (ncid, hk_energy_dot_id, hk%energy_dot, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_energy_dot_id)
       status = nf90_put_var (ncid, hk_antenna_id, hk%antenna, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_antenna_id)
       status = nf90_put_var (ncid, hk_eapar_id, hk%eapar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_eapar_id)
       status = nf90_put_var (ncid, hk_ebpar_id, hk%ebpar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_ebpar_id)

       call hk_repack (hk, 1, tmps)
       status = nf90_put_var (ncid, hk_hypervisc_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hypervisc_id)

       call hk_repack (hk, 2, tmps)
       status = nf90_put_var (ncid, hk_hyperres_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hyperres_id)

       call hk_repack (hk, 3, tmps)
       status = nf90_put_var (ncid, hk_collisions_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_collisions_id)

       call hk_repack (hk, 4, tmps)
       status = nf90_put_var (ncid, hk_gradients_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_gradients_id)

       call hk_repack (hk, 5, tmps)
       status = nf90_put_var (ncid, hk_heating_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_heating_id)

       call hk_repack (hk, 6, tmps)
       status = nf90_put_var (ncid, hk_hypercoll_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hypercoll_id)

       call hk_repack (hk, 7, tmps)
       status = nf90_put_var (ncid, hk_delfs2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_delfs2_id)

       call hk_repack (hk, 8, tmps)
       status = nf90_put_var (ncid, hk_hs2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hs2_id)

       call hk_repack (hk, 9, tmps)
       status = nf90_put_var (ncid, hk_phis2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_phis2_id)

       call hk_repack (hk, 10, tmps)
       status = nf90_put_var (ncid, hk_imp_colls_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_imp_colls_id)
    end if

    if (write_omega) then
       do it = 1, ntheta0
          tmp(it, :) = omega(it, :) * woutunits
       end do
       
       call c2r (tmp, ri2)
       status = nf90_put_var (ncid, omega_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, omega_id)
       
       do it = 1, ntheta0
          tmp(it, :) = omegaavg(it, :) * woutunits
       end do

       call c2r (tmp, ri2)
       status = nf90_put_var (ncid, omegaavg_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, omegaavg_id)
    end if

    if (write_scan_parameter) then
       status = nf90_put_var (ncid, &
                              current_scan_parameter_value_id, &
                              current_scan_parameter_value, &
                              start=(/nout/))
       if (status /= NF90_NOERR) &
       call netcdf_error (status, ncid, current_scan_parameter_value_id)
     end if

    status = nf90_put_var (ncid, phtot_id, phitot, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phtot_id)
# endif
  end subroutine nc_loop

  !> Outputs full not guiding center moments 
  subroutine nc_loop_fullmom (nout, ntot0, density0, upar0, tpar0, tperp0)
    implicit none
    integer, intent (in) :: nout
    complex, dimension(:,:,:), intent(in) :: ntot0, density0, upar0, tpar0, tperp0
#ifdef NETCDF
    integer, dimension (5) :: startmom
    character(*), dimension(*), parameter :: mom_dim = [character(len=7)::"ri", "kx", "ky", "species", "t"]

    startmom = [1, 1, 1, 1, nout]

    call netcdf_write_complex(ncid, "ntot0", ntot0, mom_dim, start=startmom, &
         long_name="Total perturbed density (theta=0) over time")
    call netcdf_write_complex(ncid, "density0", density0, mom_dim, start=startmom, &
         long_name="Nonadiabatic piece of perturbed density (theta=0) over time")
    call netcdf_write_complex(ncid, "upar0", upar0, mom_dim, start=startmom, &
         long_name="Upar (theta=0) over time")
    call netcdf_write_complex(ncid, "tpar0", tpar0, mom_dim, start=startmom, &
         long_name="Tpar (theta=0) over time")
    call netcdf_write_complex(ncid, "tperp0", tperp0, mom_dim, start=startmom, &
         long_name="Tperp (theta=0) over time")
#endif
  end subroutine nc_loop_fullmom

  !> Write \(\phi, A_\parallel, B_\parallel\) as functions of real space
  subroutine nc_loop_movie (nout_movie, time, yxphi, yxapar, yxbpar)
# ifdef NETCDF
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> Current movie timestep
    integer, intent (in) :: nout_movie
    !> Current simulation time
    real, intent (in) :: time
    !> \(\phi, A_\parallel, B_\parallel\) as functions of real space
    real, dimension (:,:,:), intent (in):: yxphi, yxapar, yxbpar 
# ifdef NETCDF
    integer :: status
    integer, dimension (4) :: start, count

    status = nf90_put_var (ncid_movie, time_movie_id, time, start=(/nout_movie/))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, time_movie_id)
    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout_movie
    count(1) = yxf_lo%nx
    count(2) = yxf_lo%ny
    count(3) = 2*ntgrid+1
    count(4) = 1

    if ( fphi > zero) then
       status = nf90_put_var (ncid_movie, phi_by_xmode_id, yxphi, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, phi_by_xmode_id)
    end if

    if ( fapar > zero) then
       status = nf90_put_var (ncid_movie, apar_by_xmode_id, yxapar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, apar_by_xmode_id)
    end if

    if ( fbpar > zero) then
       status = nf90_put_var (ncid_movie, bpar_by_xmode_id, yxbpar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, bpar_by_xmode_id)
    end if
# endif
  end subroutine nc_loop_movie

  !> Sync the files every nsync writes
  subroutine nc_sync(nout, nout_movie, nsync_freq)
# ifdef NETCDF
    use netcdf, only: nf90_sync
# endif
    implicit none
    !> Current output timestep
    integer, intent(in) :: nout
    !> Current output timestep for real-space fields
    integer, intent(in) :: nout_movie
    !> Timestep period at which to actually flush files
    integer, intent(in), optional :: nsync_freq
# ifdef NETCDF
    integer, parameter :: nsync_default = 10
    integer :: status, nsync

    nsync = nsync_default
    if (present(nsync_freq)) nsync = nsync_freq

    !Main file
    if(mod(nout,nsync) == 0) then
       status = nf90_sync (ncid)
       if (status /= NF90_NOERR) call netcdf_error (status)
    endif

    !Movie file
    if(my_make_movie) then
       if(mod(nout_movie,nsync) == 0) then
          status = nf90_sync (ncid_movie)       
          if (status /= NF90_NOERR) call netcdf_error (status)
       endif
    endif
# endif
  end subroutine nc_sync

  !> Write [[species:spec]] to output netCDF file
  subroutine nc_species(file_id)
#ifdef NETCDF
    use species, only: spec
    use neasyf, only: neasyf_write
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
#ifdef NETCDF
    call neasyf_write(file_id, "charge", spec%z, dim_names=["species"], &
         long_name="Charge", units="q")
    call neasyf_write(file_id, "mass", spec%mass, dim_names=["species"], &
         long_name="Atomic mass", units="m")
    call neasyf_write(file_id, "dens", spec%dens, dim_names=["species"], &
         long_name="Density", units="n_e")
    call neasyf_write(file_id, "temp", spec%temp, dim_names=["species"], &
         long_name="Temperature", units="T")
    call neasyf_write(file_id, "tprim", spec%tprim, dim_names=["species"], &
         long_name="-1/rho dT/drho")
    call neasyf_write(file_id, "fprim", spec%fprim, dim_names=["species"], &
         long_name="-1/rho dn/drho")
    call neasyf_write(file_id, "uprim", spec%uprim, dim_names=["species"], &
         long_name="-1/v_t du_par/drho")
    call neasyf_write(file_id, "uprim2", spec%uprim2, dim_names=["species"])
    call neasyf_write(file_id, "vnewk", spec%vnewk, dim_names=["species"], &
         long_name="Collisionality", units="v_t/L")
    call neasyf_write(file_id, "type_of_species", spec%type, dim_names=["species"])
#endif
  end subroutine nc_species
  
  !> Write various geometric quantities to output netCDF file
  subroutine nc_geo(file_id)
#ifdef NETCDF
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob, &
         shat, drhodpsi, eps_trapped, cdrift, cdrift0, qval, rplot, &
         zplot, aplot, rprime, zprime, aprime, bpol !, surfarea
    use run_parameters, only: beta
    use neasyf, only: neasyf_write
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
#ifdef NETCDF
    call neasyf_write(file_id, "bmag", bmag, dim_names=["theta"], &
         long_name="Magnitude of magnetic field", units="B_0")
    call neasyf_write(file_id, "bpol", bpol, dim_names=["theta"], &
         long_name="Poloidal magnetic field", units="B_0")
    call neasyf_write(file_id, "gradpar", gradpar, dim_names=["theta"], &
         long_name="Parallel derivative multiplier")
    call neasyf_write(file_id, "gbdrift", gbdrift, dim_names=["theta"], &
         long_name="Magnetic gradient drift")
    call neasyf_write(file_id, "gbdrift0", gbdrift0, dim_names=["theta"])
    call neasyf_write(file_id, "cvdrift", cvdrift, dim_names=["theta"])
    call neasyf_write(file_id, "cvdrift0", cvdrift0, dim_names=["theta"])
    call neasyf_write(file_id, "gds2", gds2, dim_names=["theta"])
    call neasyf_write(file_id, "gds21", gds21, dim_names=["theta"])
    call neasyf_write(file_id, "gds22", gds22, dim_names=["theta"])
    call neasyf_write(file_id, "gds23", gds23, dim_names=["theta"])
    call neasyf_write(file_id, "gds24", gds24, dim_names=["theta"])
    call neasyf_write(file_id, "gds24_noq", gds24_noq, dim_names=["theta"])
    call neasyf_write(file_id, "grho", grho, dim_names=["theta"])
    call neasyf_write(file_id, "jacob", jacob, dim_names=["theta"])
    call neasyf_write(file_id, "shat", shat, &
         long_name="(rho/q) dq/drho")
    call neasyf_write(file_id, "drhodpsi", drhodpsi, &
         long_name="drho/dPsi")
    call neasyf_write(file_id, "eps_trapped", eps_trapped, &
         long_name="Trapped fraction 1 - sqrt(bmin/bmax)")
    call neasyf_write(file_id, "cdrift", cdrift, dim_names=["theta"])
    call neasyf_write(file_id, "cdrift0", cdrift0, dim_names=["theta"])
    call neasyf_write(file_id, "qval", qval, &
         long_name="Local safety factor")
    call neasyf_write(file_id, "beta", beta, &
         long_name="Reference beta")
    call neasyf_write(file_id, "rplot", Rplot, dim_names=["theta"], &
         long_name="The major radius at the centre of the flux tube", units="a")
    call neasyf_write(file_id, "zplot", Zplot, dim_names=["theta"], &
         long_name="The height above the midplane at the centre of the flux tube", units="a")
    call neasyf_write(file_id, "aplot", aplot, dim_names=["theta"], &
         long_name="The toroidal angle at the centre of the flux tube", units="rad")
    call neasyf_write(file_id, "rprime", Rprime, dim_names=["theta"], &
         long_name="d/dx of the major radius at the centre of the flux tube", units="a/rho_r")
    call neasyf_write(file_id, "zprime", Zprime, dim_names=["theta"], &
         long_name="d/dx of the height above the midplane at the centre of the flux tube", units="a/rho_r")
    call neasyf_write(file_id, "aprime", aprime, dim_names=["theta"], &
         long_name="d/dx of the toroidal angle at the centre of the flux tube", units="rad/rho_r")

    ! FIXME: check if we need surfarea?
#endif
  end subroutine nc_geo

#ifdef NETCDF
  !> Wrapper around `nf90_inquire_dimension`. Aborts if any errors are
  !> detected.
  function get_dim_length(file_id, name) result(length)
    use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension, NF90_MAX_NAME
    use neasyf, only: neasyf_error

    integer :: length
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Name of the dimension
    character(len=*), intent(in) :: name

    integer :: dim_id
    ! Temporary variable for the dimension name
    character(len=NF90_MAX_NAME) :: dim_name
    ! NetCDF error status
    integer :: status

    status = nf90_inq_dimid(file_id, name, dim_id)
    call neasyf_error(status, file_id)

    status = nf90_inquire_dimension(file_id, dim_id, dim_name, length)
    call neasyf_error(status, file_id, dimid=dim_id, dim=name)
  end function get_dim_length

  subroutine write_complex_rank1(parent_id, name, values, dim_names, units, long_name, start)
    use neasyf, only: neasyf_write
    use convert, only: c2r
    !> Name of the variable
    character(len=*), intent(in) :: name
    !> NetCDF ID of the parent group/file
    integer, intent(in) :: parent_id
    !> Value of the integer to write
    complex, dimension(:), intent(in) :: values
    !> Array of dimension names, the first of which should be `"ri"`
    character(len=*), dimension(:), intent(in) :: dim_names
    !> Units of coordinate
    character(len=*), optional, intent(in) :: units
    !> Long descriptive name
    character(len=*), optional, intent(in) :: long_name
    !> Offset (one-indexed) for writing time-varying arrays. Must include the complex dimension as the first index 
    integer, dimension(:), optional, intent(in) :: start

    real, dimension(2, size(values, 1)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
  end subroutine write_complex_rank1

  subroutine write_complex_rank2(parent_id, name, values, dim_names, units, long_name, start)
    use neasyf, only: neasyf_write
    use convert, only: c2r
    !> Name of the variable
    character(len=*), intent(in) :: name
    !> NetCDF ID of the parent group/file
    integer, intent(in) :: parent_id
    !> Value of the integer to write
    complex, dimension(:, :), intent(in) :: values
    !> Array of dimension names
    character(len=*), dimension(:), intent(in) :: dim_names
    !> Units of coordinate
    character(len=*), optional, intent(in) :: units
    !> Long descriptive name
    character(len=*), optional, intent(in) :: long_name
    integer, dimension(:), optional, intent(in) :: start

    real, dimension(2, size(values, 1), size(values, 2)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
  end subroutine write_complex_rank2

  subroutine write_complex_rank3(parent_id, name, values, dim_names, units, long_name, start)
    use neasyf, only: neasyf_write
    use convert, only: c2r
    !> Name of the variable
    character(len=*), intent(in) :: name
    !> NetCDF ID of the parent group/file
    integer, intent(in) :: parent_id
    !> Value of the integer to write
    complex, dimension(:, :, :), intent(in) :: values
    !> Array of dimension names
    character(len=*), dimension(:), intent(in) :: dim_names
    !> Units of coordinate
    character(len=*), optional, intent(in) :: units
    !> Long descriptive name
    character(len=*), optional, intent(in) :: long_name
    integer, dimension(:), optional, intent(in) :: start

    real, dimension(2, size(values, 1), size(values, 2), size(values, 3)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
  end subroutine write_complex_rank3

  subroutine write_complex_rank4(parent_id, name, values, dim_names, units, long_name, start)
    use neasyf, only: neasyf_write
    use convert, only: c2r
    !> Name of the variable
    character(len=*), intent(in) :: name
    !> NetCDF ID of the parent group/file
    integer, intent(in) :: parent_id
    !> Value of the integer to write
    complex, dimension(:, :, :, :), intent(in) :: values
    !> Array of dimension names
    character(len=*), dimension(:), intent(in) :: dim_names
    !> Units of coordinate
    character(len=*), optional, intent(in) :: units
    !> Long descriptive name
    character(len=*), optional, intent(in) :: long_name
    integer, dimension(:), optional, intent(in) :: start

    real, dimension(2, size(values, 1), size(values, 2), size(values, 3), size(values, 4)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
  end subroutine write_complex_rank4
#endif
end module gs2_io
