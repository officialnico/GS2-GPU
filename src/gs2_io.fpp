!> A module for writing to netcdf.
module gs2_io

# ifdef NETCDF
  use netcdf_utils, only: kind_nf
# endif

  implicit none

  private

  public :: init_gs2_io, nc_eigenfunc, nc_final_fields, nc_final_epar
  public :: nc_final_moments, nc_final_an, nc_finish
  public :: nc_qflux, nc_vflux, nc_pflux, nc_pflux_tormom, nc_loop, nc_loop_moments
  public :: nc_loop_movie, nc_write_moments, nc_exchange
  public :: nc_loop_fullmom, nc_loop_dist, nc_loop_sym, nc_loop_corr, nc_loop_corr_extend
  public :: nc_loop_partsym_tormom, nc_sync, get_netcdf_file_id, get_netcdf_movie_file_id
  public :: define_dims, nc_norms, nc_species, nc_geo, save_input, nc_grids_mymovie, get_dim_length
  public :: nc_write_omega, nc_write_field, nc_write_heating, starts, nc_write_ql_metric
  public :: netcdf_write_complex
  public :: flux_dims, fluxk_dims, final_field_dims, sym_dims, omega_dims, mode_dims, mom_dims, mom_t_dims
  public :: species_dim, x_dim, kx_dim, y_dim, ky_dim, theta_dim, theta_extended_dim, time_dim
  public :: time_big_dim, time_movie_dim, complex_dim, vpa_dim, lambda_dim, egrid_dim, sign_dim
  public :: nheat_dim, dim_2, dim_3, dim_4, dim_5

  !> Maximum length of dimension names. This is required to get around
  !> Fortran hating ragged arrays of characters
  integer, parameter :: max_dim_name_length = 9

  !> Name of species dimension
  character(max_dim_name_length), parameter :: species_dim = "species"
  !> Name of x dimension
  character(max_dim_name_length), parameter :: x_dim = "x"
  !> Name of kx dimension
  character(max_dim_name_length), parameter :: kx_dim = "kx"
  !> Name of y dimension
  character(max_dim_name_length), parameter :: y_dim = "y"
  !> Name of ky dimension
  character(max_dim_name_length), parameter :: ky_dim = "ky"
  !> Name of theta dimension
  character(max_dim_name_length), parameter :: theta_dim = "theta"
  !> Name of extended theta dimension
  character(max_dim_name_length), parameter :: theta_extended_dim = "theta_ext"
  !> Name of time dimension
  character(max_dim_name_length), parameter :: time_dim = "t"
  !> Name of "big time" dimension, used for some correlation functions
  character(max_dim_name_length), parameter :: time_big_dim = "tb"
  !> Name of time dimension for real space "movies"
  character(max_dim_name_length), parameter :: time_movie_dim = "tm"
  !> Name of complex dimension
  character(max_dim_name_length), parameter :: complex_dim = "ri"
  !> Name of vpa dimension
  character(max_dim_name_length), parameter :: vpa_dim = "vpa"
  !> Name of lambda dimension
  character(max_dim_name_length), parameter :: lambda_dim = "lambda"
  !> Name of egrid dimension
  character(max_dim_name_length), parameter :: egrid_dim = "egrid"
  !> Name of sign dimension
  character(max_dim_name_length), parameter :: sign_dim = "sign"
  !> Name of nheat dimension
  character(max_dim_name_length), parameter :: nheat_dim = "nheat"
  !> Name of generic dimension for arrays of length 2
  character(max_dim_name_length), parameter :: dim_2 = "2"
  !> Name of generic dimension for arrays of length 3
  character(max_dim_name_length), parameter :: dim_3 = "3"
  !> Name of generic dimension for arrays of length 4
  character(max_dim_name_length), parameter :: dim_4 = "4"
  !> Name of generic dimension for arrays of length 5
  character(max_dim_name_length), parameter :: dim_5 = "5"

  !> Names of dimensions used by flux variables
  character(*), dimension(*), parameter :: flux_dims = [species_dim, time_dim]
  !> Names of dimensions used by flux-by-mode variables
  character(*), dimension(*), parameter :: fluxk_dims = [kx_dim, ky_dim, species_dim, time_dim]
  !> Names of dimensions used by snapshots of complex fields
  character(*), dimension(*), parameter :: final_field_dims = [complex_dim, theta_dim, kx_dim, ky_dim]
  !> Names of dimensions used by certain fluxes
  character(*), dimension(*), parameter :: sym_dims = [theta_dim, vpa_dim, species_dim, time_dim]
  !> Names of dimensions used by time-dependent complex fields
  character(*), dimension(*), parameter :: omega_dims = [complex_dim, kx_dim, ky_dim, time_dim]
  !> Names of dimensions used by time-dependent real fields
  character(*), dimension(*), parameter :: mode_dims = [kx_dim, ky_dim, time_dim]
  !> Names of dimensions used by moments
  character(*), dimension(*), parameter :: mom_dims = [complex_dim, kx_dim, ky_dim, species_dim, time_dim]
  !> Names of dimensions used by theta-dependent moments
  character(*), dimension(*), parameter :: mom_t_dims = [complex_dim, theta_dim, kx_dim, ky_dim, species_dim, time_dim]

# ifdef NETCDF
  integer (kind_nf) :: ncid
  integer (kind_nf) :: char10_dim

  ! added by EAB 03/05/04 for movies
  logical :: io_write_corr_extend
  integer :: ncid_movie = -1
  logical :: write_apar_t, write_phi_t, write_bpar_t ! Should the fields be written out every nwrite?

#else
  integer, parameter :: ncid = -1
  integer, parameter :: ncid_movie = -1
#endif
  !> Write a `complex` array to netcdf
  !>
  !> Converts the `complex` array to a `real` array with an extra dimension
  interface netcdf_write_complex
    module procedure write_complex_scalar
    module procedure write_complex_rank1, write_complex_rank2
    module procedure write_complex_rank3, write_complex_rank4
  end interface netcdf_write_complex

  !> Value to disable netCDF compression
  integer, parameter :: compression_off = 0
  !> Level of compression if enabled via
  !> `GK_NETCDF_DEFAULT_COMPRESSION_ON` macro. Must be a value between
  !> 1 (faster, less compression) and 9 (slower, more compression), or
  !> 0 (off)
  integer, parameter :: default_compression_level = 1

contains

  !> Provides read-only access to the netCDF file ID
  integer function get_netcdf_file_id()
    get_netcdf_file_id = ncid
  end function get_netcdf_file_id

  !> Provides read-only access to the netCDF movie file ID
  integer function get_netcdf_movie_file_id()
    get_netcdf_movie_file_id = ncid_movie
  end function get_netcdf_movie_file_id

  !> Return an array of 1-index offsets for writing a time-dependent
  !> field to netCDF:
  !>
  !>     starts(2, nout) == [1, nout]
  !>     starts(3, nout) == [1, 1, nout]
  function starts(rank, nout)
    !> Rank of the array being written
    integer, intent(in) :: rank
    !> Current time index
    integer, intent(in) :: nout
    integer, dimension(rank) :: starts
    starts(1:rank-1) = 1
    starts(rank) = nout
  end function starts

  !> Create and initialise the output netCDF files
  !>
  !> Create the main and movie (real-space fields) netCDF files,
  !> creating the dimensions and main variables
  subroutine init_gs2_io (append_old, make_movie, write_correlation_extend, &
      write_phi_over_time, write_apar_over_time, write_bpar_over_time, nout, ob_midplane, header)
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
    use constants, only: run_name_size
    use gs2_metadata, only: create_metadata
    use run_parameters, only: user_comments
    use neasyf, only: neasyf_open, neasyf_default_compression
# endif
    implicit none
    !> See [[gs2_diagnostics_knobs]] for details of these inputs
    logical, intent(in) :: append_old
    logical, intent(inout) :: write_correlation_extend
    logical, intent(in) :: make_movie
    logical, intent(in) :: write_phi_over_time, write_apar_over_time, write_bpar_over_time
    !> If true, write out moments at the outboard midplane only,
    !> otherwise write them out along the entire flux surface
    logical, intent(in), optional :: ob_midplane
    integer, intent(out) :: nout
    !> Header for files with build and run information
    type(standard_header_type), intent(in) :: header

# ifdef NETCDF
    logical :: accelerated
    character(run_name_size) :: filename
    logical :: previous_file_exists, appending

    !################
    !# Update internal flags

    write_phi_t = write_phi_over_time
    write_apar_t = write_apar_over_time
    write_bpar_t = write_bpar_over_time

    io_write_corr_extend = write_correlation_extend

    !----------------

#ifdef GK_NETCDF_DEFAULT_COMPRESSION_ON
    neasyf_default_compression = default_compression_level
#else
    neasyf_default_compression = compression_off
#endif

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
    if (make_movie) then
      call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
       ! only proc0 opens the file:
       if (proc0) then
         ncid_movie = neasyf_open(trim(run_name)//'.movie.nc', "w")
         call create_metadata(ncid_movie, 'GS2 Simulation x,y,theta Data for Movies', header, user_comments)
       end if
    endif

    if (proc0) then
        !Setup main file by creating netCDF file
        call define_dims (ncid, io_write_corr_extend)
        call define_vars (write_correlation_extend, ob_midplane=ob_midplane)

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

        if (make_movie) then
          call nc_grids_mymovie(ncid_movie)
        end if
    end if
# endif
  end subroutine init_gs2_io

  !> Ensure the netCDF file contains all the dimensions and grids,
  !> creating them if necessary
  !>
  !> Dimension variable IDs are module level variables
  subroutine define_dims (file_id, create_extended_correlation)
# ifdef NETCDF
    use file_utils, only: num_input_lines
    use kt_grids, only: aky, akx, theta0, jtwist_out, ntheta0, naky, nx, ny
    use theta_grid, only: theta, ntgrid, nperiod, ntheta, nbset, bset
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
    !> Create the dimensions for the extended correlation diagnostics
    logical, intent(in) :: create_extended_correlation
# ifdef NETCDF
    integer :: nmesh

    call neasyf_dim(file_id, trim(kx_dim), values=akx, &
         long_name="Wavenumber perpendicular to flux surface", units="1/rho_r")
    call neasyf_dim(file_id, trim(ky_dim), values=aky, &
         long_name="Wavenumber in direction of grad alpha", units="1/rho_r")
    call neasyf_dim(file_id, trim(theta_dim), values=theta, long_name="Poloidal angle", units="rad")
    call neasyf_dim(file_id, trim(lambda_dim), values=al, long_name="Energy/magnetic moment", units="1/(2 B_a)")
    call neasyf_dim(file_id, trim(species_dim), dim_size=nspec)

    ! Energy coordinate is a function of species, but the dimension must be 1D
    call neasyf_dim(file_id, trim(egrid_dim), dim_size=negrid, long_name="See 'energy' for energy grid values")
    call neasyf_write(file_id, "energy", energy, dim_names=[egrid_dim, species_dim], &
         long_name="Energy grid values")

    call neasyf_dim(file_id, trim(sign_dim), dim_size=2)
    call neasyf_dim(file_id, trim(time_dim), unlimited=.true., long_name="Time", units="L/vt")
    call neasyf_dim(file_id, trim(complex_dim), dim_size=2, &
         long_name="Complex components", units="(real, imaginary)")
    call neasyf_dim(file_id, trim(nheat_dim), dim_size=7)
    call neasyf_dim(file_id, trim(vpa_dim), dim_size=negrid*nlambda)

    ! Dimensions for various string variables
    call neasyf_dim(file_id, "char5", dim_size=5)
    call neasyf_dim(file_id, "char10", dim_size=10, dimid=char10_dim)
    call neasyf_dim(file_id, "char200", dim_size=200)

    ! netCDF interprets zero-sized dimensions as unlimited, which is not what we want
    call neasyf_dim(file_id, "nlines", dim_size=max(1, num_input_lines), long_name="Input file line number")

    ! Dimensions for non-physical matrices etc
    call neasyf_dim(file_id, trim(dim_2), dim_size=2, long_name="Generic dim")
    call neasyf_dim(file_id, trim(dim_3), dim_size=3, long_name="Generic dim")
    call neasyf_dim(file_id, trim(dim_4), dim_size=4, long_name="Generic dim")
    call neasyf_dim(file_id, trim(dim_5), dim_size=5, long_name="Generic dim")

    call neasyf_dim(file_id, "bset", values=bset, long_name="Bmag at bounce points", units="Bref")

    if (create_extended_correlation) then
       call neasyf_dim(file_id, trim(theta_extended_dim), dim_size=(2*ntgrid+1)*((ntheta0-1)/jtwist_out+1))
       call neasyf_dim(file_id, trim(time_big_dim), unlimited=.true.)
    end if

    call neasyf_write(file_id, "layout", layout, long_name="Decomposition layout for g_lo")
    call neasyf_write(file_id, "nproc", nproc, long_name="Number of processors")
    call neasyf_write(file_id, "ntheta_tot", 2*ntgrid+1, long_name="Total number of theta points")
    call neasyf_write(file_id, "nkx", ntheta0, long_name="Number of kx points")
    call neasyf_write(file_id, "nky", naky, long_name="Number of ky points")
    call neasyf_write(file_id, "nspecies", nspec, long_name="Number of species")
    call neasyf_write(file_id, "nenergy", negrid, long_name="Number of energy points")
    call neasyf_write(file_id, "nlambda", nlambda, long_name="Number of lambda points")

    call neasyf_write(file_id, "nperiod", nperiod, long_name="Number of 2 pi theta sections in parallel domain chunk")
    call neasyf_write(file_id, "ntheta", ntheta, long_name="Number of theta points per 2 pi section")
    call neasyf_write(file_id, "nbset", nbset, long_name="Number of bounce points")

    call neasyf_write(file_id, "theta0", theta0, dim_names=[kx_dim, ky_dim], &
         long_name="Theta_0, poloidal angle where mode radial wavenumber == 0", units="rad")

    if (nonlin) then
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
    else
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
    end if

    call neasyf_write(file_id, "nmesh", nmesh, long_name="Total number of mesh points")

# endif
  end subroutine define_dims

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
  subroutine nc_grids_mymovie(file_id)
#ifdef NETCDF
    use constants, only: pi
    use kt_grids, only: akx, aky
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid, theta
    use neasyf, only: neasyf_dim, neasyf_write
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
# ifdef NETCDF
    real, dimension(:), allocatable :: x,y
    integer :: it, ik

    allocate(x(yxf_lo%nx),y(yxf_lo%ny))

    do it = 1, yxf_lo%nx
      x(it) = 2.0*pi/akx(2)*(-0.5+ real(it-1)/real(yxf_lo%nx))
    end do
    do ik = 1, yxf_lo%ny
      y(ik) = 2.0*pi/aky(2)*(-0.5+ real(ik-1)/real(yxf_lo%ny))
    end do

    call neasyf_write(file_id, "nx", yxf_lo%nx)
    call neasyf_write(file_id, "ny", yxf_lo%ny)
    call neasyf_write(file_id, "ntheta", 2*ntgrid+1)
    call neasyf_dim(file_id, x_dim, x, long_name="x / rho")
    call neasyf_dim(file_id, y_dim, y, long_name="y / rho")
    call neasyf_dim(file_id, trim(theta_dim), theta)
    call neasyf_dim(file_id, trim(time_movie_dim), unlimited=.true., long_name="Time", units="L/vt")
    deallocate(x,y)
# endif
  end subroutine nc_grids_mymovie

  !> Close the output netCDF files
  subroutine nc_finish
# ifdef NETCDF
    use mp, only: proc0
    use neasyf, only: neasyf_close
    implicit none

    if (proc0) then
       call neasyf_close(ncid)
       if (ncid_movie >= 0) then
          call neasyf_close(ncid_movie)
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

    allocate(input_file_array(max(1, num_input_lines)))
    input_file_array(1) = ""

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
  subroutine define_vars (write_correlation_extend, ob_midplane)
# ifdef NETCDF
    use file_utils, only: error_unit
    use netcdf, only : NF90_CHAR, nf90_inq_libvers
    use netcdf_utils, only: get_netcdf_dim_length
    use netcdf_utils, only: ensure_netcdf_att_exists, ensure_netcdf_var_exists
# endif
    implicit none
    logical, intent(inout) :: write_correlation_extend
    logical, optional, intent(in) :: ob_midplane
# ifdef NETCDF
    character (20) :: datestamp, timestamp, timezone
    integer :: dimsize(5)
    real    :: total
    logical :: obmid=.false.
    integer :: code_id

    !Handle optional input
    if (present(ob_midplane)) then
       obmid=ob_midplane
    endif

    ! HJL This is a hack to stop the netcdf diagnostic output failing when the variable sizes are too big
    ! The netcdf limit is 2^31-4 bytes, however this limit still causes a problem so 178180992 is used
    ! This limit was found by running the code with different sizes, it is not exact but this is a hack
    ! and is not worth pursuing at the moment
    !<DD>Added protective if statement as theta_extended_dim is only defined if io_write_corr_extend
    if(io_write_corr_extend)then
       dimsize(1) = get_dim_length(ncid, complex_dim)
       dimsize(2) = get_dim_length(ncid, theta_extended_dim)
       dimsize(3) = get_dim_length(ncid, kx_dim)
       dimsize(4) = get_dim_length(ncid, ky_dim)
       dimsize(5) = get_dim_length(ncid, time_big_dim)
       total = product(real(dimsize(1:5)))
       if(total .gt. 178180992) then ! number of reals allowed in netcdf
          write(error_unit(), *) 'WARNING: phi_corr is larger than maximum netcdf size, removing.'
          write(error_unit(), *) 'Try increasing nwrite_mult to avoid this.'
          io_write_corr_extend = .false.
          write_correlation_extend = .false.
       endif
    endif
    ! HJL End of hack

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
# endif
  end subroutine define_vars

  !> Write the normalised fields to netCDF
  subroutine nc_eigenfunc(file_id, phase)
#ifdef NETCDF
    use convert, only: c2r
    use run_parameters, only: has_phi, has_apar, has_bpar
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

    call netcdf_write_complex(file_id, "phase", phase, [complex_dim, kx_dim, ky_dim], long_name="Normalising phase")

    if (has_phi) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = phi(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "phi_norm", tmp, final_field_dims)
    end if

    if (has_apar) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = apar(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "apar_norm", tmp, final_field_dims)
    end if

    if (has_bpar) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = bpar(ig,:,:)/phase(:,:)
       end do
       call netcdf_write_complex(file_id, "bpar_norm", tmp, final_field_dims)
    end if
# endif
  end subroutine nc_eigenfunc

  !> Write various moments to netcdf
  subroutine nc_write_moments (file_id, nout, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1, ob_midplane)
    use theta_grid, only: ntgrid
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
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

    if (present(ob_midplane)) then
       if (ob_midplane ) then
          call netcdf_write_complex(file_id, "ntot0", ntot(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Total perturbed density (theta=0) over time")
          call netcdf_write_complex(file_id, "density0", density(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Nonadiabatic piece of perturbed density (theta=0) over time")
          call netcdf_write_complex(file_id, "upar0", upar(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Upar (theta=0) over time")
          call netcdf_write_complex(file_id, "tpar0", tpar(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Tpar (theta=0) over time")
          call netcdf_write_complex(file_id, "tperp0", tperp(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Tperp (theta=0) over time")
          call netcdf_write_complex(file_id, "qparflux0", qparflux(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Qparflux (theta=0) over time")
          call netcdf_write_complex(file_id, "pperpj10", pperpj1(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Pperpj1 (theta=0) over time")
          call netcdf_write_complex(file_id, "qpperpj10", qpperpj1(0, :, :, :), mom_dims, start=starts(5, nout), &
               long_name="Qpperpj1 (theta=0) over time")
          return
       endif
    endif

    call netcdf_write_complex(file_id, 'ntot_t', ntot, mom_t_dims, start=starts(6, nout), &
          long_name="Total perturbed density over time")
    call netcdf_write_complex(file_id, 'density_t', density, mom_t_dims, start=starts(6, nout), &
          long_name="Nonadiabatic piece of perturbed density over time")
    call netcdf_write_complex(file_id, 'upar_t', upar, mom_t_dims, start=starts(6, nout), &
          long_name="Upar over time")
    call netcdf_write_complex(file_id, 'tpar_t', tpar, mom_t_dims, start=starts(6, nout), &
          long_name="Tpar over time")
    call netcdf_write_complex(file_id, 'tperp_t', tperp, mom_t_dims, start=starts(6, nout), &
          long_name="Tperp over time")
    call netcdf_write_complex(file_id, 'qparflux_t', qparflux, mom_t_dims, start=starts(6, nout), &
          long_name="Qparflux over time")
    call netcdf_write_complex(file_id, 'pperpj1_t', pperpj1, mom_t_dims, start=starts(6, nout), &
          long_name="Pperpj1 over time")
    call netcdf_write_complex(file_id, 'qpperpj1_t', qpperpj1, mom_t_dims, start=starts(6, nout), &
          long_name="Qpperpj1 over time")

# endif
  end subroutine nc_write_moments

  !> Write the fields to netCDF, overwriting existing values
  subroutine nc_final_fields(file_id)
#ifdef NETCDF
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: has_phi, has_apar, has_bpar
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
#ifdef NETCDF
    if (has_phi) then
       call netcdf_write_complex(file_id, "phi", phinew, final_field_dims, &
            long_name="Electrostatic Potential", units="T_ref/q_ref rho_ref/a_ref")
    end if

    if (has_apar) then
       call netcdf_write_complex(file_id, "apar", aparnew, final_field_dims, &
            long_name="Parallel Magnetic Potential")
    end if

    if (has_bpar) then
       call netcdf_write_complex(file_id, "bpar", bparnew, final_field_dims, &
            long_name="delta B Parallel")
    end if
# endif
  end subroutine nc_final_fields

  !> Write \(E_\parallel\) to netCDF, overwriting existing values
  subroutine nc_final_epar (file_id, epar)
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    complex, dimension (:,:,:), intent (in) :: epar
# ifdef NETCDF
    character(*), dimension(*), parameter :: dim_names = [complex_dim, theta_dim, kx_dim, ky_dim]

    call netcdf_write_complex(file_id, "epar", epar, dim_names)
# endif
  end subroutine nc_final_epar

  !> Write the various low-order moments to the output netCDF file,
  !> overwriting existing values
  subroutine nc_final_moments (file_id, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
    complex, dimension (:,:,:,:), intent (in) :: qparflux, pperpj1, qpperpj1
    character(*), dimension(*), parameter :: dim_names = [complex_dim, theta_dim, kx_dim, ky_dim, species_dim]

#ifdef NETCDF
    call netcdf_write_complex(file_id, 'ntot', ntot, dim_names)
    call netcdf_write_complex(file_id, 'density', density, dim_names)
    call netcdf_write_complex(file_id, 'upar', upar, dim_names)
    call netcdf_write_complex(file_id, 'tpar', tpar, dim_names)
    call netcdf_write_complex(file_id, 'tperp', tperp, dim_names)
    call netcdf_write_complex(file_id, 'qparflux', qparflux, dim_names)
    call netcdf_write_complex(file_id, 'pperpj1', pperpj1, dim_names)
    call netcdf_write_complex(file_id, 'qpperpj1', qpperpj1, dim_names)
#endif
  end subroutine nc_final_moments

  !> Write moments as: volume averaged, field line averaged, and evaulated at theta=0
  subroutine nc_loop_moments (file_id, nout, ntot2, ntot2_by_mode, ntot20, ntot20_by_mode, &
       phi00, ntot00, density00, upar00, tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: ntot2, ntot20
    real, dimension (:,:,:), intent (in) :: ntot2_by_mode, ntot20_by_mode
    real, dimension (:,:,:), intent (in) :: tpar2_by_mode, tperp2_by_mode
    complex, dimension (:), intent (in) :: phi00
    complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
#ifdef NETCDF
    character(*), dimension(*), parameter :: loop_phi_dim = [complex_dim, kx_dim, time_dim]
    character(*), dimension(*), parameter :: loop_mom_dim = [complex_dim, kx_dim, species_dim, time_dim]
    
    ! Volume averaged function of theta
    ! <DD>Should we only do this if .not.obmid?
    call neasyf_write (file_id, "ntot2", ntot2, dim_names=flux_dims, start=starts(2, nout))
    call neasyf_write (file_id, "ntot2_by_mode", ntot2_by_mode, dim_names=fluxk_dims, start=starts(4, nout))
    call neasyf_write (file_id, "tpar2_by_mode", tpar2_by_mode, dim_names=fluxk_dims, start=starts(4, nout))
    call neasyf_write (file_id, "tperp2_by_mode", tperp2_by_mode, dim_names=fluxk_dims, start=starts(4, nout))

    ! "Volume averaged" data from outboard midplane (really just evaluated at theta=0
    ! and then weighted with Fourier coeff factor (typically 0.5).
    ! <DD>Should we only do this if obmid?
    call neasyf_write (file_id, "ntot20", ntot20, dim_names=flux_dims, start=starts(2, nout), &
         long_name='Density**2 at theta=0')
    call neasyf_write (file_id, "ntot20_by_mode", ntot20_by_mode, dim_names=fluxk_dims, start=starts(4, nout))

    ! Fieldline averaged data, note this includes phi which is not a moment
    call netcdf_write_complex(file_id, "phi00", phi00, loop_phi_dim, start=starts(3, nout))
    call netcdf_write_complex(file_id, "ntot00", ntot00, loop_mom_dim, start=starts(4, nout))
    call netcdf_write_complex(file_id, "density00", density00, loop_mom_dim, start=starts(4, nout))
    call netcdf_write_complex(file_id, "upar00", upar00, loop_mom_dim, start=starts(4, nout))
    call netcdf_write_complex(file_id, "tpar00", tpar00, loop_mom_dim, start=starts(4, nout))
    call netcdf_write_complex(file_id, "tperp00", tperp00, loop_mom_dim, start=starts(4, nout))
#endif
  end subroutine nc_loop_moments

  !> Write the right-hand sides of the field equations to netCDF,
  !> overwriting existing values
  subroutine nc_final_an (file_id, antot, antota, antotp)
# ifdef NETCDF
    use run_parameters, only: has_phi, has_apar, has_bpar
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    complex, dimension (:,:,:), intent(in) :: antot, antota, antotp
# ifdef NETCDF
    if (has_phi) then
       call netcdf_write_complex(file_id, "antot", antot, final_field_dims, &
            long_name="Velocity moment of g used in calculating phi")
    end if

    if (has_apar) then
       call netcdf_write_complex(file_id, "antota", antota, final_field_dims, &
            long_name="Velocity space moment of g used in calculating apar")
    end if

    if (has_bpar) then
       call netcdf_write_complex(file_id, "antotp", antotp, final_field_dims, &
            long_name="Velocity space moment of g used in calculating bpar")
    end if
# endif
  end subroutine nc_final_an

  !> Write the energy exchange to netCDF
  subroutine nc_exchange (file_id, nout, exchange, exchange_avg)
# ifdef NETCDF
    use run_parameters, only: has_phi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    !> Energy exchange
    real, dimension (:,:,:), intent (in) :: exchange
    !> Volume-averaged energy exchange
    real, dimension (:), intent (in) :: exchange_avg
# ifdef NETCDF
    if (has_phi) then
       call neasyf_write (file_id, "es_exchange", exchange_avg, &
            dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write (file_id, "es_exchange_by_k", exchange, &
            dim_names=fluxk_dims, start=starts(4, nout))
    end if
# endif
  end subroutine nc_exchange

  !> Write various heat fluxes to netCDF
  subroutine nc_qflux (file_id, nout, qheat, qmheat, qbheat, &
       heat_par,  mheat_par,  bheat_par, &
       heat_perp, mheat_perp, bheat_perp, &
       heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot)
# ifdef NETCDF
    use run_parameters, only: has_phi, has_apar, has_bpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: qheat, qmheat, qbheat
    real, dimension (:), intent (in) :: heat_par, mheat_par, bheat_par
    real, dimension (:), intent (in) :: heat_perp, mheat_perp, bheat_perp
    real, dimension (:), intent (in) :: heat_fluxes, mheat_fluxes, bheat_fluxes
    real, dimension (:,:), intent (in) :: x_qmflux
    real, intent (in) :: hflux_tot
# ifdef NETCDF
    character(*), dimension(*), parameter :: fluxx_dim = [kx_dim, species_dim, time_dim]

    if (has_phi) then
       call neasyf_write(file_id, "es_heat_flux", heat_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_heat_par", heat_par, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_heat_perp", heat_perp, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_heat_by_k", qheat, dim_names=fluxk_dims, start=starts(4, nout))
    end if

    if (has_apar) then
       call neasyf_write(file_id, "apar_heat_flux", mheat_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "apar_heat_par", mheat_par, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "apar_heat_perp", mheat_perp, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "apar_heat_by_k", qmheat, dim_names=fluxk_dims, start=starts(4, nout))
       call neasyf_write(file_id, "apar_heat_by_x", x_qmflux, dim_names=fluxx_dim, start=starts(3, nout))
    end if

    if (has_bpar) then
       call neasyf_write(file_id, "bpar_heat_flux", bheat_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "bpar_heat_par", bheat_par, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "bpar_heat_perp", bheat_perp, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "bpar_heat_by_k", qbheat, dim_names=fluxk_dims, start=starts(4, nout))
    end if
    
    call neasyf_write(file_id, "hflux_tot", hflux_tot, dim_names=[time_dim], start=[nout])
# endif
  end subroutine nc_qflux

  !> FIXME : Add documentation  
  subroutine nc_pflux (file_id, nout, pflux, pmflux, pbflux, &
       part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
# ifdef NETCDF
    use run_parameters, only: has_phi, has_apar, has_bpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
    real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
    real, intent (in) :: zflux_tot
# ifdef NETCDF

    if (has_phi) then
       call neasyf_write(file_id, "es_part_flux", part_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_part_by_k", pflux, dim_names=fluxk_dims, start=starts(4, nout))
    end if

    if (has_apar) then
       call neasyf_write(file_id, "apar_part_flux", mpart_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "apar_part_by_k", pmflux, dim_names=fluxk_dims, start=starts(4, nout))
    end if

    if (has_bpar) then
       call neasyf_write(file_id, "bpar_part_flux", bpart_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "bpar_part_by_k", pbflux, dim_names=fluxk_dims, start=starts(4, nout))
    end if
    
    call neasyf_write(file_id, "zflux_tot", zflux_tot, dim_names=[time_dim], start=[nout])
# endif
  end subroutine nc_pflux

  !> Write toroidal angular momentum flux to netCDF
  subroutine nc_pflux_tormom (file_id, nout, pflux_tormom, part_tormom_fluxes)
# ifdef NETCDF
    use run_parameters, only: has_phi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux_tormom
    real, dimension(:), intent (in) :: part_tormom_fluxes
# ifdef NETCDF

    if (has_phi) then
       call neasyf_write(file_id, "es_part_tormom_flux", part_tormom_fluxes, &
            dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_part_tormom_by_k", pflux_tormom, &
            dim_names=fluxk_dims, start=starts(4, nout))
    end if

# endif
  end subroutine nc_pflux_tormom

  !> FIXME : Add documentation  
  subroutine nc_vflux (file_id, nout, vflux, vmflux, vbflux, &
       mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
       vflux_par, vflux_perp, vflux0, vflux1)
# ifdef NETCDF
    use run_parameters, only: has_phi, has_apar, has_bpar
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
    real, dimension (:,:,:), intent (in) :: vflux_par, vflux_perp
    real, dimension (:,:,:), intent (in) :: vflux0, vflux1
    real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
    real, intent (in) :: vflux_tot
# ifdef NETCDF

    if (has_phi) then
       call neasyf_write(file_id, "es_mom_flux", mom_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "es_mom_by_k", vflux, dim_names=fluxk_dims, start=starts(4, nout))
       call neasyf_write(file_id, "es_parmom_by_k", vflux_par, dim_names=fluxk_dims, start=starts(4, nout))
       call neasyf_write(file_id, "es_perpmom_by_k", vflux_perp, dim_names=fluxk_dims, start=starts(4, nout))
       call neasyf_write(file_id, "es_mom0_by_k", vflux0, dim_names=fluxk_dims, start=starts(4, nout))
       call neasyf_write(file_id, "es_mom1_by_k", vflux1, dim_names=fluxk_dims, start=starts(4, nout))
    end if

    if (has_apar) then
       call neasyf_write(file_id, "apar_mom_flux", mmom_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "apar_mom_by_k", vmflux, dim_names=fluxk_dims, start=starts(4, nout))
    end if

    if (has_bpar) then
       call neasyf_write(file_id, "bpar_mom_flux", bmom_fluxes, dim_names=flux_dims, start=starts(2, nout))
       call neasyf_write(file_id, "bpar_mom_by_k", vbflux, dim_names=fluxk_dims, start=starts(4, nout))
    end if
    
    call neasyf_write(file_id, "vflux_tot", vflux_tot, dim_names=[time_dim], start=[nout])
# endif
  end subroutine nc_vflux

  !> Output the poloidal distributions of the fluxes of particles, perpendicular
  !> momentum, parallel momentum, and energy ! JRB
  subroutine nc_loop_dist (file_id, nout, part_fluxes_dist, mom_fluxes_par_dist, mom_fluxes_perp_dist, heat_fluxes_dist, phi_dist)
    use run_parameters, only: has_phi
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use neasyf, only: neasyf_write
# endif
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:), intent (in) :: part_fluxes_dist, mom_fluxes_par_dist, mom_fluxes_perp_dist, heat_fluxes_dist
    real, dimension (2,-ntgrid:ntgrid), intent (in) :: phi_dist
# ifdef NETCDF
    character(*), dimension(*), parameter :: flux_dist_dim = [theta_dim, species_dim, time_dim]
    character(*), dimension(*), parameter :: phi_dist_dim = [complex_dim, theta_dim, time_dim]

    if (has_phi) then
       call neasyf_write(file_id, "es_part_flux_dist", part_fluxes_dist, &
            dim_names=flux_dist_dim, start=starts(3, nout), &
            long_name= "Particle flux density, see Ball et al PPCF 58 045023 2016")
       call neasyf_write(file_id, "es_mom_flux_par_dist", mom_fluxes_par_dist, &
            dim_names=flux_dist_dim, start=starts(3, nout), &
            long_name= "Parallel momentum flux density, see Ball et al PPCF 58 045023 2016")
       call neasyf_write(file_id, "es_mom_flux_perp_dist", mom_fluxes_perp_dist, &
            dim_names=flux_dist_dim, start=starts(3, nout), &
            long_name= "Perpendicular momentum flux density, see Ball et al PPCF 58 045023 2016")
       call neasyf_write(file_id, "es_heat_flux_dist", heat_fluxes_dist, &
            dim_names=flux_dist_dim, start=starts(3, nout), &
            long_name= "Heat flux density, see Ball et al PPCF 58 045023 2016")
       call neasyf_write(file_id, "es_phi_dist", phi_dist, &
            dim_names=phi_dist_dim, start=starts(3, nout))
    end if

# endif
  end subroutine nc_loop_dist
  
  !> Write the momentum flux as a function of \((v_\parallel, \theta, t)\) to netCDF
  subroutine nc_loop_sym (file_id, nout, vflx_sym)
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use run_parameters, only: has_phi
    use neasyf, only: neasyf_write
    implicit none
# endif
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:,:), intent (in) :: vflx_sym
# ifdef NETCDF
    
    if (has_phi) then
       call neasyf_write(file_id, "es_mom_sym", vflx_sym, &
            dim_names=sym_dims, start=starts(4, nout), &
            long_name="Momentum flux density as a function theta and vspace, used for looking at the effect of asymmetry, &
            & see Parra et al POP 18 062501 2011")
    end if

# endif
  end subroutine nc_loop_sym
  
  !> Write the particle flux contribution to toroidal momentum flux to netCDF
  subroutine nc_loop_partsym_tormom (file_id, nout, partflx_sym)
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use run_parameters, only: has_phi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    real, dimension (-ntgrid:,:,:), intent (in) :: partflx_sym
# ifdef NETCDF

    if (has_phi) then
       call neasyf_write(file_id, "es_part_tormom_sym", partflx_sym, &
            dim_names=sym_dims, start=starts(4, nout), &
            long_name="Particle flux density as a function theta and vspace, used for looking at the effect of asymmetry, &
            & see Parra et al POP 18 062501 2011 and ask Jung-Pyo Lee")
    end if

# endif
  end subroutine nc_loop_partsym_tormom

  !> Write the correlation function over the physical domain to netCDF
  subroutine nc_loop_corr (file_id, nout, phi_2pi_corr)
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use run_parameters, only: has_phi
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    !> The correlation function over the physical domain
    complex, dimension (-ntgrid:,:), intent (in) :: phi_2pi_corr
# ifdef NETCDF

    if (has_phi) then
       call netcdf_write_complex(file_id, "phi_corr_2pi", phi_2pi_corr, &
            dim_names=[complex_dim, theta_dim, ky_dim, time_dim], start=starts(4, nout), &
            long_name="2 point correlation function calculated from the electric potential")
    end if

# endif
  end subroutine nc_loop_corr

  !> Write the correlation function over the extended domain to netCDF
  subroutine nc_loop_corr_extend (file_id, nout_big, time, phi_corr, phi2_extend)
# ifdef NETCDF
    use run_parameters, only: has_phi
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> The current output timestep for "large" diagnostics
    integer, intent (in) :: nout_big
    !> Current simulation time
    real, intent(in) :: time
    real, dimension (:,:,:), intent (in) :: phi2_extend
    complex, dimension (:,:,:), intent (in) :: phi_corr
# ifdef NETCDF
    character(*), dimension(*), parameter :: phi_corr_dims = &
         [complex_dim, theta_extended_dim, kx_dim, ky_dim, time_big_dim]
    character(*), dimension(*), parameter :: phi2_extend_dims = &
         [theta_extended_dim, kx_dim, ky_dim, time_big_dim]

    if (has_phi) then
      call neasyf_write(file_id, trim(time_big_dim), time, dim_names=[time_big_dim], start=[nout_big])

       call netcdf_write_complex(file_id, "phi_corr", phi_corr, &
            dim_names=phi_corr_dims, start=starts(5, nout_big), &
            long_name="2 point correlation function calculated from the electric potential, time averaged")
       call neasyf_write(file_id, "phi2_extend", phi2_extend, &
            dim_names=phi2_extend_dims, start=starts(4, nout_big))
    end if

# endif
  end subroutine nc_loop_corr_extend

  !> Write various fields, heating, and frequency diagnostics to netCDF
  subroutine nc_loop (file_id, nout, time, &
       phi2, apar2, bpar2, igomega, &
       h, hk, omega, omegaavg, woutunits, tunits, phitot, ql_metric, &
       write_omega, write_heating, write_ql_metric)
    use gs2_heating, only: heating_diagnostics
# ifdef NETCDF
    use run_parameters, only: has_phi, has_apar, has_bpar, wstar_units
    use fields_arrays, only: phinew, aparnew, bparnew
    use parameter_scan_arrays, only: write_scan_parameter,current_scan_parameter_value
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    !> Current simulation time
    real, intent (in) :: time
    !> Volume average of fields squared
    real, intent(out) :: phi2, apar2, bpar2
    !> Theta-index to measure fields at
    integer, intent(in) :: igomega
    !> Normalisation factor for frequency
    real, dimension (:), intent (in) :: woutunits
    !> Normalisation factor for time
    real, dimension (:), intent (in) :: tunits
    !> Fields, frequency and time-averaged frequency
    complex, dimension(:,:), intent (in) :: omega, omegaavg
    !> FIXME: add documentation
    real, dimension(:,:), intent (in) :: phitot, ql_metric
    !> Heating diagnostics
    type(heating_diagnostics), intent (in) :: h
    !> Heating per-mode
    type(heating_diagnostics), dimension(:,:), intent (in) :: hk
    logical, intent(in) :: write_omega, write_heating, write_ql_metric
# ifdef NETCDF

    call neasyf_write(file_id, time_dim, time, dim_names=[time_dim], start=[nout])

    if (wstar_units) then
       call neasyf_write(file_id, "t_wstar", time * tunits, dim_names=[ky_dim, time_dim], &
            start=[1, nout], long_name="Time (wstar)", units="L/vt")
    end if

    if (has_phi) then
       call nc_write_field(file_id, nout, phinew, "phi", "Electrostatic potential", &
            igomega, write_phi_t, phi2)
    end if

    if (has_apar) then
       call nc_write_field(file_id, nout, aparnew, "apar", "Parallel magnetic potential", &
            igomega, write_apar_t, apar2)
    end if

    if (has_bpar) then
       call nc_write_field(file_id, nout, bparnew, "bpar", "Parallel magnetic field", &
            igomega, write_bpar_t, bpar2)
    end if

    if (write_heating) call nc_write_heating(file_id, nout, h, hk)

    if (write_omega) call nc_write_omega(file_id, nout, omega, omegaavg, woutunits)

    if (write_ql_metric) call nc_write_ql_metric(file_id, nout, ql_metric)

    if (write_scan_parameter) then
       call neasyf_write(file_id, &
            "current_scan_parameter_value", &
            current_scan_parameter_value, &
            dim_names=[time_dim], start=[nout], &
            long_name="The current value of the scan parameter")
     end if

    call neasyf_write(file_id, "phtot", phitot, dim_names=mode_dims, start=starts(3, nout), &
         long_name="Roughly the field line averaged summed field amplitudes")
# endif
  end subroutine nc_loop

  !> Writing the (normalised) complex frequency to netCDF
  subroutine nc_write_omega(file_id, nout, omega, omega_average, woutunits)
#ifdef NETCDF
    use kt_grids, only: ntheta0, naky
#endif
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    !> Normalisation factor for frequency
    real, dimension (:), intent (in) :: woutunits
    !> Fields, frequency and time-averaged frequency
    complex, dimension(:,:), intent (in) :: omega, omega_average
#ifdef NETCDF
    complex, dimension (ntheta0, naky) :: tmp
    integer :: it

    do it = 1, ntheta0
      tmp(it, :) = omega(it, :) * woutunits
    end do
    call netcdf_write_complex(file_id, "omega", tmp, dim_names=omega_dims, start=starts(4, nout), &
         long_name="Complex frequency", units="aref/vtref")

    do it = 1, ntheta0
      tmp(it, :) = omega_average(it, :) * woutunits
    end do

    ! FIXME: duplicate names for old/new diagnostics
    call netcdf_write_complex(file_id, "omegaavg", tmp, dim_names=omega_dims, start=starts(4, nout), &
         long_name="Complex frequency averaged over navg timesteps", units="aref/vtref")
    call netcdf_write_complex(file_id, "omega_average", tmp, dim_names=omega_dims, start=starts(4, nout), &
         long_name="Complex frequency averaged over navg timesteps", units="aref/vtref")
#endif
  end subroutine nc_write_omega

  !> Writes the QL flux metric by mode number and summed over mode number to netcdf file
  subroutine nc_write_ql_metric(file_id, nout, ql_metric)
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent (in) :: nout
    !> The QL flux metric as a function of kx and ky
    real, dimension(:,:), intent (in) :: ql_metric
#ifdef NETCDF
    character(*), dimension(*), parameter :: ql_dims = [kx_dim, ky_dim, time_dim]

    call neasyf_write(file_id, "ql_metric_by_mode", ql_metric, dim_names=ql_dims, &
         start=starts(3, nout), long_name="QL flux metric by mode number")
#endif
  end subroutine nc_write_ql_metric

  !> Write a field ((/E, A_\Vert, B_\vert/)) to netCDF
  subroutine nc_write_field(file_id, nout, field, field_name, long_name, igomega, write_field_over_time, field2)
# ifdef NETCDF
    use kt_grids, only: naky, ntheta0
    use neasyf, only: neasyf_write
    use volume_averages, only: average_theta, average_kx, average_ky
# endif
    use theta_grid, only: ntgrid
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent(in) :: nout
    !> One of the fields (/E, A_\Vert, B_\vert/)
    complex, dimension(-ntgrid:, :, :), intent(in) :: field
    !> Name of the field
    character(*), intent(in) :: field_name
    !> Description of the field
    character(*), intent(in) :: long_name
    !> \(\theta\)-index to measure field at
    integer, intent(in) :: igomega
    !> Should entire field be written every timestep?
    logical, intent(in) :: write_field_over_time
    !> Volume-average of the field amplitude squared
    real, intent(out) :: field2
# ifdef NETCDF
    real, dimension(ntheta0, naky) :: field2_by_mode
    real, dimension(ntheta0) :: field2_by_kx
    real, dimension(naky) :: field2_by_ky

    character(*), dimension(*), parameter :: field_dims = [complex_dim, theta_dim, kx_dim, ky_dim, time_dim]
    character(*), dimension(*), parameter :: kx_dims = [kx_dim, time_dim]
    character(*), dimension(*), parameter :: ky_dims = [ky_dim, time_dim]

    ! Note that we're not using `get_vol_average` because we'll end up
    ! with a circular dependency
    call average_theta(field, field, field2_by_mode, distributed=.false.)

    if (write_field_over_time) then
      call netcdf_write_complex(file_id, field_name // "_t", field, &
           dim_names=field_dims, start=starts(5, nout), &
           long_name=long_name // " over time")
    end if

    call average_ky(field2_by_mode, field2_by_kx, distributed=.false.)
    call neasyf_write(file_id, field_name // "2_by_kx", field2_by_kx, dim_names=kx_dims, &
         start=starts(2, nout), long_name=long_name // " squared and averaged over theta and ky")

    call average_kx(field2_by_mode, field2_by_ky, distributed=.false.)
    call neasyf_write(file_id, field_name // "2_by_ky", field2_by_ky, dim_names=ky_dims, &
         start=starts(2, nout), long_name=long_name // " squared and averaged over theta and kx")

    field2 = sum(field2_by_kx)

    ! FIXME: duplicated variable, old/new diagnostics
    call netcdf_write_complex(file_id, field_name // "0", field(igomega, :, :), &
         dim_names=omega_dims, start=starts(4, nout), &
         long_name=long_name // " at ig=igomega")
    call netcdf_write_complex(file_id, field_name // "_igomega_by_mode", field(igomega, :, :), &
         dim_names=omega_dims, start=starts(4, nout), &
         long_name=long_name // " at ig=igomega")

    call neasyf_write(file_id, field_name // "2_by_mode", field2_by_mode, &
         dim_names=mode_dims, start=starts(3, nout), &
         long_name=long_name // " per Fourier mode")
    call neasyf_write(file_id, field_name // "2", field2, dim_names=[time_dim], start=[nout], &
         long_name=long_name // " volume average", units="(T_ref/q_ref rho_ref/a_ref)**2")

    call nc_final_fields(file_id)
#endif
  end subroutine nc_write_field

  !> Writing heating diagnostics to netCDF
  subroutine nc_write_heating(file_id, nout, h, hk)
    use gs2_heating, only: heating_diagnostics
#ifdef NETCDF
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_heating, only: hk_repack
    use neasyf, only: neasyf_write
#endif
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent(in) :: nout
    !> Heating diagnostics
    type(heating_diagnostics), intent (in) :: h
    !> Heating per-mode
    type(heating_diagnostics), dimension(:,:), intent (in) :: hk
#ifdef NETCDF
    real, dimension (ntheta0, naky, nspec) :: tmps
    ! Note: Here parentheses around variables passed to neasyf_write are an intentional
    ! workaround to a gfortran bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105658
    call neasyf_write(file_id, "h_energy", (h%energy), dim_names=[time_dim], start=[nout], &
         long_name="Total free (turbulent) energy")
    call neasyf_write(file_id, "h_energy_dot", (h%energy_dot), dim_names=[time_dim], start=[nout], &
         long_name="Time derivative of total free energy")
    call neasyf_write(file_id, "h_antenna", (h%antenna), dim_names=[time_dim], start=[nout], &
         long_name="Energy injection from antenna")
    call neasyf_write(file_id, "h_eapar", (h%eapar), dim_names=[time_dim], start=[nout], &
         long_name="Free energy in perpendicular magnetic field")
    call neasyf_write(file_id, "h_ebpar", (h%ebpar), dim_names=[time_dim], start=[nout], &
         long_name="Free energy in parallel magnetic field")
    call neasyf_write(file_id, "h_delfs2", (h%delfs2), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Total free energy dfs^2")
    call neasyf_write(file_id, "h_hs2", (h%hs2), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Non-adiabatic free energy hs^2")
    call neasyf_write(file_id, "h_phis2", (h%phis2), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Adiabatic free energy phi_bar^2")
    call neasyf_write(file_id, "h_hypervisc", (h%hypervisc), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Heating from hyperviscosity")
    call neasyf_write(file_id, "h_hyperres", (h%hyperres), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Heating from hyperresitivity")
    call neasyf_write(file_id, "h_collisions", (h%collisions), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Heating from collisions -[h C(h) * T_0]")
    call neasyf_write(file_id, "h_gradients", (h%gradients), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Heating from gradient")
    call neasyf_write(file_id, "h_imp_colls", (h%imp_colls), dim_names=flux_dims, start=starts(2, nout), &
         long_name="[h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0")
    call neasyf_write(file_id, "h_hypercoll", (h%hypercoll), dim_names=flux_dims, start=starts(2, nout))
    call neasyf_write(file_id, "h_heating", (h%heating), dim_names=flux_dims, start=starts(2, nout), &
         long_name="Total heating by species")

    call neasyf_write(file_id, "hk_energy", (hk%energy), dim_names=mode_dims, start=starts(3, nout))
    call neasyf_write(file_id, "hk_energy_dot", (hk%energy_dot), dim_names=mode_dims, start=starts(3, nout))
    call neasyf_write(file_id, "hk_antenna", (hk%antenna), dim_names=mode_dims, start=starts(3, nout))
    call neasyf_write(file_id, "hk_eapar", (hk%eapar), dim_names=mode_dims, start=starts(3, nout))
    call neasyf_write(file_id, "hk_ebpar", (hk%ebpar), dim_names=mode_dims, start=starts(3, nout))

    call hk_repack (hk, 1, tmps)
    call neasyf_write(file_id, "hk_hypervisc", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 2, tmps)
    call neasyf_write(file_id, "hk_hyperres", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 3, tmps)
    call neasyf_write(file_id, "hk_collisions", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 4, tmps)
    call neasyf_write(file_id, "hk_gradients", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 5, tmps)
    call neasyf_write(file_id, "hk_heating", tmps, dim_names=fluxk_dims, start=starts(4, nout), &
         long_name="Total heating by species and mode")
    call hk_repack (hk, 6, tmps)
    call neasyf_write(file_id, "hk_hypercoll", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 7, tmps)
    call neasyf_write(file_id, "hk_delfs2", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 8, tmps)
    call neasyf_write(file_id, "hk_hs2", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 9, tmps)
    call neasyf_write(file_id, "hk_phis2", tmps, dim_names=fluxk_dims, start=starts(4, nout))
    call hk_repack (hk, 10, tmps)
    call neasyf_write(file_id, "hk_imp_colls", tmps, dim_names=fluxk_dims, start=starts(4, nout))
#endif
  end subroutine nc_write_heating

  !> Outputs full not guiding center moments
  subroutine nc_loop_fullmom (file_id, nout, ntot0, density0, upar0, tpar0, tperp0)
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current timestep
    integer, intent (in) :: nout
    complex, dimension(:,:,:), intent(in) :: ntot0, density0, upar0, tpar0, tperp0
#ifdef NETCDF
    integer, dimension (5) :: startmom
    character(*), dimension(*), parameter :: mom_dim = [complex_dim, kx_dim, ky_dim, species_dim, time_dim]

    startmom = [1, 1, 1, 1, nout]

    call netcdf_write_complex(file_id, "ntot0", ntot0, mom_dim, start=startmom, &
         long_name="Total perturbed density (theta=0) over time")
    call netcdf_write_complex(file_id, "density0", density0, mom_dim, start=startmom, &
         long_name="Nonadiabatic piece of perturbed density (theta=0) over time")
    call netcdf_write_complex(file_id, "upar0", upar0, mom_dim, start=startmom, &
         long_name="Upar (theta=0) over time")
    call netcdf_write_complex(file_id, "tpar0", tpar0, mom_dim, start=startmom, &
         long_name="Tpar (theta=0) over time")
    call netcdf_write_complex(file_id, "tperp0", tperp0, mom_dim, start=startmom, &
         long_name="Tperp (theta=0) over time")
#endif
  end subroutine nc_loop_fullmom

  !> Write \(\phi, A_\parallel, B_\parallel\) as functions of real space
  subroutine nc_loop_movie(file_id, nout_movie, time)
# ifdef NETCDF
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_layouts, only: yxf_lo
    use gs2_transforms, only: transform2
    use run_parameters, only: has_phi, has_apar, has_bpar
    use theta_grid, only: ntgrid
    use neasyf, only: neasyf_write
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current movie timestep
    integer, intent(in) :: nout_movie
    !> Current simulation time
    real, intent (in) :: time
# ifdef NETCDF
    ! \(\phi, A_\parallel, B_\parallel\) as functions of real space
    real, dimension (:,:,:), allocatable:: yxphi, yxapar, yxbpar
    character(*), dimension(*), parameter :: xmode_dim = [x_dim, y_dim, theta_dim, time_movie_dim]

    call neasyf_write(file_id, trim(time_movie_dim), time, dim_names=[time_movie_dim], start=[nout_movie])

    if (has_phi) then
       allocate(yxphi(yxf_lo%nx, yxf_lo%ny, -ntgrid:ntgrid))
       call transform2 (phinew, yxphi, yxf_lo%ny, yxf_lo%nx)
       call neasyf_write(file_id, "phi_by_xmode", yxphi, &
            dim_names=xmode_dim, start=starts(4, nout_movie))
       deallocate(yxphi)
    end if

    if (has_apar) then
       allocate(yxapar(yxf_lo%nx, yxf_lo%ny, -ntgrid:ntgrid))
       call transform2 (aparnew, yxapar, yxf_lo%ny, yxf_lo%nx)
       call neasyf_write(file_id, "apar_by_xmode", yxapar, &
            dim_names=xmode_dim, start=starts(4, nout_movie))
       deallocate(yxapar)
    end if

    if (has_bpar) then
       allocate(yxbpar(yxf_lo%nx, yxf_lo%ny, -ntgrid:ntgrid))
       call transform2 (bparnew, yxbpar, yxf_lo%ny, yxf_lo%nx)
       call neasyf_write(file_id, "bpar_by_xmode", yxbpar, &
            dim_names=xmode_dim, start=starts(4, nout_movie))
       deallocate(yxbpar)
    end if
# endif
  end subroutine nc_loop_movie

  !> Sync the files every nsync writes
  subroutine nc_sync(file_id, nout, movie_file_id, nout_movie, nsync_freq)
# ifdef NETCDF
    use netcdf, only: nf90_sync, NF90_NOERR
    use netcdf_utils, only: netcdf_error
# endif
    implicit none
    !> NetCDF ID of the file to write to
    integer, intent(in) :: file_id
    !> Current output timestep
    integer, intent(in) :: nout
    !> NetCDF ID of the movie file to write to
    integer, intent(in) :: movie_file_id
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
       status = nf90_sync (file_id)
       if (status /= NF90_NOERR) call netcdf_error (status)
    endif

    !Movie file
    if (movie_file_id >= 0) then
       if (mod(nout_movie, nsync) == 0) then
          status = nf90_sync(movie_file_id)
          if (status /= NF90_NOERR) call netcdf_error(status)
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
    ! Note: Here parentheses around variables passed to neasyf_write are an intentional
    ! workaround to a gfortran bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105658
    call neasyf_write(file_id, "charge", (spec%z), dim_names=["species"], &
         long_name="Charge", units="e")
    call neasyf_write(file_id, "mass", (spec%mass), dim_names=["species"], &
         long_name="Atomic mass", units="AMU")
    call neasyf_write(file_id, "dens", (spec%dens), dim_names=["species"], &
         long_name="Normalised density", units="nref")
    call neasyf_write(file_id, "temp", (spec%temp), dim_names=["species"], &
         long_name="Normalised temperature", units="Tref")
    call neasyf_write(file_id, "tprim", (spec%tprim), dim_names=["species"], &
         long_name="Normalised temperature gradient scale length -1/rho dT/drho", units="1/aref")
    call neasyf_write(file_id, "fprim", (spec%fprim), dim_names=["species"], &
         long_name="Normalised density gradient scale length -Lref/n dn/drho", units="1/aref")
    call neasyf_write(file_id, "uprim", (spec%uprim), dim_names=["species"], &
         long_name="Normalised parallel velocity gradient scale length -1/v_t du_par/drho", units="1/aref")
    call neasyf_write(file_id, "uprim2", (spec%uprim2), dim_names=["species"])
    call neasyf_write(file_id, "vnewk", (spec%vnewk), dim_names=["species"], &
         long_name="Collisionality", units="vtref/aref")

    ! FIXME: duplicated variable, old/new diagnostics
    call neasyf_write(file_id, "type_of_species", (spec%type), dim_names=["species"], &
         long_name="Species type: 1=ion, 2=electron, 3=trace, 4=hybrid_electron")
    call neasyf_write(file_id, "type", (spec%type), dim_names=["species"], &
         long_name="Species type: 1=ion, 2=electron, 3=trace, 4=hybrid_electron")
#endif
  end subroutine nc_species
  
  !> Write various geometric quantities to output netCDF file
  subroutine nc_geo(file_id)
#ifdef NETCDF
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob, &
         shat, drhodpsi, eps_trapped, cdrift, cdrift0, qval, rplot, &
         zplot, aplot, rprime, zprime, aprime, kxfac, bpol !, surfarea
    use kt_grids, only: kperp2
    use dist_fn, only: supercell_labels
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
    call neasyf_write(file_id, "kxfac", kxfac, long_name="kxfac factor used in NL term")
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
         long_name="Reference beta", units="2.mu0.nref.Tref/B_a**2")
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

    ! Wavenumber related quantities -- not strictly geometrical information.
    call neasyf_write(file_id, "kperp2", kperp2, dim_names=[theta_dim, kx_dim, ky_dim], &
         long_name="The perpendicular wavenumber squared", units="1/rho_r^2")
    call neasyf_write(file_id, "supercell_labels", supercell_labels, dim_names=[kx_dim, ky_dim], &
         long_name="Labels which unique supercell each wavenumber pair belongs to", units="N/A")

    ! FIXME: check if we need surfarea?
#endif
  end subroutine nc_geo

  !> Wrapper around `nf90_inquire_dimension`. Aborts if any errors are
  !> detected.
  function get_dim_length(file_id, name) result(length)
#ifdef NETCDF
    use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension, NF90_MAX_NAME
    use neasyf, only: neasyf_error
#endif
    integer :: length
    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Name of the dimension
    character(len=*), intent(in) :: name

#ifdef NETCDF
    integer :: dim_id
    ! Temporary variable for the dimension name
    character(len=NF90_MAX_NAME) :: dim_name
    ! NetCDF error status
    integer :: status

    status = nf90_inq_dimid(file_id, trim(name), dim_id)
    call neasyf_error(status, file_id)

    status = nf90_inquire_dimension(file_id, dim_id, dim_name, length)
    call neasyf_error(status, file_id, dimid=dim_id, dim=name)
#endif
  end function get_dim_length

  subroutine write_complex_scalar(parent_id, name, values, dim_names, units, long_name, start)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use convert, only: c2r
#endif
    !> Name of the variable
    character(len=*), intent(in) :: name
    !> NetCDF ID of the parent group/file
    integer, intent(in) :: parent_id
    !> Value of the integer to write
    complex, intent(in) :: values
    !> Array of dimension names, the first of which should be `"ri"`
    character(len=*), dimension(:), intent(in) :: dim_names
    !> Units of coordinate
    character(len=*), optional, intent(in) :: units
    !> Long descriptive name
    character(len=*), optional, intent(in) :: long_name
    !> Offset (one-indexed) for writing time-varying arrays. Must include the complex dimension as the first index
    integer, dimension(:), optional, intent(in) :: start

#ifdef NETCDF
    call neasyf_write(parent_id, name, [real(values), aimag(values)], &
         dim_names=dim_names, units=units, long_name=long_name, start=start)
#endif
  end subroutine write_complex_scalar

  subroutine write_complex_rank1(parent_id, name, values, dim_names, units, long_name, start)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use convert, only: c2r
#endif
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

#ifdef NETCDF
    real, dimension(2, size(values, 1)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
#endif
  end subroutine write_complex_rank1

  subroutine write_complex_rank2(parent_id, name, values, dim_names, units, long_name, start)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use convert, only: c2r
#endif
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

#ifdef NETCDF
    real, dimension(2, size(values, 1), size(values, 2)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
#endif
  end subroutine write_complex_rank2

  subroutine write_complex_rank3(parent_id, name, values, dim_names, units, long_name, start)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use convert, only: c2r
#endif
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

#ifdef NETCDF
    real, dimension(2, size(values, 1), size(values, 2), size(values, 3)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
#endif
  end subroutine write_complex_rank3

  subroutine write_complex_rank4(parent_id, name, values, dim_names, units, long_name, start)
#ifdef NETCDF
    use neasyf, only: neasyf_write
    use convert, only: c2r
#endif
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

#ifdef NETCDF
    real, dimension(2, size(values, 1), size(values, 2), size(values, 3), size(values, 4)) :: real_values

    call c2r(values, real_values)
    call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
#endif
  end subroutine write_complex_rank4
end module gs2_io
