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
  public :: nc_loop_partsym_tormom, nc_sync

# ifdef NETCDF
  integer (kind_nf) :: ncid
  integer (kind_nf) :: naky_dim, nakx_dim, nttot_dim, negrid_dim, nlambda_dim, nspec_dim, ncoord_dim, ncoordt_dim
  integer (kind_nf) :: nsign_dim, time_dim, char5_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim, vpa_dim
  integer (kind_nf) :: nttotext_dim, time_big_dim
  integer (kind_nf) :: aref_id, zref_id, nref_id, tref_id, mref_id, vref_id, bref_id, rhoref_id

  integer, dimension (6) :: lp_dim, mom_t_dim
  integer, dimension (5) :: field_dim, final_mom_dim, heatk_dim, err_dim
  integer, dimension (5) :: phi_corr_dim
  integer, dimension (4) :: omega_dim, fluxk_dim, final_field_dim, loop_mom_dim, sym_dim
  integer, dimension (4) :: phi2_extend_dim, phi_corr_2pi_dim
  integer, dimension (3) :: fluxx_dim
  integer, dimension (3) :: mode_dim, phase_dim, loop_phi_dim, heat_dim, flux_dist_dim, phi_dist_dim
  integer, dimension (2) :: kx_dim, ky_dim, om_dim, flux_dim, nin_dim, fmode_dim

  ! added by EAB 03/05/04 for movies
  logical :: my_make_movie, io_write_corr_extend
  integer :: ncid_movie
  integer :: nx_dim, ny_dim, nth_dim, time_movie_dim
  integer, dimension (4) :: xmode_dim
  integer :: nx_id, ny_id, nth_id, x_id, y_id, th_id, time_movie_id
  integer :: phi_by_xmode_id, apar_by_xmode_id, bpar_by_xmode_id

  integer :: nakx_id, naky_id, nttot_id, akx_id, aky_id, theta_id, nspec_id
  integer :: negrid_id, nlambda_id, egrid_id, lambda_id
  integer :: time_id, phi2_id, apar2_id, bpar2_id, theta0_id, nproc_id, nmesh_id
  integer :: current_scan_parameter_value_id
  integer :: phi2_by_mode_id, apar2_by_mode_id, bpar2_by_mode_id, errest_by_mode_id, lpcoef_by_mode_id
  integer :: phtot_id
  integer :: phi2_by_kx_id, apar2_by_kx_id, bpar2_by_kx_id
  integer :: phi2_by_ky_id, apar2_by_ky_id, bpar2_by_ky_id
  integer :: phi0_id, apar0_id, bpar0_id
  integer :: omega_id, omegaavg_id, phase_id
  integer :: es_heat_flux_id, es_mom_flux_id, es_part_flux_id, es_exchange_id
  integer :: es_part_tormom_flux_id
  integer :: es_part_flux_dist_id, es_mom_flux_par_dist_id, es_mom_flux_perp_dist_id, es_heat_flux_dist_id, es_phi_dist_id
  integer :: es_heat_par_id, es_heat_perp_id
  integer :: apar_heat_flux_id, apar_mom_flux_id, apar_part_flux_id
  integer :: apar_heat_par_id, apar_heat_perp_id
  integer :: bpar_heat_flux_id, bpar_mom_flux_id, bpar_part_flux_id
  integer :: bpar_heat_par_id, bpar_heat_perp_id
  integer :: hflux_tot_id, zflux_tot_id, vflux_tot_id
  integer :: es_heat_by_k_id, es_mom_by_k_id, es_part_by_k_id
  integer :: es_part_tormom_by_k_id, es_exchange_by_k_id
  integer :: es_parmom_by_k_id, es_perpmom_by_k_id, es_mom0_by_k_id, es_mom1_by_k_id
  integer :: es_mom_sym_id, phi_corr_id, phi2_extend_id, phi_corr_2pi_id
  integer :: es_part_sym_id 
  integer :: apar_heat_by_k_id, apar_mom_by_k_id, apar_part_by_k_id
  integer :: apar_heat_by_x_id
  integer :: bpar_heat_by_k_id, bpar_mom_by_k_id, bpar_part_by_k_id
  integer :: phi_t_id, apar_t_id, bpar_t_id
  integer :: ntot_t_id, density_t_id, upar_t_id, tpar_t_id, tperp_t_id
  integer :: qparflux_t_id, pperpj1_t_id, qpperpj1_t_id
  integer :: ntot0_id, density0_id, upar0_id, tpar0_id, tperp0_id
  integer :: qparflux0_id, pperpj10_id, qpperpj10_id
  integer :: phi_norm_id, apar_norm_id, bpar_norm_id
  integer :: phi_id, apar_id, bpar_id, epar_id
  integer :: antot_id, antota_id, antotp_id
  integer :: ntot_id, density_id, upar_id, tpar_id, tperp_id
  integer :: qparflux_id, pperpj1_id, qpperpj1_id
  integer :: ntot2_id, ntot2_by_mode_id, ntot20_id, ntot20_by_mode_id
  integer :: tpar2_by_mode_id, tperp2_by_mode_id
  integer :: phi00_id, ntot00_id, density00_id, upar00_id, tpar00_id, tperp00_id
  integer :: input_id
  integer :: charge_id, mass_id, dens_id, temp_id, tprim_id, fprim_id
  integer :: uprim_id, uprim2_id, vnewk_id, spec_type_id
  integer :: bmag_id, gradpar_id, gbdrift_id, gbdrift0_id
  integer :: cdrift_id, cdrift0_id
  integer :: cvdrift_id, cvdrift0_id, gds2_id, gds21_id, gds22_id, gds23_id, gds24_id
  integer :: gds24_noq_id, grho_id, jacob_id, shat_id, eps_id, drhodpsi_id, q_id, surfarea_id
  integer :: rplot_id, zplot_id, aplot_id, rprime_id, zprime_id, aprime_id, bpol_id
  integer :: beta_id
  integer :: code_id, layout_id
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
  integer, dimension (5) :: mom_dim
  logical :: write_apar_t, write_phi_t, write_bpar_t ! Should the fields be written out every nwrite?

# endif
  real :: zero
  
  interface ensure_netcdf_var_exists
    module procedure ensure_netcdf_var_exists_scalar
    module procedure ensure_netcdf_var_exists_onedim
    module procedure ensure_netcdf_var_exists_manydims
  end interface ensure_netcdf_var_exists

  interface ensure_netcdf_att_exists
    module procedure ensure_netcdf_att_exists_text
  end interface ensure_netcdf_att_exists

!  include 'netcdf.inc'
  
contains
  
  !> Create and initialise the output netCDF files
  !>
  !> Create the main and movie (real-space fields) netCDF files,
  !> creating the dimensions and main variables
  subroutine init_gs2_io (append_old, write_nl_flux, write_omega, &
      write_hrate, write_final_antot, write_eigenfunc, &
      make_movie, nmovie_tot, write_verr, &
      write_moments,&
      write_full_moments_notgc, write_sym, write_pflux_sym, write_nl_flux_dist, &
      write_pflux_tormom, write_correlation, nwrite_big_tot, write_correlation_extend, & 
      write_phi_over_time, write_apar_over_time, write_bpar_over_time, &
      write_avg_moments, write_final_moments, write_final_epar, nout, &
      ob_midplane) 
    use mp, only: proc0
    use init_g, only: tstart
# ifdef NETCDF
    use mp, only: barrier
    use file_utils, only: run_name
    use gs2_transforms, only: init_transforms
    use kt_grids, only: naky, ntheta0,nx,ny
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use netcdf, only: NF90_CLOBBER, nf90_create, nf90_open, NF90_WRITE
    use netcdf_utils, only: get_netcdf_code_precision, netcdf_real
    use constants, only: run_name_size
# endif
    implicit none
    !> See [[gs2_diagnostics_knobs]] for details of these inputs
    logical, intent(in) :: append_old
    logical, intent(in) :: write_nl_flux, write_nl_flux_dist, write_omega
    logical, intent(in) :: write_hrate, make_movie, write_moments
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
# ifdef NETCDF
    logical :: accelerated
    character(run_name_size) :: filename, filename_movie
    integer :: status
    logical :: previous_file_exists, appending

    !################
    !# Update internal flags

    write_phi_t = write_phi_over_time
    write_apar_t = write_apar_over_time
    write_bpar_t = write_bpar_over_time

    my_make_movie = make_movie
    io_write_corr_extend = write_correlation_extend

    !----------------

    zero = epsilon(0.0)

    if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

    status = NF90_NOERR

    !Create file
    filename = trim(trim(run_name)//'.out.nc')

    ! only proc0 opens the file:
    if (proc0) then
       inquire(file=trim(filename), exist=previous_file_exists)
       appending = (previous_file_exists .and. append_old)
       if (appending) then
          status = nf90_open(trim(filename), NF90_WRITE, ncid)
       else
          status = nf90_create (trim(filename), NF90_CLOBBER, ncid)
       end if
       if (status /= NF90_NOERR) call netcdf_error (status, file=filename, abort=.true.)
    end if

    !Create the movie file
    if (my_make_movie) then
       !perhaps put a call to init_y_transform
       !    call init_y_transform_layouts  
       call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
       filename_movie = trim(trim(run_name)//'.movie.nc')
       ! only proc0 opens the file:
       if (proc0) then
          status = nf90_create (trim(filename_movie), NF90_CLOBBER, ncid_movie)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, file=filename_movie)
       end if
    endif

# endif

    if (proc0) then
        !Setup main file by creating netCDF file
        call define_dims (nwrite_big_tot, io_write_corr_extend)
        call define_vars (write_nl_flux, write_omega, &
                          write_hrate, write_final_antot, write_eigenfunc, write_verr, &
                          write_moments, write_full_moments_notgc, write_sym, write_pflux_sym, &
                          write_nl_flux_dist, write_pflux_tormom, write_correlation, &
                          write_correlation_extend, write_avg_moments, write_final_moments, write_final_epar, &
                          ob_midplane=ob_midplane)

        if (appending) then
            nout = get_netcdf_dim_length(ncid, time_dim)
        else
            !Write constant data
            call nc_grids
            call nc_norms
            call nc_species
            call nc_geo
            call save_input
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

  end subroutine init_gs2_io

  !> Get the netCDF ID for a dimension, creating it if it doesn't
  !> exist already
  subroutine ensure_netcdf_dim_exists(file_id, dim_name, dim_size, dim_id)
#ifdef NETCDF
    use netcdf, only: nf90_inq_dimid, nf90_def_dim, NF90_NOERR, NF90_EBADDIM
#endif
    !> ID of the file or group to look for or create the dimension under
    integer(kind_nf), intent(in) :: file_id
    !> Name of the dimension
    character(len=*), intent(in) :: dim_name
    !> The size of the dimension
    integer(kind_nf), intent(in) :: dim_size
    !> The netCDF ID of the dimension under `file_id`
    integer(kind_nf), intent(out) :: dim_id

    ! Error code of the netCDF calls
    integer(kind_nf) :: status

#ifdef NETCDF
    status = nf90_inq_dimid(file_id, dim_name, dim_id)
    ! Dimension doesn't exist, so let's create it
    if (status == NF90_EBADDIM) then
      status = nf90_def_dim(file_id, dim_name, dim_size, dim_id)
    end if
    ! Something went wrong with one of the previous two calls
    if (status /= NF90_NOERR) then
      call netcdf_error(status, dim=dim_name, message="(ensure_netcdf_dim_exists)", abort=.true.)
    end if
#endif
  end subroutine ensure_netcdf_dim_exists

  !> Ensure the netCDF file contains all the dimensions, creating them
  !> if necessary. Doesn't fill the values.
  !>
  !> Dimension variable IDs are module level variables
  subroutine define_dims (nwrite_big_tot, create_extended_correlation)
# ifdef NETCDF
    use file_utils, only: num_input_lines
    use gs2_time, only: user_time
    use kt_grids, only: naky, ntheta0, jtwist_out
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use netcdf, only: NF90_UNLIMITED
# endif
    implicit none
    !> Total number of output timesteps for "large" diagnostics
    integer, intent(in) :: nwrite_big_tot
    !> Create the dimensions for the extended correlation diagnostics
    logical, intent(in) :: create_extended_correlation
# ifdef NETCDF
    integer(kind_nf) :: status

    call ensure_netcdf_dim_exists(ncid, 'ky', naky, naky_dim)
    call ensure_netcdf_dim_exists(ncid, 'kx', ntheta0, nakx_dim)
    call ensure_netcdf_dim_exists(ncid, 'theta', 2*ntgrid+1, nttot_dim)
    call ensure_netcdf_dim_exists(ncid, 'lambda', nlambda, nlambda_dim)
    call ensure_netcdf_dim_exists(ncid, 'egrid', negrid, negrid_dim)
    call ensure_netcdf_dim_exists(ncid, 'species', nspec, nspec_dim)
    call ensure_netcdf_dim_exists(ncid, 'sign', 2, nsign_dim)
    call ensure_netcdf_dim_exists(ncid, 't', NF90_UNLIMITED, time_dim)
    call ensure_netcdf_dim_exists(ncid, 'char5', 5, char5_dim)
    call ensure_netcdf_dim_exists(ncid, 'char10', 10, char10_dim)
    call ensure_netcdf_dim_exists(ncid, 'char200', 200, char200_dim)
    call ensure_netcdf_dim_exists(ncid, 'nlines', num_input_lines, nlines_dim)
    call ensure_netcdf_dim_exists(ncid, 'ri', 2, ri_dim)
    call ensure_netcdf_dim_exists(ncid, 'nheat', 7, nheat_dim)
    call ensure_netcdf_dim_exists(ncid, 'vpa', negrid*nlambda, vpa_dim)
    if (create_extended_correlation) then
       call ensure_netcdf_dim_exists(ncid, 'theta_ext', (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1), nttotext_dim)
       call ensure_netcdf_dim_exists(ncid, 'tb', nwrite_big_tot, time_big_dim)
    end if
# endif
  end subroutine define_dims

  !> Wrapper around `nf90_inquire_dimension`, includes error handling
  function get_netcdf_dim_length(file_id, dim_id) result(length)
#ifdef NETCDF
    use netcdf, only: nf90_inquire_dimension, NF90_MAX_NAME, NF90_NOERR
#endif
    integer :: length
    !> NetCDF ID of the file
    integer(kind_nf), intent(in) :: file_id
    !> NetCDF ID of the dimension
    integer(kind_nf), intent(in) :: dim_id

    ! Temporary variable for the dimension name
    character(len=NF90_MAX_NAME) :: dim_name
    ! NetCDF error status
    integer :: status

#ifdef NETCDF
    ! This could be less reliant on module-level variables if we
    ! passed in the dimension name and used nf90_inq_dim to get the ID
    status = nf90_inquire_dimension(file_id, dim_id, dim_name, length)
    if (status /= NF90_NOERR) then
      call netcdf_error(status, &
                        dim=trim(dim_name), &
                        message="(get_netcdf_dim_length)", &
                        abort=.true.)
    end if
#endif
  end function get_netcdf_dim_length

  !> Define the dimensions used in the movie file. Some duplication with
  !! main file -- may want to merge movie data into main file
  !! Main reason for using a new file is we can only have one unlimited dimension
  !! per file (with new versions of NETCDF this is no longer true) and we
  !! probably don't want to write movie data as often as other data.
  !! However, at the moment we precalculate the number of outputs, so it
  !! isn't actually unlimited
  subroutine define_dims_mymovie(nmovie_tot)
# ifdef NETCDF
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use netcdf, only: nf90_def_dim
# endif
    implicit none
    integer, intent(in) :: nmovie_tot
# ifdef NETCDF
    integer :: nx, ny, status

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

  !> Write grid dimensions to the output netCDF file
  !>
  !> This must be called after [[define_dims]]
  subroutine nc_grids
# ifdef NETCDF
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, nx, ny, akx, aky
    use species, only: nspec
    use le_grids, only: negrid, nlambda, energy, al
    use nonlinear_terms, only: nonlin
    use netcdf, only: nf90_put_var
    use gs2_layouts, only: layout
# endif    
    implicit none
# ifdef NETCDF
    integer :: status
    real :: nmesh

    status = nf90_put_var (ncid, layout_id, layout)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, layout_id)

    status = nf90_put_var (ncid, nttot_id, 2*ntgrid+1)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nttot_id)
    status = nf90_put_var (ncid, naky_id, naky)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, naky_id)
    status = nf90_put_var (ncid, nakx_id, ntheta0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nakx_id)
    status = nf90_put_var (ncid, nspec_id, nspec)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nspec_id)
    status = nf90_put_var (ncid, negrid_id, negrid)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, negrid_id)
    status = nf90_put_var (ncid, nlambda_id, nlambda)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nlambda_id)

    status = nf90_put_var (ncid, akx_id, akx)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, akx_id)
    status = nf90_put_var (ncid, aky_id, aky)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aky_id)
    status = nf90_put_var (ncid, theta_id, theta)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, theta_id)
    status = nf90_put_var (ncid, theta0_id, theta0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, theta0_id)
    status = nf90_put_var (ncid, egrid_id, energy)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, egrid_id)
    status = nf90_put_var (ncid, lambda_id, al)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, lambda_id)

    if (nonlin) then
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
    else
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
    end if

    status = nf90_put_var (ncid, nmesh_id, nmesh)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nmesh_id)
# endif
  end subroutine nc_grids

  !> Write the physical normalisations to the output netCDF file
  subroutine nc_norms
# ifdef NETCDF
    use netcdf, only: nf90_put_var
    use normalisations, only: norms
# endif    
    implicit none
# ifdef NETCDF
    integer :: status

    status = nf90_put_var (ncid, aref_id, norms%get_value("aref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aref_id)
    status = nf90_put_var (ncid, zref_id, norms%get_value("zref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, zref_id)
    status = nf90_put_var (ncid, nref_id, norms%get_value("nref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nref_id)
    status = nf90_put_var (ncid, tref_id, norms%get_value("tref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tref_id)
    status = nf90_put_var (ncid, mref_id, norms%get_value("mref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, mref_id)
    status = nf90_put_var (ncid, vref_id, norms%get_value("vref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vref_id)
    status = nf90_put_var (ncid, bref_id, norms%get_value("bref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bref_id)
    status = nf90_put_var (ncid, rhoref_id, norms%get_value("rhoref"))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, rhoref_id)

# endif
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
  subroutine save_input
# ifdef NETCDF    
    use file_utils, only: num_input_lines, get_input_unit
    use netcdf, only: nf90_put_var
    implicit none
    character(200) line
    integer, dimension (2) :: nin_start, nin_count
    integer :: status, n, unit

    nin_start(1) = 1
    nin_start(2) = 1

    nin_count(2) = 1

    call get_input_unit (unit)
    rewind (unit=unit)
    do n = 1, num_input_lines
       read (unit=unit, fmt="(a)") line
       nin_count(1) = len(trim(line))
       status = nf90_put_var (ncid, input_id, line, start=nin_start, count=nin_count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, input_id)
       nin_start(2) = nin_start(2) + 1
    end do
# endif
  end subroutine save_input

  !> Get the netCDF ID for a variable, creating it if it doesn't exist
  !> already
  subroutine ensure_netcdf_var_exists_scalar(file_id, var_name, var_type, var_id)
#ifdef NETCDF
    use netcdf, only: nf90_inq_varid, nf90_def_var, NF90_NOERR, NF90_ENOTVAR
#endif
    !> ID of the file or group to look for or create the variable under
    integer(kind_nf), intent(in) :: file_id
    !> Name of the variable
    character(len=*), intent(in) :: var_name
    !> The netCDF type of the variable
    integer(kind_nf), intent(in) :: var_type
    !> The netCDF ID of the variable under `file_id`
    integer(kind_nf), intent(out) :: var_id

    ! Error code of the netCDF calls
    integer(kind_nf) :: status

#ifdef NETCDF
    status = nf90_inq_varid(file_id, var_name, var_id)
    ! Variable doesn't exist, so let's create it
    if (status == NF90_ENOTVAR) then
      status = nf90_def_var(file_id, var_name, var_type, var_id)
    end if
    ! Something went wrong with one of the previous two calls
    if (status /= NF90_NOERR) then
      call netcdf_error(status, var=var_name, varid=var_id, &
                        message="(ensure_netcdf_var_exists)", abort=.true.)
    end if
#endif
  end subroutine ensure_netcdf_var_exists_scalar

  !> Get the netCDF ID for a variable, creating it if it doesn't exist
  !> already
  subroutine ensure_netcdf_var_exists_onedim(file_id, var_name, var_type, dim_id, var_id)
    !> ID of the file or group to look for or create the variable under
    integer(kind_nf), intent(in) :: file_id
    !> Name of the variable
    character(len=*), intent(in) :: var_name
    !> The netCDF type of the variable
    integer(kind_nf), intent(in) :: var_type
    !> Array of dimension IDs
    integer(kind_nf), intent(in) :: dim_id
    !> The netCDF ID of the variable under `file_id`
    integer(kind_nf), intent(out) :: var_id

    call ensure_netcdf_var_exists_manydims(file_id, var_name, var_type, [dim_id], var_id)
  end subroutine ensure_netcdf_var_exists_onedim

  !> Get the netCDF ID for a variable, creating it if it doesn't exist
  !> already
  subroutine ensure_netcdf_var_exists_manydims(file_id, var_name, var_type, dim_id, var_id)
#ifdef NETCDF
    use netcdf, only: nf90_inq_varid, nf90_def_var, NF90_NOERR, NF90_ENOTVAR
#endif
    !> ID of the file or group to look for or create the variable under
    integer(kind_nf), intent(in) :: file_id
    !> Name of the variable
    character(len=*), intent(in) :: var_name
    !> The netCDF type of the variable
    integer(kind_nf), intent(in) :: var_type
    !> Array of dimension IDs
    integer(kind_nf), dimension(:), intent(in) :: dim_id
    !> The netCDF ID of the variable under `file_id`
    integer(kind_nf), intent(out) :: var_id

    ! Error code of the netCDF calls
    integer(kind_nf) :: status

#ifdef NETCDF
    status = nf90_inq_varid(file_id, var_name, var_id)
    ! Variable doesn't exist, so let's create it
    if (status == NF90_ENOTVAR) then
      status = nf90_def_var(file_id, var_name, var_type, dim_id, var_id)
    end if
    ! Something went wrong with one of the previous two calls
    if (status /= NF90_NOERR) then
      call netcdf_error(status, var=var_name, varid=var_id, &
                        message="(ensure_netcdf_var_exists)", abort=.true.)
    end if
#endif
  end subroutine ensure_netcdf_var_exists_manydims

  subroutine ensure_netcdf_att_exists_text(file_id, var_id, att_name, att_text)
#ifdef NETCDF
    use netcdf, only: nf90_inquire_attribute, nf90_put_att, NF90_NOERR, NF90_ENOTATT
#endif
    !> ID of the file or group to look for or create the variable under
    integer(kind_nf), intent(in) :: file_id
    !> The netCDF ID of the variable under `file_id`
    integer(kind_nf), intent(in) :: var_id
    !> Name of the attribute
    character(len=*), intent(in) :: att_name
    !> Attribute text to store
    character(len=*), intent(in) :: att_text

    ! Error code of the netCDF calls
    integer(kind_nf) :: status

#ifdef NETCDF
    status = nf90_inquire_attribute(file_id, var_id, att_name)
    ! Attribute doesn't exist, so let's create it
    if (status == NF90_ENOTATT) then
      status = nf90_put_att(file_id, var_id, att_name, att_text)
    end if
    ! Something went wrong with one of the previous two calls
    if (status /= NF90_NOERR) then
      call netcdf_error(status, file_id, var_id, att=att_name, &
                        message="(ensure_netcdf_att_exists)", abort=.true.)
    end if
#endif
  end subroutine ensure_netcdf_att_exists_text

  !> Create (but don't fill) the various diagnostics variables in the
  !> output netCDF file
  !>
  !> This also calls `nf90_enddef` which leaves "define" mode, so no
  !> more variables can be added to the file after this call
  subroutine define_vars (write_nl_flux, write_omega, &
       write_hrate, write_final_antot, write_eigenfunc, write_verr, &
       write_moments, write_full_moments_notgc, write_sym, write_pflux_sym, &
       write_nl_flux_dist, write_pflux_tormom, write_correlation, &
       write_correlation_extend, write_avg_moments, write_final_moments, write_final_epar, &
       ob_midplane)
# ifdef NETCDF
    use file_utils, only: error_unit
    use mp, only: nproc
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use parameter_scan_arrays, only: write_scan_parameter
    use netcdf, only: NF90_CHAR, NF90_INT, NF90_GLOBAL
    use netcdf, only: nf90_enddef, nf90_inq_libvers, nf90_put_var
    use netcdf_utils, only: netcdf_real
# endif
    implicit none
    logical, intent(in) :: write_nl_flux, write_nl_flux_dist, write_omega
    logical, intent(in) :: write_hrate, write_final_antot
    logical, intent(in) :: write_eigenfunc, write_verr, write_moments
    logical, intent(in) :: write_full_moments_notgc, write_sym, write_correlation
    logical, intent(in) :: write_pflux_sym, write_pflux_tormom
    logical, intent(inout) :: write_correlation_extend
    logical, intent(in) :: write_avg_moments, write_final_moments, write_final_epar
    logical, optional, intent(in) :: ob_midplane
# ifdef NETCDF
    character (20) :: datestamp, timestamp, timezone, tmpc
    integer :: status, ierr, dimsize(5)
    real    :: total
    logical :: obmid=.false.

    !Handle optional input
    if (present(ob_midplane)) then
       obmid=ob_midplane
    endif

    !Define dimensions/shapes for different data types
    fmode_dim(1) = nakx_dim
    fmode_dim(2) = naky_dim

    mode_dim (1) = nakx_dim
    mode_dim (2) = naky_dim
    mode_dim (3) = time_dim

    err_dim (1) = nttot_dim
    err_dim (2) = nakx_dim
    err_dim (3) = naky_dim
    err_dim (4) = ncoordt_dim
    err_dim (5) = time_dim

    lp_dim (1) = nttot_dim
    lp_dim (2) = nakx_dim
    lp_dim (3) = naky_dim
    lp_dim (4) = nspec_dim
    lp_dim (5) = ncoord_dim
    lp_dim (6) = time_dim

    kx_dim (1) = nakx_dim
    kx_dim (2) = time_dim

    ky_dim (1) = naky_dim
    ky_dim (2) = time_dim

    om_dim (1) = ri_dim
    om_dim (2) = time_dim

    omega_dim (1) = ri_dim
    omega_dim (2) = nakx_dim
    omega_dim (3) = naky_dim
    omega_dim (4) = time_dim

    phase_dim (1) = ri_dim
    phase_dim (2) = nakx_dim
    phase_dim (3) = naky_dim

    nin_dim(1) = char200_dim
    nin_dim(2) = nlines_dim

    flux_dim (1) = nspec_dim
    flux_dim (2) = time_dim

    fluxk_dim (1) = nakx_dim
    fluxk_dim (2) = naky_dim
    fluxk_dim (3) = nspec_dim
    fluxk_dim (4) = time_dim

    fluxx_dim (1) = nakx_dim
    fluxx_dim (2) = nspec_dim
    fluxx_dim (3) = time_dim

    heat_dim (1) = nspec_dim
    heat_dim (2) = nheat_dim
    heat_dim (3) = time_dim

    heatk_dim (1) = nakx_dim
    heatk_dim (2) = naky_dim
    heatk_dim (3) = nspec_dim
    heatk_dim (4) = nheat_dim
    heatk_dim (5) = time_dim

    field_dim (1) = ri_dim
    field_dim (2) = nttot_dim
    field_dim (3) = nakx_dim
    field_dim (4) = naky_dim
    field_dim (5) = time_dim

    mom_t_dim (1) = ri_dim
    mom_t_dim (2) = nttot_dim
    mom_t_dim (3) = nakx_dim
    mom_t_dim (4) = naky_dim
    mom_t_dim (5) = nspec_dim
    mom_t_dim (6) = time_dim

    final_field_dim (1) = ri_dim
    final_field_dim (2) = nttot_dim
    final_field_dim (3) = nakx_dim
    final_field_dim (4) = naky_dim

    final_mom_dim (1) = ri_dim
    final_mom_dim (2) = nttot_dim
    final_mom_dim (3) = nakx_dim
    final_mom_dim (4) = naky_dim
    final_mom_dim (5) = nspec_dim

    loop_mom_dim (1) = ri_dim
    loop_mom_dim (2) = nakx_dim
    loop_mom_dim (3) = nspec_dim
    loop_mom_dim (4) = time_dim

    loop_phi_dim (1) = ri_dim
    loop_phi_dim (2) = nakx_dim
    loop_phi_dim (3) = time_dim

    mom_dim(1) = ri_dim
    mom_dim(2) = nakx_dim
    mom_dim(3) = naky_dim
    mom_dim(4) = nspec_dim
    mom_dim(5) = time_dim

    flux_dist_dim(1) = nttot_dim
    flux_dist_dim(2) = nspec_dim
    flux_dist_dim(3) = time_dim

    phi_dist_dim(1) = ri_dim
    phi_dist_dim(2) = nttot_dim
    phi_dist_dim(3) = time_dim

    sym_dim(1) = nttot_dim
    sym_dim(2) = vpa_dim
    sym_dim(3) = nspec_dim
    sym_dim(4) = time_dim


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
       phi2_extend_dim(1) = nttotext_dim
       phi2_extend_dim(2) = nakx_dim
       phi2_extend_dim(3) = naky_dim
       phi2_extend_dim(4) = time_big_dim

       phi_corr_dim(1) = ri_dim
       phi_corr_dim(2) = nttotext_dim
       phi_corr_dim(3) = nakx_dim
       phi_corr_dim(4) = naky_dim
       phi_corr_dim(5) = time_big_dim
    end if

    phi_corr_2pi_dim(1) = ri_dim
    phi_corr_2pi_dim(2) = nttot_dim
    phi_corr_2pi_dim(3) = naky_dim
    phi_corr_2pi_dim(4) = time_dim

    !<doc> Write some useful general information such as the website,
    ! date and time into the NetCDF file </doc>
    call ensure_netcdf_att_exists(ncid, NF90_GLOBAL, "title", "GS2 Simulation Data")
    call ensure_netcdf_att_exists(ncid, NF90_GLOBAL, "Conventions", "https://gyrokinetics.gitlab.io/gs2/")

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

    !This stores a copy of the (sanitised?) input file
    call ensure_netcdf_var_exists(ncid, 'input_file', NF90_CHAR, nin_dim, input_id)
    call ensure_netcdf_att_exists(ncid, input_id, 'long_name', 'Input file')

    !################################################
    !# NOW HANDLE DATA ABOUT GRID SIZES
    ! Write lots of input variables (e.g. nproc, nkx, nky)
    ! into the NetCDF file

    call ensure_netcdf_var_exists(ncid, 'layout', NF90_CHAR, char5_dim, layout_id)
    call ensure_netcdf_var_exists(ncid, 'nproc', NF90_INT, nproc_id)
    call ensure_netcdf_att_exists(ncid, nproc_id, 'long_name', 'Number of processors')

    call ensure_netcdf_var_exists(ncid, 'nmesh', netcdf_real, nmesh_id)
    call ensure_netcdf_att_exists(ncid, nmesh_id, 'long_name', 'Number of meshpoints')

    call ensure_netcdf_var_exists(ncid, 'nkx', NF90_INT, nakx_id)
    call ensure_netcdf_var_exists(ncid, 'nky', NF90_INT, naky_id)
    call ensure_netcdf_var_exists(ncid, 'ntheta_tot', NF90_INT, nttot_id)
    call ensure_netcdf_var_exists(ncid, 'nspecies', NF90_INT, nspec_id)
    call ensure_netcdf_var_exists(ncid, 'nenergy', NF90_INT, negrid_id)
    call ensure_netcdf_var_exists(ncid, 'nlambda', NF90_INT, nlambda_id)

    !################################################
    !# NOW DO NORMALISATIONS
    ! Could use norms%check_if_set to decide what we want to write here

    call ensure_netcdf_var_exists(ncid, 'aref', netcdf_real, aref_id)
    call ensure_netcdf_var_exists(ncid, 'zref', netcdf_real, zref_id)
    call ensure_netcdf_var_exists(ncid, 'nref', netcdf_real, nref_id)
    call ensure_netcdf_var_exists(ncid, 'tref', netcdf_real, tref_id)
    call ensure_netcdf_var_exists(ncid, 'mref', netcdf_real, mref_id)
    call ensure_netcdf_var_exists(ncid, 'vref', netcdf_real, vref_id)
    call ensure_netcdf_var_exists(ncid, 'bref', netcdf_real, bref_id)
    call ensure_netcdf_var_exists(ncid, 'rhoref', netcdf_real, rhoref_id)

    !################################################
    !# NOW DO SPECIES DATA

    call ensure_netcdf_var_exists(ncid, 'charge', NF90_INT, nspec_dim, charge_id)
    call ensure_netcdf_att_exists(ncid, charge_id, 'long_name', 'Charge')
    call ensure_netcdf_att_exists(ncid, charge_id, 'units', 'q')

    call ensure_netcdf_var_exists(ncid, 'mass', netcdf_real, nspec_dim, mass_id)
    call ensure_netcdf_att_exists(ncid, mass_id, 'long_name', 'Atomic mass')
    call ensure_netcdf_att_exists(ncid, mass_id, 'units', 'm')

    call ensure_netcdf_var_exists(ncid, 'dens', netcdf_real, nspec_dim, dens_id)
    call ensure_netcdf_att_exists(ncid, dens_id, 'long_name', 'Density')
    call ensure_netcdf_att_exists(ncid, dens_id, 'units', 'n_e')

    call ensure_netcdf_var_exists(ncid, 'temp', netcdf_real, nspec_dim, temp_id)
    call ensure_netcdf_att_exists(ncid, temp_id, 'long_name', 'Temperature')
    call ensure_netcdf_att_exists(ncid, temp_id, 'units', 'T')

    call ensure_netcdf_var_exists(ncid, 'tprim', netcdf_real, nspec_dim, tprim_id)
    call ensure_netcdf_att_exists(ncid, tprim_id, 'long_name', '-1/rho dT/drho')

    call ensure_netcdf_var_exists(ncid, 'fprim', netcdf_real, nspec_dim, fprim_id)
    call ensure_netcdf_att_exists(ncid, fprim_id, 'long_name', '-1/rho dn/drho')

    call ensure_netcdf_var_exists(ncid, 'uprim', netcdf_real, nspec_dim, uprim_id)
    call ensure_netcdf_att_exists(ncid, uprim_id, 'long_name', '-1/v_t du_par/drho')

    call ensure_netcdf_var_exists(ncid, 'uprim2', netcdf_real, nspec_dim, uprim2_id)

    call ensure_netcdf_var_exists(ncid, 'vnewk', netcdf_real, nspec_dim, vnewk_id)
    call ensure_netcdf_att_exists(ncid, vnewk_id, 'long_name', 'Collisionality')
    call ensure_netcdf_att_exists(ncid, vnewk_id, 'units', 'v_t/L')

    call ensure_netcdf_var_exists(ncid, 'type_of_species', NF90_INT, nspec_dim, spec_type_id)
    !------------------------------------------------

    !################################################
    !# NOW DO GRIDS
    call ensure_netcdf_var_exists(ncid, 't', netcdf_real, time_dim, time_id)
    call ensure_netcdf_att_exists(ncid, time_id, 'long_name', 'Time')
    call ensure_netcdf_att_exists(ncid, time_id, 'units', 'L/vt')

    call ensure_netcdf_var_exists(ncid, 'theta0', netcdf_real, fmode_dim, theta0_id)
    call ensure_netcdf_att_exists(ncid, theta0_id, 'long_name', 'Theta_0')

    call ensure_netcdf_var_exists(ncid, 'kx', netcdf_real, nakx_dim, akx_id)
    call ensure_netcdf_att_exists(ncid, akx_id, 'long_name', 'kx rho')

    call ensure_netcdf_var_exists(ncid, 'ky', netcdf_real, naky_dim, aky_id)
    call ensure_netcdf_att_exists(ncid, aky_id, 'long_name', 'ky rho')

    call ensure_netcdf_var_exists(ncid, 'egrid', netcdf_real, negrid_dim, egrid_id)
    call ensure_netcdf_var_exists(ncid, 'lambda', netcdf_real, nlambda_dim, lambda_id)

    call ensure_netcdf_var_exists(ncid, 'theta', netcdf_real, nttot_dim, theta_id)

    !################################################
    !# NOW DO GEOMETRICAL/EQUILIBRIUM DATA
    call ensure_netcdf_var_exists(ncid, 'bmag', netcdf_real, nttot_dim, bmag_id)
    call ensure_netcdf_att_exists(ncid, bmag_id, 'long_name', '|B|(theta)')
    call ensure_netcdf_att_exists(ncid, bmag_id, 'units', 'B_0')

    call ensure_netcdf_var_exists(ncid, 'gradpar', netcdf_real, nttot_dim, gradpar_id)
    call ensure_netcdf_var_exists(ncid, 'gbdrift', netcdf_real, nttot_dim, gbdrift_id)
    call ensure_netcdf_var_exists(ncid, 'gbdrift0', netcdf_real, nttot_dim, gbdrift0_id)
    call ensure_netcdf_var_exists(ncid, 'cvdrift', netcdf_real, nttot_dim, cvdrift_id)
    call ensure_netcdf_var_exists(ncid, 'cvdrift0', netcdf_real, nttot_dim, cvdrift0_id)
    call ensure_netcdf_var_exists(ncid, 'cdrift', netcdf_real, nttot_dim, cdrift_id)
    call ensure_netcdf_var_exists(ncid, 'cdrift0', netcdf_real, nttot_dim, cdrift0_id)

    call ensure_netcdf_var_exists(ncid, 'gds2', netcdf_real, nttot_dim, gds2_id)
    call ensure_netcdf_var_exists(ncid, 'gds21', netcdf_real, nttot_dim, gds21_id)
    call ensure_netcdf_var_exists(ncid, 'gds22', netcdf_real, nttot_dim, gds22_id)
    call ensure_netcdf_var_exists(ncid, 'gds23', netcdf_real, nttot_dim, gds23_id)
    call ensure_netcdf_var_exists(ncid, 'gds24', netcdf_real, nttot_dim, gds24_id)
    call ensure_netcdf_var_exists(ncid, 'gds24_noq', netcdf_real, nttot_dim, gds24_noq_id)
    call ensure_netcdf_var_exists(ncid, 'grho', netcdf_real, nttot_dim, grho_id)
    call ensure_netcdf_var_exists(ncid, 'jacob', netcdf_real, nttot_dim, jacob_id)

    call ensure_netcdf_var_exists(ncid, 'rplot', netcdf_real, nttot_dim, rplot_id)
    call ensure_netcdf_var_exists(ncid, 'zplot', netcdf_real, nttot_dim, zplot_id)
    call ensure_netcdf_var_exists(ncid, 'aplot', netcdf_real, nttot_dim, aplot_id)
    call ensure_netcdf_var_exists(ncid, 'rprime', netcdf_real, nttot_dim, rprime_id)
    call ensure_netcdf_var_exists(ncid, 'zprime', netcdf_real, nttot_dim, zprime_id)
    call ensure_netcdf_var_exists(ncid, 'aprime', netcdf_real, nttot_dim, aprime_id)

    call ensure_netcdf_var_exists(ncid, 'bpol', netcdf_real, nttot_dim, bpol_id)

    !    call ensure_netcdf_var_exists(ncid, 'surfarea', netcdf_real, nttot_dim, surfarea_id)

    call ensure_netcdf_var_exists(ncid, 'q', netcdf_real, q_id)
    call ensure_netcdf_att_exists(ncid, q_id, 'long_name', 'local safety factor')
    call ensure_netcdf_var_exists(ncid, 'eps_trapped', netcdf_real, eps_id)
    call ensure_netcdf_att_exists(ncid, eps_id, 'long_name', 'Trapped fraction 1 - sqrt(bmin/bmax)')
    call ensure_netcdf_var_exists(ncid, 'beta', netcdf_real, beta_id)
    call ensure_netcdf_att_exists(ncid, beta_id, 'long_name', 'reference beta')
    call ensure_netcdf_var_exists(ncid, 'shat', netcdf_real, shat_id)
    call ensure_netcdf_att_exists(ncid, shat_id, 'long_name', '(rho/q) dq/drho')

    call ensure_netcdf_var_exists(ncid, 'drhodpsi', netcdf_real, drhodpsi_id)
    call ensure_netcdf_att_exists(ncid, drhodpsi_id, 'long_name', 'drho/dPsi')
    !------------------------------------------------

    !################################################
    !# NOW WE HANDLE THE DIFFERENT OPTIONAL OUTPUT

    !Phi dependent data
    if (fphi > zero) then

       !Field at end of run
       call ensure_netcdf_var_exists(ncid, 'phi', netcdf_real, final_field_dim, phi_id)
       call ensure_netcdf_att_exists(ncid, phi_id, 'long_name', 'Electrostatic Potential')
       call ensure_netcdf_att_exists(ncid, phi_id, 'idl_name', '!7U!6')
       call ensure_netcdf_att_exists(ncid, phi_id, 'units', 'T/q rho/L')

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

       !Electrostatic fluxes
       if (write_nl_flux) then
          call ensure_netcdf_var_exists(ncid, 'es_heat_par',  netcdf_real, flux_dim, es_heat_par_id)
          call ensure_netcdf_var_exists(ncid, 'es_heat_perp', netcdf_real, flux_dim, es_heat_perp_id)
          call ensure_netcdf_var_exists(ncid, 'es_heat_flux', netcdf_real, flux_dim, es_heat_flux_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom_flux',  netcdf_real, flux_dim, es_mom_flux_id)
          call ensure_netcdf_var_exists(ncid, 'es_part_flux', netcdf_real, flux_dim, es_part_flux_id)
          call ensure_netcdf_var_exists(ncid, 'es_exchange', netcdf_real, flux_dim, es_exchange_id)
          call ensure_netcdf_var_exists(ncid, 'es_exchange_by_k', netcdf_real, fluxk_dim, es_exchange_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_heat_by_k', netcdf_real, fluxk_dim, es_heat_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom_by_k',  netcdf_real, fluxk_dim, es_mom_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_part_by_k', netcdf_real, fluxk_dim, es_part_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_parmom_by_k', netcdf_real, fluxk_dim, es_parmom_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_perpmom_by_k', netcdf_real, fluxk_dim, es_perpmom_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom0_by_k', netcdf_real, fluxk_dim, es_mom0_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'es_mom1_by_k', netcdf_real, fluxk_dim, es_mom1_by_k_id)
       end if

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

       !Normalised phi at end of run
       if (write_eigenfunc) then
          call ensure_netcdf_var_exists(ncid, 'phi_norm',  netcdf_real, final_field_dim, phi_norm_id)
       endif

       !Velocity moment of g used in calculating phi
       if (write_final_antot) then
          call ensure_netcdf_var_exists(ncid, 'antot', netcdf_real, final_field_dim, antot_id)
       endif

       !Field as a function of time
       if (write_phi_t) then
          call ensure_netcdf_var_exists(ncid, 'phi_t', netcdf_real, field_dim, phi_t_id)
          call ensure_netcdf_att_exists(ncid, phi_t_id, 'long_name', 'Electrostatic Potential over time')
       end if
    end if

    !Apar dependent output
    if (fapar > zero) then

       !Field at end of run
       call ensure_netcdf_var_exists(ncid, 'apar', netcdf_real, final_field_dim, apar_id)
       call ensure_netcdf_att_exists(ncid, apar_id, 'long_name', 'Parallel Magnetic Potential')
       call ensure_netcdf_att_exists(ncid, apar_id, 'idl_name', '!6A!9!D#!N!6')

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

       !Apar related fluxes
       if (write_nl_flux) then
          call ensure_netcdf_var_exists(ncid,'apar_heat_flux',netcdf_real, flux_dim, apar_heat_flux_id)
          call ensure_netcdf_var_exists(ncid, 'apar_heat_par', netcdf_real, flux_dim, apar_heat_par_id)
          call ensure_netcdf_var_exists(ncid, 'apar_heat_perp', netcdf_real, flux_dim, apar_heat_perp_id)
          call ensure_netcdf_var_exists(ncid, 'apar_mom_flux',  netcdf_real, flux_dim, apar_mom_flux_id)
          call ensure_netcdf_var_exists(ncid, 'apar_part_flux', netcdf_real, flux_dim, apar_part_flux_id)
          call ensure_netcdf_var_exists(ncid, 'apar_heat_by_k', netcdf_real, fluxk_dim, apar_heat_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'apar_heat_by_x', netcdf_real, fluxx_dim, apar_heat_by_x_id)
          call ensure_netcdf_var_exists(ncid, 'apar_mom_by_k',  netcdf_real, fluxk_dim, apar_mom_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'apar_part_by_k', netcdf_real, fluxk_dim, apar_part_by_k_id)
       end if

       !Normalised eigenfunction
       if (write_eigenfunc) then
          call ensure_netcdf_var_exists(ncid, 'apar_norm', netcdf_real, final_field_dim, apar_norm_id)
       endif

       !Velocity space moment of g used in calculating apar
       if (write_final_antot) then
          call ensure_netcdf_var_exists(ncid, 'antota', netcdf_real, final_field_dim, antota_id)
       endif

       !Field as a function of time
       if (write_apar_t) then
          call ensure_netcdf_var_exists(ncid, 'apar_t', netcdf_real, field_dim, apar_t_id)
          call ensure_netcdf_att_exists(ncid, apar_t_id, 'long_name', 'Parallel Magnetic Potential over Time')
       end if
    end if

    !Bpar related data
    if (fbpar > zero) then

       !Field at end of run
       call ensure_netcdf_var_exists(ncid, 'bpar', netcdf_real, final_field_dim, bpar_id)
       call ensure_netcdf_att_exists(ncid, bpar_id, 'long_name', 'delta B Parallel')
       call ensure_netcdf_att_exists(ncid, bpar_id, 'idl_name', '!6B!9!D#!N!6')

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

       !Bpar related fluxes
       if (write_nl_flux) then
          call ensure_netcdf_var_exists(ncid, 'bpar_heat_flux', netcdf_real, flux_dim, bpar_heat_flux_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_heat_par', netcdf_real, flux_dim, bpar_heat_par_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_heat_perp', netcdf_real, flux_dim, bpar_heat_perp_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_mom_flux', netcdf_real, flux_dim, bpar_mom_flux_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_part_flux', netcdf_real, flux_dim, bpar_part_flux_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_heat_by_k', netcdf_real, fluxk_dim, bpar_heat_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_mom_by_k', netcdf_real, fluxk_dim, bpar_mom_by_k_id)
          call ensure_netcdf_var_exists(ncid, 'bpar_part_by_k', netcdf_real, fluxk_dim, bpar_part_by_k_id)
       end if

       !Normalised eigenfunction
       if (write_eigenfunc) then
          call ensure_netcdf_var_exists(ncid, 'bpar_norm', netcdf_real, final_field_dim, bpar_norm_id)
       endif

       !Velocity space moment of g used in calculating bpar
       if (write_final_antot) then
          call ensure_netcdf_var_exists(ncid, 'antotp', netcdf_real, final_field_dim, antotp_id)
       endif

       !Field as a function of time
       if (write_bpar_t) then
          call ensure_netcdf_var_exists(ncid, 'bpar_t', netcdf_real, field_dim, bpar_t_id)
          call ensure_netcdf_att_exists(ncid, bpar_t_id, 'long_name', 'delta B Parallel over time')
       end if
    end if

!CMR
    if (write_moments) then
       !Only at outboard midplane
       !<DD>Should these really have _t at the end to identify them as being a function of time?
       if (obmid) then
          call ensure_netcdf_var_exists(ncid, 'ntot0', netcdf_real, mom_dim, ntot0_id)
          call ensure_netcdf_att_exists(ncid, ntot0_id, 'long_name', 'Total perturbed density (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'density0', netcdf_real, mom_dim, density0_id)
          call ensure_netcdf_att_exists(ncid, density0_id, 'long_name', 'Nonadiabatic piece of perturbed density (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'upar0', netcdf_real, mom_dim, upar0_id)
          call ensure_netcdf_att_exists(ncid, upar0_id, 'long_name', 'Upar (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'tpar0', netcdf_real, mom_dim, tpar0_id)
          call ensure_netcdf_att_exists(ncid, tpar0_id, 'long_name', 'Tpar (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'tperp0', netcdf_real, mom_dim, tperp0_id)
          call ensure_netcdf_att_exists(ncid, tperp0_id, 'long_name', 'Tperp (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'qparflux0', netcdf_real, mom_dim, qparflux0_id)
          call ensure_netcdf_att_exists(ncid, qparflux0_id, 'long_name', 'Qparflux (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'pperpj10', netcdf_real, mom_dim, pperpj10_id)
          call ensure_netcdf_att_exists(ncid, pperpj10_id, 'long_name', 'Pperpj1 (theta=0) over time')

          call ensure_netcdf_var_exists(ncid, 'qpperpj10', netcdf_real, mom_dim, qpperpj10_id)
          call ensure_netcdf_att_exists(ncid, qpperpj10_id, 'long_name', 'Qpperpj1 (theta=0) over time')
       !Along flux_surface
       else
          call ensure_netcdf_var_exists(ncid, 'ntot_t', netcdf_real, mom_t_dim, ntot_t_id)
          call ensure_netcdf_att_exists(ncid, ntot_t_id, 'long_name', 'Total perturbed density over time')

          call ensure_netcdf_var_exists(ncid, 'density_t', netcdf_real, mom_t_dim, density_t_id)
          call ensure_netcdf_att_exists(ncid, density_t_id, 'long_name', 'Nonadiabatic piece of perturbed density over time')

          call ensure_netcdf_var_exists(ncid, 'upar_t', netcdf_real, mom_t_dim, upar_t_id)
          call ensure_netcdf_att_exists(ncid, upar_t_id, 'long_name', 'Upar over time')

          call ensure_netcdf_var_exists(ncid, 'tpar_t', netcdf_real, mom_t_dim, tpar_t_id)
          call ensure_netcdf_att_exists(ncid, tpar_t_id, 'long_name', 'Tpar over time')

          call ensure_netcdf_var_exists(ncid, 'tperp_t', netcdf_real, mom_t_dim, tperp_t_id)
          call ensure_netcdf_att_exists(ncid, tperp_t_id, 'long_name', 'Tperp over time')

          call ensure_netcdf_var_exists(ncid, 'qparflux_t', netcdf_real, mom_t_dim, qparflux_t_id)
          call ensure_netcdf_att_exists(ncid, qparflux_t_id, 'long_name', 'Qparflux over time')

          call ensure_netcdf_var_exists(ncid, 'pperpj1_t', netcdf_real, mom_t_dim, pperpj1_t_id)
          call ensure_netcdf_att_exists(ncid, pperpj1_t_id, 'long_name', 'Pperpj1 over time')

          call ensure_netcdf_var_exists(ncid, 'qpperpj1_t', netcdf_real, mom_t_dim, qpperpj1_t_id)
          call ensure_netcdf_att_exists(ncid, qpperpj1_t_id, 'long_name', 'Qpperpj1 over time')
       end if
    end if
!CMRend

    if (write_verr) then
!       call ensure_netcdf_var_exists(ncid, 'errest_by_mode', netcdf_real, err_dim, errest_by_mode_id)
!       call ensure_netcdf_var_exists(ncid, 'lpcoef_by_mode', netcdf_real, lp_dim, lpcoef_by_mode_id)
    end if

    !The current value of the scanned parameter
    if (write_scan_parameter) then
       call ensure_netcdf_var_exists(ncid,'scan_parameter_value',netcdf_real,time_dim,current_scan_parameter_value_id)
       call ensure_netcdf_att_exists(ncid,current_scan_parameter_value_id,'long_name','The current value of the scan parameter')
    end if

    !Phase used to normalise eigenfunction/field
    if(write_eigenfunc)then
       call ensure_netcdf_var_exists(ncid, 'phase', netcdf_real, phase_dim, phase_id)
       call ensure_netcdf_att_exists(ncid, phase_id, 'long_name', 'Normalizing phase')
    endif

    !Roughly the field line averaged summed field amplitudes
    call ensure_netcdf_var_exists(ncid, 'phtot', netcdf_real, mode_dim, phtot_id)
    call ensure_netcdf_att_exists(ncid, phtot_id, 'long_name', 'Field amplitude')

    !Frequency and growth rate estimates
    if (write_omega) then
       call ensure_netcdf_var_exists(ncid, 'omega', netcdf_real, omega_dim, omega_id)
       call ensure_netcdf_var_exists(ncid, 'omegaavg', netcdf_real, omega_dim, omegaavg_id)
    end if

    !Summed fluxes
    if (write_nl_flux) then
       call ensure_netcdf_var_exists(ncid, 'hflux_tot', netcdf_real, time_dim, hflux_tot_id)
       call ensure_netcdf_var_exists(ncid, 'vflux_tot', netcdf_real, time_dim, vflux_tot_id)
       call ensure_netcdf_var_exists(ncid, 'zflux_tot', netcdf_real, time_dim, zflux_tot_id)
    end if

    !Parallel electric field -- only good for linear runs
    if(write_final_epar) then
       call ensure_netcdf_var_exists(ncid, 'epar', netcdf_real, final_field_dim, epar_id)
    endif

    !Here we setup the volume/fieldline averaged moment data, currently this triggers
    !both volume averaged, field line averaged and evaulated at theta=0 versions
    if(write_avg_moments) then
       !Volume averaged function of theta
       !<DD>Should we only do this if .not.obmid?
       call ensure_netcdf_var_exists(ncid, 'ntot2', netcdf_real, flux_dim,  ntot2_id)
       call ensure_netcdf_var_exists(ncid, 'ntot2_by_mode', netcdf_real, fluxk_dim, ntot2_by_mode_id)
       call ensure_netcdf_var_exists(ncid, 'tpar2_by_mode', netcdf_real, fluxk_dim, tpar2_by_mode_id)
       call ensure_netcdf_var_exists(ncid, 'tperp2_by_mode', netcdf_real, fluxk_dim, tperp2_by_mode_id)

       !"Volume averaged" data from outboard midplane (really just evaluated at theta=0
       !and then weighted with Fourier coeff factor (typically 0.5).
       !<DD>Should we only do this if obmid?
       call ensure_netcdf_var_exists(ncid, 'ntot20', netcdf_real, flux_dim,  ntot20_id)
       call ensure_netcdf_att_exists(ncid, ntot20_id, 'long_name', 'Density**2 at theta=0')
       call ensure_netcdf_var_exists(ncid, 'ntot20_by_mode', netcdf_real, fluxk_dim, ntot20_by_mode_id)


       !Fieldline averaged data, note this includes phi which is not a moment
       call ensure_netcdf_var_exists(ncid, 'phi00', netcdf_real, loop_phi_dim, phi00_id)
       call ensure_netcdf_var_exists(ncid, 'ntot00', netcdf_real, loop_mom_dim, ntot00_id)
       call ensure_netcdf_var_exists(ncid, 'density00', netcdf_real, loop_mom_dim, density00_id)
       call ensure_netcdf_var_exists(ncid, 'upar00', netcdf_real, loop_mom_dim, upar00_id)
       call ensure_netcdf_var_exists(ncid, 'tpar00', netcdf_real, loop_mom_dim, tpar00_id)
       call ensure_netcdf_var_exists(ncid, 'tperp00', netcdf_real, loop_mom_dim, tperp00_id)
    endif

    !Setup full moment data for end of simulation
    if(write_final_moments)then
       call ensure_netcdf_var_exists(ncid, 'ntot', netcdf_real, final_mom_dim, ntot_id)
       call ensure_netcdf_var_exists(ncid, 'density', netcdf_real, final_mom_dim, density_id)
       call ensure_netcdf_var_exists(ncid, 'upar', netcdf_real, final_mom_dim, upar_id)
       call ensure_netcdf_var_exists(ncid, 'tpar', netcdf_real, final_mom_dim, tpar_id)
       call ensure_netcdf_var_exists(ncid, 'tperp', netcdf_real, final_mom_dim, tperp_id)
       call ensure_netcdf_var_exists(ncid, 'qparflux', netcdf_real, final_mom_dim, qparflux_id)
       call ensure_netcdf_var_exists(ncid, 'pperpj1', netcdf_real, final_mom_dim, pperpj1_id)
       call ensure_netcdf_var_exists(ncid, 'qpperpj1', netcdf_real, final_mom_dim, qpperpj1_id)
    endif

    !Heating rate diagnostics
    if (write_hrate) then
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

    ! RN> not guiding center moments (every nwrite steps)
    if (write_full_moments_notgc .and. .not. (write_moments .and. obmid)) then
       call ensure_netcdf_var_exists(ncid, 'ntot0', netcdf_real, mom_dim, ntot0_id)
       call ensure_netcdf_var_exists(ncid, 'density0', netcdf_real, mom_dim, density0_id)
       call ensure_netcdf_var_exists(ncid, 'upar0', netcdf_real, mom_dim, upar0_id)
       call ensure_netcdf_var_exists(ncid, 'tpar0', netcdf_real, mom_dim, tpar0_id)
       call ensure_netcdf_var_exists(ncid, 'tperp0', netcdf_real, mom_dim, tperp0_id)
    endif
    ! <RN

    !Correlation function (phi only so should probably include if fphi>0)
    if (io_write_corr_extend) then
       call ensure_netcdf_var_exists(ncid, 'phi2_extend',  netcdf_real, phi2_extend_dim, phi2_extend_id)
       call ensure_netcdf_var_exists(ncid, 'phi_corr',  netcdf_real, phi_corr_dim, phi_corr_id)
    end if
    !------------------------------------------------

    !################################################
    !# NOW EXIT DEFINITION MODE
    status = nf90_enddef (ncid)  ! out of definition mode

    !Fill the value of nproc -- helps make sure we can write vars to file
    status = nf90_put_var (ncid, nproc_id, nproc)

# endif
  end subroutine define_vars

  !> Define the vars in the movie file. A fair bit of duplication wrt main
  !! file -- may want to think about merging movie data into main file.
  !! Main reason for using a new file is we can only have one unlimited dimension
  !! per file (with new versions of NETCDF this is no longer true) and we 
  !! probably don't want to write movie data as often as other data.
  !! However, at the moment we precalculate the number of outputs, so it
  !! isn't actually unlimited
  subroutine define_vars_mymovie
# ifdef NETCDF
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: NF90_INT, NF90_GLOBAL
    use netcdf, only: nf90_def_var, nf90_put_att, nf90_enddef
    use netcdf_utils, only: netcdf_real
    implicit none
    integer :: status

    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then
       !Setup field shape
       xmode_dim (1) = nx_dim
       xmode_dim (2) = ny_dim
       xmode_dim (3) = nth_dim
       xmode_dim (4) = time_movie_dim

       !Basic information
       call ensure_netcdf_att_exists(ncid_movie, NF90_GLOBAL, 'title', 'GS2 Simulation x,y,theta Data for Movies')
       call ensure_netcdf_att_exists(ncid_movie, NF90_GLOBAL, 'Conventions','https://sourceforge.net/projects/gyrokinetics/')

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
  subroutine nc_eigenfunc (phase)
# ifdef NETCDF
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> Normalising phase factor (nominally \(\phi\) at the outboard midplane)
    complex, dimension(:,:), intent (in) :: phase
# ifdef NETCDF
    complex, dimension(-ntgrid:ntgrid, ntheta0, naky) :: tmp
    real, dimension(2, ntheta0, naky) :: ri2
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status, ig

    call c2r (phase, ri2)
    status = nf90_put_var (ncid, phase_id, ri2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phase_id)

    if (fphi > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = phi(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, phi_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_norm_id)
    end if

    if (fapar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = apar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, apar_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_norm_id)
    end if

    if (fbpar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = bpar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, bpar_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_norm_id)
    end if
# endif
  end subroutine nc_eigenfunc

  !> Write various moments to netcdf
  subroutine nc_write_moments (nout, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1, ob_midplane) 
    use theta_grid, only: ntgrid
# ifdef NETCDF
    use convert, only: c2r
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use netcdf, only: nf90_put_var
# endif
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
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
    real, dimension (2, ntheta0, naky, nspec) :: ri3
    integer, dimension (6) :: start6, count6
    integer, dimension (5) :: start5, count5
    integer :: status

    if (present(ob_midplane)) then
       if (ob_midplane ) then
          start5(1) = 1
          start5(2) = 1
          start5(3) = 1
          start5(4) = 1
          start5(5) = nout
    
          count5(1) = 2
          count5(2) = ntheta0
          count5(3) = naky
          count5(4) = nspec
          count5(5) = 1

          call c2r (reshape(ntot(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, ntot0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot0_id)
   
          call c2r (reshape(density(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, density0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, density0_id)

          call c2r (reshape(upar(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, upar0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar0_id)
      
          call c2r (reshape(tpar(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, tpar0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar0_id)
      
          call c2r (reshape(tperp(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, tperp0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp0_id)
      
          call c2r (reshape(qparflux(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, qparflux0_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, qparflux0_id)
      
          call c2r (reshape(pperpj1(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, pperpj10_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, pperpj10_id)
      
          call c2r (reshape(qpperpj1(0,:,:,:),(/ntheta0,naky,nspec/)), ri3)
          status = nf90_put_var(ncid, qpperpj10_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, qpperpj10_id)
   
          return

       endif
    endif

    start6(1) = 1
    start6(2) = 1
    start6(3) = 1
    start6(4) = 1
    start6(5) = 1
    start6(6) = nout
    
    count6(1) = 2
    count6(2) = 2*ntgrid+1
    count6(3) = ntheta0
    count6(4) = naky
    count6(5) = nspec
    count6(6) = 1

    call c2r (ntot, ri4)
    status = nf90_put_var(ncid, ntot_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot_t_id)

    call c2r (density, ri4)
    status = nf90_put_var(ncid, density_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density_t_id)

    call c2r (upar, ri4)
    status = nf90_put_var(ncid, upar_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar_t_id)
   
    call c2r (tpar, ri4)
    status = nf90_put_var(ncid, tpar_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar_t_id)
   
    call c2r (tperp, ri4)
    status = nf90_put_var(ncid, tperp_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp_t_id)
   
    call c2r (qparflux, ri4)
    status = nf90_put_var(ncid, qparflux_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, qparflux_t_id)
   
    call c2r (pperpj1, ri4)
    status = nf90_put_var(ncid, pperpj1_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, pperpj1_t_id)
   
    call c2r (qpperpj1, ri4)
    status = nf90_put_var(ncid, qpperpj1_t_id, ri4, start=start6, count=count6)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, qpperpj1_t_id)

# endif
  end subroutine nc_write_moments

  !> Write the fields to netCDF, overwriting existing values
  subroutine nc_final_fields
# ifdef NETCDF
    use convert, only: c2r
    use fields_arrays, only: phinew, aparnew, bparnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: nf90_put_var
    implicit none
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (phinew, ri3)
       status = nf90_put_var (ncid, phi_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_id)
    end if

    if (fapar > zero) then
       call c2r (aparnew, ri3)
       status = nf90_put_var (ncid, apar_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_id)
    end if

    if (fbpar > zero) then
       call c2r (bparnew, ri3)
       status = nf90_put_var (ncid, bpar_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_id)
    end if
# endif
  end subroutine nc_final_fields

  !> Write \(E_\parallel\) to netCDF, overwriting existing values
  subroutine nc_final_epar (epar)
# ifdef NETCDF
    use convert, only: c2r
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use netcdf, only: nf90_put_var
# endif
    implicit none
    complex, dimension (:,:,:), intent (in) :: epar
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    call c2r (epar, ri3)
    status = nf90_put_var (ncid, epar_id, ri3)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, epar_id)
# endif
  end subroutine nc_final_epar

  !> Write the various low-order moments to the output netCDF file,
  !> overwriting existing values
  subroutine nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
# ifdef NETCDF
    use convert, only: c2r
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use netcdf, only: nf90_put_var
# endif
    implicit none
    complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
    complex, dimension (:,:,:,:), intent (in) :: qparflux, pperpj1, qpperpj1
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
    integer :: status

    call c2r (ntot, ri4)
    status = nf90_put_var (ncid, ntot_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot_id)

    call c2r (density, ri4)
    status = nf90_put_var (ncid, density_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density_id)

    call c2r (upar, ri4)
    status = nf90_put_var (ncid, upar_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar_id)

    call c2r (tpar, ri4)
    status = nf90_put_var (ncid, tpar_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar_id)

    call c2r (tperp, ri4)
    status = nf90_put_var (ncid, tperp_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp_id)

    call c2r (qparflux, ri4)
    status = nf90_put_var (ncid, qparflux_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, qparflux_id)

    call c2r (pperpj1, ri4)
    status = nf90_put_var (ncid, pperpj1_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, pperpj1_id)

    call c2r (qpperpj1, ri4)
    status = nf90_put_var (ncid, qpperpj1_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, qpperpj1_id)

# endif
  end subroutine nc_final_moments

  !> Write averaged (volume/field line) moments.
  subroutine nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, ntot20_by_mode, &
       phi00, ntot00, density00, upar00, tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)
# ifdef NETCDF
    use convert, only: c2r
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use netcdf, only: nf90_put_var
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: ntot2, ntot20
    real, dimension (:,:,:), intent (in) :: ntot2_by_mode, ntot20_by_mode
    real, dimension (:,:,:), intent (in) :: tpar2_by_mode, tperp2_by_mode
    complex, dimension (:), intent (in) :: phi00
    complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
# ifdef NETCDF
    real, dimension (2, ntheta0, nspec) :: ri2
    real, dimension (2, ntheta0) :: ri1
    integer, dimension (2) :: start, count
    integer, dimension (3) :: start3, count3
    integer, dimension (4) :: start00, count00, start4, count4
    integer :: status

    start00(1) = 1
    start00(2) = 1
    start00(3) = 1
    start00(4) = nout
    
    count00(1) = 2
    count00(2) = ntheta0
    count00(3) = nspec
    count00(4) = 1

    start3(1) = 1
    start3(2) = 1
    start3(3) = nout
    
    count3(1) = 2
    count3(2) = ntheta0
    count3(3) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    status = nf90_put_var (ncid, ntot2_id, ntot2, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot2_id)
    status = nf90_put_var (ncid, ntot2_by_mode_id, ntot2_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot2_by_mode_id)
    status = nf90_put_var (ncid, tpar2_by_mode_id, tpar2_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar2_by_mode_id)
    status = nf90_put_var (ncid, tperp2_by_mode_id, tperp2_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp2_by_mode_id)

    status = nf90_put_var (ncid, ntot20_id, ntot20, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot20_id)
    status = nf90_put_var (ncid, ntot20_by_mode_id, ntot20_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot20_by_mode_id)

    call c2r (phi00, ri1)
    status = nf90_put_var (ncid, phi00_id, ri1, start=start3, count=count3)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi00_id)
    
    call c2r (ntot00, ri2)
    status = nf90_put_var (ncid, ntot00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot00_id)
    
    call c2r (density00, ri2)
    status = nf90_put_var (ncid, density00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density00_id)
    
    call c2r (upar00, ri2)
    status = nf90_put_var (ncid, upar00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar00_id)
    
    call c2r (tpar00, ri2)
    status = nf90_put_var (ncid, tpar00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar00_id)

    call c2r (tperp00, ri2)
    status = nf90_put_var (ncid, tperp00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp00_id)
# endif
  end subroutine nc_loop_moments

  !> Write the right-hand sides of the field equations to netCDF,
  !> overwriting existing values
  subroutine nc_final_an (antot, antota, antotp)
# ifdef NETCDF
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use netcdf, only: nf90_put_var
# endif
    implicit none
    complex, dimension (:,:,:), intent(in) :: antot, antota, antotp
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (antot, ri3)
       status = nf90_put_var(ncid, antot_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antot_id)
    end if

    if (fapar > zero) then
       call c2r (antota, ri3)
       status = nf90_put_var(ncid, antota_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antota_id)
    end if

    if (fbpar > zero) then
       call c2r (antotp, ri3)
       status = nf90_put_var(ncid, antotp_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antotp_id)
    end if
# endif
  end subroutine nc_final_an

  !> Write the energy exchange to netCDF
  subroutine nc_exchange (nout, exchange, exchange_avg)
# ifdef NETCDF
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    !> Energy exchange
    real, dimension (:,:,:), intent (in) :: exchange
    !> Volume-averaged energy exchange
    real, dimension (:), intent (in) :: exchange_avg
# ifdef NETCDF
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start2, count2
    integer :: status

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    start2(1) = 1
    start2(2) = 1
    start2(3) = 1
    start2(4) = nout

    count2(1) = ntheta0
    count2(2) = naky
    count2(3) = nspec
    count2(4) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_exchange_id, exchange_avg, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_exchange_id)
       status = nf90_put_var (ncid, es_exchange_by_k_id, exchange, start=start2, count=count2)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_exchange_by_k_id)
    end if
# endif
  end subroutine nc_exchange

  !> Write various heat fluxes to netCDF
  subroutine nc_qflux (nout, qheat, qmheat, qbheat, &
       heat_par,  mheat_par,  bheat_par, &
       heat_perp, mheat_perp, bheat_perp, &
       heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot)
# ifdef NETCDF
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: nf90_put_var
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
    integer, dimension (2) :: start, count
    integer, dimension (3) :: start3, count3
    integer, dimension (4) :: start4, count4
    integer :: status

    start3(1) = 1
    start3(2) = 1
    start3(3) = nout

    count3(1) = ntheta0
    count3(2) = nspec
    count3(3) = 1    

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_heat_flux_id, heat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_flux_id)
       status = nf90_put_var (ncid, es_heat_par_id, heat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_par_id)
       status = nf90_put_var (ncid, es_heat_perp_id, heat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_perp_id)
       status = nf90_put_var (ncid, es_heat_by_k_id, qheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_heat_flux_id, mheat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_flux_id)
       status = nf90_put_var (ncid, apar_heat_par_id, mheat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_par_id)
       status = nf90_put_var (ncid, apar_heat_perp_id, mheat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_perp_id)
       status = nf90_put_var (ncid, apar_heat_by_k_id, qmheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_by_k_id)
       status = nf90_put_var (ncid, apar_heat_by_x_id, x_qmflux, start=start3, count=count3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_by_x_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_heat_flux_id, bheat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_flux_id)
       status = nf90_put_var (ncid, bpar_heat_par_id, bheat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_par_id)
       status = nf90_put_var (ncid, bpar_heat_perp_id, bheat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_perp_id)
       status = nf90_put_var (ncid, bpar_heat_by_k_id, qbheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_by_k_id)
    end if
    
    status = nf90_put_var (ncid, hflux_tot_id, hflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, hflux_tot_id)
# endif
  end subroutine nc_qflux

  !> FIXME : Add documentation  
  subroutine nc_pflux (nout, pflux, pmflux, pbflux, &
       part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
# ifdef NETCDF
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: nf90_put_var
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
    real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
    real, intent (in) :: zflux_tot
# ifdef NETCDF
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_part_flux_id, part_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_flux_id)
       status = nf90_put_var (ncid, es_part_by_k_id, pflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_part_flux_id, mpart_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_part_flux_id)
       status = nf90_put_var (ncid, apar_part_by_k_id, pmflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_part_by_k_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_part_flux_id, bpart_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_part_flux_id)
       status = nf90_put_var (ncid, bpar_part_by_k_id, pbflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_part_by_k_id)
    end if
    
    status = nf90_put_var (ncid, zflux_tot_id, zflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, zflux_tot_id)
# endif
  end subroutine nc_pflux

  !> Write toroidal angular momentum flux to netCDF
  subroutine nc_pflux_tormom (nout, pflux_tormom, part_tormom_fluxes)
# ifdef NETCDF
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi
    use netcdf, only: nf90_put_var
# endif
    implicit none
    !> Current output timestep
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux_tormom
    real, dimension(:), intent (in) :: part_tormom_fluxes
# ifdef NETCDF
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout

    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout

    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_part_tormom_flux_id, part_tormom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_tormom_flux_id)
       status = nf90_put_var (ncid, es_part_tormom_by_k_id, pflux_tormom, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_tormom_by_k_id)
    end if

# endif
  end subroutine nc_pflux_tormom

  !> FIXME : Add documentation  
  subroutine nc_vflux (nout, vflux, vmflux, vbflux, &
       mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
       vflux_par, vflux_perp, vflux0, vflux1)
# ifdef NETCDF
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf, only: nf90_put_var
# endif
    implicit none
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
    real, dimension (:,:,:), intent (in) :: vflux_par, vflux_perp
    real, dimension (:,:,:), intent (in) :: vflux0, vflux1
    real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
    real, intent (in) :: vflux_tot
# ifdef NETCDF
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = nf90_put_var (ncid, es_mom_flux_id, mom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_flux_id)
       status = nf90_put_var (ncid, es_mom_by_k_id, vflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_by_k_id)
       status = nf90_put_var (ncid, es_parmom_by_k_id, vflux_par, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_parmom_by_k_id)
       status = nf90_put_var (ncid, es_perpmom_by_k_id, vflux_perp, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_perpmom_by_k_id)
       status = nf90_put_var (ncid, es_mom0_by_k_id, vflux0, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom0_by_k_id)
       status = nf90_put_var (ncid, es_mom1_by_k_id, vflux1, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom1_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_mom_flux_id, mmom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_mom_flux_id)
       status = nf90_put_var (ncid, apar_mom_by_k_id, vmflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_mom_by_k_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_mom_flux_id, bmom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_mom_flux_id)
       status = nf90_put_var (ncid, bpar_mom_by_k_id, vbflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_mom_by_k_id)
    end if
    
    status = nf90_put_var (ncid, vflux_tot_id, vflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vflux_tot_id)
# endif
  end subroutine nc_vflux

  !> Output the poloidal distributions of the fluxes of particles, perpendicular
  !! momentum, parallel momentum, and energy ! JRB
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
       h, hk, omega, omegaavg, woutunits, phitot, write_omega, write_hrate)
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
    logical, intent(in) :: write_omega, write_hrate
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

    status = nf90_put_var (ncid, time_id, time, start=(/ nout /))

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

    if (write_hrate) then
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
  subroutine nc_loop_fullmom (nout, time, &
       & ntot0, density0, upar0, tpar0, tperp0)
# ifdef NETCDF
    use kt_grids, only: naky, nakx => ntheta0
    use species, only: nspec
    use convert, only: c2r
    use netcdf, only: nf90_put_var
# endif
    implicit none
    integer, intent (in) :: nout
    real, intent (in) :: time
    complex, dimension(:,:,:), intent(in) :: ntot0, density0, upar0, tpar0, tperp0
# ifdef NETCDF
    real, dimension (2, nakx, naky, nspec) :: ri3
    integer, dimension (5) :: startmom, countmom
    integer :: status

    status = nf90_put_var (ncid, time_id, time, start=(/nout/))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, time_id)

    startmom(1) = 1;    countmom(1)=2
    startmom(2) = 1;    countmom(2)=nakx
    startmom(3) = 1;    countmom(3)=naky
    startmom(4) = 1;    countmom(4)=nspec
    startmom(5) = nout; countmom(5)=1

    call c2r (ntot0, ri3)
    status = nf90_put_var (ncid, ntot0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot0_id)
    call c2r (density0, ri3)
    status = nf90_put_var (ncid, density0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density0_id)
    call c2r (upar0, ri3)
    status = nf90_put_var (ncid, upar0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar0_id)
    call c2r (tpar0, ri3)
    status = nf90_put_var (ncid, tpar0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar0_id)
    call c2r (tperp0, ri3)
    status = nf90_put_var (ncid, tperp0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp0_id)
# endif
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
  subroutine nc_species
# ifdef NETCDF
    use species, only: spec
    use netcdf, only: nf90_put_var
    implicit none
    integer :: status

    status = nf90_put_var (ncid, charge_id, spec%z)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, charge_id)
    status = nf90_put_var (ncid, mass_id, spec%mass)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, mass_id)
    status = nf90_put_var (ncid, dens_id, spec%dens)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, dens_id)
    status = nf90_put_var (ncid, temp_id, spec%temp)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, temp_id)
    status = nf90_put_var (ncid, tprim_id, spec%tprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tprim_id)
    status = nf90_put_var (ncid, fprim_id, spec%fprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, fprim_id)
    status = nf90_put_var (ncid, uprim_id, spec%uprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, uprim_id)
    status = nf90_put_var (ncid, uprim2_id, spec%uprim2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, uprim2_id)
    status = nf90_put_var (ncid, vnewk_id, spec%vnewk)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vnewk_id)
    status = nf90_put_var (ncid, spec_type_id, spec%type)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, spec_type_id)

# endif
  end subroutine nc_species
  
  !> Write various geometric quantities to output netCDF file
  subroutine nc_geo
# ifdef NETCDF
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob, &
         shat, drhodpsi, eps_trapped, cdrift, cdrift0, qval, rplot, &
         zplot, aplot, rprime, zprime, aprime, bpol !, surfarea
!CMR add qval and beta here too, as they are both useful
    use run_parameters, only: beta
    use netcdf, only: nf90_put_var
    implicit none
    integer :: status

    status = nf90_put_var (ncid, bmag_id, bmag)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bmag_id)
    status = nf90_put_var (ncid, gradpar_id, gradpar)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gradpar_id)
    status = nf90_put_var (ncid, gbdrift_id, gbdrift)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gbdrift_id)
    status = nf90_put_var (ncid, gbdrift0_id, gbdrift0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gbdrift0_id)
    status = nf90_put_var (ncid, cvdrift_id, cvdrift)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cvdrift_id)
    status = nf90_put_var (ncid, cvdrift0_id, cvdrift0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cvdrift0_id)
    status = nf90_put_var (ncid, cdrift_id, cdrift)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cdrift_id)
    status = nf90_put_var (ncid, cdrift0_id, cdrift0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cdrift0_id)
    status = nf90_put_var (ncid, gds2_id, gds2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds2_id)
    status = nf90_put_var (ncid, gds21_id, gds21)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds21_id)
    status = nf90_put_var (ncid, gds22_id, gds22)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds22_id)
    status = nf90_put_var (ncid, gds23_id, gds23)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds23_id)
    status = nf90_put_var (ncid, gds24_id, gds24)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds24_id)
    status = nf90_put_var (ncid, gds24_noq_id, gds24_noq)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds24_noq_id)
    status = nf90_put_var (ncid, grho_id, grho)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, grho_id)
    status = nf90_put_var (ncid, jacob_id, jacob)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, jacob_id)

    status = nf90_put_var (ncid, rplot_id, rplot)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, rplot_id)
    status = nf90_put_var (ncid, zplot_id, zplot)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, zplot_id)
    status = nf90_put_var (ncid, aplot_id, aplot)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aplot_id)
    status = nf90_put_var (ncid, rprime_id, rprime)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, rprime_id)
    status = nf90_put_var (ncid, zprime_id, zprime)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, zprime_id)
    status = nf90_put_var (ncid, aprime_id, aprime)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aprime_id)

    status = nf90_put_var (ncid, bpol_id, bpol)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpol_id)

    !    status = nf90_put_var (ncid, surfarea_id, surfarea)
    !    if (status /= NF90_NOERR) call netcdf_error (status, ncid, surfarea_id)

    status = nf90_put_var (ncid, beta_id, beta)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, beta_id)
    status = nf90_put_var (ncid, q_id, qval)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, q_id)
    status = nf90_put_var (ncid, shat_id, shat)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, shat_id)
    status = nf90_put_var (ncid, eps_id, eps_trapped)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, eps_id)
    status = nf90_put_var (ncid, drhodpsi_id, drhodpsi)   
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, drhodpsi_id)
# endif
  end subroutine nc_geo
end module gs2_io
