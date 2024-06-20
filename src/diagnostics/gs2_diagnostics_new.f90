!> A module for calculating and writing gs2 outputs. It can write
!> these outputs both to a netcdf file <run_name>.cdf or to ascii text
!> files. It is controlled via the namelist diagnostics_config.
!> This module is intended to replace the old gs2_diagnostics module with a
!> simpler and more structured interface.
module gs2_diagnostics_new
  use diagnostics_config, only: diagnostics_type
  implicit none

  private
  
  public :: init_gs2_diagnostics_new
  public :: finish_gs2_diagnostics_new
  public :: reset_averages_and_counters
  public :: run_diagnostics
  public :: gnostics

  !> Do we have parallel IO capability in the build. This is always disabled
  !> for now whilst we actually implement this feature
  logical, parameter :: parallel_io_capable = .false.

  type(diagnostics_type) :: gnostics 

  logical, parameter :: debug = .false.

contains

  !> Read namelist diagnostics_config, initialise submodules,
  !> open output file 'run_name.cdf' and create dimensions.
  subroutine init_gs2_diagnostics_new(header)
    use kt_grids, only: nx, ny, naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use theta_grid, only: ntgrid
    use gs2_transforms, only: init_transforms
    use diagnostics_config, only: init_diagnostics_config
    use diagnostics_fluxes, only: init_diagnostics_fluxes
    use diagnostics_omega, only: init_diagnostics_omega
    use diagnostics_velocity_space, only: init_diagnostics_velocity_space
    use diagnostics_heating, only: init_diagnostics_heating
    use diagnostics_ascii, only: init_diagnostics_ascii
    use diagnostics_antenna, only: init_diagnostics_antenna
    use diagnostics_nonlinear_convergence, only: init_nonlinear_convergence
    use diagnostics_zonal_transfer, only: init_diagnostics_transfer
    use collisional_heating, only: init_collisional
    use collisions, only: heating
    use nonlinear_terms, only: nonlin
    use gs2_save, only: save_many
    use file_utils, only: run_name, error_unit
    use mp, only: proc0, broadcast, mp_abort
    use kt_grids, only: naky, aky
    use gs2_diagnostics, only: check_restart_file_writeable
    use unit_tests, only: debug_message
    use standard_header, only: standard_header_type
    use gs2_metadata, only: create_metadata
    use run_parameters, only: user_comments
    use constants, only: run_name_size
    use neasyf, only: neasyf_open
    use gs2_io, only: define_dims, nc_norms, nc_species, nc_geo, save_input, get_dim_length, nc_grids_mymovie
    implicit none
    !> Header for files with build and run information
    type(standard_header_type), intent(in) :: header
    logical :: ex, accelerated
    character(run_name_size) :: filename
    
    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new starting')
    call init_diagnostics_config(gnostics)
    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new initialized config')
    call check_parameters
    call check_restart_file_writeable(gnostics%file_safety_check, &
                                      gnostics%save_for_restart, &
                                      gnostics%save_distfn)
    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new  checked restart file')
    
    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new inialized vol avgs')
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !! Adjust other modules
    !!!!!!!!!!!!!!!!!!!!!!!
    save_many = gnostics%save_many

    if (.not. gnostics%write_any) return
    
    gnostics%parallel = .false.
    if (gnostics%enable_parallel) then
       !> We should replace parallel_io_capable below with
       !> [[build_config:gs2_has_netcdf_parallel]] once we actually
       !> want to allow parallel diagnostics.
       if (parallel_io_capable) then
          gnostics%parallel = .true.
       else
          if (proc0) write (*,*) "WARNING: you have selected &
               & enable_parallel but this build does not have &
               & parallel capability."
       end if
    end if
    
    gnostics%user_time_old = 0.0
    
    ! fluxfac is used for summing fields, fluxes etc over ky
    ! Mostly this is not needed, since the average_ky routine in 
    ! volume_averages takes care of the factor... you only need it
    ! if you are manually summing something over ky
    allocate(gnostics%fluxfac(naky))
    gnostics%fluxfac = 0.5
    !<DD>This is only correct if running in box mode surely?
    !    I think this should be if(aky(1)==0.0) fluxfac(1)=1.0 but I may be wrong
    if(aky(1)==0.0) gnostics%fluxfac(1) = 1.0

    if (gnostics%write_heating .and. .not. heating) then
       if (proc0) write(*,'("Warning: Disabling write_heating as collisions:heating is false.")')
       gnostics%write_heating = .false.
    else if (heating .and. .not. &
       (gnostics%write_heating .or. gnostics%write_collisional)) then
       heating = .false.
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Open Text Files (if required)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (proc0) call set_ascii_file_switches
    if (proc0) call init_diagnostics_ascii(gnostics%ascii_files)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise submodules
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_diagnostics_fluxes(gnostics)
    call init_diagnostics_omega(gnostics)
    !if (gnostics%write_max_verr) gnostics%write_verr = .true.
    call init_diagnostics_velocity_space(gnostics)
    call init_diagnostics_antenna(gnostics)
    call init_diagnostics_transfer(gnostics)
    call init_collisional(gnostics)

    if (gnostics%write_heating) call init_diagnostics_heating(gnostics)
    if (nonlin.and.gnostics%use_nonlin_convergence) call init_nonlinear_convergence(gnostics)
    
    filename = trim(trim(run_name)//'.out.nc')

    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new opening file')
    if (gnostics%parallel.or.proc0) then 
       inquire(file=trim(run_name)//'.out.nc', exist=ex)
       if (gnostics%append_old .and. ex) then
         gnostics%appending=.true.
         gnostics%file_id = neasyf_open(trim(filename), "rw")
         gnostics%nout = get_dim_length(gnostics%file_id, "t")
       else
         gnostics%appending=.false.
         gnostics%file_id = neasyf_open(trim(filename), "w")
         call create_metadata(gnostics%file_id, "GS2 Simulation Data", header, user_comments)
         call debug_message(gnostics%verbosity, &
           'gs2_diagnostics_new::init_gs2_diagnostics_new written metadata')

         ! Write constants/parameters
         call define_dims(gnostics%file_id, gnostics%write_correlation_extend)
         call nc_norms(gnostics%file_id)
         call nc_species(gnostics%file_id)
         call nc_geo(gnostics%file_id)
         call save_input(gnostics%file_id)
         gnostics%nout = 1

         if (gnostics%make_movie) then
           call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
           call nc_grids_mymovie(gnostics%file_id)
         end if
       end if
    end if
    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new opened file')

    ! Vital that other procs know we are appending
    call broadcast(gnostics%appending)
    call broadcast(gnostics%nout)

    call debug_message(gnostics%verbosity, &
      'gs2_diagnostics_new::init_gs2_diagnostics_new finished')
    
  end subroutine init_gs2_diagnostics_new

  !> FIXME : Add documentation  
  subroutine check_parameters
    use run_parameters, only: fapar
    use file_utils, only: error_unit
    use mp, only: proc0
    implicit none
    if ((gnostics%print_line .or. gnostics%write_line) .and. .not.(gnostics%write_fields.and.gnostics%write_omega)) then 
       if (proc0) write (error_unit(), *) 'print_line and write_line require both write_fields and write_omega... enabling'
       gnostics%write_fields = .true.
       gnostics%write_omega = .true.
    end if
    if ((gnostics%print_flux_line .or. gnostics%write_flux_line) .and. .not.gnostics%write_fields) then 
       if (proc0) write (error_unit(), *) 'print_flux_line and write_flux_line require both write_fields ... enabling'
       gnostics%write_fields = .true.
    end if
    if (gnostics%write_jext .and. .not. fapar .gt. epsilon(0.0)) then
       if (proc0) write (*,*) "ERROR: it doesn't make sense to switch on write_jext without apar"
       !SHOULD THIS BE MP_ABORT? COULD WE NOT JUST DISABLE THE DIAGNOSTIC?
       stop 1
    end if
    
  end subroutine check_parameters

  !> This subroutine determines which ascii output files are enabled
  !! (i.e., opened, flushed at each write, and then closed).
  !! If an ascii file is not enabled here, writing to it will 
  !! cause some indeterminate unpleasant behaviour
  !!
  !! Note that the .out file is always enabled
  subroutine set_ascii_file_switches
    implicit none
    gnostics%ascii_files%write_to_out   = gnostics%write_ascii
    !gnostics%ascii_files%write_to_fields = gnostics%write_fields  .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat   = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat2  = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_lpc    = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres   = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres2  = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_cres   = gnostics%write_cerr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_dist   = gnostics%write_g    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_yxdist   = gnostics%write_gyx    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_phase  = gnostics%write_cross_phase .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_jext   = gnostics%write_jext    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_parity = gnostics%write_parity  .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_eigenfunc = gnostics%write_eigenfunc .and.  gnostics%write_ascii
  end subroutine set_ascii_file_switches

  !> Close the output file and deallocate arrays
  subroutine finish_gs2_diagnostics_new
    use diagnostics_fluxes, only: finish_diagnostics_fluxes
    use diagnostics_omega, only: finish_diagnostics_omega
    use diagnostics_heating, only: finish_diagnostics_heating
    use diagnostics_ascii, only: finish_diagnostics_ascii
    use diagnostics_config, only: finish_diagnostics_config
    use diagnostics_antenna, only: finish_diagnostics_antenna
    use diagnostics_nonlinear_convergence, only: finish_nonlinear_convergence
    use diagnostics_velocity_space, only: finish_diagnostics_velocity_space
    use collisional_heating, only: finish_collisional
    use gs2_diagnostics, only: save_restart_dist_fn, do_write_geom, do_write_fyx, do_write_f
    use nonlinear_terms, only: nonlin
    use mp, only: proc0
    use fields_arrays, only: phinew, bparnew
    use unit_tests, only: debug_message
    use neasyf, only: neasyf_close
    implicit none
    integer, parameter :: verb=3
    if (.not. gnostics%write_any) return
    
    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new &
         & calling save_restart_dist_fn')
    call save_restart_dist_fn(gnostics%save_for_restart, &
         gnostics%save_distfn, &
         gnostics%save_glo_info_and_grids, &
         gnostics%save_velocities, &
         gnostics%user_time)

    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new &
         & calling run_old_final_routines')
    call run_old_final_routines
    
    deallocate(gnostics%fluxfac)
    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new finishing submodules')
    call finish_diagnostics_fluxes
    call finish_diagnostics_omega
    call finish_collisional
    call finish_diagnostics_antenna(gnostics)
    call finish_diagnostics_velocity_space()
    
    if (nonlin.and.gnostics%use_nonlin_convergence) call finish_nonlinear_convergence()
    if (gnostics%write_heating) call finish_diagnostics_heating(gnostics)
    if (gnostics%parallel .or. proc0) then
       if(proc0.and.debug) write(*,*) "Closing new diagnostics"
       call neasyf_close(gnostics%file_id)
    end if
    if (proc0) call finish_diagnostics_ascii(gnostics%ascii_files)
    
    ! Random stuff that needs to be put in properly or removed
    ! Note we don't use gnostics%ascii_files%write_to_yxdist here
    ! because do_write_fyx is collective but gnostics%ascii_files%write_to_yxdist
    ! is only set on proc0
    if (gnostics%write_ascii .and. gnostics%write_gyx) &
         call do_write_fyx (gnostics%ascii_files%yxdist, phinew, bparnew)
    if (gnostics%write_ascii .and. gnostics%write_g) &
         call do_write_f (gnostics%ascii_files%dist)
    if (gnostics%write_ascii) call do_write_geom()
    
    call finish_diagnostics_config(gnostics)
  end subroutine finish_gs2_diagnostics_new

  !> FIXME : Add documentation  
  subroutine run_diagnostics_to_be_updated
    use fields_arrays, only: phinew, bparnew
    use gs2_diagnostics, only: do_write_fyx, do_write_f
    implicit none
    integer :: nwrite_large

    nwrite_large = gnostics%nwrite*gnostics%nwrite_mult
    if (mod(gnostics%istep,nwrite_large) /= 0) return

    ! Random stuff that needs to be put in properly or removed
    ! Note we don't use gnostics%ascii_files%write_to_yxdist here
    ! because do_write_fyx is collective but gnostics%ascii_files%write_to_yxdist
    ! is only set on proc0
    if (gnostics%write_ascii .and. gnostics%write_gyx) &
         call do_write_fyx (gnostics%ascii_files%yxdist, phinew, bparnew)
    if (gnostics%write_ascii .and. gnostics%write_g) &
         call do_write_f (gnostics%ascii_files%dist)
  end subroutine run_diagnostics_to_be_updated

  !> Calculates and write the QL flux metric to netcdf
  subroutine write_ql_metric(gnostics)
    use gs2_diagnostics, only: calculate_simple_quasilinear_flux_metric_by_k
    use gs2_diagnostics, only: calculate_instantaneous_omega
    use kt_grids, only: ntheta0, naky
    use gs2_io, only: nc_write_ql_metric
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    real, dimension(ntheta0, naky) :: growth_rates
    if (.not. gnostics%writing) return

    growth_rates = aimag(calculate_instantaneous_omega(&
         ig = gnostics%igomega, tolerance = gnostics%omegatol))
    call nc_write_ql_metric(gnostics%file_id, gnostics%nout, &
         ql_metric = calculate_simple_quasilinear_flux_metric_by_k(growth_rates))
  end subroutine write_ql_metric
  
  !> Create or write all variables according to the value of istep:
  !! istep=-1 --> Create all netcdf variables
  !! istep=0 --> Write constant arrays/parameters (e.g. aky) and initial values
  !! istep>0 --> Write variables
  subroutine run_diagnostics(istep_in, exit, force)
    use gs2_time, only: user_time, tunits
    use mp, only: proc0
    use diagnostics_zonal_transfer, only: write_zonal_transfer, calculate_zonal_transfer
    use collisional_heating, only: write_collisional, calculate_collisional   
    use diagnostics_printout, only: print_flux_line, print_line
    use diagnostics_printout, only: write_flux_line, write_line
    use diagnostics_fluxes, only: calculate_fluxes
    use diagnostics_fields, only: write_fields, write_movie
    use diagnostics_fields, only: write_eigenfunc
    use diagnostics_moments, only: write_moments, write_full_moments_notgc
    use diagnostics_omega, only: calculate_omega, write_omega
    use diagnostics_velocity_space, only: write_velocity_space_checks
    use diagnostics_velocity_space, only: write_collision_error
    use diagnostics_heating, only: calculate_heating, write_heating
    use diagnostics_nonlinear_convergence, only: check_nonlin_convergence
    use diagnostics_turbulence, only: write_cross_phase
    use diagnostics_antenna, only: write_jext, write_lorentzian
    use diagnostics_ascii, only: flush_output_files
    use collisions, only: vary_vnew
    use nonlinear_terms, only: nonlin
    use species, only: spec, has_electron_species
    use unit_tests, only: debug_message
    use mp,  only: broadcast
    use neasyf, only: neasyf_write
    use optionals, only: get_option_with_default
    use run_parameters, only: nstep, wstar_units
    use gs2_io, only: get_dim_length, nc_sync, ky_dim, time_dim
    use gs2_diagnostics, only: do_dump_fields_periodically, do_write_parity, do_write_nl_flux_dist, &
         do_write_correlation_extend, do_write_correlation, do_write_symmetry, do_write_pflux_sym

    implicit none
    integer, intent(in) :: istep_in
    logical, intent(inout) :: exit
    logical, intent(in), optional :: force
    integer, parameter :: verb=3
    integer :: istep
    integer, save :: istep_last = -1
    logical :: do_force
    if (.not. gnostics%write_any) return
    
    call broadcast(exit)

    call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics starting')
    do_force = get_option_with_default(force, .false.)
    
    gnostics%exit = exit
    
    ! If parallel, then everybody writes to netcdf,
    ! otherwise, only proc0
    if (gnostics%parallel .or. proc0) then
      gnostics%create = (istep_in==-1) .and. .not. gnostics%appending
      gnostics%writing = .true.
    else
      gnostics%create = .false.
      gnostics%writing = .false.
    end if

    ! Now that we've used istep to work out what operations we want
    ! to do, ensure the istep value is valid.
    istep = istep_in
    if(istep_in == -1) istep = 0
    gnostics%istep = istep
    
    ! Sets whether field-like arrays are assumed
    ! to be distributed across processes
    ! This line is a temporary placeholder
    ! till distributed fields are up and running
    gnostics%distributed = gnostics%parallel
    
    gnostics%calculate_fluxes = (gnostics%write_fluxes &
         .or.  gnostics%print_flux_line &
         .or.  gnostics%write_flux_line)

    gnostics%user_time = user_time
    if (istep .eq. 0) gnostics%start_time = user_time
    
    if (istep > 0) then
       call calculate_omega(gnostics)
       if (gnostics%write_heating) call calculate_heating (gnostics)
    end if
    call broadcast(gnostics%exit)

    call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics calculated &
      & omega and heating')

    if ((istep /= istep_last) .and. (mod(istep, gnostics%nwrite) == 0 .or. gnostics%exit .or. do_force)) then
       ! If istep_in = -1 (to indicate we're setting things up) then we enter
       ! this block (by forcing istep=0), but only to define variables - not to write them.
       ! We intend to write them on the next call when istep_in = 0, but if we were to
       ! record istep_last = istep when istep_in = -1 then we end up with istep_last = 0
       ! and we would actually skip the real istep_in = 0 call.
       if (istep_in > -1) istep_last = istep

       if (gnostics%writing) then
         call neasyf_write(gnostics%file_id, time_dim, user_time, dim_names=[time_dim], start=[gnostics%nout])
         if (wstar_units) then
            call neasyf_write(gnostics%file_id, "t_wstar", user_time * tunits, dim_names=[ky_dim, time_dim], &
                 start=[1, gnostics%nout], long_name="Time (wstar)", units="L/vt")
         end if
       end if

       gnostics%user_time = user_time
       gnostics%vary_vnew_only = .false.
       call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics starting write sequence')
       if (gnostics%write_omega)  call write_omega(gnostics)
       call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics written omega')
       if (gnostics%write_fields) call write_fields(gnostics)
       call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics written fields')
       if (gnostics%dump_fields_periodically) call do_dump_fields_periodically(gnostics%user_time)
       if (gnostics%write_ql_metric)  call write_ql_metric(gnostics)
       if (gnostics%calculate_fluxes) call calculate_fluxes(gnostics) ! NB  also writes fluxes if on
       call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics calculated fluxes')
       if (gnostics%write_symmetry) then
         call do_write_symmetry(gnostics%file_id, gnostics%nout)
         call do_write_pflux_sym(gnostics%file_id, gnostics%nout)
       end if
       if (gnostics%write_nl_flux_dist) call do_write_nl_flux_dist(gnostics%file_id, gnostics%nout)
       if (gnostics%write_parity) call do_write_parity(gnostics%user_time, gnostics%ascii_files%parity, gnostics%write_ascii)
       if (gnostics%write_verr) call write_velocity_space_checks(gnostics)
       if (gnostics%write_cerr) call write_collision_error(gnostics) ! NB only ascii atm
       call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics writing moments')
       if (gnostics%write_moments) call write_moments(gnostics)
       if (gnostics%write_full_moments_notgc) call write_full_moments_notgc(gnostics)
       if (gnostics%make_movie) call write_movie(gnostics)
       if (gnostics%write_heating) call write_heating(gnostics)
       if (nonlin.and.gnostics%use_nonlin_convergence) call check_nonlin_convergence(gnostics)
       if (gnostics%write_cross_phase.and.has_electron_species(spec)) call write_cross_phase(gnostics)
       if (gnostics%write_jext) call write_jext(gnostics)
       if (gnostics%write_correlation) call do_write_correlation(gnostics%file_id, gnostics%nout)
       if (gnostics%write_correlation_extend &
            .and. istep > nstep/4 &
            .and. mod(istep, gnostics%nwrite_mult * gnostics%nwrite)==0) &
            call do_write_correlation_extend(gnostics%file_id, gnostics%user_time, gnostics%user_time_old)
       if (gnostics%write_lorentzian) call write_lorentzian(gnostics)
       if (gnostics%write_eigenfunc) call write_eigenfunc(gnostics)
       
       if (gnostics%print_line) call print_line(gnostics)
       if (gnostics%write_line) call write_line(gnostics)
       if (proc0) then
          if (gnostics%print_flux_line) call print_flux_line(gnostics)
          if (gnostics%write_flux_line) call write_flux_line(gnostics)
       end if

       if (gnostics%write_zonal_transfer) then
          call calculate_zonal_transfer(gnostics)
          call write_zonal_transfer(gnostics)
       end if 

       if (gnostics%write_collisional) then  !< new diagnostic, calulate and write in new variable
          call calculate_collisional()
          call write_collisional(gnostics)
       end if

       call run_diagnostics_to_be_updated

       ! Don't sync movie file because it's the same as the main file
       if (gnostics%writing) call nc_sync(gnostics%file_id, gnostics%nout, -1, -1, gnostics%nc_sync_freq)
       if (proc0 .and. gnostics%write_ascii) call flush_output_files(gnostics%ascii_files)
       
       ! Update time used for time averages
       gnostics%user_time_old = gnostics%user_time

       gnostics%nout = gnostics%nout + 1
    else if (mod(istep, gnostics%ncheck).eq.0) then
       ! These lines cause the automated checking of velocity space resolution
       ! and correction by varying collisionality
       gnostics%vary_vnew_only = .true.
       if (gnostics%write_verr .and. vary_vnew) call write_velocity_space_checks(gnostics)
    end if

    call debug_message(verb, 'gs2_diagnostics_new::run_diagnostics finished')
    exit = gnostics%exit
    call broadcast(exit)
  end subroutine run_diagnostics

  !> Reset cumulative flux and heating averages
  !! that are used, e.g. for Trinity.
  !! Does not at the moment apply to average
  !! growth rates.
  subroutine reset_averages_and_counters
    use nonlinear_terms, only: nonlin
    use diagnostics_nonlinear_convergence, only: &
      dnc_reset => reset_averages_and_counters
    use diagnostics_fluxes, only: &
      fluxes_reset => reset_averages_and_counters
    use diagnostics_heating, only: &
      heating_reset => reset_averages_and_counters
    use gs2_time, only: user_time

    if (nonlin.and.gnostics%use_nonlin_convergence) call dnc_reset()
    call fluxes_reset(gnostics)
    call heating_reset(gnostics)

    gnostics%start_time = user_time

  end subroutine reset_averages_and_counters

  !> FIXME : Add documentation  
  subroutine run_old_final_routines
    use diagnostics_final_routines,only: do_write_final_fields
    use diagnostics_final_routines,only: do_write_kpar
    use diagnostics_final_routines,only: do_write_final_epar
    use diagnostics_final_routines,only: do_write_final_db
    use diagnostics_final_routines,only: do_write_final_moments
    use diagnostics_final_routines,only: do_write_final_antot
    use diagnostics_final_routines,only: do_write_gs
    use diagnostics_final_routines,only: init_par_filter
    use diagnostics_fields, only: get_phi0
    use nonlinear_terms, only: nonlin
    use antenna, only: dump_ant_amp
    use mp, only: proc0
    use kt_grids, only: ntheta0, naky
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (ntheta0, naky) :: phi0
    integer, parameter :: verb=3
    
    if(gnostics%write_kpar.or.gnostics%write_gs) call init_par_filter

    if (proc0) then
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_final_fields')
       if (gnostics%write_final_fields) then
          call do_write_final_fields(gnostics%write_ascii, write_netcdf=.false., file_id=gnostics%file_id)
       end if
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_kpar')
       if (gnostics%write_kpar) call do_write_kpar(gnostics%write_ascii)
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_final_epar')
       if (gnostics%write_final_epar) then
          call do_write_final_epar(gnostics%write_ascii, write_netcdf=.false., file_id=gnostics%file_id)
       end if
       
       ! definition here assumes we are not using wstar_units
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_final_db')
       if (gnostics%write_final_db) call do_write_final_db(gnostics%write_ascii)
    end if

    phi0 = get_phi0()
    if (gnostics%write_final_moments) then
      call do_write_final_moments(phi0, use_normalisation=gnostics%write_eigenfunc, &
                                  write_text=gnostics%write_ascii, write_netcdf=.false., &
                                  file_id=gnostics%file_id)
    end if

    call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
      & calling do_write_final_antot')
    
    if (gnostics%write_final_antot) then
      call do_write_final_antot(gnostics%write_ascii, write_netcdf=.false., file_id=gnostics%file_id)
    end if

    if (proc0) call dump_ant_amp
    
    if (nonlin.and.gnostics%write_gs) call do_write_gs(gnostics%write_ascii)
    call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & finished')
    
  end subroutine run_old_final_routines
end module gs2_diagnostics_new
