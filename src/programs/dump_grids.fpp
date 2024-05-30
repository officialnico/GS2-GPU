!> For the given input file initialises the main grids
!> and geometry and writes these to file (netcdf or binary).
program dump_grids
  use job_manage, only: job_fork, time_message
  use mp, only: init_mp, proc0, nproc, broadcast
  use file_utils, only: init_file_utils, run_name, run_name_target
  use kt_grids, only: init_kt_grids
  use theta_grid, only: init_theta_grid
  use le_grids, only: init_le_grids
  use species, only: init_species
  implicit none
  logical :: list
  logical, parameter :: quiet=.false.
  character (len=500) :: filename
  real :: time_measure(2)

  ! Initialize message passing
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
     !  whether we are doing a run or a list of runs.
     call init_file_utils (list, name="gs")
  end if

  call broadcast (list)

  ! If given a list of jobs, fork
  if (list) call job_fork

  if (proc0) run_name_target = trim(run_name)
  call broadcast (run_name_target)
  if (.not. proc0) run_name => run_name_target

  !Initialise
  time_measure=0.
  !/Theta
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init theta_grids")',advance='no')
  endif
  call time_message(.false.,time_measure,'init-theta')
  call init_theta_grid
  call time_message(.false.,time_measure,'init-theta')
  if(proc0.and.(.not.quiet)) then
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/KT
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init kt_grids   ")',advance='no')
  endif
  call time_message(.false.,time_measure,'init-kt')
  call init_kt_grids
  call time_message(.false.,time_measure,'init-kt')
  if(proc0.and.(.not.quiet)) then
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/SPECIES
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init species   ")',advance='no')
  endif
  call time_message(.false.,time_measure,'init-species')
  call init_species
  call time_message(.false.,time_measure,'init-species')
  if(proc0.and.(.not.quiet)) then
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/LE
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init le_grids   ")',advance='no')
  endif
  call time_message(.false.,time_measure,'init-le')
  call init_le_grids
  call time_message(.false.,time_measure,'init-le')
  if(proc0.and.(.not.quiet)) then
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

  !Now write to file
  if(proc0.and.(.not.quiet)) then
     write(*,'("Write to file   ")',advance='no')
  endif
  call time_message(.false.,time_measure,'write-file')
  filename=trim(adjustl(run_name))
  if(proc0) call write_grids_to_file(filename)
  call time_message(.false.,time_measure,'write-file')
  if(proc0.and.(.not.quiet)) then
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

contains

#ifdef NETCDF
!Prefer to write to netcdf file if possible

  !> FIXME : Add documentation
  subroutine write_grids_to_file(fname)
    use neasyf, only: neasyf_open, neasyf_default_compression, neasyf_dim, neasyf_write, neasyf_close
    use kt_grids, only: aky, akx, theta0
    use theta_grid, only: theta
    use le_grids, only: al, energy, wl, w, negrid
    use species, only: nspec
    use gs2_metadata, only: create_metadata
    use gs2_io, only: nc_geo, nc_species
    use gs2_io, only: kx_dim, ky_dim, theta_dim, lambda_dim, species_dim, egrid_dim
    implicit none
    character(len=*), intent(in) :: fname
    integer :: ncid
    !> Level of compression if enabled via
    !> `GK_NETCDF_DEFAULT_COMPRESSION_ON` macro. Must be a value between
    !> 1 (faster, less compression) and 9 (slower, more compression), or
    !> 0 (off)
    integer, parameter :: default_compression_level = 1

#ifdef GK_NETCDF_DEFAULT_COMPRESSION_ON
    neasyf_default_compression = default_compression_level
#else
    neasyf_default_compression = compression_off
#endif

    !First create a file
    ncid = neasyf_open(trim(adjustl(fname))//".grids.nc", "w")
    call create_metadata(ncid, "GS2 grids")

    call neasyf_dim(ncid, trim(kx_dim), values = akx, &
         long_name = "Wavenumber perpendicular to flux surface", units = "1/rho_r")
    call neasyf_dim(ncid, trim(ky_dim), values = aky, &
         long_name = "Wavenumber in direction of grad alpha", units = "1/rho_r")
    call neasyf_dim(ncid, trim(theta_dim), values = theta, long_name = "Poloidal angle", &
         units = "rad")
    call neasyf_dim(ncid, trim(lambda_dim), values = al, &
         long_name = "Energy/magnetic moment", units = "1/(2 B_a)")
    call neasyf_dim(ncid, trim(species_dim), dim_size = nspec)
    call neasyf_dim(ncid, trim(egrid_dim), dim_size = negrid, &
         long_name="See 'energy' for energy grid values")

    call neasyf_write(ncid, 'theta0', theta0, dim_names = [kx_dim,ky_dim], &
         long_name = 'theta0')
    call neasyf_write(ncid, "energy", energy, dim_names = [egrid_dim, species_dim], &
         long_name="Energy grid values")
    call nc_geo(ncid)
    call nc_species(ncid)
    call neasyf_write(ncid, 'wl', wl, dim_names = [theta_dim, lambda_dim], &
         long_name = "Pitch angle integration weights")
    call neasyf_write(ncid, 'w', w, dim_names = [egrid_dim, species_dim], &
         long_name = "Energy grid integration weights")
    call neasyf_close(ncid)
  end subroutine write_grids_to_file
#else
  !Fall back to binary output when netcdf unavailable

  !> FIXME : Add documentation
  subroutine write_grids_to_file(fname)
    use kt_grids, only: aky, akx, theta0, kperp2
    use theta_grid, only: theta
    use theta_grid, only: bmag, bpol, gradpar, grho
    use theta_grid, only: cdrift, cvdrift, gbdrift
    use theta_grid, only: cdrift0, cvdrift0, gbdrift0
    use theta_grid, only: gds2, gds21, gds22, gds23, gds24, gds24_noq
    use theta_grid, only: rplot, zplot, aplot
    use theta_grid, only: rprime, zprime, aprime
    use le_grids, only: al, energy, wl, w
    use file_utils, only: get_unused_unit

    implicit none
    character(len=*), intent(in) :: fname
    integer :: unit

    call get_unused_unit(unit)

    open(unit=unit,file=trim(adjustl(fname))//".ky",form="unformatted")
    write(unit) aky
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".kx",form="unformatted")
    write(unit) akx
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".theta0",form="unformatted")
    write(unit) theta0
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".theta",form="unformatted")
    write(unit) theta
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".lambda",form="unformatted")
    write(unit) al
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".energy",form="unformatted")
    write(unit) energy
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".bmag",form="unformatted")
    write(unit) bmag
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".bpol",form="unformatted")
    write(unit) bpol
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gradpar",form="unformatted")
    write(unit) gradpar
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".grho",form="unformatted")
    write(unit) grho
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".cdrift",form="unformatted")
    write(unit) cdrift
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".cvdrift",form="unformatted")
    write(unit) cvdrift
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gbdrift",form="unformatted")
    write(unit) gbdrift
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".cdrift0",form="unformatted")
    write(unit) cdrift0
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".cvdrift0",form="unformatted")
    write(unit) cvdrift0
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gbdrift0",form="unformatted")
    write(unit) gbdrift0
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds2",form="unformatted")
    write(unit) gds2
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds21",form="unformatted")
    write(unit) gds21
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds22",form="unformatted")
    write(unit) gds22
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds23",form="unformatted")
    write(unit) gds23
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds24",form="unformatted")
    write(unit) gds24
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gds24_noq",form="unformatted")
    write(unit) gds24_noq
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".rplot",form="unformatted")
    write(unit) rplot
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".zplot",form="unformatted")
    write(unit) zplot
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".aplot",form="unformatted")
    write(unit) aplot
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".rprime",form="unformatted")
    write(unit) rprime
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".zprime",form="unformatted")
    write(unit) zprime
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".aprime",form="unformatted")
    write(unit) aprime
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".kperp2",form="unformatted")
    write(unit) kperp2
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".wl",form="unformatted")
    write(unit) wl
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".w",form="unformatted")
    write(unit) w
    close(unit)
  end subroutine write_grids_to_file
#endif
end program dump_grids
