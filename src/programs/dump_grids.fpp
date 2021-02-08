!> FIXME : Add documentation
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
  logical :: accel_x, accel_y
  logical, parameter :: quiet=.false.
  character (500), target :: cbuff
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
     call time_message(.false.,time_measure,'init-theta')
  endif
  call init_theta_grid
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-theta')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/KT
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init kt_grids   ")',advance='no')
     call time_message(.false.,time_measure,'init-kt')
  endif
  call init_kt_grids
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-kt')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/SPECIES
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init species   ")',advance='no')
     call time_message(.false.,time_measure,'init-species')
  endif
  call init_species
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-species')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/LE
  if(proc0.and.(.not.quiet)) then
     write(*,'("Init le_grids   ")',advance='no')
     call time_message(.false.,time_measure,'init-le')
  endif
  call init_le_grids(accel_x,accel_y)
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-le')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

  !Now write to file
  if(proc0.and.(.not.quiet)) then
     write(*,'("Write to file   ")',advance='no')
     call time_message(.false.,time_measure,'write-file')
  endif
  filename=trim(adjustl(run_name))
  if(proc0) call write_grids_to_file(filename)
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'write-file')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

contains

#ifdef NETCDF
!Prefer to write to netcdf file if possible

  !> FIXME : Add documentation
  subroutine write_grids_to_file(fname)
    use netcdf, only: NF90_CLOBBER
    use netcdf, only: nf90_create, nf90_close
    use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
    use netcdf, only: nf90_put_var

    use netcdf_utils, only: get_netcdf_code_precision
    use netcdf_utils, only: netcdf_real

    use kt_grids, only: aky, akx, theta0
    use theta_grid, only: theta
    use theta_grid, only: bmag, bpol, gradpar, grho
    use theta_grid, only: cdrift, cvdrift, gbdrift
    use theta_grid, only: cdrift0, cvdrift0, gbdrift0
    use theta_grid, only: cvdrift_th, gbdrift_th
    use theta_grid, only: gds2, gds21, gds22, gds23, gds24, gds24_noq
    use theta_grid, only: rplot, zplot, aplot
    use theta_grid, only: rprime, zprime, aprime
    use le_grids, only: al, energy, negrid, wl, w
    use species, only: nspec

    implicit none
    character(len=*), intent(in) :: fname
    integer :: ncid, ierr
    integer :: ky_dimid, kx_dimid, theta_dimid, al_dimid, energy_dimid, species_dimid
    integer :: ky_varid, kx_varid, theta_varid, al_varid, energy_varid, theta0_varid
    integer :: bmag_varid, bpol_varid, gradpar_varid, grho_varid
    integer :: cdrift_varid, cvdrift_varid, gbdrift_varid
    integer :: cdrift0_varid, cvdrift0_varid, gbdrift0_varid
    integer :: cvdrift_th_varid, gbdrift_th_varid
    integer :: gds2_varid, gds21_varid, gds22_varid, gds23_varid, gds24_varid, gds24_noq_varid
    integer :: rplot_varid, zplot_varid, aplot_varid, rprime_varid, zprime_varid, aprime_varid
    integer :: wl_varid, w_varid

    !First get precisions
    if(netcdf_real==0) netcdf_real=get_netcdf_code_precision()

    !First create a file
    ierr=nf90_create(trim(adjustl(fname))//".grids.nc",NF90_CLOBBER,ncid)

    !Define dimensions
    ierr=nf90_def_dim(ncid,'aky',size(aky),ky_dimid)
    ierr=nf90_def_dim(ncid,'akx',size(akx),kx_dimid)
    ierr=nf90_def_dim(ncid,'theta',size(theta),theta_dimid)
    ierr=nf90_def_dim(ncid,'lambda',size(al),al_dimid)
    ierr=nf90_def_dim(ncid,'energy',negrid,energy_dimid)
    ierr=nf90_def_dim(ncid,'species',nspec,species_dimid)

    !Define basic variables
    ierr=nf90_def_var(ncid,'ky',netcdf_real,[ky_dimid],ky_varid)
    ierr=nf90_def_var(ncid,'kx',netcdf_real,[kx_dimid],kx_varid)
    ierr=nf90_def_var(ncid,'theta0',netcdf_real,[kx_dimid,ky_dimid],theta0_varid)
    ierr=nf90_def_var(ncid,'theta',netcdf_real,[theta_dimid],theta_varid)
    ierr=nf90_def_var(ncid,'lambda',netcdf_real,[al_dimid],al_varid)
    ierr=nf90_def_var(ncid,'energy',netcdf_real,[energy_dimid, species_dimid],energy_varid)

    !Define associated variables
    ierr=nf90_def_var(ncid,'bmag',netcdf_real,[theta_dimid], bmag_varid)
    ierr=nf90_def_var(ncid,'bpol',netcdf_real,[theta_dimid], bpol_varid)
    ierr=nf90_def_var(ncid,'gradpar',netcdf_real,[theta_dimid], gradpar_varid)
    ierr=nf90_def_var(ncid,'grho',netcdf_real,[theta_dimid], grho_varid)
    ierr=nf90_def_var(ncid,'cdrift',netcdf_real,[theta_dimid], cdrift_varid)
    ierr=nf90_def_var(ncid,'cvdrift',netcdf_real,[theta_dimid], cvdrift_varid)
    ierr=nf90_def_var(ncid,'gbdrift',netcdf_real,[theta_dimid], gbdrift_varid)
    ierr=nf90_def_var(ncid,'cdrift0',netcdf_real,[theta_dimid], cdrift0_varid)
    ierr=nf90_def_var(ncid,'cvdrift0',netcdf_real,[theta_dimid], cvdrift0_varid)
    ierr=nf90_def_var(ncid,'gbdrift0',netcdf_real,[theta_dimid], gbdrift0_varid)
    ierr=nf90_def_var(ncid,'cvdrift_th',netcdf_real,[theta_dimid], cvdrift_th_varid)
    ierr=nf90_def_var(ncid,'gbdrift_th',netcdf_real,[theta_dimid], gbdrift_th_varid)
    ierr=nf90_def_var(ncid,'gds2',netcdf_real,[theta_dimid], gds2_varid)
    ierr=nf90_def_var(ncid,'gds21',netcdf_real,[theta_dimid], gds21_varid)
    ierr=nf90_def_var(ncid,'gds22',netcdf_real,[theta_dimid], gds22_varid)
    ierr=nf90_def_var(ncid,'gds23',netcdf_real,[theta_dimid], gds23_varid)
    ierr=nf90_def_var(ncid,'gds24',netcdf_real,[theta_dimid], gds24_varid)
    ierr=nf90_def_var(ncid,'gds24_noq',netcdf_real,[theta_dimid], gds24_noq_varid)
    ierr=nf90_def_var(ncid,'rplot',netcdf_real,[theta_dimid], rplot_varid)
    ierr=nf90_def_var(ncid,'zplot',netcdf_real,[theta_dimid], zplot_varid)
    ierr=nf90_def_var(ncid,'aplot',netcdf_real,[theta_dimid], aplot_varid)
    ierr=nf90_def_var(ncid,'rprime',netcdf_real,[theta_dimid], rprime_varid)
    ierr=nf90_def_var(ncid,'zprime',netcdf_real,[theta_dimid], zprime_varid)
    ierr=nf90_def_var(ncid,'aprime',netcdf_real,[theta_dimid], aprime_varid)
    ierr=nf90_def_var(ncid,'wl',netcdf_real,[theta_dimid, al_dimid],wl_varid)
    ierr=nf90_def_var(ncid,'w',netcdf_real,[species_dimid, energy_dimid],w_varid)

    !End definitions
    ierr=nf90_enddef(ncid)

    !Now fill in variable data
    ierr=nf90_put_var(ncid,ky_varid,aky)
    ierr=nf90_put_var(ncid,kx_varid,akx)
    ierr=nf90_put_var(ncid,theta0_varid,theta0)
    ierr=nf90_put_var(ncid,theta_varid,theta)
    ierr=nf90_put_var(ncid,al_varid,al)
    ierr=nf90_put_var(ncid,energy_varid,energy)

    ierr=nf90_put_var(ncid, bmag_varid, bmag)
    ierr=nf90_put_var(ncid, bpol_varid, bpol)
    ierr=nf90_put_var(ncid, gradpar_varid, gradpar)
    ierr=nf90_put_var(ncid, grho_varid, grho)
    ierr=nf90_put_var(ncid, cdrift_varid, cdrift)
    ierr=nf90_put_var(ncid, cvdrift_varid, cvdrift)
    ierr=nf90_put_var(ncid, gbdrift_varid, gbdrift)
    ierr=nf90_put_var(ncid, cdrift0_varid, cdrift0)
    ierr=nf90_put_var(ncid, cvdrift0_varid, cvdrift0)
    ierr=nf90_put_var(ncid, gbdrift0_varid, gbdrift0)
    ierr=nf90_put_var(ncid, cvdrift_th_varid, cvdrift_th)
    ierr=nf90_put_var(ncid, gbdrift_th_varid, gbdrift_th)
    ierr=nf90_put_var(ncid, gds2_varid, gds2)
    ierr=nf90_put_var(ncid, gds21_varid, gds21)
    ierr=nf90_put_var(ncid, gds22_varid, gds22)
    ierr=nf90_put_var(ncid, gds23_varid, gds23)
    ierr=nf90_put_var(ncid, gds24_varid, gds24)
    ierr=nf90_put_var(ncid, gds24_noq_varid, gds24_noq)
    ierr=nf90_put_var(ncid, rplot_varid, rplot)
    ierr=nf90_put_var(ncid, zplot_varid, zplot)
    ierr=nf90_put_var(ncid, aplot_varid, aplot)
    ierr=nf90_put_var(ncid, rprime_varid, rprime)
    ierr=nf90_put_var(ncid, zprime_varid, zprime)
    ierr=nf90_put_var(ncid, aprime_varid, aprime)
    ierr=nf90_put_var(ncid,wl_varid,wl)
    ierr=nf90_put_var(ncid,w_varid,w)

    !Now close file
    ierr=nf90_close(ncid)
  end subroutine write_grids_to_file
#else
  !Fall back to binary output when netcdf unavailable

  !> FIXME : Add documentation
  subroutine write_grids_to_file(fname)
    use kt_grids, only: aky, akx, theta0
    use theta_grid, only: theta
    use theta_grid, only: bmag, bpol, gradpar, grho
    use theta_grid, only: cdrift, cvdrift, gbdrift
    use theta_grid, only: cdrift0, cvdrift0, gbdrift0
    use theta_grid, only: cvdrift_th, gbdrift_th
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

    open(unit=unit,file=trim(adjustl(fname))//".cvdrift_th",form="unformatted")
    write(unit) cvdrift_th
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".gbdrift_th",form="unformatted")
    write(unit) gbdrift_th
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

    open(unit=unit,file=trim(adjustl(fname))//".wl",form="unformatted")
    write(unit) wl
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".w",form="unformatted")
    write(unit) w
    close(unit)
  end subroutine write_grids_to_file
#endif
end program dump_grids
