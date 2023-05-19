!> A module for calculating properties of the turbulence
!! such as cross phases and correlations
module diagnostics_turbulence
  implicit none
  private

  public :: write_cross_phase
contains

  !> FIXME : Add documentation  
  subroutine write_cross_phase(gnostics)
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use neasyf, only: neasyf_write
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: phase_tot, phase_theta
    real :: t
    
    t = gnostics%user_time
    
    call get_cross_phase (gnostics, phase_tot, phase_theta)
    
    if (.not. gnostics%wryte) return

    call neasyf_write(gnostics%file_id, "electron_cross_phase_theta", phase_theta, &
         dim_names=["t"], start=[gnostics%nout], &
         long_name="Cross phase between electron density and perpendicular electron temperature, &
         &  given at the outboard midplane and averaged across x and y", units="radians")
    call neasyf_write(gnostics%file_id, "electron_cross_phase_tot", phase_tot, &
         dim_names=["t"], start=[gnostics%nout], &
         long_name="Cross phase between electron density and perpendicular electron temperature, &
         & averaged across all space", units="radians")
    
    if (proc0 .and. gnostics%ascii_files%write_to_phase .and. (.not. gnostics%create)) write (unit=gnostics%ascii_files%phase, &
         fmt="('t= ',e17.10,' phase_tot= ',e11.4,' phase_theta= ',e11.4)") &
         & t, phase_tot, phase_theta
  end subroutine write_cross_phase

  !> This is a highly simplified synthetic diagnostic which 
  !! calculates the cross phase between the electron density and the 
  !! perpendicular electron temperature for comparisons with DIII-D.  
  !! Returns the value of the cross-phase at the outboard midplane and 
  !! integrated over all v and x. We can generalize this routine to 
  !! other fields at some point, but for now this is just a skeleton for 
  !! a more realistic synthetic diagnostic. 
  subroutine get_cross_phase (gnostics, phase_tot, phase_theta)
    use species, only: nspec, spec, is_electron_species
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn, only: getemoms
    use fields_arrays, only: phinew, bparnew
    use volume_averages, only: average_all
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, intent (out) :: phase_tot, phase_theta
    complex, dimension (:,:,:,:), allocatable :: ntot, tperp
    complex :: nTp
    integer ::  is
    
    allocate ( ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))
    
    call getemoms (phinew, bparnew, ntot, tperp)
    
    do is = 1,nspec
       if (is_electron_species(spec(is))) then
          ! get cross_phase at outboard midplane
          call average_all(ntot(0,:,:,is), tperp(0,:,:,is), nTp, gnostics%distributed)
          phase_theta = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
          ! get integrated cross_phase 
          call average_all(ntot(:,:,:,is), tperp(:,:,:,is), nTp, gnostics%distributed)
          phase_tot = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
       end if
    end do
    
    deallocate (ntot, tperp)
    
  end subroutine get_cross_phase
end module diagnostics_turbulence
