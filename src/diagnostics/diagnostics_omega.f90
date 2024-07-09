!> A module which calculates and writes the growth rates and frequencies
module diagnostics_omega
 
  implicit none

  private

  public :: write_omega
  public :: init_diagnostics_omega
  public :: finish_diagnostics_omega
  public :: calculate_omega
  public :: omega_average
  public :: omegahist
  public :: debug

  complex, dimension (:,:), allocatable :: omega_average
  logical :: debug =.false.

  complex, dimension (:,:,:), allocatable :: omegahist
  complex, allocatable, save, dimension (:,:,:) :: domega


contains 

  !> Allocate arrays for storing omega history and averages
  subroutine init_diagnostics_omega(gnostics)
    use kt_grids, only: ntheta0, naky
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    allocate (omegahist(0:gnostics%navg-1,ntheta0,naky))
    allocate(domega(gnostics%navg,ntheta0,naky))
    allocate(omega_average(ntheta0, naky))
    omegahist = 0.0
    omega_average = 0.0
  end subroutine init_diagnostics_omega

  !> FIXME : Add documentation  
  subroutine finish_diagnostics_omega
    implicit none
    deallocate(omegahist)
    deallocate(domega)
    deallocate(omega_average)
  end subroutine finish_diagnostics_omega
  
  !> Write omega, omega_average as functions of kx, ky and time
  !! to the new netcdf file
  subroutine write_omega(gnostics)
    use gs2_time, only: woutunits
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: nc_write_omega
    implicit none
    type(diagnostics_type), intent(in) :: gnostics

    if (.not. gnostics%writing) return

    call nc_write_omega(gnostics%file_id, gnostics%nout, &
         omegahist(mod(gnostics%istep, gnostics%navg), :, :), omega_average, woutunits)

  end subroutine write_omega

  !> Calculates omega and stores it in omegahist. Calculates omega_average, 
  !! the average of omega over navg timesteps
  subroutine calculate_omega (gnostics)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_time, only: code_dt
    use constants, only: zi
    use diagnostics_config, only: diagnostics_type
    use run_parameters, only: ieqzip
    use nonlinear_terms, only: nonlin
    use mp, only: proc0
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    real :: fac
    integer :: j

    ! Avoid calculation of omega in nonlinear calculations.
    ! Note we rely on omega_average and related arrays being
    ! initialised elsewhere.
    if (nonlin .and. .not. any(ieqzip)) then
       return
    end if

    if (debug) write(6,*) "get_omeaavg: start"

    j = gnostics%igomega
    if(gnostics%omegatol.eq.0)then
       fac=1.0
    else
       fac=1000/abs(gnostics%omegatol)
    endif
    
    !<DD>The logic below was originally
    !    where (abs(phinew(j,:,:)+aparnew(j,:,:)+bparnew(j,:,:)) < epsilon(0.0) &
    !           .or. abs(phi(j,:,:)+apar(j,:,:)+bpar(j,:,:)) < epsilon(0.0))
    !This results in not calculating the frequency whenever the fields drop below epsilon(0.0) [~1d-16 for DP]
    !What we really want to do is prevent dividing by too small a number. 
    ! Calculate omega
    where (abs(phinew(j,:,:)+aparnew(j,:,:)+bparnew(j,:,:)) < tiny(0.0)*fac &
         .or. abs(phi(j,:,:)+apar(j,:,:)+bpar(j,:,:)) < tiny(0.0)*fac)
       omegahist(mod(gnostics%istep,gnostics%navg),:,:) = 0.0
    elsewhere
       omegahist(mod(gnostics%istep,gnostics%navg),:,:) &
            = log((phinew(j,:,:) + aparnew(j,:,:) + bparnew(j,:,:)) &
            /(phi(j,:,:)   + apar(j,:,:)    + bpar(j,:,:)))*zi/code_dt
    end where
    
    ! Calculate omega_average
    omega_average = sum(omegahist/real(gnostics%navg),dim=1)
    ! Copy the results to gnostics
    gnostics%current_results%omega_average = omega_average
    
    if (debug) write(6,*) "calculate_omega: omega_average=",omega_average
    
    ! Check if converged by comparing the difference between
    ! omega_average and all values in omegahist, to the value of omega_average
    if (gnostics%istep > gnostics%navg) then
       domega = spread(omega_average,1,gnostics%navg) - omegahist
       if (all(sqrt(sum(abs(domega)**2/real(gnostics%navg),dim=1)) &
            .le. min(abs(omega_average),1.0)*gnostics%omegatol)) &
       then
          if (proc0 .and. gnostics%ascii_files%write_to_out) write (gnostics%ascii_files%out, "('*** omega converged')")
          gnostics%exit = gnostics%exit_when_converged
       end if
       
       if (any(abs(omega_average)*code_dt > gnostics%omegatinst)) then
          !In old diagnostics this would only be written by proc0 and it would
          !go to gnostics%ascii_files%out rather than screen. Should we make this
          !consistent here?
          write (*, "('*** numerical instability detected')") 
          gnostics%exit = .true.
       end if
    end if
    
    if (debug) write(6,*) "calculate_omega: done"
  end subroutine calculate_omega
end module diagnostics_omega
