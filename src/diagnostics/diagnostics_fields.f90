!> This module contains functions for writing quantities
!! that have the dimensions of the fields, including
!! volume averaged quantities, etc.
module diagnostics_fields
  implicit none
  
  private
  
  public :: write_fields, write_movie, write_eigenfunc, get_phi0
  real, parameter :: zero = epsilon(0.0)

contains
  !> FIXME : Add documentation
  subroutine write_fields(gnostics)
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: nc_write_field
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    if (gnostics%replay) return 
    if (.not. gnostics%wryte) return

    if (fphi > zero) then
      call nc_write_field(gnostics%file_id, gnostics%nout, phinew, &
           "phi", "Electrostatic potential", gnostics%igomega, gnostics%write_phi_over_time,&
           gnostics%current_results%phi2)
    end if

    if (fapar > zero) then
      call nc_write_field(gnostics%file_id, gnostics%nout, aparnew, &
           "apar", "Parallel magnetic potential", gnostics%igomega, gnostics%write_apar_over_time,&
           gnostics%current_results%apar2)
    end if

    if (fbpar > zero) then
      call nc_write_field(gnostics%file_id, gnostics%nout, bparnew, &
           "bpar", "Parallel magnetic field", gnostics%igomega, gnostics%write_bpar_over_time,&
           gnostics%current_results%bpar2)
    end if
  end subroutine write_fields

  !> FIXME : Add documentation  
  subroutine write_movie(gnostics)
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: nc_loop_movie
    implicit none
    type(diagnostics_type), intent(in) :: gnostics

    if (proc0) then
      call nc_loop_movie(gnostics%file_id, gnostics%nout, gnostics%user_time)
    end if
  end subroutine write_movie

  !> Calculate the `phi0` array, used to normalise eigenfunctions and
  !> moments.
  function get_phi0() result(phi0)
    use fields_arrays, only: phi, apar
    use kt_grids, only: ntheta0, naky
    use dist_fn, only: def_parity, even
    use run_parameters, only: fapar
    use theta_grid, only: theta
    implicit none
    complex, dimension (ntheta0, naky) :: phi0

    !Should this actually use igomega instead of 0?
    !What if fphi==0? --> See where statements below
    phi0 = phi(0,:,:)

    !This looks like a hack for the case where we know we've forced phi(theta=0) to be 0
    !this could probably be better addressed by the use of igomega above
    if (def_parity .and. fapar > 0 .and. (.not. even)) phi0 = apar(0, :, :)

    !Address locations where phi0=0 by using next point
    where (abs(phi0) < 10.0*epsilon(0.0))
       phi0 = phi(1,:,:)/(theta(1)-theta(0))
    end where

    !Check again if any locations are 0, this could be true if fphi (or fapar)
    !is zero.
    where (abs(phi0) < 10.0*epsilon(0.0))
       phi0 = 1.0
    end where

  end function get_phi0

  !> Write out the fields, normalized to their value at theta=0
  subroutine write_eigenfunc(gnostics)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    use fields_arrays, only: phi, apar, bpar
    use kt_grids, only: ntheta0, naky, theta0, aky
    use theta_grid, only: theta, nperiod, ntheta, ntgrid
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: nc_eigenfunc
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (ntheta0, naky) :: phi0
    integer :: it, ik, ig, ntg_out
    
    if(.not.proc0 .or. gnostics%replay) return

    phi0 = get_phi0()
    
    !Do ascii output
    if (proc0 .and. gnostics%ascii_files%write_to_eigenfunc .and. (.not. gnostics%create)) then
       ntg_out = ntheta/2 + (nperiod-1)*ntheta
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out
                write (unit=gnostics%ascii_files%eigenfunc, fmt="(9(1x,e12.5))") &
                     theta(ig), theta0(it,ik), aky(ik), &
                     phi(ig,it,ik)/phi0(it,ik), &
                     apar(ig,it,ik)/phi0(it,ik), &
                     bpar(ig,it,ik)/phi0(it,ik)
             end do
             write (unit=gnostics%ascii_files%eigenfunc, fmt="()")
          end do
       end do
    end if
    
    call nc_eigenfunc(gnostics%file_id, phi0)
  end subroutine write_eigenfunc

end module diagnostics_fields

