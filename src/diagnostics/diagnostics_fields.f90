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
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_transforms, only: transform2
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0, mp_abort
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use kt_grids, only: nx, ny, naky, ntheta0
    use file_utils, only: error_unit
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t
    integer :: nnx, nny
    real, dimension(:,:,:), allocatable :: yxphi, yxapar, yxbpar
    logical :: accelerated ! dummy variable here
    logical, save :: warning = .false.

    call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

    t = gnostics%user_time
    
    ! EAB 09/17/03 -- modify dump_fields_periodically to print out inverse fft of fields in x,y
    nnx = yxf_lo%nx
    nny = yxf_lo%ny
    
    if (proc0 .and. .not. warning) then
       write (error_unit(),*) "WARNING: write_movie writes arrays of size nx*ny which contain a lot of &
            & redundancy, leading to unnecessarily large output files... consider doing &
            & the Fourier transforms of the fields in post-processing instead"
       warning = .true.
    endif
    
    if (fphi > epsilon(0.0)) then
       allocate (yxphi(nnx,nny,-ntgrid:ntgrid))
       call transform2 (phinew, yxphi, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "phi_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The electric potential in real space", "Tr/e", yxphi)
       deallocate (yxphi)
    end if

    if (fapar > epsilon(0.0)) then
       allocate (yxapar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (aparnew, yxapar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "apar_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The parallel vector potential in real space", "TBC", yxapar)
       deallocate (yxapar)
    end if

    if (fbpar > epsilon(0.0)) then 
       allocate (yxbpar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (bparnew, yxbpar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "bpar_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The parallel magnetic field in real space", "TBC", yxbpar)
       deallocate (yxbpar)
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
    if (proc0 .and. gnostics%write_ascii) then
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

