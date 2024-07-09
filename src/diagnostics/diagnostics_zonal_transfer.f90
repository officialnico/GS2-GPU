!> FIXME : Add documentation
module diagnostics_zonal_transfer

!> calculate tau(kx,ky), the transfer of free energy by the non-linearity in fourier space
  implicit none

  private

!> allocate free energy transfer arrays
  public :: init_diagnostics_transfer
  public :: calculate_zonal_transfer, write_zonal_transfer

contains
  !> FIXME : Add documentation
  subroutine init_diagnostics_transfer(gnostics)
    use kt_grids, only: naky, ntheta0, nx, ny
    use theta_grid, only: ntgrid
    use species, only: nspec
    use le_grids, only: init_le_grids, nlambda, negrid
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: init_dist_fn_layouts
    use diagnostics_config, only: diagnostics_type
    use mp, only: nproc, iproc
    implicit none

    type(diagnostics_type), intent(inout) :: gnostics
    logical :: accel_x

    if(.not. gnostics%write_zonal_transfer) return

    !> initialize zonal_transfer to zero
    gnostics%current_results%zonal_transfer = 0.0
    
    ! call every initiation routine we need:
    ! le_grids for integration, init_gs2 transforms to initiate fourier transforms routines, etc.
    
    ! calls init_theta_grid, init_kt_grids, init_gs2_layouts, and init_species for later use
    call init_le_grids
    ! will need to use fft routines
    call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accel_x)
    call init_dist_fn_layouts(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nproc, iproc)

  end subroutine init_diagnostics_transfer
  
  !> calculate nonlinear term appearing in GKE
  !! for now assume electrostatic
  subroutine non_linear_term (g1)
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: yxf_lo
    use dist_fn_arrays, only: g, g_adjust, to_g_gs2, from_g_gs2
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use kt_grids, only: aky, akx
    use fields_arrays, only: phi, bpar
    use constants, only: zi

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent(out) :: g1
    integer :: i, j
    integer :: iglo, ik, it, ig, is

    real, dimension(:,:), allocatable :: banew, gbnew, bracketnew

    allocate (banew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; banew = 0.
    allocate (gbnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; gbnew = 0.
    allocate (bracketnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; bracketnew = 0.

    ! Form g1=i*kx*phi
    call load_kx_phi

    ! Transform to real space
    ! 2D fft of g1 is returned as banew
    call transform2 (g1, banew)

    ! Form g1=i*ky*g
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g1(:,:,iglo)=g(:,:,iglo)*zi*aky(ik)
    end do

    ! Transform to real space
    call transform2 (g1, gbnew)

    ! Calculate (d phi /dx).(dg/dy)
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracketnew(i,j) = banew(i,j)*gbnew(i,j)*kxfac
       end do
    end do

    ! Form g1=i*ky*chi
    call load_ky_phi

    !Transform to real space
    call transform2 (g1, banew)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       g1(:,:,iglo)=g(:,:,iglo)*zi*akx(it)
    enddo

    ! Transform to real space
    call transform2 (g1, gbnew)

    ! Calculate (d phi /dy).(dg/dx) and subtract from (d phi /dx).(dg/dy)
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracketnew(i,j) = bracketnew(i,j) - banew(i,j)*gbnew(i,j)*kxfac
       end do
    end do

    ! Transform nonlinearity back to spectral space
    ! g1 contains Fourier coefficients associated with nonlinearity
    call inverse2 (bracketnew, g1)

    call g_adjust(g,phi,bpar, direction = from_g_gs2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*conjg(g(:,:,iglo))*spec(is)%temp
    end do

    call g_adjust(g,phi,bpar, direction = to_g_gs2)

    deallocate (banew, gbnew, bracketnew)

  contains

    !> FIXME : Add documentation
    subroutine load_kx_phi
      use dist_fn_arrays, only: aj0  !< bessel function from the fourier coeff of average
      implicit none
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*akx(it)*aj0(ig,iglo)*phi(ig,it,ik)
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_phi

    !> FIXME : Add documentation    
    subroutine load_ky_phi
      use dist_fn_arrays, only: aj0
      implicit none
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*aky(ik)*aj0(ig,iglo)*phi(ig,it,ik)
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_ky_phi

  end subroutine non_linear_term

  !> Returns tau, the transfer of free energy as a function of (kx,ky)
  subroutine total_term(tau)
    use mp, only: proc0, broadcast
    use theta_grid, only: ntgrid, field_line_average
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use species, only: nspec
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment

    implicit none
    
    complex, dimension(:,:), intent(out) :: tau

    complex, dimension (:,:,:), allocatable :: g0
    complex, dimension (:,:,:), allocatable :: inter
    complex, dimension (:,:,:,:), allocatable :: inter_s

    integer :: ik, it, is
    
    allocate (g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (inter(-ntgrid:ntgrid,ntheta0,naky)) ; inter = 0.
    allocate (inter_s(-ntgrid:ntgrid,ntheta0,naky,nspec)) ; inter_s = 0.

    ! obtain the nonlinear term
    call non_linear_term(g0)

    ! integrate over velocities
    ! can only be called after w and wl are assigned in init_le_grids
    call integrate_moment(g0, inter_s)

    if (proc0) then
       ! sum over species
       do is=1,nspec
          inter = inter + inter_s(:,:,:,is)
       end do

       ! integrate over theta, must be done after velocity space integration
       do ik = 1, naky
          do it = 1, ntheta0
             tau(it, ik) = field_line_average(inter(:,it,ik))
          end do
       end do
    end if
    
    do ik = 1, naky
       call broadcast (tau(:,ik))
    end do

    deallocate (g0, inter, inter_s)
    
  end subroutine total_term

  !> FIXME : Add documentation  
  subroutine write_zonal_transfer(gnostics)
    use gs2_io, only: starts, netcdf_write_complex
    use fields_parallelization, only: field_k_local
    use diagnostics_config, only: diagnostics_type
    use mp, only: proc0
    implicit none
    
    type(diagnostics_type), intent(in out) :: gnostics
    if (.not. proc0) return
    call netcdf_write_complex(gnostics%file_id, "zonal_transfer", gnostics%current_results%zonal_transfer, &
         dim_names=["ri" ,"kx", "ky", "t "], start=starts(4, gnostics%nout), &
         long_name="Time rate of change of free energy due to nonlinear transfer, as a function of kx, ky and time", &
         units="Tr2*rhor*c*ns*vth6 / Ba2*e*a3")
    
  end subroutine write_zonal_transfer

  !> FIXME : Add documentation  
  subroutine calculate_zonal_transfer(gnostics)
    use diagnostics_config, only: diagnostics_type
    use kt_grids, only: ntheta0, naky
    implicit none
    
    type(diagnostics_type), intent(inout) :: gnostics
    
    complex, dimension (:,:), allocatable ::  zonal_transfer

    allocate (zonal_transfer(ntheta0,naky)) ; zonal_transfer = 0.

    call total_term(zonal_transfer)
    gnostics%current_results%zonal_transfer = zonal_transfer

    deallocate (zonal_transfer)

  end subroutine calculate_zonal_transfer

end module diagnostics_zonal_transfer
