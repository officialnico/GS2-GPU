module diagnostics_zonal_transfer

!!!!! calculate tau(kx,ky), transfer of free energy by non-linearity in fourier space


  implicit none
  private

!!!!!!! Allocate arrays
  public :: init_diagnostics_transfer
!!!!!!!! Deallocate arrays
  public :: finish_diagnostics_transfer

  public :: calculate_zonal_transfer, write_zonal_transfer, total_term

  complex, dimension (:,:), allocatable ::  zonal_transfer !2 dimensions, kx, ky
	
	real, dimension(:,:), allocatable :: banew,gbnew,bracketnew  !!!!! dummy arrays

	complex, dimension (:,:,:,:), allocatable :: inter_s
	complex, dimension (:,:,:), allocatable :: g0
  !!! dummy arrays as well
contains

  subroutine init_diagnostics_transfer
    use kt_grids, only: naky, ntheta0, nx, ny, akx, aky
    use theta_grid, only: ntgrid
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    use le_grids, only: init_le_grids, nlambda, negrid
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: init_dist_fn_layouts
    use gs2_layouts, only: yxf_lo, g_lo
    
    implicit none
    logical :: accel
    
    accel=.false.
    
!!!calls every initiation routine we need: le_grids for integration, init_gs2 transforms to initiate fourier transforms routines...
    
    call init_le_grids(accel,accel)  !!!! benefit: calls init_theta_grid, init_kt_grids, init_gs2_layouts and init_species for later use
    call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accel)  !!! in order to use fft transform routines
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
!!!! accel_x =false : useless
    
!!!!! allocate and initialize arrays
    allocate (zonal_transfer(ntheta0,naky))
    zonal_transfer=0.0
    
    allocate(inter_s(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate(g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(banew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
    allocate(gbnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
    allocate(bracketnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
    banew = 0.
    gbnew = 0.
    bracketnew = 0.
    
  end subroutine init_diagnostics_transfer
  
  subroutine finish_diagnostics_transfer 
    deallocate(zonal_transfer)
    deallocate(banew)
    deallocate(gbnew)
    deallocate(bracketnew)
    deallocate(inter_s)
    deallocate(g0)
  end subroutine finish_diagnostics_transfer
  
!  subroutine non_linear_term (g1, phi, apar, bpar)   !no perturbed magnetic field: apar and bpar set to zero
  subroutine non_linear_term (g1)   !no perturbed magnetic field: apar and bpar set to zero
    use mp, only: max_allreduce, iproc
    use theta_grid, only: ntgrid, kxfac, bmag !!! get magnetic field as well
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx, il_idx
    use gs2_layouts, only: yxf_lo
    use dist_fn_arrays, only: g, g_adjust
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, reset, immediate_reset  !!! for us, fphi=1 (electrostatic case), fbar=fapar=0
    use kt_grids, only: aky, akx
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use le_grids, only: negrid, nlambda
    use dist_fn_arrays, only: g
    use fields_arrays, only: phi, apar, bpar
    use constants, only: zi   !i complex
    use species, only: nspec,spec !! need to get temperature for each species
    
    implicit none
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent(out) :: g1  !contains the distribution function g (not h), indices in velocity space as well
    complex, dimension(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g2,g3    !!!! for dummy arrays, specify exact shape
!    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar  !contains the potential phi, indices for theta, kx and ky
    integer :: i, j, k, isgn
    complex, dimension(-ntgrid:ntgrid,nspec,negrid,nlambda,2) :: test_sum
    integer, parameter :: verb = 4
    real :: fphi0
    integer :: iglo, ik, it, ig, ie, il, is   !indices for velocity, mode, time, species...
    
    fphi0=1.
    
    !Form g1=i*kx*chi     
    call load_kx_phi
    
    !Transform to real space
    call transform2 (g1, banew) !!!!the fft2d of g1 goes to array ba, which needs to be specified somewhere in the beginning! (init function)
    
    !Form g1=i*ky*h
    g1=g   !!!! g is an originally empty array from dist_fn_array, assigned with g somewhere
    call g_adjust(g1,phi,bpar,fphi0,fbpar)   !!! transforms g into h (minus boltzmann response), anyway fbpar=0
    
    !g2=g1 !! assigns values of hk in g2 for further use (see end of routine)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*aky(ik) !!!!for all indices in velocity space, multiply by i*ky 
    end do
    
    !Transform to real space    
    
    call transform2 (g1, gbnew) 
    !Calculate (d Chi /dx).(d h/dy)

    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
      do i = 1, yxf_lo%ny
        bracketnew(i,j) = banew(i,j)*gbnew(i,j)*kxfac   !!!! don't know why, but kxfac is important
      end do
    end do

    !Form g1=i*ky*chi
    call load_ky_phi

    !Transform to real space    

    call transform2 (g1, banew)

    !Form g1=i*kx*h
    g1=g

    call g_adjust(g1,phi,bpar,fphi0,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*akx(it)
    enddo
       
    !Transform to real space

    call transform2 (g1, gbnew)

    !Calculate (d Chi /dy).(d h/dx) and subtract from (d Chi /dx).(d h/dy)

    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
      do i = 1, yxf_lo%ny
        bracketnew(i,j) = bracketnew(i,j) - banew(i,j)*gbnew(i,j)*kxfac
      end do
    end do

    !Transform NL source back to spectral space

    call inverse2 (bracketnew, g1)  ! fourier non linear coeff Nk in g1

    !! TEST: is the poisson bracket correct? check from non_linear 
    !write(*,*) g1(:,1,g_lo%ulim_proc-10),'diag'
    
    g2=g
    call g_adjust(g2,phi,bpar,fphi0,fbpar)
    
    !! we have to multiply by hk                      
    do ig=-ntgrid,ntgrid
       do isgn=1,2
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             g3(ig,isgn,iglo)=g1(ig,isgn,iglo)*conjg(g2(ig,isgn,iglo))
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!! NEED TO MULTIPLY BY CONJUGATE ????????
!!!!! the sum over kx,ky of g1 should be 0! we calculate it
  
  !write(*,*) g3(10,1,g_lo%ulim_proc-10)
  
!   test_sum=0.
!   do iglo = g_lo%llim_proc, g_lo%ulim_proc
!      it = it_idx(g_lo,iglo)
!      ik = ik_idx(g_lo,iglo)
!      is = is_idx(g_lo,iglo)
!      ie = ie_idx(g_lo,iglo)
!      il = il_idx(g_lo,iglo)
!      do ig = -ntgrid,ntgrid
!         do isgn=1,2
!            test_sum(ig,is,ie,il,isgn) = test_sum(ig,is,ie,il,isgn) + g3(ig,isgn,iglo)
!         end do
!      end do
!   end do
!		write(*,*) abs(test_sum(1,1,2,20,1))/maxval(abs(g3))   !!!! should be zero

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             g1(ig,isgn,iglo)=spec(is)%temp*g3(ig,isgn,iglo)
          end do
       end do
    end do
    
  contains
    
    subroutine load_kx_phi
      
      use dist_fn_arrays, only: aj0  !!!bessel function from the fourier coeff of average
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc  !!!lower limit to upper limit?
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)  !!!ig for theta, it and ik for kx and ky
         do ig = -ntgrid, ntgrid
            fac = zi*akx(it)*aj0(ig,iglo)*phi(ig,it,ik)
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_phi

    subroutine load_ky_phi

      use dist_fn_arrays, only: aj0
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

!!!!!!!!!!!!!! integration over velocity space
  
  subroutine integrate (g, total)  !!!!! integrates g with maxwellian weights, output in total array
    
    use le_grids, only: integrate_moment  !!!! Not public, until told so in le_grids!
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use species, only:nspec
    
    implicit none
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    
    call integrate_moment(g, total)  !!!! can only be called once w and wl are assigned, after init_le_grids
    
  end subroutine integrate
  

!!!!!!! average over theta

  subroutine average_theta(a, axb)
    use theta_grid, only: ntgrid, delthet, jacob
    implicit none
    complex, dimension (-ntgrid:), intent (in) :: a
    complex, intent (out) :: axb
    integer :: ng
    real, dimension (-ntgrid:ntgrid) :: wgt     !!!! in volume_averages, ntg_out=ntgrid... suppose likewise
    real :: anorm

    ng = ntgrid
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    axb = sum(a(-ng:ng)*wgt)/anorm
  end subroutine average_theta

!  subroutine total_term(g,phi,tau)   !!!!transfer of energy in terms of kx, ky, stored in tau, sums up previous routines
  subroutine total_term(tau)   !!!!transfer of energy in terms of kx, ky, stored in tau, sums up previous routines
    
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use species, only: nspec,spec
    use kt_grids, only: ntheta0, naky
    
    implicit none
    
!    complex, dimension (-ntgrid:,:,g_lo%llim_proc:) :: g
!    complex, dimension(-ntgrid:ntgrid,ntheta0,naky) :: phi, apar, bpar   !!!! for dummy arrays, specify exact shape
    complex, dimension(:,:), intent(out) :: tau   !!!!! intent statement should not be taken lightly...may raise exceptions
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky) :: inter
    integer :: ik, it, is, iglo
    
    !!calculation of non linear term
 !   apar=0.
 !   bpar=0.
    
!    call non_linear_term(g,phi,apar,bpar)
    call non_linear_term(g0)
    
!!!integrate over velocities
    inter=0.
    inter_s=0.
!    g0=g
    call integrate(g0,inter_s)
    !write(*,*) inter_s
    
!!! sum over species
    do is=1,nspec
       inter(:,:,:)=inter(:,:,:) + inter_s(:,:,:,is)
    end do
    
    !! integrate over theta, must be done after the rest (dimension of arrays)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       !write(*,*) inter(-ntgrid:,it,ik)
       call average_theta(inter(-ntgrid:,it,ik),tau(it,ik))
    end do
    
  end subroutine total_term
  
  
  subroutine write_zonal_transfer(gnostics)
    
    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use diagnostics_dimensions, only: dim_string
    use fields_parallelization, only: field_k_local
    use kt_grids, only: ntheta0, naky
    use diagnostics_config, only: diagnostics_type
    
    implicit none
    
    type(diagnostics_type), intent(in) :: gnostics
    
    call create_and_write_distributed_fieldlike_variable(gnostics, gnostics%rtype, "zonal_transfer", &
         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         "Non linear transfer of free energy from for each mode, as a function of kx, ky and time", &
         "Tr2*rhor*c*ns*vth6 / Ba2*e*a3", zonal_transfer)
    
  end subroutine write_zonal_transfer
  
  subroutine calculate_zonal_transfer(gnostics)

!   use fields_arrays, only: phi
!   use dist_fn_arrays, only: gnew                   !!!!!!!!!!!!! HAVE TO CHOOSE PHI AND GNEW, NOT ANYTHING ELSE
   use diagnostics_config, only: diagnostics_type
   
   implicit none
   
   type(diagnostics_type), intent(inout) :: gnostics
   
!   call total_term(gnew,phi,zonal_transfer)
   call total_term(zonal_transfer)
   
   gnostics%current_results%zonal_transfer = zonal_transfer !!!copy results to gnostics
    
 end subroutine calculate_zonal_transfer
 
end module diagnostics_zonal_transfer













