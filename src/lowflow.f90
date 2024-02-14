!> FIXME : Add documentation
module lowflow

  implicit none

  private

  public :: get_lowflow_terms, finish_lowflow_terms, init_lowflow
  public :: phineo, dphidth, mach_lab, wstar_neo
  public :: hneoc, vparterm, wdfac, wstarfac, wdttpfac

  real, dimension (:), allocatable :: dphidth
  real, dimension (:,:,:,:,:), allocatable :: coefs
  real, dimension (:,:), allocatable :: phineo, dcoefsdr
  real, dimension (:,:,:), allocatable :: dcoefsdth
  real, dimension (:), allocatable :: rad_neo, theta_neo
  real, dimension (:,:), allocatable :: dens_neo, temp_neo
  real, dimension (:), allocatable :: mach_lab

  real, dimension (:,:), allocatable :: wcurv
  real, dimension (:,:,:), allocatable :: wstar_neo

  ! v-dependent factors in low-flow terms
  real, dimension (:,:,:), allocatable :: hneoc, vparterm, wdfac, wstarfac
  real, dimension (:,:,:,:,:,:), allocatable :: wdttpfac

contains
  !> FIXME : Add documentation
  subroutine finish_lowflow_terms
    implicit none
    !<DD>May actually be able to deallocate most of these much
    !    earlier as most only used during initialisation
    if(allocated(dphidth)) deallocate(dphidth)
    if(allocated(coefs)) deallocate(coefs)
    if(allocated(phineo)) deallocate(phineo)
    if(allocated(dcoefsdr)) deallocate(dcoefsdr, dcoefsdth)
    if(allocated(rad_neo)) deallocate(rad_neo)
    if(allocated(theta_neo)) deallocate(theta_neo)
    if(allocated(mach_lab)) deallocate(mach_lab)
    if(allocated(dens_neo)) deallocate(dens_neo)
    if(allocated(temp_neo)) deallocate(temp_neo)
    if(allocated(wcurv)) deallocate(wcurv)
    if(allocated(wstar_neo)) deallocate(wstar_neo)
    if(allocated(vparterm)) deallocate(vparterm)
    if(allocated(hneoc)) deallocate(hneoc)
    if(allocated(wdfac)) deallocate(wdfac)
    if(allocated(wstarfac)) deallocate(wstarfac)
    if(allocated(wdttpfac)) deallocate(wdttpfac)
  end subroutine finish_lowflow_terms

  !> FIXME : Add documentation
  !!
  !! @note This routine must be called after init_collisions as it relies on the
  !! variable use_le_layout which is set in init_collisions based on input
  !! from the user provided input file.
  subroutine init_lowflow(wdrift, wdriftttp, wstar, vpar, vpac, lf_default, lf_decompose, &
       driftknob, tpdriftknob)
    use constants, only: zi
    use species, only: spec, nspec
    use theta_grid, only: theta, ntgrid, delthet, gradpar, bmag, shat
    use theta_grid, only: gds23, gds24, gds24_noq
    use le_grids, only: energy, al, negrid, nlambda, forbid, init_map
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use kt_grids, only: theta0, ntheta0, naky, aky, akx
    use gs2_time, only: tunits, wunits, code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx, it_idx
    use run_parameters, only:rhostar, neo_test
    use file_utils, only: open_output_file, close_output_file
    use collisions, only: use_le_layout
    use mp, only: proc0, mp_abort

    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in out) :: wdrift
    real, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: wdriftttp, vpar, vpac
    real, dimension(:, :, :), intent(in) :: wstar
    logical, intent(in) :: lf_default, lf_decompose
    real, intent(in) :: driftknob, tpdriftknob
    integer, save :: neo_unit, neophi_unit, neovpth_unit
    integer :: it, ik, il, ie, is, isgn, iglo, ig
    real, dimension (:,:,:,:,:), allocatable :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    real, dimension (:,:,:), allocatable :: vpadhdec, dhdec, dhdxic, cdfac, hneovpth
    real, dimension (:), allocatable :: tmp7, tmp8, tmp9

    ! Initialise the curvature drift
    call init_wcurv(driftknob, tpdriftknob)

    allocate (vpadhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (dhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (dhdxic (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (cdfac (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    vpadhdec = 0. ; dhdec = 0. ; dhdxic = 0. ; cdfac = 0.

    if (.not. allocated(vparterm)) then
       allocate (vparterm(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (hneoc(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wstarfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdttpfac(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec,2))
       allocate (wstar_neo(-ntgrid:ntgrid,ntheta0,naky))
    end if
    vparterm = 0. ; wdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; wdttpfac = 0. ; wstar_neo = 0.

    allocate (tmp1(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp1 = 0.
    allocate (tmp2(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp2 = 0.
    allocate (tmp3(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp3 = 0.
    allocate (tmp4(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp4 = 0.
    allocate (tmp5(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp5 = 0.
    allocate (tmp6(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp6 = 0.
    allocate (tmp7(-ntgrid:ntgrid), tmp8(-ntgrid:ntgrid), tmp9(-ntgrid:ntgrid))
    tmp7 = 0. ; tmp8 = 0. ; tmp9 = 0.

    allocate (hneovpth(-ntgrid:ntgrid,negrid*nlambda,nspec)) ; hneovpth = 0.

    ! tmp1 is dH^{neo}/dE, tmp2 is dH^{neo}/dxi, tmp3 is vpa*dH^{neo}/dE,
    ! tmp4 is dH^{neo}/dr, tmp5 is dH^{neo}/dtheta, tmp6 is H^{neo}
    ! tmp7 is phi^{neo}/dr, tmp8 is dphi^{neo}/dtheta, and tmp9 phi^{neo}
    call get_lowflow_terms (theta, al, energy, bmag, tmp1, tmp2, tmp3, tmp4, &
         tmp5, tmp6, tmp7, tmp8, tmp9, lf_default, lf_decompose)

    if (proc0) then
       call open_output_file (neo_unit,".neodist")
       write (neo_unit,*) "# all quantities given at theta=0 for species 1"
       write (neo_unit,fmt='(10a12)') "# 1) vpa", "2) vpe", "3) energy", "4) vpa/v", &
            "5) dH/dE", "6) dH/dxi", "7) vpa*dH/dE", "8) dH/dr", "9) dH/dtheta", "10) H"
       ! Only write species 1:
       is = 1
       do isgn = 1, 2
          do il = 1, nlambda
             do ie = 1, negrid
                if (.not. forbid(0,il)) then
                   write (neo_unit,'(10e12.4)') sign(sqrt(energy(ie,is)*(1.-al(il)*bmag(0))),1.5-real(isgn)), &
                        sqrt(energy(ie,is)*al(il)*bmag(0)), energy(ie,is), &
                        sign(sqrt(1.-al(il)*bmag(0)),1.5-real(isgn)), &
                        tmp1(0,il,ie,isgn,is), tmp2(0,il,ie,isgn,is), tmp3(0,il,ie,isgn,is), &
                        tmp4(0,il,ie,isgn,is), tmp5(0,il,ie,isgn,is), tmp6(0,il,ie,isgn,is)
                end if
             end do
             write (neo_unit,*)
          end do
       end do
       call close_output_file (neo_unit)

       call open_output_file (neovpth_unit,".neothvp")

       ! Get Fneo(theta,vpa)
       call get_flux_vs_theta_vs_vpa (tmp6, hneovpth)

       write (neovpth_unit,'(3a12)') '1) theta', '2) vpa', '3) Fneo'
       ! Only write species 1:
       is = 1
       do ie = 1, negrid*nlambda
          do ig = -ntgrid, ntgrid
             write (neovpth_unit,'(3e12.4)') theta(ig), sqrt(energy(negrid,is))*(1.-2.*(ie-1)/real(negrid*nlambda-1)), &
                  hneovpth(ig,ie,is)
          end do
          write (neovpth_unit,*)
       end do

       call close_output_file (neovpth_unit)

       call open_output_file (neophi_unit,".neophi")
       write (neophi_unit,*) "# 1) theta, 2) dphi/dr, 3) dphi/dtheta, 4) phi"
       do ig = -ntgrid, ntgrid
          write (neophi_unit,'(4e14.5)') theta(ig), tmp7(ig), tmp8(ig), tmp9(ig)
       end do
       call close_output_file (neophi_unit)
    end if

    ! if set neo_test flag to .true. in input file, GS2 exits after writing out
    ! neoclassical quantities of interest
    if (neo_test) call mp_abort('stopping as neo_test=.true.')

    ! intialize mappings from g_lo to e_lo and lz_lo or to le_lo to facilitate
    ! energy and lambda derivatives in parallel nonlinearity, etc.
    if(use_le_layout) then
       call init_map (use_lz_layout=.false., use_e_layout=.false., use_le_layout=.true., test=.false.)
    else
       call init_map (use_lz_layout=.true., use_e_layout=.true., use_le_layout=.false., test=.false.)
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1,2
          dhdec(:,isgn,iglo) = tmp1(:,il,ie,isgn,is)
          dhdxic(:,isgn,iglo) = tmp2(:,il,ie,isgn,is)
          vpadhdec(:,isgn,iglo) = tmp3(:,il,ie,isgn,is)
          hneoc(:,isgn,iglo) = 1.0+tmp6(:,il,ie,isgn,is)
       end do

       ! get cell-centered (in theta) values

       ! this is the contribution from dH^{neo}/dtheta (part of v_E dot grad F^{neo})
       ! takes care of part of Eq. 60 in MAB's GS2 notes
       if (aky(ik) == 0.0) then
          wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.25*akx(it)/shat &
               *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
          wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.25*akx(it)/shat &
               *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       else
          wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*wunits(ik) &
               *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
          wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*wunits(ik) &
               *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       end if

       ! this is the contribution from dH^{neo}/dr (part of v_E dot grad F^{neo})
       ! takes care of part of Eq. 60 in MAB's GS2 notes
       wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + tmp4(:,il,ie,:,is)*code_dt*wunits(ik)

       if (.not. lf_default) then
          ! this is the contribution from the last term of the 2nd line of Eq. 43 in
          ! MAB's GS2 notes (arises because the NEO dist. fn. is given at v/vt, and
          ! the vt varies radially, so v_E . grad F1 must take this into account)
          ! If lf_default is true this is taken care of by multiplying the vt normalization
          ! of Chebyshev polynomial argument by appropriate temperature factor
          ! when constructing neoclassical distribution function (see lowflow.f90).
          ! Otherwise, the below line of code takes care of it.
          wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + code_dt*wunits(ik) &
               *spec(is)%tprim*energy(ie,is)*dhdec(:,:,iglo)
       end if

       wstarfac(ntgrid,:,iglo) = 0.0

       ! wdfac takes care of the gbdrift part of Eq. 60 of MAB's GS2 notes, as well
       ! as part of the curvature drift term of Eq. 54.  the other part of
       ! the curvature drift term of Eq. 54 is dealt with by cdfac below.
       ! no code_dt in wdfac because it multiplies wdrift, which has code_dt in it
       wdfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/energy(ie,is)**1.5 &
            - dhdec(-ntgrid:ntgrid-1,1,iglo)
       wdfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/energy(ie,is)**1.5 &
            - dhdec(-ntgrid:ntgrid-1,2,iglo)
       wdfac(ntgrid,:,iglo) = 0.0

       ! takes care of part of curvature drift term in Eq. 54 of MAB's GS2 notes.
       ! no code_dt in cdfac because it multiples wcurv, which has code_dt in it.
       cdfac(-ntgrid:ntgrid-1,1,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/sqrt(energy(ie,is))
       cdfac(-ntgrid:ntgrid-1,2,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/sqrt(energy(ie,is))
       cdfac(ntgrid,:,iglo) = 0.0

       ! this is the first term multiplying dF_1/dE in Eq. 42 of MAB's GS2 notes
       ! i.e. Z_s * e * vpa . grad phi^{tb} d(F1/F0)/dE
       ! note there will be a correction to this added below
       ! because actual term appearing in GKE ~ (1/F0) * d(F1)/dE
       ! also note that vpadhdec is vpa*d(F1/F0)/dE at fixed mu (not xi)
       vparterm(-ntgrid:ntgrid-1,1,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
            /delthet(-ntgrid:ntgrid-1) &
            * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
            * vpadhdec(-ntgrid:ntgrid-1,1,iglo)
       vparterm(-ntgrid:ntgrid-1,2,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
            /delthet(-ntgrid:ntgrid-1) &
            * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
            * vpadhdec(-ntgrid:ntgrid-1,2,iglo)
       vparterm(ntgrid,:,iglo) = 0.0

       if (aky(ik) == 0) then
          wstar_neo(-ntgrid:ntgrid-1,it,ik) = -0.25*code_dt*akx(it) &
               * tmp8(-ntgrid:ntgrid-1)*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))
       else
          wstar_neo(-ntgrid:ntgrid-1,it,ik) = -code_dt*wunits(ik)*(tmp7(-ntgrid:ntgrid-1) &
               + 0.5*tmp8(-ntgrid:ntgrid-1)*(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid) &
               + theta0(it,ik)*shat*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))))
       end if
       wstar_neo(ntgrid,it,ik) = 0.0
    end do

    ! vparterm is -2*vpar*(1+H^{neo}) - Ze*(vpa . grad phi^{tb})*(dH^{neo}/dE)
    ! note that vpar has contribution from v_{magnetic} . grad theta in it
    ! hneoc = 1 + H^{neo} below accounts for usual parallel streaming source term,
    ! as well as first of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
    vparterm = -2.0*vpar*hneoc + vparterm

    ! hneoc below accounts for usual v_magnetic source term,
    ! as well as second of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
    ! wdfac on RHS deals with grad-B drift part of Eq. 60, as well as part of
    ! curvature drift term in Eq. 54
    wdfac = wdfac + hneoc

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       wdfac(:,1,iglo) = wdfac(:,1,iglo)*wdrift(:,1,iglo) &
            + cdfac(:,1,iglo)*wcurv(:,iglo)
       wdfac(:,2,iglo) = wdfac(:,2,iglo)*wdrift(:,2,iglo) &
            + cdfac(:,2,iglo)*wcurv(:,iglo)
       wdttpfac(:,it,ik,ie,is,1) = wdfac(:,1,iglo)*wdriftttp(:, 1, iglo) &
            + cdfac(:,1,iglo)*wcurv(:,iglo)
       wdttpfac(:,it,ik,ie,is,2) = wdfac(:,2,iglo)*wdriftttp(:, 2, iglo) &
            + cdfac(:,2,iglo)*wcurv(:,iglo)

       ! add Z^2 * sqrt(m*T^3) * vpa * (bhat . grad(theta)) * d(phineo)/dth term to wdrift
       ! this modification of wdrift takes care of -g term in Eq. 67 of MB's gs2_notes
       ! arises from pulling 1/F0 inside dg/dE term
       wdrift(-ntgrid:ntgrid-1,1,iglo) = wdrift(-ntgrid:ntgrid-1,1,iglo) &
            - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,1,iglo) &
            * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
            * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))
       wdrift(-ntgrid:ntgrid-1,2,iglo) = wdrift(-ntgrid:ntgrid-1,2,iglo) &
            - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,2,iglo) &
            * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
            * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))

       ! hneoc below accounts for usual wstar term, as well as last of three terms
       ! multiplying F_1 in Eq. 42 of MAB'S GS2 notes
       ! note that hneo is necessary here because NEO's dist fn is
       ! normalized by F0(r) instead of F0 at center of simulation domain as in GS2
       wstarfac(:,:,iglo) = wstar(ik,ie,is)*hneoc(:,:,iglo) - wstarfac(:,:,iglo)
    end do

    deallocate (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)
    deallocate (vpadhdec,dhdec,dhdxic,cdfac,hneovpth)

  end subroutine init_lowflow

  !> FIXME : Add documentation
  subroutine init_wcurv(driftknob, tpdriftknob)
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, il_is_wfb, il_is_trapped
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use array_utils, only: zero_array
    implicit none
    real, intent(in) :: driftknob, tpdriftknob
    integer :: ig, ik, it, il, ie, is, iglo
    real :: scaling_knob
    logical,parameter :: debug = .false.

    if (.not. allocated(wcurv)) allocate (wcurv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    call zero_array(wcurv)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, it, ik, is, ig, scaling_knob) &
    !$OMP SHARED( wcurv, g_lo, tpdriftknob, driftknob, ntgrid, forbid) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il=il_idx(g_lo,iglo)
       it=it_idx(g_lo,iglo)
       ik=ik_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)

       if (il_is_wfb(il) .or. il_is_trapped(il)) then
          scaling_knob = tpdriftknob
       else
          scaling_knob = driftknob
       end if

       do ig = -ntgrid, ntgrid
          if (forbid(ig, il)) cycle
          wcurv(ig,iglo) = wcurv_func(ig, it, ik) * scaling_knob
       end do
    end do
    !$OMP END PARALLEL DO

    ! This should be weighted by bakdif to be completely consistent
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ig) &
    !$OMP SHARED(g_lo, ntgrid, wcurv) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid-1
          wcurv(ig,iglo) =  0.5*(wcurv(ig,iglo) + wcurv(ig+1,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    if (debug) write(6,*) 'init_wcurv: driftknob, tpdriftknob=',driftknob,tpdriftknob

  end subroutine init_wcurv

  !> FIXME : Add documentation
  pure real function wcurv_func (ig, it, ik)
    use theta_grid, only: cvdrift, cvdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use gs2_time, only: code_dt, wunits
    implicit none
    integer, intent (in) :: ig, ik, it

    if (aky(ik) == 0.0) then
       wcurv_func = akx(it)/shat &
            * cvdrift0(ig) * code_dt/2.0
    else
       wcurv_func = (cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
            *code_dt*wunits(ik)
    end if
  end function wcurv_func

  !> FIXME : Add documentation
  subroutine get_lowflow_terms (theta, al, energy, bmag, dHdEc, dHdxic, vpadHdEc, dHdrc, &
       dHdthc, hneoc, dphidrc, dphidthc, phi_neo, lf_default, lf_decompose)

    use mp, only: proc0, broadcast
    use le_grids, only: w, wl
    use theta_grid, only: delthet, jacob, ntgrid, Rplot
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use geometry, only: rmaj
    use species, only: spec
    use integration, only: trapezoidal_integration
    implicit none

    real, dimension (:), intent (in) :: theta, bmag
    real, dimension (:), intent (in) ::  al
    real, dimension (:,:), intent (in) :: energy
    real, dimension (:,:,:,:,:), intent (out) :: dHdec, dHdxic, dHdrc, dHdthc, vpadHdEc, hneoc
    real, dimension (:), intent (out) :: dphidthc, dphidrc, phi_neo
    logical, intent (in) :: lf_default, lf_decompose

    real :: rgeo2, drlt, drln

    real, dimension (:,:,:,:,:), allocatable :: hneo
    real, dimension (:,:,:,:), allocatable :: dHdxi, dHdE, vpadHdE, dHdr, dHdth
    real, dimension (:,:,:), allocatable :: legp
    real, dimension (:,:), allocatable :: xi, emax, dphidr
    real, dimension (:,:,:,:), allocatable :: chebyp1, chebyp2
    real, dimension (:), allocatable :: dl_over_b, drad_neo
    real, dimension (:,:), allocatable :: transport
    real, dimension (:), allocatable :: phineo2
    real, dimension (:,:), allocatable :: upar2, qpar2
    real, dimension (:,:,:), allocatable :: upar, qpar, afac, bfac
    real, dimension (:,:), allocatable :: pflx, qflx, vflx, qparB, upar1, uparB, upar_over_B
    real, dimension (:), allocatable ::  radius, jboot, phi2, vtor, upar0
    real, dimension (size(energy)) :: w_r

    integer :: il, ie, is, ns, nc, nl, nr, ig, ixi, ir, ir_loc
    integer :: ntheta, nlambda, nenergy, nxi
    integer, save :: neo_unit, neot_unit
    real :: average_weight
    logical, dimension (:,:), allocatable :: forbid

    ntheta = size(theta)
    nlambda = size(al)
    nenergy = size(energy)
    nxi = 2*nlambda-1

    if (.not. allocated(xi)) allocate (xi(ntheta,nxi))
    if (.not. allocated(forbid)) allocate(forbid(ntheta,nxi))
    forbid = .false. ; xi = 0.0
    !Note, xi is not the same order as le_grids version. Should we just use le_grids
    ! one? If so would need to reorder indexing elsewhere in this module
    ! get xi=vpar/v from lambda=(mu/E)*Bnorm
    do il = 1, nlambda-1
       do ig = 1, ntheta
          if (1.0-al(il)*bmag(ig) < 0.0) then
             forbid(ig,il) = .true.
             forbid(ig,2*nlambda-il) = .true.
          else
             xi(ig,il) = -sqrt(1.0-al(il)*bmag(ig))
             xi(ig,2*nlambda-il) = -xi(ig,il)
          end if
       end do
    end do

    ! ns is nspec, nc is nenergy, nl is nlambda,
    ! nr is number of radii
    ! taken from neo
    call read_neocoefs (theta, ns, nc, nl, nr, ir_loc)

    if (.not. allocated(chebyp1)) allocate (chebyp1(nr,nenergy,0:nc-1,ns), chebyp2(nr,nenergy,0:nc-1,ns))
    if (.not. allocated(legp)) allocate (legp(ntheta,nxi,0:nl+1))
    if (.not. allocated(hneo)) allocate (hneo(nr,ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdr)) allocate (   dHdr(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdth)) allocate (  dHdth(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdxi)) allocate (  dHdxi(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdE)) allocate (   dHdE(ntheta,nxi,nenergy,ns))
    if (.not. allocated(vpadHdE)) allocate (vpadHdE(ntheta,nxi,nenergy,ns))
    if (.not. allocated (dphidr)) allocate (dphidr(ntheta,ns))
    if (.not. allocated (dl_over_b)) allocate (dl_over_b(ntheta))
    if (.not. allocated (emax)) allocate (emax(nr,ns))

    ! better to be taken from neo
    if (proc0) write (*,*) '# make sure ENERGY_MAX=16.0 in NEO INPUT file'
    ! emax is ENERGY_MAX (from NEO input file) times v_{ts}^2

    ! two ways to deal with radial variation of temperature in
    ! NEO energy variable: 1st (lf_default) is to take it into
    ! account when constructing distribution function, i.e. H(EGS2,r);
    ! 2nd is to construct H(ENEO) and deal with it later by
    ! including a temp gradient term in wstarfac in dist_fn.fpp
    if (lf_default) then
       ! v_{ts} is a function of radius, so we need to convert emax
       ! to its equivalent value using v_{ts} from the center radius.
       ! this is necessary because we will be taking radial derivatives of
       ! the distribution function with v fixed, not v/v_t(r) fixed.
       do is = 1, ns
          emax(2,is) = 16.0 ! this is EMAX for center grid point
          emax(1,is) = emax(2,is)*temp_neo(1,is)
          emax(3,is) = emax(2,is)*temp_neo(3,is)
       end do
    else
       emax = 16.0
    end if

    allocate (drad_neo(nr)) ; drad_neo = 0.
    drad_neo = rad_neo - rad_neo(2)

    ! get legendre polynomials on gs2 pitch-angle grid
    legp = 0.0
    do ixi = 1, nxi
       do ig = 1, ntheta
          if (.not. forbid(ig,ixi)) call legendre (xi(ig,ixi), legp(ig,ixi,:))
       end do
    end do

    ! get chebyshev polynomials of the first and second kinds on the gs2 energy grid
    ! note the argument of the chebyshev polynomials is z = 2 * sqrt(energy/emax) - 1
    do is = 1, ns
       do ie = 1, nenergy
          do ir = 1, nr
             call chebyshev (zfnc(energy(ie,is),emax(ir,is)), chebyp1(ir,ie,:,is), 1)
             call chebyshev (zfnc(energy(ie,is),emax(ir,is)), chebyp2(ir,ie,:,is), 2)
          end do
       end do
    end do

    ! first, get hneo = F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
    ! H_1/F_0 is constructed from Legendre and Chebyshev polynomials
    do is = 1, ns
       do ig = 1, ntheta
          do ir = 1, nr
             call get_H (coefs(ir,ig,:,:,is), legp(ig,:,:), chebyp1(ir,:,:,is), &
                  hneo(ir,ig,:,:,is), phineo(ir,ig), ig, is)
          end do
       end do
    end do

    ! calculate and write to file some useful quantities like Upar and Qpar
    ! Note dl_over_b may benefit from using midpoint jacob. Only used for rgeo2.
    dl_over_b = delthet*jacob
    dl_over_b = dl_over_b/sum(dl_over_b)
    average_weight = 1.0 / trapezoidal_integration(theta, jacob)

    if (.not. allocated(pflx)) allocate (pflx(nr,ns), qflx(nr,ns), vflx(nr,ns), upar1(nr,ns))
    if (.not. allocated(upar)) then
       allocate (phineo2(nr))
       allocate (upar2(nr,ns), qpar2(nr,ns))
       allocate (uparB(nr,ns), qparB(nr,ns))
       allocate (upar(nr,ns,ntheta), qpar(nr,ns,ntheta))
       allocate (upar_over_b(nr,ns))
       upar = 0. ; qpar = 0.
       upar2 = 0. ; qpar2 = 0. ; phineo2 = 0.
       uparB = 0. ; qparB = 0.
       upar_over_B = 0.
    end if
    if (.not. allocated(transport)) allocate (transport(nr,(5+ns*8)))  !JPL, NEO output for all radius
    if (.not. allocated(radius)) allocate (radius(nr), phi2(nr), jboot(nr), vtor(nr), upar0(nr))
    if (.not. allocated(mach_lab)) allocate (mach_lab(ns))

    if (proc0) then
       call get_unused_unit (neo_unit)
       ! read the diagnosis in NEO from "neo_transport.out"
       open (unit=neo_unit, file='neo_transport.out', status='old', action='read')
       do ir=1, nr
          read (neo_unit,*) transport(ir,:)
          radius(ir) = transport(ir,1)
          phi2(ir) = transport(ir,2)
          jboot(ir) = transport(ir,3)
          vtor(ir) = transport(ir,4)
          upar0(ir) = transport(ir,5)
          do is = 1, ns
             pflx(ir,is) = transport(ir,5+(is-1)*8+1)/sqrt(2.0)  ! NEO uses vt=sqrt(T/m) while GS2 uses vt=sqrt(2T/m)
             qflx(ir,is) = transport(ir,5+(is-1)*8+2)/(sqrt(2.0)**3)
             vflx(ir,is) = transport(ir,5+(is-1)*8+3)/(sqrt(2.0)**2)
             upar1(ir,is) = transport(ir,5+(is-1)*8+4)/sqrt(2.0)
          end do
       end do
       close (neo_unit)

       ! evaluate the neoclassical parallel heat flux and parallel velcotiy in GS2 geometry
       do ir=1, nr
          do is=1, ns
             do ie = 1, nenergy
                w_r(ie) = w(ie,is) * exp(energy(ie,is)*(1.0-1.0/temp_neo(ir,is))) / temp_neo(ir,is)**1.5
                do ixi = 1, nxi
                   il = min(ixi, nxi+1-ixi)
                   do ig = 1, ntheta
                      upar(ir,is,ig) = upar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie,is))*xi(ig,ixi) * w_r(ie) * wl(-ntgrid+ig-1,il)
                      qpar(ir,is,ig) = qpar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie,is))*xi(ig,ixi)*(energy(ie,is)-2.5*temp_neo(ir,is)) &
                           * dens_neo(ir,is) * w_r(ie) * wl(-ntgrid+ig-1,il)
                   end do
                end do
             end do
          end do
       end do

       do is = 1, ns
          do ir = 1, nr
             ! Note we could use field_line_average from theta grid if
             ! we knew that theta and jacob are always provided by
             ! theta_grid. As we pass in theta here we better not
             ! assume it comes from theta_grid.  Note however we
             ! always use jacob from theta_grid, so if we do pass in a
             ! different theta this may be inconsistent.
             uparB(ir,is) = trapezoidal_integration(theta, upar(ir,is,:)*jacob*bmag)*average_weight
             qparB(ir,is) = trapezoidal_integration(theta, qpar(ir,is,:)*jacob*bmag)*average_weight
             upar_over_B(ir,is) = trapezoidal_integration(theta, upar(ir,is,:)*jacob/bmag)*average_weight
             upar2(ir,is) = trapezoidal_integration(theta, upar(ir,is,:)**2*jacob)*average_weight
             qpar2(ir,is) = trapezoidal_integration(theta, qpar(ir,is,:)**2*jacob)*average_weight
             phineo2(ir) = trapezoidal_integration(theta, phineo(ir,:)**2*jacob)*average_weight
          end do
       end do

       rgeo2 = 0.
       do ig = 1, ntheta
          rgeo2 = rgeo2 + dl_over_b(ig)*Rplot(-ntgrid+ig-1)**2
       end do

       ! note that for total toroidal angular momentum (including diamagnetic and ExB) to be zero,
       ! upar_over_b/rgeo2 = -(omega_phi * a / v_{t,ref}) * (a/R0)

       call open_output_file (neot_unit,".neotransp")
       ! print out in ".neotransp" file
       write (neot_unit,fmt='(a14,a8,12a14)') "# 1) rad", "2) spec", "3) pflx", "4) qflx", "5) vflx", &
            "6) qparflx", "7) uparB(GS2)", "8) upar1(NEO)", "9) <phi**2>", "10) bootstrap", "11) upar_o_B", &
            "12) upar2", "13) qpar2", "14) phi2"
       do ir=1, nr
          do is = 1, ns
             write (neot_unit,fmt='(e14.5,i8,12e14.5)') radius(ir), is, pflx(ir,is), qflx(ir,is), vflx(ir,is), &
                  qparB(ir,is),uparB(ir,is), upar1(ir,is)*sqrt(temp_neo(ir,is)), &
                  phi2(ir), jboot(ir), rmaj*upar_over_B(ir,is)/rgeo2, upar2(ir,is), qpar2(ir,is), phineo2(ir)
             if (ir==2)  mach_lab(is)=rmaj*upar_over_B(ir,is)/rgeo2
          end do
       end do
       call close_output_file (neot_unit)
    end if

    call broadcast(mach_lab)

    do ir = 1, nr
       do is = 1, ns
          call broadcast (upar(ir,is,:))
          call broadcast (qpar(ir,is,:))
       end do
    end do

    allocate (afac(nr,ns,ntheta), bfac(nr,ns,ntheta))

    ! Here only to set up F_1/F_0 = xi*(v/vth)*(a + b*(v/vth)^2)
    if (lf_decompose) then
       do is = 1, ns
          do ir = 1, nr
             drlt = drad_neo(ir)*spec(is)%tprim ; drln = drad_neo(ir)*spec(is)%fprim
             afac(ir,is,:) = 2.*(upar(ir,is,:)-qpar(ir,is,:)/(1.0-drln) - drlt*upar(ir,is,:))/(1.0-drlt)**2
             bfac(ir,is,:) = 0.8*qpar(ir,is,:) / ((1.0-drln) * (1.0-drlt)**3)
          end do
       end do

       coefs = 0.0
       do ig = 1, ntheta
          coefs(:,ig,1,0,:) = sqrt(emax)*(0.5*afac(:,:,ig) + 0.3125*bfac(:,:,ig)*emax)
          coefs(:,ig,1,1,:) = 2.*sqrt(emax)*(0.25*afac(:,:,ig) + 0.234375*bfac(:,:,ig)*emax)
          coefs(:,ig,1,2,:) = 2.*0.09375*emax**1.5*bfac(:,:,ig)
          coefs(:,ig,1,3,:) = 2.*0.015625*emax**1.5*bfac(:,:,ig)
       end do

       ! first, get hneo = F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
       ! H_1/F_0 is constructed from Legendre and Chebyshev polynomials
       do is = 1, ns
          do ig = 1, ntheta
             do ir = 1, nr
                call get_H (coefs(ir,ig,:,:,is), legp(ig,:,:), chebyp1(ir,:,:,is), &
                     hneo(ir,ig,:,:,is), phineo(ir,ig), ig, is)
             end do
          end do
       end do

       upar = 0. ; qpar = 0.
       do ir = 1, nr
          do is = 1, ns
             drlt = drad_neo(ir)*spec(is)%tprim ; drln = drad_neo(ir)*spec(is)%fprim
             do ie = 1, nenergy
                w_r(ie) = w(ie,is) * exp(energy(ie,is)*(1.0-1.0/(1.0-drlt))) / (1.0-drlt)**1.5
                do ixi = 1, nxi
                   il = min(ixi, nxi+1-ixi)
                   do ig = 1, ntheta
                      upar(ir,is,ig) = upar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie,is))*xi(ig,ixi) * w_r(ie) * wl(-ntgrid+ig-1,il)
                      qpar(ir,is,ig) = qpar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie,is))*xi(ig,ixi)*(energy(ie,is)-2.5*(1.0-drlt)) &
                           * (1.0-drln) * w_r(ie) * wl(-ntgrid+ig-1,il)
                   end do
                end do
             end do
          end do
       end do

       do is = 1, ns
          do ir = 1, nr
             upar2(ir,is) = trapezoidal_integration(theta, upar(ir,is,:)**2*jacob)*average_weight
             qpar2(ir,is) = trapezoidal_integration(theta, qpar(ir,is,:)**2*jacob)*average_weight
          end do
       end do

    end if

    ! get dphi/dr at GS2 radius (center of 3 NEO radii)
    ! note that radgrad takes derivative of e*phi/Tref,
    ! not the derivative of Z_s*e*phi/T_s
    do ig = 1, ntheta
       ! Note that we are neglecting the term prop. to dTref/dr, which is fine
       ! as long as Tref in NEO is the same for all radii.
       ! This will require a changing TEMP with radius if the species temperature
       ! changes with radius.
       call get_radgrad (phineo(:,ig), rad_neo, ir_loc, dphidr(ig,1))
    end do

    ! get theta derivative of F_1/F_0 at fixed (xi,E)
    do is = 1, ns
       do ie = 0, nc-1
          do ixi = 0, nl
             ! get theta derivative of spectral coefficients of F_1/F_0
             call get_thgrad (coefs(ir_loc,:,ixi,ie,is), theta, dcoefsdth(:,ixi,ie))
          end do
       end do
       ! note that dphidth is calculated in read_neocoefs
       do ig = 1, ntheta
          call get_gradH (dcoefsdth(ig,:,:), dphidth(ig)*spec(is)%zt, legp(ig,:,:), chebyp1(ir_loc,:,:,is), dHdth(ig,:,:,is))
       end do
    end do

    do is = 1, ns
       do ig = 1, ntheta
          ! get dH/dxi at fixed E and dH/dE at fixed xi
          call get_dHdxi (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), xi(ig,:), dHdxi(ig,:,:,is))
          call get_dHdE (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), chebyp2(ir_loc,:,:,is), &
               energy(:,is), emax(ir_loc,is), dHdE(ig,:,:,is))
       end do
    end do

    ! get radial derivative of F_1/F_0 at fixed (xi,E)
    ! note that F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
    do is = 1, ns
       do ig = 1, ntheta
          do ie = 1, nenergy
             do ixi = 1, nxi
                call get_radgrad (hneo(:,ig,ixi,ie,is), rad_neo, ir_loc, dHdr(ig,ixi,ie,is))
             end do
          end do
       end do
    end do

    dphidrc(1:ntheta-1) = 0.5*(dphidr(1:ntheta-1,1) + dphidr(2:ntheta,1))
    dphidrc(ntheta) = dphidrc(1)
    dphidthc(1:ntheta-1) = 0.5*(dphidth(1:ntheta-1) + dphidth(2:ntheta))
    dphidthc(ntheta) = dphidth(1)

    phi_neo = phineo(ir_loc,:)

    ! vpadHdE is the derivative of F1/F0 with respect to energy at fixed mu (not xi)
    do ie = 1, nenergy
       do ixi = 1, nxi
          do ig = 1, ntheta
             vpadHdE(ig,ixi,ie,:) = sqrt(energy(ie,is))*xi(ig,ixi)*dHdE(ig,ixi,ie,:) &
                  + (1.-xi(ig,ixi)**2)*dHdxi(ig,ixi,ie,:)/(2.*sqrt(energy(ie,is)))
          end do
       end do
    end do

    do is = 1, ns
       do ie = 1, nenergy
          do il = 1, nlambda
             ! should be weighted with backdif for consistency -- MAB
             dHdrc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdr(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdr(2:ntheta,2*nlambda-il,ie,is))
             dHdrc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdr(1:ntheta-1,il,ie,is) &
                  + dHdr(2:ntheta,il,ie,is))
             dHdthc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdth(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdth(2:ntheta,2*nlambda-il,ie,is))
             dHdthc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdth(1:ntheta-1,il,ie,is) &
                  + dHdth(2:ntheta,il,ie,is))
             dHdEc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdE(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdE(2:ntheta,2*nlambda-il,ie,is))
             dHdEc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdE(1:ntheta-1,il,ie,is) &
                  + dHdE(2:ntheta,il,ie,is))
             dHdxic(1:ntheta-1,il,ie,1,is) = 0.5*(dHdxi(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdxi(2:ntheta,2*nlambda-il,ie,is))
             dHdxic(1:ntheta-1,il,ie,2,is) = 0.5*(dHdxi(1:ntheta-1,il,ie,is) &
                  + dHdxi(2:ntheta,il,ie,is))
             vpadHdEc(1:ntheta-1,il,ie,1,is) = 0.5*(vpadHdE(1:ntheta-1,2*nlambda-il,ie,is) &
                  + vpadHdE(2:ntheta,2*nlambda-il,ie,is))
             vpadHdEc(1:ntheta-1,il,ie,2,is) = 0.5*(vpadHdE(1:ntheta-1,il,ie,is) &
                  + vpadHdE(2:ntheta,il,ie,is))
             hneoc(1:ntheta-1,il,ie,1,is) = 0.5*(hneo(ir_loc,1:ntheta-1,2*nlambda-il,ie,is) &
                  + hneo(ir_loc,2:ntheta,2*nlambda-il,ie,is))
             hneoc(1:ntheta-1,il,ie,2,is) = 0.5*(hneo(ir_loc,1:ntheta-1,il,ie,is) &
                  + hneo(ir_loc,2:ntheta,il,ie,is))
          end do
       end do
    end do
    dHdrc(ntheta,:,:,:,:) = 0.0 ; dHdthc(ntheta,:,:,:,:) = 0.0 ; vpadHdEc(ntheta,:,:,:,:) = 0.0
    dHdEc(ntheta,:,:,:,:) = 0.0 ; dHdxic(ntheta,:,:,:,:) = 0.0 ; hneoc(ntheta,:,:,:,:) = 0.0
    dphidrc(ntheta) = 0.0 ; dphidthc(ntheta) = 0.0

    deallocate (xi, emax, chebyp1, chebyp2, legp, dHdr, dHdth, dHdxi, dHdE, vpadHdE, hneo, forbid)
    deallocate (coefs, drad_neo, dens_neo, temp_neo)
    deallocate (dphidr, dl_over_b, transport)
    deallocate (pflx, qflx, vflx, upar1, qparB)
    deallocate (upar, qpar, uparB)
    deallocate (upar2, qpar2)
    deallocate (afac, bfac)

  end subroutine get_lowflow_terms

  !> FIXME : Add documentation
  function zfnc (enrgy, enrgymax)

    implicit none

    real, intent(in) :: enrgy, enrgymax
    real :: zfnc

    zfnc = 2.*sqrt(enrgy/enrgymax)-1.

    return

  end function zfnc

  !> FIXME : Add documentation
  !!
  !! @note knd = 1 (2) for cheb polys of first (second) kind
  subroutine chebyshev (x, chebyp, knd)

    implicit none

    integer, intent (in) :: knd
    real, intent (in) :: x
    real, dimension (0:), intent (out) :: chebyp

    integer :: n, idx

    n = size(chebyp)-1

    chebyp(0) = 1.0

    if (knd == 2) then
       chebyp(1) = 2.*x
    else
       chebyp(1) = x
    end if
    do idx = 2, n
       chebyp(idx) = 2.*x*chebyp(idx-1) - chebyp(idx-2)
    end do

  end subroutine chebyshev

  !> FIXME : Add documentation
  subroutine legendre (x, legp)

    implicit none

    real, intent (in) :: x
    real, dimension (0:), intent (out) :: legp

    integer :: n, idx

    n = size(legp)-1

    legp(0) = 1.0
    legp(1) = x

    do idx = 2, n
       legp(idx) = ((2.*idx-1.)*x*legp(idx-1) + (1.-idx)*legp(idx-2))/idx
    end do

  end subroutine legendre

  !> FIXME : Add documentation
  subroutine get_H (gjk, legdre, chebyshv, h, phi, ig, is)

    use species, only: spec

    implicit none

    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre, chebyshv
    real, intent (in) :: phi
    real, dimension (:,:), intent (out) :: h
    integer, intent (in) :: ig, is

    integer :: ix, ij, ik

    h = 0.0

    do ix = 1, size(h,1)
       do ik = 0, size(gjk,2)-1
          do ij = 0, size(gjk,1)-1
             h(ix,:) = h(ix,:) + chebyshv(:,ik)*legdre(ix,ij)*gjk(ij,ik)
          end do
       end do
    end do

    ! F_1s = H_1s - Z_s*e*Phi_1/T_s
    ! want Z_s*e*Phi/T_s, but phi from neo is e*Phi/Tnorm
    ! at center (GS2) radius, this is a simple multiply by Z_s * Tref/Ts
    ! note that this assumes the Tref used in NEO is the same as that used in GS2
    h = h - spec(is)%zt*phi

  end subroutine get_H

  !> FIXME : Add documentation
  subroutine get_dHdxi (gjk, legdre, chebyshv, x, dH)

    implicit none

    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:), intent (in) :: x
    real, dimension (:,0:), intent (in) :: chebyshv
    real, dimension (:,:), intent (out) :: dH

    integer :: ij, ik, ix

    dH = 0.0

    ! calculate dH/dxi = sum_{j,k} g_{j,k} * T_{k}(z) * d(P_{j}(xi))/d(xi), where z=2*sqrt(E/EMAX)-1
    ! note that dP_{j}/dxi = (j+1)*(P_{j+1} - xi*P_{j})/(xi^2-1)
    do ix = 1, size(x)
       do ik = 0, size(gjk,2)-1
          do ij = 0, size(gjk,1)-1
             dH(ix,:) = dH(ix,:) + (chebyshv(:,ik)*(ij+1)*(legdre(ix,ij+1)-x(ix)*legdre(ix,ij))*gjk(ij,ik))/(x(ix)**2-1.)
          end do
       end do
    end do

  end subroutine get_dHdxi

  !> FIXME : Add documentation
  subroutine get_dHdE (gjk, legdre, chebyshv1, chebyshv2, x, xmax, dH)

    implicit none

    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:), intent (in) :: x
    real, intent (in) :: xmax
    real, dimension (:,0:), intent (in) :: chebyshv1, chebyshv2
    real, dimension (:,:), intent (out) :: dH

    integer :: ij, ik, ix
    real :: cfac

    dH = 0.0

    ! dH/dE at fixed xi is sum_{j,k} c_{j,k}*P_j(xi)*d(T_{n}(z))/dz / sqrt(E*EMAX)
    ! z=2*sqrt(E/EMAX)-1
    ! note that d(T_{k})/dz = k*U_{k-1}, where U_k is the Chebyshev polynomial of the 2nd kind
    do ix = 1, size(x)
       do ik = 0, size(gjk,2)-1
          if (ik==0) then
             ! cfac nonzero here if hneo = F_1 instead of F_1/F_0
             !             cfac = -chebyshv1(ix,ik)
             cfac = 0.0
          else
             ! extra chebyshv1 term not needed since hneo = F_1/F_0 instead of F_1
             !             cfac = ik/sqrt(xmax*x(ix))*chebyshv2(ix,ik-1)-chebyshv1(ix,ik)
             cfac = ik/sqrt(xmax*x(ix))*chebyshv2(ix,ik-1)
          end if
          do ij = 0, size(gjk,1)-1
             dH(:,ix) = dH(:,ix) + (legdre(:,ij)*cfac*gjk(ij,ik))
          end do
       end do
    end do

  end subroutine get_dHdE

  !> Get radial derivative at center of 3 points
  subroutine get_radgrad (h, rad, ir, dh)

    implicit none

    real, dimension (:), intent (in) :: h, rad
    integer, intent (in) :: ir
    real, intent (out) :: dh

    dh = (h(ir+1)-h(ir-1))/(rad(ir+1)-rad(ir-1))

  end subroutine get_radgrad

  !> FIXME : Add documentation
  subroutine get_thgrad (h, th, dh)

    implicit none

    real, dimension (:), intent (in) :: h, th
    real, dimension (:), intent (out) :: dh

    integer :: ig, nth

    nth = size(th)

    do ig = 2, nth-1
       dh(ig) = (h(ig+1)-h(ig-1))/(th(ig+1)-th(ig-1))
    end do
    ! note that H_neo is periodic in theta
    dh(1) = (h(2)-h(nth))/(2.*(th(2)-th(1)))
    dh(nth) = (h(1)-h(nth-1))/(2.*(th(nth)-th(nth-1)))

  end subroutine get_thgrad

  !> FIXME : Add documentation
  subroutine get_gradH (dgjk, dphi, legdre, chebyshv, dH)

    implicit none

    real, dimension (0:,0:), intent (in) :: dgjk
    real, intent (in) :: dphi
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:,0:), intent (in) :: chebyshv
    real, dimension (:,:), intent (out) :: dH

    integer :: ij, ik, ix

    dH = 0.0

    ! calculate dH/drho or dH/dtheta at fixed (xi,E) as sum_{j,k} d(g_{j,k})/drho * T_{k}(z) * P_{j}(xi)
    ! similarly for dH/dtheta
    ! where z=2*sqrt(E/EMAX)-1
    do ix = 1, size(legdre,1)
       do ik = 0, size(dgjk,2)-1
          do ij = 0, size(dgjk,1)-1
             dH(ix,:) = dH(ix,:) + chebyshv(:,ik)*legdre(ix,ij)*dgjk(ij,ik)
          end do
       end do
    end do

    dH = dH - dphi

  end subroutine get_gradH

  !> FIXME : Add documentation
  subroutine read_neocoefs (theta, nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo)

    use mp, only: proc0, broadcast
    use splines, only: lf_spline
    use file_utils, only: get_unused_unit
    use constants, only: pi,twopi

    implicit none

    real, dimension (:), intent (in) :: theta
    integer, intent (out) :: nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo

    integer :: is, ik, ij, ig, ir, idx, ntheta, ntheta_neo
    integer :: ptheta_neo, ip
    integer :: ngrid_per_rad
    integer, save :: neo_unit, neof_unit, neophi_unit

    real, dimension (:), allocatable :: tmp, neo_coefs, dum, neo_phi, dneo_phi
    real, dimension (:), allocatable :: theta_neo_ext,neo_coefs_ext, neo_phi_ext, dneo_phi_ext !JPL
    ntheta = size(theta)

    ! for now, set ir_neo by hand, but best to derive it from neo output in future
    ir_neo = 2

    if (proc0) then
       call get_unused_unit (neo_unit)

       ! read in number of grid points from neo's grid.out file
       open (unit=neo_unit, file='neo_grid.out', status="old", action="read")
       read (neo_unit,*) nspec_neo
       read (neo_unit,*) nenergy_neo
       read (neo_unit,*) nxi_neo
       read (neo_unit,*) ntheta_neo
       if (.not. allocated(theta_neo)) allocate (theta_neo(ntheta_neo))
       do ig = 1, ntheta_neo
          read (neo_unit,*) theta_neo (ig)
       end do
       read (neo_unit,*) nrad_neo
       if (.not. allocated(rad_neo)) allocate (rad_neo(nrad_neo))
       if (.not. allocated(dens_neo)) allocate (dens_neo(nrad_neo,nspec_neo))
       if (.not. allocated(temp_neo)) allocate (temp_neo(nrad_neo,nspec_neo))
       do ir = 1, nrad_neo
          read (neo_unit,*) rad_neo(ir), dens_neo(ir,:), temp_neo(ir,:)
       end do
       close (neo_unit)
       dens_neo = dens_neo / spread(dens_neo(ir_neo,:),1,nrad_neo)
       temp_neo = temp_neo / spread(temp_neo(ir_neo,:),1,nrad_neo)
    end if

    call broadcast (nspec_neo)
    call broadcast (nenergy_neo)
    call broadcast (nxi_neo)
    call broadcast (ntheta_neo)
    call broadcast (nrad_neo)
    if (.not. allocated(theta_neo)) allocate (theta_neo(ntheta_neo))
    if (.not. allocated(rad_neo)) allocate (rad_neo(nrad_neo))
    if (.not. allocated(dens_neo)) allocate (dens_neo(nrad_neo,nspec_neo))
    if (.not. allocated(temp_neo)) allocate (temp_neo(nrad_neo,nspec_neo))
    call broadcast (theta_neo)
    call broadcast (rad_neo)
    do ir = 1, nrad_neo
       call broadcast (dens_neo)
       call broadcast (temp_neo)
    end do

    ! define number of grid points per radial location to make file i/o easier
    ngrid_per_rad = ntheta_neo*(nxi_neo+1)*nenergy_neo*nspec_neo

    if (.not. allocated(tmp)) allocate (tmp(ngrid_per_rad*nrad_neo))
    if (.not. allocated(neo_coefs)) allocate (neo_coefs(ntheta_neo), neo_phi(ntheta_neo), dneo_phi(ntheta_neo))
    if (.not. allocated(dum)) allocate (dum(ntheta))
    if (.not. allocated(coefs)) allocate (coefs(nrad_neo,ntheta,0:nxi_neo,0:nenergy_neo-1,nspec_neo))
    if (.not. allocated(phineo)) allocate (phineo(nrad_neo,ntheta))
    if (.not. allocated(dphidth)) allocate (dphidth(ntheta))
    if (.not. allocated(dcoefsdr)) then
       allocate (dcoefsdr(0:nxi_neo,0:nenergy_neo-1)) ; dcoefsdr = 0.
       allocate (dcoefsdth(ntheta,0:nxi_neo,0:nenergy_neo-1)) ; dcoefsdth = 0.
    end if

    if (proc0) then

       !JPL: if the range of theta grid in GS2 is beyond [-pi:pi](e.g. "nperiod > 1"),
       !     extend the theta grid in NEO ([-pi:pi]) to the periodic one in the theta range of GS2.
       ptheta_neo=ceiling((maxval(theta)+pi)/twopi)  !the period of 2pi in theta range
       if (ptheta_neo .gt. 1) then
          if (.not. allocated(theta_neo_ext)) allocate (theta_neo_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(neo_coefs_ext)) allocate (neo_coefs_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(neo_phi_ext)) allocate (neo_phi_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(dneo_phi_ext)) allocate (dneo_phi_ext(ntheta_neo*(2*ptheta_neo-1)))
       end if

       ! read in H1^{nc} (adiabatic piece of F1^{nc}) from neo's f.out file
       call get_unused_unit (neof_unit)
       open (unit=neof_unit, file='neo_f.out', status="old", action="read")

       do ir = 1, nrad_neo
          read (neof_unit,*) tmp(ngrid_per_rad*(ir-1)+1:ngrid_per_rad*ir)
       end do

       idx = 1
       do ir = 1, nrad_neo
          do is = 1, nspec_neo
             do ik = 0, nenergy_neo-1
                do ij = 0, nxi_neo
                   do ig = 1, ntheta_neo
                      neo_coefs(ig) = tmp(idx)
                      idx = idx+1
                   end do

                   if ( ptheta_neo .gt. 1) then
                      do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                         theta_neo_ext(ip)=theta_neo(mod((ip-1),ntheta_neo)+1)+twopi*(int((ip-1)/ntheta_neo)-(ptheta_neo-1))
                         neo_coefs_ext(ip)=neo_coefs(mod((ip-1),ntheta_neo)+1)
                      end do
                      call lf_spline (theta_neo_ext, neo_coefs_ext, theta, coefs(ir,:,ij,ik,is), dum)  !JPL
                   else
                      ! need to interpolate coefficients from neo's theta grid to gs2's
                      call lf_spline (theta_neo, neo_coefs, theta, coefs(ir,:,ij,ik,is), dum)
                   end if
                end do
             end do
          end do
       end do

       close (neof_unit)
    end if

    deallocate (tmp)

    if (.not. allocated(tmp)) allocate (tmp(ntheta_neo*nrad_neo))

    if (proc0) then
       ! read in phi1^{nc} from neo's phi.out file
       call get_unused_unit (neophi_unit)
       open (unit=neophi_unit, file='neo_phi.out', status="old", action="read")

       do ir = 1, nrad_neo
          read (neophi_unit,*) tmp((ir-1)*ntheta_neo+1:ir*ntheta_neo)
       end do

       idx = 1
       do ir = 1, nrad_neo
          do ig = 1, ntheta_neo
             neo_phi(ig) = tmp(idx)
             idx = idx+1
          end do

          ! TMP FOR TESTING -- MAB
          !          neo_phi = 0.

          ! need to interpolate coefficients from neo's theta grid to gs2's
          if ( ptheta_neo .gt. 1) then
             do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                neo_phi_ext(ip)=neo_phi(mod((ip-1),ntheta_neo)+1)
             end do
             call lf_spline (theta_neo_ext, neo_phi_ext, theta, phineo(ir,:), dum)  !JPL
          else
             call lf_spline (theta_neo, neo_phi, theta, phineo(ir,:), dum)
          end if

          ! at central radius, calculate dphi/dth and interpolate onto gs2 grid
          if (ir == 2) then
             if ( ptheta_neo .gt. 1) then
                do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                   dneo_phi_ext(ip)=dneo_phi(mod((ip-1),ntheta_neo)+1)
                end do
                call get_thgrad (neo_phi_ext, theta_neo_ext, dneo_phi_ext)
                call lf_spline (theta_neo_ext, dneo_phi_ext, theta, dphidth, dum)
             else
                call get_thgrad (neo_phi, theta_neo, dneo_phi)
                call lf_spline (theta_neo, dneo_phi, theta, dphidth, dum)
             end if
          end if
       end do

       close (neophi_unit)

    end if

    call broadcast (dphidth)
    do ir = 1, nrad_neo
       call broadcast (phineo(ir,:))
    end do
    do ir = 1, nrad_neo
       do is = 1, nspec_neo
          do ik = 0, nenergy_neo-1
             do ij = 0, nxi_neo
                call broadcast (coefs(ir,:,ij,ik,is))
             end do
          end do
       end do
    end do

    deallocate (tmp, neo_coefs, theta_neo, neo_phi, dneo_phi, dum)
    if (proc0 .and. (ptheta_neo .gt. 1)) then
       deallocate (neo_coefs_ext, theta_neo_ext, neo_phi_ext, dneo_phi_ext)
    endif
  end subroutine read_neocoefs
end module lowflow
