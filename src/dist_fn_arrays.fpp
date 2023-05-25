!> A container for the arrays that are used to store the distribution
!! function among other things.
!!
!! These need to be accessible at a lower dependency level than the
!! dist_fn module itself.  These arrays are allocated in the function
!! dist_fn::allocate_arrays.
module dist_fn_arrays

  implicit none

  private

  public :: g, gnew, g_work, g_restart_tmp, kx_shift, theta0_shift, vpa, vpac
  public :: gexp_1, gexp_2, gexp_3
  public :: vperp2, vpar, aj0, aj1
!AJ
  public :: vpa_gf, vperp2_gf, aj0_gf, aj1_gf
  public :: c_rate, to_g_gs2, from_g_gs2
  public :: g_adjust, check_g_bouncepoints, set_h_zero
  public :: antot, antota, antotp, fieldeq, fieldeqa, fieldeqp

  public :: check_is_zero_in_forbidden_region
  public :: check_are_bouncepoint_values_consistent

  ! dist fn
  complex, dimension (:,:,:), allocatable :: g, gnew, g_restart_tmp
  ! (-ntgrid:ntgrid,2, -g-layout-)

#ifndef SHMEM
  complex, dimension (:,:,:), allocatable :: g_work
#else
  complex, dimension (:,:,:), allocatable, target :: g_work
#endif

#ifndef SHMEM
  complex, dimension (:,:,:), allocatable :: gexp_1, gexp_2, gexp_3
#else
  complex, dimension (:,:,:), pointer, save, contiguous :: gexp_1 => null()
  complex, dimension (:,:,:), allocatable :: gexp_2, gexp_3
#endif

  real, dimension(:), allocatable :: kx_shift, theta0_shift
  ! (naky)

  real, dimension (:,:,:), allocatable :: vpa, vpac, vpar
!AJ
  real, dimension (:,:,:,:,:), allocatable :: vpa_gf
  real, dimension (:,:,:,:), allocatable :: vperp2_gf
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: vperp2, aj0, aj1
  ! (-ntgrid:ntgrid, -g-layout-)

!AJ
  real, dimension (:,:,:,:,:), allocatable :: aj0_gf, aj1_gf

  ! collisional diagnostic of heating rate
  complex, dimension (:,:,:,:,:), allocatable :: c_rate
  ! (-ntgrid:ntgrid,ntheta0,naky,nspecies,2) replicated

  complex, dimension (:,:,:), allocatable :: antot, antota, antotp
  complex, dimension (:,:,:), allocatable :: fieldeq, fieldeqa, fieldeqp

  !> g_adjust transforms between representations of perturbed dist'n func'n.
  !>    <delta_f> = g_wesson J0(Z) - q phi/T F_m  where <> = gyroaverage
  !>        g_gs2 = g_wesson - q phi/T J0(Z) F_m - m v_||^2/T B_||/B J1(Z)/Z F_m
  !> For numerical convenience the GS2 variable g uses the form g_gs2.
  !> g_wesson (see Wesson's book, Tokamaks) is often a more convenient form:
  !>     e.g. for handling collisions, calculating v-space moments in real space.
  !>
  !> To transform gnew from g_gs2 to g_wesson form:
  !>    call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
  !> or
  !>    call g_adjust(gnew,phinew,bparnew, to_g_gs2 = .false.)
  !> or transform from gnew from g_wesson to g_gs2 form:
  !>    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
  !> or
  !>    call g_adjust(gnew,phinew,bparnew, to_g_gs2 = .true. )
  !> CMR, 17/4/2012:
  interface g_adjust
     module procedure g_adjust_floats
     module procedure g_adjust_direction
  end interface g_adjust

  type g_adjust_direction_type
     private
     logical :: to_g_gs2
  end type g_adjust_direction_type

  type(g_adjust_direction_type), parameter :: to_g_gs2 = g_adjust_direction_type(.true.)
  type(g_adjust_direction_type), parameter :: from_g_gs2 = g_adjust_direction_type(.false.)

contains

  !> Transform between g_gs2 and g_wesson with direction indicated
  !> by passed logical [to_g_gs2].
  subroutine g_adjust_direction (g, phi, bpar, direction)
    use run_parameters, only: fphi, fbpar
    implicit none
    complex, dimension (:,:,:), intent (in out) :: g
    complex, dimension (:,:,:), intent (in) :: phi, bpar
    type(g_adjust_direction_type), intent(in) :: direction

    if (direction%to_g_gs2) then
       call g_adjust(g, phi, bpar, -fphi, -fbpar)
    else
       call g_adjust(g, phi, bpar, fphi, fbpar)
    end if
  end subroutine g_adjust_direction

  !> Transform between g_gs2 and g_wesson with direction indicated
  !> by passed floats [facphi] and [facbpar].
  subroutine g_adjust_floats (g, phi, bpar, facphi, facbpar)
    use species, only: spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use le_grids, only: grid_has_trapped_particles, jend
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx, il_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar
    integer :: iglo, ig, ik, it, ie, is, il
    complex :: adj
    logical :: trapped
    logical :: with_bpar, with_phi
    real :: phi_factor, bpar_factor

    trapped = grid_has_trapped_particles()

    with_bpar = facbpar /= 0
    with_phi = facphi /= 0

    phi_factor = 0.0

    if (with_bpar) then
       bpar_factor = 2.0*facbpar
    else
       bpar_factor = 0.0
    end if

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ig, ik, it, il, ie, is, phi_factor, adj) &
    !$OMP SHARED(g_lo, with_phi, spec, nonmaxw_corr, facphi, ntgrid, trapped, &
    !$OMP jend, bpar_factor, vperp2, aj1, bpar, phi, aj0, g) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       if (with_phi) then
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          phi_factor = spec(is)%zt*nonmaxw_corr(ie,is)*facphi
       end if

       ! BD:  bpar == delta B_parallel / B_0(theta) so no extra factor of
       ! 1/bmag is needed here.
       do ig = -ntgrid, ntgrid
          ! Avoid adjusting in the forbidden region
          if ( trapped .and. il > jend(ig)) cycle

          adj = bpar_factor*vperp2(ig,iglo)*aj1(ig,iglo)*bpar(ig,it,ik) &
               + phi_factor*phi(ig,it,ik)*aj0(ig,iglo)

          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine g_adjust_floats

  !> This routine checks trapped particle bounce conditions: 
  !!     g(thetab,1:ik,it,il,ie,is)=g(thetab,2:ik,it,il,ie,is)
  !!  and flags fractional errors exceeding a threshold tolerance, tol.
  !! CMR, 3/10/2013:   
  subroutine check_g_bouncepoints(g, ik,it,il,ie,is,err,tol)
    use theta_grid, only: ntgrid, bmag
    use gs2_layouts, only: g_lo, idx
    use le_grids, only: jend, al, il_is_passing
    use mp, only: mp_abort
    
    implicit none
    integer, intent(in) :: ik, it, il, ie, is
    complex, dimension(-ntgrid:,:,g_lo%llim_proc:), intent(in) :: g
    real, intent(out):: err
    real, optional, intent(in):: tol
    real :: tolerance, dg
    integer :: iglo, ig
    logical :: started

    started=.false.
    if (present(tol)) then
       tolerance=tol 
    else 
       tolerance=1.0e-6
    endif
    iglo=idx(g_lo,ik,it,il,ie,is)
    if (iglo.lt.g_lo%llim_proc .or. iglo.gt.g_lo%ulim_proc .or. il_is_passing(il)) then
       return
    endif
    err=0.0
    do ig=-ntgrid,ntgrid
       ! if at a bounce point, check error on g
       if (il.eq.jend(ig) .and.  al(il)*bmag(ig).gt.0.999999) then
          dg=abs(g(ig,1,iglo)-g(ig,2,iglo))/max(abs(g(ig,1,iglo)),abs(g(ig,2,iglo)))
          if ( dg .gt. tolerance) then
             if (.not. started) then
                write(6,fmt='(T7,"ig",T17,"g(ig,1,iglo)",T43,"g(ig,2,iglo)",T63,"FracBP Error" )')
                started=.true.
             endif
             write(6,fmt='(i8, "  (",1pe11.4,",",e11.4,"),(",e11.4,",",e11.4,"), ", e11.4)') ig, g(ig,1,iglo),g(ig,2,iglo), dg
             err=max(err,dg)
          endif
       endif
    enddo
    write(6,fmt='(t5,"ik",t11,"it",t17,"il",t23,"ie",t29,"is",t33,"MaxFracBP Error")')
    write(6,fmt='(5i6,1pe12.4)')ik,it,il,ie,is,err
    write(6,*) "-----"
  end subroutine check_g_bouncepoints

  !> Assign $g$ a value corresponding to $h$ (g_wesson of CMR's note
  !> in [[dist_fn_arrays:g_adjust]]) of $h = 0$ in order to be
  !> consistent with g_adjust. The correct function call to set h = 0
  !> is call set_h_zero(g, phi, bpar, iglo)
  !>
  !> @note This looks like the g_adjust loop kernel, we might be able
  !> to combine these in some way to avoid duplication.
  subroutine set_h_zero (g, phi, bpar, iglo)
    use species, only: spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use le_grids, only: is_passing_hybrid_electron
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx, il_idx
    use run_parameters, only: fphi, fbpar
    implicit none
    complex, dimension (-ntgrid:, :, g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:, :, :), intent (in) :: phi, bpar
    integer, intent (in) :: iglo
    real :: facphi, facbpar
    integer :: ig, ik, it, ie, is, il
    complex :: adj
    real :: phi_factor, bpar_factor

    facphi = -fphi
    facbpar = -fbpar

    ! Skip non hybrid electrons
    is = is_idx(g_lo,iglo)
    ! Skip zonal modes
    ik = ik_idx(g_lo,iglo)
    ! Skip trapped particles
    il = il_idx(g_lo,iglo)
    if (.not. is_passing_hybrid_electron(is, ik, il)) return

    it = it_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)

    if (facphi /= 0) then
       phi_factor = spec(is)%zt*nonmaxw_corr(ie,is)*facphi
    else
       phi_factor = 0.0
    end if

    if (facbpar /= 0) then
       bpar_factor = 2.0*facbpar
    else
       bpar_factor = 0.0
    end if

    do ig = -ntgrid, ntgrid
       adj = bpar_factor*vperp2(ig,iglo)*aj1(ig,iglo)*bpar(ig,it,ik) &
            + phi_factor*phi(ig,it,ik)*aj0(ig,iglo)
       g(ig, 1, iglo) = adj
       g(ig, 2, iglo) = adj
    end do

  end subroutine set_h_zero

  !> Check if the passed distribution function shaped array is non-zero
  !> in the forbidden region. If it is then trigger an abort.
  !>
  !> This is primarily intended as a utility to aid debugging.
  subroutine check_is_zero_in_forbidden_region(g_array)
    use mp, only: mp_abort, max_allreduce
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, il_idx
    use le_grids, only: forbid, il_is_passing
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_array
    integer :: iglo, il, ig
    real :: violation_flag
    logical :: violated

    violated = .false.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ! Don't check for passing particles
       if (il_is_passing(il)) cycle

       do ig = -ntgrid, ntgrid
          if (.not. forbid(ig, il)) cycle
          violated = violated .or. &
               any(abs(g_array(ig, :, iglo)) > 0.0)
       end do
       ! If we've found a violation don't bother checking any other points
       if (violated) exit
    end do

    ! Use a real to represent the outcome so we can use existing mp routines
    violation_flag = 0.0
    if (violated) violation_flag = 1.0
    call max_allreduce(violation_flag)
    if (violation_flag > 0.0) then
       call mp_abort('The passed array is non-zero in the forbidden region', .true.)
    end if
  end subroutine check_is_zero_in_forbidden_region

  !> Check if the passed distribution function shaped array has the same
  !> value at bounce points for both signs of v||. This provides similar
  !> functionality to check_g_bouncepoints but is intended to offer a
  !> quieter approach which could be used more routinely (i.e. throughout
  !> a run).
  !>
  !> This is primarily intended as a utility to aid debugging.
  subroutine check_are_bouncepoint_values_consistent(g_array, tolerance)
    use mp, only: mp_abort, max_allreduce
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, il_idx
    use le_grids, only: il_is_passing, is_bounce_point
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_array
    real, intent(in), optional :: tolerance
    integer :: iglo, il, ig
    real :: violation_flag, local_tolerance
    logical :: violated

    local_tolerance = get_option_with_default(tolerance, epsilon(0.0) * 1.0e2)
    violated = .false.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ! Don't check for passing particles
       if (il_is_passing(il)) cycle

       do ig = -ntgrid, ntgrid
          if (.not. is_bounce_point(ig, il)) cycle
          violated = violated .or. &
               (abs(g_array(ig, 1, iglo) - g_array(ig, 2, iglo)) > local_tolerance)
       end do
       ! If we've found a violation don't bother checking any other points
       if (violated) exit
    end do

    ! Use a real to represent the outcome so we can use existing mp routines
    violation_flag = 0.0
    if (violated) violation_flag = 1.0
    call max_allreduce(violation_flag)
    if (violation_flag > 0.0) then
       call mp_abort('The passed array has inconsistent bounce points.', .true.)
    end if
  end subroutine check_are_bouncepoint_values_consistent
end module dist_fn_arrays
