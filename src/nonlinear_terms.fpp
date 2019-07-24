!> FIXME : Add documentation
module nonlinear_terms
  !
  ! missing factors of B_a/B(theta) in A_perp terms??
  !
  implicit none

  private

  public :: init_nonlinear_terms, finish_nonlinear_terms
  public :: read_parameters, wnml_nonlinear_terms, check_nonlinear_terms
  public :: add_explicit_terms, nonlinear_mode_switch
  public :: finish_init, reset_init, nonlin, accelerated
  public :: cfl
  public :: nonlinear_mode_none, nonlinear_mode_on
  public :: gryfx_zonal
  public :: cfl_req_hand, nb_check_time_step_too_large
  public :: nonlinear_terms_testing

  !> Public type allowing tests to access otherwise-private routines
  type :: nonlinear_terms_testing
    contains
      procedure, nopass :: add_nl

  end type nonlinear_terms_testing

  !> FIXME : Add documentation  
  type gs2_gryfx_zonal_type
    logical :: on = .false.
    complex*8, dimension (:), pointer :: NLdens_ky0, NLupar_ky0, NLtpar_ky0, &
                                         NLtprp_ky0, NLqpar_ky0, NLqprp_ky0
    logical :: first_half_step = .false.
  end type gs2_gryfx_zonal_type

  type(gs2_gryfx_zonal_type) :: gryfx_zonal

  integer :: istep_last = 0

  ! knobs
  integer :: nonlinear_mode_switch

  integer, parameter :: nonlinear_mode_none = 1, nonlinear_mode_on = 2

#ifndef SHMEM
  real, dimension (:,:), allocatable :: ba, gb, bracket
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, dimension (:,:,:), allocatable :: aba, agb, abracket
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc

#else
  real, save, dimension (:,:), pointer :: ba => null(), gb => null(), bracket => null()
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, save, dimension (:,:,:), pointer, contiguous :: aba => null(), &
       agb => null(), abracket => null()
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc
#endif

! CFL coefficients
  real :: cfl, cflx, cfly
  !< Maximum velocity of the nonlinear term, the reciprocal of the maximum timestep that
  !< satisfies the CFL condition. This currently has to be module level due to potential
  !< non-blocking reduction.
  real :: max_vel
  integer :: cfl_req_hand

! for gryfx
  real :: densfac, uparfac, tparfac, tprpfac, qparfac, qprpfac

  logical :: nonlin = .false.
  logical :: initialized = .false.
  logical :: initializing = .true.
  logical :: alloc = .true.
  logical :: zip = .false.
  logical :: nl_forbid_force_zero = .true.
  logical :: accelerated = .false.

  logical :: exist = .false.
  
contains

  !> FIXME : Add documentation
  subroutine check_nonlinear_terms(report_unit,delt_adj)
    use gs2_time, only: code_dt_min
    use kt_grids, only: box
    use run_parameters, only: code_delt_max, nstep, wstar_units
    use theta_grid, only: nperiod
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: delt_adj
    if (nonlin) then
       write (report_unit, *) 
       write (report_unit, fmt="('This is a nonlinear simulation.')")
       write (report_unit, *) 
       if (wstar_units) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs require wstar_units = .false. in the knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if ( .not. box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
          write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (nperiod > 1) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs usually have nperiod = 1.')") 
          write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       write (report_unit, fmt="('The minimum delt ( code_dt_min ) = ',e11.4)") code_dt_min
       write (report_unit, fmt="('The maximum delt (code_delt_max) = ',e11.4)") code_delt_max
       write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp). (cfl)')") cfl
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
    endif
  end subroutine check_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine wnml_nonlinear_terms(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    if (nonlinear_mode_switch == nonlinear_mode_on) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e17.10)") cfl
       write (unit, fmt="(' nl_forbid_force_zero = ',L1)") nl_forbid_force_zero
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine init_nonlinear_terms 
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, nx, ny, akx, aky
    use le_grids, only: init_le_grids, nlambda, negrid
    use species, only: init_species, nspec
    use gs2_layouts, only: init_dist_fn_layouts, yxf_lo, accelx_lo
    use gs2_layouts, only: init_gs2_layouts
    use gs2_transforms, only: init_transforms
#ifdef SHMEM
    use mp, only : iproc
    use shm_mpi3, only : shm_alloc
#endif
    implicit none
    logical :: dum1, dum2
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    if (debug) write(6,*) "init_nonlinear_terms: init_gs2_layouts"
    call init_gs2_layouts
    if (debug) write(6,*) "init_nonlinear_terms: init_theta_grid"
    call init_theta_grid
    if (debug) write(6,*) "init_nonlinear_terms: init_kt_grids"
    call init_kt_grids
    if (debug) write(6,*) "init_nonlinear_terms: init_le_grids"
    call init_le_grids (dum1, dum2)
    if (debug) write(6,*) "init_nonlinear_terms: init_species"
    call init_species
    if (debug) write(6,*) "init_nonlinear_terms: init_dist_fn_layouts"
    call init_dist_fn_layouts (naky, ntheta0, nlambda, negrid, nspec)
    
    call read_parameters

    if (debug) write(6,*) "init_nonlinear_terms: init_transforms"
    if (nonlinear_mode_switch /= nonlinear_mode_none) then
       call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

       if (debug) write(6,*) "init_nonlinear_terms: allocations"
       if (alloc) then
          if (accelerated) then
#ifndef SHMEM
             allocate (     aba(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (     agb(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (abracket(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
#else
             call shm_alloc(aba, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
             call shm_alloc(agb, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
             call shm_alloc(abracket, (/ 1, 2*ntgrid+1, 1, 2, accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
#endif
             aba = 0. ; agb = 0. ; abracket = 0.
          else
#ifndef SHMEM
             allocate (     ba(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (     gb(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (bracket(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
#else
             call shm_alloc(ba, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
             call shm_alloc(gb, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
             call shm_alloc(bracket, (/ 1,yxf_lo%ny,yxf_lo%llim_proc,yxf_lo%ulim_alloc /))
#endif
             ba = 0. ; gb = 0. ; bracket = 0.
          end if
          alloc = .false.
       end if

       cfly = (aky(naky)/cfl)*0.5
       cflx = (akx((ntheta0+1)/2)/cfl)*0.5
    end if

  end subroutine init_nonlinear_terms

  !> FIXME : Add documentation  
  subroutine read_parameters
    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (4), parameter :: nonlinearopts = &
         (/ text_option('default', nonlinear_mode_none), &
            text_option('none', nonlinear_mode_none), &
            text_option('off', nonlinear_mode_none), &
            text_option('on', nonlinear_mode_on) /)
    character(20) :: nonlinear_mode
    namelist /nonlinear_terms_knobs/ nonlinear_mode, cfl, &
         zip, nl_forbid_force_zero, &
         densfac, uparfac, tparfac, tprpfac, qparfac, qprpfac
    integer :: ierr, in_file

    if (proc0) then
       nonlinear_mode = 'default'
       cfl = 0.1

       densfac = 1.
       uparfac = 1./sqrt(2.)
       tparfac = 1./2.
       tprpfac = 1./2. 
       qparfac = 1./sqrt(8.)
       qprpfac = 1./sqrt(8.)

       in_file=input_unit_exist("nonlinear_terms_knobs",exist)
       if(exist) read (unit=in_file,nml=nonlinear_terms_knobs)

       ierr = error_unit()
       call get_option_value &
            (nonlinear_mode, nonlinearopts, nonlinear_mode_switch, &
            ierr, "nonlinear_mode in nonlinear_terms_knobs",.true.)
    end if

    call broadcast (nonlinear_mode_switch)
    call broadcast (cfl)
    call broadcast (nl_forbid_force_zero)
    call broadcast (zip)
    call broadcast (densfac)
    call broadcast (uparfac)
    call broadcast (tparfac)
    call broadcast (tprpfac)
    call broadcast (qparfac)
    call broadcast (qprpfac)

    istep_last = 0

    if (nonlinear_mode_switch == nonlinear_mode_on)  then
       nonlin = .true.
    end if

  end subroutine read_parameters

  !> FIXME : Add documentation  
  subroutine add_explicit_terms (g1, g2, g3, phi, apar, bpar, istep, bd)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: save_dt_cfl
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    real :: dt_cfl
    logical, parameter :: nl = .true.

    select case (nonlinear_mode_switch)
    case (nonlinear_mode_none)
!!! NEED TO DO SOMETHING HERE...  BD GGH
       dt_cfl = 1.e8
       call save_dt_cfl (dt_cfl)
#ifdef LOWFLOW
       if (istep /=0) &
            call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd)
#endif
    case (nonlinear_mode_on)
       if (istep /= 0) then
         call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, nl)
       endif
    end select
  end subroutine add_explicit_terms

  !> FIXME : Add documentation  
  subroutine add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd,  nl)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use dist_fn_arrays, only: g
    use gs2_time, only: save_dt_cfl
    use run_parameters, only: reset
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    logical, intent (in), optional :: nl
    integer, parameter :: verb = 4
    real, parameter :: zero = epsilon(0.0)
    integer :: iglo, ik, it
    real :: dt_cfl

    call debug_message(verb, 'nonlinear_terms::add_explicit starting')

    if (initializing) then
       if (present(nl)) then
          dt_cfl = 1.e8
          call save_dt_cfl (dt_cfl)
       end if
       return
    endif


    !
    ! @todo Currently not self-starting.  Need to fix this.
    !

    call debug_message(verb, 'nonlinear_terms::add_explicit copying old terms')

    ! This ensures that the nonlinear term is calculated only once per timestep.
    ! It guards against calculating the nonlinear term in two cases:
    !   1. the second call to timeadv in a usual timestep
    !   2. both calls to timeadv after the CFL condition check has triggered a restart
    if (istep /= istep_last) then

       istep_last = istep

       ! Shuffle time history along -- this is a place pointers might help
       ! avoid some copying.
       ! Should we only need to do this if we have a nonlinear simulation?
       ! It's possible we might have other explicit sources, in which case we
       ! probably always need to do this shuffling, however in that case this
       ! subroutine probably best belongs in a different module.
       g3 = g2
       g2 = g1

       call debug_message(verb, 'nonlinear_terms::add_explicit copied old terms')

       ! if running nonlinearly, then compute the nonlinear term at grid points
       ! and store it in g1
       if (present(nl)) then
          call debug_message(verb, 'nonlinear_terms::add_explicit calling add_nl')
          if(gryfx_zonal%on) then
             call add_nl_gryfx (g1)
          else
             call add_nl (g1, phi, apar, bpar)
          endif
          ! takes g1 at grid points and returns 2*g1 at cell centers

          call nl_center (g1, bd)

          if(reset) then
             return !Return if resetting
          endif

       else
          g1 = 0.
       end if

#ifdef LOWFLOW
       ! do something
#endif

    end if

    if (zip) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          if (ik /= 1) then
             g (:,1,iglo) = 0.
             g (:,2,iglo) = 0.
             g1(:,1,iglo) = 0.
             g1(:,2,iglo) = 0.
          end if
       end do
    end if

    istep_last = istep
  end subroutine add_explicit
 
  !> Takes input array evaluated at theta grid points
  !! and overwrites it with array evaluated at cell centers
  !! note that there is an extra factor of 2 in output array
  subroutine nl_center (gtmp, bd)
    !
    !CMR, 13/10/2014
    ! Fixing some (cosmetic) issues from prior to R3021.
    ! (1) forbid(-ntgrid:ntgrid) refers to grid points at cell boundaries
    !     gtmp(-ntgrid:ntgrid-1) refers to cell centers
    !     => applying forbid to output gtmp did not zero the correct things!!!
    ! (2) totally trapped particles are special, so return source without upwinding is fine
    !     BUT source for other trapped particles should NOT enter forbidden region.
    !     if ig or ig+1 forbidden (i.e. ig+1 or ig is bounce point) now return gtmp(ig,:,iglo)=0
    !
    use dist_fn_arrays, only: ittp
    use gs2_layouts, only: g_lo, il_idx
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, ng2, jend
    use mp, only: mp_abort

    implicit none
    real, intent (in) :: bd    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: gtmp
    integer :: iglo, il, ig
    logical :: trapped

    trapped = .false.

    if (minval(jend) .gt. ng2) trapped=.true.

    ! factor of one-half appears elsewhere
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)

       !
       !CMR, 7/10/2014: 
       ! Incoming gtmp SHOULD vanish in forbidden region, 
       ! New logical in nonlinear_terms_knobs namelist: nl_forbid_force_zero
       !     nl_forbid_force_zero =.t. : force zeroing    (default)
       !     nl_forbid_force_zero =.f. : NO forced zeroing
       !                                 ie assume forbidden gtmp is zero on entry 
       !

       if ( nl_forbid_force_zero ) then
          ! force spurious gtmp outside trapped boundary to be zero
          where (forbid(:,il))
             gtmp(:,1,iglo) = 0.0
             gtmp(:,2,iglo) = 0.0
          end where
       endif

       do ig = -ntgrid, ntgrid-1
          !
          !CMR, 7/10/2014: 
          ! loop sets gtmp to value at upwinded cell center RIGHT of theta(ig)
          !           except for ttp where upwinding makes no sense!
          !
          if (il >= ittp(ig)) cycle
          if ( trapped .and. ( il > jend(ig) .or. il > jend(ig+1)) ) then
             !
             !CMR, 7/10/2014: 
             !   if either ig or ig+1 is forbidden, no source possible in a cell RIGHT of theta(ig) 
             !   => gtmp(ig,1:2,iglo)=0
             !
             gtmp(ig,1:2,iglo) = 0.0
          else 
             !
             !CMR, 7/10/2014: 
             !    otherwise ig and ig+1 BOTH allowed, and upwinding in cell RIGHT of theta(ig) is fine
             !
             gtmp(ig,1,iglo) = (1.+bd)*gtmp(ig+1,1,iglo) + (1.-bd)*gtmp(ig,1,iglo)
             gtmp(ig,2,iglo) = (1.-bd)*gtmp(ig+1,2,iglo) + (1.+bd)*gtmp(ig,2,iglo)

          endif
       end do
    end do
  end subroutine nl_center

  !> Calculate the nonlinear term and corresponding CFL condition
  subroutine add_nl (g1, phi, apar, bpar)
    use mp, only: max_allreduce, nb_max_allreduce
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: g, g_adjust
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi, reset, immediate_reset
    use kt_grids, only: aky, akx
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use constants, only: zi
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer :: i, j, k
    real :: dt_cfl
    integer, parameter :: verb = 4

    integer :: iglo, ik, it, is, ig, ia
    
    call debug_message(verb, 'nonlinear_terms::add_nl starting')

    !Form g1=i*kx*chi
    call load_kx_phi(g1, phi)    
    call load_kx_bpar(g1, bpar)    
    call load_kx_apar(g1, apar)    

    call debug_message(verb, 'nonlinear_terms::add_nl calling transform2')
    !Transform to real space
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx')
    
    !Form g1=i*ky*g_wesson
    g1=g
    call g_adjust(g1,phi,bpar,fphi,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*aky(ik)
    enddo

    !Transform to real space    
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dg/dy')

    max_vel = 0.
    call get_max_vel(cfly)

    call calculate_bracket(.true.)
    
    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx.dg/dy')

    !Form g1=i*ky*chi
    call load_ky_phi(g1, phi)
    call load_ky_bpar(g1, bpar)
    call load_ky_apar(g1, apar)

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dy')

    !Transform to real space    
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if

    !Form g1=i*kx*g_wesson
    g1=g
    call g_adjust(g1,phi,bpar,fphi,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*akx(it)
    enddo
    
    !Transform to real space
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dg/dx')
    
    call get_max_vel(cflx)

    call calculate_bracket(.false.)
    
    call debug_message(verb, &
      'nonlinear_terms::add_nl calculated dphi/dx.dg/dy.dphi/dy.dg/dx')

    !Transform NL source back to spectral space
    if (accelerated) then
       call inverse2 (abracket, g1, ia)
    else
       call inverse2 (bracket, g1)
    end if

    !Estimate the global cfl limit based on max_vel
    !If immediate_reset, do a max allreduce.
    !If not, overlap comms and work, and check for
    ! reset at the end of the timestep.
    if(immediate_reset)then
       call max_allreduce(max_vel)
       dt_cfl = 1./max_vel
       call save_dt_cfl (dt_cfl)

       !Now check to see if we've violated the
       !cfl condition if requested
       call check_time_step_too_large(reset)
       !If we have violated cfl then return immediately
       if (reset) return
    else
       !Do a nonblocking reduce and check later
       call nb_max_allreduce(max_vel, cfl_req_hand)
    endif

    call debug_message(verb, &
      'nonlinear_terms::add_nl calculated nonlinear term')

  end subroutine add_nl

  !> FIXME : Add documentation    
  subroutine load_kx_phi(g1, phi)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: aj0
    use run_parameters, only: fphi
    use kt_grids, only: akx
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi      
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer :: iglo, it, ik

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          fac = zi*akx(it)*aj0(:,iglo)*phi(:,it,ik)*fphi
          g1(:,1,iglo) = fac
          g1(:,2,iglo) = fac
       end do
    else
       g1 = 0.
    endif
  end subroutine load_kx_phi

  !> FIXME : Add documentation    
  subroutine load_ky_phi(g1, phi)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: aj0
    use run_parameters, only: fphi
    use kt_grids, only: aky
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi      
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer :: iglo, it, ik

    if (fphi > epsilon(0.0)) then      
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          fac = zi*aky(ik)*aj0(:,iglo)*phi(:,it,ik)*fphi
          g1(:,1,iglo) = fac
          g1(:,2,iglo) = fac
       end do
    else
       g1 = 0.
    endif

  end subroutine load_ky_phi

  ! should I use vpa or vpac in next two routines??

  !> FIXME : Add documentation    
  subroutine load_kx_apar(g1, apar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vpa, aj0
    use species, only: spec
    use run_parameters, only: fapar
    use kt_grids, only: akx
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: apar
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if(fapar <= epsilon(0.0)) return

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       fac = zi*akx(it)*aj0(:,iglo)*spec(is)%stm*apar(:,it,ik)*fapar
       g1(:,1,iglo) = g1(:,1,iglo) - fac*vpa(:,1,iglo)
       g1(:,2,iglo) = g1(:,2,iglo) - fac*vpa(:,2,iglo)
    end do

  end subroutine load_kx_apar

  !> FIXME : Add documentation    
  subroutine load_ky_apar(g1, apar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vpa, aj0
    use species, only: spec
    use run_parameters, only: fapar
    use kt_grids, only: aky
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: apar            
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if(fapar <= epsilon(0.0)) return

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       fac = zi*aky(ik)*aj0(:,iglo)*spec(is)%stm*apar(:,it,ik)*fapar
       g1(:,1,iglo) = g1(:,1,iglo) - fac*vpa(:,1,iglo)
       g1(:,2,iglo) = g1(:,2,iglo) - fac*vpa(:,2,iglo)
    end do

  end subroutine load_ky_apar

  !> FIXME : Add documentation    
  subroutine load_kx_bpar(g1, bpar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vperp2, aj1
    use species, only: spec
    use run_parameters, only: fbpar
    use kt_grids, only: akx
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: bpar      
    complex, dimension(-ntgrid:ntgrid) :: fac    
    integer iglo, it, ik, is

    if (fbpar <= epsilon(0.0)) return

    ! Is this factor of two from the old normalization?

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       fac = zi*akx(it)*aj1(:,iglo) &
            *2.0*vperp2(:,iglo)*spec(is)%tz*bpar(:,it,ik)*fbpar
       g1(:,1,iglo) = g1(:,1,iglo) + fac
       g1(:,2,iglo) = g1(:,2,iglo) + fac
    end do

  end subroutine load_kx_bpar

  !> FIXME : Add documentation    
  subroutine load_ky_bpar(g1, bpar)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: vperp2, aj1
    use species, only: spec
    use run_parameters, only: fbpar
    use kt_grids, only: aky
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: bpar
    complex, dimension(-ntgrid:ntgrid) :: fac
    integer iglo, it, ik, is

    if (fbpar <= epsilon(0.0)) return

    ! Is this factor of two from the old normalization?

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       fac =  zi*aky(ik)*aj1(:,iglo) &
            *2.0*vperp2(:,iglo)*spec(is)%tz*bpar(:,it,ik)*fbpar
       g1(:,1,iglo) = g1(:,1,iglo) + fac
       g1(:,2,iglo) = g1(:,2,iglo) + fac
    end do

  end subroutine load_ky_bpar

  !> Calculates an estimate of the current maximum velocity due to the
  !! perturbed potentials. This is currenly a processor local estimate.
  subroutine get_max_vel(scaling)
    implicit none
    real, intent(in) :: scaling
    real :: local_max_vel
    if (accelerated) then
       local_max_vel = maxval(abs(aba))
    else
       local_max_vel = maxval(abs(ba))
    end if
    
    local_max_vel = local_max_vel * scaling

    max_vel = max(max_vel, local_max_vel)
  end subroutine get_max_vel

  !> Calculate (d Chi /dx).(d g_wesson/dy) and store in bracket if is_first_term = .true.
  !! else calculate (d Chi /dy).(d g_wesson/dx) and subtract from bracket
  subroutine calculate_bracket(is_first_term)
    use theta_grid, only: kxfac
    implicit none
    logical :: is_first_term

    if (is_first_term) then
       if (accelerated) then
          abracket = aba*agb*kxfac
       else
          bracket = ba*gb*kxfac
       endif
    else
       if (accelerated) then
          abracket = abracket - aba*agb*kxfac
       else
          bracket = bracket - ba*gb*kxfac
       endif
    endif
  end subroutine calculate_bracket
  
  !> this subroutine constructs the GK nonlinear term from the 6 GryfX
  !! gyrofluid nonlinear terms.  we first construct in gryfx units.
  !! gryfx uses vt=sqrt(T/m) normalization, 
  !! so vt_gs2 = sqrt(2T/m) = sqrt(2) vt_gryfx.
  !! since vpa and vperp2 are normalized to vt_gs2, below we multiply by
  !! factors of sqrt(2) to get vpa and vperp2 in gryfx units, i.e.
  !! vpa_gryfx = sqrt(2) vpa
  !! vperp2_gryfx = 2 vperp2
  !!
  !! the Hermite-Laguerre construction of the distribution from the gyrofluid
  !! moments is of the form
  !! df/F_M = n + vpar/vt upar + 1/2(vpar^2/vt^2-1) tpar 
  !! + (vprp^2/2vt^2-1) tprp + 1/6(vpar^3/3vt^3 - 3 vpar/vt) qpar
  !! + vpar/vt (vprp^2/2vt^2-1) qprp
  !! where vt=vt_gryfx, and moments are normalized in gryfx units
  subroutine add_nl_gryfx (g1)
    use mp, only: max_allreduce
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use dist_fn_arrays, only: g_adjust, vpa, vperp2
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi, immediate_reset
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use mp, only: broadcast
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1

    integer :: iglo, ik, it, ig, iz, is, isgn, index_gryfx

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             iz = ig + ntgrid + 1
             if(ig==ntgrid) iz = 1 ! periodic point not included in gryfx arrays
             index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                  g_lo%naky*g_lo%ntheta0*(iz-1) + &
                  (2*ntgrid)*g_lo%naky*g_lo%ntheta0*(is-1)
             
             g1(ig,isgn,iglo) =  gryfx_zonal%NLdens_ky0(index_gryfx) + &
                  
                  0.5*(2.*vpa(ig,isgn,iglo)*vpa(ig,isgn,iglo) - 1.)* &
                  gryfx_zonal%NLtpar_ky0(index_gryfx) + &
                  
                  (vperp2(ig,iglo) - 1.0)* &
                  gryfx_zonal%NLtprp_ky0(index_gryfx) + &

                  sqrt(2.)*vpa(ig,isgn,iglo)* &
                     gryfx_zonal%NLupar_ky0(index_gryfx) + &
                     
                     1./6.*( (sqrt(2.)*vpa(ig,isgn,iglo))**3. &
                     - 3*sqrt(2.)*vpa(ig,isgn,iglo) )* &
                     gryfx_zonal%NLqpar_ky0(index_gryfx) + &
                     
                     sqrt(2.)*vpa(ig,isgn,iglo)*(vperp2(ig,iglo) - 1.)* &
                     gryfx_zonal%NLqprp_ky0(index_gryfx) 
             
             ! now we must account for the fact that phi and g are normalized to a/rho_i,
             ! and grad ~ k is normalized to rho_i, and rho_gs2 = sqrt2 rho_gryfx.
             ! we have just calculated the GK nonlinear term in gryfx units:
             ! g1 = z_hat X grad_gryfx phi_gryfx . grad_gryfx g_gryfx =
             ! z_hat X (1/sqrt2 grad_gs2)(sqrt2 phi_gs2) . (1/sqrt2 grad_gs2)(sqrt2 g_gs2)
             ! = z_hat X grad_gs2 phi_gs2 . grad_gs2 g_gs2
             ! so it turns out that the GK nonlinear term is the same in gryfx units
             ! and gs2 units, so we don't need any additional factors.

          end do
       end do
    end do
    g1 = -g1 !left-handed / right-handed conversion    
  end subroutine add_nl_gryfx

  !> Nonblocking version of check_time_step_too_large
  !!
  !! This receives the mpi_iallreduce of the maximum velocity and calculates
  !! the CFL condition.
  !!
  !! Joseph Parker, STFC
  subroutine nb_check_time_step_too_large
    use mp, only: iproc
    use run_parameters, only: reset, immediate_reset
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use mp, only: wait
    implicit none
    real :: dt_cfl

    !Estimate the global cfl limit based on max_vel
    !If immediate_reset, this check has already been done during the timestep.
    if(immediate_reset) return

    !write (*,*) 'Receiving max_iallreduce', iproc, immediate_reset, max_vel
    call wait(cfl_req_hand)    
    !write (*,*) 'Received max_iallreduce', iproc, immediate_reset, max_vel
    dt_cfl = 1./max_vel
    call save_dt_cfl (dt_cfl)

    !Now check to see if we've violated the 
    !cfl condition if requested
    call check_time_step_too_large(reset)
  
  end subroutine nb_check_time_step_too_large

  !> FIXME : Add documentation
  subroutine reset_init
    implicit none
    initialized = .false.
    initializing = .true.
  end subroutine reset_init

  !> FIXME : Add documentation  
  subroutine finish_init
    implicit none
    initializing = .false.
  end subroutine finish_init

  !> FIXME : Add documentation  
  subroutine finish_nonlinear_terms
#ifndef SHMEM    
    use gs2_transforms, only: finish_transforms
#else
    use shm_mpi3, only : shm_free
#endif
    implicit none

#ifndef SHMEM
    if (allocated(aba)) deallocate (aba, agb, abracket)
    if (allocated(ba)) deallocate (ba, gb, bracket)
#else
    if (associated(aba)) call shm_free(aba)
    if (associated(agb)) call shm_free(agb)
    if (associated(abracket)) call shm_free(abracket)
    if (associated(ba)) then
       call shm_free(ba)
    endif
    if (associated(gb)) then
       call shm_free(gb)
    endif
    if (associated(bracket)) then
       call shm_free(bracket)
    endif
#endif
    nonlin = .false. ; alloc = .true. ; zip = .false. ; accelerated = .false.
    initialized = .false. ; initializing = .true. ; exist = .false.

#ifndef SHMEM
    call finish_transforms
#endif
  end subroutine finish_nonlinear_terms
end module nonlinear_terms


