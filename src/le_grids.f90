!> FIXME : Add documentation
module egrid

! By Tomo Tatsuno, Aug 2005
! Improved accuracy and speed and maximum number of energy grid points
!

  implicit none

  private

  public :: setvgrid, setvgrid_genquad, init_egrid
  public :: zeroes, x0, zeroes_maxw, x0_maxw

  real :: x0_maxw
  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable, save :: zeroes_maxw
  real, dimension(:,:), allocatable, save :: zeroes

contains

  !> FIXME : Add documentation  
  subroutine init_egrid (negrid)
    use species, only: nspec
    implicit none
    integer, intent (in) :: negrid

    if (.not. allocated(zeroes)) then
       allocate (zeroes(negrid-1, nspec)) ; zeroes = 0.
       allocate (zeroes_maxw(negrid-1)) ; zeroes_maxw = 0.
       allocate (x0(nspec))
    end if
  end subroutine init_egrid
  
  !> FIXME : Add documentation
  subroutine setvgrid (vcut_in, negrid, epts, wgts, nesub_in)
    use constants, only: pi => dpi
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, f0_sdanalytic, f0_maxwellian, calculate_f0_arrays, f0_values

    implicit none
    integer, intent (in) :: negrid
    real, intent (in) :: vcut_in
    integer, intent (in) :: nesub_in
    real, dimension(:,:), intent (out) :: epts, wgts
    integer:: nesub, is
    real:: vcut

    call init_egrid (negrid)

    do is = 1,nspec 
       select case (spec(is)%f0type)
       case(f0_maxwellian)
          vcut = vcut_in
          nesub = nesub_in
       case(f0_sdanalytic)
          vcut = 1.0
          nesub = negrid
       case default 
          vcut = vcut_in 
          nesub = nesub_in
       end select

       ! get grid points in v up to vcut (epts is not E yet)
       call get_legendre_grids_from_cheb (0., vcut, epts(:nesub,is), wgts(:nesub,is))

       ! change from v to E
       epts(:nesub,is) = epts(:nesub,is)**2

       ! absorb exponential and volume element in weights
       !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
       !> GW: F0 now gets multiplied later
       wgts(:nesub, is) = wgts(:nesub, is)*epts(:nesub,is)*pi

       if (negrid > nesub) then

          ! get grid points in y = E - vcut**2 (epts not E yet)
          call get_laguerre_grids (epts(nesub+1:,is), wgts(nesub+1:,is))

          ! change from y to E
          epts(nesub+1:,is) = epts(nesub+1:,is) + vcut**2

          ! Do not absorb exponential in weights, since we'll be multiplying by f0 later
          wgts(nesub+1:, is) = wgts(nesub+1:, is)*exp(epts(nesub+1:,is)-vcut**2)*pi*0.5*sqrt(epts(nesub+1:,is))


       end if
    end do

    call calculate_f0_arrays(epts)
    wgts = wgts * f0_values
    
    zeroes(:,:) = sqrt(epts(:negrid-1,:))
    x0(:) = sqrt(epts(negrid,:))

  end subroutine setvgrid

  !> FIXME : Add documentation  
  subroutine setvgrid_genquad (negrid, epts, wgts, nesub)
    use genquad, only: get_quadrature_rule
    use constants, only: pi => dpi
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, calculate_f0_arrays, eval_f0, set_current_f0_species, f0_maxwellian

    implicit none
    
    integer, intent (in) :: negrid
    integer, intent (in) :: nesub
    real, dimension(:,:), intent (out) :: epts
    real, dimension(:,:), intent (out) :: wgts
    integer :: is

    call init_egrid (negrid)

    do is = 1,nspec

      ! Call the function that sets the species number so that eval_f0 takes only one argument
      call set_current_f0_species(is)

      if ( spec(is)%f0type .EQ. f0_maxwellian ) then
         ! If Maxwellian, use an semi-infinite domain
         call get_quadrature_rule(eval_f0,negrid,0.0,1.0,epts(:,is),wgts(:,is),.false.,.true.)
      else 
         ! Otherwise, use a bounded domain with vmax = vref_s.
         call get_quadrature_rule(eval_f0,negrid,0.0,1.0,epts(:,is),wgts(:,is),.false.,.false.)
      end if
    
      ! change from v to E
      epts(:,is) = epts(:,is)**2

      ! absorb exponential and volume element in weights
      !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
      ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
      ! See eq. 4.12 of M. Barnes's thesis
      wgts(:, is) = wgts(:, is)*pi*epts(:,is)

    end do

    call calculate_f0_arrays(epts)

    zeroes(:,:) = sqrt(epts(:negrid-1,:))
    x0(:) = sqrt(epts(negrid,:))

  end subroutine setvgrid_genquad

end module egrid

!> FIXME : Add documentation
module le_grids
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use redistribute, only: redist_type
  implicit none

  private

  public :: init_le_grids, finish_le_grids
  public :: read_parameters, wnml_le_grids
  public :: integrate_species, write_mpdist, write_mpdist_le
  public :: energy, energy_maxw, speed_maxw, al, delal, jend, forbid, dele, wl, w, w_maxw
  public :: negrid, nlambda, ng2, nxi, lmax, integrate_moment, nesub
  public :: xloc, sgn, ixi_to_il, ixi_to_isgn, speed, xi
  public :: xx, nterp, new_trap_int, vcut, genquad
  public :: init_weights, legendre_transform, lagrange_interp, lagrange_coefs
  public :: finish_weights
  public :: eint_error, lint_error, trap_error, wdim
  public :: integrate_kysum, integrate_volume
  public :: get_flux_vs_theta_vs_vpa
  public :: lambda_map, energy_map, g2le, g2gf, init_map, init_g2gf
  public :: integrate_species_sub !<DD>
  public :: hermite_prob

  public :: le_grids_config_type
  public :: set_le_grids_config
  
  interface integrate_moment
     module procedure integrate_moment_c34
     module procedure integrate_moment_lec
     module procedure integrate_moment_le
     module procedure integrate_moment_r33
     module procedure integrate_moment_lz 
     module procedure integrate_moment_e
  end interface

  interface integrate_species
     module procedure integrate_species_master
     module procedure integrate_species_gf_nogather
     module procedure integrate_species_le
     module procedure integrate_species_lz 
     module procedure integrate_species_e
  end interface

  interface integrate_volume
     module procedure integrate_volume_c
     module procedure integrate_volume_r
  end interface

  interface get_hermite_polynomials
     module procedure get_hermite_polynomials_1d
     module procedure get_hermite_polynomials_4d
  end interface

  real, dimension (:), allocatable :: xx ! (ng2)
  real, dimension (:,:), allocatable :: wlerr, xloc ! mbmark
  real, dimension (:,:,:), allocatable :: werr
  real, dimension (:,:,:), allocatable :: wlterr
  real, dimension (:,:,:,:), allocatable :: lgrnge
  real, dimension (:,:,:), allocatable :: wlmod
  real, dimension (:,:,:), allocatable :: wtmod
  real, dimension (:,:,:), allocatable :: wmod
  real, dimension (:,:,:), allocatable :: lpe
  real, dimension (:,:), allocatable, save :: lpl
  real, dimension (:,:,:,:), allocatable, save :: lpt

  real, dimension (:), allocatable :: energy_maxw, w_maxw, speed_maxw 
  real, dimension (:,:), allocatable :: energy, w, speed, dele
  real, dimension (:), allocatable :: al, delal ! (nlambda)
  real, dimension (:,:), allocatable :: wl ! (nlambda,-ntgrid:ntgrid)
  integer, dimension (:), allocatable :: jend ! (-ntgrid:ntgrid)
  logical, dimension (:,:), allocatable :: forbid ! (-ntgrid:ntgrid,nlambda)

  real, dimension (:,:), allocatable :: xi
  integer, dimension (:,:), allocatable :: ixi_to_il, ixi_to_isgn
  integer, dimension (2) :: sgn

  integer :: wdim
  complex, dimension (:,:), allocatable :: integration_work
  ! (-ntgrid:ntgrid, -*- processor-dependent -*-)
  logical :: exist

 ! knobs
  integer :: ngauss, npassing, negrid, nesuper, nesub
  real :: bouncefuzz, vcut
  logical :: genquad

  integer :: nlambda, ng2, lmax, nxi
  logical :: accel_x = .false.
  logical :: accel_v = .false.
  logical :: test = .false.
  logical :: trapped_particles = .true.
  logical :: new_trap_int = .false.
  logical :: radau_gauss_grid =.true. ! default lambda grid changed
  ! The default lambda grid is now gauss-radau, for passing particles up to 
  ! and including the passing-trapped boundary (wfb). The grid gives finite weight to
  ! the class of particles which bounce at theta = +/- pi (wfb) in the passing integral
  ! Note that the old gauss-legendre is used in the case that trapped_particles =.false.
  logical :: intinit = .false.
  logical :: slintinit = .false.
  logical :: lintinit = .false.
  logical :: eintinit = .false.
  logical :: initialized = .false.
  logical :: leinit = .false.
  logical :: gfinit = .false.
  logical :: lzinit = .false.
  logical :: einit = .false.
  logical :: init_weights_init = .false.

  integer :: nmax = 500
  integer :: nterp = 100

  real :: wgt_fac = 10.0

  type (redist_type), save :: lambda_map
  type (redist_type), save :: energy_map
  type (redist_type), save :: g2le
  type (redist_type), save :: g2gf

  integer, dimension(:),allocatable :: recvcnts_intspec,displs_intspec
  integer :: sz_intspec, local_rank_intspec
  
  !> Used to represent the input configuration of le_grids
  type, extends(abstract_config_type) :: le_grids_config_type
     ! namelist : le_grids_knobs
     ! indexed : false
     !> Acts as a small tolerance in deciding if a particular pitch
     !> angle is forbidden at a particular theta grid point. Rather
     !> than considering particles to be forbidden if `1.0 - lambda *
     !> B(theta) < 0` we treat them as forbidden if `1.0 - lambda *
     !> B(theta) < -bouncefuzz`. This provides a small "buffer" to
     !> account for floating round off in the calculation of `lambda`
     !> and `B`.
     real :: bouncefuzz = 10*epsilon(0.0)
     !> If true use generalised quadrature scheme for velocity
     !> integrals and energy grid. See [G. Wilkie
     !> thesis](https://drum.lib.umd.edu/handle/1903/17302).
     logical :: genquad = .false.
     !> Sets the number of energy grid points to use. If specified
     !> then overrides values of `nesuper = min(negrid/10 + 1,
     !> 4)` and `nesub = negrid - nesuper`. If not set then the total
     !> number of energy grid points is just `negrid = nesub +
     !> nesuper`. The energy grid is split into two regions at a point
     !> controlled by `vcut`.
     integer :: negrid = -10
     !> Sets the number of energy grid points below the cutoff.
     integer :: nesub = 8
     !> Sets the number of energy grid points above the cutoff.
     integer :: nesuper = 2
     !> If `true` then use a more accurate integration method for the
     !> trapped pitch angle contribution to velocity integrals. The
     !> default uses a simple, but potentially less accurate, finite
     !> difference method.
     logical :: new_trap_int = .false.
     !> The number of untrapped pitch-angles moving in one direction
     !> along field line.
     integer :: npassing = -1
     !> The number of untrapped pitch-angles moving in one direction
     !> along field line is `2*ngauss` if [[le_grids_knobs::npassing]] is
     !> not set. Note that the number of trapped pitch angles is
     !> directly related to the number of theta grid points (per
     !> \(2\pi\)) and is `ntheta/2 + 1`.
     integer :: ngauss = 5
     !> Influences the minimum number of grid points in an integration
     !> subinterval used in calculating the integration grid weights.
     integer :: nmax = 500
     !> The number of nodes used in calculating lagrange coefficients
     !> is `2*nterp-1`.  This code isn't currently used in GS2 so this
     !> input should not influence anything significant.
     integer :: nterp = 100
     !> The default lambda grid is now gauss-radau, for passing
     !> particles up to and including the passing-trapped boundary
     !> (wfb). The grid gives finite weight to the class of particles
     !> which bounce at theta = +/- pi (wfb) in the passing integral
     !> Note that the old gauss-legendre is used in the case that
     !> trapped_particles =.false.
     !>
     !> @note In reference to above comment -
     !> code doesn't seem to change radau_gauss_grid when
     !> trapped_particles is false, should it?
     logical :: radau_gauss_grid = .true.
     !> If `true` then just does a few tests and writes pitch angle
     !> and energy grids to screen before aborting the simulation.
     logical :: test = .false.
     !> If set to `false` then the lambda grid weighting `wl` is set
     !> to zero for trapped pitch angles. This means that integrals
     !> over velocity space do not include a contribution from trapped
     !> particles which is equivalent to the situation where
     !> `eps<=0.0`.  Trapped particle drifts are not set to zero so
     !> "trapped" particles still enter the source term through
     !> `wdfac`. At least for s-alpha the drifts are the main
     !> difference (*please correct if not true*) between the
     !> `eps<=0.0` and the `trapped_particles = .false.` cases as
     !> `Bmag` is not a function of theta in the `eps<=0.0` case
     !> whilst it is in the `trapped_particles = .false.` case.
     logical :: trapped_particles = .true.
     !> No. of standard deviations from the standard Maxwellian beyond
     !> which the distribution function will be set to 0.
     real :: vcut = 2.5
     !> Influences the maximum integration weight allowed.
     real :: wgt_fac = 10.0
   contains
     procedure, public :: read => read_le_grids_config
     procedure, public :: write => write_le_grids_config
     procedure, public :: reset => reset_le_grids_config
     procedure, public :: broadcast => broadcast_le_grids_config
     procedure, public, nopass :: get_default_name => get_default_name_le_grids_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_le_grids_config
  end type le_grids_config_type
  
  type(le_grids_config_type) :: le_grids_config
  
contains

  !> FIXME : Add documentation  
  subroutine wnml_le_grids(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "le_grids_knobs"
    write (unit, fmt="(' nesub = ',i4)") nesub
    write (unit, fmt="(' nesuper = ',i4)") nesuper
    write (unit, fmt="(' ngauss = ',i4)") ngauss
    write (unit, fmt="(' vcut = ',e17.10)") vcut
    write (unit, fmt="(' /')")
  end subroutine wnml_le_grids

  !> FIXME : Add documentation  
  subroutine init_le_grids (accelerated_x, accelerated_v, le_grid_config_in)
    use mp, only: proc0, finish_mp
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
    use egrid, only: setvgrid
    implicit none
    logical, intent (out) :: accelerated_x, accelerated_v
    type(le_grids_config_type), intent(in), optional :: le_grid_config_in    
    integer :: il, ie

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_species
    call init_theta_grid
    call init_kt_grids

    call read_parameters(le_grid_config_in)
    call set_vgrid
    if (proc0) then
       call set_grids
    end if
    call broadcast_results
    call init_integrations

    accelerated_x = accel_x
    accelerated_v = accel_v

    if (test) then
       if (proc0) then
          do il = 1, nlambda
             write(*,*) al(il)
          end do
          write(*,*) 
          do ie = 1, negrid
             write(*,*) energy(ie,:)
          end do
       end if
       call finish_mp
       stop
    endif
    
  end subroutine init_le_grids

  !> FIXME : Add documentation  
  subroutine set_vgrid
    use species, only: nspec, init_species
    use egrid, only: setvgrid, setvgrid_genquad
    call init_species

    allocate (energy(negrid,nspec), w(negrid,nspec))
    if (genquad) then
       call setvgrid_genquad (negrid, energy, w, nesub)
    else
       call setvgrid (vcut, negrid, energy, w, nesub)
    end if
  end subroutine set_vgrid

  !> FIXME : Add documentation 
  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use species, only: nspec
    use egrid, only: init_egrid, zeroes_maxw, x0_maxw
    use theta_grid, only: ntgrid

    implicit none
    integer :: tsize

    tsize = 2*nterp-1

    call broadcast (lmax)
    call broadcast (ng2)
    call broadcast (nlambda)
    call broadcast (nxi)

    if (.not. proc0) then
       allocate(speed(negrid,nspec))
       allocate(w_maxw(negrid))
       allocate(energy_maxw(negrid))
       allocate(speed_maxw(negrid))
       allocate(dele(negrid,nspec))
       allocate (al(nlambda), delal(nlambda))
       allocate (wl(-ntgrid:ntgrid,nlambda))
       allocate (jend(-ntgrid:ntgrid))
       allocate (forbid(-ntgrid:ntgrid,nlambda))
       allocate (xx(ng2))
       allocate (lgrnge(-ntgrid:ntgrid,nlambda,tsize,2))
       allocate (xloc(-ntgrid:ntgrid,tsize))
       allocate (xi(-ntgrid:ntgrid, 2*nlambda))
       allocate (ixi_to_il(-ntgrid:ntgrid, 2*nlambda))
       allocate (ixi_to_isgn(-ntgrid:ntgrid, 2*nlambda))
    end if

    call init_egrid (negrid)
    call broadcast (xx)
    call broadcast (x0_maxw)
    call broadcast (zeroes_maxw)

    call broadcast (al)
    call broadcast (delal)
    call broadcast (jend)
 
    call broadcast (energy_maxw)
    call broadcast (speed)
    call broadcast (speed_maxw)
    call broadcast (dele)
    call broadcast (w_maxw)

    call broadcast(wl)
    call broadcast(forbid)

    call broadcast(xloc)
    call broadcast(lgrnge)

    call broadcast(xi)
    call broadcast(ixi_to_il)
    call broadcast(ixi_to_isgn)

    call broadcast (sgn)

  end subroutine broadcast_results

  !> FIXME : Add documentation  
  subroutine read_parameters(le_grid_config_in)
    implicit none
    type(le_grids_config_type), intent(in), optional :: le_grid_config_in    
    
    if (present(le_grid_config_in)) le_grids_config = le_grid_config_in

    call le_grids_config%init(name = 'le_grids_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    bouncefuzz = le_grids_config%bouncefuzz
    genquad = le_grids_config%genquad
    negrid = le_grids_config%negrid
    nesub = le_grids_config%nesub
    nesuper = le_grids_config%nesuper
    new_trap_int = le_grids_config%new_trap_int
    npassing = le_grids_config%npassing
    ngauss = le_grids_config%ngauss
    nmax = le_grids_config%nmax
    nterp = le_grids_config%nterp
    radau_gauss_grid = le_grids_config%radau_gauss_grid
    test = le_grids_config%test
    trapped_particles = le_grids_config%trapped_particles
    vcut = le_grids_config%vcut
    wgt_fac = le_grids_config%wgt_fac

    exist = le_grids_config%exist
    
! user can choose not to set negrid (preferred for old scheme)
    if (negrid == -10) then
       negrid = nesub + nesuper

! If user chose negrid, then set nesuper and nesub accordingly
    else 
       nesuper = min(negrid/10+1, 4)
       nesub = negrid - nesuper
    endif

    ! Set ng2 from ngauss if not specified in the input file
    if (npassing < 0) then
       ng2 = 2 * ngauss
    else
       ng2 = npassing
    end if

  end subroutine read_parameters

  !> FIXME : Add documentation  
  subroutine init_integrations
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: init_dist_fn_layouts, pe_layout
    use mp, only: nproc
    implicit none
    character (1) :: char

    call init_dist_fn_layouts (naky, ntheta0, nlambda, negrid, nspec)

    if (.not. intinit) then
       intinit = .true.
       call pe_layout (char)
       if (char == 'x') then
          accel_x = mod(ntheta0*naky*nspec, nproc) == 0
          accel_v = .false.
       end if
       if (char == 'v') then
          accel_x = .false.
          accel_v = mod(negrid*nlambda*nspec, nproc) == 0
       end if
    end if
          

  end subroutine init_integrations

  !> Calculates energy grid and (untrapped) pitch angle weights for
  !> integrals which drop a point from the energy/lambda
  !> dimensions. These are only used for estimating the velocity space
  !> error in [[dist_fn::get_verr]] as a part of the adaptive
  !> collisionality algorithm.
  subroutine init_weights
    use file_utils, only: open_output_file, close_output_file
    use egrid, only: zeroes
    use constants, only: pi => dpi
    use species, only: nspec, f0_values

    implicit none

    real, dimension (:), allocatable :: modzeroes, werrtmp  ! (negrid-2)
    real, dimension (:), allocatable :: lmodzeroes, wlerrtmp ! (ng2-1)
    integer :: ipt, ndiv, divmax, is
    logical :: eflag = .false.

    if (init_weights_init) return
    init_weights_init = .true.


    allocate(lmodzeroes(ng2-1), wlerrtmp(ng2-1))
    allocate(wlerr(ng2,ng2))

    wlerr = 0.0; lmodzeroes = 0.0; wlerrtmp = 0.0

    wdim = nesub
    allocate(modzeroes(nesub-1), werrtmp(nesub-1))
    allocate(werr(negrid-1,nesub,nspec))

    werr = 0.0 ; modzeroes = 0.0 ; werrtmp = 0.0
    
    ! loop to obtain weights for energy grid points.  negrid-1 sets
    ! of weights are needed because we want to compute integrals
    ! for negrid-1 sets of energy points (corresponding to negrid-1
    ! points that we can choose to drop from the guassian quadrature)
    
    do ipt=1,nesub
       do is = 1,nspec
       
         ! drops the point corresponding to ipt from the energy grid
         
         if (ipt /= 1) modzeroes(:ipt-1) = zeroes(:ipt-1,is)
         if (ipt /= nesub) modzeroes(ipt:nesub-1) = zeroes(ipt+1:nesub,is)
         
         ! get weights for energy grid points
         
         call get_weights (nmax,0.0,vcut,modzeroes,werrtmp,ndiv,divmax,eflag)
         
         ! a zero is left in the position corresponding to the dropped point
         
         if (ipt /= 1) werr(:ipt-1,ipt,is) = werrtmp(:ipt-1)
         if (ipt /= nesub) werr(ipt+1:nesub,ipt,is) = werrtmp(ipt:nesub-1)
         werr(nesub+1:,ipt,is) = w(nesub+1:negrid-1,is)
         
         ! absorbing volume element into weights
         werr(:nesub,ipt,is) = werr(:nesub,ipt,is)*energy(:nesub,is)*2.0*pi*f0_values(:nesub,is)
       end do
       
    end do

    ! same thing done here for lamdba as was
    ! done earlier for energy space

    do ipt=1,ng2

       if (ipt /= 1) lmodzeroes(:ipt-1) = xx(:ipt-1)
       if (ipt /= ng2) lmodzeroes(ipt:ng2-1) = xx(ipt+1:)

       call get_weights (nmax,1.0,0.0,lmodzeroes,wlerrtmp,ndiv,divmax,eflag)

       if (ipt /= 1) wlerr(:ipt-1,ipt) = wlerrtmp(:ipt-1)
       if (ipt /= ng2) wlerr(ipt+1:,ipt) = wlerrtmp(ipt:)

!       do il = 1, ng2
          ! TEMP FOR TESTING -- MAB
!          write (*,*) 'wlerr', ipt, il, xx(il), wlerr(il,ipt), ndiv, divmax
!       end do

    end do

    deallocate(modzeroes,werrtmp,lmodzeroes,wlerrtmp)
    eflag = .false.

  end subroutine init_weights

  !> The get_weights subroutine determines how to divide up the integral into 
  !! subintervals and how many grid points should be in each subinterval
  subroutine get_weights (maxpts_in, llim, ulim, nodes, wgts, ndiv, divmax, err_flag)

    implicit none

    integer, intent (in) :: maxpts_in
    real, intent (in) :: llim, ulim
    real, dimension (:), intent (in) :: nodes
    real, dimension (:), intent (out) :: wgts
    logical, intent (out) :: err_flag
    integer, intent (out) :: ndiv, divmax

    integer :: npts, rmndr, basepts, divrmndr, base_idx, idiv, epts, maxpts
    integer, dimension (:), allocatable :: divpts

    real :: wgt_max

! npts is the number of grid points in the integration interval
    npts = size(nodes)

    wgts = 0.0; epts = npts; basepts = nmax; divrmndr = 0; ndiv = 1; divmax = npts

! maxpts is the max number of pts in an integration subinterval
    maxpts = min(maxpts_in,npts)

    do

!       wgt_max = wgt_fac/maxpts
       wgt_max = abs(ulim-llim)*wgt_fac/npts

! only need to subdivide integration interval if maxpts < npts
       if (maxpts .ge. npts) then
          call get_intrvl_weights (llim, ulim, nodes, wgts)
       else
          rmndr = mod(npts-maxpts,maxpts-1)
          
! if rmndr is 0, then each subinterval contains maxpts pts
          if (rmndr == 0) then
! ndiv is the number of subintervals
             ndiv = (npts-maxpts)/(maxpts-1) + 1
             allocate (divpts(ndiv))
! divpts is an array containing the # of pts for each subinterval
             divpts = maxpts
          else
             ndiv = (npts-maxpts)/(maxpts-1) + 2
             allocate (divpts(ndiv))
! epts is the effective number of pts after taking into account double
! counting of some grid points (those that are boundaries of subintervals
! are used twice)
             epts = npts + ndiv - 1
             basepts = epts/ndiv
             divrmndr = mod(epts,ndiv)
             
! determines if all intervals have same # of pts
             if (divrmndr == 0) then
                divpts = basepts
             else
                divpts(:divrmndr) = basepts + 1
                divpts(divrmndr+1:) = basepts
             end if
          end if
          
          base_idx = 0
          
! loop calls subroutine to get weights for each subinterval
          do idiv=1,ndiv
             if (idiv == 1) then
                call get_intrvl_weights (llim, nodes(base_idx+divpts(idiv)), &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             else if (idiv == ndiv) then
                call get_intrvl_weights (nodes(base_idx+1), ulim, &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             else
                call get_intrvl_weights (nodes(base_idx+1), nodes(base_idx+divpts(idiv)), &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             end if
             base_idx = base_idx + divpts(idiv) - 1
          end do
          
          divmax = maxval(divpts)

          deallocate (divpts)
       end if

       ! check to make sure the weights do not get too large
       ! @Warning: This looks a bit unusual -- we're taking the magnitude of
       ! the maximum value, whilst usually we would find the maximum value of
       ! the magnitude.
       if (abs(maxval(wgts)) .gt. wgt_max) then
          if (maxpts .lt. 3) then
             err_flag = .true.
             exit
          end if
          maxpts = divmax - 1
       else
          exit
       end if

       wgts = 0.0; epts = npts; divrmndr = 0; basepts = nmax
    end do

  end subroutine get_weights

  !> Used by [[get_weights]] to find the Lagrange quadrature weights for the given grid.
  !>
  !> Suppose we wish to calculate Int(f[x] dx) over the range [llim, ulim].
  !> We could approximate f[x] by the Lagrange interpolating polynomial -->
  !> f[x] ~ Sum_i{ f[x_i] l_i(x)} where. the Lagrange polynomial is given by
  !> l_i(x) = Pi_j{ (x - x_i)/(x_i - x_j) }. Hence our integral can be written as
  !> I = Int(f[x] dx) ~ Sum_i{ f[x_i] Int(Pi_j{ (x - x_i)/(x_i - x_j) })} so
  !> I ~ Sum_i{f[x_i] c_i} with c_i = Int(Pi_j{ (x - x_i)/(x_i - x_j) }). Here
  !> we are calculating the `c_i` values which are passed out as [[wgts]] corresponding
  !> to the set of positions {x_i} set by [[nodes]]. To evaluate the integral
  !> of the Lagrange polynomial coefficients we adopt Gauss-Legendre quadrature,
  !> which allows highly accurate evaluation of polynomials of degree 2*N-1.
  !> Given the need to find M weights we find we only need Gauss-Legendre of order
  !> N = (M + 1)/2. Here we actually use degree (M/2) + 1.
  !>
  !> One limitation of Lagrange polynomials is their lack of numerical stability
  !> at high order (i.e. when M is large, typically > 20 is quoted as a limit).
  !> The routine [[get_weights]] allows the integration domain to be sub-divided
  !> in order to limit the Lagrange order and reduce the variation in weights.
  !> This is partially controlled by the input [[nmax]].
  subroutine get_intrvl_weights (llim, ulim, nodes, wgts)
    use gauss_quad, only: get_legendre_grids_from_cheb
    
    implicit none
    
    ! llim (ulim) is lower (upper) limit of integration
    real, intent (in) :: llim, ulim
    real, dimension (:), intent (in) :: nodes
    real, dimension (:), intent (in out) :: wgts
    
    ! stuff needed to do guassian quadrature 
    real, dimension (:), allocatable :: gnodes, gwgts, omprod
    integer :: ix, iw

    allocate (gnodes(size(nodes)/2+1), gwgts(size(wgts)/2+1), omprod(size(nodes)/2+1))
    
    call get_legendre_grids_from_cheb (llim, ulim, gnodes, gwgts)

    do iw=1,size(wgts)
       omprod = 1.0
       
       do ix=1,size(nodes)
          if (ix /= iw) omprod = omprod*(gnodes - nodes(ix))/(nodes(iw) - nodes(ix))
       end do
       
       do ix=1,size(gwgts)
          wgts(iw) = wgts(iw) + omprod(ix)*gwgts(ix)
       end do
    end do

    deallocate (gnodes, gwgts, omprod)
       
  end subroutine get_intrvl_weights

  !> FIXME : Add documentation  
  subroutine integrate_species_original (g, weights, total)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use kt_grids, only: kwork_filter
    use mp, only: sum_allreduce
    use species, only: tracer_species, spec

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo

    !Ensure array is zero to begin
    total = 0.

    !Performed integral (weighted sum) over local velocity space and species
    if(any(kwork_filter))then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          if(kwork_filter(it,ik)) cycle
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total(:, it, ik) = total(:, it, ik) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total(:, it, ik) = total(:, it, ik) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    endif
    !Reduce sum across all procs to make integral over all velocity space and species
    call sum_allreduce (total) 
  end subroutine integrate_species_original


  !> Integrate species on xy subcommunicator - NO GATHER
  subroutine integrate_species_sub (g, weights, total)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, intspec_sub
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce_sub, sum_allreduce
    use kt_grids, only: kwork_filter
    use species, only: spec, tracer_species
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    complex, dimension(:,:,:),allocatable :: total_small
    integer :: is, il, ie, ik, it, iglo

    !Allocate array and ensure is zero
    if(intspec_sub)then
       !       total(:,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max)=0.
       allocate(total_small(-ntgrid:ntgrid,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max))
    else
       !total=0.
       allocate(total_small(-ntgrid:ntgrid,g_lo%ntheta0,g_lo%naky))       
    endif
    total_small=0.

    !Performed integral (weighted sum) over local velocity space and species
    if(any(kwork_filter))then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          if(kwork_filter(it,ik)) cycle
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          
          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total_small(:, it, ik) = total_small(:, it, ik) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          
          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total_small(:, it, ik) = total_small(:, it, ik) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    endif

    !Reduce sum across all procs in sub communicator to make integral over all velocity space and species
    if(intspec_sub)then
       call sum_allreduce_sub(total_small,g_lo%xyblock_comm)
    else
       call sum_allreduce(total_small)
    endif

    !Copy data into output array
    !Note: When not using sub-comms this is an added expense which will mean
    !this routine is more expensive than original version just using total.
    !In practice we should have two integrate_species routines, one for sub-comms
    !and one for world-comms.
    if(intspec_sub)then
       total(:,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max)=total_small
    else
       total=total_small
    endif

    !Deallocate
    deallocate(total_small)

  end subroutine integrate_species_sub

  !>Integrate species using gf_lo data format
  subroutine integrate_species_gf (g, weights, total)
    use species, only : nspec
    use theta_grid, only: ntgrid
    use gs2_layouts, only: gf_lo, g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce_sub, sum_allreduce
    use kt_grids, only: kwork_filter
    use redistribute, only: gather
    use species, only: spec, tracer_species
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:, :, :, :, :, :), allocatable :: gf
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:ntgrid,gf_lo%ntheta0,gf_lo%naky), intent (out) :: total
    integer :: is, il, ie, igf, it, ik

    total =0.
    allocate(gf(-ntgrid:ntgrid,2,nspec,negrid,nlambda,gf_lo%llim_proc:gf_lo%ulim_alloc))
    call gather(g2gf, g, gf, ntgrid)

    !Performed integral (weighted sum) over local velocity space and species
    if(any(kwork_filter)) then
       do igf = gf_lo%llim_proc,gf_lo%ulim_proc
          it = it_idx(gf_lo,igf)
          ik = ik_idx(gf_lo,igf)
          if(kwork_filter(it,ik)) cycle
          do il = 1,gf_lo%nlambda
             do ie = 1,gf_lo%negrid
                do is = 1,gf_lo%nspec
                   if (spec(is)%type .NE. tracer_species) then
                      total(:,it,ik) = total(:,it,ik) + weights(is)*w(ie,is)*wl(:,il)*(gf(:,1,is,ie,il,igf)+gf(:,2,is,ie,il,igf))
                   end if
                end do
             end do
          end do
       end do
    else
       do igf = gf_lo%llim_proc,gf_lo%ulim_proc
          it = it_idx(gf_lo,igf)
          ik = ik_idx(gf_lo,igf)
          do il = 1,gf_lo%nlambda
             do ie = 1,gf_lo%negrid
                do is = 1,gf_lo%nspec
                   if (spec(is)%type .NE. tracer_species) then
                      total(:,it,ik) = total(:,it,ik) + weights(is)*w(ie,is)*wl(:,il)*(gf(:,1,is,ie,il,igf)+gf(:,2,is,ie,il,igf))
                   end if
                end do
             end do
          end do
       end do
    end if

    deallocate(gf)
  end subroutine integrate_species_gf

  !> Integrate species using gf_lo data format assuming that gf has already been gathered prior to being 
  !! passed to this routine.
  !! Currently this routine isn't being used as the functionality has been directly added to 
  !! getan_nogath in dist_fn.
  !!
  !! AJ
  subroutine integrate_species_gf_nogather (gf, weights, total)
    use species, only : nspec
    use theta_grid, only: ntgrid
    use gs2_layouts, only: gf_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce_sub, sum_allreduce
    use kt_grids, only: kwork_filter
    use species, only: spec, tracer_species
    implicit none

    complex, dimension(-ntgrid:ntgrid,2,nspec,negrid,nlambda,gf_lo%llim_proc:gf_lo%ulim_alloc), intent(in) :: gf
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:ntgrid,gf_lo%ntheta0,gf_lo%naky), intent (out) :: total
    integer :: is, il, ie, igf, it, ik

    total =0.

    !Performed integral (weighted sum) over local velocity space and species
    if(any(kwork_filter)) then
       do igf = gf_lo%llim_proc,gf_lo%ulim_proc
          it = it_idx(gf_lo,igf)
          ik = ik_idx(gf_lo,igf)
          total(:,it,ik) = 0.
          if(kwork_filter(it,ik)) cycle
          do il = 1,gf_lo%nlambda
             do ie = 1,gf_lo%negrid
                do is = 1,gf_lo%nspec
                   if (spec(is)%type .NE. tracer_species) then
                      total(:,it,ik) = total(:,it,ik) + weights(is)*w(ie,is)*wl(:,il)*(gf(:,1,is,ie,il,igf)+gf(:,2,is,ie,il,igf))
                   end if
                end do
             end do
          end do
       end do
    else
       do igf = gf_lo%llim_proc,gf_lo%ulim_proc
          it = it_idx(gf_lo,igf)
          ik = ik_idx(gf_lo,igf)
          total(:,it,ik) = 0.
          do il = 1,gf_lo%nlambda
             do ie = 1,gf_lo%negrid
                do is = 1,gf_lo%nspec
                   total(:,it,ik) = total(:,it,ik) + weights(is)*w(ie,is)*wl(:,il)*(gf(:,1,is,ie,il,igf)+gf(:,2,is,ie,il,igf))
                end do
             end do
          end do
       end do
    end if
  end subroutine integrate_species_gf_nogather
  
  !> Integrate_species on subcommunicator with gather
  !! Falls back to original method if not using xyblock sub comm
  subroutine integrate_species_master (g, weights, total, nogath, gf_lo)
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use gs2_layouts, only: g_lo, intspec_sub
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce_sub, allgatherv
    use mp, only: nproc_comm,rank_comm,mp_abort
    use species, only: spec, tracer_species
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (0:,:,:), intent (out) :: total
    logical, intent(in), optional :: nogath
    logical, intent(in), optional :: gf_lo
    complex, dimension (:), allocatable :: total_flat
    complex, dimension (:,:,:), allocatable :: total_transp
    integer :: nl,nr, ik, it, iglo, ip, ie,is,il, ig


    if(present(gf_lo)) then
       if(gf_lo) then
          call integrate_species_gf(g,weights,total)
          return
       end if
    end if
    !If not using sub-communicators then just use original method
    !Note that if x and y are entirely local then we force intspec_sub=.false.
    if(.not.intspec_sub) then
       call integrate_species_original(g,weights,total)
       return
    endif

    !If we don't want to gather then use integrate_species_sub
    if(present(nogath))then
       if(nogath)then
          call integrate_species_sub(g,weights,total)
          return
       endif
    endif

    !->First intialise gather vars
    !Note: We only do this on the first call !!May be better to move this to some init routine?
    if(.not.allocated(recvcnts_intspec)) then
       !Get subcomm size
       call nproc_comm(g_lo%lesblock_comm,sz_intspec)

       !Get local rank
       call rank_comm(g_lo%lesblock_comm,local_rank_intspec)

       !Create displacement and receive count arrays
       allocate(recvcnts_intspec(sz_intspec),displs_intspec(sz_intspec))

       do ip=0,sz_intspec-1
          displs_intspec(ip+1)=MIN(g_lo%les_kxky_range(1,ip)*(2*ntgrid+1),ntheta0*naky*(2*ntgrid+1)-1)
          recvcnts_intspec(ip+1)=MAX((g_lo%les_kxky_range(2,ip)-g_lo%les_kxky_range(1,ip)+1)*(2*ntgrid+1),0)
       enddo
    endif

    !Allocate array and ensure is zero
    allocate(total_flat(g_lo%les_kxky_range(1,local_rank_intspec)*&
         (2*ntgrid+1):(1+g_lo%les_kxky_range(2,local_rank_intspec))*(2*ntgrid+1)))
    total_flat=0.

    !Performed integral (weighted sum) over local velocity space and species
    if(g_lo%x_before_y) then

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          
          !Calculate extent
          nl=(2*ntgrid+1)*(it-1+ntheta0*(ik-1))
          nr=nl+(2*ntgrid)
          
          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total_flat(nl:nr) = total_flat(nl:nr) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          !Convert from iglo to the separate indices
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          
          !Calculate extent
          nl=(2*ntgrid+1)*(ik-1+naky*(it-1))
          nr=nl+(2*ntgrid)

          !Sum up weighted g
          if (spec(is)%type .NE. tracer_species) then
             total_flat(nl:nr) = total_flat(nl:nr) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
          end if
       end do
    endif

    !Reduce sum across all procs in sub communicator to make integral over all velocity space and species
    call sum_allreduce_sub(total_flat,g_lo%xyblock_comm)

    !Now gather missing xy data from other procs (only talk to procs
    !with the same piece of les)

    if(g_lo%x_before_y)then
       call allgatherv(total_flat,recvcnts_intspec(local_rank_intspec+1),total,recvcnts_intspec,displs_intspec,g_lo%lesblock_comm)
    else
       allocate(total_transp(0:2*ntgrid,naky,ntheta0))
       call allgatherv(total_flat,recvcnts_intspec(local_rank_intspec+1),total_transp,recvcnts_intspec,displs_intspec,g_lo%lesblock_comm)
       do ig=0,2*ntgrid
          total(ig,:,:)=transpose(total_transp(ig,:,:))
       enddo
       !This is pretty bad for memory access so can do this all at once
       !using reshape with a specified order :
       !total=RESHAPE(total_transp,(/2*ntgrid+1,ntheta0,naky/),ORDER=(/1,3,2/))
       !BUT timings in a simple test code suggest loop+transpose can be faster.
       !In the case where ntgrid is large and ntheta0 is small reshape can win
       !whilst in the case where ntgrid is small and ntheta0 is large transpose wins.
       !When both are large reshape seems to win.
       !In conclusion it's not clear which method is better but if we assume we care most
       !about nonlinear simulations then small ntgrid with large ntheta0 is most likely so
       !pick transpose method
       deallocate(total_transp)
    endif

    deallocate(total_flat)
  end subroutine integrate_species_master

  !> Integrate over velocity space whilst in e_lo_LAYOUT. 
  subroutine integrate_species_e (lo, g, weights, total)
    use layouts_type, only: e_layout_type
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ig_idx, il_idx
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, sum_allreduce
    use species, only: spec, tracer_species
    implicit none

    type (e_layout_type), intent (in) :: lo
    complex, dimension (:,lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, ig, ielo

    total = 0.0
    do ielo = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ielo)
       ik = ik_idx(lo,ielo)
       it = it_idx(lo,ielo)
       is = is_idx(lo,ielo)
!CMR: in e_lo il is lambda index, and isign is separate.
       il = il_idx(lo,ielo)
       !Perform local sum
       if (spec(is)%type .NE. tracer_species) then
          do ie=1, negrid
              total(ig, it, ik) = total(ig, it, ik) + &
               weights(is)*w(ie,is)*wl(ig,il)*g(ie,ielo)
          end do
       end if
    end do
!Sum over all procs to make integral over all velocity space and species
    call sum_allreduce (total)
  end subroutine integrate_species_e

  !> Integrate over velocity space whilst in lz_lo_LAYOUT. 
  subroutine integrate_species_lz (lo, g, weights, total)
    use layouts_type, only: lz_layout_type
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ig_idx, ie_idx
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, sum_allreduce
    use species, only: spec, tracer_species
    implicit none

    type (lz_layout_type), intent (in) :: lo
    complex, dimension (:,lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer :: is, ixi, il, ie, ik, it, ig, ilzlo

    total = 0.0
    do ilzlo = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ilzlo)
       ik = ik_idx(lo,ilzlo)
       it = it_idx(lo,ilzlo)
       ie = ie_idx(lo,ilzlo)
       is = is_idx(lo,ilzlo)
       !Perform local sum
!CMR: in lz_lo il muxt be obtained from local index ixi, which includes isign.
!CMR, 2/10/2013:
!   nxi+1 limit on do loop below is CRUCIAL, as its stores phase space point
!   corresponding to g_lo (il=nlambda, isign=2).
!   This MUST contribute to the v-space integral, but is NOT
!   needed in collision operator as EQUIVALENT to g_lo(il=nlambda, isign=2).
!   (In collisions at ig=0, both of these points are EXACTLY equivalent, xi=0.)
!MRH actually nxi+1 is needed in the collision operator for consistency 16/08/2018
       if (spec(is)%type .NE. tracer_species) then
          do ixi=1, nxi+1
             il = ixi_to_il(ig,ixi)
             total(ig, it, ik) = total(ig, it, ik) + &
               weights(is)*w(ie,is)*wl(ig,il)*g(il,ilzlo)
          end do
       end if
    end do
!Sum over all procs to make integral over all velocity space and species
    call sum_allreduce (total)
  end subroutine integrate_species_lz


  !> Perform an integral over velocity space whilst in the LE_LAYOUT in 
  !! which we have ensured that all of velocity space is local. As such
  !! we don't need any calls to MPI reduction routines. Note that this means
  !! the processors for different distributed spatial points (x,y) don't know
  !! the results at other points.
  subroutine integrate_species_le (lo, g, weights, total)
    use layouts_type, only: le_layout_type
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use theta_grid, only: ntgrid
    use mp, only: sum_allreduce
    use species, only: tracer_species, spec

    implicit none

    type (le_layout_type), intent (in) :: lo
    complex, dimension (:,:,lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer :: ixi, ie, il, ile, ig, it, ik, is

    total = 0.0
    do ile = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ile)
       it = it_idx (lo,ile)
       ik = ik_idx (lo,ile)
       is = is_idx (lo,ile)

       if (spec(is)%type .NE. tracer_species) then
          do ie=1, negrid
!CMR, 2/10/2013:
!   nxi+1 limit on do loop below is CRUCIAL, as its stores phase space point
!   corresponding to g_lo (il=nlambda, isign=2).
!   This MUST contribute to the v-space integral, but is NOT
!   needed in collision operator as EQUIVALENT to g_lo(il=nlambda, isign=2).
!   (In collisions at ig=0, both of these points are EXACTLY equivalent, xi=0.)
!MRH actually nxi+1 is needed in the collision operator for consistency 16/08/2018   
             do ixi=1, nxi+1
                il = ixi_to_il(ig,ixi)
                total(ig,it,ik) = total(ig,it,ik) + weights(is)*w(ie,is) * wl(ig,il) * g(ixi,ie,ile)
             end do
          end do
       end if
    end do
!Sum over all procs to make integral over all velocity space and species
    call sum_allreduce (total)

  end subroutine integrate_species_le

  !> FIXME : Add documentation
  !>
  !> @note This routine depends on the values stored in [[xx]] as input
  !> to the legendre_polynomials method. However, with Radau grids the
  !> values in [[xx]] are not Legendre roots so this routine is probably
  !> incorrect when [[le_grids_knobs::radau_gauss_grid]] is `.true.`.
  subroutine legendre_transform (g, tote, totl, tott)  
    use egrid, only: zeroes
    use mp, only: nproc
    use theta_grid, only: ntgrid, bmag, bmax
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (0:,-ntgrid:,:,:,:), intent (out) :: tote, totl
    complex, dimension (0:,-ntgrid:,:,:,:), intent (out), optional :: tott

    complex :: totfac
    real :: ulim
    integer :: is, il, ie, ik, it, iglo, ig, im, ntrap
    integer, save :: lpesize

    real, dimension (:), allocatable :: nodes
    real, dimension (:,:), allocatable :: lpltmp, lpttmp
!    real, dimension (:,:), allocatable, save :: lpe, lpl
!    real, dimension (:,:,:,:), allocatable, save :: lpt

    if (.not. allocated(lpl)) then
       allocate(lpltmp(ng2,0:ng2-1))
       allocate(lpl(nlambda,0:ng2-1))

       lpesize = nesub
       allocate(lpe(negrid,0:lpesize-1,nspec)) ; lpe = 0.0
       
       ! get value of first nesub legendre polynomials
       ! at each of the grid points on (0,vcut)
       do is = 1,nspec
          call legendre_polynomials (0.0,vcut,zeroes(:lpesize,is),lpe(:lpesize,:,is))
       end do
       ! TEMP FOR TESTING -- MAB
       !          lpe = 2.*lpe/vcut
       
       ! get value of first ng2 legendre polynomials
       ! at each of the grid points on (0,1)
       call legendre_polynomials (0.0,1.0,xx,lpltmp)

       lpl = 0.0
       lpl(1:ng2,:) = lpltmp

       if (present(tott)) then
          allocate (lpt(nlambda,0:2*(nlambda-ng2-1),-ntgrid:ntgrid,2))
          lpt = 0.0
          do ig = -ntgrid, ntgrid
             ntrap = 1
             if (jend(ig) > ng2+1) then
                ntrap = jend(ig)-ng2
                allocate (nodes(2*ntrap-1))
                allocate (lpttmp(2*ntrap-1,0:2*(ntrap-1)))
                do il = 1, ntrap
                   nodes(il) = -sqrt(max(0.0,1.0-al(ng2+il)*bmag(ig)))
                end do
                nodes(ntrap+1:) = -nodes(ntrap-1:1:-1)
!Can we remove this?
! TEMP FOR TESTING -- MAB
!                nodes = nodes + sqrt(1.0-bmag(ig)/bmax)
!                ulim = 2.*sqrt(1.0-bmag(ig)/bmax)
                ulim = sqrt(1.0-bmag(ig)/bmax)
                call legendre_polynomials (-ulim,ulim,nodes,lpttmp)
                lpt(ng2+1:jend(ig),0:2*(ntrap-1),ig,2) = lpttmp(1:ntrap,:)
                lpt(ng2+1:jend(ig)-1,0:2*(ntrap-1),ig,1) = lpttmp(2*ntrap-1:ntrap+1:-1,:)
!Can we remove this?
!                lpt(ng2+1:jend(ig),0:2*(ntrap-1),ig,1) = lpttmp(2*ntrap-1:ntrap+1:-1,:)
!                do ie = 0, 2*(ntrap-1)
!                   do il = 1, 2*ntrap-1
!                      if (proc0) write (*,*) 'lptrap', ig, ntrap, ulim, ie, il, nodes(il), lpttmp(il,ie)
!                   end do
!                end do
!                do ie = 0, 2*(ntrap-1)
!                   do il = ng2+1,jend(ig)
!                      write (*,*) 'lpt', ig, ie, il, lpt(il,ie,ig,1), lpt(il,ie,ig,2)
!                   end do
!                end do
                deallocate (nodes, lpttmp)
             end if
          end do
       end if

       deallocate (lpltmp)
    end if

    ! carry out legendre transform to get coefficients of
    ! legendre polynomial expansion of g
    totfac = 0. ; tote = 0. ; totl = 0.
    if (present(tott)) tott = 0.

    !Loop over all indices, note this loop is optimal only for layout 'xyles' (at least in terms of
    !g memory access)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky       !Swapped ik and it loop order.
                do it = 1, ntheta0
                   iglo = idx (g_lo, ik, it, il, ie, is)
                   if (idx_local (g_lo, iglo)) then
                      do ig=-ntgrid,ntgrid
                         totfac = w(ie,is)*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
                         do im=0,lpesize-1
                            tote(im, ig, it, ik, is) = tote(im, ig, it, ik, is) + totfac*lpe(ie,im,is)*(2*im+1)
                         end do
                         do im=0,ng2-1
                            totl(im, ig, it, ik, is) = totl(im, ig, it, ik, is) + totfac*lpl(il,im)*(2*im+1)
                         end do
                         if (present(tott)) then
                            do im=0,2*(jend(ig)-ng2-1)
                               tott(im, ig, it, ik, is) = tott(im, ig, it, ik, is) + &
                                    w(ie,is)*wl(ig,il)*(lpt(il,im,ig,1)*g(ig,1,iglo)+lpt(il,im,ig,2)*g(ig,2,iglo))*(2*im+1)
                            end do
                         end if
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Now complete velocity integral, bringing back results to proc0
       call sum_reduce (tote, 0)
       call sum_reduce (totl, 0)
       if (present(tott)) call sum_reduce (tott, 0)
    end if

  end subroutine legendre_transform

  !> FIXME : Add documentation  
  subroutine legendre_polynomials (llim, ulim, xptsdum, lpdum)
    implicit none
    double precision, dimension (:), allocatable :: lp1, lp2, lp3, zshift

    real, intent (in) :: ulim, llim
    real, dimension (:), intent (in)   :: xptsdum
    real, dimension (:,0:), intent(out) :: lpdum

    integer :: j, mmax

    lpdum = 0.0

!    nmax = size(lpdum(1,:))
    mmax = size(xptsdum)

    allocate(lp1(mmax),lp2(mmax),lp3(mmax),zshift(mmax))

    lp1 = real(1.0,kind(lp1(1)))
    lp2 = real(0.0,kind(lp2(1)))

    lpdum(:,0) = real(1.0,kind(lpdum))

! TEMP FOR TESTING -- MAB
!    zshift = real(2.0,kind(zshift))*xptsdum/ulim - real(1.0,kind(zshift))
    zshift = real(2.0,kind(zshift))*(xptsdum-llim)/(ulim-llim) - real(1.0,kind(zshift))

    do j=1, size(lpdum(1,:))-1
       lp3 = lp2
       lp2 = lp1
       lp1 = ((2*j-1) * zshift * lp2 - (j-1) * lp3) / j
       lpdum(:,j) = lp1
    end do

    deallocate(lp1,lp2,lp3,zshift)

  end subroutine legendre_polynomials

  !> FIXME : Add documentation
  subroutine lagrange_interp (g, poly, istep, all)
    !<DD>WARNING: THIS ROUTINE DOESN'T USE g0 FOR ANYTHING IS THIS CORRECT?
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, sum_allreduce, nproc

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: poly
    integer, optional, intent(in) :: all
    integer, intent(in) :: istep
    real :: ypt
    integer :: il, ie, ik, it, iglo, ig

    !Initialise to zero
    poly = 0.
    do ie = 1, negrid
       do il = ng2+1, nlambda
          do it = 1, ntheta0
             do ik = 1, naky
                iglo = idx (g_lo, ik, it, il, ie, 1)
                if (idx_local (g_lo, iglo)) then
                   do ig = -ntgrid, ntgrid
                      !Do we need this if?
                      if (.not. forbid(ig,il)) then
                         ypt=sqrt(max(1.0 - bmag(ig)*al(il),0.0))                      
                      else
                         ypt=0.0
                      end if
                      poly(ig, it, ik, ie, :) = poly(ig, it, ik, ie, :) + &
                           lgrnge(ig,il,:,1)*cos(0.1*istep*ypt+1.0) + &
                           lgrnge(ig,il,:,2)*cos(-0.1*istep*ypt+1.0)
                   end do
                end if
             end do
          end do
       end do
    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Reduce integral across procs
       if (present(all)) then
          call sum_allreduce (poly)
       else
          call sum_reduce (poly, 0)
       end if
    end if

  end subroutine lagrange_interp

  !> FIXME : Add documentation
  !!
  !! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
  !! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
  !! is over all velocity space
  subroutine integrate_moment_c34 (g, total, all, full_arr)
! TT>
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx, il_idx,intmom_sub
! <TT
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, sum_allreduce_sub, nproc, sum_allreduce, sum_reduce_sub

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    complex, dimension(:,:,:,:),allocatable :: total_small
    logical, optional, intent(in) :: all
    logical, optional, intent(in) :: full_arr
    logical :: local_full_arr, local_all
    integer :: is, il, ie, ik, it, iglo

    !Do we want to know the full result?
    local_full_arr=.false.
    if(present(full_arr)) local_full_arr=full_arr

    ! Do all processors need to know the result?
    local_all = .false.
    if(present(all)) local_all = all

    !NOTE: Currently we're lazy and force the full_arr approach to reduce
    !over the whole array. Really we should still use the sub-communicator
    !approach and then gather the remaining data as we do for integrate_species

    !Allocate array and ensure is zero
    if(intmom_sub.and.local_all.and.(.not.local_full_arr))then !If we're using reduce then we don't want to make array smaller
!       total(:,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max,g_lo%is_min:g_lo%is_max)=0.
       allocate(total_small(-ntgrid:ntgrid,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max,g_lo%is_min:g_lo%is_max))
    else
!       total=0.
       allocate(total_small(-ntgrid:ntgrid,g_lo%ntheta0,g_lo%naky,g_lo%nspec))       
    endif
    total_small=0.
    !Integrate over local velocity space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Perform local sum
       total_small(:, it, ik, is) = total_small(:, it, ik, is) + &
            w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
    end do

    !Not sure that we really need to limit this to nproc>1 as if
    !we run with 1 proc MPI calls should still work ok
    if (nproc > 1) then     
       if (local_all) then 
         if (local_full_arr) then
            call sum_allreduce (total_small)
         else 
           !Complete integral over distributed velocity space and ensure all procs in sub communicator know the result
           !Note: fi intmom_sub=.false. then xysblock_comm==mp_comm  | This is why total_small must be the same size on 
           !all procs in this case.
           call sum_allreduce_sub (total_small,g_lo%xysblock_comm)
         end if
       else
         !if (local_full_arr) then
            !call sum_reduce (total_small, 0)
         !else 
           !Complete integral over distributed velocity space
           !Note: fi intmom_sub=.false. then xysblock_comm==mp_comm  | This is why total_small must be the same size on 
           !all procs in this case.
           !call sum_reduce_sub (total_small,g_lo%xysblock_comm)
         !end if

          !Complete integral over distributed velocity space but only proc0 knows the answer
          call sum_reduce (total_small, 0)
       end if
    end if

    !Copy data into output array
    !Note: When not using sub-comms this is an added expense which will mean
    !this routine is more expensive than original version just using total.
    !In practice we should have two integrate_moment_c34 routines, one for sub-comms
    !and one for world-comms.
    if(intmom_sub.and.local_all.and.(.not.local_full_arr))then
       total(:,g_lo%it_min:g_lo%it_max,g_lo%ik_min:g_lo%ik_max,g_lo%is_min:g_lo%is_max)=total_small
    else
       total=total_small
    endif

    !Deallocate
    deallocate(total_small)

  end subroutine integrate_moment_c34

  !> FIXME : Add documentation
  !!
  !! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
  !! NOTE: Takes f = f(y, z, sigma, lambda, E, species) and returns int f, where the integral
  !! is over all velocity space
  subroutine integrate_moment_r33 (g, total, all)
    use mp, only: nproc
    use gs2_layouts, only: p_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: sum_reduce, sum_allreduce
    use theta_grid, only: ntgrid

    implicit none

    real, dimension (-ntgrid:,:,p_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    integer :: is, il, ie, ik, iplo

    !Ensure zero to start
    total = 0.

    !Do local velocity space integral
    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ik = ik_idx(p_lo,iplo)
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)

       total(:, ik, is) = total(:, ik, is) + &
            w(ie,is)*wl(:,il)*(g(:,1,iplo)+g(:,2,iplo))

    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Complete distributed integral
       if (present(all)) then
          !Return result to all procs
          call sum_allreduce (total)
       else
          !Only proc0 knows the result
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_moment_r33

  !> Integrate over velocity space whilst in e_lo_LAYOUT. 
  subroutine integrate_moment_e (lo, g, total, all)
    use layouts_type, only: e_layout_type
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ig_idx, il_idx
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, sum_allreduce
    implicit none

    type (e_layout_type), intent (in) :: lo
    complex, dimension (:,lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    integer :: is, il, ie, ik, it, ig, ielo

    total = 0.0
    do ielo = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ielo)
       ik = ik_idx(lo,ielo)
       it = it_idx(lo,ielo)
       is = is_idx(lo,ielo)
!CMR: in e_lo il is lambda index, and isign is separate.
       il = il_idx(lo,ielo)
       !Perform local sum
       do ie=1, negrid
           total(ig, it, ik, is) = total(ig, it, ik, is) + &
            w(ie,is)*wl(ig,il)*g(ie,ielo)
       end do
    end do

    if (present(all)) then
       !Complete vspace integral and get result to all procs
       call sum_allreduce (total)
    else
       !Complete vspace integral, sending result to proc0
       call sum_reduce (total, 0)
    end if
  end subroutine integrate_moment_e

  !> Integrate over velocity space whilst in lz_lo_LAYOUT. 
  subroutine integrate_moment_lz (lo, g, total, all)
    use layouts_type, only: lz_layout_type
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ig_idx, ie_idx
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, sum_allreduce
    implicit none

    type (lz_layout_type), intent (in) :: lo
    complex, dimension (:,lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    integer :: is, ixi, il, ie, ik, it, ig, ilzlo

    total = 0.0
    do ilzlo = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ilzlo)
       ik = ik_idx(lo,ilzlo)
       it = it_idx(lo,ilzlo)
       ie = ie_idx(lo,ilzlo)
       is = is_idx(lo,ilzlo)
       !Perform local sum
!CMR: in lz_lo il muxt be obtained from local index ixi, which includes isign.
!CMR, 2/10/2013:
!   nxi+1 limit on do loop below is CRUCIAL, as its stores phase space point
!   corresponding to g_lo (il=nlambda, isign=2).
!   This MUST contribute to the v-space integral, but is NOT
!   needed in collision operator as EQUIVALENT to g_lo(il=nlambda, isign=2).
!   (In collisions at ig=0, both of these points are EXACTLY equivalent, xi=0.)
!MRH actually nxi+1 is needed in the collision operator for consistency 16/08/2018!  
       do ixi=1, nxi+1
          il = ixi_to_il(ig,ixi)
          total(ig, it, ik, is) = total(ig, it, ik, is) + &
            w(ie,is)*wl(ig,il)*g(il,ilzlo)
       end do
    end do

    if (present(all)) then
       !Complete vspace integral and get result to all procs
       call sum_allreduce (total)
    else
       !Complete vspace integral, sending result to proc0
       call sum_reduce (total, 0)
    end if
  end subroutine integrate_moment_lz

  !> Perform an integral over velocity space whilst in the LE_LAYOUT in 
  !! which we have ensured that all of velocity space is local. As such
  !! we don't need any calls to MPI reduction routines. Note that this means
  !! the processors for different distributed spatial points (x,y) don't know
  !! the results at other points.
  subroutine integrate_moment_lec (lo, g, total)
    use layouts_type, only: le_layout_type
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use kt_grids, only: kwork_filter
    implicit none

    type (le_layout_type), intent (in) :: lo
    complex, dimension (:,:,lo%llim_proc:), intent (in) :: g
    complex, dimension (lo%llim_proc:), intent (out) :: total
    integer :: ixi, ie, il, ile, ig, it, ik, is
    integer :: nxup
    total = cmplx(0.0,0.0)
    if (nxi .gt. 2*ng2) then
       nxup=nxi+1
    else
       nxup=nxi
    endif
    do ile = lo%llim_proc, lo%ulim_proc
       is = is_idx (lo,ile)
       it = it_idx (lo,ile)
       ik = ik_idx (lo,ile)
       if(kwork_filter(it,ik)) cycle
       ig = ig_idx (lo,ile)
       do ie=1, negrid
!CMR, 2/10/2013:
!   nxi+1 limit on do loop below is CRUCIAL, as its stores phase space point
!   corresponding to g_lo (il=nlambda, isign=2).
!   This MUST contribute to the v-space integral, but is NOT
!   needed in collision operator as EQUIVALENT to g_lo(il=nlambda, isign=2).
!   (In collisions at ig=0, both of these points are EXACTLY equivalent, xi=0.)
!MRH actually nxi+1 is needed in the collision operator for consistency 16/08/2018  
          do ixi=1, nxup
             il = ixi_to_il(ig,ixi)
             total(ile) = total(ile) + w(ie,is) * wl(ig,il) * g(ixi,ie,ile)
          end do
       end do
    end do
    ! No need for communication since all velocity grid points are together
    ! and each prcessor does not touch the unset place
    ! They actually don't need to keep all 4D array
    ! Do we stay in le_layout for total?
    ! --- ile contains necessary and sufficient information for (ig,it,ik,is)

  end subroutine integrate_moment_lec

  !> Perform an integral over velocity space whilst in the LE_LAYOUT in 
  !! which we have ensured that all of velocity space is local. As such
  !! we don't need any calls to MPI reduction routines. Note that this means
  !! the processors for different distributed spatial points (x,y) don't know
  !! the results at other points.
  subroutine integrate_moment_le (lo, g, total)
    use layouts_type, only: le_layout_type
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use theta_grid, only: ntgrid
    use kt_grids, only: kwork_filter
    implicit none

    type (le_layout_type), intent (in) :: lo
    complex, dimension (:,:,lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: ixi, ie, il, ile, ig, it, ik, is
    total = cmplx(0.0,0.0)
    do ile = lo%llim_proc, lo%ulim_proc
       it = it_idx (lo,ile)
       ik = ik_idx (lo,ile)
       if(kwork_filter(it,ik)) cycle
       ig = ig_idx (lo,ile)
       is = is_idx (lo,ile)
       do ie=1, negrid
!CMR, 2/10/2013:
!   nxi+1 limit on do loop below is CRUCIAL, as its stores phase space point
!   corresponding to g_lo (il=nlambda, isign=2).
!   This MUST contribute to the v-space integral, but is NOT
!   needed in collision operator as EQUIVALENT to g_lo(il=nlambda, isign=2).
!   (In collisions at ig=0, both of these points are EXACTLY equivalent, xi=0.)
!MRH actually nxi+1 is needed in the collision operator for consistency 16/08/2018  
          do ixi=1, nxi+1
             il = ixi_to_il(ig,ixi)
             total(ig,it,ik,is) = total(ig,it,ik,is) + w(ie,is) * wl(ig,il) * g(ixi,ie,ile)
          end do
       end do
    end do
    ! No need for communication since all velocity grid points are together
    ! and each prcessor does not touch the unset place
  end subroutine integrate_moment_le

  !> FIXME : Add documentation
  !!
  !! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
  !! NOTE: Takes f = f(y, lambda, E, species) and returns int sum_{ky} f, where the integral
  !! is over energy and lambda (not sigma)
  subroutine integrate_kysum (g, ig, total, all)
    use species, only: nspec
    use kt_grids, only: aky
    use constants, only: zi
    use gs2_layouts, only: is_idx, ik_idx, ie_idx, il_idx, p_lo
    use mp, only: sum_reduce, sum_allreduce, nproc

    implicit none

    complex, dimension (p_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: ig
    complex, dimension (:), intent (out) :: total
    integer, optional, intent(in) :: all

    complex, dimension (negrid,nlambda,nspec) :: gksum
    integer :: is, il, ie, ik, iplo

    !Initialise both arrays to zero
    total = 0. ; gksum = 0.
    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ik = ik_idx(p_lo,iplo)
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)
       gksum(ie,il,is) = gksum(ie,il,is) + real(aky(ik)*g(iplo)) + zi*aimag(aky(ik)*g(iplo))
    end do
    ! real part of gksum is | sum_{ky} ky * J0 * real[ ky*(conjg(phi+)*h- + conjg(phi-)*h+ ] |**2
    ! imag part of gksum is | sum_{ky} ky * J0 * aimag[ ky*(conjg(phi+)*h- + conjg(phi-)*h+ ] |**2
    gksum = real(gksum)**2 + zi*aimag(gksum)**2

    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)

       total(is) = total(is) + w(ie,is)*wl(ig,il)*gksum(ie,il,is)
    end do

    !Do we really need this if?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_kysum

  !> FIXME : Add documentation  
  subroutine lint_error (g, weights, total)
    use theta_grid, only: ntgrid, bmag, bmax
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt

    !If the weights array hasn't been filled in then do so now
    if (.not. allocated (wlmod)) then
       !Allocate array, don't need to initialise as below loop ensures
       !all elements are assigned a value
       allocate (wlmod(-ntgrid:ntgrid,nlambda,ng2))

       if (proc0) then
          do ipt = 1, ng2
             do il = 1, ng2
                wlmod(:,il,ipt) = wlerr(il,ipt)*2.0*sqrt((bmag(:)/bmax) &
                     *((1.0/bmax-al(il))/(1.0/bmag(:)-al(il))))
             end do
             !If we have trapped particles use the precalculated weights
             !in wlmod as above is only for passing particles
             if (nlambda > ng2) wlmod(:,ng2+1:,ipt) = wl(:,ng2+1:)
          end do
       end if

       !Now send the calculated value from proc0 to all other procs
       !We could just do the above calculations on all procs?
       call broadcast (wlmod)
    end if

    !Initialise to zero
    total = 0.

    !For each (passing) lambda point do velocity space integral
    do ipt=1,ng2
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*w(ie,is)*wlmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this outside of the ipt loop above
    call sum_allreduce (total) 

  end subroutine lint_error

  !> FIXME : Add documentation  
  subroutine trap_error (g, weights, total)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt, ntrap

    !How many trapped pitch angles are there?
    ntrap = nlambda - ng2

    !If weights not calculated yet do so now
    if (.not. allocated(wtmod)) then
       !Allocate array, don't need to initialise as below loops
       !ensure every element is assigned a value
       allocate (wtmod(-ntgrid:ntgrid,nlambda,ntrap))
          
       if (proc0) then
          do ipt=1,ntrap
             wtmod(:,:ng2,ipt) = wl(:,:ng2)
          end do

!Left below comments, but are we done testing this now?
! next line only to be used when testing!!!!
!          wtmod(:,:ng2,:) = 0.

          wtmod(:,ng2+1:,:) = wlterr(:,ng2+1:,:)
       endif

       !Send from proc0 to all others | We could just do the above calculations on all procs?
       call broadcast (wtmod)
    end if


    !Initialise to zero
    total = 0.

    !Loop over number of trapped points
    do ipt=1,ntrap
       !Do local velocity integral
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*w(ie,is)*wtmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this out of ipt loop above
    call sum_allreduce (total) 

  end subroutine trap_error

  !> FIXME : Add documentation  
  subroutine eint_error (g, weights, total)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt

    !If we don't have the weights then calculate them now
    if (.not. allocated(wmod)) then
       !Allocate array, don't initialise as we fill in all values below
       allocate (wmod(negrid,wdim,nspec))

       if (proc0) then
          wmod(:negrid-1,:,:) = werr(:,:,:)
          do is = 1,nspec
            wmod(negrid,:,is) = w(negrid,is)  
          end do
       end if

       !send from proc0 to everywhere else
       call broadcast(wmod)
    end if

    !Initialise to zero
    total=0.

    !Do velocity space integral for each ipt (for all energy grid points)
    do ipt=1,wdim
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*wmod(ie,ipt,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this out of the above loop over ipt
    call sum_allreduce (total) 
  end subroutine eint_error

  !> FIXME : Add documentation  
  subroutine set_grids
    use species, only: init_species
    use egrid, only: setvgrid, x0, x0_maxw, zeroes, zeroes_maxw
    use species, only: nspec, spec, f0_maxwellian
    use theta_grid, only: init_theta_grid, ntgrid, nbset, bset, eps_trapped
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    implicit none

    integer :: tsize,is
    logical :: has_maxwellian_species
    integer :: unit

    call init_theta_grid
    call init_species

    allocate (dele(negrid,nspec),speed(negrid,nspec),speed_maxw(negrid),w_maxw(negrid),energy_maxw(negrid))

    w_maxw = 0.0
    energy_maxw = 0.0
    speed_maxw = 0.0
    has_maxwellian_species = .false. 

    speed = sqrt(energy)

    do is = 1,nspec
      if  (spec(is)%f0type .EQ. f0_maxwellian) then
         has_maxwellian_species = .true.
         w_maxw(:) = w(:,is)
         energy_maxw(:) = energy(:,is)
         speed_maxw(:) = speed(:,is)
         zeroes_maxw(:) = zeroes(:,is)
         x0_maxw = x0(is)
         exit
      end if
    end do
    if (.not. has_maxwellian_species) write(*,*) &
        'Warning; no maxwellian species; collisions will fail'
    !call setvgrid (vcut, negrid, energy, w, nesub)
    
    tsize = 2*nterp-1

    dele(1,:) = energy(1,:)
    dele(2:,:) = energy(2:,:)-energy(:negrid-1,:)

    if (trapped_particles .and. eps_trapped > epsilon(0.0)) then
       nlambda = ng2+nbset
       lmax = nlambda-1
    else
       nlambda = ng2
       lmax = nlambda
    end if
    nxi = max(2*nlambda-1, 2*ng2)
    allocate (al(nlambda), delal(nlambda))
    allocate (wl(-ntgrid:ntgrid,nlambda))
    allocate (jend(-ntgrid:ntgrid))
    allocate (forbid(-ntgrid:ntgrid,nlambda))
    allocate (lgrnge(-ntgrid:ntgrid,nlambda,tsize,2))
    allocate (xloc(-ntgrid:ntgrid,tsize))
    if (nlambda-ng2 > 0) then
       al(ng2+1:nlambda) = 1.0/bset
    end if
    call lgridset
    delal(1) = al(1)
    delal(2:) = al(2:) - al(:nlambda-1)

    if (proc0) then
       call open_output_file(unit, '.vspace_integration_error')
       call report_velocity_integration_error_estimate(unit)
       call close_output_file(unit)
    end if

  end subroutine set_grids

  !> Calculate estimates of the velocity space integration errors
  !> and report to screen / specified unit.
  !>
  !> This error estimate provides an estimate on the lower bound
  !> for the error on subsequent calls to integrate_species/integrate_moment
  !> Alternative error estimates can be calculated following the approach
  !> adopted in [[dist_fn::get_verr]].
  subroutine report_velocity_integration_error_estimate(unit_in)
    use iso_fortran_env, only: output_unit
    use theta_grid, only: ntgrid, bmag, bmax
    use species, only: nspec
    ! Note in fortran 2008 we could just use the erf intrinsic
    use spfunc, only: erf_ext
    use constants, only: pi
    implicit none

    integer, intent(in), optional :: unit_in
    real, parameter :: expected_lambda = 2.0, expected_energy = 0.25
    real, dimension(:), allocatable :: lambda_error, energy_error, b_ratio
    real, dimension(:), allocatable :: passing_error, expected_passing
    real, dimension(:), allocatable :: trapped_error, expected_trapped
    real, dimension(:), allocatable :: energy_sub_error, energy_super_error
    real, dimension(:, :), allocatable :: mixed_error, mixed_value
    real :: expected_energy_sub, expected_energy_super
    integer :: unit, ie, ig, is

    unit = output_unit
    if (present(unit_in)) unit = unit_in

    !Get estimates of the error on the integration weights.

    ! Total pitch angle grid
    allocate(lambda_error(-ntgrid:ntgrid))
    lambda_error = 0.0
    do ig = -ntgrid, ntgrid
       lambda_error(ig) = expected_lambda - sum(wl(ig, :))
    end do
    ! Introduce factor 2 here to account for sigma doubling
    lambda_error = 2*lambda_error

    ! The trapped/passing error breakdown is currently only
    ! appropriate when not using the Radau-Gauss grids
    if (.not. radau_gauss_grid) then

       ! Passing pitch angle grid
       allocate(passing_error(-ntgrid:ntgrid))
       passing_error = 0.0
       allocate(expected_passing(-ntgrid:ntgrid))
       allocate(b_ratio(-ntgrid:ntgrid))
       b_ratio = bmax / bmag
       ! This is the analytic integral of the function we weight the Legendre
       ! weights by over the passing region, i.e.
       ! 2.0*sqrt((bmag(ig)/bmax)*((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il)))
       ! Note to integrate this correctly it is useful to transform from lambda
       ! to the Legendre variable xx, through xx = sqrt(1-lambda*bmax)
       ! The weight function then becomes:
       ! 2.0*sqrt(bmag/bmax) * sqrt([xx^2]/[(bmax/bmag)-1+xx^2])
       ! and we would integrate from xx=0 to xx=1
       ! Note this analytic result and comparison is appropriate when
       ! Radau_Gauss_Grid = F.  When this is not the case (default) then
       ! the calculation of the passing and trapped errors from this
       ! analytic result is no longer appropriate as the wfb contains
       ! contributions from both the passing and trapped regions when
       ! bmag = bmax, and as such our comparison is only valid away from
       ! where bmag=bmax.
       expected_passing = 2.0 * (1 - sqrt(b_ratio - 1)/sqrt(b_ratio))
       do ig = -ntgrid, ntgrid
          passing_error(ig) = expected_passing(ig) - sum(wl(ig, :ng2))
       end do
       ! Introduce factor 2 here to account for sigma doubling
       passing_error = 2*passing_error

       ! Trapped pitch angle grid
       allocate(trapped_error(-ntgrid:ntgrid))
       trapped_error = 0.0
       ! Only calculate this if we have any trapped pitch angles
       if (nlambda > ng2) then
          allocate(expected_trapped(-ntgrid:ntgrid))
          ! We're essentially just subtracting the passing analytic
          ! result from the total expectation, i.e. effectively
          !  expected_trapped = expected_lambda - expected_passing
          expected_trapped = 2.0 * sqrt(b_ratio - 1)/sqrt(b_ratio)
          do ig = -ntgrid, ntgrid
             trapped_error(ig) = expected_trapped(ig) - sum(wl(ig, ng2+1:))
          end do
          !Introduce factor 2 here to account for sigma doubling
          trapped_error = 2*trapped_error
       end if
    end if

    ! Total energy grid
    allocate(energy_error(nspec))
    energy_error = 0.0
    do is = 1, nspec
       energy_error(is) = expected_energy - sum(w(:,is))
    end do

    ! Energy sub grid
    allocate(energy_sub_error(nspec))
    energy_sub_error = 0.0
    ! NOTE: This assumes we are considering a species with Maxwellian
    ! background such that f0_values = exp(-energy)/pi^3/2
    ! We are analytically evaluating
    ! Integral_0^vcut{v^2 * pi * f0_values dv}
    ! Assuming Maxwellian this becomes
    ! Integral_0^vcut{v^2 exp(-v^2)/sqrt(pi) dv}
    ! Which has result:
    ! (1/4) [ erf(vcut) - 2 vcut exp(-vcut^2)/sqrt(pi)]
    expected_energy_sub = 0.25 * ( erf_ext(vcut) - 2*vcut*exp(-vcut*vcut)/sqrt(pi))
    do is = 1, nspec
       energy_sub_error(is) = expected_energy_sub - sum(w(:nesub,is))
    end do

    ! Energy super grid
    allocate(energy_super_error(nspec))
    energy_super_error = 0.0
    if (negrid > nesub) then
       ! NOTE: This assumes we are considering a species with Maxwellian
       ! background such that f0_values = exp(-energy)/pi^3/2
       ! We are analytically evaluating
       ! Integral_0^infinity{sqrt(y+vcut^2) * exp(y) * pi * f0_values * W(y) dy}/2
       ! where W(y) are the Laguerre weights = exp(-y)
       ! Assuming Maxwellian this becomes
       ! Integral_0^infinity{sqrt(y+vcut^2) * exp(y) * exp(-(y+vcut^2)) * exp(-y) dy}/2
       ! which simplifies slightly to
       ! Integral_0^infinity{sqrt(y+vcut^2) * exp(-(vcut^2)) * exp(-y) dy}/2
       ! Which has result:
       ! [sqrt(pi)*erfc(vcut)/2 + vcut*exp(-vcut^2)]/(2*sqrt(pi))
       expected_energy_super =  0.5*(sqrt(pi)*(1-erf_ext(vcut))/2 + vcut*exp(-vcut*vcut))/sqrt(pi)
       do is = 1, nspec
          energy_super_error(is) = expected_energy_super - sum(w(nesub+1:,is))
       end do
    end if

    allocate(mixed_error(-ntgrid:ntgrid, nspec))
    allocate(mixed_value(size(wl(0,:)), size(w(:,1))))
    mixed_error = 0.0
    do is = 1, nspec
       do ig = -ntgrid, ntgrid
          mixed_value = 0.0
          do ie = 1, size(w(:,1))
             mixed_value(:, ie) = wl(ig,:)*w(ie, is)
          end do
          mixed_error(ig, is) = expected_lambda*expected_energy - sum(mixed_value)
       end do
    end do

    ! Introduce factor 2 here to account for sigma doubling
    mixed_error = 2 * mixed_error

    ! Report errors to unit
    call report_error(unit, "Lambda", lambda_error)
    ! The trapped/passing error breakdown is currently only
    ! appropriate when not using the Radau-Gauss grids
    if (.not. radau_gauss_grid) then
       call report_error(unit, "Passing", passing_error)
       call report_error(unit, "Trapped", trapped_error)
    end if
    call report_error(unit, "Energy", energy_error)
    call report_error(unit, "Energy_Sub", energy_sub_error)
    call report_error(unit, "Energy_Super", energy_super_error)
    call report_error_2d(unit, "Mixed", mixed_error)

  contains
    subroutine report_error(unit, name, error)
      implicit none
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, dimension(:), intent(in) :: error
      write(unit, '(A," weights errors, max/mean/mean(|error|)  :",3(E14.6E3," "))') &
           name, maxval(abs(error)), sum(error)/size(error), sum(abs(error))/size(error)
    end subroutine report_error

    subroutine report_error_2d(unit, name, error)
      implicit none
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, dimension(:, :), intent(in) :: error
      write(unit, '(A," weights errors, max/mean/mean(|error|)  :",3(E14.6E3," "))') &
           name, maxval(abs(error)), sum(error)/size(error), sum(abs(error))/size(error)
    end subroutine report_error_2d

  end subroutine report_velocity_integration_error_estimate

  !> Determines the location lambda ([[al]]) grid points in the
  !> passing domain and associated integrations weights. Also
  !> determines the integration weights for the trapped pitch angle
  !> grid (the location of trapped [[al]] is already determined
  !> through the values of [[bset]]).
  !>
  !> Alongside determining the passing locations and the lambda
  !> integration weights we also set [[forbid]], [[jend]] and
  !> call [[xigridset]] to setup the \(\xi = \v_\|/v\) grid
  !> and associated quantities.
  subroutine lgridset
    use theta_grid, only: eps_trapped
    implicit none
    integer :: ig, il
    logical :: has_trapped_points

    ! Decide if we need to calculate the trapped lambda integration weights
    has_trapped_points = trapped_particles .and. eps_trapped > epsilon(0.0)

    ! Intialise the weights to zero
    wl = 0.0

    ! Calculate the passing lambda grid location and integration weights
    call setup_passing_lambda_grids(al, wl)

    ! Calculate the trapped lambda grid integration weights
    if (has_trapped_points) call setup_trapped_lambda_grids(al, wl)

    ! Calculate the forbid flag indicating which pitch angles are allowed at
    ! various points.
    call calculate_forbidden_region(forbid)

    ! Calculate the jend values indicating the pitch angle which bounces at
    ! each theta point.
    call calculate_jend(jend)

    ! Setup the xi (v||/v) grid used in collisions
    call xigridset

  end subroutine lgridset

  !> Determines the location of lambda ([[al]]) grid points in the
  !> passing domain and the associated integrations weights.
  subroutine setup_passing_lambda_grids(al, weights)
    use gauss_quad, only: get_legendre_grids_from_cheb, get_radau_gauss_grids
    use theta_grid, only: ntgrid, bmag, bmax
    implicit none
    !> Array which stores the lambda grid points
    real, dimension(:), intent(in out) :: al
    !> Array which stores the integration weights for the lambda grid
    real, dimension(-ntgrid:, :), intent(in out) :: weights
    real, dimension (:), allocatable :: wx, xx_radau
    integer :: ig, il

    ! Note this is a module level quantity that we keep around for possible
    ! future use. For iproc /= 0 this allocation is actually done in [[broadcast_results]]
    allocate (xx(ng2))

    if(radau_gauss_grid .and. nlambda > ng2 ) then

       ! This grid uses a fixed endpoint such that il = ng2+1 has a finite weight at bounce points
       ! which contributes to the integration over passing particles
       ! which is why xx, wx are one element longer than the standard Legendre Gauss grid
       ! This grid should only be used when nlambda > ng2, i.e. there is a wfb and trapped particles

       allocate (xx_radau(ng2 +1))
       allocate (wx(ng2 +1))
       call get_radau_gauss_grids(0.,1., xx_radau, wx,report_in=.false.)
       xx = xx_radau(:ng2)
    else
       allocate (wx(ng2))
       call get_legendre_grids_from_cheb (1., 0., xx, wx)
    end if

    ! Store the location of the passing lambda grid points
    al(:ng2) = (1.0 - xx(:ng2)**2)/bmax

    ! Transform the Legendre/Radau weights to the lambda grid
    do il = 1, ng2
       do ig = -ntgrid, ntgrid
          weights(ig,il) = wx(il)*2.0*sqrt((bmag(ig)/bmax) &
               *((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il))))
       end do
    end do

    ! Assign the weight for wfb from the radau-gauss grid
    if (radau_gauss_grid .and. nlambda > ng2 ) then
       il=ng2+1
       do ig = -ntgrid, ntgrid
          if (bmag(ig) < bmax) then ! i.e. we are not the bounce point for wfb
             ! Note this will currently always be exactly 0
             ! as al(ng2+1) == 1.0/bmax
             weights(ig,il) = wx(il)*2.0*sqrt((bmag(ig)/bmax) &
                  *((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il))))
          else ! at the bounce point
             weights(ig,il) = wx(il)*2.0
          endif
       end do
    end if

    if(allocated(wx)) deallocate(wx)
    if(allocated(xx_radau)) deallocate(xx_radau)
  end subroutine setup_passing_lambda_grids

  !> Determines the lambda grid point integration weights in the
  !> trapped domain using either new or old (default) methods,
  !> determined by [[le_grids_knobs::new_trap_int]] with
  !> new=high-order interp, old=finite difference).
  !>
  !> The new method uses Lagrange polynomial interpolation to provide
  !> an accurate approximation of the integration over the trapped
  !> region. The old method effectively uses the trapezium rule to
  !> integrate. In both cases the integration variable is the parallel
  !> velocity spanning the range [-v_||_wfb, v_||_wfb].  Both methods
  !> give zero weight to the wfb at the locations where B == Bmax.
  !> Generally the old method appears to be accurate to at most single
  !> precision, whilst the new method is more accurate, reaching the
  !> double precision round off level even with just a few trapped
  !> pitch angles. It can be noted that the the old method can reach
  !> the same level of accuracy at some theta grid points but not
  !> all. Further work may be able to optimise the grid to avoid these
  !> locations.  There are a number of reasons why we may not want to
  !> make new_trap_int default to true, including :
  !>   1. It is not compatible with radau_gauss_grid = T (which is the
  !>   default).
  !>   2. One might expect high order Lagrange interpolation to be
  !>      unstable so may anticipate that the new method will not work
  !>      well at higher ntheta. Empirically it seems this does not
  !>      become an issue here. Further, the input `nmax` can be used
  !>      to limit the maximum Lagrange order used by splitting the
  !>      domain.
  !>   3. For large ntheta the new method leads to slow
  !>      initialisation. This is due to the calculation of the error
  !>      coefficients for use with get_verr. This could perhaps be
  !>      optmised.
  !> Both points 2 and 3 could be avoided to some extent if the
  !> trapped pitch angle grid resolution is decoupled from the theta
  !> grid resolution. It may be possible to make the new approach
  !> compatible with radau_gauss_grid or we may decide
  !> radau_gauss_grid should default to false.
  !>
  !> @note Here we overwrite/replace any previously calculated weights
  !> for points considered trapped. This is usually ok as we have a clear
  !> separation between passing and trapped so the trapped weights are
  !> zero on entry to this routine. In some situations, however, we
  !> may consider some points as both passing and trapped and run the
  !> risk of discarding any existing passing weights for such points.
  !> This could potentially be the case with Radau-like grids, for example
  !> (although likely isn't an issue currently as the routine here does not
  !> calculate a weight for the wfb when bmag = bmax). A simple fix is to
  !> replace `weights = ` with effectively `weights +=` in this routine.
  subroutine setup_trapped_lambda_grids(al, weights)
    use theta_grid, only: ntgrid, bmag, bmax, eps_trapped
    use file_utils, only: open_output_file, close_output_file
    implicit none
    !> Array which stores the lambda grid points
    real, dimension(:), intent(in out) :: al
    !> Array which stores the integration weights for the lambda grid
    real, dimension(-ntgrid:, :), intent(in out) :: weights
    real, dimension (:), allocatable :: ytmp, yb, yberr, wb, wberrtmp
    real, dimension (:,:), allocatable :: wberr
    integer :: npts, ix, ntrap
    real :: wwo, llim, ulim
    logical :: eflag
    integer :: ig, il, ndiv, divmax, divmaxerr, ndiverr

    if (new_trap_int) then
       eflag = .false.

       ! max number of trapped particles (occurs at outboard midplane)
       ntrap = nlambda - ng2

       ! wlterr contains weights for less accurate trapped integrals (for error estimation)
       allocate(wlterr(-ntgrid:ntgrid,nlambda,ntrap))

       wlterr = 0.0

       ! Find the integration weights for each theta grid point
       do ig = -ntgrid, ntgrid

          ! First we count how many trapped lambda bounce points we
          ! have to consider at this point. In other words how many
          ! bounce points are not forbidden for the current theta
          ! grid point.
          ! npts is the number of al values in the trapped integral (varies with theta)
          npts = 0

          do il = ng2+1, nlambda
             if (1.0 - al(il)*bmag(ig) > -bouncefuzz) then
                npts = npts + 1
             end if
          end do

          ! If there are any valid bounce points then we need to
          ! calculate the weights. If npts == 1 then we probably also expect
          ! the integration weights to be 0 so we could probably skip all of
          ! this work for npts = 1 as well.
          if (npts > 0) then

             ! ytmp is an array containing pitch angle grid points (for vpa >= 0)
             ! These are the v||/v points for each lambda between wfb (sets maximum
             ! v||/v here) and the lambda which bounces at this point (v||/v = 0)
             ! yb is an array containing the full set of [-v||/v, v||/v] grid points.
             ! Constructed directly from ytmp by copying and flipping sign+orde
             ! wb is an array containing the integration weights calculated for the
             ! given v||/v grid.
             allocate(ytmp(npts), yb(2*npts-1), wb(2*npts-1))
             ! wberr is an array used to hold the weights used in the error estimation
             ! in [[trap_error]] called by [[get_verr]]
             allocate(wberr(2*npts-1,npts))

             ytmp = 0.0; yb = 0.0; wb = 0.0; wberr = 0.0

             ! loop computes transformed variable of integration (v||/v)
             do il = ng2+1, ng2+npts
                ytmp(il-ng2) = sqrt(max(1 - al(il)*bmag(ig), 0.0))
             end do

             ! define array (yb) with pitch-angle gridpts
             ! corresponding to both positive and negative vpa
             if (npts > 1) yb(:npts-1) = -ytmp(:npts-1)
             yb(npts:) = ytmp(npts:1:-1)


             ! get grid point weights for trapped particle integral

             ! Note : Here ulim and llim are the upper and lower vpar
             ! values for the valid trapped pitch angles at this theta
             ! location. In other words ulim == v_||_wfb(theta(ig)) =
             ! sqrt(max(1.0-al(ng2+1)*bmag(ig),0.0)) = yb(2*npts-1)
             ! We could probably replace ulim and llim with yb(2*npts-1)
             ! and yb(1) respectively.
             ulim = sqrt(max(1.0-bmag(ig)/bmax,0.0))
             llim = -ulim

             ! Call get_weights to find the integration weights for the v||/v grid
             ! Note we don't call this if ulim == 0 as the grid then has no extent.
             ! This situation can arise when bmag(ig) == bmax where we would find
             ! npts = 1, corresponding to just the wfb being a valid lambda. This
             ! means that at this location and the wfb has zero integration weight.
             ! Note we currently ignore the error flag and other returned values other
             ! than the weights. We should probably check this to at least warn the
             ! user if something looks badly behaved.
             if (ulim > 0.0) call get_weights (nmax, llim, ulim, yb, wb, ndiv, divmax, eflag)

             ! Calculate the weights for use in the trapped
             ! integration error estimation code.  This can be quite
             ! expensive so we would like to to only call this if the
             ! error estimation code is active (i.e. if we're going to
             ! call [[get_verr]]).
             if (npts > 1) call get_trapped_lambda_grid_error_estimate_weights(yb, llim, ulim, npts, wberr)

             ! Convert from v||/v weights to lambda
             ! weights. Essentially just sum up the postive and
             ! negtive v|| points corresponding to the same lambda
             ! point. Note we skip the v|| == 0 point currently to ensure we
             ! don't double count.
             if (npts > 1) then
                do il = ng2+1, ng2+npts-1
                   ! take into account possible asymmetry of weights about xi = 0
                   ! due to unequal # of grid points per integration interval
                   ! Whilst the v||/v grid should always be symmetric
                   ! the algorithm in [[get_weights]] could
                   ! potentially lead to asymmetry in the resulting
                   ! weights. This may be an indication that things
                   ! aren't behaving well.
                   weights(ig,il) = wb(il-ng2) + wb(2*npts-il+ng2)
                   wlterr(ig,il,:npts) = wberr(il-ng2,:) + wberr(2*npts-il+ng2,:)
                end do
             end if

             ! avoid double counting of gridpoint at vpa=0
             weights(ig,ng2+npts) = wb(npts)
             wlterr(ig,ng2+npts,:npts) = wberr(npts,:)

             deallocate(ytmp, yb, wb, wberr)
          end if
       end do

    else
       ! old (finite-difference) integration scheme

       ! Here we find the trapped pitch angle weights corresponding to a trapezium
       ! rule integration in v||/v. With this method we can write:
       !   Int(f) ~ Sum_i((f(l_{i+1}) - f(l_{i}))*(l_{i+1}-l_{i}))/2
       ! with some handling required for the upper and lower boundaries. Here we find
       ! the v||/v spacing. Note there is no factor 1/2 due to the fact the lambda
       ! weight is the sum of the (identical) positive and negative v||/v weights, which
       ! cancels this factor of 1/2. Also note that there is no weight for the bouncing
       ! lambda.
       do il = ng2+1, nlambda-1
          do ig = -ntgrid, ntgrid
             if (1.0-al(il)*bmag(ig) > -bouncefuzz &
                  .and. 1.0-al(il+1)*bmag(ig) > -bouncefuzz) &
                  then
                wwo = sqrt(max(1.0 -   al(il)*bmag(ig),0.0)) - &
                     sqrt(max(1.0 - al(il+1)*bmag(ig),0.0))
                weights(ig,il)   = weights(ig,il)   + wwo
                weights(ig,il+1) = weights(ig,il+1) + wwo
             end if
          end do
       end do
    end if
  end subroutine setup_trapped_lambda_grids

  !> Routine for getting the weights used in estimating the error on
  !> trapped particle integrals as a part of [[trap_error]]
  !>
  !> This routine is all to calculate wberr/wlterr which is just used
  !> in [[trap_error]], which in turn is only used in
  !> [[get_verr]]. This is used as a part of the [[vary_vnew]]
  !> adaptive collisionality code triggered through the
  !> diagnostics. This routine is rather expensive for high theta
  !> resolution - measured in the order of 4-5 minutes for ntheta=256.
  !> As such we might want to consider skipping this if the variable
  !> vnewk code is not active. As we do this for each npt and each
  !> theta, and npt ~ntheta/2 we might expect this to scale
  !> quadratically with ntheta. In practice it seems to scale somewhat
  !> worse than this. This may in part be due to an increase in the
  !> number of iterations required in get_weights such that this may
  !> scale closer to ntheta cubed.
  subroutine get_trapped_lambda_grid_error_estimate_weights(yb, llim, ulim, npts, wberr)
    implicit none
    real, dimension(:), intent(in) :: yb
    real, intent(in) :: llim, ulim
    integer, intent(in) :: npts
    real, dimension(:, :), intent(in out) :: wberr
    real, dimension (:), allocatable :: yberr, wberrtmp
    integer :: ix
    logical :: eflag
    integer :: ndiv, divmax, divmaxerr, ndiverr

    ! Can't get error estimate for npts = 0 or 1 as we don't have any points
    ! that we can drop.
    if (npts > 1) then
       do ix=1,npts
          if (ix == 1) then
             ! drop the first and last grid points from the integral
             allocate (yberr(2*npts-3),wberrtmp(2*npts-3))
             yberr = 0.0; wberrtmp = 0.0
             yberr = yb(2:2*npts-2)
          else if (ix == npts) then
             ! drop the vpa=0 grid point from the integral
             allocate (yberr(2*npts-2),wberrtmp(2*npts-2))
             yberr = 0.0; wberrtmp = 0.0
             yberr(:npts-1) = yb(:npts-1)
             yberr(npts:) = yb(npts+1:)
          else
             ! drop the grid points corresponding to ix and its negative from the integral
             allocate (yberr(2*npts-3),wberrtmp(2*npts-3))
             yberr = 0.0; wberrtmp = 0.0
             yberr(:ix-1) = yb(:ix-1)
             yberr(ix:2*npts-ix-2) = yb(ix+1:2*npts-ix-1)
             yberr(2*npts-ix-1:) = yb(2*npts-ix+1:)
          end if

          call get_weights (nmax, llim, ulim, yberr, wberrtmp, ndiverr, divmaxerr, eflag)

          ! insert a weight of zero into indices corresponding to ix and its conjugate
          if (ix == 1) then
             wberr(2:2*npts-2,1) = wberrtmp
          else if (ix == npts) then
             wberr(:npts-1,npts) = wberrtmp(:npts-1)
             wberr(npts+1:,npts) = wberrtmp(npts:)
          else
             wberr(:ix-1,ix) = wberrtmp(:ix-1)
             wberr(ix+1:2*npts-ix-1,ix) = wberrtmp(ix:2*npts-ix-2)
             wberr(2*npts-ix+1,ix) = wberrtmp(2*npts-ix-1)
          end if

          deallocate (yberr,wberrtmp)
       end do
    end if
  end subroutine get_trapped_lambda_grid_error_estimate_weights

  !> Determine which lambda grid points are forbidden at each
  !> theta/bmag value.
  subroutine calculate_forbidden_region(forbid)
    use theta_grid, only: ntgrid, bmag
    use mp, only: mp_abort
    implicit none
    logical, dimension(-ntgrid:, :), intent(out) :: forbid
    integer :: ig, il

    forbid = .false.

    ! Set the forbid flag
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          forbid(ig,il) = 1.0 - al(il)*bmag(ig) < -bouncefuzz
       end do
    end do

    ! Check that none of our supposedly passing particles are forbidden
    ! at any point in our domain
    do il = 1, ng2
       if ( any(forbid(:,il)) ) then
          call mp_abort("Fatal error: supposedly passing particle was trapped, in calculate_forbidden_region, in le_grids.f90", .true.)
       end if
    end do
  end subroutine calculate_forbidden_region

  !> Determine which lambda grid point bounces at this
  !> theta/bmag value.
  subroutine calculate_jend(jend)
    use theta_grid, only: ntgrid, bmag
    use mp, only: mp_abort
    implicit none
    integer, dimension(-ntgrid:), intent(out) :: jend
    integer :: ig, il

    jend = 0

    ! jend(ig) is total # of valid al grid points at each theta value

    !CMR, 1/11/2013:
    ! Above, with no trapped particles, we set: jend(ig)=   0
    ! Here, with trapped particles, we set:     jend(ig)=  il
    !  where il is the lambda index of the trapped particle bouncing at theta(ig)

    do ig = -ntgrid, ntgrid
       ! We could initialise jend to ng2 and start this loop at ng2 + 1
       ! as we assume that all lambda up to and including with ng2 satisfy
       ! 1-lambda*bmag > -bouncefuzz for all theta. This is enforced/checked
       ! in [[calculate_forbidden_region]].
       do il = 1, nlambda
          if (1.0 - al(il)*bmag(ig) > -bouncefuzz) jend(ig) = jend(ig) + 1
       end do
    end do
  end subroutine calculate_jend

  !> FIXME : Add documentation  
  subroutine xigridset
    use theta_grid, only: ntgrid, bmag
    implicit none

    integer :: ig, je, ixi

    if (.not. allocated(xi)) allocate (xi(-ntgrid:ntgrid,2*nlambda))
    if (.not. allocated(ixi_to_il)) allocate (ixi_to_il(-ntgrid:ntgrid,2*nlambda))
    if (.not. allocated(ixi_to_isgn)) allocate (ixi_to_isgn(-ntgrid:ntgrid,2*nlambda))

    ! define array 'sgn' that returns sign of vpa associated with isgn=1,2
    sgn(1) = 1
    sgn(2) = -1

    ! xi is vpa / v and goes from 1 --> -1
    xi = 0.0
    do ig = -ntgrid, ntgrid
       xi(ig,:jend(ig)) = sqrt(max(1.0 - al(:jend(ig))*spread(bmag(ig),1,jend(ig)),0.0))
       xi(ig,jend(ig)+1:2*jend(ig)-1) = -xi(ig,jend(ig)-1:1:-1)
!       do ixi = 1, nxi
!          write (*,*) 'xi', xi(ig,ixi)
!       end do
    end do

!    xi = max(xi,0.0)
!    xi = sqrt(xi)
    
    ! get mapping from ixi (which runs between 1 and 2*nlambda) and il (runs between 1 and nlambda)
    do ig = -ntgrid, ntgrid
       je = jend(ig)
       ! if no trapped particles
       if (je == 0) then
          do ixi = 1, 2*nlambda
             if (ixi > nlambda) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*nlambda - ixi + 1
             else
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi
             end if
          end do
       else
!CMR, 1/11/2013:
! Sketch of how ixi=>il mapping is arranged
!===============================================================================================
!  ixi=   1, ... , je-1, je, je+1, ... , 2je-1, || 2je, 2je+1, ... , nl+je, nl+je+1, ... ,  2nl
!   il=   1, ... , je-1, je, je-1, ... ,     1, ||  je,  je+1, ... ,    nl,      nl, ... , je+1
! isgn=   1, ... ,    1,  2,    2, ... ,     2, ||   1,     1, ... ,     1,       2, ... ,    2
!         particles passing through             ||  je,   + forbidden trapped velocity space
!         nb only need one isigma for je, as v||=0 at the bounce point
!===============================================================================================

          do ixi = 1, 2*nlambda
             if (ixi >= nlambda + je + 1) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*nlambda + je + 1 - ixi
             else if (ixi >= 2*je) then
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi - je
             else if (ixi >= je) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*je - ixi
             else
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi
             end if
          end do
       end if
    end do

  end subroutine xigridset

  !> FIXME : Add documentation  
  subroutine lagrange_coefs (ig, nodes, lfac, xloc)
    use theta_grid, only: bmag, bmax
    implicit none

    integer, intent (in) :: ig
    real, dimension (:), intent (in) :: nodes
    real, dimension (:,:,:), intent (out) :: lfac
    real, dimension (:), intent (out) :: xloc

    integer :: ix, il, j, icnt, tsize

    xloc(1) = -sqrt(max(1.0-bmag(ig)/bmax,0.0))
    xloc(2*nterp-1) = sqrt(max(1.0-bmag(ig)/bmax,0.0))
    do ix=2,2*nterp-1
       xloc(ix) = xloc(1) + (ix-1)*(xloc(2*nterp-1)-xloc(1))/(2*nterp-2)
    end do

    lfac = 0.
    icnt = 0

    tsize = size(nodes)

    do il=ng2+1,nlambda
       if (.not. forbid(ig,il)) then
          lfac(il,:,:) = 1.0
          icnt = icnt + 1
          if (icnt .eq. tsize-icnt+1) then
             lfac(il,:,1) = 0.
             do j=1,tsize
                if (j /= icnt) lfac(il,:,2) = lfac(il,:,2)*(xloc - nodes(j))/(nodes(icnt)-nodes(j))
             end do
          else
             do j=1,tsize
                if (j /= icnt) lfac(il,:,2) = lfac(il,:,2)*(xloc - nodes(j))/(nodes(icnt)-nodes(j))
                if (j /= tsize-icnt+1) lfac(il,:,1) = lfac(il,:,1)*(xloc - nodes(j))/(nodes(tsize-icnt+1)-nodes(j))
             end do
          end if
       end if
    end do

  end subroutine lagrange_coefs

  !> FIXME : Add documentation  
  subroutine stop_message (message)
!<DD>WARNING THIS COULD CAUSE A HANG IF NOT CALLED BY ALL PROCS
!JAB: print an error message and end program
    use file_utils, only: error_unit
    use mp, only: proc0, finish_mp
    implicit none
    character(*), intent (in) :: message
    integer :: ierr

    if (proc0) then
       ierr = error_unit()
       write (unit=ierr, fmt='(a)') message
    end if
    call finish_mp
    stop
  end subroutine stop_message

  !> FIXME : Add documentation
  !!
  !! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
  !! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
  !! is over x-y space
  subroutine integrate_volume_c (g, total, all)
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: nproc, sum_reduce, sum_allreduce

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    real :: fac
    integer :: is, il, ie, ik, iglo, isgn

    !Initialise to zero
    total = 0.

    !Do integral over local x-y space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Pick the weighting factor
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       !For both signs of vpar do sum
       !May be more efficient to move ign loop above iglo loop (good for total but bad for g memory access)
       do isgn = 1, 2
          total(:, il, ie, isgn, is) = total(:, il, ie, isgn, is) + &
               fac*g(:,isgn,iglo)
       end do
    end do

    !Do we really need this if statement?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_volume_c

  !> FIXME : Add documentation
  !!
  !! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
  !! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
  !! is over x-y space  
  subroutine integrate_volume_r (g, total, all)
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: nproc,sum_reduce, sum_allreduce

    implicit none

    real, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    real :: fac
    integer :: is, il, ie, ik, iglo, isgn

    !Initialise to zero
    total = 0.

    !Do integral over local x-y space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Pick the weighting factor
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       !For both signs of vpar do sum
       !May be more efficient to move ign loop above iglo loop (good for total but bad for g memory access)
       do isgn = 1, 2
          total(:, il, ie, isgn, is) = total(:, il, ie, isgn, is) + &
               fac*g(:,isgn,iglo)
       end do
    end do

    !Do we really need this if statement?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_volume_r

  !> Calculates and returns toroidal momentum flux as a function
  !! of vpar and theta
  subroutine get_flux_vs_theta_vs_vpa (f, vflx, dealloc)

    use theta_grid, only: ntgrid, bmag
    use species, only: nspec

    implicit none
    logical, intent(in), optional :: dealloc
    real, dimension (-ntgrid:,:,:,:,:), intent (in) :: f
    real, dimension (-ntgrid:,:,:), intent (out) :: vflx

    real, dimension (:,:,:), allocatable :: favg
    real, dimension (:,:), allocatable, save :: vpa1d
    real, dimension (:,:,:), allocatable, save :: hermp1d
    real, dimension (:,:,:,:,:), allocatable, save :: vpapts
    real, dimension (:,:,:,:,:,:), allocatable, save :: hermp

    integer :: is, il, ie, ig, iv
    integer :: norder

    if(present(dealloc))then
       if(allocated(vpa1d)) deallocate(vpa1d)
       if(allocated(hermp1d)) deallocate(hermp1d)
       if(allocated(vpapts)) deallocate(vpapts)
       if(allocated(hermp)) deallocate(hermp)
       return
    endif

    norder = min(negrid, nlambda)/2

    allocate (favg(-ntgrid:ntgrid,nspec,0:norder-1))

    if (.not. allocated(vpapts)) then
       allocate (vpa1d(negrid*nlambda,nspec))
       allocate (hermp1d(negrid*nlambda,0:norder-1,nspec))
       allocate (vpapts(-ntgrid:ntgrid,nlambda,negrid,2,nspec))
       allocate (hermp(-ntgrid:ntgrid,nlambda,negrid,2,0:norder-1,nspec))
       vpapts = 0.0 ; hermp = 0.0 ; vpa1d = 0.0 ; hermp1d = 0.0

       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                vpapts(ig,il,ie,1,:) = sqrt(energy(ie,:)*max(0.0, 1.0-al(il)*bmag(ig)))
                vpapts(ig,il,ie,2,:) = -vpapts(ig,il,ie,1,:)
             end do
          end do
       end do

       do iv = 1, negrid*nlambda
          vpa1d(iv,:) = sqrt(energy(negrid,:))*(1. - 2.*(iv-1)/real(negrid*nlambda-1))
       end do

       do is = 1,nspec
         call get_hermite_polynomials (vpa1d(:,is), hermp1d(:,:,is))
         call get_hermite_polynomials (vpapts(:,:,:,:,is), hermp(:,:,:,:,:,is))
       end do
    end if

    favg = 0.
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                favg(ig,is,:) = favg(ig,is,:) &
                     +w(ie,is)*wl(ig,il)*(hermp(ig,il,ie,1,:,is)*f(ig,il,ie,1,is) &
                     +hermp(ig,il,ie,2,:,is)*f(ig,il,ie,2,is))
             end do
          end do
       end do
    end do

    do is = 1, nspec
       do iv = 1, negrid*nlambda
          do ig = -ntgrid, ntgrid
             vflx(ig,iv,is) = sum(favg(ig,is,:)*hermp1d(iv,:,is))*exp(-vpa1d(iv,is)**2)
          end do
       end do
    end do

    deallocate (favg)
    
  end subroutine get_flux_vs_theta_vs_vpa

  !> Returns Gn = Hn / sqrt(2^n n!) / pi^(1/4),
  !! where Hn are the hermite polynomials
  !! i.e. int dx Gm * Gn exp(-x^2) = 1
  subroutine get_hermite_polynomials_4d (xptsdum, hpdum)
    use constants, only: pi
    use theta_grid, only: ntgrid

    implicit none

    real, dimension (-ntgrid:,:,:,:), intent (in)   :: xptsdum
    real, dimension (-ntgrid:,:,:,:,0:), intent (out) :: hpdum

    integer :: j
    double precision, dimension (:,:,:,:), allocatable :: hp1, hp2, hp3

    hpdum = 0.0

    allocate (hp1(-ntgrid:ntgrid,nlambda,negrid,2))
    allocate (hp2(-ntgrid:ntgrid,nlambda,negrid,2))
    allocate (hp3(-ntgrid:ntgrid,nlambda,negrid,2))

    hp1 = real(1.0,kind(hp1(0,1,1,1)))
    hp2 = real(0.0,kind(hp2(0,1,1,1)))

    hpdum(:,:,:,:,0) = 1.0

    do j=1, size(hpdum,5)-1
       hp3 = hp2
       hp2 = hp1
       hp1 = sqrt(2./j)*xptsdum*hp2 - sqrt(real(j-1)/j)*hp3
       hpdum(:,:,:,:,j) = hp1
    end do

    hpdum = hpdum/pi**(0.25)

    deallocate (hp1,hp2,hp3)

  end subroutine get_hermite_polynomials_4d

  !> Returns Gn = Hn / sqrt(2^n n!) / pi^(1/4),
  !! where Hn are the hermite polynomials
  !! i.e. int dx Gm * Gn exp(-x^2) = 1
  subroutine get_hermite_polynomials_1d (xptsdum, hpdum)
    use constants, only: pi
    implicit none

    real, dimension (:), intent (in)   :: xptsdum
    real, dimension (:,0:), intent (out) :: hpdum

    integer :: j
    double precision, dimension (:), allocatable :: hp1, hp2, hp3

    hpdum = 0.0

    allocate (hp1(size(xptsdum)))
    allocate (hp2(size(xptsdum)))
    allocate (hp3(size(xptsdum)))

    hp1 = real(1.0,kind(hp1(1)))
    hp2 = real(0.0,kind(hp2(1)))

    hpdum(:,0) = 1.0

    do j=1, size(hpdum,2)-1
       hp3 = hp2
       hp2 = hp1
       hp1 = sqrt(2./j)*xptsdum*hp2 - sqrt(real(j-1)/j)*hp3
       hpdum(:,j) = hp1
    end do

    hpdum = hpdum/pi**(0.25)

    deallocate (hp1,hp2,hp3)

  end subroutine get_hermite_polynomials_1d

  !> Returns probabilist's Hermite polynomials, He_n(x)
  elemental function hermite_prob (n, x)
    implicit none
    integer, intent (in) :: n
    real, intent (in) :: x
    integer :: k
    double precision :: hermite_prob, p, p1, p2

    p1 = x
    p2 = dble(1.0)

    if (n==0) then
       hermite_prob = p2
       return
    else if (n==1) then
       hermite_prob = p1
       return
    end if

    do k=2, n
       p = x*p1 - (k-1)*p2
       p2 = p1
       p1 = p
    end do

    hermite_prob = p

  end function hermite_prob

  !> FIXME : Add documentation
  subroutine init_map (use_lz_layout, use_e_layout, use_le_layout, test)
    use mp, only: finish_mp, proc0
    use redistribute, only: report_map_property
    implicit none

    logical, intent (in) :: use_lz_layout, use_e_layout, use_le_layout, test

    ! initialize maps from g_lo to lz_lo, e_lo, and/or le_lo

    if (use_lz_layout) then
       ! init_lambda_layout is called in redistribute
       call init_lambda_redistribute_local

       if (test) then
          if (proc0) print *, '=== Lambda map property ==='
          call report_map_property (lambda_map)
       end if
    end if

    if (use_e_layout) then
       ! init_energy_layout is called in redistribute
       call init_energy_redistribute_local

       if (test) then
          if (proc0) print *, '=== Energy map property ==='
          call report_map_property (energy_map)
       end if
    end if

    if (use_le_layout) then
       call init_g2le_redistribute_local
       if (test) call check_g2le

    end if

    call init_g2gf(test)

    if (test) then
       if (proc0) print *, 'init_map done'
    end if

  end subroutine init_map

  !> FIXME : Add documentation  
  subroutine init_g2gf(test)
    implicit none

    logical, intent (in) :: test

    call init_g2gf_redistribute
    
    if (test) call check_g2gf

  end subroutine init_g2gf

  !> Constructs the redistribute mapping from the g_lo data
  !> decomposition to the le_lo decomposition.
  subroutine init_g2le_redistribute_local
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_le_layouts
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ig_idx, isign_idx
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
    use sorting, only: quicksort
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list, sort_list
!    type (index_list_type), dimension(0:nproc-1) :: old_list !Testing
    integer, dimension (0:nproc-1) :: nn_to, nn_from
!    integer, dimension (0:nproc-1) :: nn_to_old, nn_from_old !Testing
    integer, dimension (3) :: from_low, from_high
    integer, dimension (3) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ile
    integer :: ie
    integer :: n, ip, je
    integer :: ile_bak, il0

    !Early exit if possible
    if (leinit) return
    leinit = .true.

    !Setup the le layout object (le_lo)
    call init_le_layouts (ntgrid, naky, ntheta0, nspec)

    !Initialise the data counters
    nn_to = 0
    nn_from = 0

    !First count the data to be sent | g_lo-->le_lo
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo = g_lo%llim_proc,g_lo%ulim_alloc
          !Get le_lo idx for ig=-ntgrid
          ile=idx(le_lo,-ntgrid,ik_idx(g_lo,iglo),&
               it_idx(g_lo,iglo),is_idx(g_lo,iglo))

          !Loop over remaining dimensions, note ile is independent of isign
          !so add two to count
          do ig=-ntgrid,ntgrid
             !Increment the data sent counter for this proc
             nn_from(proc_id(le_lo,ile))=nn_from(proc_id(le_lo,ile))+2

             !Increment ile
             ile=ile+1
          enddo
       enddo
    endif

    !Now count how much data to receive | le_lo<--g_lo
    !Protect against procs with no data
    if(le_lo%ulim_proc.ge.le_lo%llim_proc)then
       do ile=le_lo%llim_proc,le_lo%ulim_alloc
          !Loop over local dimensions, adding 2 to account for each sign
          do ie=1,g_lo%negrid
             do il=1,g_lo%nlambda
                !Get index
                iglo=idx(g_lo,ik_idx(le_lo,ile),it_idx(le_lo,ile),&
                     il,ie,is_idx(le_lo,ile))

                !Increment the data to receive counter
                nn_to(proc_id(g_lo,iglo))=nn_to(proc_id(g_lo,iglo))+2
             enddo
          enddo
       enddo
    endif

    !Now allocate storage for index arrays
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip)))
          !For sorting message order later
          allocate (sort_list(ip)%first(nn_to(ip)))
       end if
    end do

    !Reinitialise counters
    nn_to = 0
    nn_from = 0

    !First fill in sending indices, these define the message order
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          !Get ile for ig=-ntgrid
          ile=idx(le_lo,-ntgrid,ik_idx(g_lo,iglo),&
               it_idx(g_lo,iglo),is_idx(g_lo,iglo))
          ile_bak=ile

          !Loop over sign
          do isign=1,2
             do ig=-ntgrid,ntgrid
                !Get proc id
                ip=proc_id(le_lo,ile)

                !Increment procs message counter
                n=nn_from(ip)+1
                nn_from(ip)=n

                !Store indices
                from_list(ip)%first(n)=ig
                from_list(ip)%second(n)=isign
                from_list(ip)%third(n)=iglo

                !Increment ile
                ile=ile+1
             enddo

             !Restore ile
             ile=ile_bak
          enddo
       enddo
    endif

    !Now fill in the receiving indices, these must match message data order
    !Protect against procs with no data
    if(le_lo%ulim_proc.ge.le_lo%llim_proc)then
       do ile=le_lo%llim_proc,le_lo%ulim_alloc
          !Get ig index
          ig=ig_idx(le_lo,ile)
 
          !Loop over local dimensions,Whilst ile is independent of sign this information
          !is in lambda so loop over sign included here
          do ie=1,g_lo%negrid
             do isign=1,2
                do il0=1,g_lo%nlambda
                   !Pick correct extended lambda value
                   je=jend(ig)
                   il=il0
                   if (je.eq.0) then
                      if (isign.eq.2) il=2*g_lo%nlambda+1-il
                   else
                      if(il.eq.je) then
                         if(isign.eq.1) il=2*je
                      else if(il.gt.je) then
                         if(isign.eq.1) then
                            il=il+je
                         else
                            il=2*g_lo%nlambda+1-il+je
                         endif
                      else
                         if(isign.eq.2) il=2*je-il
                      endif
                   endif
                   
                   !Get iglo value
                   iglo=idx(g_lo,ik_idx(le_lo,ile),it_idx(le_lo,ile),&
                        il0,ie,is_idx(le_lo,ile))

                   !Get proc_id
                   ip=proc_id(g_lo,iglo)

                   !Increment counter
                   n=nn_to(ip)+1
                   nn_to(ip)=n

                   !Store indices
                   to_list(ip)%first(n)=il
                   to_list(ip)%second(n)=ie
                   to_list(ip)%third(n)=ile

                   !Store sorting index
                   sort_list(ip)%first(n)=ig+ntgrid-1+(2*ntgrid+1)*(isign-1+2*(iglo-g_lo%llim_world))
                enddo
             enddo
          enddo
       enddo
    endif

    !Now sort receive indices into message order
    do ip=0,nproc-1
       if(nn_to(ip)>0) then
          !Apply quicksort
          call quicksort(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second,to_list(ip)%third)
       endif
    enddo

    !Now setup array range values
    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    to_low = le_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = negrid + 1  ! TT: just followed convention with +1.
    ! TT: It may be good to avoid bank conflict.
    to_high(3) = le_lo%ulim_alloc

    !Create g2le redist object
    call init_redist (g2le, 'c', to_low, to_high, to_list, from_low, from_high, from_list)

    !Deallocate lists
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list (sort_list)
  end subroutine init_g2le_redistribute_local

  !> FIXME : Add documentation
  subroutine check_g2le
    use file_utils, only: error_unit
    use mp, only: finish_mp, iproc, proc0
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use redistribute, only: gather, scatter, report_map_property

    implicit none

    integer :: iglo, ile, ig, isgn, ik, it, il, ie, is, ierr
    integer :: ixi, je
    complex, dimension (:,:,:), allocatable :: gtmp, letmp

    if (proc0) then
       ierr = error_unit()
    else
       ierr = 6
    end if

    ! report the map property
    if (proc0) print *, '=== g2le map property ==='
    call report_map_property (g2le)

    allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (letmp(nlambda*2+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    letmp = 0.

    ! gather check
    gtmp = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = rule(ig,isgn,ik,it,il,ie,is)
          end do
       end do
    end do
    call gather (g2le, gtmp, letmp)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ig = ig_idx(le_lo,ile)
       je = jend(ig)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       is = is_idx(le_lo,ile)
       do ie = 1, negrid
          do ixi = 1, 2*nlambda
             isgn = ixi_to_isgn(ig,ixi)
             il = ixi_to_il(ig,ixi)
             if (int(real(letmp(ixi,ie,ile))) /= rule(ig,isgn,ik,it,il,ie,is)) &
                  write (ierr,'(a,8i6)') 'ERROR: gather by g2le broken!', iproc
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2le gather check done'

    ! scatter check
    gtmp = 0.0
    call scatter (g2le, letmp, gtmp)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             if (gtmp(ig,isgn,iglo) /= rule(ig,isgn,ik,it,il,ie,is)) &
                  write (ierr,'(a,i6)') 'ERROR: scatter by g2le broken!', iproc
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2le scatter check done'

    deallocate (gtmp,letmp)

!    call finish_mp
!    stop

  contains

    function rule (ig, isgn, ik, it, il, ie, is)
      integer, intent (in) :: ig, isgn, ik, it, il, ie, is
      integer :: rule
      rule = ig + isgn + ik + it + il + ie + is  ! make whatever you want
    end function rule

  end subroutine check_g2le

  !> FIXME : Add documentation  
  subroutine init_g2gf_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_gf_layouts
    use gs2_layouts, only: g_lo, gf_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ig_idx, isign_idx
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
    use redistribute, only: index_list_type, init_redist, delete_list
    use sorting, only: quicksort

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list, sort_list, bak_sort_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (6) :: to_low, to_high
    integer :: iglo, il, igf
    integer :: ie, is
    integer :: n, ip

    !Early exit if possible
    if (gfinit) return
    gfinit = .true.

    !Setup the le layout object (gf_lo)
    call init_gf_layouts (ntgrid, naky, ntheta0, negrid, nlambda, nspec)

    !Initialise the data counters
    nn_to = 0
    nn_from = 0

    !First count the data to be sent | g_lo-->gf_lo
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo = g_lo%llim_proc,g_lo%ulim_alloc
          
          igf = idx(gf_lo, ik_idx(g_lo, iglo), it_idx(g_lo, iglo))
          nn_from(proc_id(gf_lo,igf))=nn_from(proc_id(gf_lo,igf))+1
       enddo
    endif

    !Now count how much data to receive | gf_lo<--g_lo
    !Protect against procs with no data
    if(gf_lo%ulim_proc.ge.gf_lo%llim_proc)then
       do igf=gf_lo%llim_proc,gf_lo%ulim_alloc
          do is=1,g_lo%nspec
             do ie=1,g_lo%negrid
                do il=1,g_lo%nlambda
                   iglo = idx(g_lo, ik_idx(gf_lo, igf), it_idx(gf_lo, igf), il, ie, is)
                   !Increment the data to receive counter
                   nn_to(proc_id(g_lo,iglo))=nn_to(proc_id(g_lo,iglo))+1
                enddo
             enddo
          enddo
       enddo
    endif


    !Now allocate storage for index arrays
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip)))
          allocate (to_list(ip)%fourth (nn_to(ip)))
          allocate (to_list(ip)%fifth (nn_to(ip)))
          allocate (to_list(ip)%sixth (nn_to(ip)))
          allocate (sort_list(ip)%first (nn_to(ip)))
          allocate (bak_sort_list(ip)%first (nn_to(ip)))
       end if
    end do

    !Reinitialise counters
    nn_to = 0
    nn_from = 0

    !First fill in sending indices, these define the message order
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          igf=idx(gf_lo,ik_idx(g_lo,iglo),it_idx(g_lo,iglo))

          !Get proc id
          ip=proc_id(gf_lo,igf)

          n=nn_from(ip)+1
          nn_from(ip)=n
          
          !Store indices (we are sending ntgridtotal*isign sized messages each time at a time)
          from_list(ip)%first(n)=-gf_lo%ntgrid
          from_list(ip)%second(n)=1
          from_list(ip)%third(n)=iglo
          
       enddo
    endif

    !Now fill in the receiving indices, these must match message data order
    !Protect against procs with no data
    if(gf_lo%ulim_proc.ge.gf_lo%llim_proc)then
       do igf=gf_lo%llim_proc,gf_lo%ulim_alloc
          do il=1,g_lo%nlambda
             do ie=1,g_lo%negrid
                do is=1,g_lo%nspec
                   !Get iglo value
                   iglo = idx(g_lo,ik_idx(gf_lo,igf),it_idx(gf_lo,igf),il,ie,is)
                   !Get proc_id
                   ip=proc_id(g_lo,iglo)
                                         
                   !Increment counter
                   n=nn_to(ip)+1
                   nn_to(ip)=n
                   
                   !Store indices
                   !We are sending messages of size ntgridtotal*isign
                   to_list(ip)%first(n)=-gf_lo%ntgrid
                   to_list(ip)%second(n)=1
                   to_list(ip)%third(n)=is
                   to_list(ip)%fourth(n)=ie
                   to_list(ip)%fifth(n)=il
                   to_list(ip)%sixth(n)=igf

                   sort_list(ip)%first(n) = iglo-g_lo%llim_world
                enddo
             enddo
          enddo
       enddo
    endif

    do ip=0,nproc-1
       if(associated(sort_list(ip)%first)) then
          bak_sort_list(ip)%first(:) = sort_list(ip)%first(:)
       end if
    end do

    !Now sort receive indices into message order
    do ip=0,nproc-1
       if(nn_to(ip)>0) then
          !Apply quicksort
          !AJ I'm calling quicksort twice to save me writing a six array version of quicksort and insertsort.
          !AJ This is just lazyness and probably should be corrected at some point.
          call quicksort(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second,to_list(ip)%third)
          call quicksort(nn_to(ip),bak_sort_list(ip)%first,to_list(ip)%fourth,to_list(ip)%fifth,to_list(ip)%sixth)
       endif
    enddo


    !Now setup array range values
    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    to_low(1) = -ntgrid
    to_low(2) = 1
    to_low(3) = 1
    to_low(4) = 1
    to_low(5) = 1
    to_low(6) = gf_lo%llim_proc

    to_high(1) = ntgrid
    to_high(2) = 2
    to_high(3) = nspec
    to_high(4) = negrid
    to_high(5) = nlambda
    to_high(6) = gf_lo%ulim_alloc

    !Create g2gf redist object
    call init_redist (g2gf, 'c', to_low, to_high, to_list, from_low, from_high, from_list)

    !Deallocate lists
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list (sort_list)
    call delete_list (bak_sort_list)

  end subroutine init_g2gf_redistribute

  !> FIXME : Add documentation  
  subroutine check_g2gf
    use file_utils, only: error_unit
    use mp, only: finish_mp, iproc, proc0
    use species, only : nspec
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gf_lo
    use gs2_layouts, only: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use redistribute, only: gather, scatter, report_map_property

    implicit none

    integer :: iglo, ig, isgn, ik, it, il, ie, igf, is, ierr
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:,:,:), allocatable :: gftmp

    if (proc0) then
       ierr = error_unit()
    else
       ierr = 6
    end if

    allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gftmp(-ntgrid:ntgrid, 2, nspec, negrid, nlambda, gf_lo%llim_proc:gf_lo%ulim_alloc))

    ! gather check
    gtmp = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = rule(ig,isgn,ik,it,il,ie,is)
          end do
       end do
    end do
    call gather (g2gf, gtmp, gftmp, ntgrid)
    do igf = gf_lo%llim_proc, gf_lo%ulim_proc
       ik = ik_idx(gf_lo,igf)
       it = it_idx(gf_lo,igf)
       do is = 1, nspec
          do il = 1,nlambda
             do ie = 1,negrid
                do isgn = 1,2
                   do ig=-ntgrid,ntgrid
                      if (int(real(gftmp(ig,isgn,is,ie,il,igf))) /= rule(ig,isgn,ik,it,il,ie,is)) &
                           write (ierr,'(a,8i6)') 'ERROR: gather by g2gf broken!', iproc
                   end do
                end do
             end do
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2gf gather check done'

    ! scatter check
    gtmp = 0.0
    call scatter (g2gf, gftmp, gtmp, ntgrid)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             if (gtmp(ig,isgn,iglo) /= rule(ig,isgn,ik,it,il,ie,is)) &
                  write (ierr,'(a,i6)') 'ERROR: scatter by g2gf broken!', iproc
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2gf scatter check done'

    deallocate (gtmp,gftmp)

!    call finish_mp
!    stop

  contains
    !> FIXME : Add documentation
    function rule (ig, isgn, ik, it, il, ie, is)
      integer, intent (in) :: ig, isgn, ik, it, il, ie, is
      integer :: rule
      rule = ig + isgn + ik + it + il + ie + is  ! make whatever you want
    end function rule
  end subroutine check_g2gf

  !> Computes the redistribute mapping from the g_lo data
  !> decomposition to the lz_lo decomposition
  subroutine init_lambda_redistribute_local
    use mp, only: nproc, barrier
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_lambda_layouts
    use gs2_layouts, only: g_lo, lz_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx, ig_idx,isign_idx
    use redistribute, only: index_list_type, init_redist, delete_list
    use sorting, only: quicksort
    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list, sort_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ilz, il0
    integer :: je, ilz_bak
    integer :: n, ip

    !Early exit if possible
    if (lzinit) return
    lzinit = .true.

    !Setup the lambda layout object (lz_lo)
    call init_lambda_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)

    !Initialise data counters
    nn_to = 0
    nn_from = 0

    !First count the data to be send | g_lo-->lz_lo
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo = g_lo%llim_proc, g_lo%ulim_alloc
          !Get lz_lo idx for ig=-ntgrid
          ilz=idx(lz_lo,-ntgrid,ik_idx(g_lo,iglo), &
               it_idx(g_lo,iglo),ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
          
          !Loop over other local dimensions, noting that ig->ig+1 ==> ilz->ilz+1
          !Note that ilz is independent ofi sign so we just add two pieces of data per point instead
          do ig=-ntgrid,ntgrid
             !Increment the data sent counter for the processor with ilz
             nn_from(proc_id(lz_lo,ilz))=nn_from(proc_id(lz_lo,ilz))+2
             
             !Increment ilz
             ilz=ilz+1
          enddo
       enddo
    endif

    !Now count how much data to receive | lz_lo<--g_lo
    !Protect against procs with no data
    if(lz_lo%ulim_proc.ge.lz_lo%llim_proc)then
       do ilz = lz_lo%llim_proc, lz_lo%ulim_alloc
          do il=1,g_lo%nlambda
             !Get iglo
             iglo=idx(g_lo,ik_idx(lz_lo,ilz),it_idx(lz_lo,ilz),&
                  il,ie_idx(lz_lo,ilz),is_idx(lz_lo,ilz))
             
             !Increment the data to receive count for this proc
             !Note we increment by two due to independence of isign
             nn_to(proc_id(g_lo,iglo))=nn_to(proc_id(g_lo,iglo))+2
          enddo
       enddo
    endif

   !Now allocate storage for index arrays
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          !For sorting messsages later
          allocate (sort_list(ip)%first(nn_to(ip)))
       end if
    end do

    !Reinitialise counters
    nn_to = 0
    nn_from = 0

    !First fill in the sending indices, these define the message order
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          !Convert to ilz for ig=-ntgrid
          ilz=idx(lz_lo,-ntgrid,ik_idx(g_lo,iglo), &
            it_idx(g_lo,iglo),ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
          
          !Store backup of ilz value
          ilz_bak=ilz
          
          !Loop over other local dimensions
          do isign=1,2
             do ig=-ntgrid,ntgrid
                !Get proc id
                ip=proc_id(lz_lo,ilz)
                
                !Increment procs message counter
                n=nn_from(ip)+1
                nn_from(ip)=n
                
                !Store indices
                from_list(ip)%first(n)=ig
                from_list(ip)%second(n)=isign
                from_list(ip)%third(n)=iglo
                
                !Increment ilz
                ilz=ilz+1
             enddo
             
          !Restore ilz
             ilz=ilz_bak
          enddo
       enddo
    endif

    !Now fill in the receiving indices, these must match the message order (through sorting later)
    !NOTE: Not all procs have data in lz_lo so protect against this 
    if(lz_lo%ulim_proc.ge.lz_lo%llim_proc)then
       do ilz=lz_lo%llim_proc,lz_lo%ulim_alloc
          ig=ig_idx(lz_lo,ilz)
          !Whilst lz_lo is independent of sign we actually have a lambda dimension double that of other layouts, which is how the sign dependent information is stored, so loop over sign here.
          do isign=1,2
             !Should the upper limit actually be max(nlambda,ng2+1)?
             do il0=1,g_lo%nlambda
                je=jend(ig)
                !Pick the correct value of il
                il=il0
                if (je.eq.0) then
                   if (isign.eq.2) il=2*g_lo%nlambda+1-il
                else
                   if(il.eq.je) then
                      if(isign.eq.1) il=2*je
                   else if(il.gt.je) then
                      if(isign.eq.1) then
                         il=il+je
                      else
                         il=2*g_lo%nlambda+1-il+je
                      endif
                   else
                      if(isign.eq.2) il=2*je-il
                   endif
                endif
                
                !Get iglo value. Note we use il0 and not il
                iglo=idx(g_lo,ik_idx(lz_lo,ilz),it_idx(lz_lo,ilz),&
                     il0,ie_idx(lz_lo,ilz),is_idx(lz_lo,ilz))
                
                !Get proc id
                ip=proc_id(g_lo,iglo)
                
                !Increment counter
                n=nn_to(ip)+1
                nn_to(ip)=n
                
                !Store indices
                to_list(ip)%first(n)=il
                to_list(ip)%second(n)=ilz
                
                !Store sorting index
                sort_list(ip)%first(n)=ig+ntgrid-1+(2*ntgrid+1)*(isign-1+2*(iglo-g_lo%llim_world))
             enddo
          enddo
       enddo
    endif
    
    !Now sort receive indices into message order
    do ip=0,nproc-1
       !Only worry about the cases where we're receiving data
       if(nn_to(ip)>0) then
          !Apply quicksort
          call quicksort(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second)
       endif
    enddo

    !Now setup array range values
    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = lz_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = lz_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    !Create lambda map redistribute objects
    call init_redist (lambda_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    !Deallocate the list objects
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list (sort_list)

  end subroutine init_lambda_redistribute_local

  !> Constructs the redistribute mapping from the g_lo data
  !> decomposition to the e_lo decomposition
  subroutine init_energy_redistribute_local
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky,ntheta0
    use gs2_layouts, only: init_energy_layouts
    use gs2_layouts, only: g_lo, e_lo, ie_idx, ik_idx, it_idx, il_idx, is_idx, ig_idx, isign_idx
    use gs2_layouts, only:proc_id, idx, idx_local
    use redistribute, only: index_list_type, init_redist, delete_list
    use sorting, only: quicksort

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list,sort_list  
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, ie, ielo
    integer :: n, ip

    !Early exit if possible
    if (einit) return
    einit = .true.

    !Initialise e_lo layout object
    call init_energy_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)

    !Initialise counters to zero
    nn_to = 0
    nn_from = 0

    !First count how much data to send | g_lo-->e_lo
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo = g_lo%llim_proc, g_lo%ulim_alloc
          !Get e_lo idx for ig=-ntgrid and isign=1
          ielo=idx(e_lo,-ntgrid,1,ik_idx(g_lo,iglo), &
               it_idx(g_lo,iglo),il_idx(g_lo,iglo),is_idx(g_lo,iglo))

          !Loop over other local dimensions, noting that ig->ig+1 ==> ielo->ielo+1
          do isign = 1,2
             do ig=-ntgrid,ntgrid
                !Increment the data sent counter for the processor with ielo
                nn_from(proc_id(e_lo,ielo))=nn_from(proc_id(e_lo,ielo))+1

                !Increment ielo
                ielo=ielo+1
             enddo
          enddo
       enddo
    endif

    !Now count how much data to receive | e_lo<--g_lo
    !Protect against procs with no data
    if(e_lo%ulim_proc.ge.e_lo%llim_proc)then
       do ielo = e_lo%llim_proc, e_lo%ulim_alloc
          do ie=1,g_lo%negrid
             !Get iglo
             iglo=idx(g_lo,ik_idx(e_lo,ielo),it_idx(e_lo,ielo),&
                  il_idx(e_lo,ielo),ie,is_idx(e_lo,ielo))

             !Increment the data to receive count for this proc
             nn_to(proc_id(g_lo,iglo))=nn_to(proc_id(g_lo,iglo))+1
          enddo
       enddo
    endif

    !Now we've done counting allocate index arrays
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          !For sorting messages later
          allocate (sort_list(ip)%first(nn_to(ip)))
       end if
    end do

    !Reinitialise counters to zero
    nn_to = 0
    nn_from = 0

    !First fill in the sending indices, these define the message order
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          !Convert to ielo for ig=-ntgrid and isign=1
          ielo=idx(e_lo,-ntgrid,1,ik_idx(g_lo,iglo), &
               it_idx(g_lo,iglo),il_idx(g_lo,iglo),is_idx(g_lo,iglo))

          !Loop over other local dimensions
          do isign=1,2
             do ig=-ntgrid,ntgrid
                !Get proc id
                ip=proc_id(e_lo,ielo)

                !Increment procs message counter
                n=nn_from(ip)+1
                nn_from(ip)=n

                !Store indices
                from_list(ip)%first(n)=ig
                from_list(ip)%second(n)=isign
                from_list(ip)%third(n)=iglo

                !Increment ielo
                ielo=ielo+1
             enddo
          enddo
       enddo
    endif

    !Now fill in receive indices, these must match message data order (achieved through sorting later)
    !Protect against procs with no data
    if(e_lo%ulim_proc.ge.e_lo%llim_proc)then
       do ielo=e_lo%llim_proc,e_lo%ulim_alloc
          !Get indices used for creating sort index
          ig=ig_idx(e_lo,ielo)
          isign=isign_idx(e_lo,ielo)

          do ie=1,g_lo%negrid
             !Get iglo index
             iglo=idx(g_lo,ik_idx(e_lo,ielo),it_idx(e_lo,ielo),&
                  il_idx(e_lo,ielo),ie,is_idx(e_lo,ielo))

             !Get proc id
             ip=proc_id(g_lo,iglo)

             !Increment procs data counter
             n=nn_to(ip)+1
             nn_to(ip)=n

             !Store message indices
             to_list(ip)%first(n)=ie
             to_list(ip)%second(n)=ielo

             !Store index for sorting
             sort_list(ip)%first(n)=ig+ntgrid-1+(2*ntgrid+1)*(isign-1+(iglo-g_lo%llim_world)*2)
          enddo
       enddo
    endif

    !Now sort receive indices into message order
    do ip=0,nproc-1
       !Only worry about cases where we're receiving data
       if(nn_to(ip)>0) then
          !Sort based on quicksort
          call quicksort(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second)
       endif
    enddo

    !Now setup array range values
    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = e_lo%llim_proc

    to_high(1) = negrid+1
    to_high(2) = e_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    !Create energy map redistribute object
    call init_redist (energy_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    !Deallocate the list objects 
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list (sort_list)

  end subroutine init_energy_redistribute_local

  !> Used for testing
  !! takes as input an array using g_lo and
  !! writes it to a .distmp output file
  subroutine write_mpdist (dist, extension, last)
    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: ie_idx, il_idx, idx_local, proc_id
    use theta_grid, only: ntgrid, bmag, theta
    use kt_grids, only: aky, akx ! MRH added to index with ky, kx

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: dist
!    real, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last

    integer :: iglo, ik, it, is, ie, il, ig
    integer, save :: unit
    logical :: done = .false.
    real :: vpa1
    complex, dimension (2) :: gtmp
!    real, dimension (2) :: gtmp

     if (.not. done) then
!        if (proc0) call open_output_file (unit, ".distmp")
        if (proc0) call open_output_file (unit, trim(extension))
        do iglo=g_lo%llim_world, g_lo%ulim_world
           ik = ik_idx(g_lo, iglo)
           it = it_idx(g_lo, iglo)
           is = is_idx(g_lo, iglo)
           ie = ie_idx(g_lo, iglo)
           il = il_idx(g_lo, iglo)
           do ig = -ntgrid, ntgrid
              if (idx_local (g_lo, ik, it, il, ie, is)) then
                 if (proc0) then
                    gtmp = dist(ig,:,iglo)
                 else
                    call send (dist(ig,:,iglo), 0)
                 end if
              else if (proc0) then
                 call receive (gtmp, proc_id(g_lo, iglo))
              end if
              if (proc0) then
                 if (.not. forbid(ig,il)) then
                    vpa1 = sqrt(energy_maxw(ie)*max(0.0,1.0-al(il)*bmag(ig)))
                    write (unit,'(9e14.5)') aky(ik),akx(it),theta(ig), vpa1, energy_maxw(ie), &
                         real(gtmp(1)), aimag(gtmp(1)), real(gtmp(2)), aimag(gtmp(2))
                    ! write (unit,'(7e14.5)') theta(ig), vpa1, energy(ie), &
                         ! real(gtmp(1)), aimag(gtmp(1)), real(gtmp(2)), aimag(gtmp(2))
!                    write (unit,'(5e14.5)') theta(ig), vpa1, energy(ie), gtmp(1), gtmp(2)
                 end if
              end if
           end do
           if (proc0) then
              write (unit,*)
              write (unit,*)
           end if
        end do
        if (proc0) call close_output_file (unit)
        if (present(last)) done = .true.
     end if

   end subroutine write_mpdist

  !> Used for testing
  !! takes as input an array using le_lo and
  !! writes it to a .distmp output file
  subroutine write_mpdist_le (dist, extension, last)
    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: le_lo, ig_idx, ik_idx, it_idx, is_idx
    use gs2_layouts, only: idx_local, proc_id
    use theta_grid, only: bmag, theta

    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in) :: dist
!    real, dimension (:,:,le_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last
    
    integer :: ile, ixi, nxi, ik, it, is, ie, il, ig, isgn
    integer, save :: unit
    logical :: done = .false.
    real :: vpa1
    complex :: gtmp
!    real :: gtmp

    if (.not. done) then
       if (proc0) call open_output_file (unit, trim(extension))
       nxi = max(2*nlambda-1, 2*ng2)
       do ile=le_lo%llim_world, le_lo%ulim_world
          ig = ig_idx(le_lo, ile)
          ik = ik_idx(le_lo, ile)
          it = it_idx(le_lo, ile)
          is = is_idx(le_lo, ile)
          do ie = 1, negrid
             do ixi = 1, nxi
                il = ixi_to_il(ig,ixi)
                isgn = ixi_to_isgn(ig,ixi)
                if (idx_local (le_lo, ig, ik, it, is)) then
                   if (proc0) then
                      gtmp = dist(ixi,ie,ile)
                   else
                      call send (dist(ixi,ie,ile), 0)
                   end if
                else if (proc0) then
                   call receive (gtmp, proc_id(le_lo, ile))
                end if
                if (proc0) then
                   if (.not. forbid(ig,il)) then
                      vpa1 = sqrt(energy_maxw(ie)*max(0.0,1.0-al(il)*bmag(ig)))
                      write (unit,'(5e14.5)') theta(ig), sgn(isgn)*vpa1, energy_maxw(ie), &
                           real(gtmp), aimag(gtmp)
                      !                    write (unit,'(5e14.5)') theta(ig), vpa1, energy(ie), gtmp(1), gtmp(2)
                   end if
                end if
             end do
          end do
       end do
       if (proc0) call close_output_file (unit)
       if (present(last)) done = .true.
    end if
    
  end subroutine write_mpdist_le

  !> FIXME : Add documentation  
  subroutine finish_le_grids
    use redistribute, only: delete_redist
    use egrid, only: zeroes, x0, zeroes_maxw
    use species, only: nonmaxw_corr, f0_values, dlnf0drho

    implicit none
    real, dimension(1,1,1,1,1) :: tmpf
    real, dimension(1,1,1) :: tmpvflx

    if (allocated(zeroes)) deallocate (zeroes,x0,zeroes_maxw)
    if (allocated(energy)) deallocate (energy, energy_maxw,dele, al, wl, jend, forbid, xx, lgrnge, xloc, speed, speed_maxw)
    if (allocated(integration_work)) deallocate (integration_work)
    if (allocated(wlterr)) deallocate (wlterr)
    if (allocated(w)) deallocate (w,w_maxw)
    if (allocated(delal)) deallocate (delal)
    if (allocated(lpl)) deallocate (lpl, lpe)
    if (allocated(lpt)) deallocate (lpt)
    if (allocated(wlmod)) deallocate (wlmod)
    if (allocated(wmod)) deallocate (wmod)
    if (allocated(wtmod)) deallocate (wtmod)
    if (allocated(xi)) deallocate (xi)
    if (allocated(ixi_to_il)) deallocate (ixi_to_il)
    if (allocated(ixi_to_isgn)) deallocate (ixi_to_isgn)
    if (allocated(recvcnts_intspec)) deallocate (recvcnts_intspec)
    if (allocated(displs_intspec)) deallocate (displs_intspec)

    if (allocated(nonmaxw_corr)) deallocate(nonmaxw_corr)
    if (allocated(f0_values)) deallocate(f0_values)
    if (allocated(dlnf0drho)) deallocate(dlnf0drho)

    !Deallocate redistribute maps
    call delete_redist(lambda_map)
    call delete_redist(energy_map)
    call delete_redist(g2le)
    call delete_redist(g2gf)

    !Free internal saved arrays
    call get_flux_vs_theta_vs_vpa(tmpf,tmpvflx,.true.)

    accel_x = .false. ; accel_v = .false.
    test = .false. ; trapped_particles = .true.
    new_trap_int = .false. 

    intinit = .false. ; slintinit = .false. ; lintinit = .false. ; eintinit = .false.

    leinit = .false.
    gfinit = .false.
    lzinit = .false.
    einit = .false.
    initialized = .false.

    call le_grids_config%reset()
  end subroutine finish_le_grids

  !> Deallocate the weights used to provide error estimates
  subroutine finish_weights
    if (allocated(werr)) deallocate (werr)
    if (allocated(wlerr)) deallocate (wlerr)
    init_weights_init = .false.
  end subroutine finish_weights

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_le_grids_config(le_grids_config_in)
    use mp, only: mp_abort
    type(le_grids_config_type), intent(in), optional :: le_grids_config_in
    if (initialized) then
       call mp_abort("Trying to set le_grids_config when already initialized.", to_screen = .true.)
    end if
    if (present(le_grids_config_in)) then
       le_grids_config = le_grids_config_in
    end if
  end subroutine set_le_grids_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the le_grids_knobs namelist and populates the member variables
  subroutine read_le_grids_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(le_grids_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    integer :: negrid, nesub, nesuper, ngauss, nmax, nterp
    logical :: genquad, new_trap_int, radau_gauss_grid, test, trapped_particles
    real :: bouncefuzz, vcut, wgt_fac

    namelist /le_grids_knobs/ bouncefuzz, genquad, negrid, nesub, nesuper, new_trap_int, npassing, ngauss, nmax, nterp, radau_gauss_grid, test, trapped_particles, vcut, wgt_fac

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    bouncefuzz = self%bouncefuzz
    genquad = self%genquad
    negrid = self%negrid
    nesub = self%nesub
    nesuper = self%nesuper
    new_trap_int = self%new_trap_int
    npassing = self%npassing
    ngauss = self%ngauss
    nmax = self%nmax
    nterp = self%nterp
    radau_gauss_grid = self%radau_gauss_grid
    test = self%test
    trapped_particles = self%trapped_particles
    vcut = self%vcut
    wgt_fac = self%wgt_fac

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = le_grids_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = le_grids_knobs)
    endif

    ! Now copy from local variables into type members
    self%bouncefuzz = bouncefuzz
    self%genquad = genquad
    self%negrid = negrid
    self%nesub = nesub
    self%nesuper = nesuper
    self%new_trap_int = new_trap_int
    self%npassing = npassing
    self%ngauss = ngauss
    self%nmax = nmax
    self%nterp = nterp
    self%radau_gauss_grid = radau_gauss_grid
    self%test = test
    self%trapped_particles = trapped_particles
    self%vcut = vcut
    self%wgt_fac = wgt_fac

    self%exist = exist
  end subroutine read_le_grids_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_le_grids_config(self, unit)
    implicit none
    class(le_grids_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("bouncefuzz", self%bouncefuzz, unit_internal)
    call self%write_key_val("genquad", self%genquad, unit_internal)
    call self%write_key_val("negrid", self%negrid, unit_internal)
    call self%write_key_val("nesub", self%nesub, unit_internal)
    call self%write_key_val("nesuper", self%nesuper, unit_internal)
    call self%write_key_val("new_trap_int", self%new_trap_int, unit_internal)
    call self%write_key_val("npassing", self%npassing, unit_internal)
    call self%write_key_val("ngauss", self%ngauss, unit_internal)
    call self%write_key_val("nmax", self%nmax, unit_internal)
    call self%write_key_val("nterp", self%nterp, unit_internal)
    call self%write_key_val("radau_gauss_grid", self%radau_gauss_grid, unit_internal)
    call self%write_key_val("test", self%test, unit_internal)
    call self%write_key_val("trapped_particles", self%trapped_particles, unit_internal)
    call self%write_key_val("vcut", self%vcut, unit_internal)
    call self%write_key_val("wgt_fac", self%wgt_fac, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_le_grids_config

  !> Resets the config object to the initial empty state
  subroutine reset_le_grids_config(self)
    class(le_grids_config_type), intent(in out) :: self
    type(le_grids_config_type) :: empty
    select type (self)
    type is (le_grids_config_type)
       self = empty
    end select
  end subroutine reset_le_grids_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_le_grids_config(self)
    use mp, only: broadcast
    implicit none
    class(le_grids_config_type), intent(in out) :: self
    call broadcast(self%bouncefuzz)
    call broadcast(self%genquad)
    call broadcast(self%negrid)
    call broadcast(self%nesub)
    call broadcast(self%nesuper)
    call broadcast(self%new_trap_int)
    call broadcast(self%npassing)
    call broadcast(self%ngauss)
    call broadcast(self%nmax)
    call broadcast(self%nterp)
    call broadcast(self%radau_gauss_grid)
    call broadcast(self%test)
    call broadcast(self%trapped_particles)
    call broadcast(self%vcut)
    call broadcast(self%wgt_fac)

    call broadcast(self%exist)
  end subroutine broadcast_le_grids_config

  !> Gets the default name for this namelist
  function get_default_name_le_grids_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_le_grids_config
    get_default_name_le_grids_config = "le_grids_knobs"
  end function get_default_name_le_grids_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_le_grids_config()
    implicit none
    logical :: get_default_requires_index_le_grids_config
    get_default_requires_index_le_grids_config = .false.
  end function get_default_requires_index_le_grids_config
end module le_grids



