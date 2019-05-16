!> Analytic local equilibria via generalised forms of the Miller model
!!
!! Justin Ball's MIT Masters thesis and/or "GS2 analytic geometry
!! specification" on the wiki provide a good overview of how this
!! module works.
!!
!! There are four equilibrium models in this module:
!!
!! 1. generalised Miller
!! 2. global MHD
!! 3. generalised elongation
!! 4. Fourier
!!
!! These can be chosen via the `geoType` input parameter (`0-4`,
!! respectively)
!!
!! All of these models allow you to specify the amplitudes of two
!! shaping modes, along with their phases and radial derivatives. The
!! generalised elongation and Fourier models also allow you to specify
!! the exact mode numbers.
module leq
  implicit none

  private
  
  public :: leq_init, leqin, gradient, eqitem, bgradient, leq_finish
  public :: invR, Rpos, Zpos, diameter, btori, dbtori, qfun, pfun
  public :: dpfun, betafun, psi, rcenter, dpdrhofun

  public :: leq_testing_type

  !> Number of radial grid points
  !!
  !! Unconditionally always set to 3?
  integer :: nr
  !> Number of theta grid points
  integer :: nt

  real, allocatable, dimension (:)     :: eqpsi, fp, beta, pressure
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi

  ! The following quantities represent changes between adjacent grid
  ! points. The first two dimensions are minor radius and theta,
  ! respectively. The last dimension is the difference along either
  ! the (minor) radial or theta directions, respectively

  !> Major radius differences on theta, minor radius grids
  real, allocatable, dimension (:,:,:) :: drm
  !> Vertical differences on theta, minor radius grids
  real, allocatable, dimension (:,:,:) :: dzm
  !> Magnetic field magnitude differences on theta, minor radius grids
  real, allocatable, dimension (:,:,:) :: dbtm
  !> Minor radius differences on theta, minor radius grids, note: `dpm(:,:,2) == 0.0`
  real, allocatable, dimension (:,:,:) :: dpm
  !> Theta differences on theta, minor radius grids, note: `dtm(:,:,1) == 0.0`
  real, allocatable, dimension (:,:,:) :: dtm
  !> Minor radius gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dpcart
  !> Theta gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dtcart
  !> Magnetic field magnitude gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dbtcart
  !> Minor radius gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dpbish
  !> Theta gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dtbish
  !> Magnetic field magnitude gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
  real, allocatable, dimension (:,:,:) :: dbtbish

  !> FIXME : Add documentation
  real :: beta_0

  !> Used to translate definition of the poloidal angle such that the location
  !> of the maximium magnetic field is at pi (see section 3.2.3 of Ball MIT
  !> Masters thesis or section 2.4 of "GS2 analytic geometry specification")
  real :: thetaShift

  !> Private helper type to wrap up some module-level variables
  type :: flux_surface
     !> Flux surface major radius, \(R_{0N}\), input variable: `Rmaj`
     real :: Rmaj
     !> Location of reference magnetic field, \(R_{geoN}\), input variable: `R_geo`
     real :: Rgeo
     !> Minor radius, \(r_{\psi N}\), input variable: `rhoc`
     real :: r
     !> Step size for radial derivatives, \(\Delta r_{\psi N}\), input variable: `delrho`
     real :: dr
     !> Minor radius of shaping surface, \(a_{\psi N}\), input variable: `aSurf`
     real :: aSurf
     !> Horizontal Shafranov shift, \(dR_{0N}/dr_{\psi N}\), input variable: `shift`, `Raxis`
     real :: sHorz
     !> Vertical Shafranov shift, \(dZ_{0N}/dr_{\psi N}\), input variable: `shiftVert`, `Zaxis`
     real :: sVert
     !> First shaping mode strength, \(\Delta_m\), input variable: `akappa`, `delta2`, `deltam`
     real :: delm
     !> Second shaping mode strength, \(\Delta_n\), input variable: `tri`, `delta3`, `deltan`
     real :: deln
     !> First shaping mode derivative, \(d\Delta_m/dr_{\psi N}\), input variable: `akappri`, `deltampri`
     real :: delmp
     !> Second shaping mode derivative, \(d\Delta_n/dr_{\psi N}\), input variable: `tripri`, `deltanpri`
     real :: delnp
     !> First shaping mode tilt angle, \(\theta_m\), input variable: `thetak`, `theta2`, `thetam`
     real :: thm
     !> Second shaping mode tilt angle, \(\theta_n\), input variable: `thetad`, `theta2`, `thetan`
     real :: thn
     !> Safety factor, \(q\), input variable: `qinp`
     real :: q
     !> Magnetic shear, \(\hat{s}\), input variable: `s_hat_input`
     real :: shat
     !> Number of points in theta, input variable: `ntheta(nperiod - 1/2)`
     integer :: nt
     !> Geometry specification type, input variable: `geoType`
     integer :: geoType
     !> The first poloidal mode number, input variable: `mMode`
     integer :: mMode
     !> The second poloidal mode number, input variable: `nMode`
     integer :: nMode
  end type flux_surface

  type (flux_surface) :: surf

  !> Expose private routines for testing only!
  type :: leq_testing_type
  contains
    procedure, nopass :: mod2pi
    procedure, nopass :: derm
    procedure, nopass :: set_ntnr
    procedure, nopass :: set_drm
    procedure, nopass :: set_dzm
    procedure, nopass :: set_dpcart
    procedure, nopass :: set_cart_arrays
    procedure, nopass :: set_bish_arrays
    procedure, nopass :: eqdcart
    procedure, nopass :: eqdbish
    procedure, nopass :: alloc_arrays, dealloc_arrays
    procedure, nopass :: set_surf
  end type leq_testing_type

contains
  
  !> Set module-level `surf`
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_surf(Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, &
       delm, deln, delmp, delnp, thm, thn, q, shat, nt, geoType, mMode, nMode)
     real, intent(in) :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert,&
          delm, deln, delmp, delnp, thm, thn, q, shat
     integer, intent(in) :: nt, geoType, mMode, nMode

     surf = flux_surface(Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, &
       delm, deln, delmp, delnp, thm, thn, q, shat, nt, geoType, mMode, nMode)
  end subroutine set_surf

  !> Set module-level `nt`, `nr`
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_ntnr(ntheta, nradius)
    integer, intent(in) :: ntheta, nradius
    nt = ntheta
    nr = nradius
  end subroutine set_ntnr

  !> Set module-level `drm`
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_drm(drm_in)
    real, dimension(nr,-nt:nt,2), intent(in) :: drm_in
    drm = drm_in
  end subroutine set_drm

  !> Set module-level `dzm`
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_dzm(dzm_in)
    real, dimension(nr,-nt:nt,2), intent(in) :: dzm_in
    dzm = dzm_in
  end subroutine set_dzm

  !> Set module-level dpcart
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_dpcart(dpcart_in)
    real, dimension(nr,-nt:nt,2), intent(in) :: dpcart_in
    dpcart = dpcart_in
  end subroutine set_dpcart

  !> Set module-level dpcart, dbtcart and dtcart
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_cart_arrays(dpcart_in, dbtcart_in, dtcart_in)
    real, dimension(nr,-nt:nt,2), intent(in) :: dpcart_in
    real, dimension(nr,-nt:nt,2), intent(in) :: dbtcart_in
    real, dimension(nr,-nt:nt,2), intent(in) :: dtcart_in
    dpcart = dpcart_in
    dbtcart = dbtcart_in
    dtcart = dtcart_in
  end subroutine set_cart_arrays

  !> Set module-level dpbish, dbtbish and dtbish
  !! PURELY FOR TESTING
  !! FIXME: Remove when we have more control over module-level variables for testing
  subroutine set_bish_arrays(dpbish_in, dbtbish_in, dtbish_in)
    real, dimension(nr,-nt:nt,2), intent(in) :: dpbish_in
    real, dimension(nr,-nt:nt,2), intent(in) :: dbtbish_in
    real, dimension(nr,-nt:nt,2), intent(in) :: dtbish_in
    dpbish = dpbish_in
    dbtbish = dbtbish_in
    dtbish = dtbish_in
  end subroutine set_bish_arrays

  !> Set module-level parameters and initialise the leq module
  !!
  !! Note: in contrast to the other geometry modules, this also calls
  !! `leq_init`, whereas others have separate `Xeqin` and `Xeq_init`
  !! procedures
  subroutine leqin(geoType, Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, mMode, nMode, &
       delm, deln, delmp, delnp, thm, thn, q, shat, nt_used, thShift)
    implicit none
    real, intent(in) :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, delm, deln, &
         delmp, delnp, thm, thn, q, shat
    real, intent(in) :: thShift
    integer, intent(in) :: nt_used, geoType, mMode, nMode

    surf%geoType = geoType
    surf%Rmaj = Rmaj
    surf%Rgeo = Rgeo
    surf%r = r
    surf%dr = dr 
    surf%aSurf = aSurf
    surf%sHorz = sHorz
    surf%sVert = sVert
    surf%mMode = mMode
    surf%nMode = nMode
    surf%delm = delm
    surf%deln = deln
    surf%delmp = delmp
    surf%delnp = delnp
    surf%thm = thm
    surf%thn = thn
    surf%q = q
    surf%shat = shat
 
    beta_0 = 1.

    thetaShift = thShift;

    nr = 3
    nt = nt_used
    if(.not.allocated(beta)) call alloc_arrays(3, nt)
    surf%nt = nt
    call leq_init
  end subroutine leqin

  !> Allocate internal module-level arrays
  !!
  !! Note that this doesn't use the module-level `nr`, `nt`
  subroutine alloc_arrays(nr, nt)
    implicit none
    !> Number of radial points
    integer, intent(in) :: nr
    !> Number of points in theta
    integer, intent(in) :: nt

    allocate(eqpsi(nr), fp(nr), beta(nr), pressure(nr))
    allocate(R_psi(nr, -nt:nt), Z_psi(nr, -nt:nt))
    allocate(drm(nr, -nt:nt, 2), dzm(nr, -nt:nt, 2), dbtm(nr, -nt:nt, 2), &
         dpm(nr, -nt:nt, 2), dtm(nr, -nt:nt, 2))
    allocate(dpcart(nr, -nt:nt, 2), dtcart(nr, -nt:nt, 2), dbtcart(nr, -nt:nt,2))
    allocate(dpbish(nr, -nt:nt, 2), dtbish(nr, -nt:nt, 2), dbtbish(nr, -nt:nt,2))

  end subroutine alloc_arrays

  !> Deallocate internal module-level arrays
  subroutine dealloc_arrays
    implicit none
    if(allocated(eqpsi)) deallocate(eqpsi, fp, beta, pressure)
    if(allocated(R_psi)) deallocate(R_psi, Z_psi)
    if(allocated(drm)) deallocate(drm,dzm,dbtm,dpm,dtm)
    if(allocated(dpcart)) deallocate(dpcart,dtcart,dbtcart)
    if(allocated(dpbish)) deallocate(dpbish,dtbish,dbtbish)
  end subroutine dealloc_arrays

  !> Finalise leq module
  subroutine leq_finish
    implicit none
    call dealloc_arrays
  end subroutine leq_finish

  !> Initialise leq module
  subroutine leq_init
    use constants, only: pi
    implicit none
    real, dimension(nr, -nt:nt) :: eqpsi1, eqth, eqbtor
    real :: dr(3)
    real :: t, r
    integer :: i, j
   
    dr(1) = -surf%dr
    dr(2) = 0.
    dr(3) = surf%dr
   
    ! All calculation must be extended over the full [-pi,pi] range in order to
    ! support up-down asymmetric geometries (see section 3.2.2 of Ball MIT
    ! Masters thesis or section 2.2 of "GS2 analytic geometry specification")
    !
    ! NOTE: actual theta grid used here includes one extra point outside [-pi, pi]!
    do j=-nt,nt
       do i=1,nr
          r = surf%r + dr(i)
          t = (j)*pi/real(nt-1)
          R_psi(i,j) = Rpos(r, t) 
          Z_psi(i,j) = Zpos(r, t)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
          eqbtor(i,j) = surf%Rgeo/R_psi(i,j)
       enddo
    enddo

    do i=1,nr
       pressure(i) = -dr(i)
    enddo

    eqpsi(:) = eqpsi1(:,1)

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')

    ! below is actually grad(rho) instead of grad(psi),
    ! and 'cartesian' refers to (R,Z) coordinates -- MAB
    ! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
    ! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

    ! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
    ! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

    ! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
    ! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

  end subroutine leq_init

  !> Index-space derivatives of `f` in each direction. Gradient in
  !> first direction is in `dfm(:,:,1)`, second is in `dfm(:,:,2)`.
  !> `char` controls the boundary conditions in the second direction.
  !!
  !! Acceptable values for `char`:
  !!   - `E` means continuous and even at theta = +/- pi -- FIXME: confirm
  !!   - `O` means continuous and odd at theta = +/- pi -- FIXME: confirm
  !!   - `T` means a jump of 2 pi at theta = +/- pi
  subroutine derm(f, dfm, char)
    use constants, only: pi
    implicit none
    real, dimension(:,:), intent(in) :: f(:,-nt:)
    real, dimension(:,:,:), intent(out) :: dfm(:,-nt:,:)
    character(1), intent(in) :: char
    integer :: i, j
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    select case (char)
    case ('E')
       j=-nt
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=-nt
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=-nt
       dfm(:,j,2) = f(:,j+1) + pi

       j=nt
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo

    do j=-nt+1,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo

  end subroutine derm

  !> Calculate the derivative of \(r_\psi\) w.r.t. R and Z
  !!
  !! \(r_\psi\) is the "minor radius flux function"
  subroutine gradient(theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array is -ntm
    integer, intent(in) :: ntm
    !> If char = 'R', return the gradient of pressure instead.
    !> Note that this is useless here!
    character(1), intent(in) :: char
    !> Theta grid
    real, dimension(-ntm:), intent(in) :: theta
    !> grad(:,1): derivative w.r.t. R as function of theta
    !> grad(:,2): derivative w.r.t  Z as function of theta
    real, dimension(-ntm:,:), intent(out) :: grad
    !> Value of psi/(a^2 B_a) on the flux surface where we want the gradient.
    !> Note: never used!
    real, intent(in) :: rp
    real :: tmp(2), aa(1), daa(1), rpt(1)
    real, dimension(nr,-nt:nt,2) :: dcart
    integer :: i
    
    select case(char)
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R' -- never reached: other calculation done instead
    case('T')
       dcart = dtcart  ! Never reached: bgradient always used
    end select
    
    do i=-nth_used,-1
       call eqitem(theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(theta(i), dcart(:,:,2), tmp(2), 'Z')
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       call eqitem(theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(theta(i), dcart(:,:,2), tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1)*0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1)*0.5*beta_0
       enddo
    endif
  end subroutine gradient

  !> Calculate the derivative of \(r_\psi\) in the Bishop coordinate system.
  !> If `char=='T'`, calculate the gradient of \(\theta\) instead
  subroutine bgradient(theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array is -ntm
    integer, intent(in) :: ntm
    !> If char = 'R', return the gradient of pressure instead.
    !> Note that this is useless here!
    character(1), intent(in) :: char
    !> Theta grid
    real, dimension(-ntm:), intent(in) :: theta
    !> grad(:,1): derivative w.r.t. R as function of theta
    !> grad(:,2): derivative w.r.t  Z as function of theta
    real, dimension(-ntm:,:), intent(out) :: grad
    !> Value of psi/(a^2 B_a) on the flux surface where we want the gradient.
    !> Note: never used!
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(nr,-nt:nt,2) ::  dbish
    integer :: i

    select case(char)
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R' -- never reached: other calculation done instead
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(theta(i), dbish(:,:,1), grad(i,1), 'R')
       call eqitem(theta(i), dbish(:,:,2), grad(i,2), 'Z')
    enddo

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif
  end subroutine bgradient

  !> Calculates `fstar` which is `f` interpolated at the location `theta_in`
  subroutine eqitem(theta_in, f, fstar, char)
    use constants, only: pi
    implicit none
    !> Desired theta location
    real, intent(in) :: theta_in
    !> Array to interpolate
    real, dimension(:,-nt:), intent(in) :: f
    !> Interpolation result
    real, intent(out) :: fstar
    !> FIXME: Unused
    character(1), intent(in) :: char
    !> Internal loop variable
    integer :: j
    !> Found theta index
    integer :: jstar
    !> r index, always centre of grid
    integer, parameter :: istar = 2
    !> theta point in [-pi, pi]
    real :: tp
    !> linear interpolation weights
    real :: st, dt
    !> Internal theta grid
    real, dimension(-nt:nt) :: mtheta

    tp = mod2pi(theta_in)

    ! get thet on theta mesh

    ! NOTE: actual theta grid used here includes one extra point outside [-pi, pi]!
    mtheta = (/ ( real(j)*pi/real(nt-1), j=-nt,nt) /)
  
    ! note that theta(1)=0 for local_eq theta 

    jstar=-nt-1
    do j=-nt,nt
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
          exit
       endif
       if(jstar /= -nt-1) write(*,*) 'exit error j'
    enddo
      
    ! treat theta = pi separately
  
    if(jstar == -nt-1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

    ! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt
    fstar=fstar &
         /(mtheta(jstar+1)-mtheta(jstar))
  end subroutine eqitem

  !> Convert `dfm` from index-space derivative to derivative on
  !> cylindrical coordinates
  !!
  !! - Gradient in the first coordinate is in `(:,:,1)`
  !! - Gradient in the second coordinate is in `(:,:,2)`
  subroutine eqdcart(dfm, dfcart)
    !> Index-space derivative
    real, dimension(nr,-nt:nt,2), intent(in)  :: dfm
    !> Cylindrical coordinates derivative
    real, dimension(nr,-nt:nt,2), intent(out) :: dfcart
    real, dimension(nr,-nt:nt) :: denom

    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

    dfcart = 0.
    
    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)
    
    ! FIXME: document why this is only on the range [2,nr], or
    ! simplify
    do j=-nt,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo
  end subroutine eqdcart

  !> Convert `dcart` from cylindrical coordinates derivative to
  !> derivative on Bishop coordinates
  !!
  !! - Gradient in the first coordinate is in `(:,:,1)`
  !! - Gradient in the second coordinate is in `(:,:,2)`
  subroutine eqdbish(dcart, dbish)
    !> Cylindrical coordinates derivative
    real, dimension(nr,-nt:nt,2), intent(in) :: dcart
    !> Bishop coordinates derivative
    real, dimension(nr,-nt:nt,2), intent(out) :: dbish
    real, dimension(nr,-nt:nt) :: denom

    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    ! FIXME: document why this is only on the range [2,nr], or
    ! simplify
    do j=-nt,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
       enddo
    enddo
  end subroutine eqdbish

  !> Return the inverse of the major radius at \((r, \theta)\)
  function invR (r, theta)
    implicit none
    real, intent (in) :: r, theta
    real :: invR
    
    invR=1./Rpos(r, theta)
  end function invR

  !> Return the major radius of the centre of the flux surface
  function rcenter()
    implicit none
    real :: rcenter
    
    ! If using the global geometry, Rmaj is the center of the shaped flux surface
    ! (i.e. the surface at rhoc=aSurf), so here we calculate the major radial
    ! location of the center of the flux surface of interest ! JRB
    if (surf%geoType==1) then
       rcenter = surf%Rmaj+(1-(surf%r/surf%aSurf)**2)*surf%sHorz
    else
       rcenter = surf%Rmaj
    end if
  end function rcenter

  !> Return the major radius at \((r, \theta)\)
  function Rpos (r, theta)
    implicit none
    real, intent (in) :: r, theta
    real :: Rpos, dummy

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (surf%geoType)
       case (0)
          call geo_Miller(r, theta, Rpos, dummy)
       case (1)
          call geo_global(r, theta, Rpos, dummy)
       case (2)
          call geo_generalizedEllipticity(r, theta, Rpos, dummy)
       case (3)
          call geo_FourierSeries(r, theta, Rpos, dummy)
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select
  end function Rpos

  !> Return the axial coordinate at \((r, \theta)\)
  function Zpos (r, theta)
    implicit none
    real, intent (in) :: r, theta
    real :: dummy, Zpos

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (surf%geoType)
       case (0)
          call geo_Miller(r, theta, dummy, Zpos)
       case (1)
          call geo_global(r, theta, dummy, Zpos)
       case (2)
          call geo_generalizedEllipticity(r, theta, dummy, Zpos)
       case (3)
          call geo_FourierSeries(r, theta, dummy, Zpos)
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select
  end function Zpos

  !> Calculate the cylindrical \((R, Z)$$ at $$(r, \theta)\) for the
  !> generalised Miller model
  !!
  !! See section 3.2.1 of Ball MIT Masters thesis or section 2.1.1 of
  !! "GS2 analytic geometry specification"
  subroutine geo_Miller(r, theta, Rpos, Zpos)
    implicit none
    real, intent(in) :: r, theta
    real, intent(out) :: Rpos, Zpos
    real :: dr, thAdj
    real :: Rcirc, Rcircp, Relong, Relongp, RelongTilt, RelongTiltp, Rtri, Rtrip, RtriTilt, RtriTiltp, Rfinal, Rfinalp
    real :: Zcirc, Zcircp, Zelong, Zelongp, ZelongTilt, ZelongTiltp, Ztri, Ztrip, ZtriTilt, ZtriTiltp, Zfinal, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    Rcirc = surf%r*cos(surf%thn - surf%thm - thAdj)
    Zcirc = -surf%r*sin(surf%thn - surf%thm - thAdj)
    Rcircp = cos(surf%thn - surf%thm - thAdj)
    Zcircp = -sin(surf%thn - surf%thm - thAdj)

    Relong = Rcirc
    Zelong = -surf%r*(surf%delm - 1)*sin(surf%thn - surf%thm - thAdj) + Zcirc
    Relongp = Rcircp
    Zelongp = -((surf%delm - 1) + surf%r*surf%delmp)*sin(surf%thn - surf%thm - thAdj) + Zcircp

    RelongTilt = cos(surf%thn - surf%thm)*Relong - sin(surf%thn - surf%thm)*Zelong
    ZelongTilt = sin(surf%thn - surf%thm)*Relong + cos(surf%thn - surf%thm)*Zelong
    RelongTiltp = cos(surf%thn - surf%thm)*Relongp - sin(surf%thn - surf%thm)*Zelongp
    ZelongTiltp = sin(surf%thn - surf%thm)*Relongp + cos(surf%thn - surf%thm)*Zelongp

    Rtri = -surf%r*cos(thAdj) + surf%r*cos(sin(thAdj)*surf%deln + thAdj) + RelongTilt
    Ztri = ZelongTilt
    Rtrip = -cos(thAdj) + cos(sin(thAdj)*surf%deln + thAdj) &
            - surf%r*sin(sin(thAdj)*surf%deln + thAdj)*sin(thAdj)*surf%delnp + RelongTiltp
    Ztrip = ZelongTiltp

    RtriTilt = cos(surf%thn)*Rtri + sin(surf%thn)*Ztri
    ZtriTilt = -sin(surf%thn)*Rtri + cos(surf%thn)*Ztri
    RtriTiltp = cos(surf%thn)*Rtrip + sin(surf%thn)*Ztrip
    ZtriTiltp = -sin(surf%thn)*Rtrip + cos(surf%thn)*Ztrip

    Rfinal = RtriTilt + surf%Rmaj
    Zfinal = ZtriTilt
    Rfinalp = RtriTiltp + surf%sHorz
    Zfinalp = ZtriTiltp + surf%sVert

    Rpos = Rfinal + dr*Rfinalp
    Zpos = Zfinal + dr*Zfinalp

  end subroutine geo_Miller

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> global MHD model
  !!
  !! See section 2.1.2 of "GS2 analytic geometry specification"
  subroutine geo_global(r, theta, Rpos, Zpos)
    implicit none
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: rho0, thAdj, rCyl, Rc, Zc
    real :: Cm, Cn, psiN, xAng, yAng, rAng   

    rho0 = r/surf%aSurf

    thAdj = theta + thetaShift

    Cm=(surf%delm**2-1)/(surf%delm**2+1)
    Cn=(surf%deln**2-1)/(surf%deln**3+1)

    psiN=rho0**2+Cm*rho0**2+Cn*rho0**3

    if (abs(Cn*cos(3*(thAdj+surf%thn)))<10e-6) then
      rCyl=surf%aSurf*sqrt(psiN/(1+Cm*cos(2*(thAdj+surf%thm))))
    else
      xAng=2*(1+Cm*cos(2*(thAdj+surf%thm)))**3-27*Cn**2*psiN*cos(3*(thAdj+surf%thn))**2
      yAng=3*Cn*sqrt(3*psiN)*cos(3*(thAdj+surf%thn))*sqrt(xAng+2*(1+Cm*cos(2*(thAdj+surf%thm)))**3)

      rAng=(1/3.0)*atan2(yAng,xAng)

      rCyl=surf%aSurf*(1+Cm*cos(2*(thAdj+surf%thm)))/(3*Cn*cos(3*(thAdj+surf%thn))) &
           *(cos(rAng)+sqrt(3.0)*sin(rAng)-1)
    end if

    Rc=surf%sHorz-(surf%sHorz-surf%Rmaj)*rho0**2
    Zc=surf%sVert-(surf%sVert)*rho0**2

    Rpos=Rc+rCyl*cos(thAdj)
    Zpos=Zc+rCyl*sin(thAdj)

  end subroutine geo_global

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> generalised elongation model
  !!
  !! See section 5.1.2 of Ball Oxford PhD thesis or section 2.1.3 of
  !! "GS2 analytic geometry specification"
  subroutine geo_generalizedEllipticity(r, theta, Rpos, Zpos)
    implicit none
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    rCyl=surf%r * (1 &
        + (-1+surf%delm/sqrt(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)) &
        + (-1+surf%deln/sqrt(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2.0)**2)))
    rCylp=1 &
     +(-1+(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)**(-0.5)*(surf%delm+surf%r*surf%delmp &
      *(1-(surf%delm*cos(surf%mMode*(thAdj+surf%thm)/2.0))**2/(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)))) &
     +(-1+(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2)**2)**(-0.5)*(surf%deln+surf%r*surf%delnp &
      *(1-(surf%deln*cos(surf%nMode*(thAdj+surf%thn)/2.0))**2/(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2.0)**2))))

    Rfinal=rCyl*cos(thAdj)+surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_generalizedEllipticity

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> Fourier model
  !!
  !! See section 5.1.1 of Ball Oxford PhD thesis or section 2.1.4 of
  !! "GS2 analytic geometry specification"
  subroutine geo_FourierSeries(r, theta, Rpos, Zpos)
    implicit none
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    rCyl=surf%r * (1+((surf%delm-1)/2) * (1-cos(surf%mMode*(thAdj+surf%thm))) + ((surf%deln-1)/2) * (1-cos(surf%nMode*(thAdj+surf%thn))))
    rCylp=1 + (((surf%delm-1)/2) + (surf%r*surf%delmp/2)) * (1-cos(surf%mMode*(thAdj+surf%thm))) + (((surf%deln-1)/2) + (surf%r*surf%delnp/2)) * (1-cos(surf%nMode*(thAdj+surf%thn)))

    Rfinal=rCyl*cos(thAdj)+surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_FourierSeries

  !> FIXME : Add documentation  
  function psi (r)
    implicit none
    real, intent (in) :: r
    real :: psi

    psi = r - surf%r
  end function psi

  !> Put `theta` into the range `[-pi, pi]`
  !!
  !! This version also sets `theta == +(-)pi +(-) theta_tol` to be
  !! exactly `+(-)pi`, where `theta_tol == 1.e-6`
  !!
  !! General purpose, could be moved into utils? Could also drop the
  !! tolerance?
  function mod2pi (theta)
    use constants, only: pi
    real, intent(in) :: theta
    real :: mod2pi
    real, parameter :: theta_tol = 1.e-6
    if(theta <= pi .and. theta >= -pi) then
      mod2pi = theta
      return
    else if(theta - theta_tol <= pi .and. theta >= -pi) then
      mod2pi = pi
      return
    else if(theta <= pi .and. theta + theta_tol >= -pi) then
      mod2pi = -pi
      return
    endif
    mod2pi = modulo(theta + pi, 2.*pi) - pi
  end function mod2pi

  !> FIXME : Add documentation  
  function diameter (rp)
    implicit none
    real, intent(in) :: rp
    real :: diameter

    diameter = 2.*rp
  end function diameter

  !> FIXME : Add documentation  
  function dbtori ()
    implicit none
    real :: dbtori
    dbtori = 1.
  end function dbtori

  !> FIXME : Add documentation  
  function btori ()
    implicit none
    real :: btori
    btori = surf%Rgeo
  end function btori

  !> FIXME : Add documentation  
  function qfun ()
    implicit none
    real :: qfun
    qfun = surf%q
  end function qfun

  !> FIXME : Add documentation
  function pfun ()
    implicit none
    real :: pfun
    pfun = 0.5*beta_0
  end function pfun

  !> FIXME : Add documentation  
  function dpfun ()
    implicit none
    real :: dpfun    
    dpfun = -1.
  end function dpfun

  !> FIXME : Add documentation  
  function dpdrhofun()
    implicit none
    real :: dpdrhofun
    dpdrhofun = 0 ! surf%pp
  end function dpdrhofun

  !> FIXME : Add documentation  
  function betafun ()
    implicit none
    real :: betafun
    betafun = beta_0
  end function betafun
end module leq
