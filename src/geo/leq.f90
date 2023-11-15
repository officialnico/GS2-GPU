!> Analytic local equilibria via generalised forms of the Miller model
!>
!> Justin Ball's MIT Masters thesis and/or ["GS2 analytic geometry
!> specification"](https://bitbucket.org/gyrokinetics/wikifiles/raw/master/JRB/GS2GeoDoc.pdf)
!> provide a good overview of how this module works.
!>
!> There are four equilibrium models in this module:
!>
!> 1. generalised Miller
!> 2. global MHD
!> 3. generalised elongation
!> 4. Fourier
!>
!> These can be chosen via the `geoType` input parameter (`0-4`,
!> respectively)
!>
!> All of these models allow you to specify the amplitudes of two
!> shaping modes, along with their phases and radial derivatives. The
!> generalised elongation and Fourier models also allow you to specify
!> the exact mode numbers.
module leq
  use geo_utils, only: abstract_geo_type, flux_surface_type, EQFILE_LENGTH
  implicit none

  private

  public :: leq_type, dpdrhofun

  type, extends(abstract_geo_type) :: leq_type
     private
     !> Used to translate definition of the poloidal angle such that the location
     !> of the maximium magnetic field is at pi (see section 3.2.3 of Ball MIT
     !> Masters thesis or section 2.4 of "GS2 analytic geometry specification")
     real :: thetaShift
     type(flux_surface_type) :: surf
     ! Regular variables
     real, allocatable, dimension(:) :: eqpsi, pressure
   contains
     procedure :: initialise => leqin
     procedure :: calculate_gradients => leq_init
     procedure :: finish_setup => leq_finish_setup
     procedure :: finalise => leq_finish
     procedure :: rpos => rpos
     procedure :: zpos => zpos
     procedure :: invr => invr
     procedure :: rcenter => rcenter
     procedure :: diameter => diameter
     procedure :: btori => btori
     procedure :: dbtori => dbtori
     procedure :: qfun => qfun
     procedure :: pfun => pfun
     procedure :: dpfun => dpfun
     procedure :: betafun => betafun
     procedure :: psi => psi
     procedure :: rfun => rfun
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     procedure :: eqdcart => eqdcart
     procedure :: eqdbish => eqdbish
     procedure :: alloc_arrays => alloc_arrays
     procedure :: dealloc_arrays => dealloc_arrays
     procedure :: geo_Miller => geo_Miller
     procedure :: geo_global => geo_global
     procedure :: geo_generalizedEllipticity => geo_generalizedEllipticity
     procedure :: geo_FourierSeries => geo_FourierSeries
  end type leq_type

contains

  !> Set module-level parameters and initialise the leq module
  subroutine leqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use mp, only: mp_abort
    use geo_utils, only: geo_input_type
    implicit none
    class(leq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(out) :: nthg
    if (inputs%surf%geoType > 4 .or. inputs%surf%geoType < 0) &
       call mp_abort("ERROR: invalid analytic geometry specification", .true.)

    if (.not. self%initialised) then
       self%surf = inputs%surf
       self%beta_0 = 1.
       self%thetaShift = inputs%thShift
       self%nr = 3
       self%nt = self%surf%nt
       call self%alloc_arrays()
    end if

    !Set outputs
    psi_0_out = 0.0
    psi_a_out = 1.0
    B_T0 = -1.0
    avgrmid = 1.0
    rmaj = self%surf%Rmaj
    nthg = self%nt * 2 !Output is equivalent to ntheta_input
    self%initialised = .true.
  end subroutine leqin

  !> Allocate internal module-level arrays
  !>
  !> Note that this doesn't use the module-level `nr`, `nt`
  pure subroutine alloc_arrays(self)
    implicit none
    class(leq_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%eqpsi(nr), self%pressure(nr))
    allocate(self%R_psi(nr, -nt:nt), self%Z_psi(nr, -nt:nt))
    allocate(self%drm(nr, -nt:nt, 2), self%dzm(nr, -nt:nt, 2), &
         self%dpm(nr, -nt:nt, 2), self%dtm(nr, -nt:nt, 2))
    allocate(self%dpcart(nr, -nt:nt, 2), self%dtcart(nr, -nt:nt, 2))
    allocate(self%dpbish(nr, -nt:nt, 2), self%dtbish(nr, -nt:nt, 2))
  end subroutine alloc_arrays

  !> Deallocate internal module-level arrays
  pure subroutine dealloc_arrays(self)
    implicit none
    class(leq_type), intent(in out) :: self
    if (allocated(self%eqpsi)) deallocate(self%eqpsi, self%pressure)
    if (allocated(self%R_psi)) deallocate(self%R_psi, self%Z_psi)
    if (allocated(self%drm)) deallocate(self%drm, self%dzm, self%dpm, self%dtm)
    if (allocated(self%dpcart)) deallocate(self%dpcart, self%dtcart)
    if (allocated(self%dpbish)) deallocate(self%dpbish, self%dtbish)
  end subroutine dealloc_arrays

  !> Finalise leq module
  pure subroutine leq_finish(self)
    implicit none
    class(leq_type), intent(in out) :: self
    self%initialised = .false.
    call self%dealloc_arrays
  end subroutine leq_finish

  !> Initialise leq module
  subroutine leq_init(self)
    use constants, only: pi
    implicit none
    class(leq_type), intent(in out) :: self
    real, dimension(self%nr, -self%nt:self%nt) :: eqpsi1, eqth
    real :: t, r, dr(3)
    integer :: i, j

    dr(1) = -self%surf%dr
    dr(2) = 0.
    dr(3) = self%surf%dr

    ! All calculation must be extended over the full [-pi,pi] range in order to
    ! support up-down asymmetric geometries (see section 3.2.2 of Ball MIT
    ! Masters thesis or section 2.2 of "GS2 analytic geometry specification")
    !
    do j = -self%nt, self%nt
       do i = 1, self%nr
          r = self%surf%r + dr(i)
          t = (j) * pi / real(self%nt)
          self%R_psi(i, j) = self%Rpos(r, t)
          self%Z_psi(i, j) = self%Zpos(r, t)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
       end do
    end do

    do i = 1, self%nr
       self%pressure(i) = -dr(i)
    end do

    self%eqpsi = eqpsi1(:,1)

    call self%derm(eqth,   self%dtm,  'T')
    ! char arguments E and O are redundant, but are retained so that `derm` interface
    ! is the same in all geometry modules.
    call self%derm(self%R_psi,  self%drm,  'E')
    call self%derm(self%Z_psi,  self%dzm,  'O')
    call self%derm(eqpsi1, self%dpm,  'E')

    ! below is actually grad(rho) instead of grad(psi),
    ! and 'cartesian' refers to (R,Z) coordinates -- MAB
    ! grad(psi) in cartesian form
    call self%eqdcart(self%dpm, self%dpcart)
    ! grad(psi) in Bishop form
    call self%eqdbish(self%dpcart, self%dpbish)

    ! grad(theta) in cartesian form
    call self%eqdcart(self%dtm, self%dtcart)
    ! grad(theta) in Bishop form
    call self%eqdbish(self%dtcart, self%dtbish)
  end subroutine leq_init

  !> Finish the leq setup
  subroutine leq_finish_setup(self)
    class(leq_type), intent(in out) :: self
    call self%calculate_gradients
  end subroutine leq_finish_setup

  !> `derm` computes the index-space derivatives of `f` in radius and theta directions,
  !> assuming periodicity in theta, `f(:,-nt) = f(:,nt)`.
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(leq_type), intent(in) :: self
    !> 2D array to be differentiated, indices radius, theta
    !> @note Theta grid end points assumed to be \[\pm \pi\]
    real, dimension(:,:), intent(in) :: f(:,-self%nt:)
    !> Derivatives of f:
    !>   dfm(nr,-nt:nt,1) = df/dr
    !>   dfm(nr,-nt:nt,2) = df/dtheta
    real, dimension(:,:,:), intent(out) :: dfm(:,-self%nt:,:)
    !> Specifies symmetry constraint at end of grid
    !>   "O" is odd in theta at \(\pm \pi\) (deprecated)
    !>   "E" is even in theta at \(\pm \pi\) (deprecated)
    !>   "T" denotes theta which trivially increases in theta at \(\pm \pi\)
    character(1), intent(in) :: char
    integer :: i, j, nr, nt
    nr = self%nr
    nt = self%nt
    ! Use second-order one-sided finite-differences for extreme radial points.
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    ! Compute df/dr
    do i = 2, nr-1
       dfm(i,:,1) = 0.5*(f(i+1,:)-f(i-1,:))
    end do

    ! Compute df/dtheta (except at endpoints -pi,pi)
    do j = -nt+1, nt-1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j-1))
    end do

    ! Compute df/dtheta at -pi,+pi assuming
    !   (i) periodicity in theta,
    !  (ii) endpoints at -pi,+pi
    ! (iii) equispaced grid in theta, and
    !  (iv) -nt and +nt correspond to -pi, +pi and are at the same place
    dfm(:,-nt,2) = 0.5*(f(:,-nt+1)-f(:,nt-1))
    dfm(:,nt,2) = dfm(:,-nt,2)

    ! `T` flag to denotes trivial derivative when f=theta,
    ! where care is needed at end points.
    if (char .eq. 'T') then
       dfm(:,-nt,2) = f(:,-nt+1) + pi
       dfm(:,nt,2) = pi - f(:,nt-1)
    end if
  end subroutine derm

  !> Calculate the derivative of \(r_\psi\) w.r.t. R and Z
  !>
  !> \(r_\psi\) is the "minor radius flux function"
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(leq_type), intent(in) :: self
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array is -ntm
    integer, intent(in) :: ntm
    !> If char = 'R', return the gradient of pressure instead.
    !> Note that this is useless here!
    character(1), intent(in) :: char
    !> Theta grid and rgrid (unused)
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    !> grad(:,1): derivative w.r.t. R as function of theta
    !> grad(:,2): derivative w.r.t  Z as function of theta
    real, dimension(-ntm:,:), intent(out) :: grad
    !> Value of psi/(a^2 B_a) on the flux surface where we want the gradient.
    !> Note: never used!
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, -self%nt:self%nt, 2) :: dcart
    integer :: i

    select case(char)
    case('P', 'R')
       dcart = self%dpcart
    case('T')
       dcart = self%dtcart  ! Never reached: bgradient always used
    end select

    do i = -nth_used, -1
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,2), 'Z')
       if (char == 'T') then
          grad(i, 1) = -grad(i, 1)
          grad(i, 2) = -grad(i, 2)
       end if
    end do

    do i = 0, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nr, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i = -nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5*self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5*self%beta_0
       end do
    end if
  end subroutine gradient

  !> Calculate the derivative of \(r_\psi\) in the Bishop coordinate system.
  !> If `char=='T'`, calculate the gradient of \(\theta\) instead
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(leq_type), intent(in) :: self
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array is -ntm
    integer, intent(in) :: ntm
    !> If char = 'R', return the gradient of pressure instead.
    !> Note that this is useless here!
    character(1), intent(in) :: char
    !> Theta grid and rgrid (unused)
    real, dimension(-ntm:), intent(in) :: theta, rgrid
    !> grad(:,1): derivative w.r.t. R as function of theta
    !> grad(:,2): derivative w.r.t  Z as function of theta
    real, dimension(-ntm:,:), intent(out) :: grad
    !> Value of psi/(a^2 B_a) on the flux surface where we want the gradient.
    !> Note: never used!
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, -self%nt:self%nt, 2) ::  dbish
    integer :: i

    select case(char)
    case('P')
       dbish = self%dpbish
    case('R')
       dbish = self%dpbish  ! dpcart is correct for 'R' -- never reached: other calculation done instead
    case('T')
       dbish = self%dtbish
    end select

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dbish(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dbish(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nr, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i = -nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5*self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5*self%beta_0
       enddo
    endif
  end subroutine bgradient

  !> Calculates `fstar` which is `f` interpolated at the location `theta_in`
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use constants, only: pi
    implicit none
    class(leq_type), intent(in) :: self
    real, intent(in) :: r_in !Ignored
    !> Desired theta location
    real, intent(in) :: theta_in
    !> Array to interpolate
    real, dimension(:,-self%nt:), intent(in) :: f
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
    real, dimension(-self%nt:self%nt) :: mtheta

    tp = mod2pi(theta_in)

    ! get theta points on the local theta mesh
    ! Note: mtheta must be the same points as the `t` variable in leq_init
    mtheta = [( real(j)*pi/real(self%nt), j=-self%nt,self%nt)]

    ! note that theta(1)=0 for local_eq theta
    jstar = -self%nt-1
    do j = -self%nt, self%nt
       if (tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar = j-1
          exit
       end if
       if(jstar /= -self%nt-1) write(*,*) 'exit error j'
    end do

    ! treat theta = pi separately
    if (jstar == -self%nt-1) then
       jstar = self%nt-1
       dt = mtheta(jstar+1) - mtheta(jstar)
       st = 0.
    end if

    ! use opposite area stencil to interpolate
    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt
    fstar = fstar / (mtheta(jstar+1)-mtheta(jstar))
  end function eqitem

  !> Convert `dfm` from index-space derivative to derivative on
  !> cylindrical coordinates
  !>
  !> - Gradient in the first coordinate is in `(:,:,1)`
  !> - Gradient in the second coordinate is in `(:,:,2)`
  pure subroutine eqdcart(self, dfm, dfcart)
    implicit none
    class(leq_type), intent(in) :: self
    !> Index-space derivative
    real, dimension(self%nr,-self%nt:self%nt,2), intent(in)  :: dfm
    !> Cylindrical coordinates derivative
    real, dimension(self%nr,-self%nt:self%nt,2), intent(out) :: dfcart
    real, dimension(self%nr,-self%nt:self%nt) :: denom
    integer :: i, j
    denom = self%drm(:,:,1) * self%dzm(:,:,2) - self%drm(:,:,2) * self%dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * self%dzm(:,:,2) - self%dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * self%drm(:,:,2) + self%drm(:,:,1) * dfm(:,:,2)

    ! FIXME: document why this is only on the range [2,nr], or
    ! simplify. Likely copy-paste from geq etc. See ceq for note
    do j = -self%nt, self%nt
       do i = 2, self%nr
          dfcart(i,j,:) = dfcart(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdcart

  !> Convert `dcart` from cylindrical coordinates derivative to
  !> derivative on Bishop coordinates
  !>
  !> - Gradient in the first coordinate is in `(:,:,1)`
  !> - Gradient in the second coordinate is in `(:,:,2)`
  pure subroutine eqdbish(self, dcart, dbish)
    implicit none
    class(leq_type), intent(in) :: self
    !> Cylindrical coordinates derivative
    real, dimension(self%nr, -self%nt:self%nt, 2), intent(in) :: dcart
    !> Bishop coordinates derivative
    real, dimension(self%nr, -self%nt:self%nt, 2), intent(out) :: dbish
    real, dimension(self%nr, -self%nt:self%nt) :: denom
    integer :: i, j
    denom = sqrt(self%dpcart(:,:,1)**2 + self%dpcart(:,:,2)**2)
    dbish(:,:,1) = dcart(:,:,1) * self%dpcart(:,:,1) + dcart(:,:,2) * self%dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * self%dpcart(:,:,2) + dcart(:,:,2) * self%dpcart(:,:,1)

    do j = -self%nt, self%nt
       do i = 2, self%nr
          dbish(i,j,:) = dbish(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdbish

  !> Return the inverse of the major radius at \((r, \theta)\)
  real function invR (self, r, theta)
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%Rpos(r, theta)
  end function invR

  !> Return the major radius at \((r, \theta)\)
  real function Rpos (self, r, theta)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: Zpos

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (self%surf%geoType)
       case (0)
          call self%geo_Miller(r, theta, Rpos, Zpos)
       case (1)
          call self%geo_global(r, theta, Rpos, Zpos)
       case (2)
          call self%geo_generalizedEllipticity(r, theta, Rpos, Zpos)
       case (3)
          call self%geo_FourierSeries(r, theta, Rpos, Zpos)
    end select
  end function Rpos

  !> Return the axial coordinate at \((r, \theta)\)
  real function Zpos (self, r, theta)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: Rpos

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (self%surf%geoType)
       case (0)
          call self%geo_Miller(r, theta, Rpos, Zpos)
       case (1)
          call self%geo_global(r, theta, Rpos, Zpos)
       case (2)
          call self%geo_generalizedEllipticity(r, theta, Rpos, Zpos)
       case (3)
          call self%geo_FourierSeries(r, theta, Rpos, Zpos)
    end select
  end function Zpos

  !> Calculate the cylindrical \((R, Z)$$ at $$(r, \theta)\) for the
  !> generalised Miller model
  !>
  !> See section 3.2.1 of Ball MIT Masters thesis or section 2.1.1 of
  !> "GS2 analytic geometry specification"
  pure subroutine geo_Miller(self, r, theta, Rpos, Zpos)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent(in) :: r, theta
    real, intent(out) :: Rpos, Zpos
    real :: dr, thAdj
    real :: Rcirc, Rcircp, Relong, Relongp, RelongTilt, RelongTiltp, Rtri, Rtrip, RtriTilt, RtriTiltp, Rfinal, Rfinalp
    real :: Zcirc, Zcircp, Zelong, Zelongp, ZelongTilt, ZelongTiltp, Ztri, Ztrip, ZtriTilt, ZtriTiltp, Zfinal, Zfinalp

    dr = r - self%surf%r

    thAdj = theta + self%thetaShift

    Rcirc = self%surf%r*cos(self%surf%thn - self%surf%thm - thAdj)
    Zcirc = -self%surf%r*sin(self%surf%thn - self%surf%thm - thAdj)
    Rcircp = cos(self%surf%thn - self%surf%thm - thAdj)
    Zcircp = -sin(self%surf%thn - self%surf%thm - thAdj)

    Relong = Rcirc
    Zelong = -self%surf%r*(self%surf%delm - 1)*sin(self%surf%thn - self%surf%thm - thAdj) + Zcirc
    Relongp = Rcircp
    Zelongp = -((self%surf%delm - 1) + self%surf%r*self%surf%delmp)*sin(self%surf%thn - self%surf%thm - thAdj) + Zcircp

    RelongTilt = cos(self%surf%thn - self%surf%thm)*Relong - sin(self%surf%thn - self%surf%thm)*Zelong
    ZelongTilt = sin(self%surf%thn - self%surf%thm)*Relong + cos(self%surf%thn - self%surf%thm)*Zelong
    RelongTiltp = cos(self%surf%thn - self%surf%thm)*Relongp - sin(self%surf%thn - self%surf%thm)*Zelongp
    ZelongTiltp = sin(self%surf%thn - self%surf%thm)*Relongp + cos(self%surf%thn - self%surf%thm)*Zelongp

    Rtri = -self%surf%r*cos(thAdj) + self%surf%r*cos(sin(thAdj)*self%surf%deln + thAdj) + RelongTilt
    Ztri = ZelongTilt
    Rtrip = -cos(thAdj) + cos(sin(thAdj)*self%surf%deln + thAdj) &
            - self%surf%r*sin(sin(thAdj)*self%surf%deln + thAdj)*sin(thAdj)*self%surf%delnp + RelongTiltp
    Ztrip = ZelongTiltp

    RtriTilt = cos(self%surf%thn)*Rtri + sin(self%surf%thn)*Ztri
    ZtriTilt = -sin(self%surf%thn)*Rtri + cos(self%surf%thn)*Ztri
    RtriTiltp = cos(self%surf%thn)*Rtrip + sin(self%surf%thn)*Ztrip
    ZtriTiltp = -sin(self%surf%thn)*Rtrip + cos(self%surf%thn)*Ztrip

    Rfinal = RtriTilt + self%surf%Rmaj
    Zfinal = ZtriTilt
    Rfinalp = RtriTiltp + self%surf%sHorz
    Zfinalp = ZtriTiltp + self%surf%sVert

    Rpos = Rfinal + dr*Rfinalp
    Zpos = Zfinal + dr*Zfinalp

  end subroutine geo_Miller

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> global MHD model
  !>
  !> See section 2.1.2 of "GS2 analytic geometry specification"
  pure subroutine geo_global(self, r, theta, Rpos, Zpos)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: rho0, thAdj, rCyl, Rc, Zc
    real :: Cm, Cn, psiN, xAng, yAng, rAng

    rho0 = r/self%surf%aSurf

    thAdj = theta + self%thetaShift

    Cm=(self%surf%delm**2-1)/(self%surf%delm**2+1)
    Cn=(self%surf%deln**2-1)/(self%surf%deln**3+1)

    psiN=rho0**2+Cm*rho0**2+Cn*rho0**3

    if (abs(Cn*cos(3*(thAdj+self%surf%thn)))<10e-6) then
      rCyl=self%surf%aSurf*sqrt(psiN/(1+Cm*cos(2*(thAdj+self%surf%thm))))
    else
      xAng=2*(1+Cm*cos(2*(thAdj+self%surf%thm)))**3-27*Cn**2*psiN*cos(3*(thAdj+self%surf%thn))**2
      yAng=3*Cn*sqrt(3*psiN)*cos(3*(thAdj+self%surf%thn))*sqrt(xAng+2*(1+Cm*cos(2*(thAdj+self%surf%thm)))**3)

      rAng=(1/3.0)*atan2(yAng,xAng)

      rCyl=self%surf%aSurf*(1+Cm*cos(2*(thAdj+self%surf%thm)))/(3*Cn*cos(3*(thAdj+self%surf%thn))) &
           *(cos(rAng)+sqrt(3.0)*sin(rAng)-1)
    end if

    Rc=self%surf%sHorz-(self%surf%sHorz-self%surf%Rmaj)*rho0**2
    Zc=self%surf%sVert-(self%surf%sVert)*rho0**2

    Rpos=Rc+rCyl*cos(thAdj)
    Zpos=Zc+rCyl*sin(thAdj)

  end subroutine geo_global

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> generalised elongation model
  !>
  !> See section 5.1.2 of Ball Oxford PhD thesis or section 2.1.3 of
  !> "GS2 analytic geometry specification"
  pure subroutine geo_generalizedEllipticity(self, r, theta, Rpos, Zpos)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - self%surf%r

    thAdj = theta + self%thetaShift

    rCyl=self%surf%r * (1 &
        + (-1+self%surf%delm/sqrt(1+(self%surf%delm**2-1)*cos(self%surf%mMode*(thAdj+self%surf%thm)/2.0)**2)) &
        + (-1+self%surf%deln/sqrt(1+(self%surf%deln**2-1)*cos(self%surf%nMode*(thAdj+self%surf%thn)/2.0)**2)))
    rCylp=1 &
     +(-1+(1+(self%surf%delm**2-1)*cos(self%surf%mMode*(thAdj+self%surf%thm)/2.0)**2)**(-0.5)*(self%surf%delm+self%surf%r*self%surf%delmp &
      *(1-(self%surf%delm*cos(self%surf%mMode*(thAdj+self%surf%thm)/2.0))**2/(1+(self%surf%delm**2-1)*cos(self%surf%mMode*(thAdj+self%surf%thm)/2.0)**2)))) &
     +(-1+(1+(self%surf%deln**2-1)*cos(self%surf%nMode*(thAdj+self%surf%thn)/2)**2)**(-0.5)*(self%surf%deln+self%surf%r*self%surf%delnp &
      *(1-(self%surf%deln*cos(self%surf%nMode*(thAdj+self%surf%thn)/2.0))**2/(1+(self%surf%deln**2-1)*cos(self%surf%nMode*(thAdj+self%surf%thn)/2.0)**2))))

    Rfinal=rCyl*cos(thAdj)+self%surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+self%surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+self%surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_generalizedEllipticity

  !> Calculate the cylindrical \((R, Z)\) at \((r, \theta)\) for the
  !> Fourier model
  !>
  !> See section 5.1.1 of Ball Oxford PhD thesis or section 2.1.4 of
  !> "GS2 analytic geometry specification"
  pure subroutine geo_FourierSeries(self, r, theta, Rpos, Zpos)
    implicit none
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - self%surf%r

    thAdj = theta + self%thetaShift

    rCyl=self%surf%r * (1+((self%surf%delm-1)/2) * (1-cos(self%surf%mMode*(thAdj+self%surf%thm))) + ((self%surf%deln-1)/2) * (1-cos(self%surf%nMode*(thAdj+self%surf%thn))))
    rCylp=1 + (((self%surf%delm-1)/2) + (self%surf%r*self%surf%delmp/2)) * (1-cos(self%surf%mMode*(thAdj+self%surf%thm))) + (((self%surf%deln-1)/2) + (self%surf%r*self%surf%delnp/2)) * (1-cos(self%surf%nMode*(thAdj+self%surf%thn)))

    Rfinal=rCyl*cos(thAdj)+self%surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+self%surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+self%surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_FourierSeries

  !> FIXME : Add documentation
  real function psi (self, r, theta)
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    if (r==0.) then
       psi = self%psi_0
    else if (r==1. .and. theta == 0.0) then
       psi = self%psi_a
    else
       psi = r - self%surf%r
    end if
  end function psi

  !> FIXME : Add documentation | Should this be equivalent to psi as in others?
  real function rfun (self, r, theta)
    class(leq_type), intent(in) :: self
    real, intent (in) :: r, theta
    rfun = r
  end function rfun

  !> Put `theta` into the range `[-pi, pi]`
  !>
  !> This version also sets `theta == +(-)pi +(-) theta_tol` to be
  !> exactly `+(-)pi`, where `theta_tol == 1.e-6`
  !>
  !> General purpose, could be moved into utils? Could also drop the
  !> tolerance?
  !>
  !> @note Currently different to the form in geo_utils
  pure real function mod2pi (theta)
    use constants, only: pi
    real, intent(in) :: theta
    real, parameter :: theta_tol = 1.e-6
    if(theta <= pi .and. theta >= -pi) then
      mod2pi = theta
    else if(theta - theta_tol <= pi .and. theta >= -pi) then
      mod2pi = pi
    else if(theta <= pi .and. theta + theta_tol >= -pi) then
      mod2pi = -pi
    else
      mod2pi = modulo(theta + pi, 2.*pi) - pi
    end if
  end function mod2pi

  !> FIXME : Add documentation
  real function diameter (self, rp)
    class(leq_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = 2. * rp
  end function diameter

  !> Return the major radius of the centre of the flux surface
  real function rcenter(self, rp)
    class(leq_type), intent(in) :: self
    real, intent(in) :: rp
    ! If using the global geometry, Rmaj is the center of the shaped flux surface
    ! (i.e. the surface at rhoc=aSurf), so here we calculate the major radial
    ! location of the center of the flux surface of interest ! JRB
    if (self%surf%geoType==1) then
       rcenter = self%surf%Rmaj+(1-(self%surf%r/self%surf%aSurf)**2)*self%surf%sHorz
    else
       rcenter = self%surf%Rmaj
    end if
  end function rcenter

  !> FIXME : Add documentation
  real function dbtori (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    dbtori = 1.
  end function dbtori

  !> FIXME : Add documentation
  real function btori (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    btori = self%surf%Rgeo
  end function btori

  !> FIXME : Add documentation
  real function qfun (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    qfun = self%surf%q
  end function qfun

  !> FIXME : Add documentation
  real function pfun (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    pfun = 0.5 * self%beta_0
  end function pfun

  !> FIXME : Add documentation
  real function dpfun (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    dpfun = -1.
  end function dpfun

  !> FIXME : Add documentation
  real function betafun (self, pbar)
    class(leq_type), intent(in) :: self
    real, intent(in) :: pbar
    betafun = self%beta_0
  end function betafun

  !> FIXME : Add documentation
  real function dpdrhofun()
    dpdrhofun = 0
  end function dpdrhofun
end module leq
