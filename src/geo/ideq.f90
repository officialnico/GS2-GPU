!> FIXME : Add documentation
module ideq
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi, EQFILE_LENGTH
  use splines, only: spline

  implicit none

  private

  public :: ideq_type

  type, extends(abstract_eqfile_geo_type) :: ideq_type
     private
     real, allocatable, dimension (:,:)   :: B_psi
     real, allocatable, dimension (:,:,:) :: dbm, dbcart, dbbish
     character(len = EQFILE_LENGTH) :: filename
     integer :: nr, nt
     real :: psi_0, psi_a, B_T, beta_0
     real :: R_mag, Z_mag, aminor
     real, allocatable, dimension (:)     :: eqpsi, psi_bar, fp, qsf
     real, allocatable, dimension (:)     :: beta, pressure, diam, rc
     real, allocatable, dimension (:,:)   :: R_psi, Z_psi
     real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
     real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
     real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish
     type(spline) :: diam_spl, rcenter_spl, btori_spl, qsf_spl, pressure_spl, beta_spl
   contains
     procedure :: initialise => idfitin
     procedure, private :: calculate_gradients => idfit_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
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
     procedure :: rfun => psi
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure, private :: derm => derm
     procedure, private :: eqitem => eqitem
     procedure, private :: eqdcart => eqdcart
     procedure, private :: eqdbish => eqdbish
     procedure, private :: alloc_arrays => alloc_arrays
     procedure, private :: dealloc_arrays => dealloc_arrays
  end type ideq_type

contains

  !> FIXME : Add documentation
  subroutine idfitin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use constants, only: twopi, pi
    use splines, only: new_periodic_spline, delete_periodic_spline, &
         periodic_splint, periodic_spline
    use geo_utils, only: sort, geo_input_type
    implicit none
    class(ideq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    !> FIXME : Add documentation
    type :: grid_type
       integer :: nt
       real, dimension (:), allocatable :: R
       real, dimension (:), allocatable :: Z
       real, dimension (:), allocatable :: B
       real, dimension (:), allocatable :: theta
       type (periodic_spline) :: spl
    end type grid_type

    !> FIXME : Add documentation
    type :: eq_type
       integer :: nr
       real, dimension (:), allocatable :: pbar
       real, dimension (:), allocatable :: psi
       real, dimension (:), allocatable :: pressure
       real, dimension (:), allocatable :: volume
       type (grid_type), dimension (:), allocatable :: grid
    end type eq_type
    type (eq_type) :: d
    real, intent (out) :: rmaj, avgrmid, psi_a_out, psi_0_out, B_T0
    integer, intent (out) :: nthg
    real :: xdum, I_ring, dpsidr, dpsidz, psi_N
    integer :: i, j, idum, jj, jmax
    character (200) :: line
    logical, dimension(:), allocatable :: has_duplicate_point
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = self%nt
       return
    end if

    ! Need to generalize initialization condition if equilibrium changes
    self%nt = size(inputs%theta)

    self%filename = trim(adjustl(inputs%eqfile))

    open(unit = 5, file = self%filename, status = 'old', form = 'formatted')

    ! Read the data

    read (5, fmt="(a)") line
    read (5, fmt="(a)") line
    read (5,*) self%R_mag, self%Z_mag, I_ring
    read (5, fmt="(a)") line
    read (5, *) self%nr
    read (5, fmt="(a)") line
    read (5, fmt="(a)") line

    d%nr = self%nr
    allocate (d%pbar     (d%nr))
    allocate (d%psi      (d%nr))
    allocate (d%pressure (d%nr))
    allocate (d%volume   (d%nr))
    allocate (d%grid     (d%nr))

    ! Read in each surface
    do i = 1, d%nr
       read(5,*) idum, d%pbar(i), d%psi(i), d%pressure(i), xdum, d%volume(i), d%grid(i)%nt
       d%grid(i)%nt = d%grid(i)%nt-1
       jmax = d%grid(i)%nt
       allocate (d%grid(i)%R     (jmax))
       allocate (d%grid(i)%Z     (jmax))
       allocate (d%grid(i)%B     (jmax))
       allocate (d%grid(i)%theta (jmax))
       do j = 1, jmax
          read (5,*) d%grid(i)%R(j), d%grid(i)%Z(j), dpsidR, dpsidZ
          ! make B:
          d%grid(i)%B(j) = sqrt(dpsidR**2 + dpsidZ**2) / (d%grid(i)%R(j)*twopi)
          d%grid(i)%theta(j) = atan2 (d%grid(i)%Z(j) - self%Z_mag, &
               (d%grid(i)%R(j) - self%R_mag))
       end do
       read (5,*) xdum, xdum, xdum, xdum
    end do

    ! different from tokamak case here, psi_0 is the value of psi on
    ! the innermost flux surface that is kept in equilibrium and psi_a
    ! is the last flux surface kept.
    self%psi_0 = d%psi(1)
    self%psi_a = d%psi(d%nr)

    do i = 1, d%nr
       call sort (d%grid(i)%theta, d%grid(i)%R, d%grid(i)%Z, d%grid(i)%B)
    end do

    ! Spline the data onto a regular grid for later use.

    !     R_psi(1:nr, 1:nt)
    !     Z_psi(1:nr, 1:nt)
    !     B_psi(1:nr, 1:nt)
    !
    ! nt determined from input file by user

    ! need to define our theta grid.  still not right
    call self%alloc_arrays

    allocate(has_duplicate_point(d%nr))

    ! Need to determine if we have a duplicate grid point as the spline routines
    ! require that we do not include the duplicate point
    do i = 1, d%nr
       has_duplicate_point(i) = twopi - &
            (d%grid(i)%theta(d%grid(i)%nt) - d%grid(i)%theta(1)) <= epsilon(twopi)
    end do

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%R(:), &
            twopi, d%grid(i)%spl, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%R_psi (i, jj) = periodic_splint (inputs%theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%Z(:), &
            twopi, d%grid(i)%spl, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%Z_psi (i, jj) = periodic_splint (inputs%theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%B(:), &
            twopi, d%grid(i)%spl, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%B_psi (i, jj) = periodic_splint (inputs%theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    rmaj = 1.
    avgrmid = -1.0 !< Historically no value given

    ! B at center of ring is mu_0 I / (2 a)
    self%B_T = 4.*pi*1.e-7*I_ring*0.5 / self%R_mag

    write(*,*) 'B_T0 = ',self%B_T
    B_T0 = self%B_T

    psi_N = self%B_T * self%R_mag**2
    ! Why do we only normalise the outputs and not the stored values?
    psi_a_out = self%psi_a / psi_N
    psi_0_out = self%psi_0 / psi_N
    self%eqpsi = d%psi / psi_N

    self%R_psi = self%R_psi / self%R_mag
    self%Z_psi = self%Z_psi / self%R_mag
    self%B_psi = self%B_psi/B_T0
    self%rc = 0.5*(self%R_psi(:, 1) + self%R_psi(:, (self%nt-1)/2))
    ! need to define the rho_d grid (called diam)
    self%diam = 1. - self%R_psi(:,1)


    self%qsf = 0.
    self%fp = 0.
    ! Note beta and pressure were allocated but not set historically.
    ! As such now forcing these to zero, but we _do_ read in pressure
    ! so we could probably actually set this?
    self%beta = 0.
    self%pressure= 0.
    self%beta_0 = self%beta(1)

    nthg = self%nt

    ! recover some memory:
    do i = 1, d%nr
       deallocate (d%grid(i)%R)
       deallocate (d%grid(i)%Z)
       deallocate (d%grid(i)%B)
       deallocate (d%grid(i)%theta)
    end do
    deallocate (d%grid)

    self%initialised = .true.
  end subroutine idfitin

  !> FIXME : Add documentation
  pure subroutine alloc_arrays(self)
    implicit none
    class(ideq_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%eqpsi(nr), self%psi_bar(nr), self%fp(nr), self%qsf(nr), &
         self%beta(nr), self%pressure(nr), self%rc(nr), self%diam(nr))
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt), self%B_psi(nr, nt))
    allocate(self%drm(nr, nt, 2), self%dzm(nr, nt, 2), self%dbm(nr, nt, 2), &
         self%dbtm(nr, nt, 2), self%dpm(nr, nt, 2), self%dtm(nr, nt, 2))
    allocate(self%dpcart(nr, nt, 2), self%dbcart(nr, nt, 2), self%dtcart(nr, nt, 2))
    allocate(self%dpbish(nr, nt, 2), self%dbbish(nr, nt, 2), self%dtbish(nr, nt, 2))
    allocate(self%dbtcart(nr, nt, 2), self%dbtbish(nr, nt, 2))
  end subroutine alloc_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_arrays(self)
    implicit none
    class(ideq_type), intent(in out) :: self
    deallocate(self%eqpsi, self%psi_bar, self%fp, self%qsf, self%beta, &
         self%pressure, self%rc, self%diam)
    deallocate(self%R_psi, self%Z_psi, self%B_psi)
    deallocate(self%drm, self%dzm, self%dbm, self%dbtm, self%dpm, self%dtm)
    deallocate(self%dpcart, self%dbcart, self%dtcart, self%dbtcart)
    deallocate(self%dpbish, self%dbbish, self%dtbish, self%dbtbish)
  end subroutine dealloc_arrays

  !> FIXME : Add documentation
  subroutine idfit_init(self)
    use constants, only: pi
    implicit none
    class(ideq_type), intent(in out) :: self
    real, dimension(self%nr, self%nt) :: eqpsi1, eqth
    integer :: j

    do j = 1, self%nt
       eqpsi1(:,j) = self%eqpsi
       eqth(:,j) = (j-1) * pi / float(self%nt-1)
    end do

    call self%derm(eqth, self%dtm, 'T')
    call self%derm(self%R_psi, self%drm, 'E')
    call self%derm(self%Z_psi, self%dzm, 'O')
    call self%derm(self%B_psi, self%dbm, 'E')
    call self%derm(eqpsi1, self%dpm, 'E')

    ! grad(psi) in cartesian form
    call self%eqdcart(self%dpm, self%dpcart)
    ! grad(psi) in Bishop form
    call self%eqdbish(self%dpcart, self%dpbish)

    ! grad(B) in cartesian form
    call self%eqdcart(self%dbm, self%dbcart)
    ! grad(B) in Bishop form
    call self%eqdbish(self%dbcart, self%dbbish)

    ! grad(BT) in cartesian form
    self%dbtcart = 0.
    ! grad(BT) in Bishop form
    self%dbtbish = 0.

    ! grad(theta) in cartesian form
    call self%eqdcart(self%dtm, self%dtcart)
    ! grad(theta) in Bishop form
    call self%eqdbish(self%dtcart, self%dtbish)
  end subroutine idfit_init

  subroutine setup_splines(self)
    use splines, only: new_spline
    implicit none
    class(ideq_type), intent(in out) :: self
    call new_spline(self%nr, self%eqpsi, self%diam, self%diam_spl)
    call new_spline(self%nr, self%eqpsi, self%rc, self%rcenter_spl)
    call new_spline(self%nr, self%psi_bar, self%fp, self%btori_spl)
    call new_spline(self%nr, self%psi_bar, self%qsf, self%qsf_spl)
    call new_spline(self%nr, self%psi_bar, self%pressure, self%pressure_spl)
    call new_spline(self%nr, self%psi_bar, self%beta, self%beta_spl)
  end subroutine setup_splines

  subroutine delete_splines(self)
    use splines, only: delete_spline
    implicit none
    class(ideq_type), intent(in out) :: self
    call delete_spline(self%diam_spl)
    call delete_spline(self%rcenter_spl)
    call delete_spline(self%btori_spl)
    call delete_spline(self%qsf_spl)
    call delete_spline(self%pressure_spl)
    call delete_spline(self%beta_spl)
  end subroutine delete_splines

  !> FIXME : Add documentation
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(ideq_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nr, nt
    nr = self%nr ; nt = self% nt

    i = 1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i = nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

!!! this is a bad assumption
    ! assume up-down symmetry for now:

    select case (char)
    case ('E')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=1
       dfm(:,j,2) = f(:,j+1)

       j=nt
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i = 2, nr-1
       dfm(i,:,1) = 0.5*(f(i+1,:)-f(i-1,:))
    end do

    do j = 2, nt-1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j-1))
    end do
  end subroutine derm

  !> FIXME : Add documentation
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(ideq_type), intent(in) :: self
    integer, intent (in) :: nth_used, ntm
    character(1), intent (in) :: char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real, dimension (self%nr, self%nt, 2) :: dcart
    real, dimension (1) :: aa, daa, rpt
    integer :: i

    select case(char)
    case('B')
       dcart = self%dbcart
    case('D')  ! diagnostic
       dcart = self%dbtcart
    case('P')
       dcart = self%dpcart
    case('R')
       dcart = self%dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = self%dtcart
    end select

    do i = -nth_used, -1
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), grad(i, 1), 'R')
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,2), grad(i, 2), 'Z')
       if(char == 'T') then
          grad(i,1) = -grad(i, 1)
          grad(i,2) = -grad(i, 2)
       end if
    end do

    do i = 0, nth_used
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), grad(i, 1), 'R')
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,2), grad(i, 2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nr, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i = -nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5 * self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5 * self%beta_0
       end do
    end if
  end subroutine gradient

  !> FIXME : Add documentation
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(ideq_type), intent(in) :: self
    integer, intent (in) :: nth_used, ntm
    character(1), intent (in) :: char
    real, dimension (-ntm:), intent (in)  ::  rgrid, theta
    real, dimension (-ntm:,:), intent (out) :: grad
    real, intent(in) :: rp
    real, dimension (1) :: aa, daa, rpt
    real, dimension(self%nr, self%nt, 2) ::  dbish
    integer :: i

    select case(char)
    case('B')
       dbish = self%dbbish
    case('D')  ! diagnostic
       dbish = self%dbtbish
    case('P')
       dbish = self%dpbish
    case('R')
       dbish = self%dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = self%dtbish
    end select

    do i = -nth_used, nth_used
       call self%eqitem(rgrid(i), theta(i), dbish(:,:,1), grad(i,1), 'R')
       call self%eqitem(rgrid(i), theta(i), dbish(:,:,2), grad(i,2), 'Z')
    end do

    if (char == 'T') then
       where (theta(-nth_used:nth_used) < 0.0)
          grad(-nth_used:nth_used,1) = -grad(-nth_used:nth_used,1)
          grad(-nth_used:nth_used,2) = -grad(-nth_used:nth_used,2)
       end where
    end if

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nr, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i =-nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5 * self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5 * self%beta_0
       enddo
    endif
  end subroutine bgradient

  !> FIXME : Add documentation
  subroutine eqitem(self, r_in, theta_in, f, fstar, char)
    use constants, only: pi
    use mp, only: mp_abort
    class(ideq_type), intent(in) :: self
    real, intent(in) :: r_in
    real, intent(in) :: theta_in
    real, dimension(:,:), intent(in) :: f
    real, intent(out) :: fstar
    character(1), intent(in) :: char
    integer :: i, j, istar, jstar
    real :: thet, sign, tp, tps, r
    real :: st, dt, sr, dr, rt
    real, dimension(size(f,2)) :: mtheta

    !Copy value to avoid any modification
    r=r_in

    ! check for axis evaluation
    if (r == self%eqpsi(1)) then
       write(*,*) 'no evaluation at axis allowed in eqitem'
       write(*,*) r, theta_in, self%eqpsi(1)
       call mp_abort('no evaluation at axis allowed in eqitem')
    end if

    ! allow psi(r) to be a decreasing function
    sign = 1.
    if (self%eqpsi(2) < self%eqpsi(1)) sign = -1.

    if (r < sign * self%eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r, sign, self%eqpsi(1)
       call mp_abort('r < Psi_0 in eqitem')
    end if

    ! find r on psi mesh

    ! disallow evaluations outside the plasma surface for now
    if (r == self%eqpsi(self%nr)) then
       rt = 0.9999999999 * r
       if (rt > self%eqpsi(self%nr)) rt = 1.0000000001*r
       r = rt
    end if

    if (r >= self%eqpsi(self%nr)) then
       write(*,*) 'No evaluation of eqitem allowed outside surface'
       write(*,*) r, theta_in, self%eqpsi(self%nr), sign
       call mp_abort('No evaluation of eqitem allowed outside surface')
    end if

    istar = 0
    do i = 2, self%nr
       if (r < sign * self%eqpsi(i)) then
          dr = r - sign * self%eqpsi(i-1)
          sr = sign * self%eqpsi(i) - r
          istar = i-1
          exit
       end if
    end do

    ! no gradients available at axis, so do not get too close
    if (istar == 1) then
       write(*,*) 'Too close to axis in eqitem'
       write(*,*) r, theta_in, self%eqpsi(1), self%eqpsi(2)
       call mp_abort('Too close to axis in eqitem')
    end if

    ! Now do theta direction
    thet = mod2pi(theta_in)

    ! assume up-down symmetry
    tp =abs(thet)
    tps = 1.
    if (char == 'Z' .and. tp >= 1.e-10) tps = thet/tp

    ! get thet on theta mesh
    do j = 1, self%nt
       mtheta(j) = (j-1)*pi / float(self%nt-1)
    end do

    ! note that theta(1)=0 for gen_eq theta
    jstar = -1
    do j = 1, self%nt
       if (jstar /= -1) cycle
       if (tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar = j-1
       end if
    end do

    ! treat theta = pi separately
    if (jstar == -1) then
       jstar = self%nt-1
       dt = mtheta(jstar+1) - mtheta(jstar)
       st = 0.
    end if

    ! use opposite area stencil to interpolate
    fstar = f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar * tps / (abs(self%eqpsi(istar+1) - self%eqpsi(istar))  * &
         (mtheta(jstar+1) - mtheta(jstar)))
  end subroutine eqitem

  !> FIXME : Add documentation
  pure subroutine eqdcart(self, dfm, dfcart)
    implicit none
    class(ideq_type), intent(in) :: self
    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (self%nr, self%nt) :: denom
    integer :: i, j
    denom = self%drm(:,:,1) * self%dzm(:,:,2) - self%drm(:,:,2) * self%dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * self%dzm(:,:,2) - self%dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * self%drm(:,:,2) + self%drm(:,:,1) * dfm(:,:,2)

    do j = 1, self%nt
       do i = 2, self%nr !Note this means we don't normalise i = 1 to denom?
          dfcart(i,j,:) = dfcart(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdcart

  !> FIXME : Add documentation
  pure subroutine eqdbish(self, dcart, dbish)
    implicit none
    class(ideq_type), intent(in) :: self
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j
    denom = sqrt(self%dpcart(:,:,1)**2 + self%dpcart(:,:,2)**2)
    dbish(:,:,1) = dcart(:,:,1) * self%dpcart(:,:,1) + dcart(:,:,2) * self%dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * self%dpcart(:,:,2) + dcart(:,:,2) * self%dpcart(:,:,1)

    do j = 1, self%nt
       do i = 2, self%nr !Note this means we don't normalise i = 1 to denom?
          dbish(i,j,:) = dbish(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdbish

  !> FIXME : Add documentation
  real function invR (self, r, theta)
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%Rpos(r, theta)
  end function invR

  !> FIXME : Add documentation
  real function Rpos (self, r, theta)
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%R_psi, Rpos, 'R')
  end function Rpos

  !> FIXME : Add documentation
  real function Zpos (self, r, theta)
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%Z_psi, Zpos, 'Z')
  end function Zpos

  !> FIXME : Add documentation
  real function psi (self, r, theta)
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent (in) :: r, theta
    if (r==0.) then
       psi = self%psi_0
    else if (r==1. .and. theta == 0.0) then
       psi = self%psi_a
    else
       psi = r
    end if
  end function psi

  !> FIXME : Add documentation
  !>
  !> Not really the diameter in this case.  Instead, return the
  !> normalized minor radius, measured inside the ring, in the plane of
  !> the ring, starting at the ring and going inward.
  real function diameter (self, rp)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = splint(rp, self%diam_spl)
  end function diameter

  !> FIXME : Add documentation
  real function rcenter (self, rp)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: rp
    rcenter = splint(rp, self%rcenter_spl)
  end function rcenter

  !> FIXME : Add documentation | Always zero.
  real function dbtori (self, pbar)
    use splines, only: dsplint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    dbtori = dsplint(pbar, self%btori_spl)/(self%psi_a-self%psi_0)
  end function dbtori

  !> FIXME : Add documentation | Always zero
  real function btori (self, pbar)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    btori = splint(pbar, self%btori_spl)
  end function btori

  !> FIXME : Add documentation | Always zero
  real function qfun (self, pbar)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    qfun = splint(pbar, self%qsf_spl)
  end function qfun

  !> FIXME : Add documentation | Always zero
  real function pfun (self, pbar)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    pfun = 0.5 * self%beta_0 * splint(pbar, self%pressure_spl)
  end function pfun

  !> FIXME : Add documentation
  real function dpfun (self, pbar)
    use splines, only: dsplint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    dpfun = 0.5 * self%beta_0 * dsplint(pbar, self%pressure_spl) / (self%psi_a-self%psi_0)
  end function dpfun

  !> FIXME : Add documentation
  real function betafun (self, pbar)
    use splines, only: splint
    implicit none
    class(ideq_type), intent(in) :: self
    real, intent(in) :: pbar
    if (pbar == 0.) then
       betafun = self%beta_0
       return
    end if
    betafun = splint(pbar, self%beta_spl)
  end function betafun
end module ideq
