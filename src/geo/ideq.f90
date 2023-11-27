!> FIXME : Add documentation
module ideq
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi

  implicit none

  private

  public :: ideq_type

  type, extends(abstract_eqfile_geo_type) :: ideq_type
     private
   contains
     procedure :: initialise => idfitin
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     procedure :: alloc_special_arrays => alloc_special_arrays
     procedure :: dealloc_special_arrays => dealloc_special_arrays
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
    !> Not really the diameter in this case.  Instead, return the
    !> normalized minor radius, measured inside the ring, in the plane of
    !> the ring, starting at the ring and going inward.
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

    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * pi / float(self%nt-1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    ! recover some memory:
    do i = 1, d%nr
       deallocate (d%grid(i)%R)
       deallocate (d%grid(i)%Z)
       deallocate (d%grid(i)%B)
       deallocate (d%grid(i)%theta)
    end do
    deallocate (d%grid)

    self%initialised = .true.
    self%has_full_theta_range = .false.
  end subroutine idfitin

  !> FIXME : Add documentation
  pure subroutine alloc_special_arrays(self)
    implicit none
    class(ideq_type), intent (in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%psi_bar(nr), self%fp(nr), self%qsf(nr), self%beta(nr))
    allocate(self%rc(nr), self%diam(nr))
  end subroutine alloc_special_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(ideq_type), intent(in out) :: self
    deallocate(self%psi_bar, self%fp, self%qsf, self%beta, self%rc, self%diam)
  end subroutine dealloc_special_arrays

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
    case('P', 'R')
       dcart = self%dpcart
    case('T')
       dcart = self%dtcart
    end select

    do i = -nth_used, -1
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,2), 'Z')
       if(char == 'T') then
          grad(i,1) = -grad(i, 1)
          grad(i,2) = -grad(i, 2)
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
    case('P', 'R')
       dbish = self%dpbish
    case('T')
       dbish = self%dtbish
    end select

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dbish(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dbish(:,:,2), 'Z')
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
       do i = -nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5 * self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5 * self%beta_0
       end do
    end if
  end subroutine bgradient

  !> FIXME : Add documentation
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    class(ideq_type), intent(in) :: self
    real, intent(in) :: r_in
    real, intent(in) :: theta_in
    real, dimension(:,:), intent(in) :: f
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
    mtheta = self%eqth(1, j)

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
  end function eqitem
end module ideq
