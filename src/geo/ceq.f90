!> Reads in the geometry using a CHEASE output file.
!> The CHEASE output file is read using the helper module read_chease.
module ceq
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi, EQFILE_LENGTH
  use splines, only: spline

  implicit none

  public :: Hahm_Burrell, ceq_type

  private

  type, extends(abstract_eqfile_geo_type) :: ceq_type
     private
     character(len = EQFILE_LENGTH) :: filename
     integer :: nr, nt
     real :: psi_0, psi_a, B_T, I_0, beta_0
     real :: R_mag, Z_mag, aminor
     real, allocatable, dimension (:)     :: rho_d, eqpsi, psi_bar, fp, qsf
     real, allocatable, dimension (:)     :: beta, pressure, diam, rc
     real, allocatable, dimension (:,:)   :: R_psi, Z_psi !, B_psi
     real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm!, dbm
     real, allocatable, dimension (:,:,:) :: dpcart, dbcart, dtcart, dbtcart
     real, allocatable, dimension (:,:,:) :: dpbish, dbbish, dtbish, dbtbish
     type(spline) :: diam_spl, rcenter_spl, btori_spl, qsf_spl, pressure_spl, beta_spl
   contains
     procedure :: initialise_eqfile => ceqin
     procedure :: calculate_gradients => ceq_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: finalise => ceq_finish
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
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure, private :: derm => derm
     procedure, private :: eqitem => eqitem
     procedure, private :: eqdcart => eqdcart
     procedure, private :: eqdbish => eqdbish
     procedure, private :: alloc_arrays => alloc_arrays
     procedure, private :: dealloc_arrays => dealloc_arrays
     procedure, private :: set_ceq_from_chease => set_ceq_from_chease
  end type ceq_type

  logical, parameter :: debug = .true.
  logical :: skip_file_read = .false.
contains

  !> Convert from theta-index which is `-pi` to `pi`, to chi-index, which
  !> is `0` to `2*pi.`
  !> Assumes `nchi` is even?
  elemental integer function chease_chi_index(nchi, itheta)
    integer, intent(in) :: nchi, itheta
    ! when ichi = 0, itheta = (ntheta-1)/2 + 1
    ! when itheta = 0, ichi = ntheta / 2
    ! itheta  1 2 3 4 5 6 7 8
    ! ichi    5 6 7 8 1 2 3 4
    if (itheta > nchi/2) then
       chease_chi_index = itheta - nchi/2
    else
       chease_chi_index = itheta + nchi/2
    end if
  end function chease_chi_index

  !> FIXME : Add documentation
  subroutine ceqin(self, eqfile, psi_0_out, psi_a_out, rmaj, B_T0, &
       avgrmid, initeq, nthg, big_in)
    use read_chease, only: read_infile
    implicit none
    class(ceq_type), intent(in out) :: self
    character(len = EQFILE_LENGTH), intent(in) :: eqfile
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(in) :: initeq
    integer, intent(out) :: nthg
    integer, intent(in), optional :: big_in
    real :: f_N, psi_N

    ! Exit early if already initialised
    if(initeq == 0) then
      ! Ensure `nthg` gets set if exiting early or it may be uninitialised in
      ! the calling procedure. The value it gets set to must be consistent
      ! with the value set at the end of this subroutine when not exiting
      ! early or the calling procedure may end up with inconsistently sized
      ! arrays.
      nthg = self%nt
      return
    endif

    self%filename = trim(adjustl(eqfile))
    if (.not. skip_file_read) then ! Primarily for testing
       write (*,*)  'Reading CHEASE input file: ', self%filename
       !    Read the data
       call read_infile(self%filename)
       call self%set_ceq_from_chease
    end if
    !
    !     Normalize, rename quantities
    !

    avgrmid = self%aminor
    ! rmaj changes depending on iflux.. FIX!!
    ! see geometr... for Miller is the rgeo for the flux surface of
    ! interest... rgeo is rgeo of LCFS... check!!!
    ! actually when iflux = 1, it is the major radius of the magnetic axis,
    ! so below is OK.
    rmaj = self%R_mag
    B_T0 = abs(self%B_T)

    psi_N = B_T0 * avgrmid**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    psi_a_out = self%psi_a
    psi_0_out = self%psi_0
    self%eqpsi = self%eqpsi / psi_N

    f_N = B_T0*avgrmid
    self%fp = self%fp / f_N

    nthg = self%nt

    if (debug) then
       write (*,*) "Finished ceqin... imported CHEASE equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", self%aminor
       write (*,*) 'R_mag', self%R_mag
       write (*,*) 'B_T0', B_T0
       write (*,*) 'f_N', f_N
       write (*,*) 'nthg', nthg
       write (*,*) 'beta', self%beta_0
    end if

    self%diam = abs(self%R_psi(:, self%nt/2 + 1) - self%R_psi(:, 1))
    self%rc = 0.5*(self%R_psi(:, 1) + self%R_psi(:, self%nt/2 + 1))

    self%initialised = .true.
  end subroutine ceqin

  !> Set the ceq module-level variables from the read_chease module
  !> variables.
  !>
  !> This involves changing the theta grid from 0 to 2pi => -pi to pi,
  !> adding a point at theta == 0, and some normalisations
  pure subroutine set_ceq_from_chease(self)
    use constants, only: pi
    use read_chease, only: npsi_chease, nchi_chease, b0exp_chease
    use read_chease, only: psi_chease, f_chease, q_chease, p_chease
    use read_chease, only: r_chease, z_chease
    implicit none
    class(ceq_type), intent(in out) :: self
    real :: d, R_geo
    integer :: i, j

    !     nr == number of actual grid points in radial array
    self%nr = npsi_chease

    !     nt == number of theta grid points in theta eq grid
    self%nt = nchi_chease + 1

    call self%alloc_arrays()

    !     eqpsi(1:nr) == values of psi on the radial grid
    self%eqpsi = psi_chease

    !     fp(1:nr) == the function that satisfies
    !              B = fp grad zeta + grad zeta x grad psi
    self%fp = f_chease

    self%qsf = q_chease

    !     pressure(1:nr) == pressure profile on the radial grid,
    !     normalized to the value at the magnetic axis.
    self%pressure = p_chease

    !     psi_0 == value of psi at the magnetic axis
    !     psi_a == value of psi at the boundary of the plasma
    self%psi_0 = self%eqpsi(1)
    self%psi_a = self%eqpsi(self%nr)
    self%psi_bar = (self%eqpsi - self%psi_0) / (self%psi_a - self%psi_0)

    do j = 1, self%nt
      do i = 1, self%nr
        self%R_psi(i, j) = r_chease(i, chease_chi_index(nchi_chease, j))
      end do
    end do

    do j = 1, self%nt
      do i = 1, self%nr
        self%Z_psi(i,j) = z_chease(i, chease_chi_index(nchi_chease, j))
      end do
    end do

    !     aminor    == half diameter of last flux surface at elevation of Z_mag
    d = self%R_psi(self%nr, self%nt/2+1) - self%R_psi(self%nr, 1)
    R_geo = (self%R_psi(self%nr, self%nt/2+1) + self%R_psi(self%nr, 1))/2.
    self%aminor = d/2.

    self%B_T = self%fp(self%nr)/R_geo

    self%R_psi = self%R_psi / self%aminor
    self%Z_psi = self%Z_psi / self%aminor

    !     R_mag == R position of magnetic axis
    !     Z_mag == Z position of magnetic axis
    self%R_mag = self%R_psi(1,1)
    self%Z_mag = self%Z_psi(1,1)

    self%beta = 8. * pi * 1.e-7 * self%pressure / b0exp_chease**2  !not correct definition yet
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)
  end subroutine set_ceq_from_chease

  !> Allocate module-level arrays
  pure subroutine alloc_arrays(self)
    implicit none
    class(ceq_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%rho_d(nr), self%eqpsi(nr), self%psi_bar(nr), self%fp(nr), self%beta(nr), &
         self%pressure(nr), self%rc(nr), self%diam(nr), self%qsf(nr))
    ! rho_d = 0. ; eqpsi = 0. ; psi_bar = 0. ; fp = 0. ; beta = 0. ; pressure = 0.
    ! rc = 0. ; diam = 0. ; qsf = 0.
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt))
    ! R_psi = 0.  ; Z_psi = 0.
    allocate(self%drm(nr, nt, 2), self%dzm(nr, nt, 2), self%dbtm(nr, nt, 2), &
         self%dpm(nr, nt, 2), self%dtm(nr, nt, 2))
    ! drm = 0. ; dzm = 0. ; dbtm = 0. ; dpm = 0. ; dtm = 0.
    allocate(self%dpcart(nr, nt, 2), self%dtcart(nr, nt, 2), self%dbtcart(nr, nt, 2))
    allocate(self%dpbish(nr, nt, 2), self%dtbish(nr, nt, 2), self%dbtbish(nr, nt, 2))
    ! dpcart = 0. ; dtcart = 0. ; dbtcart = 0.
    ! dpbish = 0. ; dtbish = 0. ; dbtbish = 0.
  end subroutine alloc_arrays

  !> Deallocate module-level arrays
  pure subroutine dealloc_arrays(self)
    implicit none
    class(ceq_type), intent(in out) :: self
    if (allocated(self%rho_d)) deallocate(self%rho_d, self%eqpsi, self%psi_bar, self%fp, &
         self%beta, self%pressure, self%rc, self%diam, self%qsf)
    if (allocated(self%R_psi)) deallocate(self%R_psi, self%Z_psi)
    if (allocated(self%drm)) deallocate(self%drm, self%dzm, self%dbtm, self%dpm, self%dtm)
    if (allocated(self%dpcart)) deallocate(self%dpcart, self%dtcart, self%dbtcart)
    if (allocated(self%dpbish)) deallocate(self%dpbish, self%dtbish, self%dbtbish)
  end subroutine dealloc_arrays

  !> Finish this module; deallocate arrays
  subroutine ceq_finish(self)
    use read_chease, only : finish
    implicit none
    class(ceq_type), intent(in out) :: self
    self%initialised = .false.
    call self%dealloc_arrays
    call self%delete_splines
    call finish
  end subroutine ceq_finish

  !> Initialise CHEASE equilibrium module
  !>
  !> Assumes data has already been read in
  pure subroutine ceq_init(self)
    use constants, only: pi
    implicit none
    class(ceq_type), intent(in out) :: self
    real, dimension(self%nr, self%nt) :: eqpsi1, eqth, eqbtor
    integer :: j

    do j = 1, self%nt
       eqbtor(:, j) = self%fp / self%R_psi(:, j)
       eqpsi1(:, j) = self%eqpsi
       eqth(:, j) = (j-1)*2.*pi / real(self%nt-1) - pi
    end do

    call self%derm(eqth, self%dtm, 'T')
    call self%derm(self%R_psi, self%drm, 'E')
    call self%derm(self%Z_psi,  self%dzm, 'O')
    call self%derm(eqbtor, self%dbtm, 'E')
    call self%derm(eqpsi1, self%dpm, 'E')

    ! grad(psi) in cartesian form
    call self%eqdcart(self%dpm, self%dpcart)
    ! grad(psi) in Bishop form
    call self%eqdbish(self%dpcart, self%dpbish)

    ! grad(BT) in cartesian form
    call self%eqdcart(self%dbtm, self%dbtcart)
    ! grad(BT) in Bishop form
    call self%eqdbish(self%dbtcart, self%dbtbish)

    ! grad(theta) in cartesian form
    call self%eqdcart(self%dtm, self%dtcart)
    ! grad(theta) in Bishop form
    call self%eqdbish(self%dtcart, self%dtbish)
  end subroutine ceq_init

  subroutine setup_splines(self) !Identical to geq
    use splines, only: new_spline
    implicit none
    class(ceq_type), intent(in out) :: self
    call new_spline(self%nr, self%eqpsi, self%diam, self%diam_spl)
    call new_spline(self%nr, self%eqpsi, self%rc, self%rcenter_spl)
    call new_spline(self%nr, self%psi_bar, self%fp, self%btori_spl)
    call new_spline(self%nr, self%psi_bar, self%qsf, self%qsf_spl)
    call new_spline(self%nr, self%psi_bar, self%pressure, self%pressure_spl)
    call new_spline(self%nr, self%psi_bar, self%beta, self%beta_spl)
  end subroutine setup_splines

  subroutine delete_splines(self) !Identical to geq
    use splines, only: delete_spline
    implicit none
    class(ceq_type), intent(in out) :: self
    call delete_spline(self%diam_spl)
    call delete_spline(self%rcenter_spl)
    call delete_spline(self%btori_spl)
    call delete_spline(self%qsf_spl)
    call delete_spline(self%pressure_spl)
    call delete_spline(self%beta_spl)
  end subroutine delete_splines

  !> Calculate the derivative of f w.r.t to the radial
  !> and poloidal indexes (i.e. calculate the finite
  !> differences without dividing by
  !> delta psi and delta theta).
  !> - dfm(:,:,1) is the psi derivative
  !> - dfm(:,:,2) is the theta derivative
  !> - char specifies the periodicity in theta
  !>   - 'E', 'O' mean continuous at theta = +/- pi
  !>   - 'T' means a jump of 2 pi at theta = +/- pi
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(ceq_type), intent(in) :: self
    character(1), intent(in) :: char
    real, intent(in) :: f(:,:)
    real, intent(out) :: dfm(:,:,:)
    integer :: i, j, nr, nt
    nr = self%nr
    nt = self%nt

    i = 1
    ! Radial derivative
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i = nr
    ! Radial derivative
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    select case (char)
    case ('E')
       j = 1
       ! theta derivative, periodic
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))

       j = nt
       dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
    case ('O')
       j = 1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))

       j = nt
       dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
    case ('T')
       j = 1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1)+2.*pi)

       j = nt
       dfm(:,j,2) = 0.5*(2.*pi + f(:,2) - f(:,j-1))
    end select

    do i = 2, nr-1
       dfm(i,:,1) = 0.5*(f(i+1,:)-f(i-1,:))
    end do

    do j = 2, nt-1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j-1))
    end do
  end subroutine derm

  !> Calculate the derivative of a selected field w.r.t. R and Z
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(ceq_type), intent(in) :: self
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array
    integer, intent(in) :: ntm
    !> Select field to take gradient of:
    !>   - 'D': Bt
    !>   - 'P': psi
    !>   - 'R': pressure
    !>   - 'T': theta
    !>   - 'B': for debugging with `bishop==0` only
    character(1), intent(in) :: char
    !> Minor radius grid of flux surface (?)
    real, intent(in) :: rgrid(-ntm:)
    !> Theta grid
    real, intent(in) :: theta(-ntm:)
    !> grad(:,1): derivative w.r.t. R as function of theta
    !> grad(:,2): derivative w.r.t  Z as function of theta
    real, intent(out) :: grad(-ntm:,:)
    !> Radius for taking pressure derivative (FIXME: double-check)
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, self%nt, 2) :: dcart
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

    do i = -nth_used, nth_used
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), grad(i,1), 'R')
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,2), grad(i,2), 'Z')
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

  !> Calculate the derivative of a selected field w.r.t. the Bishop
  !> coordinate system
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use mp, only: mp_abort
    use splines, only: inter_d_cspl
    implicit none
    class(ceq_type), intent(in) :: self
    !> Number of theta points
    integer, intent(in) :: nth_used
    !> Lower index of theta array
    integer, intent(in) :: ntm
    !> Select field to take gradient of:
    !>   - 'D': Bt
    !>   - 'P': psi
    !>   - 'R': pressure
    !>   - 'T': theta
    !>   - 'B': for debugging with `bishop==0` only
    character(1), intent(in) :: char
    !> Minor radius grid of flux surface (?)
    real, intent(in) :: rgrid(-ntm:)
    !> Theta grid
    real, intent(in) :: theta(-ntm:)
    !> Computed derivative, second dimension is the Bishop coordinate
    real, intent(out) :: grad(-ntm:,:)
    !> Radius for taking pressure derivative (?)
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, self%nt, 2) ::  dbish
    integer :: i

    select case(char)
    case('B')
       call mp_abort('error: bishop = 1 not allowed with ceq. (2)',.true.)
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

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nr, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i = -nth_used, nth_used
          grad(i,1) = grad(i,1)*daa(1) * 0.5*self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5*self%beta_0
       end do
    end if
  end subroutine bgradient

  !> Calculates fstar which is f interpolated at the location (r,theta).
  !> Here r is the normalised poloidal flux coordinate rp (= psi_pN + psi_0N)
  !> and theta_in is theta. f is a grid of values of f as a function of psi_p,theta
  subroutine eqitem(self, r_in, theta_in, f, fstar, char)
    use mp, only: mp_abort
    use constants, only: pi
    class(ceq_type), intent(in) :: self
    real, intent(in) :: r_in, theta_in
    real, intent(out) :: fstar
    real, dimension(:,:), intent(in) :: f
    character(len = 1), intent(in) :: char
    integer :: i, j, istar, jstar
    real :: thet, sign, r
    real :: st, dt, sr, dr
    real, dimension(size(f,2)) :: mtheta

    r = r_in

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
    endif

    ! find r on psi mesh

    ! disallow evaluations on or outside the surface for now
    if (r >= self%eqpsi(self%nr)) then
       write(*,*) 'No evaluation of eqitem allowed on or outside surface'
       write(*,*) '(Could this be relaxed a bit?)'
       write(*,*) r, theta_in, self%eqpsi(self%nr), sign
       call mp_abort('No evaluation of eqitem allowed on or outside surface')
    endif

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
    endif

    ! Now do theta direction
    thet = mod2pi(theta_in)

    ! get thet on theta mesh
    do j = 1, self%nt
       mtheta(j) = (j-1) * 2. * pi / real(self%nt-1) - pi
    end do

    ! note that theta(1)=0 for gen_eq theta
    jstar = -1
    do j = 1, self%nt
       if (jstar /= -1) cycle
       if (thet < mtheta(j)) then
          dt = thet - mtheta(j-1)
          st = mtheta(j) - thet
          jstar = j-1
       end if
    end do

    ! treat pi separately
    if (jstar == -1) then
       jstar = self%nt-1
       dt = mtheta(jstar+1) - mtheta(jstar)
       st = 0.
    end if

    ! use opposite area stencil to interpolate
    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar / (abs(self%eqpsi(istar+1) - self%eqpsi(istar)) * &
         (mtheta(jstar+1) - mtheta(jstar)))
  end subroutine eqitem

  !> Converts derivatives w.r.t. (psi_index,theta_index) to derivatives
  !> w.r.t. (R,Z).
  !> Note that `dfcart(1, :, :)` is not valid
  pure subroutine eqdcart(self, dfm, dfcart)
    implicit none
    class(ceq_type), intent(in) :: self
    !> dfm(:,:,1) is deriv w.r.t. psi_index (i.e. (df/dpsi)_theta * delta psi);
    !> dfm(:,:,2) is deriv w.r.t. theta_index
    real, dimension (:,:,:), intent(in)  :: dfm
    !> dfcart(:,:,1) is deriv w.r.t. R;
    !> dfcart(:,:,2) is deriv w.r.t. Z
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (self%nr, self%nt) :: denom
    integer :: i, j
    denom = self%drm(:,:,1) * self%dzm(:,:,2) - self%drm(:,:,2) * self%dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * self%dzm(:,:,2) - self%dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * self%drm(:,:,2) + self%drm(:,:,1) * dfm(:,:,2)

    ! Inner loop skips first point because it's likely exactly on the
    ! magnetic axis, which means the Jacobian is zero. Do we need
    ! need a better way to either calculate dfcart there, or to ensure
    ! we ignore it in other places?
    do j = 1, self%nt
       do i = 2, self%nr
          dfcart(i,j,:) = dfcart(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdcart

  !> Convert gradients of a function f w.r.t. R,Z into bishop form
  !> Note that `dbish(1, :, :)` is not valid
  pure subroutine eqdbish(self, dcart, dbish)
    implicit none
    class(ceq_type), intent(in) :: self
    !> dcart(:,:,1) is gradient of f w.r.t. R;
    !> dcart(:,:,2) is gradient of f w.r.t. Z
    real, dimension(:, :, :), intent (in) :: dcart
    !> dbish(:,:,1) is set to (df/dR dpsi/dR + df/dZ dpsi/dZ)/|grad psi|;
    !> dbish(:,:,2) is set to (-df/dR dpsi/dZ + df/dZ dpsi/dR)/|grad psi|.
    !> Note that in the special case where f=psi:
    !> dbish(:,:,1) is |grad psi|;
    !> dbish(:,:,2) is 0
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j
    ! denom is |grad psi|
    denom = sqrt(self%dpcart(:,:,1)**2 + self%dpcart(:,:,2)**2)
    dbish(:,:,1) = dcart(:,:,1) * self%dpcart(:,:,1) + dcart(:,:,2) * self%dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * self%dpcart(:,:,2) + dcart(:,:,2) * self%dpcart(:,:,1)

    do j = 1, self%nt
       do i = 2, self%nr
          dbish(i,j,:) = dbish(i,j,:)/ denom(i,j)
       end do
    end do
  end subroutine eqdbish

  !> Return the inverse of the major radius at \((r, \theta)\)
  real function invR (self, r, theta)
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%Rpos(r, theta)
  end function invR

  !> Return the major radius at \((r, \theta)\)
  real function Rpos (self, r, theta)
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%R_psi, Rpos, 'R')
  end function Rpos

  !> Return the vertical coordinate at \((r, \theta)\)
  real function Zpos (self, r, theta)
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%Z_psi, Zpos, 'Z')
  end function Zpos

  !> Return the flux surface label `psi`
  real function psi (self, r, theta)
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent (in) :: r, theta
    psi = r
  end function psi

  !> Return the diameter of the flux surface at a given major radius
  real function diameter (self, rp)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = splint(rp, self%diam_spl)
  end function diameter

  !> Return the major radius of the centre of a given flux surface
  real function rcenter (self, rp)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: rp
    rcenter = splint(rp, self%rcenter_spl)
  end function rcenter

  !> \(\frac{dI}{d\psi}\) at the given normalised \(\psi\)
  real function dbtori (self, pbar)
    use splines, only: dsplint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    dbtori = dsplint(pbar, self%btori_spl)/(self%psi_a - self%psi_0)
  end function dbtori

  !> FIXME : Add documentation
  real function btori (self, pbar)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    btori = splint(pbar, self%btori_spl)
  end function btori

  !> FIXME : Add documentation
  real function qfun (self, pbar)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    qfun = splint(pbar, self%qsf_spl)
  end function qfun

  !> FIXME : Add documentation
  real function pfun (self, pbar)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    pfun = 0.5*self%beta_0*splint(pbar, self%pressure_spl)
  end function pfun

  !> FIXME : Add documentation
  real function dpfun (self, pbar)
    use splines, only: dsplint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    dpfun = dsplint(pbar, self%pressure_spl) / (self%psi_a-self%psi_0) * self%beta_0/2.
  end function dpfun

  !> FIXME : Add documentation
  real function betafun (self, pbar)
    use splines, only: splint
    implicit none
    class(ceq_type), intent(in) :: self
    real, intent(in) :: pbar
    if (pbar == 0.) then
       betafun = self%beta(1)
       return
    end if
    betafun = splint(pbar, self%beta_spl)
  end function betafun

  !> FIXME : Add documentation
  subroutine Hahm_Burrell(geom, irho, a)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    implicit none
    type(ceq_type), intent(in) :: geom
    real, intent(in) :: a
    integer, intent(in) :: irho
    integer :: i, nr
    real :: gradpsi, mag_B, rho_eq, rp1, rp2, rho1, rho2, drhodpsiq
    real, dimension(geom%nr) :: gamma, pbar, dp, d2p, pres
    integer :: fort24_unit
    if (.not. proc0) return
    nr = geom%nr

    gamma = 0.

    pbar = (geom%eqpsi-geom%eqpsi(1))/(geom%eqpsi(nr)-geom%eqpsi(1))

    do i=2, nr-1
       dp(i)  = geom%dpfun(pbar(i))
    enddo

    do i=3,nr-2
       d2p(i) = (dp(i+1)-dp(i-1))/(geom%eqpsi(i+1)-geom%eqpsi(i-1))
    enddo

    pres=0.
    do i=3, nr-2
       rp1=geom%eqpsi(i+1)
       rp2=geom%eqpsi(i-1)
       rho1=0.5*geom%diameter(rp1)
       rho2=0.5*geom%diameter(rp2)
       drhodpsiq=(rho1-rho2)/(rp1-rp2)

       call geom%eqitem(geom%eqpsi(i), 0., geom%dpbish(:,:,1), gradpsi, 'R')

       mag_B = sqrt( (geom%btori(pbar(i))**2 + gradpsi**2))/geom%Rpos(geom%eqpsi(i),0.)

       pres(i) = geom%pfun(pbar(i))

       gamma(i) = 0.01*gradpsi**2*(d2p(i)/pres(i)-a*(dp(i)/pres(i))**2) &
            /mag_B*(2.*pres(i)/geom%beta_0)**((1-a)/2.) &
            *(-pres(i)/(dp(i)/drhodpsiq))
    enddo

    call open_output_file(fort24_unit,".fort.24")
    do i=3,nr-2
       if(irho == 1) then
          rho_eq = geom%eqpsi(i)
       else if(irho == 2) then
          rho_eq = 0.5 * geom%diameter(geom%eqpsi(i))
       else if(irho == 3) then
          rho_eq = pbar(i)
       endif
       write(fort24_unit,1000) i, pbar(i), rho_eq, pres(i), gamma(i)
    enddo

    if(irho == 1) write(fort24_unit,*) '# irho = 1 produces psi instead of rho_eq'

    call close_output_file(fort24_unit)

1000 format(i5,11(1x,e16.9))
  end subroutine Hahm_Burrell
end module ceq
