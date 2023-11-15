!> FIXME : Add documentation
module geq
# ifdef NETCDF
  use netcdf, only: NF90_FLOAT, NF90_DOUBLE
  use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
  use netcdf, only: nf90_create, nf90_open, nf90_sync, nf90_close
  use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
  use netcdf, only: nf90_put_var, nf90_get_var, nf90_strerror
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid, nf90_inquire_variable
  use netcdf_utils, only: netcdf_error
# endif
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi, EQFILE_LENGTH
  use splines, only: spline
  implicit none

  private

  public :: Hahm_Burrell, geq_type

  type, extends(abstract_eqfile_geo_type) :: geq_type
     private
     real, allocatable, dimension (:,:)   :: B_psi
     real, allocatable, dimension (:,:,:) :: dbm, dbcart, dbbish
     real, allocatable, dimension (:)     :: diam, rc
     type(spline) :: diam_spl, rcenter_spl
   contains
     procedure :: initialise => eqin
     procedure :: calculate_gradients => geq_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: rcenter => rcenter
     procedure :: diameter => diameter
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     procedure :: alloc_arrays => alloc_arrays
     procedure :: dealloc_arrays => dealloc_arrays
  end type geq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file
  !! containing the axisymmetric magnetic field geometry in flux
  !! coordinates
  subroutine eqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use geo_utils, only: geo_input_type
# ifndef NETCDF
    use mp, only: mp_abort
# endif
    implicit none
    class(geq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(out) :: nthg
    real :: I_0
    real, dimension(:), allocatable :: rho_d
    integer :: istatus
    integer :: ncid, id
    ! what is the best way to handle the netcdf single/double problem?
    real :: f_N, psi_N

    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = self%nt
       return
    end if

    !     read the data
    self%filename = trim(adjustl(inputs%eqfile))

# ifdef NETCDF
    !     netcdf open file
    istatus = nf90_open(self%filename, NF90_NOWRITE, ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file = self%filename)

    !     netcdf read scalar: nr
    !
    !     nr == number of radial grid points in radial eq grid
    istatus = nf90_inq_dimid (ncid, 'psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='psi')

    istatus = nf90_inquire_dimension (ncid, id, len = self%nr)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    !     netcdf read scalar: nt
    !
    !     nt == number of theta grid points in theta eq grid

    istatus = nf90_inq_dimid (ncid, 'theta', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='theta')
    istatus = nf90_inquire_dimension (ncid, id, len = self%nt)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    call self%alloc_arrays()
    allocate(rho_d(self%nr))

    !     netcdf read scalars: psi_0,psi_a,B_T,I_0
    !
    !     psi_0 == value of psi at the magnetic axis
    !     psi_a == value of psi at the boundary of the plasma
    !     B_T == vacuum toroidal magnetic field at R_center
    !     I_0 == total plasma current

    istatus = nf90_inq_varid (ncid, 'psi_0', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi_0')
    istatus = nf90_get_var (ncid, id, self%psi_0)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%psi_0 = self%psi_0 * 1.e-8

    istatus = nf90_inq_varid (ncid, 'psi_a', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi_a')
    istatus = nf90_get_var (ncid, id, self%psi_a)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%psi_a = self%psi_a * 1.e-8

    istatus = nf90_inq_varid (ncid, 'B_T', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='B_T')
    istatus = nf90_get_var (ncid, id, self%B_T)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%B_T = self%B_T * 1.e-4  ! cgs to MKS

    istatus = nf90_inq_varid (ncid, 'I_0', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='I_0')
    istatus = nf90_get_var (ncid, id, I_0)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     netcdf read scalars: R_mag,Z_mag,aminor
    !
    !     R_mag == R position of magnetic axis
    !     Z_mag == Z position of magnetic axis
    !     aminor    == half diameter of last flux surface at elevation of Z_mag

    istatus = nf90_inq_varid (ncid, 'Rmag', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Rmag')
    istatus = nf90_get_var (ncid, id, self%R_mag)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%R_mag = self%R_mag / 100.

    istatus = nf90_inq_varid (ncid, 'Zmag', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Zmag')
    istatus = nf90_get_var (ncid, id, self%Z_mag)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%Z_mag = self%Z_mag / 100.

    istatus = nf90_inq_varid (ncid, 'a', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a')
    istatus = nf90_get_var (ncid, id, self%aminor)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%aminor = self%aminor / 100.

    !     netcdf read vectors: rho_d,eqpsi,psi_bar,fp,qsf,beta,pressure
    !
    !     rho_d(1:nr) == half diameters of flux surfaces at elevation
    !                             of Z_mag on the radial grid
    !     psi is the poloidal flux
    !     eqpsi(1:nr) == values of psi on the radial grid
    !     psi_bar(1:nr) == values of psi_bar on the radial grid
    !     [psi_bar == (eqpsi - psi_0)/(psi_a - psi_0) if not available]
    !     fp(1:nr) == the function that satisfies
    !              B = fp grad zeta + grad zeta x grad psi
    !     qsf(1:nr) == q profile on the radial grid
    !     beta(1:nr) == local beta, with the magnetic field defined
    !     to be vacuum magnetic field on axis
    !     pressure(1:nr) == pressure profile on the radial grid,
    !     normalized to the value at the magnetic axis.

    istatus = nf90_inq_varid (ncid, 'rho', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='rho')
    istatus = nf90_get_var (ncid, id, rho_d)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi')
    istatus = nf90_get_var (ncid, id, self%eqpsi)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%eqpsi = self%eqpsi * 1.e-8

    istatus = nf90_inq_varid (ncid, 'psibar', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psibar')
    istatus = nf90_get_var (ncid, id, self%psi_bar)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'fp', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='fp')
    istatus = nf90_get_var (ncid, id, self%fp)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%fp = self%fp * 1.e-6

    istatus = nf90_inq_varid (ncid, 'q', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='q')
    istatus = nf90_get_var (ncid, id, self%qsf)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'beta', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='beta')
    istatus = nf90_get_var (ncid, id, self%beta)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    self%beta_0 = self%beta(1)

    istatus = nf90_inq_varid (ncid, 'pressure', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='pressure')
    istatus = nf90_get_var (ncid, id, self%pressure)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     netcdf read 2d field: R_psi,Z_psi and B_psi (mod(B))
    !     eq_Rpsi(1:nr, 1:nt)
    !     eq_Zpsi(1:nr, 1:nt)
    !     eq_Bpsi(1:nr, 1:nt)
    !

    istatus = nf90_inq_varid (ncid, 'R_psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='R_psi')
    istatus = nf90_get_var (ncid, id, self%R_psi, count = [self%nr, self%nt])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%R_psi = self%R_psi / (100. * self%aminor)

    istatus = nf90_inq_varid (ncid, 'Z_psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Z_psi')
    istatus = nf90_get_var (ncid, id, self%Z_psi, count = [self%nr, self%nt])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%Z_psi = self%Z_psi / (100. * self%aminor)

    istatus = nf90_inq_varid (ncid, 'B_psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='B_psi')
    istatus = nf90_get_var (ncid, id, self%B_psi, count = [self%nr, self%nt])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%B_psi = self%B_psi * 1.e-4 / self%B_T

    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)

    !
    !     Normalize, rename quantities
    !

    avgrmid = self%aminor
    rmaj = self%R_mag / self%aminor   ! used to reference the grid
    B_T0 = self%B_T

    psi_N = B_T0 * avgrmid**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    psi_a_out = self%psi_a
    psi_0_out = self%psi_0
    self%eqpsi = self%eqpsi / psi_N

    f_N = B_T0 * avgrmid
    self%fp = self%fp / f_N

    self%diam = self%R_psi(:, 1) - self%R_psi(:, self%nt)
    self%rc = 0.5 * (self%R_psi(:, 1) + self%R_psi(:, self%nt))

    nthg = self%nt
# else
    call mp_abort('error: geq eqin is called without netcdf',.true.)
#endif

    self%initialised = .true.
  end subroutine eqin

  !> FIXME : Add documentation
  pure subroutine alloc_arrays(self)
    implicit none
    class(geq_type), intent (in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%eqpsi(nr), self%psi_bar(nr), self%fp(nr), self%qsf(nr))
    allocate(self%beta(nr), self%pressure(nr), self%rc(nr), self%diam(nr))
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt), self%B_psi(nr, nt))
    allocate(self%drm(nr, nt, 2), self%dzm(nr, nt, 2), self%dbm(nr, nt, 2))
    allocate(self%dpm(nr, nt, 2), self%dtm(nr, nt, 2))
    allocate(self%dpcart(nr, nt, 2), self%dbcart(nr, nt, 2), self%dtcart(nr, nt, 2))
    allocate(self%dpbish(nr, nt, 2), self%dbbish(nr, nt, 2), self%dtbish(nr, nt, 2))
  end subroutine alloc_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_arrays(self)
    implicit none
    class(geq_type), intent(in out) :: self
    deallocate(self%eqpsi, self%psi_bar, self%fp, &
         self%qsf, self%beta, self%pressure, self%rc, self%diam)
    deallocate(self%R_psi, self%Z_psi, self%B_psi)
    deallocate(self%drm, self%dzm, self%dbm, self%dpm, self%dtm)
    deallocate(self%dpcart, self%dbcart, self%dtcart, self%dtbish)
    deallocate(self%dpbish, self%dbbish)
  end subroutine dealloc_arrays

  !> FIXME : Add documentation
  subroutine geq_init(self)
    use constants, only: pi
    implicit none
    class(geq_type), intent(in out) :: self
    real, dimension(self%nr, self%nt) :: eqpsi1, eqth
    integer j

    do j = 1, self%nt
       eqpsi1(:, j) = self%eqpsi
       eqth(:, j) = (j-1) * pi / real(self%nt-1)
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

    ! grad(theta) in cartesian form
    call self%eqdcart(self%dtm, self%dtcart)
    ! grad(theta) in Bishop form
    call self%eqdbish(self%dtcart, self%dtbish)
  end subroutine geq_init

  subroutine setup_splines(self)
    use splines, only: new_spline
    implicit none
    class(geq_type), intent(in out) :: self
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
    class(geq_type), intent(in out) :: self
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
    class(geq_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nr, nt
    nr = self%nr ; nt = self%nt

    i = 1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i = nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

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

       j = nt
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j = 1
       dfm(:,j,2) = f(:,j+1)

       j = nt
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i = 2, nr - 1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    end do

    do j = 2, nt - 1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j-1))
    end do
  end subroutine derm

  !> FIXME : Add documentation
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(geq_type), intent(in) :: self
    integer, intent (in) :: nth_used, ntm
    character(1), intent (in) :: char
    real, dimension (-ntm:), intent (in)  :: rgrid, theta
    real, dimension (-ntm:,:), intent (out) :: grad
    real, intent(in) :: rp
    real, dimension(self%nr, self%nt, 2) :: dcart
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
          grad(i,1)=grad(i,1)*daa(1) * 0.5 * self%beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5 * self%beta_0
       end do
    end if
  end subroutine gradient

  !> FIXME : Add documentation
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(geq_type), intent(in) :: self
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
    case('P')
       dbish = self%dpbish
    case('R')
       dbish = self%dpbish  ! dpcart is correct for 'R'
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
    use constants, only: pi
    use mp, only: mp_abort
    implicit none
    class(geq_type), intent(in) :: self
    real, intent(in) :: r_in
    real, intent(in) :: theta_in
    real, dimension(:,:), intent(in) :: f
    character(len = 1), intent(in) :: char
    integer :: i, j, istar, jstar
    real :: thet, sign, tp, tps, r
    real :: st, dt, sr, dr, rt
    real, dimension(size(f,2)) :: mtheta

    !Copy input r to avoid any modification
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

    if (r < sign*self%eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r, sign, self%eqpsi(1)
       call mp_abort('r < Psi_0 in eqitem')
    end if

    ! find r on psi mesh

    ! disallow evaluations outside the plasma surface for now
    if (r == self%eqpsi(self%nr)) then
       rt = 0.9999999999*r
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
       if (r < sign*self%eqpsi(i)) then
          dr = r - sign*self%eqpsi(i-1)
          sr = sign*self%eqpsi(i) - r
          istar = i-1
          exit
       end if
    end do

    ! Now do theta direction
    thet = mod2pi(theta_in)

    ! assume up-down symmetry
    tp = abs(thet)
    tps = 1.
    if (char == 'Z' .and. tp >= 1.e-10) tps = thet / tp

    ! get thet on theta mesh
    do j = 1, self%nt
       mtheta(j) = (j-1) * pi / real(self%nt-1)
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
    fstar = fstar * tps / (abs(self%eqpsi(istar+1) - self%eqpsi(istar)) * &
         (mtheta(jstar+1) - mtheta(jstar)))
  end function eqitem

  !> FIXME : Add documentation
  real function diameter (self, rp)
    use splines, only: splint
    implicit none
    class(geq_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = splint(rp, self%diam_spl)
  end function diameter

  !> FIXME : Add documentation
  real function rcenter (self, rp)
    use splines, only: splint
    implicit none
    class(geq_type), intent(in) :: self
    real, intent(in) :: rp
    rcenter = splint(rp, self%rcenter_spl)
  end function rcenter

  !> FIXME : Add documentation
  subroutine Hahm_Burrell(geom, irho, a)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    implicit none
    type(geq_type), intent(in) :: geom
    real, intent(in) :: a
    integer, intent(in) :: irho
    integer :: i, nr
    real :: gradpsi, mag_B, rho_eq, rp1, rp2, rho1, rho2, drhodpsiq
    real, dimension(geom%nr) :: gamma, pbar, dp, d2p, pres, bs_coll, s__hat, bs_sh
    integer :: fort24_unit
    if (.not. proc0) return
    nr = geom%nr
    if(irho /= 2) then
       write(*,*) 'use irho=2 to get correct shearing rate.'
    endif

    gamma = 0.

    pbar = (geom%eqpsi-geom%eqpsi(1))/(geom%eqpsi(nr)-geom%eqpsi(1))

    do i=2, nr-1
       dp(i)  = geom%dpfun(pbar(i))
    enddo

    do i=3,nr-2
       d2p(i) = (dp(i+1)-dp(i-1))/(geom%eqpsi(i+1)-geom%eqpsi(i-1))
    enddo

    pres=0.

    do i=3,nr-2
       rp1=geom%eqpsi(i+1)
       rp2=geom%eqpsi(i-1)
       rho1=0.5*geom%diameter(rp1)
       rho2=0.5*geom%diameter(rp2)
       drhodpsiq=(rho1-rho2)/(rp1-rp2)

       s__hat(i) = (geom%qfun(pbar(i+1))-geom%qfun(pbar(i-1)))/(rho2-rho1) &
            * 0.5*(rho1+rho2)/geom%qfun(pbar(i))

       pres(i) = geom%pfun(pbar(i))

       gradpsi = geom%eqitem(geom%eqpsi(i), 0., geom%dpbish(:,:,1), 'R')
       mag_B = geom%eqitem(geom%eqpsi(i), 0., geom%B_psi, 'R')

       gamma(i) = (d2p(i)/pres(i)-a*(dp(i)/pres(i))**2)
       gamma(i) = 0.01*gradpsi**2*gamma(i) &
            /mag_B*(2.*pres(i)/geom%beta_0)**((1-a)/2.) &
            *(-pres(i)/(dp(i)/drhodpsiq))
    enddo
    !
    ! Assume rho_*0 = 0.01, nu_*0 = 1.0
    !
    do i=3,nr-2
       bs_coll(i)=0.01*((2.*pres(i)/geom%beta_0)**(1-a))**2.5 &
            / (2.*pres(i)/geom%beta_0)**a * geom%qfun(0.) &
            * 0.5*geom%diameter(geom%eqpsi(i)) &
            * (-dp(i)/pres(i)) * geom%R_psi(1,1)**1.5
       bs_sh(i)=(-pres(i)/(dp(i)/drhodpsiq))*s__hat(i)/geom%R_psi(1,1) &
            / geom%qfun(pbar(i)) /sqrt(pres(i))
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
       write(fort24_unit,1000) i, pbar(i), rho_eq, pres(i), gamma(i), bs_coll(i), bs_sh(i)
    enddo

    if(irho == 1) write(fort24_unit,*) '# irho = 1 produces psi instead of rho_eq'

    call close_output_file(fort24_unit)

1000 format(i5,11(1x,e16.9))
  end subroutine Hahm_Burrell
end module geq
