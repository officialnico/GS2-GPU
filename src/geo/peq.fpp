!> FIXME : Add documentation
module peq
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

  public :: Hahm_Burrell, peq_type, teq_type

  type, abstract, extends(abstract_eqfile_geo_type) :: peq_base_type
     private
     logical :: transp
     real, allocatable, dimension(:) :: rho_b
     character(len = EQFILE_LENGTH) :: filename
     integer :: nr, nt
     real :: psi_0, psi_a, B_T, beta_0
     real :: R_mag, Z_mag, aminor
     real, allocatable, dimension (:)     :: rho_d, eqpsi, psi_bar, fp, qsf
     real, allocatable, dimension (:)     :: beta, pressure, diam, rc
     real, allocatable, dimension (:,:)   :: R_psi, Z_psi
     real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
     real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
     real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish
     type(spline) :: diam_spl, rcenter_spl, btori_spl, qsf_spl, pressure_spl, beta_spl
   contains
     procedure :: calculate_gradients => peq_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: finalise => peq_finish
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
  end type peq_base_type

  type, extends(peq_base_type) :: peq_type
   contains
     procedure :: initialise_eqfile => eqin
  end type peq_type

  type, extends(peq_base_type) :: teq_type
   contains
     procedure :: initialise_eqfile => teqin
  end type teq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file
  !> containing the axisymmetric magnetic field geometry in flux
  !> coordinates
  subroutine eqin(self, eqfile, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, initeq, nthg, &
       big_in)
    use mp, only: mp_abort
    use constants, only: pi
    implicit none
    class(peq_type), intent(in out) :: self
    character (len = EQFILE_LENGTH), intent(in) :: eqfile
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(in) :: initeq
    integer, intent(out) :: nthg
    integer, intent(in), optional :: big_in
    real :: d, R_geo
    integer :: istatus
    integer :: ncid, id, i, j
    integer :: nchar, i_sym
    integer nz1, nz2
    real, allocatable, dimension(:) :: work
    real :: f_N, psi_N
    !     read the data

# ifdef NETCDF
    ! WARNING: Exiting early without setting `nthg` may result in `nthg`
    ! being uninitialised in the calling procedure. The value it gets set to
    ! must be consistent with the value set at the end of this subroutine
    ! when not exiting early or the calling procedure may end up with
    ! inconsistently sized arrays. This is expected to affect `ideal_ball`
    ! (unconfirmed) but not `gs2` (also unconfirmed).
    if(initeq == 0) return

    self%filename = trim(adjustl(eqfile))

    !     netcdf open file
    istatus = nf90_open(self%filename, NF90_NOWRITE, ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file = self%filename)

    !     netcdf read scalar: nr
    !
    !     nz2 == number of points in radial array
    !     nr == number of actual grid points in radial array
    istatus = nf90_inq_dimid (ncid, 'z2', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='z2')
    istatus = nf90_inquire_dimension (ncid, id, len=nz2)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    istatus = nf90_inq_dimid (ncid, 'z1', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='z1')
    istatus = nf90_inquire_dimension (ncid, id, len=nz1)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    istatus = nf90_inq_dimid (ncid, 'npsi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='npsi')
    istatus = nf90_inquire_dimension (ncid, id, len = self%nr)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    istatus = nf90_inq_varid (ncid, 'nxy', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='nxy')
    istatus = nf90_get_var (ncid, id, i_sym, start=[8])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     netcdf read scalar: nt
    !
    !     nt == number of theta grid points in theta eq grid
    istatus = nf90_inq_dimid (ncid, 'nthe', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='nthe')
    istatus = nf90_inquire_dimension (ncid, id, len = self%nt)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

    call self%alloc_arrays()

    !     netcdf read vectors: rho_d, eqpsi, psi_bar, fp, beta, pressure
    !
    !     rho_d(1:nr) == half diameters of flux surfaces at elevation
    !                             of Z_mag on the radial grid
    !     psi is the poloidal flux
    !     eqpsi(1:nr) == values of psi on the radial grid
    !     psi_bar(1:nr) == values of psi_bar on the radial grid
    !     [psi_bar == (eqpsi - psi_0)/(psi_a - psi_0) if not available]
    !     fp(1:nr) == the function that satisfies
    !              B = fp grad zeta + grad zeta x grad psi
    !     beta(1:nr) == local beta, with the magnetic field defined
    !     to be vacuum magnetic field on axis
    !     pressure(1:nr) == pressure profile on the radial grid,
    !     normalized to the value at the magnetic axis.

    istatus = nf90_inq_varid (ncid, 'psivec', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psivec')
    istatus = nf90_get_var (ncid, id, self%eqpsi, count = [self%nr])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'fvec', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='fvec')
    istatus = nf90_get_var (ncid, id, self%fp, count = [self%nr])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'pvec', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='pvec')
    istatus = nf90_get_var (ncid, id, self%eqpsi, count = [self%nr])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     assign scalars: psi_0, psi_a, B_T
    !
    !     psi_0 == value of psi at the magnetic axis
    !     psi_a == value of psi at the boundary of the plasma
    !     B_T == vacuum toroidal magnetic field at R_center
    !
    self%psi_0 = self%eqpsi(1)
    self%psi_a = self%eqpsi(self%nr)

    self%psi_bar = (self%eqpsi - self%psi_0) / (self%psi_a - self%psi_0)

    !     netcdf read 2d field: R_psi,Z_psi and B_psi (mod(B))
    !     R_psi(1:nr, 1:nt) -- transpose of Menard's storage ...
    !     Z_psi(1:nr, 1:nt)
    !     B_psi(1:nr, 1:nt)
    allocate(work(nz1*nz2))

    istatus = nf90_inq_varid (ncid, 'x', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='x')
    istatus = nf90_get_var (ncid, id, work, count = [nz1*nz2])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    if (i_sym == 0) then
       write(*,*) 'Up-down asymmetries not yet handled correctly.'
    end if

    do j = 1, self%nt
       do i = 1, self%nr
          self%R_psi(i,j) = work(3+j-1+nz1*(i-1))
       end do
    end do

    istatus = nf90_inq_varid (ncid, 'z', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='z')
    istatus = nf90_get_var (ncid, id, work, count=[nz1, nz2]) !Should this be nz1*nz2?
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    if (i_sym == 0) then
       write(*,*) 'Up-down asymmetries not yet handled correctly.'
    end if

    do j = 1, self%nt
       do i = 1, self%nr
          self%Z_psi(i,j) = work(3+j-1+nz1*(i-1))
       end do
    end do


    !     assign scalars: R_mag, Z_mag, aminor
    !
    !     R_mag == R position of magnetic axis
    !     Z_mag == Z position of magnetic axis
    !     aminor    == half diameter of last flux surface at elevation of Z_mag

    d = self%R_psi(self%nr, 1) - self%R_psi(self%nr, self%nt)
    R_geo = (self%R_psi(self%nr, 1) + self%R_psi(self%nr, self%nt))/2.
    self%aminor = d/2.

    self%B_T = self%fp(self%nr) / R_geo

    self%R_psi = self%R_psi / self%aminor
    self%Z_psi = self%Z_psi / self%aminor

    self%R_mag = self%R_psi(1,1)
    self%Z_mag = self%Z_psi(1,1)

    self%beta = 8. * pi * 1.e-7 * self%pressure / self%B_T**2  ! not correct definition yet
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    istatus = nf90_close(ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)

    deallocate(work)

    !
    !     Normalize, rename quantities
    !

    avgrmid = self%aminor
    rmaj = self%R_mag
    B_T0 = abs(self%B_T)

    psi_N = B_T0 * avgrmid**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    psi_a_out = self%psi_a
    psi_0_out = self%psi_0
    self%eqpsi = self%eqpsi / psi_N

    f_N = B_T0 * avgrmid
    self%fp = self%fp / f_N

    nthg = self%nt
    self%diam = abs(self%R_psi(:, self%nt/2) - self%R_psi(:, 1))
    self%rc = 0.5 * (self%R_psi(:, 1) + self%R_psi(:, self%nt/2))

    self%transp = .false.
    self%initialised = .true.
# else
    call mp_abort('error: peq eqin is called without netcdf',.true.)
# endif
  end subroutine eqin

  !> This subroutine reads a generic NetCDF equilibrium file
  !> containing the axisymmetric magnetic field geometry in flux
  !> coordinates
  subroutine teqin(self, eqfile, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, initeq, nthg, &
       big_in)
    use mp, only: mp_abort
    use constants, only: pi
    implicit none
    class(teq_type), intent(in out) :: self
    character (len = EQFILE_LENGTH), intent(in) :: eqfile
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(in) :: initeq
    integer, intent(out) :: nthg
    integer, intent(in), optional :: big_in
    real :: R_geo
    integer :: istatus
    integer :: ncid, id, i, j
    integer :: nchar
    ! what is the best way to handle the netcdf single/double problem?
    real, allocatable, dimension(:,:) :: work
    real :: f_N, psi_N
    !     read the data

    if (initeq == 0) then
       nthg = self%nt/2
       return
    end if

    self%filename = trim(adjustl(eqfile))

# ifdef NETCDF
    !     netcdf open file
    istatus = nf90_open(self%filename, NF90_NOWRITE, ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file = self%filename)

    !     netcdf read scalar: nr
    !
    !     nr == number of radial grid points in radial eq grid

    istatus = nf90_inq_varid (ncid, 'ns', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='ns')
    istatus = nf90_get_var (ncid, id, self%nr)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     netcdf read scalar: nt
    !
    !     nt == number of theta grid points in theta eq grid
    istatus = nf90_inq_varid (ncid, 'nt1', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='nt1')
    istatus = nf90_get_var (ncid, id, self%nt)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    call self%alloc_arrays()

    !     netcdf read scalars: psi_0,psi_a,B_T
    !
    !     psi_0 == value of psi at the magnetic axis
    !     psi_a == value of psi at the boundary of the plasma
    !     B_T == vacuum toroidal magnetic field at R_center

    !       id = ncvid (ncid, 'B_T', ifail)
    !       call ncvgt1 (ncid, id, start, work1, ifail)
    !       B_T = work1*1.e-4  ! cgs to MKS

    !     netcdf read vectors: eqpsi,fp,qsf,beta,pressure
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
    istatus = nf90_get_var (ncid, id, self%rho_b)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi')
    istatus = nf90_get_var (ncid, id, self%eqpsi)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%eqpsi = abs(self%eqpsi)

    self%psi_0 = self%eqpsi(1)
    self%psi_a = self%eqpsi(self%nr)

    self%psi_bar = (self%eqpsi - self%psi_0) / (self%psi_a - self%psi_0)

    istatus = nf90_inq_varid (ncid, 'g', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='g')
    istatus = nf90_get_var (ncid, id, self%fp)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%fp = abs(self%fp)

    istatus = nf90_inq_varid (ncid, 'q', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='q')
    istatus = nf90_get_var (ncid, id, self%qsf)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
    self%qsf = abs(self%qsf)

    istatus = nf90_inq_varid (ncid, 'p', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='p')
    istatus = nf90_get_var (ncid, id, self%pressure)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     netcdf read 2d field: R_psi,Z_psi and B_psi (mod(B))
    !     eq_Rpsi(1:nr, 1:nt)
    !     eq_Zpsi(1:nr, 1:nt)
    !     eq_Bpsi(1:nr, 1:nt)

    allocate(work(self%nt, self%nr))

    istatus = nf90_inq_varid (ncid, 'R', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='R')
    istatus = nf90_get_var (ncid, id, work)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    do j = 1, self%nt
       do i = 1, self%nr
          self%R_psi(i, j) = work(j, i)
       end do
    end do

    istatus = nf90_inq_varid (ncid, 'Z', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Z')
    istatus = nf90_get_var (ncid, id, work)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    do j = 1, self%nt
       do i = 1, self%nr
          self%Z_psi(i, j) = work(j, i)
       end do
    end do

    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)

    deallocate(work)

    ! mildly gimpy at the moment, because the TRANSP output does not
    ! regularly include the point at theta=0.

    ! find aminor
    self%aminor = 0.5 * abs(self%R_psi(self%nr, self%nt/2) - self%R_psi(self%nr, 1))

    ! find R_geo
    R_geo = 0.5 * (self%R_psi(self%nr, self%nt/2) + self%R_psi(self%nr, 1))

    ! find rho_d
    self%rho_d = 0.5 * (self%R_psi(:, self%nt/2) - self%R_psi(:, 1)) / self%aminor

    ! find magnetic axis
    self%R_mag = self%R_psi(1, 1)
    self%Z_mag = self%Z_psi(1, 1)

    ! Find normalizing magnetic field
    self%B_T = self%fp(self%nr) / R_geo

    !
    !     Normalize, rename quantities
    !
    R_geo = R_geo / self%aminor
    self%R_psi = self%R_psi / self%aminor
    self%Z_psi = self%Z_psi / self%aminor

    self%beta = 8. * pi * 1.e-7 * self%pressure / self%B_T**2  ! not correct definition yet?
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    avgrmid = self%aminor
    rmaj = self%R_mag / self%aminor   ! used to reference the grid
    B_T0 = abs(self%B_T)

    psi_N = B_T0 * avgrmid**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    psi_a_out = self%psi_a
    psi_0_out = self%psi_0
    self%eqpsi = self%eqpsi / psi_N

    f_N = B_T0*avgrmid
    self%fp = self%fp / f_N

    !CMR 28/6/2011: why divide the number of theta points by two here?
    nthg = self%nt/2
    self%diam = abs(self%R_psi(:, self%nt/2) - self%R_psi(:, 1))
    self%rc = 0.5 * (self%R_psi(:, 1) + self%R_psi(:, self%nt/2))

    self%transp = .true.
    self%initialised = .true.
# else
    call mp_abort('error: peq teqin is called without netcdf',.true.)
# endif
  end subroutine teqin

  !> FIXME : Add documentation
  pure subroutine alloc_arrays(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt

    allocate(self%rho_d(nr), self%eqpsi(nr), self%psi_bar(nr), self%fp(nr), self%beta(nr), &
         self%pressure(nr), self%rc(nr), self%diam(nr), self%qsf(nr), self%rho_b(nr))
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt))
    allocate(self%drm(nr, nt, 2), self%dzm(nr, nt, 2), self%dbtm(nr, nt, 2), &
         self%dpm(nr, nt, 2), self%dtm(nr, nt, 2))
    allocate(self%dpcart(nr, nt, 2), self%dtcart(nr, nt, 2), self%dbtcart(nr, nt, 2))
    allocate(self%dpbish(nr, nt, 2), self%dtbish(nr, nt, 2), self%dbtbish(nr, nt, 2))
  end subroutine alloc_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_arrays(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    if (allocated(self%rho_d)) deallocate(self%rho_d, self%eqpsi, self%psi_bar, self%fp, &
         self%beta, self%pressure, self%rc, self%diam, self%qsf, self%rho_b)
    if (allocated(self%R_psi)) deallocate(self%R_psi, self%Z_psi)
    if (allocated(self%drm)) deallocate(self%drm, self%dzm, self%dbtm, self%dpm, self%dtm)
    if (allocated(self%dpcart)) deallocate(self%dpcart, self%dtcart, self%dbtcart)
    if (allocated(self%dpbish)) deallocate(self%dpbish, self%dtbish, self%dbtbish)
  end subroutine dealloc_arrays

  !> FIXME : Add documentation
  subroutine peq_finish(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    self%initialised = .false.
    call self%dealloc_arrays
    call self%delete_splines
  end subroutine peq_finish

  !> FIXME : Add documentation
  subroutine peq_init(self)
    use constants, only: pi
    implicit none
    class(peq_base_type), intent(in out) :: self
    real, dimension(self%nr, self%nt) :: eqpsi1, eqth, eqbtor
    integer :: i, j

    do j = 1, self%nt
       eqbtor(:, j) = self%fp / self%R_psi(:, j)
       eqpsi1(:, j) = self%eqpsi
    end do

    if (self%transp) then
       do j = 1, self%nt
          eqth(:, j) = (j-1) * 2.*pi / float(self%nt-1) - pi
       end do
    else
       do j = 1, self%nt
          eqth(:, j) = (j-1) *pi / float(self%nt-1)
       end do
    end if

    call self%derm(eqth, self%dtm, 'T')
    call self%derm(self%R_psi, self%drm, 'E')
    call self%derm(self%Z_psi, self%dzm, 'O')
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
  end subroutine peq_init

  subroutine setup_splines(self)
    use splines, only: new_spline
    implicit none
    class(peq_base_type), intent(in out) :: self
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
    class(peq_base_type), intent(in out) :: self
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
    class(peq_base_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nr, nt
    nr = self%nr
    nt = self%nt

    i = 1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i = nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    if (.not. self%transp) then
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

    else

       select case (char)
       case ('E')
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))

          j=nt
          dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
       case ('O')
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))

          j=nt
          dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
       case ('T')
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1)+2.*pi)

          j=nt
          dfm(:,j,2) = 0.5*(2.*pi + f(:,2) - f(:,j-1))
       end select

    end if

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
    class(peq_base_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real tmp(2), aa(1), daa(1), rpt(1)
    real, dimension(self%nr, self%nt, 2) :: dcart
    integer :: i

    select case(char)
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
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       if (char == 'T' .and. .not. self%transp) then
          grad(i,1) = -tmp(1)
          grad(i,2) = -tmp(2)
       else
          grad(i,1) = tmp(1)
          grad(i,2) = tmp(2)
       end if
    end do

    do i = 0, nth_used
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call self%eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(2), 'Z')
       grad(i,1) = tmp(1)
       grad(i,2) = tmp(2)
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

  !> FIXME : Add documentation
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(peq_base_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, self%nt, 2) ::  dbish
    integer :: i

    select case(char)
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
    enddo

    if (char == 'T' .and. .not. self%transp) then
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
          grad(i,1) = grad(i,1)*daa(1) * 0.5*self%beta_0
          grad(i,2) = grad(i,2)*daa(1) * 0.5*self%beta_0
       end do
    end if
  end subroutine bgradient

  !> FIXME : Add documentation
  subroutine eqitem(self, r_in, theta_in, f, fstar, char)
    use mp, only: mp_abort
    use constants, only: pi
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: r_in, theta_in
    real, dimension(:,:), intent(in) :: f
    character(1), intent(in) :: char
    integer :: i, j, istar, jstar
    real, intent(out) :: fstar
    real :: thet, sign, tp, tps
    real :: st, dt, sr, dr, r
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

    if (r < sign * self%eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r, sign, self%eqpsi(1)
       call mp_abort('r < Psi_0 in eqitem')
    end if

    ! find r on psi mesh

    ! disallow evaluations on or outside the surface for now
    if (r >= self%eqpsi(self%nr)) then
       write(*,*) 'No evaluation of eqitem allowed on or outside surface'
       write(*,*) '(Could this be relaxed a bit?)'
       write(*,*) r, theta_in, self%eqpsi(self%nr), sign
       call mp_abort('No evaluation of eqitem allowed on or outside surface')
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
    if (.not. self%transp) then
       tp = abs(thet)
       tps = 1.
       if (char == 'Z' .and. tp >= 1.e-10) tps = thet / tp

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
          jstar = self%nt - 1
          dt = mtheta(jstar+1) - mtheta(jstar)
          st = 0.
       end if

       ! use opposite area stencil to interpolate
       fstar=f(istar    , jstar    ) * sr * st &
            +f(istar + 1, jstar    ) * dr * st &
            +f(istar    , jstar + 1) * sr * dt &
            +f(istar + 1, jstar + 1) * dr * dt
       fstar = fstar * tps / (abs(self%eqpsi(istar+1) - self%eqpsi(istar)) * &
            (mtheta(jstar+1) - mtheta(jstar)))

    else

       ! get thet on theta mesh
       do j = 1, self%nt
          mtheta(j) = (j-1)*2.*pi / float(self%nt-1) - pi
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
          jstar = self%nt - 1
          dt = mtheta(jstar+1) - mtheta(jstar)
          st = 0.
       end if

       ! use opposite area stencil to interpolate
       fstar=f(istar    , jstar    ) * sr * st &
            +f(istar + 1, jstar    ) * dr * st &
            +f(istar    , jstar + 1) * sr * dt &
            +f(istar + 1, jstar + 1) * dr * dt
       fstar=fstar / (abs(self%eqpsi(istar+1) - self%eqpsi(istar)) * &
            (mtheta(jstar+1) - mtheta(jstar)))
    end if
  end subroutine eqitem

  !> FIXME : Add documentation
  pure subroutine eqdcart(self, dfm, dfcart)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (self%nr, self%nt) :: denom
    integer :: i, j
    denom = self%drm(:,:,1) * self%dzm(:,:,2) - self%drm(:,:,2) * self%dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * self%dzm(:,:,2) - self%dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * self%drm(:,:,2) + self%drm(:,:,1) * dfm(:,:,2)

    do j = 1, self%nt
       do i = 2, self%nr
          dfcart(i,j,:) = dfcart(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdcart

  !> FIXME : Add documentation
  subroutine eqdbish(self, dcart, dbish)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j
    denom = sqrt(self%dpcart(:,:,1)**2 + self%dpcart(:,:,2)**2)
    dbish(:,:,1) = dcart(:,:,1) * self%dpcart(:,:,1) + dcart(:,:,2) * self%dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * self%dpcart(:,:,2) + dcart(:,:,2) * self%dpcart(:,:,1)

    do j = 1, self%nt
       do i = 2, self%nr
          dbish(i,j,:) = dbish(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdbish

  !> FIXME : Add documentation
  real function invR (self, r, theta)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%Rpos(r, theta)
  end function invR

  !> FIXME : Add documentation
  real function Rpos (self, r, theta)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%R_psi, Rpos, 'R')
  end function Rpos

  !> FIXME : Add documentation
  real function Zpos (self, r, theta)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    call self%eqitem(r, th, self%Z_psi, Zpos, 'Z')
  end function Zpos

  !> FIXME : Add documentation
  real function psi (self, r, theta)
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent (in) :: r, theta
    psi = r
  end function psi

  !> FIXME : Add documentation
  real function diameter (self, rp)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = splint(rp, self%diam_spl)
  end function diameter

  !> FIXME : Add documentation
  real function rcenter (self, rp)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: rp
    rcenter = splint(rp, self%rcenter_spl)
  end function rcenter

  !> FIXME : Add documentation
  real function dbtori (self, pbar)
    use splines, only: dsplint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    dbtori = dsplint(pbar, self%btori_spl) / (self%psi_a - self%psi_0)
  end function dbtori

  !> FIXME : Add documentation
  real function btori (self, pbar)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    btori = splint(pbar, self%btori_spl)
  end function btori

  !> FIXME : Add documentation
  real function qfun (self, pbar)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    qfun = splint(pbar, self%qsf_spl)
  end function qfun

  !> FIXME : Add documentation
  real function pfun (self, pbar)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    pfun = 0.5 * self%beta_0 * splint(pbar, self%pressure_spl)
  end function pfun

  !> FIXME : Add documentation
  real function dpfun (self, pbar)
    use splines, only: dsplint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    dpfun =  0.5 * self%beta_0 *dsplint(pbar, self%pressure_spl) / (self%psi_a-self%psi_0)
  end function dpfun

  !> FIXME : Add documentation
  real function betafun (self, pbar)
    use splines, only: splint
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    if(pbar == 0.) then
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
    type(peq_type), intent(in) :: geom
    real, intent(in) :: a
    integer, intent(in) :: irho
    integer :: i, fort24_unit, nr
    real :: gradpsi, mag_B, rho_eq, rp1, rp2, rho1, rho2, drhodpsiq
    real, dimension(geom%nr) :: gamma, pbar, dp, d2p, pres
    if (.not. proc0) return
    nr = geom%nr

    gamma = 0.

    pbar = (geom%eqpsi - geom%eqpsi(1)) / (geom%eqpsi(geom%nr) - geom%eqpsi(1))

    do i = 2, nr-1
       dp(i)  = geom%dpfun(pbar(i))
    end do

    do i = 3, nr-2
       d2p(i) = (dp(i+1)-dp(i-1))/(geom%eqpsi(i+1)-geom%eqpsi(i-1))
    enddo

    pres=0.
    do i = 3, nr-2
       rp1=geom%eqpsi(i+1)
       rp2=geom%eqpsi(i-1)
       rho1=0.5*geom%diameter(rp1)
       rho2=0.5*geom%diameter(rp2)
       drhodpsiq=(rho1-rho2)/(rp1-rp2)

       call geom%eqitem(geom%eqpsi(i), 0., geom%dpbish(:,:,1), gradpsi, 'R')

       mag_B=sqrt( (geom%btori(pbar(i))**2 + gradpsi**2))/geom%Rpos(geom%eqpsi(i),0.)

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
end module peq
