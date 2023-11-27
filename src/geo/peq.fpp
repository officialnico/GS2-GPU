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
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi
  implicit none

  private

  public :: Hahm_Burrell, peq_type, teq_type

  type, abstract, extends(abstract_eqfile_geo_type) :: peq_base_type
     private
     logical :: transp
   contains
     procedure :: calculate_gradients => peq_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     procedure :: alloc_special_arrays => alloc_special_arrays
     procedure :: dealloc_special_arrays => dealloc_special_arrays
  end type peq_base_type

  type, extends(peq_base_type) :: peq_type
   contains
     procedure :: initialise => eqin
  end type peq_type

  type, extends(peq_base_type) :: teq_type
   contains
     procedure :: initialise => teqin
  end type teq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file
  !> containing the axisymmetric magnetic field geometry in flux
  !> coordinates
  subroutine eqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use mp, only: mp_abort
    use constants, only: pi
    use geo_utils, only: geo_input_type
    implicit none
    class(peq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(out) :: nthg
    real :: d, R_geo
    integer :: istatus
    integer :: ncid, id, i, j
    integer :: nchar, i_sym
    integer nz1, nz2
    real, allocatable, dimension(:) :: work
    real :: f_N, psi_N
    !     read the data
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = self%nt
       return
    end if

# ifdef NETCDF
    self%filename = trim(adjustl(inputs%eqfile))

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
    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * pi / float(self%nt - 1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    self%B_psi = 0.0
    self%transp = .false.
    self%initialised = .true.
    self%has_full_theta_range = .false.
# else
    call mp_abort('error: peq eqin is called without netcdf',.true.)
# endif
  end subroutine eqin

  !> This subroutine reads a generic NetCDF equilibrium file
  !> containing the axisymmetric magnetic field geometry in flux
  !> coordinates
  subroutine teqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use mp, only: mp_abort
    use constants, only: pi, twopi
    use geo_utils, only: geo_input_type
    implicit none
    class(teq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(out) :: nthg
    real :: R_geo
    integer :: istatus
    integer :: ncid, id, i, j
    integer :: nchar
    ! what is the best way to handle the netcdf single/double problem?
    real, allocatable, dimension(:,:) :: work
    real :: f_N, psi_N
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag / self%aminor
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

    ! find magnetic axis -- Note currently we don't normalise these
    ! whilst we _do_ for eqin (transp = .false.). They are not used
    ! so this is OK for now.
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
    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * twopi / float(self%nt - 1) - pi
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    self%B_psi = 0.0
    self%transp = .true.
    self%initialised = .true.
    self%has_full_theta_range = .true.
# else
    call mp_abort('error: peq teqin is called without netcdf',.true.)
# endif
  end subroutine teqin

  !> Allocate module-level arrays
  pure subroutine alloc_special_arrays(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    integer :: nr
    nr = self%nr
    allocate(self%psi_bar(nr), self%fp(nr), self%qsf(nr), self%beta(nr))
    allocate(self%rc(nr), self%diam(nr))
  end subroutine alloc_special_arrays

  !> Deallocate module-level arrays
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    deallocate(self%psi_bar, self%fp, self%qsf, self%beta, self%rc, self%diam)
  end subroutine dealloc_special_arrays

  !> FIXME : Add documentation
  subroutine peq_init(self)
    implicit none
    class(peq_base_type), intent(in out) :: self
    call self%derm(self%eqth, self%dtm, 'T')
    call self%derm(self%R_psi, self%drm, 'E')
    call self%derm(self%Z_psi, self%dzm, 'O')
    call self%derm(self%eqpsi_2d, self%dpm, 'E')

    ! grad(psi) in cartesian form
    call self%eqdcart(self%dpm, self%dpcart)
    ! grad(psi) in Bishop form
    call self%eqdbish(self%dpcart, self%dpbish)

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
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nr, self%nt, 2) :: dcart
    integer :: i

    select case(char)
    case('P', 'R')
       dcart = self%dpcart
    case('T')
       dcart = self%dtcart
    end select

    do i = -nth_used, -1
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,2), 'Z')
       if (char == 'T' .and. .not. self%transp) then
          grad(i, 1) = -grad(i, 1)
          grad(i, 2) = -grad(i, 2)
       end if
    end do

    do i = 0, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'Z')
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
    case('P', 'R')
       dbish = self%dpbish
    case('T')
       dbish = self%dtbish
    end select

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dbish(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dbish(:,:,2), 'Z')
    end do

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
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    implicit none
    class(peq_base_type), intent(in) :: self
    real, intent(in) :: r_in, theta_in
    real, dimension(:,:), intent(in) :: f
    character(1), intent(in) :: char
    integer :: i, j, istar, jstar
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
    mtheta = self%eqth(1, :)

    ! assume up-down symmetry
    if (.not. self%transp) then
       tp = abs(thet)
       tps = 1.
       if (char == 'Z' .and. tp >= 1.e-10) tps = thet / tp

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
  end function eqitem

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

       gradpsi = geom%eqitem(geom%eqpsi(i), 0., geom%dpbish(:,:,1), 'R')

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
