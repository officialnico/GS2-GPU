!> FIXME : Add documentation
module geq
# ifdef NETCDF
  use netcdf, only: NF90_NOWRITE, NF90_NOERR, nf90_open, nf90_close, nf90_get_var
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid
  use netcdf_utils, only: netcdf_error
# else
  use mp, only: mp_abort
# endif
  use geo_utils, only: abstract_eqfile_geo_type, geo_input_type
  implicit none
  private

  public :: geq_type

  type, extends(abstract_eqfile_geo_type) :: geq_type
     private
   contains
     procedure :: read_and_set
  end type geq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file containing
  !> the axisymmetric magnetic field geometry in flux coordinates
  subroutine read_and_set(self, inputs)
    use constants, only: pi
    implicit none
    class(geq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, parameter :: psi_scale = 1.0e-8, b_scale = 1.0e-4, l_scale = 1.0e-2
    real :: f_N, psi_N
    integer :: istatus, ncid, id, j
    self%type_name = 'geq'
    self%filename = trim(adjustl(inputs%eqfile))
# ifdef NETCDF
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

    istatus = nf90_inq_varid (ncid, 'psi_a', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi_a')
    istatus = nf90_get_var (ncid, id, self%psi_a)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'B_T', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='B_T')
    istatus = nf90_get_var (ncid, id, self%B_T)
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

    istatus = nf90_inq_varid (ncid, 'Zmag', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Zmag')
    istatus = nf90_get_var (ncid, id, self%Z_mag)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'a', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a')
    istatus = nf90_get_var (ncid, id, self%aminor)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

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
    istatus = nf90_inq_varid (ncid, 'psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psi')
    istatus = nf90_get_var (ncid, id, self%eqpsi)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'psibar', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psibar')
    istatus = nf90_get_var (ncid, id, self%psi_bar)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'fp', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='fp')
    istatus = nf90_get_var (ncid, id, self%fp)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'q', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='q')
    istatus = nf90_get_var (ncid, id, self%qsf)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'beta', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='beta')
    istatus = nf90_get_var (ncid, id, self%beta)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

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

    istatus = nf90_inq_varid (ncid, 'Z_psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Z_psi')
    istatus = nf90_get_var (ncid, id, self%Z_psi, count = [self%nr, self%nt])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'B_psi', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='B_psi')
    istatus = nf90_get_var (ncid, id, self%B_psi, count = [self%nr, self%nt])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)

    ! Convert units
    self%psi_0 = self%psi_0 * psi_scale
    self%psi_a = self%psi_a * psi_scale
    self%eqpsi = self%eqpsi * psi_scale
    self%B_T = self%B_T * b_scale  ! cgs to MKS
    self%B_psi = self%B_psi * b_scale / self%B_T
    self%R_psi = self%R_psi / self%aminor
    self%Z_psi = self%Z_psi / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    self%aminor = self%aminor * l_scale

    self%fp = self%fp * b_scale * l_scale
    self%beta_0 = self%beta(1)

    !     Normalize, rename quantities
    psi_N = self%B_T * self%aminor**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%eqpsi = self%eqpsi / psi_N
    f_N = self%B_T * self%aminor
    self%fp = self%fp / f_N

    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * pi / real(self%nt-1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    self%diam = self%R_psi(:, 1) - self%R_psi(:, self%nt)
    self%rc = 0.5 * (self%R_psi(:, 1) + self%R_psi(:, self%nt))
    self%has_full_theta_range = .false.
# else
    call mp_abort('error: geq eqin is called without netcdf',.true.)
#endif
  end subroutine read_and_set
end module geq
