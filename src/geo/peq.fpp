!> FIXME : Add documentation
module peq
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

  public :: Hahm_Burrell, peq_type, teq_type

  type, abstract, extends(abstract_eqfile_geo_type) :: peq_base_type
   contains
     procedure :: shared_setup
  end type peq_base_type

  type, extends(peq_base_type) :: peq_type
   contains
     procedure :: read_and_set => read_and_set_peq
  end type peq_type

  type, extends(peq_base_type) :: teq_type
   contains
     procedure :: read_and_set => read_and_set_teq
  end type teq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file containing
  !> the axisymmetric magnetic field geometry in flux coordinates
  subroutine read_and_set_peq(self, inputs)
    use constants, only: pi, twopi
    implicit none
    class(peq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    integer :: istatus, ncid, id, i, j, i_sym, nz1, nz2
    real, allocatable, dimension(:) :: work
    self%filename = trim(adjustl(inputs%eqfile))

# ifdef NETCDF
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

    if (i_sym == 0) write(*,*) 'Up-down asymmetries not yet handled correctly.'

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

    !     netcdf read 2d field: R_psi,Z_psi and B_psi (mod(B))
    !     R_psi(1:nr, 1:nt) -- transpose of Menard's storage ...
    !     Z_psi(1:nr, 1:nt)
    !     B_psi(1:nr, 1:nt)
    allocate(work(nz1*nz2))

    istatus = nf90_inq_varid (ncid, 'x', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='x')
    istatus = nf90_get_var (ncid, id, work, count = [nz1*nz2])
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    do j = 1, self%nt
       do i = 1, self%nr
          self%R_psi(i,j) = work(3+j-1+nz1*(i-1))
       end do
    end do

    istatus = nf90_inq_varid (ncid, 'z', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='z')
    istatus = nf90_get_var (ncid, id, work, count=[nz1, nz2]) !Should this be nz1*nz2?
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    do j = 1, self%nt
       do i = 1, self%nr
          self%Z_psi(i,j) = work(3+j-1+nz1*(i-1))
       end do
    end do

    istatus = nf90_close(ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)
    deallocate(work)

    self%has_full_theta_range = .false.
    call self%shared_setup(self%nt, 1)
# else
    call mp_abort('error: peq eqin is called without netcdf',.true.)
# endif
  end subroutine read_and_set_peq

  !> This subroutine reads a generic NetCDF equilibrium file
  !> containing the axisymmetric magnetic field geometry in flux
  !> coordinates
  subroutine read_and_set_teq(self, inputs)
    use constants, only: pi, twopi
    implicit none
    class(teq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    integer :: istatus, ncid, id
    real, allocatable, dimension(:,:) :: work
    self%filename = trim(adjustl(inputs%eqfile))

# ifdef NETCDF
    istatus = nf90_open(self%filename, NF90_NOWRITE, ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file = self%filename)

    !     nr == number of radial grid points in radial eq grid
    istatus = nf90_inq_varid (ncid, 'ns', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='ns')
    istatus = nf90_get_var (ncid, id, self%nr)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    !     nt == number of theta grid points in theta eq grid
    istatus = nf90_inq_varid (ncid, 'nt1', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='nt1')
    istatus = nf90_get_var (ncid, id, self%nt)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    call self%alloc_arrays()

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

    istatus = nf90_inq_varid (ncid, 'g', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='g')
    istatus = nf90_get_var (ncid, id, self%fp)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'q', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='q')
    istatus = nf90_get_var (ncid, id, self%qsf)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    istatus = nf90_inq_varid (ncid, 'p', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='p')
    istatus = nf90_get_var (ncid, id, self%pressure)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    allocate(work(self%nt, self%nr))

    !     netcdf read 2d field: R_psi and Z_psi
    istatus = nf90_inq_varid (ncid, 'R', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='R')
    istatus = nf90_get_var (ncid, id, work)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    self%R_psi = transpose(work)

    istatus = nf90_inq_varid (ncid, 'Z', id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='Z')
    istatus = nf90_get_var (ncid, id, work)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

    self%Z_psi = transpose(work)

    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus)
    deallocate(work)

    self%eqpsi = abs(self%eqpsi) ; self%fp = abs(self%fp) ; self%qsf = abs(self%qsf)
    self%has_full_theta_range = .true.
    call self%shared_setup(self%nt, self%nt/2)
# else
    call mp_abort('error: peq teqin is called without netcdf',.true.)
# endif
  end subroutine read_and_set_teq

  subroutine shared_setup(self, tind_in, tind_out)
    use constants, only: twopi, pi
    class(peq_base_type), intent(in out) :: self
    integer, intent(in) :: tind_in, tind_out !< Theta index of inboard and outboard
    real :: R_geo, psi_N, f_n
    integer :: j
    self%psi_0 = self%eqpsi(1) ; self%psi_a = self%eqpsi(self%nr)
    self%psi_bar = (self%eqpsi - self%psi_0) / (self%psi_a - self%psi_0)
    self%aminor = 0.5 * abs(self%R_psi(self%nr, tind_out) - self%R_psi(self%nr, tind_in))
    R_geo = 0.5 * (self%R_psi(self%nr, tind_out) + self%R_psi(self%nr, tind_in))
    self%B_T = self%fp(self%nr) / R_geo ; self%B_psi = 0.0
    self%R_psi = self%R_psi / self%aminor ; self%Z_psi = self%Z_psi / self%aminor
    self%R_mag = self%R_psi(1, 1) ; self%Z_mag = self%Z_psi(1, 1)
    self%beta = 4. * twopi * 1.e-7 * self%pressure / self%B_T**2 ! not correct definition yet?
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    psi_N = abs(self%B_T) * self%aminor**2
    self%psi_a = self%psi_a / psi_N ; self%psi_0 = self%psi_0 / psi_N
    self%eqpsi = self%eqpsi / psi_N
    f_N = abs(self%B_T) * self%aminor
    self%fp = self%fp / f_N

    self%diam = abs(self%R_psi(:, tind_out) - self%R_psi(:, tind_in))
    self%rc = 0.5 * (self%R_psi(:, tind_out) + self%R_psi(:, tind_in))
    do j = 1, self%nt
       self%eqth(:, j) = (j - 1) * pi / float(self%nt - 1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do
    if (self%has_full_theta_range) self%eqth = 2 * self%eqth - pi
  end subroutine shared_setup

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
