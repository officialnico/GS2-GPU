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
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi
  implicit none

  private

  public :: Hahm_Burrell, geq_type

  type, extends(abstract_eqfile_geo_type) :: geq_type
     private
   contains
     procedure :: initialise => eqin
  end type geq_type

contains

  !> This subroutine reads a generic NetCDF equilibrium file
  !! containing the axisymmetric magnetic field geometry in flux
  !! coordinates
  subroutine eqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use geo_utils, only: geo_input_type
    use constants, only: pi
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
    integer :: istatus, ncid, id, j
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

    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * pi / real(self%nt-1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    f_N = B_T0 * avgrmid
    self%fp = self%fp / f_N

    self%diam = self%R_psi(:, 1) - self%R_psi(:, self%nt)
    self%rc = 0.5 * (self%R_psi(:, 1) + self%R_psi(:, self%nt))

    nthg = self%nt
# else
    call mp_abort('error: geq eqin is called without netcdf',.true.)
#endif

    self%initialised = .true.
    self%has_full_theta_range = .false.
  end subroutine eqin

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
