!> Reads in the geometry using a CHEASE output file.
!> The CHEASE output file is read using the helper module read_chease.
module ceq
  use geo_utils, only: abstract_eqfile_geo_type, mod2pi

  implicit none

  public :: Hahm_Burrell, ceq_type

  private

  type, extends(abstract_eqfile_geo_type) :: ceq_type
     private
   contains
     procedure :: initialise => ceqin
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: eqitem => eqitem
     procedure :: alloc_special_arrays => alloc_special_arrays
     procedure :: dealloc_special_arrays => dealloc_special_arrays
     procedure :: set_ceq_from_chease => set_ceq_from_chease
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
  subroutine ceqin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use read_chease, only: read_infile, finish
    use geo_utils, only: geo_input_type
    use constants, only: pi
    implicit none
    class(ceq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(out) :: nthg
    real :: f_N, psi_N
    integer :: j
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = self%nt
       return
    end if

    self%filename = trim(adjustl(inputs%eqfile))
    if (.not. skip_file_read) then ! Primarily for testing
       write (*,*)  'Reading CHEASE input file: ', self%filename
       !    Read the data
       call read_infile(self%filename)
       call self%set_ceq_from_chease
       call finish
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

    do j = 1, self%nt
       self%eqth(:, j) = (j-1)*2.*pi / real(self%nt-1) - pi
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

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

    self%B_psi = 0.0
    self%initialised = .true.
    self%has_full_theta_range = .true.
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
  pure subroutine alloc_special_arrays(self)
    implicit none
    class(ceq_type), intent(in out) :: self
    integer :: nr
    nr = self%nr
    allocate(self%psi_bar(nr), self%fp(nr), self%qsf(nr), self%beta(nr))
    allocate(self%rc(nr), self%diam(nr))
  end subroutine alloc_special_arrays

  !> Deallocate module-level arrays
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(ceq_type), intent(in out) :: self
    deallocate(self%psi_bar, self%fp, self%qsf, self%beta, self%rc, self%diam)
  end subroutine dealloc_special_arrays

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

  !> Calculates fstar which is f interpolated at the location (r,theta).
  !> Here r is the normalised poloidal flux coordinate rp (= psi_pN + psi_0N)
  !> and theta_in is theta. f is a grid of values of f as a function of psi_p,theta
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    class(ceq_type), intent(in) :: self
    real, intent(in) :: r_in, theta_in
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
    mtheta = self%eqth(1, :)

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
  end function eqitem

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

       gradpsi = geom%eqitem(geom%eqpsi(i), 0., geom%dpbish(:,:,1), 'R')

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
