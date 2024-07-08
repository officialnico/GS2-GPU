!> Reads in the geometry using a CHEASE output file.
!> The CHEASE output file is read using the helper module read_chease.
module ceq
  use geo_utils, only: abstract_eqfile_geo_type, geo_input_type
  implicit none
  public :: ceq_type, chease_chi_index

  private

  type, extends(abstract_eqfile_geo_type) :: ceq_type
     private
   contains
     procedure :: read_and_set, set_ceq_from_chease
  end type ceq_type

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

  subroutine read_and_set(self, inputs)
    use read_chease, only: read_infile, finish
    use constants, only: pi, twopi
    implicit none
    class(ceq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real :: f_N, psi_N
    integer :: j
    logical, parameter :: debug = .true.
    self%type_name = 'ceq'
    self%filename = trim(adjustl(inputs%eqfile))
    if (.not. skip_file_read) then ! Primarily for testing
       write (*,*)  'Reading CHEASE input file: ', trim(self%filename)
       !    Read the data
       call read_infile(self%filename)
       call self%set_ceq_from_chease
       call finish
    end if

    !     Normalize, rename quantities
    psi_N = abs(self%B_T) * self%aminor**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%eqpsi = self%eqpsi / psi_N
    f_N = abs(self%B_T) * self%aminor
    self%fp = self%fp / f_N

    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * twopi / real(self%nt-1) - pi
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    self%diam = abs(self%R_psi(:, self%nt/2 + 1) - self%R_psi(:, 1))
    self%rc = 0.5*(self%R_psi(:, 1) + self%R_psi(:, self%nt/2 + 1))
    self%has_full_theta_range = .true.

    if (debug) then
       write (*,*) "Finished ceqin... imported CHEASE equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", self%aminor
       write (*,*) 'R_mag', self%R_mag
       write (*,*) 'B_T0', abs(self%B_T)
       write (*,*) 'f_N', abs(self%B_T) * self%aminor
       write (*,*) 'nthg', self%nt
       write (*,*) 'beta', self%beta_0
    end if
  end subroutine read_and_set

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
    real :: R_geo
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
    self%psi_0 = self%eqpsi(1) ; self%psi_a = self%eqpsi(self%nr)
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
    R_geo = (self%R_psi(self%nr, self%nt/2+1) + self%R_psi(self%nr, 1))/2.
    self%aminor = (self%R_psi(self%nr, self%nt/2+1) - self%R_psi(self%nr, 1))/2.
    self%B_T = self%fp(self%nr)/R_geo
    self%R_psi = self%R_psi / self%aminor
    self%Z_psi = self%Z_psi / self%aminor
    self%B_psi = 0.0

    !     R_mag == R position of magnetic axis
    !     Z_mag == Z position of magnetic axis
    self%R_mag = self%R_psi(1,1)
    self%Z_mag = self%Z_psi(1,1)
    self%beta = 8. * pi * 1.e-7 * self%pressure / b0exp_chease**2  !not correct definition yet
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)
  end subroutine set_ceq_from_chease
end module ceq
