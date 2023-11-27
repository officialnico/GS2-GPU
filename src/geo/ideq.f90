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
  end subroutine alloc_special_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(ideq_type), intent(in out) :: self
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
end module ideq
