!> FIXME : Add documentation
module ideq
  use geo_utils, only: abstract_eqfile_geo_type, geo_input_type
  implicit none
  private

  public :: ideq_type

  type, extends(abstract_eqfile_geo_type) :: ideq_type
     private
   contains
     procedure :: read_and_set
  end type ideq_type

contains
  !> FIXME : Add documentation
  subroutine read_and_set(self, inputs)
    use constants, only: twopi, pi
    use splines, only: new_periodic_spline, delete_periodic_spline, &
         periodic_splint, periodic_spline
    use geo_utils, only: sort
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
    end type grid_type

    !> FIXME : Add documentation
    type :: eq_type
       integer :: nr
       real, dimension (:), allocatable :: pbar
       real, dimension (:), allocatable :: pressure
       real, dimension (:), allocatable :: volume
       type (grid_type), dimension (:), allocatable :: grid
    end type eq_type
    type (periodic_spline) :: spl
    type (eq_type) :: d
    real :: xdum, I_ring, dpsidr, dpsidz, psi_N
    integer :: i, j, idum, jj, jmax, funit
    character (200) :: line
    logical, dimension(:), allocatable :: has_duplicate_point
    self%type_name = 'ideq'
    ! Need to generalize initialization condition if equilibrium changes
    self%nt = size(inputs%theta)

    self%filename = trim(adjustl(inputs%eqfile))

    open(newunit = funit, file = self%filename, status = 'old', form = 'formatted')

    ! Read the data
    read (funit, fmt="(a)") line
    read (funit, fmt="(a)") line
    read (funit,*) self%R_mag, self%Z_mag, I_ring
    read (funit, fmt="(a)") line
    read (funit, *) self%nr
    read (funit, fmt="(a)") line
    read (funit, fmt="(a)") line
    call self%alloc_arrays

    d%nr = self%nr
    allocate (d%pbar     (d%nr))
    allocate (d%pressure (d%nr))
    allocate (d%volume   (d%nr))
    allocate (d%grid     (d%nr))

    ! Read in each surface
    do i = 1, d%nr
       read(funit,*) idum, d%pbar(i), self%eqpsi(i), d%pressure(i), xdum, d%volume(i), d%grid(i)%nt
       d%grid(i)%nt = d%grid(i)%nt-1
       jmax = d%grid(i)%nt
       allocate (d%grid(i)%R     (jmax))
       allocate (d%grid(i)%Z     (jmax))
       allocate (d%grid(i)%B     (jmax))
       allocate (d%grid(i)%theta (jmax))
       do j = 1, jmax
          read (funit,*) d%grid(i)%R(j), d%grid(i)%Z(j), dpsidR, dpsidZ
          ! make B:
          d%grid(i)%B(j) = hypot(dpsidR, dpsidZ) / (d%grid(i)%R(j)*twopi)
          d%grid(i)%theta(j) = atan2 (d%grid(i)%Z(j) - self%Z_mag, &
               (d%grid(i)%R(j) - self%R_mag))
       end do
       read (funit,*) xdum, xdum, xdum, xdum
    end do
    close(funit)

    ! different from tokamak case here, psi_0 is the value of psi on
    ! the innermost flux surface that is kept in equilibrium and psi_a
    ! is the last flux surface kept.
    self%psi_0 = self%eqpsi(1) ; self%psi_a = self%eqpsi(d%nr)

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
    allocate(has_duplicate_point(d%nr))

    ! Need to determine if we have a duplicate grid point as the spline routines
    ! require that we do not include the duplicate point
    do i = 1, d%nr
       has_duplicate_point(i) = twopi - &
            (d%grid(i)%theta(d%grid(i)%nt) - d%grid(i)%theta(1)) <= epsilon(twopi)
    end do

    do i = 1, d%nr
       spl = new_periodic_spline (d%grid(i)%theta(:), d%grid(i)%R(:), &
            twopi, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%R_psi (i, jj) = periodic_splint (inputs%theta(jj), spl)
       end do
       call delete_periodic_spline (spl)
    end do

    do i = 1, d%nr
       spl = new_periodic_spline (d%grid(i)%theta(:), d%grid(i)%Z(:), &
            twopi, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%Z_psi (i, jj) = periodic_splint (inputs%theta(jj), spl)
       end do
       call delete_periodic_spline (spl)
    end do

    do i = 1, d%nr
       spl = new_periodic_spline (d%grid(i)%theta(:), d%grid(i)%B(:), &
            twopi, drop_last_point = has_duplicate_point(i))
       do jj = 1, self%nt
          self%B_psi (i, jj) = periodic_splint (inputs%theta(jj), spl)
       end do
       call delete_periodic_spline (spl)
    end do

    ! recover some memory:
    deallocate (d%grid, d%pbar, d%pressure, d%volume)

    ! B at center of ring is mu_0 I / (2 a)
    self%B_T = 4.*pi*1.e-7*I_ring*0.5 / self%R_mag

    write(*,*) 'B_T0 = ',self%B_T
    self%aminor = self%R_mag
    psi_N = self%B_T * self%aminor**2
    self%eqpsi = self%eqpsi / psi_N
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%R_psi = self%R_psi / self%R_mag
    self%Z_psi = self%Z_psi / self%R_mag
    self%R_mag = self%R_mag / self%aminor !Normalised to itself
    self%B_psi = self%B_psi / self%B_T
    self%rc = 0.5*(self%R_psi(:, 1) + self%R_psi(:, (self%nt-1)/2))
    ! need to define the rho_d grid (called diam)
    !> Not really the diameter in this case.  Instead, return the
    !> normalized minor radius, measured inside the ring, in the plane of
    !> the ring, starting at the ring and going inward.
    self%diam = 1. - self%R_psi(:,1)

    do j = 1, self%nt
       self%eqth(:, j) = (j-1) * pi / float(self%nt-1)
       self%eqpsi_2d(:, j) = self%eqpsi
    end do

    self%qsf = 0.
    self%fp = 0.
    ! Note beta and pressure were allocated but not set historically.
    ! As such now forcing these to zero, but we _do_ read in pressure
    ! so we could probably actually set this?
    self%beta = 0.
    self%pressure= 0.
    self%beta_0 = self%beta(1)
    self%has_full_theta_range = .false.
  end subroutine read_and_set
end module ideq
