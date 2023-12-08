!> FIXME : Add documentation
module deq
  use geo_utils, only: abstract_eqfile_cart_geo_type
  use splines, only: spline

  implicit none

  private

  public :: deq_type

  type, extends(abstract_eqfile_cart_geo_type) :: deq_type
     private
     integer :: nrho
     real, allocatable, dimension(:) :: rho_mid, psi_mid
     type(spline) :: rho_spl
   contains
     procedure :: initialise => dfitin
     procedure :: setup_special_splines => setup_special_splines
     procedure :: delete_special_splines => delete_special_splines
     procedure :: rcenter => rcenter
     procedure :: rhofun => rhofun
     procedure :: alloc_special_arrays => alloc_special_arrays
     procedure :: dealloc_special_arrays => dealloc_special_arrays
  end type deq_type

contains

  !> This subroutine reads an DFIT output file containing
  !> the axisymmetric magnetic field geometry on a rectangular
  !> domain defined by the coordinates (R,Z).
  !
  !> eq_R     R grid
  !> eq_Z     Z grid
  !> fp    F on psibar grid
  !> eqpsi_2d   psi on (R,Z) grid
  !> R_mag        major radius of the magnetic axis
  !> Z_mag        elevation of the magnetic axis
  !> rwid       total width of the domain
  !> rleft      position of leftmost point of domain
  !> zhei       total height of domain
  subroutine dfitin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use mp, only: mp_abort
    use constants, only: pi, twopi
    use splines, only: inter_cspl
    use geo_utils, only: sort, geo_input_type
    implicit none
    class(deq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, B_T0, avgrmid, rmaj
    integer, intent(out) :: nthg
    real, dimension(:), allocatable :: spsi_bar, sdfit_R, sdfit_Z, rbbbs, zbbbs, dummy, eq_R, eq_Z
    real :: p_0, rwid, zhei, bcentr, delta_R, delta_Z, lref, bref
    integer :: i, j, jmin, jmax, nw_in, nh_in, funit, nbbbs
    logical :: eqnorm_exists
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = -1
       return
    end if
    ! Need to generalize initialization condition if equilibrium changes
    lref = 0.
    bref = 0.
    ! read eqnorm file if exists
    inquire(file = trim(inputs%eqnormfile), exist = eqnorm_exists)
    if (eqnorm_exists) then
       open(newunit = funit, file = trim(inputs%eqnormfile), status='old', form='formatted')
       read(funit, *) lref, bref
       close(funit)
    end if

    self%filename = trim(adjustl(inputs%eqfile))
    open(newunit = funit, file = self%filename, status='old', form='formatted')

    ! Read the data
    read(funit, *) nw_in, nh_in, self%nrho
    ! Note here we don't allow for rescaling the final grids with big
    ! unlike in eeq. This makes many of the spline interpolation calls
    ! in this routine effectively just copies
    self%nw = nw_in ; self%nh = nh_in
    self%nr = self%nw ; self%nt = self%nh
    call self%alloc_arrays
    allocate(spsi_bar(nw_in), sdfit_R(nw_in), sdfit_Z(nh_in), dummy(nw_in))
    allocate(eq_R(nw_in), eq_Z(nh_in))
    read(funit, 2020) rwid, zhei
    read(funit, 2020) self%R_mag, self%Z_mag
    read(funit, 2020) self%psi_0, self%psi_a, bcentr

    ! pbar is defined by
    ! pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    do i = 1, nw_in
       spsi_bar(i) = float(i) / float(nw_in)
       sdfit_R(i) = rwid * float(i) / float(nw_in)
    end do

    do j = 1, nh_in
       sdfit_Z(j) = zhei*(float(j-1)/float(nh_in-1)-0.5)
    end do

    do i = 1, self%nw
       self%psi_bar(i) = float(i) / float(self%nw)
       eq_R(i) = rwid * float(i) / float(self%nw)
    end do

    do j = 1, self%nh
       eq_Z(j) = zhei*(float(j-1)/float(self%nh-1)-0.5)
    end do

    read(funit, 2020) (dummy(j), j = 1, nw_in)  ! pressure read
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%pressure)
    ! Note we read directly into the eqpsi_2d array here so self%nw == nw_in etc.
    ! has to be true for this to work. See efit_in for a more general approach.
    read(funit, 2020) ((self%eqpsi_2d(i,j) , i = 1, nw_in) , j = 1, nh_in)

    ! Only needed for irho = 4
    read (funit, 2020) (self%rho_mid(i), i = 1, self%nrho)
    read (funit, 2020) (self%psi_mid(i), i = 1, self%nrho)
    close(funit)

    ! stay away from edge of box
    delta_R = 1.02 ; delta_Z = 1.02

    nbbbs = 2 * (self%nw + self%nh) - 5
    allocate(rbbbs(nbbbs), zbbbs(nbbbs))
    allocate(self%thetab(nbbbs), self%r_bound(nbbbs))
    jmin = 1
    jmax = self%nh / 2
    do j = jmin, jmax
       rbbbs(j) = (eq_R(1) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = eq_Z(j + self%nh / 2) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nw - 2
    do j = jmin, jmax
       rbbbs(j) = (eq_R(j - jmin + 2) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = eq_Z(self%nh) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nh - 1
    do j = jmin, jmax
       rbbbs(j) = (eq_R(self%nw) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = eq_Z(jmax - j + 2) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nw - 1
    do j = jmin, jmax
       rbbbs(j) = (eq_R(jmax - j + 2) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = eq_Z(1) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nh/2
    do j = jmin, jmax
       rbbbs(j) = (eq_R(1) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = eq_Z(j - jmin + 2) / delta_Z
    end do

    do j = 1, nbbbs-1
       if ((rbbbs(j) == rbbbs(j+1)) .and. (zbbbs(j) == zbbbs(j+1))) then
          write(*,*) 'duplicates: ',rbbbs(j),zbbbs(j)
       end if
    end do

    ! get r_boundary(theta)
    self%thetab = atan2 ((zbbbs - self%Z_mag), (rbbbs - self%R_mag))
    self%r_bound = sqrt( (rbbbs - self%R_mag)**2 + (zbbbs - self%Z_mag)**2 )

    call sort(self%thetab, self%r_bound, zbbbs, rbbbs)

    ! Allow for duplicated points near +- pi:
    if (self%thetab(1) == self%thetab(2)) then
       self%thetab(1) = self%thetab(1) + twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    if (self%thetab(nbbbs-1) == self%thetab(nbbbs)) then
       self%thetab(nbbbs) = self%thetab(nbbbs) - twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    ! It isn't likely that a duplicate point would exist near theta = 0,
    ! so I am not allowing this possibility for now.
    do i = 1, nbbbs-1
       if (self%thetab(i) == self%thetab(i+1)) then
          write(*,*) 'Duplicates near theta = 0 not allowed.'
          write(*,*) i, i+1, ' Stopping.'
          call mp_abort('Duplicates near theta = 0 not allowed.')
       end if
    end do
    deallocate (rbbbs, zbbbs)

    self%B_psi = 0.0 ; self%diam = 0.0 ; self%rc = 0.0
    self%fp = 0. ; self%qsf = 0.

    if (lref /= 0.) then
       self%aminor = lref
       write(6,*) 'aminor = ', self%aminor, lref
    else
       self%aminor = self%R_mag
    end if

    self%r_bound = self%r_bound / self%aminor

    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor

    eq_R = eq_R / self%aminor
    eq_Z = eq_Z / self%aminor

    do i = 1, self%nw
       self%R_psi(i, :) = eq_R(i)
       self%Z_psi(i, :) = eq_Z
    end do

    if (bref /= 0.) then
       self%B_T = bref
       write(6,*) 'B_T0 = ', self%B_T, bref
    else
       self%B_T = abs(bcentr)
    end if
    self%psi_a = self%psi_a / (self%B_T * self%aminor**2)
    self%psi_0 = self%psi_0 / (self%B_T * self%aminor**2)
    p_0 = self%pressure(1)

    ! Pointless as fp == 0 here
    self%fp = self%fp / (self%B_T * self%aminor)
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / p_0

    self%eqpsi_2d = self%eqpsi_2d / (self%B_T * self%aminor**2)
    do i = 1, self%nw
       do j = 1, self%nh
          if (eq_Z(j) == self%Z_mag .and. eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (eq_Z(j) - self%Z_mag), (eq_R(i) - self%R_mag))
          end if
       end do
    end do

    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0

    psi_0_out = self%psi_0 ; psi_a_out = self%psi_a
    rmaj = self%R_mag ; B_T0 = self%B_T  ! should rmaj be R_mag? use R_mag for now.
    avgrmid = self%aminor ; nthg = -1

    self%initialised = .true.
    self%has_full_theta_range = .true.
2020 format (5e16.8)
  end subroutine dfitin

  pure subroutine alloc_special_arrays(self)
    class(deq_type), intent(in out) :: self
    allocate (self%rho_mid(self%nrho), self%psi_mid(self%nrho))
  end subroutine alloc_special_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(deq_type), intent(in out) :: self
    deallocate(self%rho_mid, self%psi_mid, self%thetab, self%r_bound)
  end subroutine dealloc_special_arrays

  subroutine setup_special_splines(self)
    use splines, only: new_spline, new_periodic_spline
    use constants, only: twopi
    implicit none
    class(deq_type), intent(in out) :: self
    call new_periodic_spline(size(self%thetab), self%thetab, self%r_bound, twopi, &
         self%rbound_spl, drop_last_point = .true.)
    call new_spline(self%nrho, self%psi_mid, self%rho_mid, self%rho_spl)
  end subroutine setup_special_splines

  subroutine delete_special_splines(self)
    use splines, only: delete_spline, delete_periodic_spline
    implicit none
    class(deq_type), intent(in out) :: self
    call delete_periodic_spline(self%rbound_spl)
    call delete_spline(self%rho_spl)
  end subroutine delete_special_splines

  !> FIXME : Add documentation   | Not implemented
  real function rcenter (self, rp)
    implicit none
    class(deq_type), intent(in) :: self
    real, intent (in) :: rp
    rcenter = self%R_mag !Historic fallback
  end function rcenter

  !> FIXME : Add documentation
  real function rhofun (self, pbar)
    use splines, only: splint
    implicit none
    class(deq_type), intent(in) :: self
    real, intent(in) :: pbar
    rhofun = splint(pbar, self%rho_spl)
  end function rhofun

end module deq
