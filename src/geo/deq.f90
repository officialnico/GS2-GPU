!> FIXME : Add documentation
module deq
  use geo_utils, only: abstract_eqfile_cart_geo_type, geo_input_type
  use splines, only: spline
  implicit none

  private
  public :: deq_type

  type, extends(abstract_eqfile_cart_geo_type) :: deq_type
     private
     real, allocatable, dimension(:) :: rho_mid, psi_mid
     type(spline) :: rho_spl
   contains
     procedure :: read_and_set, setup_special_splines, delete_special_splines
     procedure :: rcenter, rhofun, dealloc_special_arrays
  end type deq_type

contains

  !> This subroutine reads an DFIT output file containing
  !> the axisymmetric magnetic field geometry on a rectangular
  !> domain defined by the coordinates (R,Z).
  subroutine read_and_set(self, inputs)
    use mp, only: mp_abort
    use splines, only: inter_cspl
    use geo_utils, only: sort
    implicit none
    class(deq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, dimension(:), allocatable :: spsi_bar, sdfit_R, sdfit_Z, rbbbs, zbbbs, eq_R, eq_Z, f, p, q
    real, dimension(:, :), allocatable :: sdfit_psi
    real :: psi_N, f_N, rwid, zhei, bcentr, delta_R, delta_Z, lref, bref
    integer :: i, j, jmin, jmax, nw_in, nh_in, funit, nbbbs, nrho
    logical :: eqnorm_exists
    self%type_name = 'deq'
    ! Need to generalize initialization condition if equilibrium changes
    lref = 0. ; bref = 0.
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
    read(funit, *) nw_in, nh_in, nrho
    allocate (self%rho_mid(nrho), self%psi_mid(nrho))
    allocate(spsi_bar(nw_in), sdfit_R(nw_in), sdfit_Z(nh_in), sdfit_psi(nw_in, nh_in))
    allocate(eq_R(self%nw), eq_Z(self%nh), f(nw_in), p(nw_in), q(nw_in))
    read(funit, 2020) rwid, zhei
    read(funit, 2020) self%R_mag, self%Z_mag
    read(funit, 2020) self%psi_0, self%psi_a, bcentr
    f = 0 ; q = 0
    read(funit, 2020) (p(j), j = 1, nw_in)  ! pressure read
    read(funit, 2020) ((sdfit_psi(i,j) , i = 1, nw_in) , j = 1, nh_in)

    ! Only needed for irho = 4
    read (funit, 2020) (self%rho_mid(i), i = 1, nrho)
    read (funit, 2020) (self%psi_mid(i), i = 1, nrho)
    close(funit)

    ! pbar is defined by pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    do i = 1, nw_in
       spsi_bar(i) = float(i - 1) / float(nw_in - 1)
       sdfit_R(i) = rwid * float(i - 1) / float(nw_in - 1)
    end do

    do j = 1, nh_in
       sdfit_Z(j) = (float(j - 1) / float(nh_in - 1) - 0.5) * zhei
    end do

    ! stay away from edge of box
    delta_R = 1.02 ; delta_Z = 1.02

    nbbbs = 2 * (nw_in + nh_in) - 5
    allocate(rbbbs(nbbbs), zbbbs(nbbbs))
    jmin = 1 ; jmax = nh_in / 2
    do j = jmin, jmax
       rbbbs(j) = (sdfit_R(1) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = sdfit_Z(j + nh_in / 2) / delta_Z
    end do

    jmin = jmax + 1 ; jmax = jmax + nw_in - 2
    do j = jmin, jmax
       rbbbs(j) = (sdfit_R(j - jmin + 2) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = sdfit_Z(nh_in) / delta_Z
    end do

    jmin = jmax + 1 ; jmax = jmax + nh_in - 1
    do j = jmin, jmax
       rbbbs(j) = (sdfit_R(nw_in) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = sdfit_Z(jmax - j + 2) / delta_Z
    end do

    jmin = jmax + 1 ; jmax = jmax + nw_in - 1
    do j = jmin, jmax
       rbbbs(j) = (sdfit_R(jmax - j + 2) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = sdfit_Z(1) / delta_Z
    end do

    jmin = jmax + 1 ; jmax = jmax + nh_in / 2
    do j = jmin, jmax
       rbbbs(j) = (sdfit_R(1) - 0.5 * rwid) / delta_R + 0.5 * rwid
       zbbbs(j) = sdfit_Z(j - jmin + 2) / delta_Z
    end do

    do j = 1, nbbbs-1
       if ((rbbbs(j) == rbbbs(j+1)) .and. (zbbbs(j) == zbbbs(j+1))) then
          write(*,*) 'duplicates: ',rbbbs(j),zbbbs(j)
       end if
    end do

    ! Overwrite inputs with normalisations if provided
    if (bref /= 0.) then
       bcentr = bref
       write(6,*) 'B_T0 = ', bref
    end if

    if (lref /= 0.) then
       self%aminor = lref
       write(6,*) 'aminor = ', self%aminor, lref
    else
       self%aminor = self%R_mag
    end if

    call self%shared_setup(nw_in, nh_in, spsi_bar, sdfit_R, sdfit_Z, f, p, q, sdfit_psi, &
         nbbbs, rbbbs, zbbbs, bcentr, 0.0, rwid, -0.5*zhei, zhei, 1, .false.)
2020 format (5e16.8)
  end subroutine read_and_set

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
    call new_spline(size(self%psi_mid), self%psi_mid, self%rho_mid, self%rho_spl)
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
