!> This module contains routines to read from CMR's GS2D equilibrium solver
module gs2d
  private

  public :: write_gs2d, read_gs2d, psi, b0, f, ippsi, nr, nz, p, ps
  public :: psip, psmin, q, r0, rgrid, rmag, rsep, zgrid, zmag, zsep

  real, dimension(:), allocatable :: ps, amin_gs2d, q, f, p, pp
  real, dimension(:), allocatable :: rsep, zsep, rgrid, zgrid
  real, dimension(:,:), allocatable :: psi
  real :: r0, a, rmag, zmag, psmin, psip, b0, ippsi
  integer :: nsurf, nsep, nr, nz

contains

  !> FIXME : Add documentation
  subroutine read_gs2d(filename)
    implicit none
    character(len = *), intent(in) :: filename
    character(len = 80) :: line
    integer :: i, unit
    open(newunit = unit, file = trim(filename) , status = 'unknown')
    do i = 1, 7
       read(unit, fmt='(a80)') line
    end do
    read(unit, *) r0,a,rmag,zmag
    read(unit, fmt='(a80)') line ; read(unit, *) psmin, psip, b0, ippsi
    read(unit, fmt='(a80)') line ; read(unit, *) nsurf
    allocate(ps(nsurf),amin_gs2d(nsurf),q(nsurf),f(nsurf),p(nsurf),pp(nsurf))
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') ps
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') amin_gs2d
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') q
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') f
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') p
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') pp
    read(unit, fmt='(T22,I6)') nsep
    allocate(rsep(nsep),zsep(nsep))
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') rsep
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') zsep
    read(unit, fmt='(a80)') line ; read(unit, *) nr,nz
    allocate(rgrid(nr), zgrid(nz), psi(nr, nz))
    read(unit, fmt='(a80)') line ; read(unit, *) rgrid
    read(unit, fmt='(a80)') line ; read(unit, *) zgrid
    read(unit, fmt='(a80)') line ; read(unit, fmt='(1p,8e16.8)') psi
    close(unit)
  end subroutine read_gs2d

  !> FIXME : Add documentation
  subroutine write_gs2d(filename)
    implicit none
    character(len = *), intent(in) :: filename
    integer :: unit
    open(newunit = unit, file = filename, status = 'unknown')
    write(unit, fmt='("GS2 input file",T30,"Produced by GS2D at:",a40)') ' '
    write(unit, fmt='(A80/,"GS2D Equilibrium Boundary description:",3(/A80))') repeat('-',80),' ',' ',repeat('-',80)
    write(unit, fmt='(T2,"r0",T15,"a",T27,"rmag",T39,"zmag (m)"/1p,4e16.8)') r0,a,rmag,zmag
    write(unit, fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi
    write(unit, fmt='("nfs"/I6)') nsurf
    write(unit, fmt='("Psi on 1d grid (for FS quantities) (T)")')
    write(unit, fmt='(1p,8e16.8)') ps
    write(unit, fmt='("amin_gs2d (m)")')
    write(unit, fmt='(1p,8e16.8)') amin_gs2d
    write(unit, fmt='("q")')
    write(unit, fmt='(1p,8e16.8)') q
    write(unit, fmt='("f =r B_phi (Tm)")')
    write(unit, fmt='(1p,8e16.8)') f
    write(unit, fmt='("p (Pa)")')
    write(unit, fmt='(1p,8e16.8)') p
    write(unit, fmt='("dp/dpsi (Pa/Wb)")')
    write(unit, fmt='(1p,8e16.8)') pp
    write(unit, fmt='("No of points on LCFS=",I6)') nsep
    write(unit, fmt='("r(j) (m) on LCFS")')
    write(unit, fmt='(1p,8e16.8)') rsep
    write(unit, fmt='("z(j) (m) on LCFS")')
    write(unit, fmt='(1p,8e16.8)') zsep
    write(unit, fmt='("NR",T14,"NZ"/2I6)') NR, NZ
    write(unit, fmt='("rgrid (m)")')
    write(unit, fmt='(1p,8e16.8)') rgrid
    write(unit, fmt='("zgrid (m)")')
    write(unit, fmt='(1p,8e16.8)') zgrid
    write(unit, fmt='("Psi on grid (Wb) : NB Br=(1/2pi r)*dpsi/dz")')
    write(unit, fmt='(1p,8e16.8)') psi
    close(unit)
  end subroutine write_gs2d
end module gs2d

!> This module is a submodule of geometry which handles reading from the
!> ascii EQDSK format output by EFIT, but also now by other codes. This
!> file contains psi on an R,Z grid, as well as other quantities such as
!> q (safety factor), I (a.k.a. f) and p (pressure) on a psi grid.
!>
!> The normalising field is set to the field on axis; the normalising length
!> is set to the half-diameter of the LCFS.
module eeq
  use geo_utils, only: abstract_eqfile_cart_geo_type, geo_input_type, sort

  implicit none

  private

  public :: eeq_type, gs2d_type

  type, abstract, extends(abstract_eqfile_cart_geo_type) :: eeq_base_type
     private
   contains
     procedure :: setup_special_splines, delete_special_splines
     procedure :: rcenter, dealloc_special_arrays
  end type eeq_base_type

  type, extends(eeq_base_type) :: eeq_type
   contains
     procedure :: read_and_set => efit_read_and_set
  end type eeq_type

  type, extends(eeq_base_type) :: gs2d_type
   contains
     procedure :: read_and_set => gs2d_read_and_set
  end type gs2d_type

contains

  !> This subroutine reads a GS2D output file containing
  !> the axisymmetric magnetic field geometry on a rectangular
  !> domain defined by the coordinates (R,Z).
  !>
  !> eq_R     R grid
  !> eq_Z     Z grid
  !> fp    F on psibar grid
  !> eqpsi_2d   psi on (R,Z) grid
  !> R_mag        major radius of the magnetic axis
  !> Z_mag        elevation of the magnetic axis
  !> rwid       total width of the domain
  !> rleft      position of leftmost point of domain
  !> zhei       total height of domain
  subroutine gs2d_read_and_set(self, inputs)
    use constants, only: twopi
    use gs2d, only: read_gs2d, psi, b0, f, ippsi, nr, nz, p, ps
    use gs2d, only: psip, psmin, q, r0, rgrid, rmag, rsep, zgrid, zmag, zsep
    use mp, only: mp_abort
    implicit none
    class(gs2d_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real :: rwid, rleft, zhei, rcentr
    real, dimension(:), allocatable :: spsi_bar
    real, dimension(:, :), allocatable :: sefit_psi
    integer :: nw_in, nh_in
    logical, parameter :: debug =.false.
    self%type_name = 'GS2DEQ'
    self%filename = trim(adjustl(inputs%eqfile))
    if (debug) write(6,*) "gs2din: filename=",self%filename
    ! read GS2D datafile
    call read_gs2d(self%filename)
    if (debug) write(6,*) "gs2din: read_gs2d done, "
    if (debug) write(6,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi

    nw_in = nr ; nh_in = nz
    if (size(ps) /= nw_in) then
       write(6,*) 'gs2din: size(ps) (',size(ps),') /= nw (',nw_in,')'
       write(6,*) 'gs2din: => should fix the GS2D output file and try again'
       call mp_abort('gs2din: => should fix the GS2D output file and try again')
    end if
    allocate(spsi_bar(nw_in), sefit_psi(nw_in, nh_in))
    rwid = maxval(rgrid) - minval(rgrid) ; zhei = maxval(zgrid) - minval(zgrid)
    rcentr = r0 ; rleft = minval(rgrid)
    self%R_mag = rmag ; self%Z_mag = zmag
    self%psi_0 = psmin ; self%psi_a = psip
    self%psi_0 = self%psi_0 / twopi ; self%psi_a = self%psi_a / twopi
    if (debug) write(6,*) "gs2din: psi_0, psi_a=", self%psi_0, self%psi_a
    ! divide GS2D psi by 2pi to get poloidal flux in (Wb/rad)
    sefit_psi = psi / twopi

    ! pbar is defined by pbar = (psi-psi_0)/(psi_a-psi_0). fp and q are
    ! functions of pbar
    spsi_bar = (ps - minval(ps)) / (maxval(ps) - minval(ps))
    call self%shared_setup(nw_in, nh_in, spsi_bar, rgrid, zgrid, f, p, q, sefit_psi, &
         size(rsep), rsep, zsep, b0, rleft, rwid, minval(zgrid), zhei, inputs%big, .true.)
  end subroutine gs2d_read_and_set

  !> This subroutine reads an EFIT output file containing
  !> the axisymmetric magnetic field geometry on a rectangular
  !> domain defined by the coordinates (R,Z). It reads and stores
  !> the following quantities (among others).
  !>
  !> - eq_R:     R grid
  !> - eq_Z:     Z grid
  !> - fp:    F on psibar grid
  !> - eqpsi_2d:   psi on (R,Z) grid
  !> - R_mag:        major radius of the magnetic axis
  !> - Z_mag:        elevation of the magnetic axis
  !> - rwid:       total width of the domain
  !> - rleft:      position of leftmost point of domain
  !> - zhei:       total height of domain
  subroutine efit_read_and_set(self, inputs)
    implicit none
    class(eeq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real :: xdum, rwid, rleft, zhei, rcentr, bcentr, zmid
    real, dimension(:), allocatable :: dummy, spsi_bar, sefit_R, sefit_Z, rbbbs, zbbbs, f, p, q
    real, dimension(:, :), allocatable :: sefit_psi
    character(10) ::  char
    integer :: i, j, ndum, nw_in, nh_in, funit, nbbbs
    self%type_name = 'EFIT'
    self%filename = trim(adjustl(inputs%eqfile))
    open(newunit = funit,  file = self%filename, status='old',form='formatted')
    read(funit, '(5(a10),i2,i4,i4)') char, char, char, char, char, i, nw_in, nh_in
    allocate(spsi_bar(nw_in), sefit_R(nw_in), sefit_Z(nh_in), dummy(nw_in))
    allocate(sefit_psi(nw_in, nh_in), f(nw_in), p(nw_in), q(nw_in))

    read(funit, 2020) rwid, zhei, rcentr, rleft, zmid
    read(funit, 2020) self%R_mag, self%Z_mag, self%psi_0, self%psi_a, bcentr
    read(funit, 2020) xdum,xdum,xdum,xdum,xdum
    read(funit, 2020) xdum,xdum,xdum,xdum,xdum

    read(funit, 2020) (f(j), j = 1, nw_in)
    read(funit, 2020) (p(j), j = 1, nw_in)
    read(funit, 2020) (dummy(j), j = 1, nw_in)
    read(funit, 2020) (dummy(j), j = 1 ,nw_in)
    read(funit, 2020) ((sefit_psi(i,j) , i = 1, nw_in) , j = 1, nh_in)
    read(funit, 2020) (q(j), j = 1, nw_in)

    read(funit, '(2i5)') nbbbs, ndum
    allocate(rbbbs(nbbbs), zbbbs(nbbbs))
    read(funit, 2020) (rbbbs(i), zbbbs(i) , i = 1, nbbbs)
    close(funit)

    ! pbar is defined by pbar = (psi-psi_0)/(psi_a-psi_0) fp and q are
    ! functions of pbar
    do i = 1, nw_in
       spsi_bar(i) = float(i - 1) / float(nw_in - 1)
       sefit_R(i) = rleft + rwid * float(i - 1) / float(nw_in - 1)
    end do

    do j = 1, nh_in
       sefit_Z(j) = ((float(j - 1) / float(nh_in - 1) ) - 0.5) * zhei !Should we add zmid?
    end do

    call self%shared_setup(nw_in, nh_in, spsi_bar, sefit_R, sefit_Z, f, p, q, sefit_psi, &
         nbbbs, rbbbs, zbbbs, bcentr, rleft, rwid, -0.5*zhei, zhei, inputs%big, .true.) !Should -0.5*zhei be zmid
2020 format (5e16.9)
  end subroutine efit_read_and_set

  !> FIXME : Add documentation
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(eeq_base_type), intent(in out) :: self
    deallocate(self%thetab, self%r_bound)
  end subroutine dealloc_special_arrays

  subroutine setup_special_splines(self)
    use splines, only: new_periodic_spline
    use constants, only: twopi
    implicit none
    class(eeq_base_type), intent(in out) :: self
    call new_periodic_spline(size(self%thetab), self%thetab, self%r_bound, twopi, &
         self%rbound_spl, drop_last_point = .true.)
  end subroutine setup_special_splines

  subroutine delete_special_splines(self)
    use splines, only: delete_periodic_spline
    implicit none
    class(eeq_base_type), intent(in out) :: self
    call delete_periodic_spline(self%rbound_spl)
  end subroutine delete_special_splines

  !> FIXME : Add documentation
  real function rcenter (self, rp)
    use constants, only: pi
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, intent (in) :: rp
    rcenter = self%R_mag + 0.5 * (self%rfun(rp, 0.) - self%rfun(rp, pi))
  end function rcenter
end module eeq
