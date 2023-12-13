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
  use geo_utils, only: abstract_eqfile_cart_geo_type, geo_input_type

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
    ! It may be clearer if we use our spline/periodic_spline wrappers?
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    use constants, only: pi, twopi
    use gs2d, only: read_gs2d, psi, b0, f, ippsi, nr, nz, p, ps
    use gs2d, only: psip, psmin, q, r0, rgrid, rmag, rsep, zgrid, zmag, zsep
    use geo_utils, only: sort, geo_input_type
    use mp, only: mp_abort
    implicit none
    class(gs2d_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real :: psi_N, f_N, rwid, rleft, zhei, rcentr, bcentr
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn
    real, dimension(:), allocatable :: spsi_bar, sefit_R, sefit_Z, rbbbs, zbbbs, eq_R, eq_Z
    real, dimension(:, :), allocatable :: sefit_psi
    real :: zxy11, zxym1, zxy1n, zxymn
    integer :: i, j, ierr, nw_in, nh_in, nbbbs, big
    logical, parameter :: debug =.false.
    self%filename = trim(adjustl(inputs%eqfile))
    big = inputs%big
    if (debug) write(6,*) "gs2din: filename=",self%filename
    ! read GS2D datafile
    call read_gs2d(self%filename)
    if (debug) write(6,*) "gs2din: read_gs2d done, "
    if (debug) write(6,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi

    nw_in = nr ; nh_in = nz
    self%nw = nw_in * big ; self%nr = self%nw
    self%nh = nh_in * big ; self%nt = self%nh
    call self%alloc_arrays()
    allocate(spsi_bar(nw_in), sefit_R(nw_in), sefit_Z(nh_in))
    allocate(sefit_psi(nw_in, nh_in))
    allocate(eq_R(self%nw), eq_Z(self%nh))
    rwid = maxval(rgrid) - minval(rgrid) ; zhei = maxval(zgrid) - minval(zgrid)
    rcentr = r0 ; rleft = minval(rgrid)
    self%R_mag = rmag ; self%Z_mag = zmag
    self%psi_0 = psmin ; self%psi_a = psip ; bcentr = b0

    self%psi_0 = self%psi_0 / twopi
    self%psi_a = self%psi_a / twopi
    if (debug) write(6,*) "gs2din: psi_0, psi_a=", self%psi_0, self%psi_a

    ! pbar is defined by
    ! pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    sefit_R = rgrid ; sefit_Z = zgrid
    if (size(ps) /= nw_in) then
       write(6,*) 'gs2din: size(ps) (',size(ps),') /= nw (',nw_in,')'
       write(6,*) 'gs2din: => should fix the GS2D output file and try again'
       call mp_abort('gs2din: => should fix the GS2D output file and try again')
    end if
    spsi_bar = (ps - minval(ps)) / (maxval(ps) - minval(ps))

    do i = 1, self%nw
       self%psi_bar(i) = float(i-1) / float(self%nw-1)
    end do

    call inter_cspl(nw_in, spsi_bar, f, self%nw, self%psi_bar, self%fp)
    call inter_cspl(nw_in, spsi_bar, p, self%nw, self%psi_bar, self%pressure)
    call inter_cspl(nw_in, spsi_bar, q, self%nw, self%psi_bar, self%qsf)
    ! divide GS2D psi by 2pi to get poloidal flux in (Wb/rad)
    sefit_psi = psi / twopi

    nbbbs = size(rsep)
    allocate(self%thetab(nbbbs), self%r_bound(nbbbs))
    rbbbs = rsep ; zbbbs = zsep
    !Not quite identical to efit unless minval(sefit_Z)==0
    do j = 1, self%nh ! nb Zgrid is not necessarily UpDown symmetric about zero
       eq_Z(j) = ((float(j-1) / float(self%nh-1) )) * zhei + minval(sefit_Z)
    end do

    ! Identical to efit_read_and_set (except prints) from here
    do i = 1, self%nw
       eq_R(i) = rleft + rwid * float(i-1) / float(self%nw-1)
    end do

    allocate(zp(3*nw_in*nh_in), temp(nw_in+2*nh_in))
    allocate(zx1(nh_in), zxm(nh_in), zy1(nw_in), zyn(nw_in))

    ! Setup 2D interpolation of psi(R, Z)
    call fitp_surf1(nw_in, nh_in, sefit_R, sefit_Z, sefit_psi, &
         nw_in, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, self%nh
       do i = 1, self%nw
          self%eqpsi_2d(i,j) = fitp_surf2(eq_R(i), eq_Z(j), nw_in, nh_in, &
               sefit_R, sefit_Z, sefit_psi, nw_in, zp, 1.)
       end do
    end do

    deallocate(zp, temp)
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
       if (self%thetab(i+1) == self%thetab(i)) then
          ! put in kluge for duplicate points, which happens near theta=0:
          self%thetab(i+1) = self%thetab(i+1) + 1.e-8
       end if
    end do

    if (debug) write(6,*) "gs2din:  rbbbs=",rbbbs
    if (debug) write(6,*) "gs2din:  zbbbs=",zbbbs
    call a_minor(rbbbs, zbbbs, self%Z_mag, self%aminor)
    if (debug) write(6,*) "gs2din:  aminor=", self%aminor
    deallocate (rbbbs, zbbbs)

    ! Normalise
    self%B_psi = 0.0 ; self%diam = 0.0 ; self%rc = 0.0
    self%r_bound = self%r_bound / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    eq_R = eq_R / self%aminor ; eq_Z = eq_Z / self%aminor

    do i = 1, self%nw
       self%R_psi(i, :) = eq_R(i)
       self%Z_psi(i, :) = eq_Z
    end do

    self%B_T = abs(bcentr)
    psi_N = self%B_T * self%aminor**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%eqpsi_2d = self%eqpsi_2d / psi_N
    f_N = self%B_T * self%aminor
    self%fp = self%fp / f_N
    if (debug) write(6,*) "gs2din: B_T0, aminor, psi_0, psi_a=", self%B_T, self%aminor, &
         self%psi_0, self%psi_a
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    do i = 1, self%nw
       do j = 1, self%nh
          if (eq_Z(j) == self%Z_mag .and. eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (eq_Z(j) - self%Z_mag), &
                  (eq_R(i) - self%R_mag))
          end if
       end do
    end do
    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0
    self%has_full_theta_range = .true.
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
    use constants, only: pi, twopi
    ! It may be clearer if we use our spline/periodic_spline wrappers?
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    use geo_utils, only: sort
    implicit none
    class(eeq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real :: xdum, psi_N, f_N, rwid, rleft, zhei, rcentr, bcentr, zmid
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn, dummy, eq_R, eq_Z
    real, dimension(:), allocatable :: spsi_bar, sefit_R, sefit_Z, rbbbs, zbbbs
    real, dimension(:, :), allocatable :: sefit_psi
    real :: zxy11, zxym1, zxy1n, zxymn
    character(10) ::  char
    integer :: i, j, ndum, ierr, nw_in, nh_in, funit, nbbbs, big
    self%filename = trim(adjustl(inputs%eqfile))
    big = inputs%big
    open(newunit = funit,  file = self%filename, status='old',form='formatted')
    read(funit, 1000) char, char, char, char, char, i, nw_in, nh_in
    self%nw = nw_in * big ; self%nr = self%nw
    self%nh = nh_in * big ; self%nt = self%nh

    call self%alloc_arrays()
    allocate(spsi_bar(nw_in), sefit_R(nw_in), sefit_Z(nh_in), dummy(nw_in))
    allocate(sefit_psi(nw_in, nh_in))
    allocate(eq_R(self%nw), eq_Z(self%nh))

    read(funit, 2020) rwid, zhei, rcentr, rleft, zmid
    read(funit, 2020) self%R_mag, self%Z_mag, self%psi_0, self%psi_a, bcentr

    ! pbar is defined by pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    do i = 1, 2
       read(funit, 2020) xdum,xdum,xdum,xdum,xdum
    end do

    do i = 1, nw_in
       spsi_bar(i) = float(i-1) / float(nw_in-1)
       sefit_R(i) = rleft + rwid * float(i-1) / float(nw_in-1)
    end do

    do j = 1, nh_in
       sefit_Z(j) = ((float(j-1) / float(nh_in-1) ) - 0.5) * zhei
    end do

    do i = 1, self%nw
       self%psi_bar(i) = float(i-1) / float(self%nw-1)
    end do

    read(funit, 2020) (dummy(j),   j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%fp)
    read(funit, 2020) (dummy(j), j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%pressure)
    read(funit, 2020) (dummy(j),     j = 1, nw_in)
    read(funit, 2020) (dummy(j),     j = 1 ,nw_in)
    read(funit, 2020) ((sefit_psi(i,j) , i = 1, nw_in) , j = 1, nh_in)
    read(funit, 2020) (dummy(j) ,   j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%qsf)

    read(funit, 2022) nbbbs, ndum
    allocate(rbbbs(nbbbs), zbbbs(nbbbs))
    allocate(self%thetab(nbbbs), self%r_bound(nbbbs))
    read(funit, 2020) (rbbbs(i), zbbbs(i) , i = 1, nbbbs)
    close(funit)
    !Not quite identical to gs2d unless minval(sefit_Z)==0 - Should we add zmid?
    do j = 1, self%nh
       eq_Z(j) = ((float(j-1) / float(self%nh-1) ) - 0.5) * zhei
    end do

    ! Identical to gs2d_read_and_set (except prints) from here
    do i = 1, self%nw
       eq_R(i) = rleft + rwid * float(i-1) / float(self%nw-1)
    end do

    allocate(zp(3*nw_in*nh_in), temp(nw_in+2*nh_in))
    allocate(zx1(nh_in), zxm(nh_in), zy1(nw_in), zyn(nw_in))

    ! Setup 2D interpolation of psi(R, Z)
    call fitp_surf1(nw_in, nh_in, sefit_R, sefit_Z, sefit_psi, &
         nw_in, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, self%nh
       do i = 1, self%nw
          self%eqpsi_2d(i,j) = fitp_surf2(eq_R(i), eq_Z(j), nw_in, nh_in, &
               sefit_R, sefit_Z, sefit_psi, nw_in, zp, 1.)
       end do
    end do

    deallocate(zp, temp)
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
       if (self%thetab(i+1) == self%thetab(i)) then
          ! put in kluge for duplicate points, which happens near theta=0:
          self%thetab(i+1) = self%thetab(i+1) +  1.e-8
       end if
    end do

    call a_minor(rbbbs, zbbbs, self%Z_mag, self%aminor)

    deallocate (rbbbs, zbbbs)

    self%B_psi = 0.0 ; self%diam = 0.0 ; self%rc = 0.0
    self%r_bound = self%r_bound / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    eq_R = eq_R / self%aminor ; eq_Z = eq_Z / self%aminor

    do i = 1, self%nw
       self%R_psi(i, :) = eq_R(i)
       self%Z_psi(i, :) = eq_Z
    end do

    self%B_T = abs(bcentr)
    psi_N = self%B_T * self%aminor**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%eqpsi_2d = self%eqpsi_2d / psi_N
    f_N = self%B_T * self%aminor
    self%fp = self%fp / f_N
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    do i = 1, self%nw
       do j = 1, self%nh
          if (eq_Z(j) == self%Z_mag .and. eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (eq_Z(j) - self%Z_mag), &
                  (eq_R(i) - self%R_mag))
          end if
       end do
    end do
    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0
    self%has_full_theta_range = .true.

    if (.true.) then
       write (*,*) "Finished efitin... imported EFIT equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", self%aminor
       write (*,*) 'R_mag', self%R_mag
       write (*,*) 'B_T0', self%B_T
       write (*,*) 'beta', self%beta_0
    end if
1000 format(5(a10),i2,i4,i4)
2020 format (5e16.9)
2022 format (2i5)
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

  !> FIXME : Add documentation
  subroutine a_minor(r, z, Z_mag, a)
    use mp, only: mp_abort
    use splines, only: new_spline, splint, delete_spline, spline
    implicit none
    real, dimension(:), intent (in) :: r, z
    real, intent(in) :: Z_mag
    real, intent(out) :: a
    real :: r1, r2
    integer, parameter :: nz = 5, half_nz = int(nz/2.0), invalidPoint = -nz
    real, dimension(nz) :: rtmp, ztmp
    integer :: i, j, i1, n, k, index
    type (spline) :: spl
    !CMR, 28/10/08: add code to avoid duplicate points in 5 point spline
    !               to determine r2 on inboard mid-plane
    logical, parameter :: debug=.false.

    k = 0
    n = size(r)

    if (debug) write(6,*) "aminor:"
    if (debug) write(6,fmt='("r=",10f8.4)') r
    if (debug) write(6,fmt='("z=",10f8.4)') z

    if (n < nz) then
       write(*,*) 'nbbbs < nz -- very strange.  Stopping.'
       write(*,*) 'Look in eeq.f90.'
!      call mp_abort('nbbbs < nz -- very strange.  Stopping : Look in eeq.f90')
    end if

    j = 0
    do i = half_nz + 1, 1, -1
       j = j + 1
       ztmp(j) = z(i)
       rtmp(j) = r(i)
    end do

    if (r(1) == r(n) .and. z(1) == z(n)) then
       do i = n-1, n-half_nz, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       end do
    else
       do i = n, n-half_nz+1, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       end do
    end if

    if (debug) write(6,fmt='("rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("Z_mag=",f8.4)') Z_mag

    call new_spline(nz, ztmp, rtmp, spl)
    r1 = splint(Z_mag, spl)
    call delete_spline(spl)

    ! find point near magnetic axis elevation on low field side

    ! Initialise i1 to some invalid value so we can check if an appropriate
    ! point is not found.
    i1 = invalidPoint

    ! Note there appears to be an assumption here about the order
    ! of z which may not be true
    do i = nz, n
       if (z(i) - Z_mag > 0.) then
          i1 = i - 1
          exit
       end if
    end do

    if (i1 == invalidPoint) then
       write(*, '("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")')
       call mp_abort("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")
    end if

    !CMR, 28/10/08: modify this code to avoid duplicate points in 5 point spline
    !               to determine r2 on inboard mid-plane
    do i = 1, nz
       index = i1 - half_nz + i - 1 + k
       rtmp(i) = r(index)
       ztmp(i) = z(index)
       if ( i > 1 ) then
          ! If we find a duplicate point then we increase the
          ! duplicate point counter in k, used to offset index and
          ! store the next value instead Note, here we assume we can't
          ! get three duplicate values. A better approach may be to
          ! instead use a while loop where we keep going until we've
          ! found enough points.
          if ((rtmp(i)-rtmp(i-1))**2+(ztmp(i)-ztmp(i-1))**2 < 1.0e-7) then
             k = k + 1
             if (debug) write(6,fmt='("a_minor: duplicate pt, set k=",i2)') k
             rtmp(i) = r(index + 1)
             ztmp(i) = z(index + 1)
          end if
       end if
    end do

    if (debug) write(6,fmt='("a_minor: rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("a_minor: ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("a_minor: Z_mag=",f8.4)') Z_mag

    call new_spline(nz, ztmp, rtmp, spl)
    r2 = splint(Z_mag, spl)
    call delete_spline(spl)

    a = (r2 - r1) / 2.
    if (debug) write(6,*) "a_minor: r1,r2=",r1,r2
    if (debug) write(6,*) "a_minor: a=",a
  end subroutine a_minor
end module eeq
