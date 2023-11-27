!> This module contains routines to read from CMR's GS2D equilibrium solver
module gs2d
  private

  public :: write_gs2d, read_gs2d, psi, b0, f, ippsi, nr, nz, p, ps
  public :: psip, psmin, q, r0, rgrid, rmag, rsep, zgrid, zmag, zsep

  real, dimension(:), allocatable :: ps,amin_gs2d,q,f,p,pp
  real, dimension(:), allocatable :: rsep,zsep,rgrid,zgrid
  real, dimension(:,:), allocatable :: psi
  real :: r0,a,rmag,zmag
  real :: psmin,psip,b0,ippsi
  integer :: nsurf,nsep,nr,nz

contains

  !> FIXME : Add documentation
  subroutine read_gs2d(filename)
    implicit none
!CMR, 2/3/2006: lengthen the string filename to 80
    character(len=80), intent(in) :: filename
    character(len=80) :: line
    integer :: i
    open(1,file=trim(filename),status='unknown')
    do i=1,7
       read(1,fmt='(a80)') line
    enddo
    read(1,*) r0,a,rmag,zmag
    read(1,fmt='(a80)') line ; read(1,*) psmin,psip,b0,ippsi
    read(1,fmt='(a80)') line ; read(1,*) nsurf
    allocate(ps(nsurf),amin_gs2d(nsurf),q(nsurf),f(nsurf),p(nsurf),pp(nsurf))
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') ps
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') amin_gs2d
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') q
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') f
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') p
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') pp
    read(1,fmt='(T22,I6)') nsep
    allocate(rsep(nsep),zsep(nsep))
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') rsep
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') zsep
    read(1,fmt='(a80)') line ; read(1,*) nr,nz
    allocate(rgrid(nr),zgrid(nz),psi(nr,nz))
    read(1,fmt='(a80)') line ; read(1,*) rgrid
    read(1,fmt='(a80)') line ; read(1,*) zgrid

    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') psi
    close(1)
  end subroutine read_gs2d

  !> FIXME : Add documentation
  subroutine write_gs2d(filename)
    implicit none
    character(len=40), intent(in) :: filename
    open(1,file=filename,status='unknown')
    write(1,fmt='("GS2 input file",T30,"Produced by GS2D at:",a40)') ' '
    write(1,fmt='(A80/,"GS2D Equilibrium Boundary description:",3(/A80))') repeat('-',80),' ',' ',repeat('-',80)
    write(1,fmt='(T2,"r0",T15,"a",T27,"rmag",T39,"zmag (m)"/1p,4e16.8)') r0,a,rmag,zmag
    write(1,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi
    write(1,fmt='("nfs"/I6)') nsurf
    write(1,fmt='("Psi on 1d grid (for FS quantities) (T)")')
    write(1,fmt='(1p,8e16.8)') ps
    write(1,fmt='("amin_gs2d (m)")')
    write(1,fmt='(1p,8e16.8)') amin_gs2d
    write(1,fmt='("q")')
    write(1,fmt='(1p,8e16.8)') q
    write(1,fmt='("f =r B_phi (Tm)")')
    write(1,fmt='(1p,8e16.8)') f
    write(1,fmt='("p (Pa)")')
    write(1,fmt='(1p,8e16.8)') p
    write(1,fmt='("dp/dpsi (Pa/Wb)")')
    write(1,fmt='(1p,8e16.8)') pp
    write(1,fmt='("No of points on LCFS=",I6)') nsep
    write(1,fmt='("r(j) (m) on LCFS")')
    write(1,fmt='(1p,8e16.8)') rsep
    write(1,fmt='("z(j) (m) on LCFS")')
    write(1,fmt='(1p,8e16.8)') zsep
    write(1,fmt='("NR",T14,"NZ"/2I6)') NR, NZ
    write(1,fmt='("rgrid (m)")')
    write(1,fmt='(1p,8e16.8)') rgrid
    write(1,fmt='("zgrid (m)")')
    write(1,fmt='(1p,8e16.8)') zgrid
    write(1,fmt='("Psi on grid (Wb) : NB Br=(1/2pi r)*dpsi/dz")')
    write(1,fmt='(1p,8e16.8)') psi
    close(1)
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
  use geo_utils, only: abstract_eqfile_cart_geo_type

  implicit none

  private

  public :: eeq_type, gs2d_type
  public :: verbosity  !Should this be replaced by usage of runtime_test:verbosity?

  integer :: verbosity = 0  !Should this be replaced by usage of runtime_test:verbosity?

  type, abstract, extends(abstract_eqfile_cart_geo_type) :: eeq_base_type
     private
   contains
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: rcenter => rcenter
     procedure :: diameter => diameter
     procedure :: dpfun => dpfun
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     procedure :: alloc_special_arrays => alloc_special_arrays
     procedure :: dealloc_special_arrays => dealloc_special_arrays
  end type eeq_base_type

  type, extends(eeq_base_type) :: eeq_type
   contains
     procedure :: initialise => efitin
  end type eeq_type

  type, extends(eeq_base_type) :: gs2d_type
   contains
     procedure :: initialise => gs2din
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
  subroutine gs2din(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
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
    real, intent(out) :: psi_0_out, psi_a_out, B_T0, avgrmid, rmaj
    integer, intent(out) :: nthg
    real :: p_0, rwid, rleft, zhei, rcentr, bcentr
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn
    real, dimension(:), allocatable :: spsi_bar, sefit_R, sefit_Z, rbbbs, zbbbs
    real, dimension(:, :), allocatable :: sefit_psi
    real :: zxy11, zxym1, zxy1n, zxymn
    integer :: i, j, ierr, big, nw_in, nh_in
    logical, parameter :: debug =.false.
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = -1
       return
    end if

    nthg = -1
    big = inputs%big
    self%filename = trim(adjustl(inputs%eqfile))
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
       self%eq_R(i) = rleft +rwid * float(i-1) / float(self%nw-1)
    end do

    ! nb Zgrid is not necessarily UpDown symmetric
    do j = 1, self%nh
       self%eq_Z(j) = ((float(j-1) / float(self%nh-1) ))*zhei+minval(zgrid)
    end do

    call inter_cspl(nw_in, spsi_bar, f, self%nw, self%psi_bar, self%fp)
    call inter_cspl(nw_in, spsi_bar, p, self%nw, self%psi_bar, self%pressure)
    ! divide GS2D psi by 2pi to get poloidal flux in (Wb/rad)
    sefit_psi = psi / twopi

    allocate(zp(3*nw_in*nh_in), temp(nw_in+2*nh_in))
    allocate(zx1(nh_in), zxm(nh_in), zy1(nw_in), zyn(nw_in))

    ! Setup 2D interpolation of psi(R, Z)
    call fitp_surf1(nw_in, nh_in, sefit_R, sefit_Z, sefit_psi, &
         nw_in, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, self%nh
       do i = 1, self%nw
          self%eqpsi_2d(i,j) = fitp_surf2(self%eq_R(i), self%eq_Z(j), nw_in, nh_in, &
               sefit_R, sefit_Z, sefit_psi, nw_in, zp, 1.)
       end do
    end do

    deallocate(zp, temp)

    call inter_cspl(nw_in, spsi_bar, q, self%nw, self%psi_bar, self%qsf)

    self%nbbbs = size(rsep)
    allocate(rbbbs(self%nbbbs), zbbbs(self%nbbbs))
    allocate(self%thetab(self%nbbbs), self%r_bound(self%nbbbs))
    rbbbs = rsep ; zbbbs = zsep

    ! get r_boundary(theta)
    self%thetab = atan2 ((zbbbs - self%Z_mag), (rbbbs - self%R_mag))
    self%r_bound = sqrt( (rbbbs - self%R_mag)**2 + (zbbbs - self%Z_mag)**2 )

    call sort(self%thetab, self%r_bound, zbbbs, rbbbs)

    ! Allow for duplicated points near +- pi:
    if (self%thetab(1) == self%thetab(2)) then
       self%thetab(1) = self%thetab(1) + twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    if (self%thetab(self%nbbbs-1) == self%thetab(self%nbbbs)) then
       self%thetab(self%nbbbs) = self%thetab(self%nbbbs) - twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    ! It isn't likely that a duplicate point would exist near theta = 0,
    ! so I am not allowing this possibility for now.
    do i = 1, self%nbbbs-1
       if (self%thetab(i+1) == self%thetab(i)) then
          ! put in kluge for duplicate points, which happens near theta=0:
          self%thetab(i+1) = self%thetab(i+1)+1.e-8
       end if
    end do

    if (debug) write(6,*) "gs2din:  rbbbs=",rbbbs
    if (debug) write(6,*) "gs2din:  zbbbs=",zbbbs
    call a_minor(rbbbs, zbbbs, self%Z_mag, avgrmid)
    self%aminor = avgrmid
    if (debug) write(6,*) "gs2din:  aminor=", self%aminor
    deallocate (rbbbs, zbbbs)

    ! Normalise
    self%r_bound = self%r_bound / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    self%eq_R = self%eq_R / self%aminor
    self%eq_Z = self%eq_Z / self%aminor

    ! should rmaj be R_mag or rcentr? use R_mag for now.
    rmaj = self%R_mag
    B_T0 = abs(bcentr)
    self%B_T = B_T0
    self%psi_a = self%psi_a / (B_T0*self%aminor**2)
    self%psi_0 = self%psi_0 / (B_T0*self%aminor**2)
    psi_0_out = self%psi_0 ; psi_a_out = self%psi_a
    if (debug) write(6,*) "gs2din: B_T0, aminor, psi_0, psi_a=", self%B_T, self%aminor, &
         self%psi_0, self%psi_a

    p_0 = self%pressure(1)
    self%fp = self%fp / (self%B_T * self%aminor)
    ! MKS: beta = 2 mu_0 p / B**2

    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / p_0
    self%eqpsi_2d = self%eqpsi_2d / (self%B_T * self%aminor**2)
    do i = 1, self%nw
       do j = 1, self%nh
          if (self%eq_Z(j) == self%Z_mag .and. self%eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (self%eq_Z(j) - self%Z_mag), &
                  (self%eq_R(i) - self%R_mag))
          end if
       end do
    end do
    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0
    self%B_psi = 0.0
    self%initialised = .true.
    self%has_full_theta_range = .true.
  end subroutine gs2din

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
  subroutine efitin(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid, nthg)
    use constants, only: pi, twopi
    ! It may be clearer if we use our spline/periodic_spline wrappers?
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    use geo_utils, only: sort, geo_input_type
    implicit none
    class(eeq_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, B_T0, avgrmid, rmaj
    integer, intent(out) :: nthg
    real :: xdum, p_0, rwid, rleft, zhei, rcentr, bcentr
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn, dummy
    real, dimension(:), allocatable :: spsi_bar, sefit_R, sefit_Z, rbbbs, zbbbs
    real, dimension(:, :), allocatable :: sefit_psi
    real :: zxy11, zxym1, zxy1n, zxymn
    character(10) ::  char
    integer :: i, j, init, ndum, ierr, big, nw_in, nh_in
    if (self%initialised) then
       psi_0_out = self%psi_0
       psi_a_out = self%psi_a
       rmaj = self%R_mag
       B_T0 = self%B_T
       avgrmid = self%aminor
       nthg = -1
       return
    end if

    nthg = -1
    big = inputs%big
    ! Need to generalize initialization condition if equilibrium changes
    if (verbosity>2) write (*,*) 'Starting efitin'
    self%filename = trim(adjustl(inputs%eqfile))
    if (verbosity>1) write (*,*) 'Opening file ', self%filename
    open(unit=5, file = self%filename, status='old',form='formatted')

    ! Read the data
    read(5,1000) char, char, char, char, char, i, nw_in, nh_in

    if (verbosity>1) write (*,*) 'EFIT dimensions are: ', nw_in, nh_in

    self%nw = nw_in * big ; self%nr = self%nw
    self%nh = nh_in * big ; self%nt = self%nh

    call self%alloc_arrays()
    allocate(spsi_bar(nw_in), sefit_R(nw_in), sefit_Z(nh_in), dummy(nw_in))
    allocate(sefit_psi(nw_in, nh_in))

    read(5,2020) rwid, zhei, rcentr, rleft, xdum
    read(5,2020) self%R_mag, self%Z_mag, self%psi_0, self%psi_a, bcentr

    ! pbar is defined by
    ! pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    do i = 1, 2
       read(5,2020) xdum,xdum,xdum,xdum,xdum
    end do

    do i = 1, nw_in
       spsi_bar(i) = float(i-1) / float(nw_in-1)
       sefit_R(i) = rleft +rwid * float(i-1) / float(nw_in-1)
    end do

    do j = 1, nh_in
       sefit_Z(j) = ((float(j-1) / float(nh_in-1) ) - 0.5)*zhei
    end do

    do i = 1, self%nw
       self%psi_bar(i) = float(i-1) / float(self%nw-1)
       self%eq_R(i) = rleft +rwid * float(i-1) / float(self%nw-1)
    end do

    do j = 1, self%nh
       self%eq_Z(j) = ((float(j-1) / float(self%nh-1) ) - 0.5)*zhei
    end do

    read(5,2020) (dummy(j),   j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%fp)
    read(5,2020) (dummy(j), j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%pressure)
    read(5,2020) (dummy(j),     j = 1, nw_in)
    read(5,2020) (dummy(j),     j = 1 ,nw_in)
    read(5,2020) ((sefit_psi(i,j) , i = 1, nw_in) , j = 1, nh_in)

    allocate(zp(3*nw_in*nh_in), temp(nw_in+2*nh_in))
    allocate(zx1(nh_in), zxm(nh_in), zy1(nw_in), zyn(nw_in))

    call fitp_surf1(nw_in, nh_in, sefit_R, sefit_Z, sefit_psi, &
         nw_in, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, self%nh
       do i = 1, self%nw
          self%eqpsi_2d(i,j) = fitp_surf2(self%eq_R(i), self%eq_Z(j), nw_in, nh_in, &
               sefit_R, sefit_Z, sefit_psi, nw_in, zp, 1.)
       end do
    end do

    deallocate(zp, temp)

    read(5,2020) (dummy(j) ,   j = 1, nw_in)
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%qsf)

    read(5,2022) self%nbbbs, ndum
    allocate(rbbbs(self%nbbbs), zbbbs(self%nbbbs))
    allocate(self%thetab(self%nbbbs), self%r_bound(self%nbbbs))
    read(5,2020) (rbbbs(i), zbbbs(i) , i = 1, self%nbbbs)
    close(5)

    self%thetab = atan2 ((zbbbs - self%Z_mag), (rbbbs - self%R_mag))
    self%r_bound = sqrt( (rbbbs - self%R_mag)**2 + (zbbbs - self%Z_mag)**2 )

    call sort(self%thetab, self%r_bound, zbbbs, rbbbs)

    ! Allow for duplicated points near +- pi:
    if (self%thetab(1) == self%thetab(2)) then
       self%thetab(1) = self%thetab(1) + twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    if (self%thetab(self%nbbbs-1) == self%thetab(self%nbbbs)) then
       self%thetab(self%nbbbs) = self%thetab(self%nbbbs) - twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    ! It isn't likely that a duplicate point would exist near theta = 0,
    ! so I am not allowing this possibility for now.
    do i = 1, self%nbbbs-1
       if (self%thetab(i+1) == self%thetab(i)) then
          ! put in kluge for duplicate points, which happens near theta=0:
          self%thetab(i+1) = self%thetab(i+1) +  1.e-8
       end if
    end do

    call a_minor(rbbbs, zbbbs, self%Z_mag, avgrmid)
    self%aminor = avgrmid

    deallocate (rbbbs, zbbbs)

    self%r_bound = self%r_bound / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    self%eq_R = self%eq_R / self%aminor
    self%eq_Z = self%eq_Z / self%aminor

    do i = 1, self%nw
       self%R_psi(i, :) = self%eq_R(i)
       self%Z_psi(i, :) = self%eq_Z
    end do

    ! should rmaj be R_mag or rcentr? use R_mag for now.
    rmaj = self%R_mag
    B_T0 = abs(bcentr)
    self%B_T = B_T0
    self%psi_a = self%psi_a / (self%B_T * self%aminor**2)
    self%psi_0 = self%psi_0 / (self%B_T * self%aminor**2)
    psi_0_out = self%psi_0 ; psi_a_out = self%psi_a

    p_0 = self%pressure(1)
    self%fp = self%fp / (self%B_T * self%aminor)
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / p_0
    self%eqpsi_2d = self%eqpsi_2d / (self%B_T * self%aminor**2)

    do i = 1, self%nw
       do j = 1, self%nh
          if (self%eq_Z(j) == self%Z_mag .and. self%eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (self%eq_Z(j) - self%Z_mag), &
                  (self%eq_R(i) - self%R_mag))
          end if
       end do
    end do

    if (.true.) then
       write (*,*) "Finished efitin... imported EFIT equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", self%aminor
       write (*,*) 'R_mag', self%R_mag
       write (*,*) 'B_T0', self%B_T
       write (*,*) 'beta', self%beta_0
    end if

    ! To allow consistency with other impls in gradient etc.
    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0
    self%B_psi = 0.0
    self%initialised = .true.
    self%has_full_theta_range = .true.
1000 format(5(a10),i2,i4,i4)
2020 format (5e16.9)
2022 format (2i5)
  end subroutine efitin

  pure subroutine alloc_special_arrays(self)
    class(eeq_base_type), intent(in out) :: self
    integer :: nw, nh
    nw = self%nw ; nh = self%nh
    allocate (self%psi_bar(nw), self%fp(nw), self%qsf(nw), self%beta(nw))
    allocate (self%eq_R(nw), self%eq_Z(nh))
  end subroutine alloc_special_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_special_arrays(self)
    implicit none
    class(eeq_base_type), intent(in out) :: self
    deallocate(self%psi_bar, self%fp, self%qsf, self%beta)
    deallocate(self%eq_R, self%eq_Z)
    deallocate(self%thetab, self%r_bound)
  end subroutine dealloc_special_arrays

  subroutine setup_splines(self)
    use splines, only: new_spline, new_periodic_spline
    use constants, only: twopi
    implicit none
    class(eeq_base_type), intent(in out) :: self
    call new_spline(self%nw, self%psi_bar, self%fp, self%btori_spl)
    call new_spline(self%nw, self%psi_bar, self%qsf, self%qsf_spl)
    call new_spline(self%nw, self%psi_bar, self%pressure, self%pressure_spl)
    call new_spline(self%nw, self%psi_bar, self%beta, self%beta_spl)
    call new_periodic_spline(self%nbbbs, self%thetab, self%r_bound, twopi, &
         self%rbound_spl, drop_last_point = .true.)
  end subroutine setup_splines

  subroutine delete_splines(self)
    use splines, only: delete_spline, delete_periodic_spline
    implicit none
    class(eeq_base_type), intent(in out) :: self
    call delete_spline(self%btori_spl)
    call delete_spline(self%qsf_spl)
    call delete_spline(self%pressure_spl)
    call delete_spline(self%beta_spl)
    call delete_periodic_spline(self%rbound_spl)
  end subroutine delete_splines

  !> Calculate the derivative of f w.r.t. R, Z
  !> - dfm(:,:,1) is deriv w.r.t. R
  !> - def(:,:,2) is deriv w.r.t. Z
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nw, nh
    nw = self%nw ; nh = self%nh
    i = 1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))

    i = nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    j = 1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2))

    j = nh
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2))

    do i = 2, nw-1
       dfm(i,:,1) = 0.5*(f(i+1,:)-f(i-1,:))
    end do

    select case (char)
    case ('T')
       do j = 2, nh-1
          do i = 1, nw
             if(f(i,j+1)-f(i,j-1) > pi) then
                dfm(i,j,2) = 0.5*(f(i,j+1)-f(i,j-1)-2.*pi)
             else
                dfm(i,j,2) = 0.5*(f(i,j+1)-f(i,j-1))
             end if
          end do
       end do
    case('E','O')
       do j = 2, nh-1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j-1))
       end do
    end select
  end subroutine derm

  !> Calculate the derivative of psi/(a^2 B_a) w.r.t. R and Z. Parameters
  !> are:
  !> - rgrid: distance to mag axis as a fn of theta
  !> - theta: grid of theta
  !> - grad(:,1): deriv w.r.t. R as fn of theta
  !> - grad(:,2): deriv w.r.t  Z as fn of theta
  !> - char: if char = 'R', return the gradient of pressure instead
  !> - rp: value of psi/(a^2 B_a) on the flux surface where we want the grad
  !> - nth: number of theta points
  !> - ntm: lower index of theta array is -ntm
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    ! It may be clearer if we use our spline/periodic_spline wrappers?
    use splines, only: inter_d_cspl
    implicit none
    class(eeq_base_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    character(len=1), intent(in) ::  char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nw, self%nh, 2) :: dcart
    integer :: i

    select case(char)
    case('P', 'R')
       dcart = self%dpcart
    case('T')
       dcart = self%dtcart
    end select

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dcart(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dcart(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       ! This returns d (pressure/pressure(1)) / d psibar
       ! multiply by psi_N to get d (pressure/pressure(1)) / d (psi/ a^2 B_a)
       ! multiply by beta_0 to get d (pressure) / d (psi)
       ! However, I think there is a problem because dpm is in units of a^2 B_a ?
       call inter_d_cspl(self%nw, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       grad = grad*daa(1) * 0.5 * self%beta_0
    end if
  end subroutine gradient

  !> Calculate the derivative of psi/(a^2 B_a) w.r.t. R and Z and return
  !> the modulus sqrt(dpsi/dR ^ 2 + dpsi/dZ^2). I.e. return |grad psi|.
  !> Parameters are:
  !> - rgrid: distance to mag axis as a fn of theta
  !> - theta: grid of theta
  !> - grad(:,1): |grad psi|
  !> - grad(:,2): 0.0
  !> - char: if char = 'R', return |grad pressure| instead
  !> - char: if char = 'T', then return theta gradient in bishop form
  !>  - grad(:,1): (dtheta/dZ d psi/dR - dtheta/dR dpsi/dZ)/ |grad psi|
  !>  - grad(:,2): (dtheta/dR d psi/dR + dtheta/dZ dpsi/dZ)/ |grad psi|
  !> - rp: value of psi/(a^2 B_a) on the flux surface where we want the grad
  !> - nth: number of theta points
  !> - ntm: lower index of theta array is -ntm
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    ! It may be clearer if we use our spline/periodic_spline wrappers?
    use splines, only: inter_d_cspl
    implicit none
    class(eeq_base_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    character(len=1), intent(in) :: char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nw, self%nh, 2) ::  dbish
    integer :: i

    select case(char)
    case('P', 'R')
       dbish = self%dpbish
    case('T')
       dbish = self%dtbish
    end select

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dbish(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dbish(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       ! This returns d (pressure/pressure(1)) / d psibar
       ! multiply by psi_N to get d (pressure/pressure(1)) / d (psi/ a^2 B_a)
       ! multiply by beta_0 to get d (pressure) / d (psi)
       ! However, I think there is a problem because dpm is in units of a^2 B_a ?
       call inter_d_cspl(self%nw, self%eqpsi, self%pressure, 1, rpt, aa, daa)
       do i = -nth_used, nth_used
          grad = grad*daa(1) * 0.5 * self%beta_0 / (self%psi_a - self%psi_0)
       end do
    end if
  end subroutine bgradient

  !> fstar is f(R,Z) interpolated at the values (r,thetain). The parameter
  !> r is the distance to the magnetic axis.
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, intent (in) :: r_in, theta_in
    real, dimension(:, :), intent(in) :: f
    character(1), intent(in) :: char
    integer :: i, j, istar, jstar
    real :: st, dt, sr, dr, r
    real :: r_pos, z_pos
    !Copy input r to avoid any modification
    r = r_in

    r_pos = self%Rpos(r, theta_in)
    z_pos = self%Zpos(r, theta_in)

    ! find point on R mesh
    if (r_pos >= self%eq_R(self%nw) .or. r_pos <= self%eq_R(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of R domain'
       write(*,*) r, theta_in, self%eq_R(self%nw), r_pos
       call mp_abort('No evaluation of eqitem allowed outside or on edge of R domain')
    end if

    ! ensure point is on Z mesh
    if (z_pos >= self%eq_Z(self%nh) .or. z_pos <= self%eq_Z(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of Z domain'
       write(*,*) r, theta_in, self%eq_Z(1), self%eq_Z(self%nh), z_pos
       call mp_abort('No evaluation of eqitem allowed outside or on edge of Z domain')
    end if

    istar = 0
    do i = 2, self%nw
       if (istar /= 0) cycle
       if (r_pos < self%eq_R(i)) then
          dr = r_pos - self%eq_R(i-1)
          sr = self%eq_R(i) - r_pos
          istar = i-1
       end if
    end do

    ! Now do Z direction
    jstar = 0
    do j = 1, self%nh
       if (jstar /= 0) cycle
       if (z_pos < self%eq_Z(j)) then
          dt = z_pos - self%eq_Z(j-1)
          st = self%eq_Z(j) - z_pos
          jstar = j-1
       end if
    end do

    ! use opposite area stencil to interpolate
    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar / (abs(self%eq_R(istar+1)-self%eq_R(istar)) * &
         (self%eq_Z(jstar+1)-self%eq_Z(jstar)))
  end function eqitem

  !> FIXME : Add documentation
  real function diameter (self, rp)
    use constants, only: pi
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, intent (in) :: rp
    diameter = self%rfun(rp, 0.) + self%rfun(rp, pi)
  end function diameter

  !> FIXME : Add documentation
  real function rcenter (self, rp)
    use constants, only: pi
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, intent (in) :: rp
    rcenter = self%R_mag + &
         0.5 * (self%rfun(rp, 0.) - self%rfun(rp, pi))
  end function rcenter

  !> FIXME : Add documentation | Note this differs from all other types
  real function dpfun (self, pbar)
    use splines, only: dsplint
    implicit none
    class(eeq_base_type), intent(in) :: self
    real, intent(in) :: pbar
    dpfun = 0.5 * self%beta_0 * dsplint(pbar, self%pressure_spl)
  end function dpfun

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
    integer :: index
    real, dimension(nz) :: rtmp, ztmp
    integer i, j, i1, n
    type (spline) :: spl
!CMR, 28/10/08: add code to avoid duplicate points in 5 point spline
!               to determine r2 on inboard mid-plane
    integer :: k
    logical, parameter :: debug=.false.

    k = 0
    n = size(r)

    if (debug) write(6,*) "aminor:"
    if (debug) write(6,fmt='("r=",10f8.4)') r
    if (debug) write(6,fmt='("z=",10f8.4)') z

    if(n < nz) then
       write(*,*) 'nbbbs < nz -- very strange.  Stopping.'
       write(*,*) 'Look in eeq.f90.'
!      call mp_abort('nbbbs < nz -- very strange.  Stopping : Look in eeq.f90')
    endif

    j = 0
    do i = half_nz+1,1,-1
       j = j + 1
       ztmp(j) = z(i)
       rtmp(j) = r(i)
    enddo

    if(r(1) == r(n) .and. z(1) == z(n)) then
       do i = n-1, n-half_nz, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       enddo
    else
       do i = n, n-half_nz+1, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       enddo
    endif

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
       if(z(i)-Z_mag > 0.) then
          i1 = i - 1
          exit
       endif
    enddo

    if (i1 == invalidPoint) then
       write(*, '("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")')
       call mp_abort("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")
    endif

!CMR, 28/10/08: modify this code to avoid duplicate points in 5 point spline
!               to determine r2 on inboard mid-plane
    do i = 1, nz
       index = i1 - half_nz + i - 1 + k
       rtmp(i) = r(index)
       ztmp(i) = z(index)
       if ( i.gt.1 ) then
          ! If we find a duplicate point then we increase the
          ! duplicate point counter in k, used to offset index and
          ! store the next value instead Note, here we assume we can't
          ! get three duplicate values. A better approach may be to
          ! instead use a while loop where we keep going until we've
          ! found enough points.
          if ((rtmp(i)-rtmp(i-1))**2+(ztmp(i)-ztmp(i-1))**2 .lt. 1.0e-7) then
             k=k+1
             if (debug) write(6,fmt='("a_minor: duplicate pt, set k=",i2)') k
             rtmp(i) = r(index + 1)
             ztmp(i) = z(index + 1)
          endif
       endif
    enddo

    if (debug) write(6,fmt='("a_minor: rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("a_minor: ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("a_minor: Z_mag=",f8.4)') Z_mag

    call new_spline(nz, ztmp, rtmp, spl)
    r2 = splint(Z_mag, spl)
    call delete_spline(spl)

    a = (r2 - r1)/2.
    if (debug) write(6,*) "a_minor: r1,r2=",r1,r2
    if (debug) write(6,*) "a_minor: a=",a

  end subroutine a_minor
end module eeq
