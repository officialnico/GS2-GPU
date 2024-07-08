!> This module contains routines to read from CMR's GS2D equilibrium solver
module gs2d
  private
  public :: write_gs2d, read_gs2d, psi, b0, f, ippsi, nr, nz, p, ps
  public :: psip, psmin, q, r0, rgrid, rmag, rsep, zgrid, zmag, zsep

  real, dimension(:), allocatable :: ps, amin_gs2d, q, f, p, pp, rsep, zsep, rgrid, zgrid
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
