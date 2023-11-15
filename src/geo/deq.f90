!> FIXME : Add documentation
module deq
  use geo_utils, only: abstract_eqfile_cart_geo_type, EQFILE_LENGTH
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
     procedure :: calculate_gradients => dfit_init
     procedure :: setup_splines => setup_splines
     procedure :: delete_splines => delete_splines
     procedure :: rcenter => rcenter
     procedure :: diameter => diameter
     procedure :: gradient => gradient
     procedure :: bgradient => bgradient
     procedure :: rhofun => rhofun
     procedure :: derm => derm
     procedure :: eqitem => eqitem
     !procedure :: eqdcart => eqdcart
     !procedure :: eqdbish => eqdbish
     procedure :: alloc_arrays => alloc_arrays
     procedure :: dealloc_arrays => dealloc_arrays
  end type deq_type

contains

  !> This subroutine reads an DFIT output file containing
  !> the axisymmetric magnetic field geometry on a rectangular
  !> domain defined by the coordinates (R,Z).
  !
  !> eq_R     R grid
  !> eq_Z     Z grid
  !> fp    F on psibar grid
  !> eq_psi   psi on (R,Z) grid
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
    real, dimension(:), allocatable :: spsi_bar, sdfit_R, sdfit_Z, rbbbs, zbbbs, dummy
    real :: p_0, rwid, zhei, bcentr
    real :: delta_R, delta_Z, lref, bref
    integer :: i, j, jmin, jmax, nw_in, nh_in
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
       open(unit=5, file = trim(inputs%eqnormfile), status='old', form='formatted')
       read(5,*) lref, bref
       close(5)
    end if

    self%filename = trim(adjustl(inputs%eqfile))
    open(unit = 5, file = self%filename, status='old', form='formatted')

    ! Read the data
    read(5,*) nw_in, nh_in, self%nrho
    ! Note here we don't allow for rescaling the final grids with big
    ! unlike in eeq. This makes many of the spline interpolation calls
    ! in this routine effectively just copies
    self%nw = nw_in ; self%nh = nh_in
    call self%alloc_arrays
    allocate(spsi_bar(nw_in), sdfit_R(nw_in), sdfit_Z(nh_in), dummy(nw_in))

    read(5,2020) rwid, zhei
    read(5,2020) self%R_mag, self%Z_mag
    read(5,2020) self%psi_0, self%psi_a, bcentr

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
       self%eq_R(i) = rwid * float(i) / float(self%nw)
    end do

    do j = 1, self%nh
       self%eq_Z(j) = zhei*(float(j-1)/float(self%nh-1)-0.5)
    end do

    read(5,2020) (dummy(j), j = 1, nw_in)  ! pressure read
    call inter_cspl(nw_in, spsi_bar, dummy, self%nw, self%psi_bar, self%pressure)
    ! Note we read directly into the eq_psi array here so self%nw == nw_in etc.
    ! has to be true for this to work. See efit_in for a more general approach.
    read(5,2020) ((self%eq_psi(i,j) , i = 1, nw_in) , j = 1, nh_in)

    ! Only needed for irho = 4
    read (5,2020) (self%rho_mid(i), i = 1, self%nrho)
    read (5,2020) (self%psi_mid(i), i = 1, self%nrho)

    close(5)

    self%fp = 0.
    self%qsf = 0.
    self%nbbbs = 2*(self%nw+self%nh)-5

    ! stay away from edge of box
    delta_R = 1.02 ; delta_Z = 1.02

    allocate(rbbbs(self%nbbbs), zbbbs(self%nbbbs))
    allocate(self%thetab(self%nbbbs), self%r_bound(self%nbbbs))
    jmin = 1
    jmax = self%nh / 2
    do j = jmin, jmax
       rbbbs(j) = (self%eq_R(1)-0.5*rwid) / delta_R + 0.5*rwid
       zbbbs(j) = self%eq_Z(j+self%nh/2) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nw - 2
    do j = jmin, jmax
       rbbbs(j) = (self%eq_R(j-jmin+2)-0.5*rwid) / delta_R + 0.5*rwid
       zbbbs(j) = self%eq_Z(self%nh) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nh - 1
    do j = jmin, jmax
       rbbbs(j) = (self%eq_R(self%nw)-0.5*rwid) / delta_R + 0.5*rwid
       zbbbs(j) = self%eq_Z(jmax-j+2) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nw - 1
    do j = jmin, jmax
       rbbbs(j) = (self%eq_R(jmax-j+2)-0.5*rwid) / delta_R + 0.5*rwid
       zbbbs(j) = self%eq_Z(1) / delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + self%nh/2
    do j = jmin, jmax
       rbbbs(j) = (self%eq_R(1)-0.5*rwid) / delta_R + 0.5*rwid
       zbbbs(j) = self%eq_Z(j-jmin+2) / delta_Z
    end do

    do j = 1, self%nbbbs-1
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

    if (self%thetab(self%nbbbs-1) == self%thetab(self%nbbbs)) then
       self%thetab(self%nbbbs) = self%thetab(self%nbbbs) - twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    ! It isn't likely that a duplicate point would exist near theta = 0,
    ! so I am not allowing this possibility for now.
    do i = 1, self%nbbbs-1
       if (self%thetab(i) == self%thetab(i+1)) then
          write(*,*) 'Duplicates near theta = 0 not allowed.'
          write(*,*) i, i+1, ' Stopping.'
          call mp_abort('Duplicates near theta = 0 not allowed.')
       end if
    end do
    deallocate (rbbbs, zbbbs)

    if (lref /= 0.) then
       self%aminor = lref
       write(6,*) 'aminor = ', self%aminor, lref
    else
       self%aminor = self%R_mag
    end if
    avgrmid = self%aminor

    self%r_bound = self%r_bound / self%aminor

    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor

    self%eq_R = self%eq_R / self%aminor
    self%eq_Z = self%eq_Z / self%aminor

    ! should rmaj be R_mag? use R_mag for now.
    rmaj = self%R_mag
    if (bref /= 0.) then
       B_T0 = bref
       write(6,*) 'B_T0 = ', B_T0, bref
    else
       B_T0 = abs(bcentr)
    end if
    self%B_T = B_T0
    self%psi_a = self%psi_a / (self%B_T * self%aminor**2)
    self%psi_0 = self%psi_0 / (self%B_T * self%aminor**2)
    psi_0_out = self%psi_0 ; psi_a_out = self%psi_a
    p_0 = self%pressure(1)

    ! Pointless as fp == 0 here
    self%fp = self%fp / (self%B_T * self%aminor)
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 8. * pi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / p_0

    self%eq_psi = self%eq_psi / (self%B_T * self%aminor**2)

    self%eq_dR = self%eq_R(2) - self%eq_R(1)
    self%eq_dZ = self%eq_Z(2) - self%eq_Z(1)

    nthg = -1
    self%initialised = .true.
2020 format (5e16.8)

  end subroutine dfitin

  pure subroutine alloc_arrays(self)
    class(deq_type), intent(in out) :: self
    integer :: nw, nh, nrho
    nw = self%nw ; nh = self%nh ; nrho = self%nrho
    allocate (self%rho_mid(nrho), self%psi_mid(nrho))
    allocate (self%psi_bar(nw), self%fp(nw), self%qsf(nw), self%pressure(nw), self%beta(nw))
    allocate (self%eq_R(nw), self%eq_Z(nh))
    allocate (self%eq_psi(nw, nh))
    allocate (self%dpm(nw, nh, 2), self%dtm(nw, nh, 2))
  end subroutine alloc_arrays

  !> FIXME : Add documentation
  pure subroutine dealloc_arrays(self)
    implicit none
    class(deq_type), intent(in out) :: self
    if (allocated(self%rho_mid)) deallocate(self%rho_mid, self%psi_mid, self%beta)
    if (allocated(self%psi_bar)) deallocate(self%psi_bar, self%fp, self%qsf, self%pressure)
    if (allocated(self%eq_R)) deallocate(self%eq_R, self%eq_Z, self%eq_psi)
    if (allocated(self%dpm)) deallocate(self%dpm, self%dtm)
    if (allocated(self%thetab)) deallocate(self%thetab, self%r_bound)
  end subroutine dealloc_arrays

  !> FIXME : Add documentation
  subroutine dfit_init(self)
    implicit none
    class(deq_type), intent(in out) :: self
    real, dimension(self%nw, self%nh) :: eqth
    integer :: i, j
    do i = 1, self%nw
       do j = 1, self%nh
          if (self%eq_Z(j) == self%Z_mag .and. self%eq_R(i) == self%R_mag) then
             eqth(i, j) = 0.  ! value should not matter
          else
             eqth(i, j) = atan2( (self%eq_Z(j) - self%Z_mag), (self%eq_R(i) - self%R_mag))
          end if
       end do
    end do

    call self%derm(self%eq_psi, self%dpm, 'E')
    call self%derm(eqth, self%dtm, 'T')
  end subroutine dfit_init


  subroutine setup_splines(self)
    use splines, only: new_spline, new_periodic_spline
    use constants, only: twopi
    implicit none
    class(deq_type), intent(in out) :: self
    call new_spline(self%nw, self%psi_bar, self%fp, self%btori_spl)
    call new_spline(self%nw, self%psi_bar, self%qsf, self%qsf_spl)
    call new_spline(self%nw, self%psi_bar, self%pressure, self%pressure_spl)
    call new_spline(self%nw, self%psi_bar, self%beta, self%beta_spl)
    call new_periodic_spline(self%nbbbs, self%thetab, self%r_bound, twopi, &
         self%rbound_spl, drop_last_point = .true.)
    call new_spline(self%nrho, self%psi_mid, self%rho_mid, self%rho_spl)
  end subroutine setup_splines

  subroutine delete_splines(self)
    use splines, only: delete_spline, delete_periodic_spline
    implicit none
    class(deq_type), intent(in out) :: self
    call delete_spline(self%btori_spl)
    call delete_spline(self%qsf_spl)
    call delete_spline(self%pressure_spl)
    call delete_spline(self%beta_spl)
    call delete_periodic_spline(self%rbound_spl)
    call delete_spline(self%rho_spl)
  end subroutine delete_splines

  !> FIXME : Add documentation
  !> DFIT grid is equally spaced in R, Z -- this routine uses that fact and
  !> is therefore not completely general.  It is fine for DFIT output.
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(deq_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(len=1), intent(in) :: char
    integer :: i, j, nw, nh
    nw = self%nw ; nh = self%nh
    i = 1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:)) / self%eq_dR

    i = nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:)) / self%eq_dR

    j = 1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2)) / self%eq_dZ

    j = nh
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2)) / self%eq_dZ

    do i = 2, nw-1
       dfm(i,:,1) = 0.5*(f(i+1,:)-f(i-1,:)) / self%eq_dR
    end do

    select case (char)
    case ('T')
       do j = 2, nh-1
          do i = 1, nw
             if (f(i,j+1)-f(i,j-1) > pi) then
                dfm(i,j,2) = 0.5*(f(i,j+1)-f(i,j-1)-2.*pi) / self%eq_dZ
             else
                dfm(i,j,2) = 0.5*(f(i,j+1)-f(i,j-1)) / self%eq_dZ
             end if
          end do
       end do
    case('E')
       do j =2, nh-1
          dfm (:,j,2) = 0.5*(f(:,j+1)-f(:,j-1)) / self%eq_dZ
       end do
    end select

  end subroutine derm

  !> FIXME : Add documentation
  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(deq_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    character(len=1), intent(in) :: char
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    integer :: i
    grad = 0.
    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), self%dpm(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), self%dpm(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nw, self%psi_bar, self%pressure, 1, rpt, aa, daa)
       grad = grad * daa(1) * 0.5 * self%beta_0 / (self%psi_a - self%psi_0)
    end if
  end subroutine gradient

  !> FIXME : Add documentation
  subroutine bgradient(self, rgrid, theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    class(deq_type), intent(in) :: self
    integer, intent(in) :: nth_used, ntm
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    character(len=1), intent(in) :: char
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(self%nw, self%nh, 2) ::  dbish
    integer :: i
    ! Other eq types precalculates this in the init stage, e.g. calling eqdcart etc.
    dbish(:,:,1) = sqrt(self%dpm(:,:,1)**2 + self%dpm(:,:,2)**2)
    dbish(:,:,2) = 0.

    if (char == 'T') then  ! the order of the next two statements matters
       dbish(:,:,2) = (self%dtm(:,:,2) * self%dpm(:,:,1) - &
            self%dtm(:,:,1) * self%dpm(:,:,2)) / dbish(:,:,1)
       dbish(:,:,1) = (self%dtm(:,:,1) * self%dpm(:,:,1) + &
            self%dtm(:,:,2) * self%dpm(:,:,2)) / dbish(:,:,1)
    end if

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), dbish(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), dbish(:,:,2), 'Z')
    end do

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    if (char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(self%nw, self%psi_bar, self%pressure, 1, rpt, aa, daa)
       grad = grad * daa(1) * 0.5 * self%beta_0 / (self%psi_a - self%psi_0)
    end if
  end subroutine bgradient

  !> FIXME : Add documentation
  real function eqitem(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    implicit none
    class(deq_type), intent(in) :: self
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
    fstar = fstar / (abs(self%eq_R(istar+1) - self%eq_R(istar)) * &
         (self%eq_Z(jstar+1) - self%eq_Z(jstar)))
  end function eqitem

  !> FIXME : Add documentation
  real function diameter (self, rp)
    use constants, only: pi
    implicit none
    class(deq_type), intent(in) :: self
    real, intent (in) :: rp
    diameter = self%rfun(rp, 0.) + self%rfun(rp, pi)
  end function diameter

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
