!> FIXME : Add documentation
module gs2_heating

  implicit none
  private

  public :: heating_diagnostics
  public :: init_htype, zero_htype, del_htype
  public :: avg_h, avg_hk
  public :: hk_repack
  public :: dens_vel_diagnostics
  public :: init_dvtype, zero_dvtype, del_dvtype
  public :: avg_dv, avg_dvk, get_heat
   
  !> Allocate the array components of [[heating_diagnostics]]
  !>
  !> FIXME: Could be replaced with `elemental` constructor
  interface init_htype
     module procedure init_htype_0, init_htype_1, init_htype_2, init_htype_3
  end interface

  interface zero_htype
     module procedure zero_htype_0, zero_htype_1, zero_htype_2, zero_htype_3
  end interface

  interface del_htype
     module procedure del_htype_0, del_htype_1, del_htype_2, del_htype_3
  end interface

  !GGH
  interface init_dvtype
     module procedure init_dvtype_0, init_dvtype_1, init_dvtype_2, init_dvtype_3
  end interface

  interface zero_dvtype
     module procedure zero_dvtype_0, zero_dvtype_1, zero_dvtype_2, zero_dvtype_3
  end interface

  interface del_dvtype
     module procedure del_dvtype_0, del_dvtype_1, del_dvtype_2, del_dvtype_3
  end interface

  !> FIXME : Add documentation
  type :: heating_diagnostics
! total quantities:
     real :: energy
     real :: energy_dot
     real :: antenna
     real :: eapar                                       !int k_perp^2 A_par^2/8 pi
     real :: ebpar                                       !int B_par^2/8 pi
! species by species:
     real, dimension(:), allocatable :: delfs2      !int T/F0 dfs^2/2
     real, dimension(:), allocatable :: hs2         !int T/F0 hs^2/2
     real, dimension(:), allocatable :: phis2       !int q^2 n/T phi^2/2
     real, dimension(:), allocatable :: hypervisc
     real, dimension(:), allocatable :: hyperres
     real, dimension(:), allocatable :: hypercoll
     real, dimension(:), allocatable :: collisions
     real, dimension(:), allocatable :: imp_colls
     real, dimension(:), allocatable :: gradients
!     real, dimension(:), allocatable :: curvature
     real, dimension(:), allocatable :: heating
  end type heating_diagnostics

  !> Density-velocity perturbations
  type :: dens_vel_diagnostics
     !GGH NOTE: Dimension here is for species
     real, dimension(:), allocatable :: dvpar
     real, dimension(:), allocatable :: dvperp
     real, dimension(:), allocatable :: dn
  end type dens_vel_diagnostics

contains

  !> FIXME : Add documentation  
  subroutine init_htype_0 (h, nspec)
    implicit none
    type (heating_diagnostics), intent(in out) :: h
    integer, intent (in) :: nspec
    
    allocate (h % delfs2(nspec))      
    allocate (h % hs2(nspec))      
    allocate (h % phis2(nspec))      
    allocate (h % hypervisc(nspec))   
    allocate (h % hyperres(nspec))    
    allocate (h % hypercoll(nspec))  
    allocate (h % collisions(nspec))  
    allocate (h % imp_colls(nspec))  
    allocate (h % gradients(nspec))   
!    allocate (h % curvature(nspec))   
    allocate (h % heating(nspec))     

    call zero_htype (h)

  end subroutine init_htype_0

  !> FIXME : Add documentation
  subroutine init_htype_1 (h, nspec)
    implicit none
    type (heating_diagnostics), dimension(:), intent(in out) :: h
    integer, intent (in) :: nspec
    integer :: n, nmax

    nmax = size(h)
    
    do n=1,nmax
       allocate (h(n) % delfs2(nspec))      
       allocate (h(n) % hs2(nspec))      
       allocate (h(n) % phis2(nspec))      
       allocate (h(n) % hypervisc(nspec))   
       allocate (h(n) % hyperres(nspec))    
       allocate (h(n) % hypercoll(nspec))  
       allocate (h(n) % collisions(nspec))  
       allocate (h(n) % imp_colls(nspec))  
       allocate (h(n) % gradients(nspec))   
       !    allocate (h(n) % curvature(nspec))   
       allocate (h(n) % heating(nspec))     
    end do

    call zero_htype (h)

  end subroutine init_htype_1

  !> FIXME : Add documentation  
  subroutine init_htype_2 (h, nspec)
    implicit none
    type (heating_diagnostics), dimension(:,:), intent(in out) :: h
    integer, intent (in) :: nspec
    integer :: m, n, mmax, nmax

    mmax = size(h, 1)
    nmax = size(h, 2)

    do n = 1, nmax
       do m = 1, mmax
          allocate (h(m,n) % delfs2(nspec))      
          allocate (h(m,n) % hs2(nspec))      
          allocate (h(m,n) % phis2(nspec))      
          allocate (h(m,n) % hypervisc(nspec))   
          allocate (h(m,n) % hyperres(nspec))    
          allocate (h(m,n) % hypercoll(nspec))  
          allocate (h(m,n) % collisions(nspec))  
          allocate (h(m,n) % imp_colls(nspec))  
          allocate (h(m,n) % gradients(nspec))   
          !    allocate (h(m,n) % curvature(nspec))   
          allocate (h(m,n) % heating(nspec))     
       end do
    end do

    call zero_htype (h)

  end subroutine init_htype_2

  !> FIXME : Add documentation
  subroutine init_htype_3 (h, nspec)
    implicit none
    type (heating_diagnostics), dimension(:,:,:), intent(in out) :: h
    integer, intent (in) :: nspec
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             allocate (h(l,m,n) % delfs2(nspec))      
             allocate (h(l,m,n) % hs2(nspec))      
             allocate (h(l,m,n) % phis2(nspec))      
             allocate (h(l,m,n) % hypervisc(nspec))   
             allocate (h(l,m,n) % hyperres(nspec))    
             allocate (h(l,m,n) % hypercoll(nspec))  
             allocate (h(l,m,n) % collisions(nspec))  
             allocate (h(l,m,n) % imp_colls(nspec))  
             allocate (h(l,m,n) % gradients(nspec))   
             !    allocate (h(l,m,n) % curvature(nspec))   
             allocate (h(l,m,n) % heating(nspec))     
          end do
       end do 
    end do

    call zero_htype (h)

  end subroutine init_htype_3

  !> FIXME : Add documentation
  subroutine zero_htype_0 (h)
    implicit none
    type (heating_diagnostics), intent(in out) :: h
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    h % delfs2 = 0. 
    h % hs2 = 0. 
    h % phis2 = 0. 
    h % hypervisc = 0. 
    h % hyperres = 0.  
    h % hypercoll = 0.  
    h % collisions = 0.
    h % imp_colls = 0.
    h % gradients = 0. 
!    h % curvature = 0. 
    h % heating = 0.   

  end subroutine zero_htype_0

  !> FIXME : Add documentation
  subroutine zero_htype_1 (h)
    implicit none
    type (heating_diagnostics), dimension(:), intent(in out) :: h
    integer :: n, nmax

    nmax = size(h)
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       h(n) % delfs2 = 0. 
       h(n) % hs2 = 0. 
       h(n) % phis2 = 0. 
       h(n) % hypervisc = 0. 
       h(n) % hyperres = 0.  
       h(n) % hypercoll = 0.  
       h(n) % collisions = 0.
       h(n) % imp_colls = 0.
       h(n) % gradients = 0. 
       !    h(n) % curvature = 0. 
       h(n) % heating = 0.   
    end do

  end subroutine zero_htype_1

  !> FIXME : Add documentation
  subroutine zero_htype_2 (h)
    implicit none
    type (heating_diagnostics), dimension(:,:), intent(in out) :: h
    integer :: n, nmax, m, mmax

    mmax = size(h, 1)
    nmax = size(h, 2)
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       do m=1,mmax
          h(m,n) % delfs2 = 0. 
          h(m,n) % hs2 = 0. 
          h(m,n) % phis2 = 0. 
          h(m,n) % hypervisc = 0. 
          h(m,n) % hyperres = 0.  
          h(m,n) % hypercoll = 0.  
          h(m,n) % collisions = 0.
          h(m,n) % imp_colls = 0.
          h(m,n) % gradients = 0. 
          !    h(m,n) % curvature = 0. 
          h(m,n) % heating = 0.   
       end do
    end do

  end subroutine zero_htype_2

  !> FIXME : Add documentation
  subroutine zero_htype_3 (h)
    implicit none
    type (heating_diagnostics), dimension(:,:,:), intent(in out) :: h
    integer :: n, nmax, m, mmax, l, lmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)
        
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       do m=1,mmax
          do l=1,lmax
             h(l,m,n) % delfs2 = 0. 
             h(l,m,n) % hs2 = 0. 
             h(l,m,n) % phis2 = 0. 
             h(l,m,n) % hypervisc = 0. 
             h(l,m,n) % hyperres = 0.  
             h(l,m,n) % hypercoll = 0.  
             h(l,m,n) % collisions = 0.
             h(l,m,n) % imp_colls = 0.
             h(l,m,n) % gradients = 0. 
             !    h(l,m,n) % curvature = 0. 
             h(l,m,n) % heating = 0.   
          end do
       end do
    end do

  end subroutine zero_htype_3

  !> FIXME : Add documentation
  subroutine del_htype_0 (h)
    implicit none
    type (heating_diagnostics), intent(in out) :: h
    
    deallocate (h % delfs2)
    deallocate (h % hs2)
    deallocate (h % phis2)
    deallocate (h % hypervisc)
    deallocate (h % hyperres)
    deallocate (h % hypercoll)
    deallocate (h % collisions)
    deallocate (h % imp_colls)
    deallocate (h % gradients)
!    deallocate (h % curvature)
    deallocate (h % heating)

  end subroutine del_htype_0

  !> FIXME : Add documentation  
  subroutine del_htype_1 (h)
    implicit none
    type (heating_diagnostics), dimension(:), intent(in out) :: h
    integer :: m, mmax
    
    mmax = size (h, 1)
    
    do m=1,mmax
       deallocate (h(m) % delfs2)
       deallocate (h(m) % hs2)
       deallocate (h(m) % phis2)
       deallocate (h(m) % hypervisc)
       deallocate (h(m) % hyperres)
       deallocate (h(m) % hypercoll)
       deallocate (h(m) % collisions)
       deallocate (h(m) % imp_colls)
       deallocate (h(m) % gradients)
       !       deallocate (h(m,n) % curvature)
       deallocate (h(m) % heating)
    end do

  end subroutine del_htype_1

  !> FIXME : Add documentation
  subroutine del_htype_2 (h)
    implicit none
    type (heating_diagnostics), dimension(:,:), intent(in out) :: h
    integer :: m, mmax, n, nmax
    
    mmax = size (h, 1)
    nmax = size (h, 2)
    
    do n=1,nmax
       do m=1,mmax
          deallocate (h(m,n) % delfs2)
          deallocate (h(m,n) % hs2)
          deallocate (h(m,n) % phis2)
          deallocate (h(m,n) % hypervisc)
          deallocate (h(m,n) % hyperres)
          deallocate (h(m,n) % hypercoll)
          deallocate (h(m,n) % collisions)
          deallocate (h(m,n) % imp_colls)
          deallocate (h(m,n) % gradients)
          !    deallocate (h(m,n) % curvature)
          deallocate (h(m,n) % heating)
       end do
    end do

  end subroutine del_htype_2

  !> FIXME : Add documentation
  subroutine del_htype_3 (h)
    implicit none
    type (heating_diagnostics), dimension(:,:,:), intent(in out) :: h
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             deallocate (h(l,m,n) % delfs2)      
             deallocate (h(l,m,n) % hs2)      
             deallocate (h(l,m,n) % phis2)      
             deallocate (h(l,m,n) % hypervisc)   
             deallocate (h(l,m,n) % hyperres)    
             deallocate (h(l,m,n) % hypercoll)  
             deallocate (h(l,m,n) % collisions)  
             deallocate (h(l,m,n) % imp_colls)  
             deallocate (h(l,m,n) % gradients)   
             !    deallocate (h(l,m,n) % curvature)   
             deallocate (h(l,m,n) % heating)     
          end do
       end do 
    end do

  end subroutine del_htype_3

  !> Calculate the average of various heating quantities
  subroutine avg_h (h, h_hist, istep, navg)
    use mp, only: proc0
    use species, only: nspec
    implicit none
    !> The heating diagnostics at the current timestep
    type (heating_diagnostics), intent(in out) :: h
    !> The last [[gs2_diagnostics:navg]] timesteps of heating diagnostics
    type (heating_diagnostics), dimension (0:), intent(in out) :: h_hist
    integer, intent (in) :: istep, navg
    integer :: is, i

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             h_hist(mod(istep,navg)) % energy     = h % energy
             h_hist(mod(istep,navg)) % energy_dot = h % energy_dot
             h_hist(mod(istep,navg)) % antenna    = h % antenna
             h_hist(mod(istep,navg)) % eapar      = h % eapar
             h_hist(mod(istep,navg)) % ebpar      = h % ebpar
             h_hist(mod(istep,navg)) % delfs2     = h % delfs2
             h_hist(mod(istep,navg)) % hs2        = h % hs2
             h_hist(mod(istep,navg)) % phis2      = h % phis2
             h_hist(mod(istep,navg)) % hypervisc  = h % hypervisc
             h_hist(mod(istep,navg)) % hyperres   = h % hyperres
             h_hist(mod(istep,navg)) % hypercoll  = h % hypercoll
             h_hist(mod(istep,navg)) % collisions = h % collisions
             h_hist(mod(istep,navg)) % imp_colls  = h % imp_colls
             h_hist(mod(istep,navg)) % gradients  = h % gradients
!             h_hist(mod(istep,navg)) % curvature  = h % curvature
             h_hist(mod(istep,navg)) % heating    = h % heating
          end if
          
          if (istep >= navg) then
             call zero_htype(h)
             do i=0,navg-1
                h % energy     = h%energy     + h_hist(i) % energy / real(navg)
                h % energy_dot = h%energy_dot + h_hist(i) % energy_dot / real(navg)
                h % antenna    = h%antenna    + h_hist(i) % antenna / real(navg)
                h % eapar      = h%eapar      + h_hist(i) % eapar    / real(navg)
                h % ebpar      = h%ebpar      + h_hist(i) % ebpar    / real(navg)
                do is = 1,nspec
                   h % delfs2(is)     = h%delfs2(is)     + h_hist(i) % delfs2(is) / real(navg)
                   h % hs2(is)        = h%hs2(is)        + h_hist(i) % hs2(is) / real(navg)
                   h % phis2(is)      = h%phis2(is)      + h_hist(i) % phis2(is) / real(navg)
                   h % hypervisc(is)  = h%hypervisc(is)  + h_hist(i) % hypervisc(is) / real(navg)
                   h % hyperres(is)   = h%hyperres(is)   + h_hist(i) % hyperres(is) / real(navg)
                   h % hypercoll(is)  = h%hypercoll(is)  + h_hist(i) % hypercoll(is) / real(navg)
                   h % collisions(is) = h%collisions(is) + h_hist(i) % collisions(is) / real(navg)
                   h % imp_colls(is) = h%imp_colls(is) + h_hist(i) % imp_colls(is) / real(navg)
                   h % gradients(is)  = h%gradients(is)  + h_hist(i) % gradients(is) / real(navg)
!                  h % curvature(is)  = h%curvature(is)  + h_hist(i) % curvature(is) / real(navg)
                   h % heating(is)    = h%heating(is)    + h_hist(i) % heating(is) / real(navg)
                end do
             end do
          end if
       end if
    end if

  end subroutine avg_h

  !> Calculate the average of various heating quantities as a function of \((\theta, k_y)\)
  subroutine avg_hk (hk, hk_hist, istep, navg)
    use mp, only: proc0
    use species, only: nspec
    implicit none
    !> Heating diagnostics as a function of \((\theta, k_y)\) at the current timestep
    type (heating_diagnostics), dimension(:,:), intent(in out) :: hk
    !> Heating diagnostics as a function of \((\theta, k_y)\) over the last [[gs2_diagnostics:navg]] timesteps
    type (heating_diagnostics), dimension(:,:,0:), intent(in out) :: hk_hist
    integer, intent (in) :: istep, navg
    integer :: is, i, m, n, mmax, nmax

    mmax = size(hk, 1)
    nmax = size(hk, 2)

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             hk_hist(:,:,mod(istep,navg)) % energy     = hk % energy
             hk_hist(:,:,mod(istep,navg)) % energy_dot = hk % energy_dot
             hk_hist(:,:,mod(istep,navg)) % antenna    = hk % antenna
             hk_hist(:,:,mod(istep,navg)) % eapar      = hk % eapar
             hk_hist(:,:,mod(istep,navg)) % ebpar      = hk % ebpar
            do n=1,nmax
                do m=1,mmax
                   do is = 1,nspec
                      hk_hist(m,n,mod(istep,navg)) % delfs2(is)     = hk(m,n) % delfs2(is)
                      hk_hist(m,n,mod(istep,navg)) % hs2(is)        = hk(m,n) % hs2(is)
                      hk_hist(m,n,mod(istep,navg)) % phis2(is)      = hk(m,n) % phis2(is)
                      hk_hist(m,n,mod(istep,navg)) % hypervisc(is)  = hk(m,n) % hypervisc(is)
                      hk_hist(m,n,mod(istep,navg)) % hyperres(is)   = hk(m,n) % hyperres(is)
                      hk_hist(m,n,mod(istep,navg)) % hypercoll(is)  = hk(m,n) % hypercoll(is)
                      hk_hist(m,n,mod(istep,navg)) % collisions(is) = hk(m,n) % collisions(is)
                      hk_hist(m,n,mod(istep,navg)) % imp_colls(is)  = hk(m,n) % imp_colls(is)
                      hk_hist(m,n,mod(istep,navg)) % gradients(is)  = hk(m,n) % gradients(is)
!                     hk_hist(m,n,mod(istep,navg)) % curvature(is)  = hk(m,n) % curvature(is)
                      hk_hist(m,n,mod(istep,navg)) % heating(is)    = hk(m,n) % heating(is)
                   end do
                end do
             end do
          end if
          
          if (istep >= navg) then
             call zero_htype (hk)
             do n=1,nmax
                do m=1,mmax
                   do i=0,navg-1
                      hk(m,n) % energy     = hk(m,n)%energy     + hk_hist(m,n,i) % energy / real(navg)
                      hk(m,n) % energy_dot = hk(m,n)%energy_dot + hk_hist(m,n,i) % energy_dot / real(navg)
                      hk(m,n) % antenna    = hk(m,n)%antenna    + hk_hist(m,n,i) % antenna / real(navg)
                      hk(m,n) % eapar      = hk(m,n)%eapar      + hk_hist(m,n,i) % eapar / real(navg)
                      hk(m,n) % ebpar      = hk(m,n)%ebpar      + hk_hist(m,n,i) % ebpar / real(navg)
                      do is = 1,nspec
                         hk(m,n)%delfs2(is)     = hk(m,n)%delfs2(is)    + hk_hist(m,n,i) % delfs2(is) / real(navg)
                         hk(m,n)%hs2(is)        = hk(m,n)%hs2(is)       + hk_hist(m,n,i) % hs2(is) / real(navg)
                         hk(m,n)%phis2(is)      = hk(m,n)%phis2(is)     + hk_hist(m,n,i) % phis2(is) / real(navg)
                         hk(m,n)%hypervisc(is)  = hk(m,n)%hypervisc(is) + hk_hist(m,n,i) % hypervisc(is) / real(navg)
                         hk(m,n)%hyperres(is)   = hk(m,n)%hyperres(is)  + hk_hist(m,n,i) % hyperres(is) / real(navg)
                         hk(m,n)%hypercoll(is)  = hk(m,n)%hypercoll(is) + hk_hist(m,n,i) % hypercoll(is) / real(navg)
                         hk(m,n)%collisions(is) = hk(m,n)%collisions(is)+ hk_hist(m,n,i) % collisions(is) / real(navg)
                         hk(m,n)%imp_colls(is)  = hk(m,n)%imp_colls(is) + hk_hist(m,n,i) % imp_colls(is) / real(navg)
                         hk(m,n)%gradients(is)  = hk(m,n)%gradients(is) + hk_hist(m,n,i) % gradients(is) / real(navg)
!                        hk(m,n)%curvature(is)  = hk(m,n)%curvature(is) + hk_hist(m,n,i) % curvature(is) / real(navg)
                         hk(m,n)%heating(is)    = hk(m,n)%heating(is)   + hk_hist(m,n,i) % heating(is) / real(navg)
                      end do
                   end do
                end do
             end do
          end if
       end if
    end if

  end subroutine avg_hk

  !> FIXME : Add documentation
  subroutine hk_repack (hk, i, tmp)
    use species, only: nspec
    implicit none    
    type (heating_diagnostics), dimension (:,:), intent(in) :: hk
    real, dimension(:,:,:), intent (out) :: tmp
    integer, intent (in) :: i
    integer :: is, m, n, mmax, nmax

    mmax = size(tmp, 1)
    nmax = size(tmp, 2)

    select case (i)
    case (1) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hypervisc(is)
             end do
          end do
       end do
    case (2) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hyperres(is)
             end do
          end do
       end do
    case (3) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%collisions(is)
             end do
          end do
       end do
    case (4) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%gradients(is)
             end do
          end do
       end do
    case (5) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%heating(is)
             end do
          end do
       end do
    case (6) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hypercoll(is)
             end do
          end do
       end do
    case (7) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%delfs2(is)
             end do
          end do
       end do
    case (8) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hs2(is)
             end do
          end do
       end do
    case (9) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%phis2(is)
             end do
          end do
       end do
    case (10) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%imp_colls(is)
             end do
          end do
       end do
    end select
    
  end subroutine hk_repack

  !=============================================================================
  !<GGH
  ! Density velocity perturbation routines
  !=============================================================================

  !> FIXME : Add documentation
  subroutine init_dvtype_0 (dv, nspec)
    implicit none
    type (dens_vel_diagnostics), intent(in out) :: dv
    integer, intent (in) :: nspec

       allocate (dv % dvpar(nspec))
       allocate (dv % dvperp(nspec))    
       allocate (dv % dn(nspec))  

    call zero_dvtype (dv)

  end subroutine init_dvtype_0
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine init_dvtype_1 (dv, nspec)
    implicit none
    type (dens_vel_diagnostics), dimension(:), intent(in out) :: dv
    integer, intent (in) :: nspec
    integer :: n, nmax

    nmax = size(dv)
    
    do n=1,nmax
       allocate (dv(n) % dvpar(nspec))
       allocate (dv(n) % dvperp(nspec))    
       allocate (dv(n) % dn(nspec))  
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_1
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine init_dvtype_2 (dv, nspec)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:), intent(in out) :: dv
    integer, intent (in) :: nspec
    integer :: m, n, mmax, nmax

    mmax = size(dv, 1)
    nmax = size(dv, 2)

    do n = 1, nmax
       do m = 1, mmax
          allocate (dv(m,n) % dvpar(nspec))
          allocate (dv(m,n) % dvperp(nspec))    
          allocate (dv(m,n) % dn(nspec))  
       end do
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_2

  !=============================================================================
  !> FIXME : Add documentation
  subroutine init_dvtype_3 (dv, nspec)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:,:), intent(in out) :: dv
    integer, intent (in) :: nspec
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(dv, 1)
    mmax = size(dv, 2)
    nmax = size(dv, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             allocate (dv(l,m,n) % dvpar(nspec))
             allocate (dv(l,m,n) % dvperp(nspec))    
             allocate (dv(l,m,n) % dn(nspec))  
          end do
       end do 
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_3

  !=============================================================================
  !> FIXME : Add documentation
  subroutine zero_dvtype_0 (dv)
    implicit none
    type (dens_vel_diagnostics), intent(in out) :: dv

    dv % dvpar = 0. 
    dv % dvperp = 0.  
    dv % dn = 0.

  end subroutine zero_dvtype_0
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine zero_dvtype_1 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:), intent(in out) :: dv
    integer :: n, nmax

    nmax = size(dv)
    
    do n=1,nmax
       dv(n) % dvpar = 0. 
       dv(n) % dvperp = 0.  
       dv(n) % dn = 0.
    end do

  end subroutine zero_dvtype_1

  !=============================================================================
  !> FIXME : Add documentation
  subroutine zero_dvtype_2 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:), intent(in out) :: dv
    integer :: n, nmax, m, mmax

    mmax = size(dv, 1)
    nmax = size(dv, 2)
        
    do n=1,nmax
       do m=1,mmax
          dv(m,n) % dvpar = 0. 
          dv(m,n) % dvperp = 0.  
          dv(m,n) % dn = 0.
       end do
    end do

  end subroutine zero_dvtype_2
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine zero_dvtype_3 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:,:), intent(in out) :: dv
    integer :: n, nmax, m, mmax, l, lmax

    lmax = size(dv, 1)
    mmax = size(dv, 2)
    nmax = size(dv, 3)
        
    do n=1,nmax
       do m=1,mmax
          do l=1,lmax
             dv(l,m,n) % dvpar = 0. 
             dv(l,m,n) % dvperp = 0.  
             dv(l,m,n) % dn = 0.
          end do
       end do
    end do

  end subroutine zero_dvtype_3

  !=============================================================================
  !> FIXME : Add documentation
  subroutine del_dvtype_0 (dv)
    implicit none
    type (dens_vel_diagnostics), intent(in out) :: dv
    
    deallocate (dv % dvpar)
    deallocate (dv % dvperp)
    deallocate (dv % dn)

  end subroutine del_dvtype_0
  
  !=============================================================================
  !> FIXME : Add documentation  
  subroutine del_dvtype_1 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:), intent(in out) :: dv
    integer :: m, mmax

    mmax = size (dv, 1)

    do m=1,mmax
       deallocate (dv(m) % dvpar)
       deallocate (dv(m) % dvperp)
       deallocate (dv(m) % dn)
    end do

  end subroutine del_dvtype_1

  !=============================================================================
  !> FIXME : Add documentation
  subroutine del_dvtype_2 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:), intent(in out) :: dv
    integer :: m, mmax, n, nmax
    
    mmax = size (dv, 1)
    nmax = size (dv, 2)
    
    do n=1,nmax
       do m=1,mmax
          deallocate (dv(m,n) % dvpar)
          deallocate (dv(m,n) % dvperp)
          deallocate (dv(m,n) % dn)
       end do
    end do

  end subroutine del_dvtype_2
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine del_dvtype_3 (dv)
    implicit none
    type (dens_vel_diagnostics), dimension(:,:,:), intent(in out) :: dv
    integer :: l, lmax, m, mmax, n, nmax
    
    lmax = size (dv, 1)
    mmax = size (dv, 2)
    nmax = size (dv, 3)
    
    do l=1,lmax
       do n=1,nmax
          do m=1,mmax
             deallocate (dv(l,m,n) % dvpar)
             deallocate (dv(l,m,n) % dvperp)
             deallocate (dv(l,m,n) % dn)
          end do
       end do
    end do

  end subroutine del_dvtype_3
  
  !=============================================================================
  !> FIXME : Add documentation
  subroutine avg_dv (dv, dv_hist, istep, navg)  
    use mp, only: proc0
    implicit none
    type (dens_vel_diagnostics), intent(in out) :: dv
    type (dens_vel_diagnostics), dimension (0:), intent(in out) :: dv_hist
    integer, intent (in) :: istep, navg
    integer :: i

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             dv_hist(mod(istep,navg)) % dvpar(:)     = dv % dvpar(:)
             dv_hist(mod(istep,navg)) % dvperp(:)    = dv % dvperp(:)
             dv_hist(mod(istep,navg)) % dn(:)        = dv % dn(:)
          end if
          
          if (istep >= navg) then
             call zero_dvtype(dv)
             do i=0,navg-1
                dv % dvpar(:)  = dv % dvpar(:)  + dv_hist(i) % dvpar(:)  / real(navg)
                dv % dvperp(:) = dv % dvperp(:) + dv_hist(i) % dvperp(:) / real(navg)
                dv % dn(:)     = dv % dn(:)     + dv_hist(i) % dn (:)    / real(navg)
             end do
          end if
       end if
    end if

  end subroutine avg_dv

  !=============================================================================
  !> FIXME : Add documentation
  subroutine avg_dvk (dvk, dvk_hist, istep, navg)    
    use mp, only: proc0
    implicit none
    type (dens_vel_diagnostics), dimension(:,:), intent(in out) :: dvk
    type (dens_vel_diagnostics), dimension (:,:,0:), intent(in out) :: dvk_hist
    integer, intent (in) :: istep, navg
    integer :: i, m, n, mmax, nmax

    mmax = size(dvk, 1)
    nmax = size(dvk, 2)

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             do n=1,nmax
                do m=1,mmax
                   dvk_hist(m,n,mod(istep,navg)) % dvpar(:)     = dvk(m,n) % dvpar(:)
                   dvk_hist(m,n,mod(istep,navg)) % dvperp(:)    = dvk(m,n) % dvperp(:)
                   dvk_hist(m,n,mod(istep,navg)) % dn(:)        = dvk(m,n) % dn(:)
             enddo
          enddo
          end if
          
          if (istep >= navg) then
             call zero_dvtype (dvk)
             do n=1,nmax
                do m=1,mmax
                   do i=0,navg-1
                      dvk(m,n) % dvpar(:)  = dvk(m,n) % dvpar(:)  + dvk_hist(m,n,i) % dvpar(:)  / real(navg)
                      dvk(m,n) % dvperp(:) = dvk(m,n) % dvperp(:) + dvk_hist(m,n,i) % dvperp(:) / real(navg)
                      dvk(m,n) % dn(:)     = dvk(m,n) % dn(:)     + dvk_hist(m,n,i) % dn (:)    / real(navg)
                   end do
                enddo
             enddo
          end if
       end if
    end if

  end subroutine avg_dvk

!=============================================================================
! Density: Calculate Density perturbations
 !=============================================================================

  !> Calculate some heating quantities:
  !> - ion/electron heating
  !> - antenna power and B-field contributions to E and E_dot
  !> - gradient contributions to heating
  !> - hyperviscosity
  !> - hyperresistivity
  !>
  !> FIXME: This is too large and does too much, split into smaller routines
  subroutine get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0
    use constants, only: zi
    use kt_grids, only: ntheta0, naky, aky, akx, kperp2
#ifdef LOWFLOW
    use lowflow, only: hneoc
#endif
    use dist_fn_arrays, only: vpac, aj0, aj1, vperp2, g, gnew, g_adjust, g_work, to_g_gs2, from_g_gs2, wstar
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: jacob, delthet, ntgrid
    use run_parameters, only: fphi, fapar, fbpar, beta, tite
    use gs2_time, only: code_dt, tunits
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp, hypervisc_filter

    implicit none
    type (heating_diagnostics), intent(in out) :: h
    type (heating_diagnostics), dimension(:,:), intent(in out) :: hk
    complex, dimension (-ntgrid:,:,:), intent(in) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension(:,:,:,:), allocatable :: tot
!    complex, dimension (:,:,:), allocatable :: epar
    complex, dimension(:,:,:), allocatable :: bpardot, apardot, phidot, j_ext
    complex :: chi, havg
    complex :: chidot, j0phiavg, j1bparavg, j0aparavg
!    complex :: pstar, pstardot, gdot
    complex :: phi_m, apar_m, bpar_m, hdot
    complex :: phi_avg, bpar_avg, bperp_m, bperp_avg
    complex :: de, denew
    complex :: dgdt_hypervisc
    real, dimension (:), allocatable :: wgt
    real :: fac2, dtinv, akperp4
    integer :: isgn, iglo, ig, is, ik, it, ie

    g_work(ntgrid,:,:) = 0.

! ==========================================================================
! Ion/Electron heating------------------------------------------------------
! ==========================================================================

    allocate ( phidot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (apardot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (bpardot(-ntgrid:ntgrid, ntheta0, naky))

    call dot ( phi,  phinew,  phidot, fphi)
    call dot (apar, aparnew, apardot, fapar)
    call dot (bpar, bparnew, bpardot, fbpar)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Next two calls make the variables g, gnew = h, hnew
!!! until the end of this procedure!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call g_adjust (g,    phi,    bpar,    direction = from_g_gs2)
    call g_adjust (gnew, phinew, bparnew, direction = from_g_gs2)

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1

             chidot = aj0(ig,iglo)*(phidot(ig,it,ik) &
                  - vpac(ig,isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
                  + aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpardot(ig,it,ik)*spec(is)%tz

             hdot = fdot (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo), dtinv)

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

! First term on RHS and LHS of Eq B-10 of H1:

             g_work(ig,isgn,iglo) = spec(is)%dens*conjg(havg)* &
                  (chidot*spec(is)%z-hdot*spec(is)%temp)

          end do
       end do
    end do

    deallocate (phidot, apardot, bpardot)

    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

    call integrate_moment (g_work, tot)

    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
!          wgt(ig) = delthet(ig)*jacob(ig)
! delthet is cell-centered, but jacob is given on grid
          wgt(ig) = delthet(ig)*(jacob(ig)+jacob(ig+1))*0.5
       end do
       wgt = wgt/sum(wgt)

       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                    hk(it,ik) % heating(is) = hk(it,ik) % heating(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % heating(is) = h % heating(is) + hk(it,ik) % heating(is)
             end do
          end do
       end do
    end if

! ==========================================================================
! Antenna Power and B-field contribution to E and E_dot---------------------
! ==========================================================================
    if (proc0) then
!       allocate (epar (-ntgrid:ntgrid, ntheta0, naky)) ; epar = 0.
       allocate (j_ext(-ntgrid:ntgrid, ntheta0, naky)) ; j_ext = 0.
       call antenna_apar (kperp2, j_ext)

       if (beta > epsilon(0.)) then
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   !GGH Time and space averaged estimate of d/dt(apar)
                   apar_m = fdot (apar   (ig  ,it,ik), &
                                  apar   (ig+1,it,ik), &
                                  aparnew(ig  ,it,ik), &
                                  aparnew(ig+1,it,ik), dtinv)*fapar

! J_ext.E when driving antenna only includes A_parallel:

                   hk(it,ik) % antenna = hk(it, ik) % antenna + real(conjg(j_ext(ig,it,ik))*apar_m)*wgt(ig)*fac2

                   !GGH Time and space averaged estimate of d/dt(bperp)
                   bperp_m = fdot (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                   aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), dtinv) * fapar

                   !GGH Time and space averaged estimate of d/dt(bpar)
                   bpar_m = fdot (bpar   (ig  ,it,ik), &
                                  bpar   (ig+1,it,ik), &
                                  bparnew(ig  ,it,ik), &
                                  bparnew(ig+1,it,ik), dtinv)*fbpar

                   !GGH Time and space averaged estimate of bperp
                   bperp_avg = favg (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                     aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik))) * fapar

                   !GGH Time and space averaged estimate of bpar
                   bpar_avg = favg (bpar   (ig  ,it,ik), &
                                    bpar   (ig+1,it,ik), &
                                    bparnew(ig  ,it,ik), &
                                    bparnew(ig+1,it,ik)) * fbpar

! 1/2 * d/dt B**2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        real(0.25 * conjg(bperp_m)*bperp_avg + conjg(bpar_m)*bpar_avg) &
                        * wgt(ig)*fac2*(2.0/beta)

! B**2/2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2
                   hk(it,ik) % energy = hk(it,ik) % energy &
                        + 0.5*real((0.25*conjg(bperp_avg)*bperp_avg + conjg(bpar_avg)*bpar_avg)) &
                        * wgt(ig)*fac2*(2.0/beta)

                   !Eapar = int k_perp^2 A_par^2/(8 pi)
                   hk(it,ik) % eapar = hk(it,ik) % eapar &
                        + 0.5*real(0.25*conjg(bperp_avg)*bperp_avg) * wgt(ig)*fac2*(2.0/beta)

                   !Ebpar = int B_par^2/(8 pi)
                   hk(it,ik) % ebpar = hk(it,ik) % ebpar &
                        + 0.5*real(conjg(bpar_avg)*bpar_avg) * wgt(ig)*fac2*(2.0/beta)

                end do
                h % antenna = h % antenna + hk(it, ik) % antenna
                h % eapar = h % eapar + hk(it, ik) % eapar
                h % ebpar = h % ebpar + hk(it, ik) % ebpar
             end do
          end do
       else
          hk % antenna = 0.
          h  % antenna = 0.
          hk % energy_dot = 0.
          hk % energy = 0.
          hk % eapar = 0.
          h  % eapar = 0.
          hk % ebpar = 0.
          h  % ebpar = 0.
       end if
       deallocate (j_ext)
    end if

! ==========================================================================
! Finish E_dot--------------------------------------------------------------
! ==========================================================================

!GGH Include response of Boltzmann species for single-species runs

    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   phi_m   = fdot (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik), dtinv)

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        fphi * real(conjg(phi_avg)*phi_m) &
                        * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to E_dot for single species runs---------------------

!GGH New E_dot calc
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1

             !Calculate old fluctuating energy de
             havg = favg_x (g(ig  ,isgn,iglo), &
                            g(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phi(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phi(ig  ,it,ik), &
                               phi(ig+1,it,ik)) * fphi * spec(is)%zt

             de=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)

            !Calculate new fluctuating energy denew
             havg = favg_x (gnew(ig  ,isgn,iglo), &
                            gnew(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phinew(ig  ,it,ik), &
                               phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             denew=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)

             !Set g_work as the change of energy (denew-de)/dt
             g_work(ig,isgn,iglo) = fdot_t(de,denew,dtinv)

          end do
       end do
    end do
    !GGH -END new e_dot calc

   call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy_dot = h % energy_dot + hk(it,ik) % energy_dot
          end do
       end do
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ==========================================================================
! Gradient Contributions to Heating-----------------------------------------
! ==========================================================================

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle

       do isgn=1,2
          do ig=-ntgrid, ntgrid-1

             chi = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phi   (ig+1,it,ik),  &
                         aj0(ig  ,iglo)*phinew(ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phinew(ig+1,it,ik)) &
                         *fphi

!!GGH Bug fix: The apar part should be subtracted (because chi= phi - v|| A|| + B||)
             chi = chi - &
                  favg (aj0(ig  ,iglo)*apar   (ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*apar   (ig+1,it,ik)*vpac(ig+1,isgn,iglo),  &
                        aj0(ig  ,iglo)*aparnew(ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*aparnew(ig+1,it,ik)*vpac(ig+1,isgn,iglo)) &
                        *spec(is)%stm*fapar

             chi = chi + &
                  favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                        aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                        *spec(is)%tz*fbpar

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

#ifdef LOWFLOW
             g_work(ig,isgn,iglo) = zi * wstar(ik,ie,is)*hneoc(ig,isgn,iglo)/code_dt * conjg(havg)*chi &
#else
             g_work(ig,isgn,iglo) = zi * wstar(ik,ie,is)/code_dt * conjg(havg)*chi &
#endif
                  * spec(is)%dens

          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % gradients(is) = hk(it,ik) % gradients(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % gradients(is) = h % gradients(is) + hk(it,ik) % gradients(is)
             end do
          end do
       end do
    end if
! ==========================================================================
! Hyperviscosity------------------------------------------------------------
! ==========================================================================

    if (D_v > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
!          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1

                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo))

                j0phiavg = favg (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phi(ig+1,it,ik), &
                                 aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                j1bparavg= favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                                 aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                                 *fbpar

                dgdt_hypervisc = 0.5*((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo) &
                     + (1.0-1./hypervisc_filter(ig+1,it,ik))*gnew(ig+1,isgn,iglo))/code_dt

!Set g_work for hyperviscous heating
!                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_v*akperp4* &
!                     ( conjg(havg)*havg - conjg(havg)*j0phiavg - &
!                     conjg(havg)*j1bparavg)
                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*conjg(havg)*dgdt_hypervisc

             end do
          end do
       end do

       call integrate_moment (g_work, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hypervisc(is) = hk(it,ik) % hypervisc(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hypervisc(is) = h % hypervisc(is) + hk(it,ik) % hypervisc(is)
                end do
             end do
          end do
       end if

    end if !End Hyperviscous Heating Calculation


! ==========================================================================
! Hyperresistivity------------------------------------------------------------
! ==========================================================================

    if (D_eta > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1

                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo))

                j0aparavg = favg (aj0(ig  ,iglo) * apar(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * apar(ig+1,it,ik), &
                                 aj0(ig  ,iglo)  * aparnew(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * aparnew(ig+1,it,ik)) &
                                 * fapar * spec(is)%zstm * vpac(ig,isgn,iglo)

!Set g_work for hyperresistive heating
                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_eta*akperp4* &
                     conjg(havg)*j0aparavg

             end do
          end do
       end do

       call integrate_moment (g_work, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hyperres(is) = hk(it,ik) % hyperres(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hyperres(is) = h % hyperres(is) + hk(it,ik) % hyperres(is)
                end do
             end do
          end do
       end if

    end if !End Hyperresistivity Heating Calculation

!==========================================================================
!Finish Energy-------------------------------------------------------------
!==========================================================================

!GGH Calculate hs2-------------------------------------------------------------
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g_work(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg)
          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                end do
             end do
             h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
          end do
       end do
    end if

!Calculate phis2-------------------------------------------------------------
    if (proc0) then
       do ik=1,naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1,ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % phis2(is) = hk(it,ik) % phis2(is)  &
                        +0.5*spec(is)%temp*spec(is)%dens*real(conjg(phi_avg)*phi_avg) &
                        * wgt(ig) * fac2
                enddo
             end do
             h % phis2(:) = h % phis2(:) + hk(it,ik) % phis2(:)
          end do
       end do
    endif

! Calculate delfs2 (rest of energy)-----------------------------------------------

!GGH  Include response of Boltzmann species for single species runs
    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy = hk(it,ik) % energy + &
                        fphi * real(conjg(phi_avg)*phi_avg) &
                        * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to energy for single species runs---------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             j0phiavg = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik), &
                              aj0(ig+1,iglo)*phi   (ig+1,it,ik), &
                              aj0(ig  ,iglo)*phinew(ig  ,it,ik), &
                              aj0(ig+1,iglo)*phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg (phi   (ig  ,it,ik), &
                             phi   (ig+1,it,ik), &
                             phinew(ig  ,it,ik), &
                             phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             g_work(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)
          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy = hk(it,ik) % energy  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                   !Delfs2 = int_r int_v T/F0 dfs^2/2
                   hk(it,ik) % delfs2(is) = hk(it,ik) % delfs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy = h % energy + hk(it,ik) % energy
             h % delfs2(:) = h % delfs2(:) + hk(it,ik) % delfs2(:)
          end do
       end do
       deallocate (wgt)
    end if

    deallocate (tot)

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    direction = to_g_gs2)
    call g_adjust (gnew, phinew, bparnew, direction = to_g_gs2)

  end subroutine get_heat


  !> Get a theta-centered and time-centered estimate of the time derivative
  !! of a field.
  !!
  !! tunits(ky) == 1. unless the "omega_*" units are chosen.
  !! omega_* units normalize time by an additional factor of ky.
  pure subroutine dot (a, anew, adot, fac)
    use gs2_time, only: code_dt, tunits
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid

    implicit none
    complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
    complex, intent (out), dimension (-ntgrid:,:,:) :: adot
    real, intent (in) :: fac
    real :: dtinv
    integer :: ig, it, ik

    do ik=1,naky
       dtinv = 1./(code_dt*tunits(ik))
       do it=1,ntheta0
          do ig=-ntgrid,ntgrid-1
             adot(ig,it,ik) = 0.5*fac*(anew(ig+1,it,ik)+anew(ig,it,ik) - &
                  (a(ig+1,it,ik)+a(ig,it,ik)))*dtinv
          end do
       end do
    end do

  end subroutine dot

  !> FIXME : Add documentation
  pure complex function fdot (fl, fr, fnewl, fnewr, dtinv)
    complex, intent (in) :: fl, fr, fnewl, fnewr
    real, intent (in) :: dtinv
    fdot = 0.5*(fnewl+fnewr-fl-fr)*dtinv
  end function fdot

  !> Construct time and space-centered quantities
  !! (should use actual bakdif and fexpr values?)
  pure complex function favg (fl, fr, fnewl, fnewr)
    complex, intent (in) :: fl, fr, fnewl, fnewr
    favg = 0.25*(fnewl+fnewr+fl+fr)
  end function favg

  !> FIXME : Add documentation
  pure complex function fdot_t (f,fnew, dtinv)
    complex, intent (in) :: f, fnew
    real, intent (in) :: dtinv
    fdot_t = (fnew-f)*dtinv
  end function fdot_t

  !> FIXME : Add documentation
  pure complex function favg_x (fl, fr)
    complex, intent (in) :: fl, fr
    favg_x = 0.5*(fl+fr)
  end function favg_x
end module gs2_heating
