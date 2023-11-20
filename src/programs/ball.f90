!> FIXME : Add documentation
program ball
  use geometry, only: eikcoefs, eikcoefs_output
  !> Inputs
  use geometry, only: nperiod, geoType, rmaj, shift, eqfile
  use geometry, only: equal_arc, rhoc, itor, qinp, iflux, delrho
  use geometry, only: bishop, irho, isym, efit_eq, writelots, R_geo
  use geometry, only: gen_eq, ppl_eq, local_eq, chs_eq, s_hat_input
  use geometry, only: p_prime_input, invLp_input, beta_prime_input, big, gs2d_eq
  use geometry, only: transp_eq, alpha_input, dp_mult
  use geometry, only: deltam, deltan, deltampri, deltanpri
  implicit none

  real, dimension(:), pointer :: psi

  real, dimension(:), allocatable :: theta
  real :: dbetadrho
  integer :: nbeta, ntgrid, ntheta, iunstable, ibeta, ishat
  integer :: nshat, order

  real :: beta_p1, beta_p2, diffscheme, psiend, beta_prime, cflmax, cflmin
  real :: beta_prime_times, beta_prime_over
  real :: bprev, pprev, stabmhd, bcrit, s_hat_max, shat_0
  real :: akappa, tri, akappri, tripri, rmin, rmax
  logical :: diagram
  type(eikcoefs_output) :: eikcoefs_results

  namelist/stuff/ntheta,nperiod,geoType,rmaj,akappri,akappa,shift,equal_arc, &
       rhoc,rmin,rmax,itor,qinp,iflux,delrho,tri,bishop, irho,isym, &
       tripri,writelots,R_geo, gen_eq,efit_eq,local_eq,eqfile, &
       chs_eq, s_hat_input,p_prime_input,invLp_input,ppl_eq,&
       diffscheme,nbeta,beta_p1,beta_p2,beta_prime_input, &
       alpha_input, beta_prime_times, beta_prime_over, big, gs2d_eq, &
       transp_eq, diagram, nshat, s_hat_max


! set defaults

  ntheta=32; nperiod=1
  equal_arc = .false. ; bishop = 0
  itor=1; iflux=0

  eqfile='dskeq.cdf'
  rhoc=0.5

  rmin = -1.0 ; rmax = -1.0

  diagram = .false.
  nshat = 1 ; s_hat_max = 2.0

  geoType = 0

  rmaj = 3   ;     akappa = 1   ;     tri = 0   ;     qinp = 2
  R_geo = 3  ;     akappri = 0  ;     tripri=0  ;    shift = 0.
   
  s_hat_input = 1
  p_prime_input = 1
  invLp_input = 5

  beta_prime_times = -1. ; beta_prime_over = -1.

  irho = 2 ;      isym = 0

  delrho = 0.01

  diffscheme=0.33; nbeta = 1; beta_p1 = 0.; beta_p2 = 0.
      
!     read in data
  open(10,file='eik.in')
  read(10,stuff)
  close(10)

  if (rmin /= -1) write(*,*) "Option deprecated : rmin set but not used"
  if (rmax /= -1) write(*,*) "Option deprecated : rmax set but not used"

  writelots = .false.

  ! Copy shaping inputs to generic names
  deltam = akappa
  deltan = tri
  deltampri = akappri
  deltanpri = tripri

  if(.not.local_eq) then
     if(beta_prime_times >= 0) then
        if(beta_prime_over == -1) beta_prime_over = beta_prime_times
     endif     
  endif

  equal_arc = .false.   ! this should not be required
  
  !if(k1.lt.0) k1=(ntheta/abs(k1))**2
  !if(k2.lt.0) k2=(ntheta/abs(k2))**2
  
  !if(vmom_eq .and. gen_eq) then
     !write(*,*) 'Choose either vmom_eq or gen_eq, not both'               
     !write(*,*) 'Stopping'
     !stop
  !endif
  
!
! Note that if iflux=1 then always choose itor=1
!
  if(iflux.eq.1 .and. itor.ne.1) then
     itor=1
     write(*,*) 'Forcing itor=1'
  endif
  
  if(iflux /= 0 .and. iflux /=0) then
     write(*,*) 'iflux=0 or 1 is the only option.'
     stop
  endif
!
! Catch if no q profile 
!
  if(iflux.ne.1 .and. irho.eq.1) then
     irho=2
     write(*,*) 'Forcing irho=2'
  endif
  
  if(iflux.ne.1) then
     !if(vmom_eq) write(*,*) 'Forcing vmom_eq to be false'
     if(gen_eq) write(*,*) 'Forcing gen_eq to be false'
     if(gen_eq) write(*,*) 'because iflux.ne.1'
     !vmom_eq=.false.
     gen_eq=.false.
  endif
  
!     compute the theta grid

  if( (.not. gen_eq) .and. (.not. ppl_eq) &
       .and. (.not. transp_eq) .and. (.not. chs_eq)) then
     call eikcoefs(ntheta, eikcoefs_results)
     theta = eikcoefs_results%theta
     dbetadrho = eikcoefs_results%dbetadrho
     ntgrid = (2*nperiod - 1)*(ntheta/2)
     allocate(psi(-ntgrid:ntgrid))
  else 
     allocate(psi(-1:1))
  endif 

  open(11,file='ball.out')
  open(29,file='psi.out')

  if (.not. diagram) then
     order = 1
     call beta_loop (bcrit, order, ntheta)
     stop
  endif

! Attempt to sketch out an s-alpha type diagram

  if (bishop /= 8) stop  ! just do one case for now

  shat_0 = s_hat_input

  do order = 1, 2
     do ishat = 1, nshat
        if (nshat > 1) &
             s_hat_input = shat_0 + (ishat-1)*(s_hat_max-shat_0)/(nshat-1)
        
        bcrit = 0.
        call beta_loop (bcrit, order, ntheta)
        
        if (bcrit /= 0.) then
           write (*,*) 'shat= ',s_hat_input,' alpha= ', -Rmaj*qinp**2*bcrit
        end if
     end do
     write (*,*) 
  end do

contains
  
  subroutine beta_loop (bcrit, order, ntheta)
    integer, intent(inout) :: ntheta
    integer :: order
    real :: bcrit
    logical :: done 
    integer :: ibmin, ibmax, ibstep

    done = .false.

    if (order == 1) then
       ibmin = 1
       ibmax = nbeta
       ibstep = 1
    else
       ibmin = nbeta
       ibmax = 1
       ibstep = -1
    end if

    do ibeta=ibmin, ibmax, ibstep
       
       if(bishop >= 7) then
          if(nbeta == 1) then
             dp_mult = beta_prime_times
          else
             dp_mult = 1./beta_prime_over+(ibeta-1) &
                  *(beta_prime_times-1./beta_prime_over)/(nbeta-1)
          endif
       else
          if(nbeta == 1) then
             beta_prime=beta_p1
          else
             beta_prime=beta_p1+(ibeta-1)*(beta_p2-beta_p1)/(nbeta-1)
          endif
          beta_prime_input = beta_prime
       endif
       
       call eikcoefs(ntheta, eikcoefs_results)
       theta = eikcoefs_results%theta
       dbetadrho = eikcoefs_results%dbetadrho
       
       if(.not.local_eq) then
          ntgrid = (size(theta)-1)/2
          ntheta = (size(theta)-1)/(2*nperiod-1)
          deallocate(psi)
          allocate(psi(-ntgrid:ntgrid))
       endif
       
       if(bishop>= 7) beta_prime = dbetadrho
       
       call mhdballn(ntgrid, beta_prime, diffscheme, iunstable, psiend, &
            cflmax, cflmin, psi, eikcoefs_results%bmag, eikcoefs_results%cvdrift, &
            eikcoefs_results%gradpar, eikcoefs_results%gds2, theta)
       
       if(ibeta.eq.1) then
          stabmhd=0.
       else
          stabmhd=0.
          if(psiend /= pprev) then
             stabmhd=beta_prime-psiend*(beta_prime-bprev)/(psiend-pprev)
          end if

          if (diagram) then
             if (iunstable == 1 .and. .not. done) then
                bcrit = stabmhd
                done = .true.
             end if
          end if

       endif
       bprev=beta_prime
       pprev=psiend

       if (done) cycle
       
       if (.not. diagram) then
        write(6,998) beta_prime,iunstable,stabmhd,max(abs(cflmax),abs(cflmin)),dp_mult
       end if
    enddo

998  format(1x,g13.6,2x,i2,3x,g13.6,2x,2(1pe11.4,1x))
  end subroutine beta_loop
  
  subroutine mhdballn(ntgrid, beta_prime, diffscheme, iunstable, psiend, cflmax, &
       cflmin, psi, bmag, cvdrift, gradpar, gds2, theta)
    implicit none
    integer, intent(in) :: ntgrid
    real, intent(in) :: beta_prime, diffscheme
    integer, intent(out) :: iunstable
    real, intent(out) :: psiend, cflmax, cflmin
    real, dimension(-ntgrid:ntgrid), intent(out) :: psi
    real, dimension(-ntgrid:ntgrid), intent(in) :: bmag, cvdrift, gradpar, gds2, theta
    integer :: i
    real, dimension(-ntgrid:ntgrid) :: psi_p, c, g, c1, c2, ch, g1, g2, gh, delx
    real :: alpha,  b_prime

    !       include 'ball.inc'

    !       input arguments:
    !       npt=total number of gridpoints with theta>0 = nperiod*ntheta+1
    !       beta_prime=d beta/ d rho
    !       alpha=numerical knob (0 to .333 are good if cflmax,min<1)
    !
    !       output arguments:
    !       iunstable=number of times the solution crossed the axis. If
    !                 this is greater than 0 it is unstable.
    !
    !       psi is the solution of the Euler Lagrange equation for delta W.
    !       At the first theta point, psi is zero with slope positive 1.
    !
    !       psiend=psi at the last gridpoint
    !       psi is psi at every grid point
    !       cflmax,min are indicators of whether resolution is adequate.
    !       diffscheme is a numerical difference scheme knob.

    !       For greatest accuracy diffscheme=0. or .333333 .
    !       If the absolute value of cflmax,min approaches one, inaccuracy and/or
    !       numerical instability are possible.  If cflmax,min are >1-3, increase
    !       diffscheme to .666 or 1.0 for stability, but probably resolution is
    !       inadequate.
    !

    b_prime = -beta_prime
    iunstable=0
    alpha=1.-diffscheme

    g = gds2*gradpar/bmag
    c = b_prime*cvdrift/(2.*bmag*gradpar)

    i=-ntgrid
    cflmax=(theta(i+1)-theta(i))**2*(c(i+1)+c(i))/(g(i+1)+g(i))
    cflmin=cflmax
    do i=-ntgrid+1, ntgrid
       delx(i)=theta(i)-theta(i-1)
       ch(i)=.5*(c(i)+c(i-1))
       gh(i)=.5*(g(i)+g(i-1))
       cflmax=max(cflmax,delx(i)**2*ch(i)/gh(i))
       cflmin=min(cflmin,delx(i)**2*ch(i)/gh(i))
    enddo


    do i=-ntgrid+1, ntgrid-1
       c1(i)=-delx(i+1)*(alpha*c(i)+.5*diffscheme*ch(i+1)) &
            -delx(i  )*(alpha*c(i)+.5*diffscheme*ch(i  )) &
            -delx(i)*.5*diffscheme*ch(i)
       c1(i)=.5*c1(i)
    enddo

    do i=-ntgrid+1, ntgrid
       c2(i)=-.25*diffscheme*ch(i)*delx(i)
       g1(i)=gh(i)/delx(i)
       g2(i)=1./(gh(i)/delx(i)+.25*diffscheme*ch(i)*delx(i))
    enddo

    i=-ntgrid
    psi(i)=0.
    i=i+1
    psi(i)=delx(i)
    psi_p(i)=psi(i)/g2(i)-g1(i)*psi(i-1)

    do i=-ntgrid+1,ntgrid-1
       psi_p(i+1)=psi_p(i)+c1(i)*psi(i)+c2(i)*psi(i-1)
       psi(i+1)=(g1(i+1)*psi(i)+psi_p(i+1))*g2(i+1)
    enddo
    psiend=psi(ntgrid)

    do i=-ntgrid+1,ntgrid-1
       if(psi(i)*psi(i+1) <= 0.) iunstable=iunstable+1
    enddo

  end subroutine mhdballn
  
end program ball

