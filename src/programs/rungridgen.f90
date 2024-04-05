!> FIXME : Add documentation  
program rungridgen
  use file_utils
  use gridgen4mod, only: gridgen4_1
  use constants, only: pi, twopi
  use mp, only: init_mp, finish_mp
  implicit none

  ! input parameters
  character(200) :: source, gingrid, gsource
  integer :: nthetaout, nlambdaout, nperiodout
  integer :: npadd, iperiod
  real :: alknob, epsknob, extrknob, bpknob, smoothknob
  integer :: nfinegrid
  real :: thetamax, deltaw, widthw, tension
  logical :: screenout
  logical :: auto_width, three_dim
  real :: cv_fraction, delth_max
  integer :: max_autoiter
  logical :: list = .false.
  logical :: no_end_point
  character(2000) :: ncsource

  ! work variables
  integer :: ntheta, nlambda, ntgrid
  integer :: ntgridin, nperiodin, nthetain
  real, dimension (:), allocatable :: thetain, bmagin, bmagsm
  real, dimension (:), allocatable :: thetaout, bmagout, alambdaout

  real :: drhodpsi, rmaj, shat, kxfac, qval
  real, dimension (:), allocatable :: thetagrid
  real, dimension (:), allocatable :: gbdrift, gradpar, grho
  real, dimension (:), allocatable :: cvdrift, gds2, bmag
  real, dimension (:), allocatable :: gds21, gds22
  real, dimension (:), allocatable :: cvdrift0, gbdrift0
  real, dimension (:), allocatable :: Rplot, Rprime
  real, dimension (:), allocatable :: Zplot, Zprime
  real, dimension (:), allocatable :: aplot, aprime
  
  call init_mp
  call init_file_utils (list, name="grid")
  call read_parameters
  if (three_dim) then
     call allocate_arrays_3d
     call get_initial_grids_3d
  else
     call allocate_arrays
     call get_initial_grids
  end if
  call do_smoothing
  call generate_grids
  call write_output
  call finish_file_utils
  call finish_mp
  stop

contains
  !> FIXME : Add documentation  
  subroutine read_parameters
    implicit none
    namelist /testgridgen/ source, gsource, &
         nthetaout,nlambdaout,nperiodout, &
         npadd,alknob,epsknob,bpknob,extrknob,smoothknob, nfinegrid, &
         thetamax, deltaw, widthw, tension, gingrid, screenout, &
         auto_width, cv_fraction, delth_max, max_autoiter, three_dim, &
         iperiod, &
         no_end_point, ncsource
    gsource = "eik6.out"
    source = "eik.out"
    nthetaout = 32
    nlambdaout = 20
    nperiodout = 2
    npadd = 2
    alknob = 0.1
    bpknob = 1.e-8
    epsknob = 1e-9
    extrknob = 10.0
    smoothknob = 0.0
    nfinegrid=200
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    tension = 1.0
    screenout = .false.
    auto_width = .false.
    cv_fraction = 0.6
    delth_max = 0.5
    max_autoiter = 3
    gingrid = "gingrid"
    three_dim = .false.
    iperiod = 1
    no_end_point = .false.
    ncsource = ""
    read (unit=input_unit("testgridgen"), nml=testgridgen)
    nthetaout = nthetaout + 1
  end subroutine read_parameters

  !> FIXME : Add documentation  
  subroutine allocate_arrays
    use gs2_io_grids, only:nc_grid_file_open
    use gs2_io_grids, only: nc_get_grid_sizes, nc_get_grid_scalars
    implicit none
    integer :: unit
    character(200) :: line
    logical :: ncsource_exist

    inquire(file=ncsource, exist=ncsource_exist)
    if (ncsource_exist) then
       call nc_grid_file_open(ncsource, "r")
       call nc_get_grid_sizes(nthetain, ntgridin, nperiodin)
       call nc_get_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj)
    else
       call get_unused_unit(unit)
       open (unit=unit, file=trim(source), status="old")
       read (unit=unit, fmt="(a)") line
       read (unit=unit, fmt=*,err=10) ntgridin, nperiodin, nthetain, drhodpsi, rmaj, shat, kxfac, qval
10     continue
       close (unit=unit)
    end if

    allocate (thetain(nthetain+1),bmagin(nthetain+1),bmagsm(nthetain+1))
    allocate (thetaout(nthetaout),bmagout(nthetaout),alambdaout(nlambdaout))
    allocate (thetagrid(-ntgridin:ntgridin))
    allocate (gbdrift(-ntgridin:ntgridin))
    allocate (gradpar(-ntgridin:ntgridin))
    allocate (grho(-ntgridin:ntgridin))
    allocate (cvdrift(-ntgridin:ntgridin))
    allocate (gds2(-ntgridin:ntgridin))
    allocate (bmag(-ntgridin:ntgridin))
    allocate (gds21(-ntgridin:ntgridin))
    allocate (gds22(-ntgridin:ntgridin))
    allocate (cvdrift0(-ntgridin:ntgridin))
    allocate (gbdrift0(-ntgridin:ntgridin))
    allocate (Rplot(-ntgridin:ntgridin))
    allocate (Rprime(-ntgridin:ntgridin))
    allocate (Zplot(-ntgridin:ntgridin))
    allocate (Zprime(-ntgridin:ntgridin))
    allocate (aplot(-ntgridin:ntgridin))
    allocate (aprime(-ntgridin:ntgridin))
  end subroutine allocate_arrays

  !> FIXME : Add documentation    
  subroutine allocate_arrays_3d
    implicit none
    integer :: unit, ntor
!    character(200) :: line

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old")
    read (unit=unit, fmt=*) ntgridin, nthetain, ntor
    close (unit=unit)
    
    drhodpsi = 1.0  ! assumes our radial variable matches VVBAL

    nthetain = nthetain * iperiod

    allocate (thetain(nthetain+1),bmagin(nthetain+1),bmagsm(nthetain+1))
    allocate (thetaout(nthetaout),bmagout(nthetaout),alambdaout(nlambdaout))
    allocate (thetagrid(-ntgridin:ntgridin))
    allocate (gbdrift(-ntgridin:ntgridin))
    allocate (gradpar(-ntgridin:ntgridin))
    allocate (grho(-ntgridin:ntgridin))
    allocate (cvdrift(-ntgridin:ntgridin))
    allocate (gds2(-ntgridin:ntgridin))
    allocate (bmag(-ntgridin:ntgridin))
    allocate (gds21(-ntgridin:ntgridin))
    allocate (gds22(-ntgridin:ntgridin))
    allocate (cvdrift0(-ntgridin:ntgridin))
    allocate (gbdrift0(-ntgridin:ntgridin))
    allocate (Rplot(-ntgridin:ntgridin))
    allocate (Rprime(-ntgridin:ntgridin))
    allocate (Zplot(-ntgridin:ntgridin))
    allocate (Zprime(-ntgridin:ntgridin))
    allocate (aplot(-ntgridin:ntgridin))
    allocate (aprime(-ntgridin:ntgridin))

    cvdrift0 = 0.   ! assumes theta0 = 0.
    gbdrift0 = 0.   ! assumes theta0 = 0.
    gds21 = 0.      ! assumes theta0 = 0.
    gds22 = 0.      ! assumes theta0 = 0.
    grho  = 1.      ! assumes linear calculation
    Rplot = 1.      ! pure fiction
    Rprime = 1.     ! pure fiction
    Zplot = 1.      ! pure fiction
    Zprime = 1.     ! pure fiction
    aplot = 1.      ! pure fiction
    aprime = 1.     ! pure fiction
    
  end subroutine allocate_arrays_3d

  !> FIXME : Add documentation  
  subroutine get_initial_grids
    use file_utils, only: get_unused_unit
    use gs2_io_grids, only: nc_get_grids
    use gs2_io_grids, only: nc_grid_file_close
    implicit none
    integer :: unit
    integer :: i
!    real :: discard
    character(200) :: line
    logical :: ncsource_exist

    inquire(file=ncsource, exist=ncsource_exist)
    if (ncsource_exist) then
       call nc_get_grids(ntgridin, &
            bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
            gds2, gds21, gds22, grho, thetagrid, &
            Rplot=Rplot, Rprime=Rprime, Zplot=Zplot, Zprime=Zprime, aplot=aplot, aprime=aprime)
       call nc_grid_file_close()
    else
       call get_unused_unit(unit)
       open (unit=unit, file=trim(source), status="old")
       read (unit=unit, fmt="(a)") line
       read (unit=unit, fmt="(a)") line

       ! gbdrift gradpar grho theta
       read (unit=unit, fmt="(a)") line
       do i = -ntgridin,ntgridin
          read (unit=unit, fmt=*) gbdrift(i),gradpar(i),grho(i),thetagrid(i)
       end do

       ! cvdrift gds2 bmag theta
       read (unit=unit, fmt="(a)") line
       do i = -ntgridin,ntgridin
          read (unit=unit, fmt=*) cvdrift(i),gds2(i),bmag(i),thetagrid(i)
       end do

       ! gds21 gds22 theta
       read (unit=unit, fmt="(a)") line
       do i = -ntgridin,ntgridin
          read (unit=unit, fmt=*) gds21(i),gds22(i),thetagrid(i)
       end do

       ! cvdrift0 gbdrift0 theta
       read (unit=unit, fmt="(a)") line
       do i = -ntgridin,ntgridin
          read (unit=unit, fmt=*) cvdrift0(i),gbdrift0(i),thetagrid(i)
       end do

       close (unit=unit)

       call get_plotdata
    end if

    thetain = thetagrid(-nthetain/2:nthetain/2)
    bmagin = bmag(-nthetain/2:nthetain/2)

  end subroutine get_initial_grids

  subroutine get_plotdata
    implicit none
    integer :: unit
    integer :: i
    character(200) :: line
    real, allocatable :: dummy(:)

    allocate (dummy(-ntgridin:ntgridin))

    Rplot = 0.; Rprime = 0.
    Zplot = 0.; Zprime = 0.
    aplot = 0.; aprime = 0.

    call get_unused_unit(unit)
    open (unit=unit, file=trim(gsource), status="old", err=100)

    read (unit=unit, fmt="(a)", err=100) line
    read (unit=unit, fmt="(a)", err=100) line
    read (unit=unit, fmt="(a)", err=100) line

    do i = -ntgridin, ntgridin
       read (unit=unit, fmt=*, err=100) dummy(i), Rplot(i), Zplot(i), aplot(i), Rprime(i), Zprime(i), aprime(i), dummy(i)
    end do

    close(unit=unit)

    return

100 continue
    close(unit=unit)
    write(6,*) 'read error: ', gsource

    return

  end subroutine get_plotdata

  !> FIXME : Add documentation  
  subroutine get_initial_grids_3d
    use file_utils, only: get_unused_unit
    implicit none
    integer :: unit
    integer :: i, ntor
    real :: dpdpsi, pres, bavg, cvavg

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old", action="read")

    read (unit=unit, fmt=*) ntgridin, nthetain, ntor

    nthetain = nthetain * iperiod
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) thetagrid(i),bmag(i),gradpar(i),gds2(i),cvdrift(i),gbdrift(i),dpdpsi,pres
    end do

    close (unit=unit)

    bavg = sum(bmag(-nthetain/2:nthetain/2))/(nthetain+1)
    cvavg = sum(cvdrift)/size(cvdrift)
    
    write(*,*) 'bavg = ',bavg
    write(*,*) 'cvavg = ',cvavg,' and remember that cvdrift assumed beta_prime=0'

!    bmag = bmag / bavg

    gds2 = gds2 * ntor**2 
    gbdrift = gbdrift * ntor * 2. / bmag**2
!    cvdrift = cvdrift * ntor * 2. / bmag**3    ! not correct for beta_prime term, so use: 
    cvdrift = gbdrift 

!    alp = -1./pres*dpdpsi

!
! these are theta and mod(b) grids that are periodic
!
    thetain = thetagrid(-nthetain/2:nthetain/2)
    bmagin = bmag(-nthetain/2:nthetain/2)
  end subroutine get_initial_grids_3d

  !> FIXME : Add documentation  
  subroutine do_smoothing
    implicit none
    real :: var
    integer :: ifail
    if (smoothknob == 0.0) then
       bmagsm = bmagin
    else
       var = smoothknob
       ifail = 0
       call smooth (nthetain+1,thetain,bmagin,var,bmagsm,ifail)
       if (ifail /= 0) then
          print *, "smooth failed: ",ifail
          select case (ifail)
          case (129)
             print *, "IC IS LESS THAN N-1."
          case (130)
             print *, "N IS lESS THAN 3."
          case (131)
             print *, "INPUT ABSCISSAE ARE NOT ORDERED SO THAT X(I).LT.X(I+1)."
          case (132)
             print *, "DF(I) IS NOT POSITIVE FOR SOME I."
          case (133)
             print *, "JOB IS NOT 0 OR 1."
          case default
          end select
          stop
       end if
    end if
  end subroutine do_smoothing

  !> FIXME : Add documentation  
  subroutine generate_grids
    implicit none
    integer :: i, iter, nin
    real :: d, deltaw_ok

    if (.not. auto_width) then
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       return
    end if

    widthw = pi
    do i = 0, nthetain/2
       if (cvdrift(i) < 0.0) then
          widthw = thetagrid(i)
          exit
       end if
    end do

    print *, "widthw: ", widthw
    deltaw_ok = 0.0
    do iter = 1, max_autoiter
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       print *, "iter: ", iter
       d = maxval(thetaout(2:ntheta) - thetaout(1:ntheta-1))
       print *, "max deltheta: ", d
       nin = count(abs(thetaout(1:ntheta) - thetamax) < widthw &
                   .or. abs(thetaout(1:ntheta) + thetamax) < widthw)
       print *, "count in +ve cvdrift: ", nin
       print *, "fraction in +ve cvdrift: ", real(nin)/real(ntheta)
       print *, "deltaw: ", deltaw
       if (d > delth_max) then
          if (deltaw_ok /= 0.0) then
             deltaw = sqrt(deltaw*deltaw_ok)
          else
             deltaw = 0.5*deltaw
          end if
          print *, "max deltheta too large: ", d
          if (deltaw < deltaw_ok) exit
       else if (real(nin)/real(ntheta) < cv_fraction) then
          print *, "fraction in +ve cvdrift too small: ", &
               real(nin)/real(ntheta)
          deltaw_ok = deltaw
          deltaw = deltaw*(cv_fraction/(real(nin)/real(ntheta)))**2
       end if
    end do

    if (deltaw_ok /= 0.0 .and. d > delth_max) then
       deltaw = deltaw_ok
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       d = maxval(thetaout(2:ntheta) - thetaout(1:ntheta-1))
       print *, "max deltheta: ", d
       nin = count(abs(thetaout(1:ntheta) - thetamax) < widthw &
                   .or. abs(thetaout(1:ntheta) + thetamax) < widthw)
       print *, "count in +ve cvdrift: ", nin
       print *, "fraction in +ve cvdrift: ", real(nin)/real(ntheta)
       print *, "deltaw: ", deltaw
    end if

  end subroutine generate_grids

  !> FIXME : Add documentation    
  subroutine write_output
    use splines, only: inter_cspl
    use gs2_io_grids, only: nc_write_grid_sizes, nc_write_grid_scalars, nc_write_grids
    use gs2_io_grids, only: nc_write_lambda_grid
    use gs2_io_grids, only: nc_grid_file_open, nc_grid_file_close
    use constants, only: run_name_size
    use mp, only: proc0
    implicit none
    integer :: unit
    integer :: i
!    real, allocatable, dimension (:) :: bmaginaux, bmagsmaux, tmp
    real, allocatable, dimension (:) :: thetagridout
    real, allocatable, dimension (:) :: gbdriftout, gradparout, grhoout
    real, allocatable, dimension (:) :: cvdriftout, gds2out, bmaggridout
    real, allocatable, dimension (:) :: gds21out, gds22out
    real, allocatable, dimension (:) :: cvdrift0out, gbdrift0out
    real, allocatable, dimension (:) :: Rplotout, Rprimeout
    real, allocatable, dimension (:) :: Zplotout, Zprimeout
    real, allocatable, dimension (:) :: aplotout, aprimeout
!    real :: th, bmin
    real :: bmin
!    integer :: ierr
    integer, dimension (1) :: minloca
    character(len=run_name_size) :: ncout

    call open_output_file (unit, ".input.out")
    write (unit=unit, fmt="('#',i5)") nthetain+1
    do i = 1, nthetain+1
       write (unit=unit, fmt="(3(1x,g19.12))") thetain(i), bmagin(i), bmagsm(i)
    enddo
    call close_output_file (unit)

!    allocate (bmaginaux(nthetain+1), bmagsmaux(nthetain+1), tmp(2*nthetain))
!    call fitp_curvp1 (nthetain,thetain,bmagin,twopi,bmaginaux,tmp,1.0,ierr)
!    call fitp_curvp1 (nthetain,thetain,bmagsm,twopi,bmagsmaux,tmp,1.0,ierr)
!    call open_output_file (unit, ".input.fine")
!    do i = -nfinegrid, nfinegrid
!       th = pi*real(i)/real(nfinegrid)
!       write (unit=unit, fmt="(3(x,g19.12))") th, &
!            fitp_curvp2(th,nthetain,thetain,bmagin,twopi,bmaginaux,1.0), &
!            fitp_curvp2(th,nthetain,thetain,bmagsm,twopi,bmagsmaux,1.0)
!    end do
!    call close_output_file (unit)
!    deallocate (bmaginaux, bmagsmaux, tmp)

    call open_output_file (unit, ".bmag.out")
    write (unit=unit, fmt="('#',i5)") ntheta
    do i = 1, ntheta
       write (unit=unit, fmt="(2(1x,g19.12))") thetaout(i), bmagout(i)
    enddo
    call close_output_file (unit)

    call open_output_file (unit, ".lambda.out")
    write (unit=unit, fmt="('#',i5)") nlambda
    do i = 1, nlambda
       write (unit=unit, fmt=*) alambdaout(i)
    enddo
    call close_output_file (unit)

    ntgrid = ntheta/2 + (nperiodout-1)*ntheta
    allocate (thetagridout(-ntgrid:ntgrid))
    allocate (gbdriftout(-ntgrid:ntgrid))
    allocate (gradparout(-ntgrid:ntgrid))
    allocate (grhoout(-ntgrid:ntgrid))
    allocate (cvdriftout(-ntgrid:ntgrid))
    allocate (gds2out(-ntgrid:ntgrid))
    allocate (bmaggridout(-ntgrid:ntgrid))
    allocate (gds21out(-ntgrid:ntgrid))
    allocate (gds22out(-ntgrid:ntgrid))
    allocate (cvdrift0out(-ntgrid:ntgrid))
    allocate (gbdrift0out(-ntgrid:ntgrid))
    allocate (Rplotout(-ntgrid:ntgrid))
    allocate (Rprimeout(-ntgrid:ntgrid))
    allocate (Zplotout(-ntgrid:ntgrid))
    allocate (Zprimeout(-ntgrid:ntgrid))
    allocate (aplotout(-ntgrid:ntgrid))
    allocate (aprimeout(-ntgrid:ntgrid))

    thetagridout(-ntheta/2:ntheta/2-1) = thetaout(:ntheta)
    bmaggridout(-ntheta/2:ntheta/2-1) = bmagout(:ntheta)
    thetagridout(ntheta/2) = thetagridout(-ntheta/2) + twopi*iperiod
    bmaggridout(ntheta/2) = bmaggridout(-ntheta/2)
    do i = 1, nperiodout - 1
       thetagridout(-ntheta/2-ntheta*i:ntheta/2-ntheta*i-1) &
            = thetagridout(-ntheta/2:ntheta/2-1) - twopi*i
       thetagridout(-ntheta/2+ntheta*i+1:ntheta/2+ntheta*i) &
            = thetagridout(-ntheta/2+1:ntheta/2) + twopi*i
       bmaggridout(-ntheta/2-ntheta*i:ntheta/2-ntheta*i-1) &
            = bmaggridout(-ntheta/2:ntheta/2-1)
       bmaggridout(-ntheta/2+ntheta*i+1:ntheta/2+ntheta*i) &
            = bmaggridout(-ntheta/2+1:ntheta/2)
    end do

    call inter_cspl(thetagrid, gbdrift, thetagridout, gbdriftout)
    call inter_cspl(thetagrid, gradpar, thetagridout, gradparout)
    call inter_cspl(thetagrid, grho, thetagridout, grhoout)
    call inter_cspl(thetagrid, cvdrift, thetagridout, cvdriftout)
    call inter_cspl(thetagrid, gds2, thetagridout, gds2out)
    call inter_cspl(thetagrid, bmag, thetagridout, bmagout)
    call inter_cspl(thetagrid, gds21, thetagridout, gds21out)
    call inter_cspl(thetagrid, gds22, thetagridout, gds22out)
    call inter_cspl(thetagrid, cvdrift0, thetagridout, cvdrift0out)
    call inter_cspl(thetagrid, gbdrift0, thetagridout, gbdrift0out)
    call inter_cspl(thetagrid, Rplot, thetagridout, Rplotout)
    call inter_cspl(thetagrid, Rprime, thetagridout, Rprimeout)
    call inter_cspl(thetagrid, Zplot, thetagridout, Zplotout)
    call inter_cspl(thetagrid, Zprime, thetagridout, Zprimeout)
    call inter_cspl(thetagrid, aplot, thetagridout, aplotout)
    call inter_cspl(thetagrid, aprime, thetagridout, aprimeout)

    call open_output_file (unit, ".out")
    write (unit=unit, fmt=*) 'nlambda'
    write (unit=unit, fmt=*) nlambda
    write (unit=unit, fmt=*) 'lambda'
    do i = 1, nlambda
       write (unit=unit, fmt=*) alambdaout(i)
    enddo

    write (unit=unit, fmt=*) 'ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q'
    write (unit=unit, fmt="(3i4,5(1x,g19.10))") ntgrid, nperiodout, ntheta, drhodpsi, rmaj, shat, kxfac, qval

    write (unit=unit, fmt=*) 'gbdrift gradpar grho tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            gbdriftout(i),gradparout(i),grhoout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'cvdrift gds2 bmag tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            cvdriftout(i),gds2out(i),bmaggridout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'gds21 gds22 tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            gds21out(i),gds22out(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'cvdrift0 gbdrift0 tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            cvdrift0out(i),gbdrift0out(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'Rplot Rprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            Rplotout(i),Rprimeout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'Zplot Zprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            Zplotout(i),Zprimeout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'aplot aprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            aplotout(i),aprimeout(i),thetagridout(i)
    end do

    call close_output_file (unit)

    if (screenout) then
       write (*, *) 'cvdrift       gds2          bmag          theta'
       do i = -ntheta/2,ntheta/2
          write (unit=*, fmt="(4(1x,g13.5))") &
               cvdriftout(i),gds2out(i),bmaggridout(i),thetagridout(i)
       end do
    end if

    minloca = minloc(bmagout(:ntheta))
    bmin = bmagout(minloca(1))
    print *, "theta,bmag minimum: ", thetaout(minloca(1)), bmin

    minloca = minloc(bmagin)
    print *, "theta,bmag input minimum: ", &
         thetain(minloca(1)), bmagin(minloca(1))

    if (bmin < bmagin(minloca(1))) print *, "WARNING: interpolated new minimum"

    if (proc0) then
       ncout = trim(run_name)//".out.nc"
       call nc_grid_file_open(ncout, "w")
       call nc_write_grid_sizes(ntheta, ntgrid, nperiodout)
       call nc_write_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj)
       call nc_write_grids(ntgrid, bmaggridout, gradparout, gbdriftout, gbdrift0out, &
            cvdriftout, cvdrift0out, gds2out, gds21out, gds22out, grhoout, thetagridout, &
            Rplot=Rplotout, Rprime=Rprimeout, Zplot=Zplotout, Zprime=Zprimeout, &
            aplot=aplotout, aprime=aprimeout, no_end_point_in=no_end_point)
       call nc_write_lambda_grid(nlambda, alambdaout(:nlambda))
       call nc_grid_file_close()
    end if

  end subroutine write_output

  !> FIXME : Add documentation
  subroutine smooth (n, xin, yin, var, yout, ifail)
    implicit none

    integer, intent(in) :: n
    real, dimension (:), intent (in) :: xin, yin
    real, dimension (:), intent (out) :: yout
    real, intent(in out) :: var
    integer, intent(out) :: ifail

! these next arrays should be double precision, which we usually get with
! compiler options
    real, dimension (1) :: se
    real, dimension (:), allocatable :: x, f, y, df
    real, dimension (:,:), allocatable :: wk, c
    real :: dvar

    allocate (x(n), f(n), y(n), df(n), c(n, 3))
    allocate (wk(n+2,7))

    ifail = 0
    x = xin
    f = yin
    df = 1.
    dvar=var
    call cubgcv (x,f,df,n,y,c,n,dvar,0,se,wk,ifail)
    var=dvar
    yout = y

    deallocate (x, f, y, df, c, wk)

  end subroutine smooth

!     algorithm 642 collected algorithms from acm.
!     algorithm appeared in acm-trans. math. software, vol.12, no. 2,
!     jun., 1986, p. 150.
!   subroutine name     - cubgcv
!
!--------------------------------------------------------------------------
!
!   computer            - vax/double
!
!   author              - m.f.hutchinson
!                         csiro division of mathematics and statistics
!                         p.o. box 1965
!                         canberra, act 2601
!                         australia
!
!   latest revision     - 15 august 1985
!
!   purpose             - cubic spline data smoother
!
!   usage               - call cubgcv (x,f,df,n,y,c,ic,var,job,se,wk,ier)
!
!   arguments    x      - vector of length n containing the
!                           abscissae of the n data points
!                           (x(i),f(i)) i=1..n. (input) x
!                           must be ordered so that
!                           x(i) .lt. x(i+1).
!                f      - vector of length n containing the
!                           ordinates (or function values)
!                           of the n data points (input).
!                df     - vector of length n. (input/output)
!                           df(i) is the relative standard deviation
!                           of the error associated with data point i.
!                           each df(i) must be positive.  the values in
!                           df are scaled by the subroutine so that
!                           their mean square value is 1, and unscaled
!                           again on normal exit.
!                           the mean square value of the df(i) is returned
!                           in wk(7) on normal exit.
!                           if the absolute standard deviations are known,
!                           these should be provided in df and the error
!                           variance parameter var (see below) should then
!                           be set to 1.
!                           if the relative standard deviations are unknown,
!                           set each df(i)=1.
!                n      - number of data points (input).
!                           n must be .ge. 3.
!                y,c    - spline coefficients. (output) y
!                           is a vector of length n. c is
!                           an n-1 by 3 matrix. the value
!                           of the spline approximation at t is
!                           s(t)=((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!                           where x(i).le.t.lt.x(i+1) and
!                           d = t-x(i).
!                ic     - row dimension of matrix c exactly
!                           as specified in the dimension
!                           statement in the calling program. (input)
!                var    - error variance. (input/output)
!                           if var is negative (i.e. unknown) then
!                           the smoothing parameter is determined
!                           by minimizing the generalized cross validation
!                           and an estimate of the error variance is
!                           returned in var.
!                           if var is non-negative (i.e. known) then the
!                           smoothing parameter is determined to minimize
!                           an estimate, which depends on var, of the true
!                           mean square error, and var is unchanged.
!                           in particular, if var is zero, then an
!                           interpolating natural cubic spline is calculated.
!                           var should be set to 1 if absolute standard
!                           deviations have been provided in df (see above).
!                job    - job selection parameter. (input)
!                         job = 0 should be selected if point standard error
!                           estimates are not required in se.
!                         job = 1 should be selected if point standard error
!                           estimates are required in se.
!                se     - vector of length n containing bayesian standard
!                           error estimates of the fitted spline values in y.
!                           se is not referenced if job=0. (output)
!                wk     - work vector of length 7*(n + 2). on normal exit the
!                           first 7 values of wk are assigned as follows:-
!
!                           wk(1) = smoothing parameter (= rho/(rho + 1))
!                           wk(2) = estimate of the number of degrees of
!                                   freedom of the residual sum of squares
!                           wk(3) = generalized cross validation
!                           wk(4) = mean square residual
!                           wk(5) = estimate of the true mean square error
!                                   at the data points
!                           wk(6) = estimate of the error variance
!                           wk(7) = mean square value of the df(i)
!
!                           if wk(1)=0 (rho=0) an interpolating natural cubic
!                           spline has been calculated.
!                           if wk(1)=1 (rho=infinite) a least squares
!                           regression line has been calculated.
!                           wk(2) is an estimate of the number of degrees of
!                           freedom of the residual which reduces to the
!                           usual value of n-2 when a least squares regression
!                           line is calculated.
!                           wk(3),wk(4),wk(5) are calculated with the df(i)
!                           scaled to have mean square value 1.  the
!                           unscaled values of wk(3),wk(4),wk(5) may be
!                           calculated by dividing by wk(7).
!                           wk(6) coincides with the output value of var if
!                           var is negative on input.  it is calculated with
!                           the unscaled values of the df(i) to facilitate
!                           comparisons with a priori variance estimates.
!
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 129, ic is less than n-1.
!                           ier = 130, n is less than 3.
!                           ier = 131, input abscissae are not
!                             ordered so that x(i).lt.x(i+1).
!                           ier = 132, df(i) is not positive for some i.
!                           ier = 133, job is not 0 or 1.
!
!   precision/hardware  - double
!
!   required routines   - spint1,spfit1,spcof1,sperr1
!
!   remarks      the number of arithmetic operations required by the
!                subroutine is proportional to n.  the subroutine
!                uses an algorithm developed by m.f. hutchinson and
!                f.r. de hoog, 'smoothing noisy data with spline
!                functions', numer. math. (in press)
!
!-----------------------------------------------------------------------
!
  !> FIXME : Add documentation
  subroutine cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!
!---specifications for arguments---
    integer n,ic,job,ier
!
! All real variables should be double precision:
!
    real, dimension (:) :: x, f, df, y, se
    real, dimension (:,:) :: c
    real, dimension (0:n+1,7) :: wk
    real :: var
!
!---specifications for local variables---
    real :: delta,err,gf1,gf2,gf3,gf4,r1,r2,r3,r4,avh,avdf,avar,stat(6),p,q
    real :: ratio = 2.
    real :: tau = 1.618033989
    real :: zero = 0.
    real :: one = 1.

    integer :: i
!
!---initialize---
    ier = 133
    if (job.lt.0 .or. job.gt.1) go to 140
    call spint1(x,avh,f,df,avdf,n,y,c,ic,wk,wk(0,4),ier)
    if (ier.ne.0) go to 140
    avar = var
    if (var.gt.zero) avar = var*avdf*avdf
!
!---check for zero variance---
    if (var.ne.zero) go to 10
    r1 = zero
    go to 90
!
!---find local minimum of gcv or the expected mean square error---
10  r1 = one
    r2 = ratio*r1
    call spfit1(x,avh,df,n,r2,p,q,gf2,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
20  call spfit1(x,avh,df,n,r1,p,q,gf1,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    if (gf1.gt.gf2) go to 30
!
!---exit if p zero---
    if (p.le.zero) go to 100
    r2 = r1
    gf2 = gf1
    r1 = r1/ratio
    go to 20

30  r3 = ratio*r2
40  call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    if (gf3.gt.gf2) go to 50
!
!---exit if q zero---
    if (q.le.zero) go to 100
    r2 = r3
    gf2 = gf3
    r3 = ratio*r3
    go to 40

50  r2 = r3
    gf2 = gf3
    delta = (r2-r1)/tau
    r4 = r1 + delta
    r3 = r2 - delta
    call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    call spfit1(x,avh,df,n,r4,p,q,gf4,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
!
!---golden section search for local minimum---
60  if (gf3.gt.gf4) go to 70
    r2 = r4
    gf2 = gf4
    r4 = r3
    gf4 = gf3
    delta = delta/tau
    r3 = r2 - delta
    call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    go to 80

70  r1 = r3
    gf1 = gf3
    r3 = r4
    gf3 = gf4
    delta = delta/tau
    r4 = r1 + delta
    call spfit1(x,avh,df,n,r4,p,q,gf4,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
80  err = (r2-r1)/ (r1+r2)
    if (err*err+one.gt.one .and. err.gt.1.0e-6) go to 60
    r1 = (r1+r2)*0.5
!
!---calculate spline coefficients---
90  call spfit1(x,avh,df,n,r1,p,q,gf1,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
100 call spcof1(x,avh,f,df,n,p,q,y,c,ic,wk(:,6),wk(:,7))
!
!---optionally calculate standard error estimates---
    if (var.ge.zero) go to 110
    avar = stat(6)
    var = avar/ (avdf*avdf)
110 if (job.eq.1) call sperr1(x,avh,df,n,wk,p,avar,se)
!
!---unscale df---
    df = df * avdf
!
!--put statistics in wk---
    do i = 0,5
       wk(i,1) = stat(i+1)
    end do
    wk(5,1) = stat(6)/ (avdf*avdf)
    wk(6,1) = avdf*avdf
    go to 150
!
!---check for error condition---
140 continue
!     if (ier.ne.0) continue
150 return
  end subroutine cubgcv

  !> FIXME : Add documentation
  subroutine spint1(x,avh,y,dy,avdy,n,a,c,ic,r,t,ier)
!
! initializes the arrays c, r and t for one dimensional cubic
! smoothing spline fitting by subroutine spfit1.  the values
! df(i) are scaled so that the sum of their squares is n
! and the average of the differences x(i+1) - x(i) is calculated
! in avh in order to avoid underflow and overflow problems in
! spfit1.
!
! subroutine sets ier if elements of x are non-increasing,
! if n is less than 3, if ic is less than n-1 or if dy(i) is
! not positive for some i.
!
!---specifications for arguments---
    integer n,ic,ier
! All variables should be double precision
    real, dimension (:) :: x, y, dy, a
    real, dimension (:,:) :: c
    real, dimension (0:n+1,3) :: r
    real, dimension (0:n+1,2) :: t
    real :: avh, avdy

!
!---specifications for local variables---
    integer i
    real :: e,f,g,h
    real :: zero = 0.

!
!---initialization and input checking---
    ier = 0
    if (n.lt.3) go to 60
    if (ic.lt.n-1) go to 70
!
!---get average x spacing in avh---
    g = zero
    do i = 1,n - 1
       h = x(i+1) - x(i)
       if (h.le.zero) go to 80
       g = g + h
    end do
    avh = g/ (n-1)
!
!---scale relative weights---
    g = zero
    do i = 1,n
       if (dy(i).le.zero) go to 90
       g = g + dy(i)*dy(i)
    end do
    avdy = sqrt(g/n)

    dy = dy / avdy

!
!---initialize h,f---
    h = (x(2)-x(1))/avh
    f = (y(2)-y(1))/h
!
!---calculate a,t,r---
    do i = 2,n - 1
       g = h
       h = (x(i+1)-x(i))/avh
       e = f
       f = (y(i+1)-y(i))/h
       a(i) = f - e
       t(i,1) = 2.0d0* (g+h)/3.0d0
       t(i,2) = h/3.0d0
       r(i,3) = dy(i-1)/g
       r(i,1) = dy(i+1)/h
       r(i,2) = -dy(i)/g - dy(i)/h
    end do
!
!---calculate c = r'*r---
    r(n,2) = zero
    r(n,3) = zero
    r(n+1,3) = zero
    do i = 2,n - 1
       c(i,1) = r(i,1)*r(i,1) + r(i,2)*r(i,2) + r(i,3)*r(i,3)
       c(i,2) = r(i,1)*r(i+1,2) + r(i,2)*r(i+1,3)
       c(i,3) = r(i,1)*r(i+2,3)
    end do
    return
!
!---error conditions---
60  ier = 130
    return

70  ier = 129
    return

80  ier = 131
    return

   90 ier = 132
    return
  end subroutine spint1

  !> FIXME : Add documentation
  subroutine spfit1(x,avh,dy,n,rho,p,q,fun,var,stat,a,c,ic,r,t,u,v)
!
! fits a cubic smoothing spline to data with relative
! weighting dy for a given value of the smoothing parameter
! rho using an algorithm based on that of c.h. reinsch (1967),
! numer. math. 10, 177-183.
!
! the trace of the influence matrix is calculated using an
! algorithm developed by m.f.hutchinson and f.r.de hoog (numer.
! math., in press), enabling the generalized cross validation
! and related statistics to be calculated in order n operations.
!
! the arrays a, c, r and t are assumed to have been initialized
! by the subroutine spint1.  overflow and underflow problems are
! avoided by using p=rho/(1 + rho) and q=1/(1 + rho) instead of
! rho and by scaling the differences x(i+1) - x(i) by avh.
!
! the values in df are assumed to have been scaled so that the
! sum of their squared values is n.  the value in var, when it is
! non-negative, is assumed to have been scaled to compensate for
! the scaling of the values in df.
!
! the value returned in fun is an estimate of the true mean square
! when var is non-negative, and is the generalized cross validation
! when var is negative.
!
!---specifications for arguments---
    integer ic,n
! all variables should be double precision:
    real, dimension (:) :: x, dy, a
    real, dimension (6) :: stat
    real, dimension (ic,3) :: c
    real, dimension (0:n+1,3) :: r
    real, dimension (0:n+1,2) :: t
    real, dimension (0:n+1) :: u, v
    real :: fun, var, avh, p, q, rho
!
!---local variables---
    integer i
    real :: e,f,g,h,rho1
    real :: zero = 0.
    real :: one = 1.
    real :: two = 2.

!
!---use p and q instead of rho to prevent overflow or underflow---
    rho1 = one + rho
    p = rho/rho1
    q = one/rho1
    if (rho1.eq.one) p = zero
    if (rho1.eq.rho) q = zero
!
!---rational cholesky decomposition of p*c + q*t---
    f = zero
    g = zero
    h = zero
    do i = 0,1
       r(i,1) = zero
    end do
    do i = 2,n - 1
       r(i-2,3) = g*r(i-2,1)
       r(i-1,2) = f*r(i-1,1)
       r(i,1) = one/ (p*c(i,1)+q*t(i,1)-f*r(i-1,2)-g*r(i-2,3))
       f = p*c(i,2) + q*t(i,2) - h*r(i-1,2)
       g = h
       h = p*c(i,3)
    end do
!
!---solve for u---
    u(0) = zero
    u(1) = zero
    do i = 2,n - 1
       u(i) = a(i) - r(i-1,2)*u(i-1) - r(i-2,3)*u(i-2)
    end do
    u(n) = zero
    u(n+1) = zero
    do i = n - 1,2,-1
       u(i) = r(i,1)*u(i) - r(i,2)*u(i+1) - r(i,3)*u(i+2)
    end do
!
!---calculate residual vector v---
    e = zero
    h = zero
    do i = 1,n - 1
       g = h
       h = (u(i+1)-u(i))/ ((x(i+1)-x(i))/avh)
       v(i) = dy(i)* (h-g)
       e = e + v(i)*v(i)
    end do
    v(n) = dy(n)* (-h)
    e = e + v(n)*v(n)
!
!---calculate upper three bands of inverse matrix---
    r(n,1) = zero
    r(n,2) = zero
    r(n+1,1) = zero
    do i = n - 1,2,-1
       g = r(i,2)
       h = r(i,3)
       r(i,2) = -g*r(i+1,1) - h*r(i+1,2)
       r(i,3) = -g*r(i+1,2) - h*r(i+2,1)
       r(i,1) = r(i,1) - g*r(i,2) - h*r(i,3)
    end do
!
!---calculate trace---
    f = zero
    g = zero
    h = zero
    do i = 2,n - 1
       f = f + r(i,1)*c(i,1)
       g = g + r(i,2)*c(i,2)
       h = h + r(i,3)*c(i,3)
    end do
    f = f + two* (g+h)
!
!---calculate statistics---
    stat(1) = p
    stat(2) = f*p
    stat(3) = n*e/ (f*f)
    stat(4) = e*p*p/n
    stat(6) = e*p/f
    if (var.ge.zero) go to 80
    stat(5) = stat(6) - stat(4)
    fun = stat(3)
    go to 90

80  stat(5) = amax1(stat(4)-two*var*stat(2)/n+var,zero)
    fun = stat(5)
90  return
  end subroutine spfit1

  !> FIXME : Add documentation
  subroutine sperr1(x,avh,dy,n,r,p,var,se)
!
! calculates bayesian estimates of the standard errors of the fitted
! values of a cubic smoothing spline by calculating the diagonal elements
! of the influence matrix.
!
!---specifications for arguments---
    integer n
    real, dimension (:) :: x, dy, se
    real, dimension (0:n+1, 3) :: r
    real :: avh, p, var
!
!---specifications for local variables---
    integer i
    real :: f,g,h,f1,g1,h1
    real :: zero = 0.
    real :: one = 1.
!
!---initialize---
    h = avh/ (x(2)-x(1))
    se(1) = one - p*dy(1)*dy(1)*h*h*r(2,1)
    r(1,1) = zero
    r(1,2) = zero
    r(1,3) = zero
!
!---calculate diagonal elements---
    do i = 2,n - 1
       f = h
       h = avh/ (x(i+1)-x(i))
       g = -f - h
       f1 = f*r(i-1,1) + g*r(i-1,2) + h*r(i-1,3)
       g1 = f*r(i-1,2) + g*r(i,1) + h*r(i,2)
       h1 = f*r(i-1,3) + g*r(i,2) + h*r(i+1,1)
       se(i) = one - p*dy(i)*dy(i)* (f*f1+g*g1+h*h1)
    end do
    se(n) = one - p*dy(n)*dy(n)*h*h*r(n-1,1)
!
!---calculate standard error estimates---
    do i = 1,n
       se(i) = sqrt(amax1(se(i)*var,zero))*dy(i)
    end do

  end subroutine sperr1

  !> FIXME : Add documentation
  subroutine spcof1(x,avh,y,dy,n,p,q,a,c,ic,u,v)
!
! calculates coefficients of a cubic smoothing spline from
! parameters calculated by subroutine spfit1.
!
!---specifications for arguments---
    integer ic,n
    real, dimension (:) :: x, y, dy, a
    real, dimension (ic, 3) :: c
    real, dimension (0:n+1) :: u, v
    real :: p, q, avh
!
!---specifications for local variables---
    integer i
    real :: h,qh
!
!---calculate a---
    qh = q/ (avh*avh)
    do i = 1,n
       a(i) = y(i) - p*dy(i)*v(i)
       u(i) = qh*u(i)
    end do
!
!---calculate c---
    do i = 1,n - 1
       h = x(i+1) - x(i)
       c(i,3) = (u(i+1)-u(i))/ (3.0d0*h)
       c(i,1) = (a(i+1)-a(i))/h - (h*c(i,3)+u(i))*h
       c(i,2) = u(i)
    end do
  end subroutine spcof1

!---c cubgcv test driver
!---c ------------------
!---c
!---c author          - m.f.hutchinson
!---c                   csiro division of water and land resources
!---c                   gpo box 1666
!---c                   canberra act 2601
!---c                   australia
!---c
!---c latest revision - 7 august 1986
!---c
!---c computer        - vax/double
!---c
!---c usage           - main program
!---c
!---c required routines - cubgcv,spint1,spfit1,spcof1,sperr1,ggrand
!---c
!---c remarks   uses subroutine cubgcv to fit a cubic smoothing spline
!---c           to 50 data points which are generated by adding a random
!---c           variable with uniform density in the interval [-0.3,0.3]
!---c           to 50 points sampled from the curve  y=sin(3*pi*x/2).
!---c           random deviates in the interval [0,1] are generated by the
!---c           double precision function ggrand (similar to imsl function
!---c           ggubfs).  the abscissae are unequally spaced in [0,1].
!---c
!---c           point standard error estimates are returned in se by
!---c           setting job=1.  the error variance estimate is returned
!---c           in var.  it compares favourably with the true value of 0.03.
!---c           summary statistics from the array wk are written to
!---c           unit 6.  data values and fitted values with estimated
!---c           standard errors are also written to unit 6.
!---c
!---      parameter (n=50, ic=49)
!---c
!---      integer            job,ier
!---      double precision   x(n),f(n),y(n),df(n),c(ic,3),wk(7*(n+2)),
!---     *                   var,se(n)
!---      double precision   ggrand,dseed
!---c
!---c---initialize---
!---      dseed=1.2345d4
!---      job=1
!---      var=-1.0d0
!---c
!---c---calculate data points---
!---      do 10 i=1,n
!---      x(i)=(i - 0.5)/n + (2.0*ggrand(dseed) - 1.0)/(3.0*n)
!---      f(i)=dsin(4.71238*x(i)) + (2.0*ggrand(dseed) - 1.0)*0.3
!---      df(i)=1.0d0
!---  10  continue
!---c
!---c---fit cubic spline---
!---      call cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---c
!---c---write out results---
!---      write(6,20)
!---  20  format(' cubgcv test driver results:')
!---      write(6,30) ier,var,wk(3),wk(4),wk(2)
!---  30  format(/' ier =',i4/' var =',f7.4/
!---     *        ' generalized cross validation =',f7.4/
!---     *        ' mean square residual         =',f7.4/
!---     *        ' residual degrees of freedom  =',f7.2)
!---      write(6,40)
!---  40  format(/' input data',17x,'output results'//
!---     *         '   i    x(i)    f(i)',6x,'    y(i)   se(i)',
!---     *          '      c(i,1)      c(i,2)      c(i,3)')
!---      do 60 i=1,n-1
!---      write(6,50) i,x(i),f(i),y(i),se(i),(c(i,j),j=1,3)
!---  50  format(i4,2f8.4,6x,2f8.4,3e12.4)
!---  60  continue
!---      write(6,50) n,x(n),f(n),y(n),se(n)
!---      stop
!---      end
!---      double precision function ggrand(dseed)
!---c
!---c double precision uniform random number generator
!---c
!---c constants: a = 7**5
!---c            b = 2**31 - 1
!---c            c = 2**31
!---c
!---c reference: imsl manual, chapter g - generation and testing of
!---c                                     random numbers
!---c
!---c---specifications for arguments---
!---      double precision dseed
!---c
!---c---specifications for local variables---
!---      double precision a,b,c,s
!---c
!---      data a,b,c/16807.0d0, 2147483647.0d0, 2147483648.0d0/
!---c
!---      s=dseed
!---      s=dmod(a*s, b)
!---      ggrand=s/c
!---      dseed=s
!---      return
!---      end
!---
!---cubgcv test driver results:
!---
!---ier =   0
!---var = 0.0279
!---generalized cross validation = 0.0318
!---mean square residual         = 0.0246
!---residual degrees of freedom  =  43.97
!---
!---input data                 output results
!---
!---  i    x(i)    f(i)          y(i)   se(i)      c(i,1)      c(i,2)      c(i,3)
!---  1  0.0046  0.2222        0.0342  0.1004  0.3630e+01  0.0000e+00  0.2542e+02
!---  2  0.0360 -0.1098        0.1488  0.0750  0.3705e+01  0.2391e+01 -0.9537e+01
!---  3  0.0435 -0.0658        0.1767  0.0707  0.3740e+01  0.2175e+01 -0.4233e+02
!---  4  0.0735  0.3906        0.2900  0.0594  0.3756e+01 -0.1642e+01 -0.2872e+02
!---  5  0.0955  0.6054        0.3714  0.0558  0.3642e+01 -0.3535e+01  0.2911e+01
!---  6  0.1078  0.3034        0.4155  0.0549  0.3557e+01 -0.3428e+01 -0.1225e+02
!---  7  0.1269  0.7386        0.4822  0.0544  0.3412e+01 -0.4131e+01  0.2242e+02
!---  8  0.1565  0.4616        0.5800  0.0543  0.3227e+01 -0.2143e+01  0.6415e+01
!---  9  0.1679  0.4315        0.6165  0.0543  0.3180e+01 -0.1923e+01 -0.1860e+02
!--- 10  0.1869  0.5716        0.6762  0.0544  0.3087e+01 -0.2985e+01 -0.3274e+02
!--- 11  0.2149  0.6736        0.7595  0.0542  0.2843e+01 -0.5733e+01 -0.4435e+02
!--- 12  0.2356  0.7388        0.8155  0.0539  0.2549e+01 -0.8486e+01 -0.5472e+02
!--- 13  0.2557  1.1953        0.8630  0.0537  0.2139e+01 -0.1180e+02 -0.9784e+01
!--- 14  0.2674  1.0299        0.8864  0.0536  0.1860e+01 -0.1214e+02  0.9619e+01
!--- 15  0.2902  0.7981        0.9225  0.0534  0.1322e+01 -0.1149e+02 -0.7202e+01
!--- 16  0.3155  0.8973        0.9485  0.0532  0.7269e+00 -0.1203e+02 -0.1412e+02
!--- 17  0.3364  1.2695        0.9583  0.0530  0.2040e+00 -0.1292e+02  0.2796e+02
!--- 18  0.3557  0.7253        0.9577  0.0527 -0.2638e+00 -0.1130e+02 -0.3453e+01
!--- 19  0.3756  1.2127        0.9479  0.0526 -0.7176e+00 -0.1151e+02  0.3235e+02
!--- 20  0.3881  0.7304        0.9373  0.0525 -0.9889e+00 -0.1030e+02  0.4381e+01
!--- 21  0.4126  0.9810        0.9069  0.0525 -0.1486e+01 -0.9977e+01  0.1440e+02
!--- 22  0.4266  0.7117        0.8842  0.0525 -0.1756e+01 -0.9373e+01 -0.8925e+01
!--- 23  0.4566  0.7203        0.8227  0.0524 -0.2344e+01 -0.1018e+02 -0.2278e+02
!--- 24  0.4704  0.9242        0.7884  0.0524 -0.2637e+01 -0.1112e+02 -0.4419e+01
!--- 25  0.4914  0.7345        0.7281  0.0523 -0.3110e+01 -0.1140e+02 -0.3562e+01
!--- 26  0.5084  0.7378        0.6720  0.0524 -0.3500e+01 -0.1158e+02  0.5336e+01
!--- 27  0.5277  0.7441        0.6002  0.0525 -0.3941e+01 -0.1127e+02  0.2479e+02
!--- 28  0.5450  0.5612        0.5286  0.0527 -0.4310e+01 -0.9980e+01  0.2920e+02
!--- 29  0.5641  0.5049        0.4429  0.0529 -0.4659e+01 -0.8309e+01  0.3758e+02
!--- 30  0.5857  0.4725        0.3390  0.0531 -0.4964e+01 -0.5878e+01  0.5563e+02
!--- 31  0.6159  0.1380        0.1850  0.0531 -0.5167e+01 -0.8307e+00  0.4928e+02
!--- 32  0.6317  0.1412        0.1036  0.0531 -0.5157e+01  0.1499e+01  0.5437e+02
!--- 33  0.6446 -0.1110        0.0371  0.0531 -0.5091e+01  0.3614e+01  0.3434e+02
!--- 34  0.6707 -0.2605       -0.0927  0.0532 -0.4832e+01  0.6302e+01  0.1164e+02
!--- 35  0.6853 -0.1284       -0.1619  0.0533 -0.4640e+01  0.6812e+01  0.1617e+02
!--- 36  0.7064 -0.3452       -0.2564  0.0536 -0.4332e+01  0.7834e+01  0.4164e+01
!--- 37  0.7310 -0.5527       -0.3582  0.0538 -0.3939e+01  0.8141e+01 -0.2214e+02
!--- 38  0.7531 -0.3459       -0.4415  0.0540 -0.3611e+01  0.6674e+01 -0.9205e+01
!--- 39  0.7686 -0.5902       -0.4961  0.0541 -0.3410e+01  0.6245e+01 -0.2193e+02
!--- 40  0.7952 -0.7644       -0.5828  0.0541 -0.3125e+01  0.4494e+01 -0.4649e+02
!--- 41  0.8087 -0.5392       -0.6242  0.0541 -0.3029e+01  0.2614e+01 -0.3499e+02
!--- 42  0.8352 -0.4247       -0.7031  0.0539 -0.2964e+01 -0.1603e+00  0.2646e+01
!--- 43  0.8501 -0.6327       -0.7476  0.0538 -0.2967e+01 -0.4132e-01  0.1817e+02
!--- 44  0.8726 -0.9983       -0.8139  0.0538 -0.2942e+01  0.1180e+01 -0.6774e+01
!--- 45  0.8874 -0.9082       -0.8574  0.0542 -0.2911e+01  0.8778e+00 -0.1364e+02
!--- 46  0.9139 -0.8930       -0.9340  0.0566 -0.2893e+01 -0.2044e+00 -0.8094e+01
!--- 47  0.9271 -1.0233       -0.9723  0.0593 -0.2903e+01 -0.5258e+00 -0.1498e+02
!--- 48  0.9473 -0.8839       -1.0313  0.0665 -0.2942e+01 -0.1433e+01  0.4945e+01
!--- 49  0.9652 -1.0172       -1.0843  0.0766 -0.2989e+01 -0.1168e+01  0.1401e+02
!--- 50  0.9930 -1.2715       -1.1679  0.0998
!---documentation:
!---c   computer            - vax/double
!---c
!---c   author              - m.f.hutchinson
!---c                         csiro division of mathematics and statistics
!---c                         p.o. box 1965
!---c                         canberra, act 2601
!---c                         australia
!---c
!---c   latest revision     - 15 august 1985
!---c
!---c   purpose             - cubic spline data smoother
!---c
!---c   usage               - call cubgcv (x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---c
!---c   arguments    x      - vector of length n containing the
!---c                           abscissae of the n data points
!---c                           (x(i),f(i)) i=1..n. (input) x
!---c                           must be ordered so that
!---c                           x(i) .lt. x(i+1).
!---c                f      - vector of length n containing the
!---c                           ordinates (or function values)
!---c                           of the n data points (input).
!---c                df     - vector of length n. (input/output)
!---c                           df(i) is the relative standard deviation
!---c                           of the error associated with data point i.
!---c                           each df(i) must be positive.  the values in
!---c                           df are scaled by the subroutine so that
!---c                           their mean square value is 1, and unscaled
!---c                           again on normal exit.
!---c                           the mean square value of the df(i) is returned
!---c                           in wk(7) on normal exit.
!---c                           if the absolute standard deviations are known,
!---c                           these should be provided in df and the error
!---c                           variance parameter var (see below) should then
!---c                           be set to 1.
!---c                           if the relative standard deviations are unknown,
!---c                           set each df(i)=1.
!---c                n      - number of data points (input).
!---c                           n must be .ge. 3.
!---c                y,c    - spline coefficients. (output) y
!---c                           is a vector of length n. c is
!---c                           an n-1 by 3 matrix. the value
!---c                           of the spline approximation at t is
!---c                           s(t)=((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!---c                           where x(i).le.t.lt.x(i+1) and
!---c                           d = t-x(i).
!---c                ic     - row dimension of matrix c exactly
!---c                           as specified in the dimension
!---c                           statement in the calling program. (input)
!---c                var    - error variance. (input/output)
!---c                           if var is negative (i.e. unknown) then
!---c                           the smoothing parameter is determined
!---c                           by minimizing the generalized cross validation
!---c                           and an estimate of the error variance is
!---c                           returned in var.
!---c                           if var is non-negative (i.e. known) then the
!---c                           smoothing parameter is determined to minimize
!---c                           an estimate, which depends on var, of the true
!---c                           mean square error, and var is unchanged.
!---c                           in particular, if var is zero, then an
!---c                           interpolating natural cubic spline is calculated.
!---c                           var should be set to 1 if absolute standard
!---c                           deviations have been provided in df (see above).
!---c                job    - job selection parameter. (input)
!---c                         job = 0 should be selected if point standard error
!---c                           estimates are not required in se.
!---c                         job = 1 should be selected if point standard error
!---c                           estimates are required in se.
!---c                se     - vector of length n containing bayesian standard
!---c                           error estimates of the fitted spline values in y.
!---c                           se is not referenced if job=0. (output)
!---c                wk     - work vector of length 7*(n + 2). on normal exit the
!---c                           first 7 values of wk are assigned as follows:-
!---c
!---c                           wk(1) = smoothing parameter (= rho/(rho + 1))
!---c                           wk(2) = estimate of the number of degrees of
!---c                                   freedom of the residual sum of squares
!---c                           wk(3) = generalized cross validation
!---c                           wk(4) = mean square residual
!---c                           wk(5) = estimate of the true mean square error
!---c                                   at the data points
!---c                           wk(6) = estimate of the error variance
!---c                           wk(7) = mean square value of the df(i)
!---c
!---c                           if wk(1)=0 (rho=0) an interpolating natural cubic
!---c                           spline has been calculated.
!---c                           if wk(1)=1 (rho=infinite) a least squares
!---c                           regression line has been calculated.
!---c                           wk(2) is an estimate of the number of degrees of
!---c                           freedom of the residual which reduces to the
!---c                           usual value of n-2 when a least squares regression
!---c                           line is calculated.
!---c                           wk(3),wk(4),wk(5) are calculated with the df(i)
!---c                           scaled to have mean square value 1.  the
!---c                           unscaled values of wk(3),wk(4),wk(5) may be
!---c                           calculated by dividing by wk(7).
!---c                           wk(6) coincides with the output value of var if
!---c                           var is negative on input.  it is calculated with
!---c                           the unscaled values of the df(i) to facilitate
!---c                           comparisons with a priori variance estimates.
!---c
!---c                ier    - error parameter. (output)
!---c                         terminal error
!---c                           ier = 129, ic is less than n-1.
!---c                           ier = 130, n is less than 3.
!---c                           ier = 131, input abscissae are not
!---c                             ordered so that x(i).lt.x(i+1).
!---c                           ier = 132, df(i) is not positive for some i.
!---c                           ier = 133, job is not 0 or 1.
!---c
!---c   precision/hardware  - double
!---c
!---c   required routines   - spint1,spfit1,spcof1,sperr1
!---c
!---c   remarks      the number of arithmetic operations required by the
!---c                subroutine is proportional to n.  the subroutine
!---c                uses an algorithm developed by m.f. hutchinson and
!---c                f.r. de hoog, 'smoothing noisy data with spline
!---c                functions', numer. math. (in press)
!---c
!---c-----------------------------------------------------------------------
!---c
!---
!---
!---algorithm
!---
!---cubgcv calculates a natural cubic spline curve which smoothes a given set
!---of data points, using statistical considerations to determine the amount
!---of smoothing required, as described in reference 2.  if the error variance
!---is known, it should be supplied to the routine in var.  the degree of
!---smoothing is then determined by minimizing an unbiased estimate of the true
!---mean square error.  on the other hand, if the error variance is not known, var
!---should be set to -1.0.  the routine then determines the degree of smoothing
!---by minimizing the generalized cross validation.  this is asymptotically the
!---same as minimizing the true mean square error (see reference 1).  in this
!---case, an estimate of the error variance is returned in var which may be
!---compared with any a priori approximate estimates.  in either case, an
!---estimate of the true mean square error is returned in wk(5).  this estimate,
!---however, depends on the error variance estimate, and should only be accepted
!---if the error variance estimate is reckoned to be correct.
!---
!---if job is set to 1, bayesian estimates of the standard error of each
!---smoothed data value are returned in the array se.  these also depend on
!---the error variance estimate and should only be accepted if the error
!---variance estimate is reckoned to be correct.  see reference 4.
!---
!---the number of arithmetic operations and the amount of storage required by
!---the routine are both proportional to n, so that very large data sets may be
!---analysed.  the data points do not have to be equally spaced or uniformly
!---weighted.  the residual and the spline coefficients are calculated in the
!---manner described in reference 3, while the trace and various statistics,
!---including the generalized cross validation, are calculated in the manner
!---described in reference 2.
!---
!---when var is known, any value of n greater than 2 is acceptable.  it is
!---advisable, however, for n to be greater than about 20 when var is unknown.
!---if the degree of smoothing done by cubgcv when var is unknown is not
!---satisfactory, the user should try specifying the degree of smoothing by
!---setting var to a reasonable value.
!---
!---references:
!---
!---1.  craven, peter and wahba, grace, "smoothing noisy data with spline
!---    functions", numer. math. 31, 377-403 (1979).
!---2.  hutchinson, m.f. and de hoog, f.r., "smoothing noisy data with spline
!---    functions", numer. math. (in press).
!---3.  reinsch, c.h., "smoothing by spline functions", numer. math. 10,
!---    177-183 (1967).
!---4.  wahba, grace, "bayesian 'confidence intervals' for the cross-validated
!---    smoothing spline", j.r.statist. soc. b 45, 133-150 (1983).
!---
!---
!---example
!---
!---a sequence of 50 data points are generated by adding a random variable with
!---uniform density in the interval [-0.3,0.3] to the curve y=sin(3*pi*x/2).
!---the abscissae are unequally spaced in [0,1].  point standard error estimates
!---are returned in se by setting job to 1.  the error variance estimate is
!---returned in var.  it compares favourably with the true value of 0.03.
!---the imsl function ggubfs is used to generate sample values of a uniform
!---variable on [0,1].
!---
!---
!---input:
!---
!---      integer          n,ic,job,ier
!---      double precision x(50),f(50),y(50),df(50),c(49,3),wk(400),
!---     *                 var,se(50)
!---      double precision ggubfs,dseed,dn
!---      data dseed/1.2345d4/
!---c
!---      n=50
!---      ic=49
!---      job=1
!---      var=-1.0d0
!---      dn=n
!---c
!---      do 10 i=1,n
!---      x(i)=(i - 0.5)/dn + (2.0*ggubfs(dseed) - 1.0)/(3.0*dn)
!---      f(i)=dsin(4.71238*x(i)) + (2.0*ggubfs(dseed) - 1.0)*0.3
!---      df(i)=1.0d0
!---  10  continue
!---      call cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---       .
!---       .
!---       .
!---      end
!---
!---output:
!---
!---ier = 0
!---var = 0.0279
!---generalized cross validation = 0.0318
!---mean square residual = 0.0246
!---residual degrees of freedom = 43.97
!---for checking purposes the following output is given:
!---
!---x(1)  = 0.0046    f(1)  =  0.2222     y(1)  =  0.0342     se(1)  = 0.1004
!---x(21) = 0.4126    f(21) =  0.9810     y(21) =  0.9069     se(21) = 0.0525
!---x(41) = 0.8087    f(41) = -0.5392     y(41) = -0.6242     se(41) = 0.0541
!---
end program rungridgen
