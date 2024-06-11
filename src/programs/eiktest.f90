!> FIXME : Add documentation
program eiktest
  use mp, only: init_mp, finish_mp, proc0
  use constants, only: pi, run_name_size
  !> Inputs
  use geometry, only: surf, eqfile, equal_arc, use_large_aspect
  use geometry, only: bishop, irho, force_sym, efit_eq, dfit_eq, writelots
  use geometry, only: gen_eq, ppl_eq, local_eq, idfit_eq, chs_eq
  use geometry, only: p_prime_input, invLp_input, beta_prime_input, big, gs2d_eq
  use geometry, only: transp_eq, eqnormfile, alpha_input, dp_mult
  !> Methods from geometry
  use geometry, only: eikcoefs_output_type, eikcoefs, geom
  use geometry, only: pbarfun, pbarofrho, rpofrho, diameter
  use file_utils, only: init_file_utils, run_name
  use gs2_io_grids, only: nc_write_grid_sizes, nc_write_grid_scalars, nc_write_grids
  use gs2_io_grids, only: nc_grid_file_open, nc_grid_file_close
  implicit none

  integer :: i, nth, ntgrid, ntheta, shotnum, geoType
  real :: q, qq, qplus, qmnus
  real :: bb, profile_fac
  real :: rk, rkm, rkp, rkpri
  real :: rt, rtm, rtp, rtpri
  real :: rd, rdm, rdp, rdpri
  real :: rhocm, rhocp
  real :: bbar, rmin, rmax, pbar
  real :: akappa, akappri, tri, tripri
  real :: rhoc, rmaj, r_geo, qinp, s_hat_input, shift

  real :: diffscheme, beta_p1, beta_p2, beta_prime_times, beta_prime_over
  real :: tstar
  integer :: nbeta, funit
  logical :: no_end_point
  logical :: fast, dipole, list, Xanthopoulos, in_nt
  character(len=:), allocatable :: in_stem, out_stem

  character(len=run_name_size) :: ncout
  type(eikcoefs_output_type) :: eikcoefs_results
  real, dimension(:), allocatable :: theta, bmag, aplot, aprime, bpol, cvdrift, grho
  real, dimension(:), allocatable :: cvdrift0, gbdrift, gbdrift0, gds2, gds21, gds22
  real, dimension(:), allocatable :: gradpar, zplot, zprime, rprime, rplot
  real :: qsf, surfarea, kxfac, dvdrhon, shat, drhodpsi, dbetadrho, delrho
  integer :: fort11_unit, isym, iflux, itor, nperiod
  namelist/stuff/ntheta,nperiod,geoType,rmaj,akappri,akappa,shift,equal_arc, &
       rhoc,rmin,rmax,itor,qinp,iflux,delrho,tri,bishop, &
       irho,isym,tripri,efit_eq,dfit_eq,writelots,R_geo, &
       gen_eq, ppl_eq, eqfile, local_eq, idfit_eq,&
       chs_eq, force_sym, s_hat_input,p_prime_input,invLp_input,beta_prime_input, &
       diffscheme,nbeta,beta_p1,beta_p2,alpha_input,big, &
       beta_prime_times, beta_prime_over, fast, profile_fac, &
       tstar, shotnum, gs2d_eq, transp_eq, Xanthopoulos, dipole, &
       eqnormfile, no_end_point, in_nt

  call init_mp

  ! set defaults
  gen_eq = .false.
  ppl_eq = .false.
  chs_eq = .false.
  idfit_eq = .false.
  dfit_eq = .false.
  efit_eq = .false.
  gs2d_eq = .false.
  transp_eq = .false.
  Xanthopoulos = .false.
  dipole = .false.
  shotnum=1001109020
  tstar = 1.0
  ntheta=32; nperiod=1
  equal_arc = .false. ; bishop = 0 ; fast = .true.
  big = 8
  eqfile='dskeq.cdf'
  rhoc=0.5
  rmin=-1.0 ; rmax=-1.0
  beta_prime_times = -1.
  profile_fac = 0.5
  geoType = 0
  rmaj = 3   ;     akappa = 1   ;     tri = 0   ;     qinp = 2
  R_geo = 3  ;     akappri = 0  ;     tripri=0  ;     shat = 1
  shift = 0
  s_hat_input = shat
  beta_prime_input = -1
  p_prime_input = -2
  invLp_input = 5
  irho = 2 ;      isym = -1 ; itor = -1 ; iflux=-1
  delrho = 0.01
  no_end_point = .false. ; in_nt = .false.

  if (command_argument_count() == 0) then
     allocate(character(len=3)::in_stem)
     in_stem = 'eik'
     allocate(character(len=0)::out_stem)
     out_stem = ''
  else
     call init_file_utils(list)
     allocate(character(len=len_trim(run_name))::in_stem)
     allocate(character(len=len_trim(run_name)+1)::out_stem)
     in_stem = trim(run_name)
     out_stem = trim(in_stem)//'.'
  endif

  !     read in data
  open(newunit=funit,file=in_stem//'.in')
  read(funit,nml=stuff)
  close(funit)

  if (local_eq .and. (geoType /= 0)) print*,'Warning : geoType /= 0 not supported'
  if (Xanthopoulos) write(*,*) "Option deprecated : Xanthopoulos set but unused"
  if (in_nt) write(*,*) "Option deprecated : in_nt set but unused"
  if (.not. fast) write(*,*) "Option deprecated : fast set but doesn't do anything"
  if (rmin /= -1) write(*,*) "Option deprecated : rmin set but not used"
  if (rmax /= -1) write(*,*) "Option deprecated : rmax set but not used"
  if (isym /= -1) write(*,*) "Option deprecated : isym set but not used (see force_sym)"
  if (iflux /= -1) write(*,*) "Option deprecated : iflux set but not used"
  if (itor /= -1) write(*,*) "Option deprecated : itor set but not used (see use_large_aspect)"

  surf%geoType = geoType
  surf%mMode = 2 ; surf%nMode = 3
  surf%rmaj = rmaj ; surf%R_geo = R_geo ; surf%r = rhoc ; surf%dr = delrho
  surf%sHorz = shift ; surf%delm = akappa
  surf%deln = tri ; surf%delmp = akappri
  surf%delnp = tripri ; surf%thm = 0
  surf%thn = 0 ; surf%q = qinp ; surf%shat = s_hat_input

  if(bishop >= 7) then
     if(beta_prime_times == -1) then
        write(*,*) 'beta_prime_times should be set for bishop = 7 or 8'
        beta_prime_times = 1.
     endif
     dp_mult = beta_prime_times
  endif
  fort11_unit = 12

  ! Note that if .not. local_eq then always choose use_large_aspect = F
  if ((.not. local_eq) .and. use_large_aspect) then
     use_large_aspect = .false.
     write(*,*) 'Forcing use_large_aspect to false'
  end if

  ! Catch if no q profile
  if (local_eq .and. irho /= 2) then
     irho = 2
     write(*,*) 'Forcing irho=2'
  end if

  if(local_eq) then
     if (gen_eq) write(*,*) 'Forcing gen_eq to be false'
     if (efit_eq) write(*,*) 'Forcing efit_eq to be false'
     if (dfit_eq) write(*,*) 'Forcing dfit_eq to be false'
     if (gs2d_eq) write(*,*) 'Forcing gs2d_eq to be false'
     if (idfit_eq) write(*,*) 'Forcing idfit_eq to be false'
     if (ppl_eq) write(*,*) 'Forcing ppl_eq to be false'
     if (chs_eq) write(*,*) 'Forcing chs_eq to be false'
     if (transp_eq) write(*,*) 'Forcing transp_eq to be false'
     if (any([gen_eq, ppl_eq, transp_eq, chs_eq, efit_eq, dfit_eq, gs2d_eq, idfit_eq])) &
          write(*,*) 'because local_eq'
     gen_eq = .false. ; ppl_eq = .false. ; transp_eq = .false. ; chs_eq = .false.
     efit_eq = .false. ; dfit_eq = .false. ; gs2d_eq = .false. ; idfit_eq = .false.
  endif

  open(newunit=funit,file=out_stem//'eik.out',status='unknown')

  call eikcoefs(ntheta, nperiod, eikcoefs_results)

  ! Unpack results for ease of use
  allocate(theta, source = eikcoefs_results%theta)
  allocate(bmag, source = eikcoefs_results%bmag)
  allocate(gradpar, source = eikcoefs_results%gradpar)
  allocate(grho, source = eikcoefs_results%grho)
  allocate(cvdrift, source = eikcoefs_results%cvdrift)
  allocate(gbdrift, source = eikcoefs_results%gbdrift)
  allocate(cvdrift0, source = eikcoefs_results%cvdrift0)
  allocate(gbdrift0, source = eikcoefs_results%gbdrift0)
  allocate(gds2, source = eikcoefs_results%gds2)
  allocate(gds21, source = eikcoefs_results%gds21)
  allocate(gds22, source = eikcoefs_results%gds22)
  allocate(rplot, source = eikcoefs_results%rplot)
  allocate(zplot, source = eikcoefs_results%zplot)
  allocate(aplot, source = eikcoefs_results%aplot)
  allocate(rprime, source = eikcoefs_results%rprime)
  allocate(zprime, source = eikcoefs_results%zprime)
  allocate(aprime, source = eikcoefs_results%aprime)
  allocate(bpol, source = eikcoefs_results%bpol)
  allocate(kxfac, source = eikcoefs_results%kxfac)
  allocate(qsf, source = eikcoefs_results%qsf)
  allocate(surfarea, source = eikcoefs_results%surfarea)
  allocate(dvdrhon, source = eikcoefs_results%dvdrhon)
  allocate(shat, source = eikcoefs_results%shat)
  allocate(drhodpsi, source = eikcoefs_results%drhodpsi)

  ! note that ntheta may have changed
  ntgrid = (size(theta)-1)/2 ; nth = ntheta / 2

  if (writelots) then
     bbar = 0.
     do i=-ntgrid,ntgrid-1
        bbar = bbar + bmag(i)*(theta(i+1)-theta(i))
     end do
     bbar = bbar/(2.*pi)

     write(*,*) 'B_bar = ',bbar
  end if

  if (dipole) then
     shat = max(1.e-7, shat)
  end if

  write(funit,*) ' ntgrid  nperiod  ntheta, drhodpsi, rmaj, shat, kxfac, q'
  !CMR 14/01/05 change format to prevent line breaks with intel compiler and allow runngridgen to read the output files
  write(funit,fmt='(3i8,1p5e12.4)') ntgrid, nperiod, ntheta, drhodpsi, rmaj, shat, kxfac, qsf
  !CMR

  write(funit,*) '    gbdrift     gradpar       grho    tgrid'
  do i= -ntgrid,ntgrid
     write(funit,1000) gbdrift(i), gradpar(i),grho(i),theta(i)
  enddo

  write(funit,*) '    cvdrift            gds2            bmag   tgrid'
  do i= -ntgrid,ntgrid
     write(funit,1000) cvdrift(i), gds2(i),bmag(i),theta(i)
  enddo

  write(funit,*) '    gds21             gds22  tgrid'
  do i= -ntgrid,ntgrid
     write(funit,1000) gds21(i), gds22(i),theta(i)
  enddo

  write(funit,*) '    cvdrift0          gbdrift0   tgrid'
  do i= -ntgrid,ntgrid
     write(funit,1000) cvdrift0(i),gbdrift0(i),theta(i)
  enddo

  call geom%hahm_burrell(profile_fac)
  dipole = dfit_eq .or. idfit_eq

  if (.not. dipole) then
     pbar = pbarfun(geom%rfun(rpofrho(rhoc),0.),0.)
     bb = geom%betafun(pbar)
     write(funit,*) '-8*pi/B_0**2 * dp/drho= '
     if(bb.ne.0.) then
        write(funit,1000) -dbetadrho
        write(funit,1000) -dbetadrho/bb, -0.5*dbetadrho/bb
        write(funit,1000) bb, bb/2.
     endif
  end if
  close(funit)

  open(newunit=funit,file=out_stem//'eik2.out',status='unknown')
  write(funit,*) 'dV/drhon= ',dvdrhon
  pbar = pbarofrho(rhoc)
  if (.not. dipole) then
     write (funit,*)
     write (funit,*) 'The following q and shat values'
     write (funit,*) 'were calculated directly'
     write (funit,*) 'from the equilibrium.'
     write (funit,*) 'They should match the ones above' !EGH
     q = geom%qfun(0.)
     write(funit,*) 'q_0= ',q
     q = geom%qfun(pbar)
     write(funit,*) 'rho= ',rhoc,' q(rho)= ',q
     q = geom%qfun(0.95)
     write(funit,*) 'q_95= ',q

     call geofax(rhoc, theta, rt, rk, rd, nth, ntgrid)
     if (.not. local_eq) then
        rhocp=rhoc+delrho
        rhocm=rhoc-delrho
        call geofax(rhocp, theta, rtp, rkp, rdp, nth, ntgrid)
        call geofax(rhocm, theta, rtm, rkm, rdm, nth, ntgrid)
        qplus = geom%qfun(pbarofrho(rhocp))
        qmnus = geom%qfun(pbarofrho(rhocm))
        rkpri=0.5*(rkp-rkm)/delrho
        rtpri=0.5*(rtp-rtm)/delrho
        rdpri=0.5*(rdp-rdm)/delrho
        if (dipole) then
           shat = 0.
        else
           shat = 0.5 * rhoc * (qplus - qmnus) / geom%qfun(pbar) / delrho
        endif
     else
        rkpri=akappri
        rtpri=tripri
        rdpri=shift
     endif

     qq = geom%qfun(pbar)
     write(funit,*) 'shat= ',shat
     write(funit,*) 'qinp= ',qq
     write(funit,*) 'shift= ',rdpri
     write(funit,*) 'akappa= ',rk
     write(funit,*) 'akappri= ',rkpri
     write(funit,*) 'tri= ',rt
     write(funit,*) 'tripri= ',rtpri
     write(funit,*) 'surface area= ',surfarea,' avgrmid**2'
     write(funit, *) 'diameter/a = ', diameter(rpofrho(rhoc))
     write(funit, *) 'diameter(rhoc=1.0)/a/2.0 = ', diameter(rpofrho(1.0))/2.0, ' ...  should be 1'
     write(funit, *) 'r/a = ', diameter(rpofrho(rhoc))/diameter(rpofrho(1.0))
  end if

  close(funit)

  open (newunit=funit,file='eik5.out',status='unknown')
  write (funit,'(''#theta1 gradpar2 bmag3 grho4 '', &
       &   ''gbdrift5 gbdrift06 cvdrift7 cvdrift08 '', &
       &    ''gds29 gds2110 gds2211'')')
  do i = -ntgrid, ntgrid
     write (funit,1000) theta(i),gradpar(i),bmag(i),grho(i), &
          gbdrift(i),gbdrift0(i),cvdrift(i),cvdrift0(i), &
          gds2(i),gds21(i),gds22(i)
  end do
  close (funit)

  open (newunit=funit,file='eik6.out',status='unknown')
  write (funit,fmt="('# shape: torus')")
  write (funit,fmt="('# q = ',e10.4,' drhodpsi = ',e10.4)") qsf, drhodpsi
  write (funit,fmt="('# theta1             R2                  Z3               alpha4      ', &
            &   '       Rprime5              Zprime6           alpha_prime7         bpol8')")
  do i = -ntgrid, ntgrid
     write (funit,1000) theta(i),Rplot(i),Zplot(i),aplot(i), &
          Rprime(i),Zprime(i),aprime(i),Bpol(i)
  end do
  close (funit)
1000 format(20(1x,1pg18.11))

  if (proc0) then
     ncout = trim(out_stem)//"out.nc"
     call nc_grid_file_open(ncout, "w")
     call nc_write_grid_sizes(ntheta, ntgrid, nperiod)
     call nc_write_grid_scalars(shat, drhodpsi, kxfac, qsf, rmaj, geom%B_T, geom%aminor, surfarea/dvdrhon, surfarea)
     call nc_write_grids(ntgrid, bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
          gds2, gds21, gds22, grho, theta, &
          Rplot=Rplot, Rprime=Rprime, Zplot=Zplot, Zprime=Zprime, aplot=aplot, aprime=aprime, &
          Bpol=Bpol, no_end_point_in=no_end_point)
     call nc_grid_file_close()
  end if

  call finish_mp

contains

  !> FIXME : Add documentation
  subroutine geofax(rho, theta, t, e, d, nth, ntgrid)
    use geometry, only: rpofrho, geom
    implicit none
    integer, intent(in) :: nth, ntgrid
    real, intent (in) :: rho
    real, dimension(-ntgrid:ntgrid), intent(in) :: theta
    real, intent (out) :: t, e, d
    real, dimension(-nth:nth) :: rgrid
    real :: a, b, aa, bb, rp, hmaxu, hmaxl, rmaxu, rmaxl, rmid, Rlo, Rhi, Rinnr, Routr
    integer :: i
    logical, save :: first=.true.
    if (first) then
       write(*,'("WARNING: GEOFAX USES UNINITIALISED VALUES -- RESULT MAY NOT BE SENSIBLE")')
       first = .false.
    end if
    rp = rpofrho(rho)
    do i = -nth, nth
       rgrid(i) = geom%rfun(rp, theta(i))
    end do

    ! do not treat up-down asymmetry completely for now
    hmaxu = -1.e6 ; hmaxl = 1.e6
    rmaxu = 0. ; rmaxl = 0.

    do i = 0, nth
       a = geom%Rpos(rgrid(i), theta(i)) ; b = geom%Zpos(rgrid(i), theta(i))

       aa = geom%Rpos(rgrid(-i), theta(-i)) ; bb = geom%Zpos(rgrid(-i), theta(-i))

       if (hmaxu < b ) then
          hmaxu = b ; rmaxu = a
       end if

       if (hmaxl > bb ) then
          hmaxl = bb ; rmaxl = aa
       end if

       Rlo = min(a, aa) ; Rhi = max(a, aa)
       !<DD>WARNING : Rinnr has not been defined on first loop through
       if (Rinnr > Rlo) Rinnr = Rlo
       !<DD>WARNING : Routr has not been defined on first loop through
       if (Routr < Rhi) Routr = Rhi
    end do
    a = (geom%Rpos(rgrid(0), theta(0)) - geom%Rpos(rgrid(nth), theta(nth)))
    e = (hmaxu - hmaxl) / a
    ! Ro=(Routr+Rinnr)/2. should = rcenter(rp)
    rmid = (Routr - Rinnr) / 2. ! midplane minor radius of this surface
    ! Using only the upper triangularity:
    ! wasn't properly normalized before 29 Feb 08, and don't use asin:
    !    t = asin(rcenter(rp)-rmax)
    t = (geom%rcenter(rp) - rmaxu) / rmid
    d = geom%rcenter(rp) - geom%rcenter(geom%psi_a)
  end subroutine geofax
end program eiktest
