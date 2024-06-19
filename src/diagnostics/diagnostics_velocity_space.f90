!> A module which writes out quantities which writes out 
!! diagnostic quantities which assess whether the velocity 
!! space resolution is sufficient.
module diagnostics_velocity_space
  implicit none
  private

  public :: init_diagnostics_velocity_space, write_velocity_space_checks, get_gtran
  public :: write_collision_error, finish_diagnostics_velocity_space, get_verr
  public :: collision_error

contains

  !> FIXME : Add documentation  
  subroutine init_diagnostics_velocity_space(gnostics)
    use le_grids, only: init_weights, wdim
    use mp, only: proc0, broadcast
    use diagnostics_config, only: diagnostics_type
    use collisions, only: collision_model_switch
    use collisions, only: init_lorentz_error
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics 
    
    if (.not. gnostics%write_verr) return
    
    ! initialize weights for less accurate integrals used
    ! to provide an error estimate for v-space integrals (energy and untrapped)
    if (proc0) call init_weights
    call broadcast(wdim)
    
    if (gnostics%write_cerr) then
       if (collision_model_switch == 1 .or. collision_model_switch == 5) then
          call init_lorentz_error
       else
          gnostics%write_cerr = .false.
       end if
    end if
  end subroutine init_diagnostics_velocity_space

  !> FIXME : Add documentation  
  subroutine finish_diagnostics_velocity_space()
    use le_grids, only: finish_weights
    use mp, only: proc0
    implicit none
    if (proc0) call finish_weights
  end subroutine finish_diagnostics_velocity_space

  !> FIXME : Add documentation  
  subroutine write_velocity_space_checks(gnostics)
    use mp, only: proc0
    use le_grids, only: grid_has_trapped_particles
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_time
    use collisions, only: vnmult
    use species, only: spec
    use diagnostics_config, only: diagnostics_type
    use neasyf, only: neasyf_write
    use gs2_io, only: starts, time_dim, dim_2, dim_3, dim_5
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, dimension (5,2) :: errest
    integer, dimension (5,3) :: erridx
    real :: geavg, glavg, gtavg

    errest = 0.0; erridx = 0
    geavg = 0.0 ; glavg = 0.0 ; gtavg = 0.0

    ! error estimate obtained by comparing standard integral with less-accurate integral
    call get_verr (errest, erridx, phinew, bparnew)

    ! error estimate based on monitoring amplitudes of legendre polynomial coefficients
    call get_gtran (geavg, glavg, gtavg, phinew, bparnew)

    if (.not. gnostics%vary_vnew_only .and. gnostics%writing) then
       call neasyf_write(gnostics%file_id, "vspace_lpcfrac", [geavg, glavg, gtavg], &
            dim_names=[dim_3, time_dim], start=starts(2, gnostics%nout), &
            long_name="Fraction of free energy contained in the high order coefficients of &
            & the Legendre polynomial transform of (1) energy space, (2) untrapped &
            & pitch angles and (3) trapped pitch angles (each should ideally be < 0.1).  &
            & Note that there are no trapped pitch angles for certain geometries")
       call neasyf_write(gnostics%file_id, "vspace_err", errest, &
            dim_names=[dim_5, dim_2, time_dim], start=starts(3, gnostics%nout), &
            long_name="Estimate of the (1) absolute and (2) relative errors resulting from &
            & velocity space integrals in the calculation of the following quantities &
            & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
            & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
            & angles. Relative errors should be < 0.1. ")
       call neasyf_write(gnostics%file_id, "vspace_vnewk", &
            [vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk], &
            dim_names=[dim_2, time_dim], start=starts(2, gnostics%nout), &
            long_name="If the simulation is set to vary the collisionality in order to keep &
            & error in velocity integrals to acceptable levels, contains species 1 &
            & collisionality in (1) pitch angle and (2) energy  ")

       if (gnostics%write_max_verr) then
          call neasyf_write(gnostics%file_id, "vspace_err_maxindex", erridx, &
               dim_names=[dim_5, dim_3, time_dim], start=starts(3, gnostics%nout), &
               long_name="Gives the (1) theta index, (2) ky index and (3) kx index of the maximum &
               & error resulting from the &
               & velocity space integrals in the calculation of the following quantities &
               & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
               & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
               & angles. Relative errors should be < 0.1. ")
       end if
    end if

    if (proc0 .and. gnostics%ascii_files%write_to_vres .and. (.not. gnostics%create)) call write_ascii
  contains
    !> FIXME : Add documentation
    subroutine write_ascii
      if (grid_has_trapped_particles()) then
         write(gnostics%ascii_files%lpc,"(4(1x,e13.6))") user_time, geavg, glavg, gtavg
      else
         write(gnostics%ascii_files%lpc,"(3(1x,e13.6))") user_time, geavg, glavg
      end if
      write(gnostics%ascii_files%vres,"(8(1x,e13.6))") user_time, errest(1,2), errest(2,2), errest(3,2), &
           errest(4,2), errest(5,2), vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk
      if (gnostics%write_max_verr) then
         write(gnostics%ascii_files%vres2,"(3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6))") &
              erridx(1,1), erridx(1,2), erridx(1,3), errest(1,1), &
              erridx(2,1), erridx(2,2), erridx(2,3), errest(2,1), &
              erridx(3,1), erridx(3,2), erridx(3,3), errest(3,1), &
              erridx(4,1), erridx(4,2), erridx(4,3), errest(4,1), &
              erridx(5,1), erridx(5,2), erridx(5,3), errest(5,1)
      end if
    end subroutine write_ascii
  end subroutine write_velocity_space_checks

  !> Calculate and write the collision error
  subroutine write_collision_error (gnostics)
    use diagnostics_config, only: diagnostics_type
    use fields_arrays, only: phinew, bparnew
    type(diagnostics_type), intent(in) :: gnostics
    if (gnostics%write_ascii) call collision_error(gnostics%ascii_files%cres, phinew, bparnew)
  end subroutine write_collision_error

  !> Calculate the collision error and write to text file `<runname>.cres`
  subroutine collision_error(unit, phi, bpar)
    use mp, only: proc0, send, receive, barrier
    use le_grids, only: ng2, jend, nlambda, lambda_map
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_layouts, only: lz_lo, ig_idx, idx_local, proc_id
    use gs2_layouts, only: ik_idx, ie_idx, is_idx, it_idx, il_idx, g_lo
    use collisions, only: dtot, fdf, fdb
    use redistribute, only: gather, scatter
    use gs2_time, only: user_time
    use array_utils, only: zero_array
    implicit none
    integer, intent(in) :: unit
    complex, dimension(-ntgrid:, :, :), intent(in) :: phi, bpar
    integer :: je, te, ig, il, ip, ilz, ie, is, ik, it
    integer :: igmax, ikmax, itmax, iemax, ilmax, ismax
    complex, dimension (:), allocatable :: ltmp, ftmp
    complex, dimension (:,:), allocatable :: lcoll, fdcoll, glze
    real :: etmp, emax, etot, eavg, edenom, ltmax

    allocate (ltmp(2*nlambda), ftmp(2*nlambda))
    allocate (lcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (fdcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glze(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))

    call zero_array(lcoll); call zero_array(fdcoll); call zero_array(glze)
    ltmp = 0.0; ftmp = 0.0
    etmp = 0.0; emax = 0.0; etot = 0.0; eavg = 0.0; edenom = 1.0; ltmax = 1.0

    ! convert gnew from g to h
    call g_adjust(gnew, phi, bpar, direction = from_g_gs2)
    g_work = gnew
    ! convert gnew from h back to g
    call g_adjust(gnew, phi, bpar, direction = to_g_gs2)

    ! map from g_work(ig,isgn,iglo) to glze(il,ilz)
    call gather (lambda_map, g_work, glze)

    ! loop over ig, isign, ik, it, ie, is
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx(lz_lo,ilz)

       if (jend(ig) == 0) then      ! no trapped particles
          ! je = number of + vpa grid pts
          ! te = number of + and - vpa grid pts
          je = ng2
          te = 2*ng2

          ! find d/d(xi) ((1+xi**2)( d g(xi)/ d(xi) )) at each xi
          ! using lagrange (lcoll) and finite difference (fdcoll)
          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,ng2
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=ng2+1, 2*ng2-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2*ng2-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = 2*ng2-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = 2*ng2
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

          ! deal with xi from 1-eps -> eps
          do il=2,ng2
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          ! deal with xi from -eps -> -1+eps
          do il=ng2+1, 2*ng2-1
             fdcoll(il,ilz) = fdb(ig,2*ng2-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,2*ng2-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do

          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(2*ng2,ilz) = -fdf(ig,1)*(glze(2*ng2,ilz) - glze(2*ng2-1,ilz))

       else       ! trapped particle runs          
          je = jend(ig)
          te = 2*je - 1

          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,je
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=je+1, te-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,te-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = te-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = te
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

          ! is il=je handled correctly here?
          do il=2,je
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          do il=je+1,te-1
             fdcoll(il,ilz) = fdb(ig,te-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,te-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do

          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(te,ilz) = -fdf(ig,1)*(glze(te,ilz) - glze(te-1,ilz))

       end if
    end do

    do ilz=lz_lo%llim_world, lz_lo%ulim_world
       ig = ig_idx(lz_lo, ilz)
       ik = ik_idx(lz_lo, ilz)
       it = it_idx(lz_lo, ilz)
       ie = ie_idx(lz_lo, ilz)
       is = is_idx(lz_lo, ilz)
       je = jend(ig)

       if (je == 0) then
          te = 2*ng2
       else
          te = 2*je-1
       end if

       if (idx_local (lz_lo, ilz)) then
          if (proc0) then 
             ltmp = lcoll(:,ilz)
             ftmp = fdcoll(:,ilz)
          else
             call send (lcoll(:,ilz), 0)
             call send (fdcoll(:,ilz), 0)
          endif
       else if (proc0) then
          call receive (ltmp, proc_id(lz_lo, ilz))
          call receive (ftmp, proc_id(lz_lo, ilz))
       endif
       call barrier

       do il=1,te
          etmp = abs(ltmp(il) - ftmp(il))

          if (etmp > emax) then
             emax = etmp
             ltmax = cabs(ltmp(il))
             ikmax = ik
             itmax = it
             iemax = ie
             ismax = is
             ilmax = il
             igmax = ig
          end if

          etot = etot + etmp
          edenom = edenom + abs(ltmp(il))
       end do
    end do

    eavg = etot/edenom
    emax = emax/ltmax

    if (proc0) then
       write(unit,"((1x,e13.6),6(i8),2(1x,e13.6))") user_time, &
            igmax, ikmax, itmax, iemax, ilmax, ismax, emax, eavg
    end if

    deallocate (lcoll, fdcoll, glze, ltmp, ftmp)

  end subroutine collision_error

  !> Error estimate obtained by comparing standard integral with less-accurate integral
  !>
  !> Estimate of the (1) absolute and (2) relative errors resulting from velocity space
  !> integrals in the calculation of the following quantities in the given dimensions: (1)
  !> k phi, energy (2) k phi, untrapped pitch angles (3) k phi, trapped pitch angles, (4)
  !> k apar, energy, (5) k apar, untrapped angles. Relative errors should be < 0.1.
  !>
  !> @note Should this be extended to include error estimates on bpar as well?
  subroutine get_verr (errest, erridx, phi, bpar)
    use le_grids, only: integrate_species, eint_error, trap_error, lint_error, wdim
    use le_grids, only: ng2, nlambda, new_trap_int, grid_has_trapped_particles
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky, aky, akx
    use species, only: nspec, spec
    use dist_fn_arrays, only: gnew, aj0, vpa, g_adjust, g_work, to_g_gs2, from_g_gs2
    use run_parameters, only: has_phi, has_apar, beta
    use gs2_layouts, only: g_lo
    use collisions, only: adjust_vnmult
    use array_utils, only: zero_array
    implicit none
    !> Indices of maximum error.
    integer, dimension (5,3), intent (out) :: erridx
    !> Estimated error.
    real, dimension (5,2), intent (out) :: errest
    !> Electrostatic potential and parallel magnetic field
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    integer :: ig, it, ik, iglo, isgn, ntrap
    complex, dimension (:,:,:), allocatable :: phi_app, apar_app
    complex, dimension (:,:,:,:), allocatable :: phi_e, phi_l, phi_t, apar_e, apar_l, apar_t
    real, dimension (:,:), allocatable :: kmax
    real, dimension (:), allocatable :: wgt
    real, dimension(2) :: errtmp
    integer, dimension(3) :: idxtmp
    logical, parameter :: trap_flag = .true.
    logical :: compute_trapped_error

    allocate(wgt(nspec))

    if (has_phi) then
       allocate(phi_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(phi_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(phi_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

    if (has_apar) then
       allocate(apar_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(apar_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(apar_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

    compute_trapped_error = grid_has_trapped_particles() .and. new_trap_int

    ! first call to g_adjust converts gyro-averaged dist. fn. (g)
    ! into nonadiabatic part of dist. fn. (h)
    call g_adjust (gnew, phi, bpar, direction = from_g_gs2)

    ! take gyro-average of h at fixed total position (not g.c. position)
    ! Note here phi_app and apar_app are effectively antot and antota that
    ! would be returned from a call to getan_from_dnf(gnew, antot, antota,...)
    if (has_phi) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig = -ntgrid, ntgrid
                g_work(ig, isgn, iglo) = aj0(ig, iglo) * gnew(ig, isgn, iglo)
             end do
          end do
       end do

       wgt = spec%z*spec%dens
       call integrate_species (g_work, wgt, phi_app)

       ! integrates dist fn of each species over v-space
       ! after dropping an energy grid point and returns
       ! phi_e, which contains the integral approximations
       ! to phi for each point dropped
       call eint_error (g_work, wgt, phi_e)

       ! integrates dist fn of each species over v-space
       ! after dropping an untrapped lambda grid point and returns phi_l.
       ! phi_l contains ng2 approximations for the integral over lambda that
       ! come from dropping different pts from the gaussian quadrature grid
       call lint_error (g_work, wgt, phi_l)

       ! next loop gets error estimate for trapped particles, if there are any
       if (compute_trapped_error) then
          ntrap = nlambda - ng2
          allocate(phi_t(-ntgrid:ntgrid, ntheta0, naky, ntrap))
          call zero_array(phi_t)
          call trap_error (g_work, wgt, phi_t)
       end if

    end if

    if (has_apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig = -ntgrid, ntgrid
                g_work(ig, isgn, iglo) = aj0(ig, iglo) * vpa(ig, isgn, iglo) * &
                     gnew(ig, isgn, iglo)
             end do
          end do
       end do

       wgt = 2.0 *beta * spec%z * spec%dens * sqrt(spec%temp / spec%mass)
       call integrate_species (g_work, wgt, apar_app)

       call eint_error (g_work, wgt, apar_e)
       call lint_error (g_work, wgt, apar_l)
       if (compute_trapped_error) then
          ntrap = nlambda - ng2
          allocate(apar_t(-ntgrid:ntgrid, ntheta0, naky, ntrap))
          call zero_array(apar_t)
          call trap_error (g_work, wgt, apar_t)
       end if
    end if
    deallocate (wgt)

    ! second call to g_adjust converts from h back to g
    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    allocate (kmax(ntheta0, naky))
    do ik = 1, naky
       do it = 1, ntheta0
          kmax(it,ik) = max(akx(it),aky(ik))
       end do
    end do

    errest = 0.0
    erridx = 0

    if (has_phi) then

       call estimate_error (phi_app, phi_e, kmax, errtmp, idxtmp)
       errest(1,:) = errtmp
       erridx(1,:) = idxtmp

       call estimate_error (phi_app, phi_l, kmax, errtmp, idxtmp)
       errest(2,:) = errtmp
       erridx(2,:) = idxtmp

       if (compute_trapped_error) then
          call estimate_error (phi_app, phi_t, kmax, errtmp, idxtmp, trap_flag)
          errest(3,:) = errtmp
          erridx(3,:) = idxtmp
          deallocate (phi_t)
       end if
       deallocate(phi_app, phi_e, phi_l)
    end if

    ! No trapped errors for apar?
    if (has_apar) then

       call estimate_error (apar_app, apar_e, kmax, errtmp, idxtmp)
       errest(4,:) = errtmp
       erridx(4,:) = idxtmp

       call estimate_error (apar_app, apar_l, kmax, errtmp, idxtmp)
       errest(5,:) = errtmp
       erridx(5,:) = idxtmp

       deallocate(apar_app, apar_e, apar_l)
    end if

    call adjust_vnmult(errest, compute_trapped_error)
  end subroutine get_verr

  !> FIXME : Add documentation
  subroutine estimate_error (app1, app2, kmax_local, errest, erridx, trap)
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use le_grids, only: ng2, jend, il_is_trapped
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: app1
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: app2
    real, dimension (:, :), intent (in) :: kmax_local
    logical, intent (in), optional :: trap
    real, dimension (:), intent (out) :: errest
    integer, dimension (:), intent (out) :: erridx
    integer :: ik, it, ig, ipt, end_idx
    integer :: igmax, ikmax, itmax, gpcnt
    real :: gdsum, gdmax, gpavg, gnsum, gsmax, gpsum, gptmp

    igmax = 0; ikmax = 0; itmax = 0
    gdsum = 0.0; gdmax = 0.0; gpavg = 0.0; gnsum = 0.0; gsmax = 1.0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             if (.not. present(trap) .or. il_is_trapped(jend(ig))) then
                gpcnt = 0; gpsum = 0.0

                if (present(trap)) then
                   end_idx = jend(ig)-ng2
                else
                   end_idx = size(app2, 4)
                end if

                do ipt=1,end_idx
                   gptmp = kmax_local(it, ik) * abs(app1(ig, it, ik) - app2(ig, it, ik, ipt))
                   gpsum = gpsum + gptmp
                   gpcnt = gpcnt + 1
                end do

                gpavg = gpsum/gpcnt

                if (gpavg > gdmax) then
                   igmax = ig
                   ikmax = ik
                   itmax = it
                   gdmax = gpavg
                   gsmax = kmax_local(it, ik) * abs(app1(ig, it, ik))
                end if

                gnsum = gnsum + gpavg
                gdsum = gdsum + kmax_local(it, ik) * abs(app1(ig, it, ik))
             end if
          end do
       end do
    end do

    gdmax = gdmax / gsmax

    erridx(1) = igmax
    erridx(2) = ikmax
    erridx(3) = itmax
    errest(1) = gdmax

    ! Guard against divide by zero. This can commonly occur on the first time step
    ! in runs with fapar /= 0 due to our initial conditions typically ensuring that
    ! the initial apar is identically zero.
    if (abs(gdsum) < epsilon(0.0)) then
       errest(2) = 0.0
    else
       errest(2) = gnsum / gdsum
    end if
  end subroutine estimate_error

  !> Error estimate based on monitoring amplitudes of legendre polynomial coefficients.
  !>
  !> FIXME: finish documentation
  subroutine get_gtran (geavg, glavg, gtavg, phi, bpar)
    use le_grids, only: legendre_transform, nlambda, ng2, nesub, jend, grid_has_trapped_particles, il_is_trapped
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: gnew, aj0, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_layouts, only: g_lo
    use mp, only: proc0, broadcast
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (out) :: geavg, glavg, gtavg
    real, dimension (:), allocatable :: gne2, gnl2, gnt2
    complex, dimension (:,:,:,:,:), allocatable :: getran, gltran, gttran
    real :: genorm, gemax, genum, gedenom
    real :: glnorm, glmax, glnum, gldenom
    real :: gtnorm, gtmax, gtnum, gtdenom
    integer :: ig, it, ik, is, ie, il, iglo, isgn

    allocate(gne2(0:nesub-1))
    allocate(gnl2(0:ng2-1))

    genorm = 0.0 ; glnorm = 0.0
    gne2  = 0.0 ; gnl2 = 0.0
    gemax = 0.0 ; glmax = 0.0
    genum = 0.0 ; gedenom = 0.0
    glnum = 0.0 ; gldenom = 0.0
    gtnum = 0.0 ; gtdenom = 0.0

    allocate(getran(0:nesub-1, -ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate(gltran(0:ng2-1, -ntgrid:ntgrid, ntheta0, naky, nspec))

    call zero_array(getran); call zero_array(gltran)

    if (grid_has_trapped_particles()) then
       allocate(gnt2(0:2*(nlambda-ng2-1)))
       allocate(gttran(0:2*(nlambda-ng2-1), -ntgrid:ntgrid, ntheta0, naky, nspec))
       gtnorm = 0.0 ; gnt2 = 0.0 ; gtmax = 0.0 ; call zero_array(gttran)
    end if

    ! transform from g to h
    call g_adjust (gnew, phi, bpar, direction = from_g_gs2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g_work(ig, isgn, iglo) = aj0(ig, iglo) * gnew(ig, isgn, iglo)
          end do
       end do
    end do

    ! transform from h back to g
    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    ! perform legendre transform on dist. fn. to obtain coefficients
    ! used when expanding dist. fn. in legendre polynomials
    if (grid_has_trapped_particles()) then
       call legendre_transform (g_work, getran, gltran, gttran)
    else
       call legendre_transform (g_work, getran, gltran)
    end if

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   do ie = 0, nesub-1
                      gne2(ie) = real(getran(ie,ig,it,ik,is)*conjg(getran(ie,ig,it,ik,is)))
                   end do

                   do il=0,ng2-1
                      gnl2(il) = real(gltran(il,ig,it,ik,is)*conjg(gltran(il,ig,it,ik,is)))
                   end do
                   genorm = maxval(gne2)
                   if (nesub < 3) then
                      gemax = gne2(size(gne2)-1)
                   else
                      gemax = maxval(gne2(nesub-3:nesub-1))
                   end if
                   glnorm = maxval(gnl2)
                   glmax = maxval(gnl2(ng2-3:ng2-1))

                   genum = genum + gemax
                   gedenom = gedenom + genorm
                   glnum = glnum + glmax
                   gldenom = gldenom + glnorm

                   if (grid_has_trapped_particles()) then
                      do il=0, 2*(jend(ig)-ng2-1)
                         gnt2(il) = real(gttran(il,ig,it,ik,is)*conjg(gttran(il,ig,it,ik,is)))
                      end do
                      gtnorm = maxval(gnt2(0:2*(jend(ig)-ng2-1)))
                      if (il_is_trapped(jend(ig))) then
                         gtmax = maxval(gnt2(2*(jend(ig)-ng2-1)-2:2*(jend(ig)-ng2-1)))
                      else
                         gtmax = gnt2(0)
                      end if

                      gtnum = gtnum + gtmax
                      gtdenom = gtdenom + gtnorm
                   end if
                end do
             end do
          end do
       end do
       geavg = genum / gedenom
       glavg = glnum / gldenom
       if (grid_has_trapped_particles()) gtavg = gtnum/gtdenom
    end if

    call broadcast (geavg)
    call broadcast (glavg)
    if (grid_has_trapped_particles()) then
       call broadcast (gtavg)
       deallocate (gnt2, gttran)
    end if

    deallocate(gne2, gnl2)
    deallocate(getran, gltran)

  end subroutine get_gtran
end module diagnostics_velocity_space
