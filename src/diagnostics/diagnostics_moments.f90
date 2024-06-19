!> A module for writing out moments of the distribution function,
!! for example density, temperature, parallel flow
module diagnostics_moments
  implicit none
  private
  public :: write_moments, write_full_moments_notgc, getmoms, getmoms_notgc

!#ifdef NETCDF_PARALLEL
!  logical, parameter :: moment_to_allprocs = .true.
!#else
  logical, parameter :: moment_to_allprocs = .false.
!#endif

contains

  !> FIXME : Add documentation  
  subroutine write_moments(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, qparflux, pperpj1, qpperpj1

    call getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    call write_standard_moment_properties(gnostics, &
         'ntot',  'The total perturbed species density', 'n_ref', ntot, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'density',  'The non-adiabatic part of the perturbed species density', 'n_ref', density, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'upar',  'The perturbed parallel flow', 'v_ths', upar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tpar',  'The perturbed parallel temperature', 'T_ref', tpar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tperp',  'The perturbed perpendicular temperature', 'T_ref', tperp, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qparflux',  'The parallel heat flux', 'n_ref T_ref v_thr', qparflux, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'pperpj1',  'A modified perpendicular pressure which give the particle flux from bpar',&
         'n_ref T_ref / q_ref', pperpj1, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qpperpj1',  'A modified perpendicular pressure * energy which give the heat flux from bpar',&
         'n_ref T_ref^2 / q_ref', qpperpj1, gnostics%distributed)
  end subroutine write_moments

  !> FIXME : Add documentation  
  subroutine write_full_moments_notgc(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, upar, tpar, tperp

    call getmoms_notgc(phinew, bparnew, density,upar,tpar,tperp,ntot)

    call write_standard_moment_properties(gnostics, &
         'ntot_notgc',  'The total perturbed species density &
         & in non-guiding centre coordinates', 'n_r', ntot, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'density_notgc',  'The non-adiabatic part of the perturbed species density &
         & in non-guiding centre coordinates', 'n_r', density, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'upar_notgc',  'The perturbed parallel flow &
         & in non-guiding centre coordinates', 'v_ths', upar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tpar_notgc',  'The perturbed parallel temperature &
         & in non-guiding centre coordinates', 'T_r', tpar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tperp_notgc',  'The perturbed perpendicular temperature &
         & in non-guiding centre coordinates', 'T_r', tperp, gnostics%distributed)
  end subroutine write_full_moments_notgc

  !> FIXME : Add documentation  
  subroutine write_standard_moment_properties(gnostics, moment_name, moment_description, &
       moment_units, moment_value, distributed)
    use volume_averages, only: average_theta, average_kx, average_ky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: starts, mom_dims, mom_t_dims, fluxk_dims, netcdf_write_complex, complex_dim, &
         kx_dim, ky_dim, species_dim, time_dim, theta_dim
    use neasyf, only: neasyf_write
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: moment_name, moment_description, moment_units
    complex, dimension(-ntgrid:,:,:,:), intent(inout) :: moment_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky, nspec) :: moment2_by_mode
    real, dimension(naky,nspec) :: moment2_by_ky
    real, dimension(ntheta0,nspec) :: moment2_by_kx
    complex, dimension(ntheta0, naky, nspec) :: moment_by_mode
    complex, dimension(ntheta0, naky, nspec) :: moment_igomega_by_mode
    complex, dimension(ntheta0, nspec) :: moment_flx_surfavg
    logical :: write_moment_by_time
    integer :: it !, ik

    if (.not. gnostics%writing) return

    call average_theta(moment_value, moment_value, moment2_by_mode, distributed)

    call average_kx(moment2_by_mode, moment2_by_ky, distributed)
    call neasyf_write(gnostics%file_id, moment_name//"2_by_ky", moment2_by_ky, &
         dim_names=[ky_dim, species_dim, time_dim], start=starts(3, gnostics%nout), &
         long_name=moment_description//" squared and averaged over theta and kx", &
         units="("//moment_units//")^2")

    call average_ky(moment2_by_mode, moment2_by_kx, distributed)
    call neasyf_write(gnostics%file_id, moment_name//"2_by_kx", moment2_by_kx, &
         dim_names=[kx_dim, species_dim, time_dim], start=starts(3, gnostics%nout), &
         long_name=moment_description//" squared and averaged over theta and ky", &
         units="("//moment_units//")^2")

    call neasyf_write(gnostics%file_id, moment_name//"2", sum(moment2_by_kx, 1), &
         dim_names=[species_dim, time_dim], start=starts(2, gnostics%nout), &
         long_name=moment_description//" squared and averaged over theta, kx and ky", &
         units="("//moment_units//")^2")

    call average_theta(moment_value, moment_by_mode, gnostics%distributed)
    moment_igomega_by_mode(:,:,:) = moment_value(gnostics%igomega,:,:,:)

    ! moment_by_mode could be distributed so we have to be careful here
    moment_flx_surfavg(:,:) = 0.0
    do it = 1,ntheta0
       if ((.not. gnostics%distributed).or.field_k_local(it, 1))&
            moment_flx_surfavg(it, :) = moment_by_mode(it, 1, :)
    end do
    if (gnostics%distributed) call sum_allreduce(moment_flx_surfavg)
    
    call netcdf_write_complex(gnostics%file_id, moment_name//"_flxsurf_avg", moment_flx_surfavg, &
         dim_names=[complex_dim, kx_dim, species_dim, time_dim], start=starts(4, gnostics%nout), &
         long_name=moment_description//" flux surface averaged over theta, at ky=0 (actally ik==1)", &
         units=moment_units)
    
    call netcdf_write_complex(gnostics%file_id, moment_name, moment_value, &
         dim_names=[complex_dim, theta_dim, kx_dim, ky_dim, species_dim], &
         long_name=moment_description, units=moment_units)
    call netcdf_write_complex( &
         gnostics%file_id, moment_name//"_igomega_by_mode", moment_igomega_by_mode, &
         dim_names=mom_dims, start=starts(5, gnostics%nout), &
         long_name=moment_description//" at ig=igomega" , units=moment_units)
    call netcdf_write_complex(gnostics%file_id, moment_name//"_by_mode", moment_by_mode, &
         dim_names=mom_dims, start=starts(5, gnostics%nout), &
         long_name=moment_description//"  averaged over theta" , &
         units=moment_units)
    call neasyf_write(gnostics%file_id, moment_name//"2_by_mode", moment2_by_mode, &
         dim_names=fluxk_dims, start=starts(4, gnostics%nout), &
         long_name=moment_description//" squared and averaged over theta" , &
         units=moment_units)

    write_moment_by_time = .false.
    if (moment_name .eq. 'ntot'  .and. gnostics%write_ntot_over_time ) write_moment_by_time = .true.
    if (moment_name .eq. 'density' .and. gnostics%write_density_over_time) write_moment_by_time = .true.
    if (moment_name .eq. 'upar' .and. gnostics%write_upar_over_time) write_moment_by_time = .true.
    if (moment_name .eq. 'tperp' .and. gnostics%write_tperp_over_time) write_moment_by_time = .true.

    if (write_moment_by_time) then
      call netcdf_write_complex(gnostics%file_id, moment_name//"_t", moment_value, &
           dim_names=mom_t_dims, start=starts(6, gnostics%nout), &
           long_name=moment_description//": the whole moment" , &
           units=moment_units)
    end if
  end subroutine write_standard_moment_properties

  !> Calculate various moments of the distribution function
  subroutine getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_adjust, to_g_gs2, from_g_gs2, g_work
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, energy
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    !> Total perturbed density
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: ntot
    !> Non-adiabatic part of the perturbed density. Normalised to \(n_s n_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density
    !> Parallel velocity. Normalised to \(v_{t,s} = v_{t,ref}\sqrt{T_s / m_s}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: upar
    !> Parallel temperature. Normalised to \(T_s T_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tpar
    !> Perpendicular temperature. Normalised to \(T_s T_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp
    !> Parallel heat flux. Normalised to \(n_s n_{ref} T_s T_{ref} v_{t,s}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: qparflux
    !> Pressure part of the particle flux. Normalised to \(\frac{n_s T_s}{q_s} \frac{n_{ref} T_{ref}}{q_{ref}}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: pperpj1
    !> Pressure part of the heat flux. Normalised to \(\frac{n_s T_s^2}{q_s} \frac{n_{ref} T_{ref}^2}{q_{ref}}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: qpperpj1
    !> Electrostatic potential, \(\phi\)
    complex, dimension (-ntgrid:,:,:), intent(in) :: phinew
    !> Parallel magnetic field, \(B_\parallel\)
    complex, dimension (-ntgrid:,:,:), intent(in) :: bparnew

    integer :: ik, it, isgn, ie, is, iglo

! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
!                   from gnew, phinew and bparnew.
!           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
!           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
!
! use g_work as dist_fn dimensioned working space for all moments
! (avoid making more copies of gnew to save memory!)
!
! set gnew = g_wesson, but return gnew to entry state prior to exit
    call g_adjust(gnew,phinew,bparnew, direction = from_g_gs2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate moments over the nonadiabatic part of <delta_f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g_work= J0 g_wesson = nonadiabatic piece of <delta_f>/F_m
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

! CMR: density is the nonadiabatic piece of perturbed density
! NB normalised wrt equ'm density for species s: n_s n_ref
!    ie multiply by (n_s n_ref) to get abs density pert'n
    call integrate_moment (g_work, density,  moment_to_allprocs, full_arr)

! DJA/CMR: upar and tpar moments
! (nb adiabatic part of <delta f> does not contribute to upar, tpar or tperp)
! NB UPAR is normalised to vt_s = sqrt(T_s/m_s) vt_ref
!    ie multiply by spec(is)%stm * vt_ref to get abs upar
    g_work = vpa*g_work
    call integrate_moment (g_work, upar,  moment_to_allprocs, full_arr)

    g_work = 2.*vpa*g_work
    call integrate_moment (g_work, tpar,  moment_to_allprocs, full_arr)
! tpar transiently stores ppar, nonadiabatic perturbed par pressure
!      vpa normalised to: sqrt(2 T_s T_ref/m_s m_ref)
!  including factor 2 in g_work product ensures
!     ppar normalised to: n_s T_s n_ref T_ref
!                         ppar = tpar + density, and so:
    tpar = tpar - density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = vperp2(:,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)
       end do
    end do
    call integrate_moment (g_work, tperp,  moment_to_allprocs, full_arr)
! tperp transiently stores pperp, nonadiabatic perturbed perp pressure
!                          pperp = tperp + density, and so:
    tperp = tperp - density
! NB TPAR, and TPERP are normalised by T_s T_ref
!    ie multiply by T_s T_ref to get abs TPAR, TPERP

! Now compute QPARFLUX
! NB QPARFLUX is normalised to n_s n_ref T_s T_ref v_ts
!    ie multiply by (n_s T_s spec(is)%stm) n_ref T_ref vt_ref to get abs qparflux
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = vpa(:,isgn,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
    end do
    call integrate_moment (g_work, qparflux,  moment_to_allprocs, full_arr)

! Now compute PPERPJ1, a modified p_perp which gives particle flux from Bpar
! NB PPERPJ1 is normalised to (n_s T_s/q_s)  n_ref T_ref/q_ref
!    ie multiply by (n_s spec(is)%tz) n_ref T_ref/q_ref to get abs PPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) &
               = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
       end do
    end do
    call integrate_moment (g_work, pperpj1,  moment_to_allprocs, full_arr)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
    end do
    call integrate_moment (g_work, qpperpj1, moment_to_allprocs, full_arr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g_work = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - phinew(:,it,ik)*spec(is)%zt*nonmaxw_corr(ie,is)
       end do
    end do

! total perturbed density
    call integrate_moment (g_work, ntot,  moment_to_allprocs, full_arr)

!CMR, now multiply by species dependent normalisations to leave only reference normalisations
    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       density(:,:,:,is)=density(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
       qparflux(:,:,:,is)=qparflux(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%stm
       pperpj1(:,:,:,is)=pperpj1(:,:,:,is)*spec(is)%dens*spec(is)%tz
       qpperpj1(:,:,:,is)=qpperpj1(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%tz
    end do

! return gnew to its initial state, the variable evolved in GS2
    call g_adjust(gnew,phinew,bparnew, direction = to_g_gs2)
  end subroutine getmoms

  !> Calculates moments at not guiding center coordinate
  subroutine getmoms_notgc (phinew, bparnew, dens, upar, tpar, tper, ntot, jpar, reset_coefficients)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_work
    use dist_fn, only: init_mom_coeff
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx
    use species, only: nspec, spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx => ntheta0, naky
    use le_grids, only: integrate_moment
    use optionals, only: get_option_with_default
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, intent (out) :: &
         & dens(-ntgrid:,:,:,:), upar(-ntgrid:,:,:,:), &
         & tpar(-ntgrid:,:,:,:), tper(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ntot(-ntgrid:,:,:,:)
    complex, intent (out), optional :: jpar(-ntgrid:,:,:)
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew, bparnew
    logical, intent(in), optional :: reset_coefficients
    real, dimension(:, :, :, :), allocatable, save :: mom_coeff
    real, dimension(:, :, :), allocatable, save :: mom_coeff_npara, mom_coeff_nperp
    real, dimension(:, :, :), allocatable, save :: mom_coeff_tpara, mom_coeff_tperp
    real, dimension(:, :, :), allocatable, save :: mom_shift_para, mom_shift_perp

    integer :: isgn, iglo, is, ie, it, ik, ig
    complex :: tpar2, tper2
    real :: a, b

    ! returns moment integrals to PE 0

    if (get_option_with_default(reset_coefficients, .false.)) then
       if (allocated(mom_coeff)) deallocate(mom_coeff)
       if (allocated(mom_coeff_npara)) deallocate(mom_coeff_npara)
       if (allocated(mom_coeff_nperp)) deallocate(mom_coeff_nperp)
       if (allocated(mom_coeff_tpara)) deallocate(mom_coeff_tpara)
       if (allocated(mom_coeff_tperp)) deallocate(mom_coeff_tperp)
       if (allocated(mom_shift_para)) deallocate(mom_shift_para)
       if (allocated(mom_shift_perp)) deallocate(mom_shift_perp)
    end if

    ! Note we mark mom_coeff* as 'save' so that the result of this call is cached
    ! and not recalculated each time. This can help save overhead and should be
    ! fine. The result of this call should only depend on the velocity grids, the
    ! wavenumbers and the species parameters (and the derived Bessel functions).
    ! These should not change during a single standard run, but if we wish to enable
    ! more exotic behaviour we might want to remove the save from above declarations and
    ! recalculate these coefficients each time. For now just add an optional flag to
    ! request recalculation.
    call init_mom_coeff(mom_coeff, mom_coeff_npara, mom_coeff_nperp, mom_coeff_tpara, &
       mom_coeff_tperp, mom_shift_para, mom_shift_perp)

    ! not guiding center n_total
    if(present(ntot)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)

          do isgn = 1, 2
             g_work(:,isgn,iglo) = aj0(:,iglo) * gnew(:,isgn,iglo)
          end do
          do isgn = 1, 2
             g_work(:,isgn,iglo) = g_work(:,isgn,iglo) + phinew(:,it,ik) &
                  & *(aj0(:,iglo)**2-1.0) * spec(is)%zt * nonmaxw_corr(ie,is)
          end do
          do isgn = 1, 2
             g_work(:,isgn,iglo) = g_work(:,isgn,iglo) &
                  & + 2.*vperp2(:,iglo)*aj1(:,iglo)*aj0(:,iglo) &
                  & * bparnew(:,it,ik) * nonmaxw_corr(ie,is)
          end do
       end do
       call integrate_moment (g_work, ntot,  moment_to_allprocs,full_arr)
    endif

! not guiding center density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, dens,  moment_to_allprocs,full_arr)

! not guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, upar,  moment_to_allprocs,full_arr)

! not guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = 2.*aj0(:,iglo)*vpa(:,isgn,iglo)**2*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, tpar,  moment_to_allprocs,full_arr)

! not guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = 2.*vperp2(:,iglo)*aj1(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, tper,  moment_to_allprocs,full_arr)

    do ig=-ntgrid,ntgrid
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                tpar2=tpar(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_npara(it,ik,is)
                tper2=tper(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_nperp(it,ik,is)

                a=mom_coeff_tperp(it,ik,is)
                b=mom_coeff_tpara(it,ik,is)

                tpar(ig,it,ik,is)=(   tpar2-a*tper2)/(1.-a*b)
                tper(ig,it,ik,is)=(-b*tpar2+  tper2)/(1.-a*b)
             end do
          end do
       end do
    end do

    do is=1,nspec
       dens(:,:,:,is)=dens(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tper(:,:,:,is)=tper(:,:,:,is)*spec(is)%temp
    end do

    if(present(jpar)) then
       jpar(:,:,:)=cmplx(0.,0.)
       do is=1,nspec
          jpar(:,:,:)=jpar(:,:,:)+spec(is)%z*spec(is)%dens*upar(:,:,:,is)
       end do
    endif
  end subroutine getmoms_notgc
end module diagnostics_moments
