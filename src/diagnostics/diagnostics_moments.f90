!> A module for writing out moments of the distribution function,
!! for example density, temperature, parallel flow
module diagnostics_moments
  implicit none
  private
  public :: write_moments, write_full_moments_notgc
contains

  !> FIXME : Add documentation  
  subroutine write_moments(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use dist_fn, only: getmoms
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, qparflux, pperpj1, qpperpj1

    if (.not. gnostics%replay) then
       call getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    else
       ntot = 0. ; density = 0. ; upar = 0. ; tpar = 0. ; tperp = 0.; qparflux = 0. ; pperpj1 = 0. ; qpperpj1 = 0.
    end if

    call write_standard_moment_properties(gnostics, &
         'ntot',  'The total perturbed species density', 'n_r', ntot, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'density',  'The non-adiabatic part of the perturbed species density', 'n_r', density, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'upar',  'The perturbed parallel flow', 'v_ths', upar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tpar',  'The perturbed parallel temperature', 'T_r', tpar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tperp',  'The perturbed perpendicular temperature', 'T_r', tperp, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qparflux',  'The parallel heat flux', 'n_r T_r v_thr', qparflux, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'pperpj1',  'A modified perpendicular pressure which give the particle flux from bpar',&
         'n_r T_r / e', pperpj1, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qpperpj1',  'A modified perpendicular pressure * energy which give the heat flux from bpar',&
         'n_r T_r^2 / e', qpperpj1, gnostics%distributed)
  end subroutine write_moments

  !> FIXME : Add documentation  
  subroutine write_full_moments_notgc(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use dist_fn, only: getmoms_notgc
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp !, qparflux, pperpj1, qpperpj1

    if (.not. gnostics%replay) then
       call getmoms_notgc(phinew, bparnew, density,upar,tpar,tperp,ntot)
    else
       ntot = 0. ; density = 0. ; upar = 0. ; tpar = 0. ; tperp = 0.
    end if

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
end module diagnostics_moments
