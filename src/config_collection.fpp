!> Provides a defined type that contains all possible config objects.
!> This can be used to represent the configuration used for a particular
!> run etc.
module config_collection
  use antenna, only: driver_config_type
  use antenna, only: stir_config_type
  use ballstab, only: ballstab_config_type
  use collisions, only: collisions_config_type
  use dist_fn, only: dist_fn_config_type
  use dist_fn, only: dist_fn_species_config_type
#ifdef WITH_EIG
  use eigval, only: eigval_config_type
#endif
  use fields, only: fields_config_type
  use gs2_diagnostics, only: gs2_diagnostics_config_type
  use gs2_init, only: init_config_type
  use gs2_layouts, only: layouts_config_type
  use hyper, only: hyper_config_type
  use ingen_mod, only: ingen_config_type
  use init_g, only: init_g_config_type
  use kt_grids, only: kt_grids_config_type
  use kt_grids_box, only: kt_grids_box_config_type
  use kt_grids_range, only: kt_grids_range_config_type
  use kt_grids_single, only: kt_grids_single_config_type
  use kt_grids_specified, only: kt_grids_specified_config_type
  use kt_grids_specified, only: kt_grids_specified_element_config_type
  use le_grids, only: le_grids_config_type
#ifdef NEW_DIAG
  use diagnostics_config, only: new_gs2_diagnostics_config_type
#endif
  use nonlinear_terms, only: nonlinear_terms_config_type
  use normalisations, only: normalisations_config_type
  use optimisation_config, only: optimisation_config_type
  use run_parameters, only: knobs_config_type
  use run_parameters, only: parameters_config_type
  use parameter_scan, only: parameter_scan_config_type
  use gs2_reinit, only: reinit_config_type
  use dist_fn, only: source_config_type
  use species, only: species_config_type
  use species, only: species_element_config_type
  use theta_grid, only: theta_grid_config_type
  use theta_grid_eik, only: theta_grid_eik_config_type
  use theta_grid_file, only: theta_grid_file_config_type
  use theta_grid_gridgen, only: theta_grid_gridgen_config_type
  use theta_grid_params, only: theta_grid_parameters_config_type
  use theta_grid_salpha, only: theta_grid_salpha_config_type

  implicit none
  private

  public :: gs2_config_type

  type :: gs2_config_type
     type(ballstab_config_type) :: ballstab_config
     type(collisions_config_type) :: collisions_config
     type(dist_fn_config_type) :: dist_fn_config
     type(dist_fn_species_config_type), dimension(:), allocatable :: dist_fn_species_config
     type(driver_config_type) :: driver_config
#ifdef WITH_EIG
     type(eigval_config_type) :: eigval_config
#endif
     type(fields_config_type) :: fields_config
     type(gs2_diagnostics_config_type) :: gs2_diagnostics_config
     type(hyper_config_type) :: hyper_config
     type(ingen_config_type) :: ingen_config
     type(init_config_type) :: init_config
     type(init_g_config_type) :: init_g_config
     type(knobs_config_type) :: knobs_config
     type(kt_grids_config_type) :: kt_grids_config
     type(kt_grids_box_config_type) :: kt_grids_box_config
     type(kt_grids_range_config_type) :: kt_grids_range_config
     type(kt_grids_single_config_type) :: kt_grids_single_config
     type(kt_grids_specified_config_type) :: kt_grids_specified_config
     type(kt_grids_specified_element_config_type), dimension(:), allocatable :: kt_grids_specified_element_config
     type(layouts_config_type) :: layouts_config
     type(le_grids_config_type) :: le_grids_config
#ifdef NEW_DIAG
     type(new_gs2_diagnostics_config_type) :: new_gs2_diagnostics_config
#endif
     type(nonlinear_terms_config_type) :: nonlinear_terms_config
     type(normalisations_config_type) :: normalisations_config
     type(optimisation_config_type) :: optimisation_config
     type(parameters_config_type) :: parameters_config
     type(parameter_scan_config_type) :: parameter_scan_config
     type(reinit_config_type) :: reinit_config
     type(source_config_type) :: source_config
     type(species_config_type) :: species_config
     type(species_element_config_type), dimension(:), allocatable :: species_element_config
     type(stir_config_type), dimension(:), allocatable :: stir_config
     type(theta_grid_config_type) :: theta_grid_config
     type(theta_grid_eik_config_type) :: theta_grid_eik_config
     type(theta_grid_file_config_type) :: theta_grid_file_config
     type(theta_grid_gridgen_config_type) :: theta_grid_gridgen_config
     type(theta_grid_parameters_config_type) :: theta_grid_parameters_config
     type(theta_grid_salpha_config_type) :: theta_grid_salpha_config
   contains
     !> See [[write_to_unit]] for details.
     procedure :: write_to_unit
     !> See [[get_configs]] for details
     procedure :: get_configs
  end type gs2_config_type

contains

  !> Writes out the current configuration state to an optionally provided unit.
  !> If no unit is provided then we write to `output_unit`.
  subroutine write_to_unit(self, unit)
    use iso_fortran_env, only: output_unit
    implicit none
    ! Note we have to make this intent out as we may want/need to allocate
    ! some components that have not been allocated. If so we do deallocate
    ! so self should be unchanged on exit.
    class(gs2_config_type), intent(in out) :: self
    ! The unit to write to, if not given then will use output_unit. Note we
    ! assume this unit has been opened ready for writing.
    integer, intent(in), optional :: unit
    integer :: local_output_unit, i
    logical :: allocated_flag

    ! Handle optional argument
    local_output_unit = output_unit
    if (present(unit)) local_output_unit = unit

    ! Take care of writing
    call self%ballstab_config%write(local_output_unit)
    call self%collisions_config%write(local_output_unit)
    call self%dist_fn_config%write(local_output_unit)

    allocated_flag = allocated(self%dist_fn_species_config)
    if (.not.allocated_flag) allocate(self%dist_fn_species_config(self%species_config%nspec))
    do i = 1, size(self%dist_fn_species_config)
       call self%dist_fn_species_config(i)%write(local_output_unit)
    end do
    if (.not.allocated_flag) deallocate(self%dist_fn_species_config)

    call self%driver_config%write(local_output_unit)
#ifdef WITH_EIG
    call self%eigval_config%write(local_output_unit)
#endif
    call self%fields_config%write(local_output_unit)
#ifdef NEW_DIAG
    if (self%knobs_config%use_old_diagnostics) then
       call self%gs2_diagnostics_config%write(local_output_unit)
    else
       call self%new_gs2_diagnostics_config%write(local_output_unit)
    end if
#else
    call self%gs2_diagnostics_config%write(local_output_unit)
#endif
    call self%hyper_config%write(local_output_unit)
    call self%ingen_config%write(local_output_unit)
    call self%init_config%write(local_output_unit)
    call self%init_g_config%write(local_output_unit)
    call self%knobs_config%write(local_output_unit)
    call self%kt_grids_config%write(local_output_unit)
    call self%kt_grids_box_config%write(local_output_unit)
    call self%kt_grids_range_config%write(local_output_unit)
    call self%kt_grids_single_config%write(local_output_unit)
    call self%kt_grids_specified_config%write(local_output_unit)

    allocated_flag = allocated(self%kt_grids_specified_element_config)
    if (.not.allocated_flag) allocate(self%kt_grids_specified_element_config(max(self%kt_grids_specified_config%naky,self%kt_grids_specified_config%ntheta0)))
    do i = 1, size(self%kt_grids_specified_element_config)
       call self%kt_grids_specified_element_config(i)%write(local_output_unit)
    end do
    if (.not.allocated_flag) deallocate(self%kt_grids_specified_element_config)

    call self%layouts_config%write(local_output_unit)
    call self%le_grids_config%write(local_output_unit)
    call self%nonlinear_terms_config%write(local_output_unit)
    call self%normalisations_config%write(local_output_unit)
    call self%optimisation_config%write(local_output_unit)
    call self%parameters_config%write(local_output_unit)
    call self%parameter_scan_config%write(local_output_unit)
    call self%reinit_config%write(local_output_unit)
    call self%source_config%write(local_output_unit)
    call self%species_config%write(local_output_unit)

    allocated_flag = allocated(self%species_element_config)
    if (.not.allocated_flag) allocate(self%species_element_config(self%species_config%nspec))
    do i = 1, size(self%species_element_config)
       call self%species_element_config(i)%write(local_output_unit)
    end do
    if (.not.allocated_flag) deallocate(self%species_element_config)

    allocated_flag = allocated(self%stir_config)
    if (.not.allocated_flag) allocate(self%stir_config(self%driver_config%nk_stir))
    do i = 1, size(self%stir_config)
       call self%stir_config(i)%write(local_output_unit)
    end do
    if (.not.allocated_flag) deallocate(self%stir_config)

    call self%theta_grid_config%write(local_output_unit)
    call self%theta_grid_eik_config%write(local_output_unit)
    call self%theta_grid_file_config%write(local_output_unit)
    call self%theta_grid_gridgen_config%write(local_output_unit)
    call self%theta_grid_parameters_config%write(local_output_unit)
    call self%theta_grid_salpha_config%write(local_output_unit)

  end subroutine write_to_unit

  !> Copy all the module-level configs into \p self.
  !>
  !> Can be used to get a "snapshot" of the current state of all
  !> module's configs.
  subroutine get_configs(self)
    use antenna, only : get_antenna_driver_config, get_antenna_stir_config
    use ballstab, only: get_ballstab_config
    use collisions, only: get_collisions_config
    use dist_fn, only: get_dist_fn_config, get_dist_fn_species_config, get_source_config
#ifdef WITH_EIG
    use eigval, only: get_eigval_config
#endif
    use fields, only: get_fields_config
    use gs2_diagnostics, only: get_gs2_diagnostics_config
    use gs2_init, only: get_gs2_init_config
    use gs2_layouts, only: get_gs2_layouts_config
    use hyper, only: get_hyper_config
    use ingen_mod, only: get_ingen_mod_config
    use init_g, only: get_init_g_config
    use kt_grids, only: get_kt_grids_config
    use kt_grids_box, only: get_kt_grids_box_config
    use kt_grids_range, only: get_kt_grids_range_config
    use kt_grids_single, only: get_kt_grids_single_config
    use kt_grids_specified, only: get_kt_grids_specified_config, get_kt_grids_specified_element_config
    use le_grids, only: get_le_grids_config
#ifdef NEW_DIAG
    use diagnostics_config, only: get_new_gs2_diagnostics_config => get_diagnostics_config_config
#endif
    use nonlinear_terms, only: get_nonlinear_terms_config
    use normalisations, only: get_normalisations_config
    use optimisation_config, only: get_optimisation_config_config
    use run_parameters, only: get_knobs_config, get_parameters_config
    use parameter_scan, only: get_parameter_scan_config
    use gs2_reinit, only: get_gs2_reinit_config
    use species, only: get_species_config, get_species_element_config
    use theta_grid, only: get_theta_grid_config
    use theta_grid_eik, only: get_theta_grid_eik_config
    use theta_grid_file, only: get_theta_grid_file_config
    use theta_grid_gridgen, only: get_theta_grid_gridgen_config
    use theta_grid_params, only: get_theta_grid_params_config
    use theta_grid_salpha, only: get_theta_grid_salpha_config

    class(gs2_config_type), intent(in out) :: self

    self%driver_config = get_antenna_driver_config()
    self%stir_config = get_antenna_stir_config()

    self%ballstab_config = get_ballstab_config()
    self%collisions_config = get_collisions_config()
    self%dist_fn_config = get_dist_fn_config()
    self%dist_fn_species_config = get_dist_fn_species_config()
#ifdef WITH_EIG
    self%eigval_config = get_eigval_config()
#endif
    self%fields_config = get_fields_config()
    self%gs2_diagnostics_config = get_gs2_diagnostics_config()
    self%init_config = get_gs2_init_config()
    self%layouts_config = get_gs2_layouts_config()
    self%hyper_config = get_hyper_config()
    self%ingen_config = get_ingen_mod_config()
    self%init_g_config = get_init_g_config()
    self%kt_grids_config = get_kt_grids_config()
    self%kt_grids_box_config = get_kt_grids_box_config()
    self%kt_grids_range_config = get_kt_grids_range_config()
    self%kt_grids_single_config = get_kt_grids_single_config()
    self%kt_grids_specified_config = get_kt_grids_specified_config()
    self%kt_grids_specified_element_config = get_kt_grids_specified_element_config()
    self%le_grids_config = get_le_grids_config()
#ifdef NEW_DIAG
    self%new_gs2_diagnostics_config = get_new_gs2_diagnostics_config()
#endif
    self%nonlinear_terms_config = get_nonlinear_terms_config()
    self%normalisations_config = get_normalisations_config()
    self%optimisation_config = get_optimisation_config_config()
    self%knobs_config = get_knobs_config()
    self%parameters_config = get_parameters_config()
    self%parameter_scan_config = get_parameter_scan_config()
    self%reinit_config = get_gs2_reinit_config()
    self%source_config = get_source_config()
    self%species_config = get_species_config()
    self%species_element_config = get_species_element_config()
    self%theta_grid_config = get_theta_grid_config()
    self%theta_grid_eik_config = get_theta_grid_eik_config()
    self%theta_grid_file_config = get_theta_grid_file_config()
    self%theta_grid_gridgen_config = get_theta_grid_gridgen_config()
    self%theta_grid_parameters_config = get_theta_grid_params_config()
    self%theta_grid_salpha_config = get_theta_grid_salpha_config()

  end subroutine get_configs

end module config_collection
