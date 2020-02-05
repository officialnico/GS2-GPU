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

  public :: config_collection_type

  type :: config_collection_type
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
  end type config_collection_type
contains

  !> Writes out the current configuration state to an optionally provided unit.
  !> If no unit is provided then we write to `output_unit`.
  subroutine write_to_unit(self, unit)
    use iso_fortran_env, only: output_unit
    implicit none
    ! Note we have to make this intent out as we may want/need to allocate
    ! some components that have not been allocated. If so we do deallocate
    ! so self should be unchanged on exit.
    class(config_collection_type), intent(in out) :: self
    ! The unit to write to, if not given then will use output_unit. Note we
    ! assume this unit has been opened ready for writing.
    integer, intent(in), optional :: unit
    integer :: the_unit, i
    logical :: allocated_flag

    ! Handle optional argument
    the_unit = output_unit
    if (present(unit)) the_unit = unit

    ! Take care of writing
    call self%ballstab_config%write(the_unit)
    call self%collisions_config%write(the_unit)
    call self%dist_fn_config%write(the_unit)

    allocated_flag = allocated(self%dist_fn_species_config)
    if (.not.allocated_flag) allocate(self%dist_fn_species_config(self%species_config%nspec))
    do i = 1, size(self%dist_fn_species_config)
       call self%dist_fn_species_config(i)%write(the_unit)
    end do
    if (.not.allocated_flag) deallocate(self%dist_fn_species_config)

    call self%driver_config%write(the_unit)
#ifdef WITH_EIG
    call self%eigval_config%write(the_unit)
#endif
    call self%fields_config%write(the_unit)
    call self%gs2_diagnostics_config%write(the_unit)
    call self%hyper_config%write(the_unit)
    call self%ingen_config%write(the_unit)
    call self%init_config%write(the_unit)
    call self%init_g_config%write(the_unit)
    call self%knobs_config%write(the_unit)
    call self%kt_grids_config%write(the_unit)
    call self%kt_grids_box_config%write(the_unit)
    call self%kt_grids_range_config%write(the_unit)
    call self%kt_grids_single_config%write(the_unit)
    call self%kt_grids_specified_config%write(the_unit)

    allocated_flag = allocated(self%kt_grids_specified_element_config)
    if (.not.allocated_flag) allocate(self%kt_grids_specified_element_config(max(self%kt_grids_specified_config%naky,self%kt_grids_specified_config%ntheta0)))
    do i = 1, size(self%kt_grids_specified_element_config)
       call self%kt_grids_specified_element_config(i)%write(the_unit)
    end do
    if (.not.allocated_flag) deallocate(self%kt_grids_specified_element_config)

    call self%layouts_config%write(the_unit)
    call self%le_grids_config%write(the_unit)
#ifdef NEW_DIAG
    ! Note if we have new diagnostics then we currently write out both
    ! versions of the namelist -- we probably/may want to check the choice
    ! of `use_old_diagnostics` and then write the appropriate one based on that
    ! ideally there would only be one diagnostics config shared between the two
    ! approaches.
    call self%new_gs2_diagnostics_config%write(the_unit)
#endif
    call self%nonlinear_terms_config%write(the_unit)
    call self%normalisations_config%write(the_unit)
    call self%optimisation_config%write(the_unit)
    call self%parameters_config%write(the_unit)
    call self%parameter_scan_config%write(the_unit)
    call self%reinit_config%write(the_unit)
    call self%source_config%write(the_unit)
    call self%species_config%write(the_unit)

    allocated_flag = allocated(self%species_element_config)
    if (.not.allocated_flag) allocate(self%species_element_config(self%species_config%nspec))
    do i = 1, size(self%species_element_config)
       call self%species_element_config(i)%write(the_unit)
    end do
    if (.not.allocated_flag) deallocate(self%species_element_config)

    allocated_flag = allocated(self%stir_config)
    if (.not.allocated_flag) allocate(self%stir_config(self%driver_config%nk_stir))
    do i = 1, size(self%stir_config)
       call self%stir_config(i)%write(the_unit)
    end do
    if (.not.allocated_flag) deallocate(self%stir_config)

    call self%theta_grid_config%write(the_unit)
    call self%theta_grid_eik_config%write(the_unit)
    call self%theta_grid_file_config%write(the_unit)
    call self%theta_grid_gridgen_config%write(the_unit)
    call self%theta_grid_parameters_config%write(the_unit)
    call self%theta_grid_salpha_config%write(the_unit)

  end subroutine write_to_unit
end module config_collection
