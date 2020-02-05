program write_default_input_file
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: get_unused_unit, open_output_file, close_output_file

  use config_collection, only: config_collection_type
  implicit none

  type(config_collection_type) :: the_config

  logical :: list
  integer :: output_unit
  integer :: i

  call init_mp

  if (proc0) then
     call get_unused_unit(output_unit)
     open(unit = output_unit, action = "write", file = "defaults.in")
     call the_config%ballstab_config%write(output_unit)
     call the_config%collisions_config%write(output_unit)
     call the_config%dist_fn_config%write(output_unit)
     if (allocated(the_config%dist_fn_species_config)) then
        do i = 1, size(the_config%dist_fn_species_config)
           call the_config%dist_fn_species_config(i)%write(output_unit)
        end do
     else
        allocate(the_config%dist_fn_species_config(1))
        call the_config%dist_fn_species_config(1)%write(output_unit)
        deallocate(the_config%dist_fn_species_config)
     end if
     call the_config%driver_config%write(output_unit)
#ifdef WITH_EIG
     call the_config%eigval_config%write(output_unit)
#endif
     call the_config%fields_config%write(output_unit)
     call the_config%gs2_diagnostics_config%write(output_unit)
     call the_config%hyper_config%write(output_unit)
     call the_config%ingen_config%write(output_unit)
     call the_config%init_config%write(output_unit)
     call the_config%init_g_config%write(output_unit)
     call the_config%knobs_config%write(output_unit)
     call the_config%kt_grids_config%write(output_unit)
     call the_config%kt_grids_box_config%write(output_unit)
     call the_config%kt_grids_range_config%write(output_unit)
     call the_config%kt_grids_single_config%write(output_unit)
     call the_config%kt_grids_specified_config%write(output_unit)
     if (allocated(the_config%kt_grids_specified_element_config)) then
        do i = 1, size(the_config%kt_grids_specified_element_config)
           call the_config%kt_grids_specified_element_config(i)%write(output_unit)
        end do
     else
        allocate(the_config%kt_grids_specified_element_config(1))
        call the_config%kt_grids_specified_element_config(1)%write(output_unit)
        deallocate(the_config%kt_grids_specified_element_config)
     end if
     call the_config%layouts_config%write(output_unit)
     call the_config%le_grids_config%write(output_unit)
#ifdef NEW_DIAG
     ! Note if we have new diagnostics then we currently write out both
     ! versions of the namelist -- we probably/may want to check the choice
     ! of `use_old_diagnostics` and then write the appropriate one based on that
     ! ideally there would only be one diagnostics config shared between the two
     ! approaches.
     call the_config%new_gs2_diagnostics_config%write(output_unit)
#endif
     call the_config%nonlinear_terms_config%write(output_unit)
     call the_config%normalisations_config%write(output_unit)
     call the_config%optimisation_config%write(output_unit)
     call the_config%parameters_config%write(output_unit)
     call the_config%parameter_scan_config%write(output_unit)
     call the_config%reinit_config%write(output_unit)
     call the_config%source_config%write(output_unit)
     call the_config%species_config%write(output_unit)
     if (allocated(the_config%species_element_config)) then
        do i = 1, size(the_config%species_element_config)
           call the_config%species_element_config(i)%write(output_unit)
        end do
     else
        allocate(the_config%species_element_config(1))
        call the_config%species_element_config(1)%write(output_unit)
        deallocate(the_config%species_element_config)
     end if
     if (allocated(the_config%stir_config)) then
        do i = 1, size(the_config%stir_config)
           call the_config%stir_config(i)%write(output_unit)
        end do
     else
        allocate(the_config%stir_config(1))
        call the_config%stir_config(1)%write(output_unit)
        deallocate(the_config%stir_config)
     end if
     call the_config%theta_grid_config%write(output_unit)
     call the_config%theta_grid_eik_config%write(output_unit)
     call the_config%theta_grid_file_config%write(output_unit)
     call the_config%theta_grid_gridgen_config%write(output_unit)
     call the_config%theta_grid_parameters_config%write(output_unit)
     call the_config%theta_grid_salpha_config%write(output_unit)
     call close_output_file(output_unit)
  end if

  call finish_mp
end program write_default_input_file
