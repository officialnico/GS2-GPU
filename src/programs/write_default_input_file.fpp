program write_default_input_file
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: get_unused_unit, open_output_file, close_output_file
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
#ifdef NEW_DIAG
  use diagnostics_config, only: gs2_diagnostics_config_type => new_gs2_diagnostics_config_type
#else
  use gs2_diagnostics, only: gs2_diagnostics_config_type
#endif
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

  type(ballstab_config_type) :: ballstab_config
  type(collisions_config_type) :: collisions_config
  type(dist_fn_config_type) :: dist_fn_config
  type(dist_fn_species_config_type) :: dist_fn_species_config
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
  type(kt_grids_specified_element_config_type) :: kt_grids_specified_element_config
  type(layouts_config_type) :: layouts_config
  type(le_grids_config_type) :: le_grids_config
  type(nonlinear_terms_config_type) :: nonlinear_terms_config
  type(normalisations_config_type) :: normalisations_config
  type(optimisation_config_type) :: optimisation_config
  type(parameters_config_type) :: parameters_config
  type(parameter_scan_config_type) :: parameter_scan_config
  type(reinit_config_type) :: reinit_config
  type(source_config_type) :: source_config
  type(species_config_type) :: species_config
  type(species_element_config_type) :: species_element_config
  type(stir_config_type) :: stir_config
  type(theta_grid_config_type) :: theta_grid_config
  type(theta_grid_eik_config_type) :: theta_grid_eik_config
  type(theta_grid_file_config_type) :: theta_grid_file_config
  type(theta_grid_gridgen_config_type) :: theta_grid_gridgen_config
  type(theta_grid_parameters_config_type) :: theta_grid_parameters_config
  type(theta_grid_salpha_config_type) :: theta_grid_salpha_config

  logical :: list
  integer :: output_unit

  call init_mp

  if (proc0) then

     ballstab_config = ballstab_config_type()
     collisions_config = collisions_config_type()
     dist_fn_config = dist_fn_config_type()
     dist_fn_species_config = dist_fn_species_config_type()
     driver_config = driver_config_type()
#ifdef WITH_EIG
     eigval_config = eigval_config_type()
#endif
     fields_config = fields_config_type()
     gs2_diagnostics_config   = gs2_diagnostics_config_type()
     hyper_config = hyper_config_type()
     ingen_config = ingen_config_type()
     init_config = init_config_type()
     init_g_config = init_g_config_type()
     knobs_config = knobs_config_type()
     kt_grids_config = kt_grids_config_type()
     kt_grids_box_config = kt_grids_box_config_type()
     kt_grids_range_config = kt_grids_range_config_type()
     kt_grids_single_config = kt_grids_single_config_type()
     kt_grids_specified_config = kt_grids_specified_config_type()
     kt_grids_specified_element_config = kt_grids_specified_element_config_type()
     layouts_config = layouts_config_type()
     le_grids_config = le_grids_config_type()
     nonlinear_terms_config = nonlinear_terms_config_type()
     normalisations_config = normalisations_config_type()
     optimisation_config = optimisation_config_type()
     parameters_config = parameters_config_type()
     parameter_scan_config = parameter_scan_config_type()
     reinit_config = reinit_config_type()
     source_config = source_config_type()
     species_config = species_config_type()
     species_element_config = species_element_config_type()
     stir_config = stir_config_type()
     theta_grid_config = theta_grid_config_type()
     theta_grid_eik_config = theta_grid_eik_config_type()
     theta_grid_file_config = theta_grid_file_config_type()
     theta_grid_gridgen_config = theta_grid_gridgen_config_type()
     theta_grid_parameters_config = theta_grid_parameters_config_type()
     theta_grid_salpha_config = theta_grid_salpha_config_type()

     call get_unused_unit(output_unit)
     open(unit = output_unit, action = "write", file = "defaults.in")
     call ballstab_config%write(output_unit)
     call collisions_config%write(output_unit)
     call dist_fn_config%write(output_unit)
     call dist_fn_species_config%write(output_unit)
     call driver_config%write(output_unit)
#ifdef WITH_EIG
     call eigval_config%write(output_unit)
#endif
     call fields_config%write(output_unit)
     call gs2_diagnostics_config%write(output_unit)
     call hyper_config%write(output_unit)
     call ingen_config%write(output_unit)
     call init_config%write(output_unit)
     call init_g_config%write(output_unit)
     call knobs_config%write(output_unit)
     call kt_grids_config%write(output_unit)
     call kt_grids_box_config%write(output_unit)
     call kt_grids_range_config%write(output_unit)
     call kt_grids_single_config%write(output_unit)
     call kt_grids_specified_config%write(output_unit)
     call kt_grids_specified_element_config%write(output_unit)
     call layouts_config%write(output_unit)
     call le_grids_config%write(output_unit)
     call nonlinear_terms_config%write(output_unit)
     call normalisations_config%write(output_unit)
     call optimisation_config%write(output_unit)
     call parameters_config%write(output_unit)
     call parameter_scan_config%write(output_unit)
     call reinit_config%write(output_unit)
     call source_config%write(output_unit)
     call species_config%write(output_unit)
     call species_element_config%write(output_unit)
     call stir_config%write(output_unit)
     call theta_grid_config%write(output_unit)
     call theta_grid_eik_config%write(output_unit)
     call theta_grid_file_config%write(output_unit)
     call theta_grid_gridgen_config%write(output_unit)
     call theta_grid_parameters_config%write(output_unit)
     call theta_grid_salpha_config%write(output_unit)
     call close_output_file(output_unit)
  end if

  call finish_mp
end program write_default_input_file
