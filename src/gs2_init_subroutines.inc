subroutine gs2_layouts_subroutine(going_up)
  use unit_tests, only: debug_message
  use gs2_layouts, only: init_gs2_layouts, finish_gs2_layouts
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- gs2_layouts (2)')
    call init_gs2_layouts
    current%level = init_level_list%gs2_layouts
  else
    call debug_message(1, 'gs2_init:init down reached init level -- gs2_layouts (2)')
    call finish_gs2_layouts
    current%level = init_level_list%gs2_layouts - 1
  end if
end subroutine gs2_layouts_subroutine

subroutine normalisations_subroutine(going_up)
  use unit_tests, only: debug_message
  use normalisations, only: init_normalisations, finish_normalisations
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- normalisations (3)')
    call init_normalisations
    current%level = init_level_list%normalisations
  else
    call debug_message(1, 'gs2_init:init down reached init level -- normalisations (3)')
    call finish_normalisations
    current%level = init_level_list%normalisations - 1
  end if
end subroutine normalisations_subroutine

subroutine theta_grid_params_subroutine(going_up)
  use unit_tests, only: debug_message
  use theta_grid_params, only: init_theta_grid_params, finish_theta_grid_params
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- theta_grid_params (4)')
    call init_theta_grid_params
    current%level = init_level_list%theta_grid_params
  else
    call debug_message(1, 'gs2_init:init down reached init level -- theta_grid_params (4)')
    call finish_theta_grid_params
    current%level = init_level_list%theta_grid_params - 1
  end if
end subroutine theta_grid_params_subroutine

subroutine fields_parameters_subroutine(going_up)
  use unit_tests, only: debug_message
  use fields, only: init_fields_parameters, finish_fields_parameters
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- fields_parameters (5)')
    call init_fields_parameters
    current%level = init_level_list%fields_parameters
  else
    call debug_message(1, 'gs2_init:init down reached init level -- fields_parameters (5)')
    call finish_fields_parameters
    current%level = init_level_list%fields_parameters - 1
  end if
end subroutine fields_parameters_subroutine

subroutine gs2_save_subroutine(going_up)
  use unit_tests, only: debug_message
  use gs2_save, only: init_gs2_save, finish_gs2_save
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- gs2_save (6)')
    call init_gs2_save
    current%level = init_level_list%gs2_save
  else
    call debug_message(1, 'gs2_init:init down reached init level -- gs2_save (6)')
    call finish_gs2_save
    current%level = init_level_list%gs2_save - 1
  end if
end subroutine gs2_save_subroutine

subroutine init_g_subroutine(going_up)
  use unit_tests, only: debug_message
  use init_g, only: init_init_g, finish_init_g
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- init_g (7)')
    call init_init_g
    current%level = init_level_list%init_g
  else
    call debug_message(1, 'gs2_init:init down reached init level -- init_g (7)')
    call finish_init_g
    current%level = init_level_list%init_g - 1
  end if
end subroutine init_g_subroutine

subroutine override_optimisations_subroutine(going_up)
  use unit_tests, only: debug_message
  use gs2_layouts, only: lso=>set_overrides
  use fields, only: fso=>set_overrides
  use dist_fn, only: dso=>set_overrides
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_optimisations (8)')
    if (current%opt_ov%is_initialised()) call lso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call fso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call dso(current%opt_ov)
    current%level = init_level_list%override_optimisations
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_optimisations (8)')

    current%level = init_level_list%override_optimisations - 1
  end if
end subroutine override_optimisations_subroutine

subroutine override_miller_geometry_subroutine(going_up)
  use unit_tests, only: debug_message
  use theta_grid_params, only: tgpso => set_overrides
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_miller_geometry (9)')
    if (current%mgeo_ov%is_initialised()) call tgpso(current%mgeo_ov)
    current%level = init_level_list%override_miller_geometry
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_miller_geometry (9)')

    current%level = init_level_list%override_miller_geometry - 1
  end if
end subroutine override_miller_geometry_subroutine

subroutine theta_grid_subroutine(going_up)
  use unit_tests, only: debug_message
  use theta_grid, only: init_theta_grid, finish_theta_grid
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- theta_grid (10)')
    call init_theta_grid
    current%level = init_level_list%theta_grid
  else
    call debug_message(1, 'gs2_init:init down reached init level -- theta_grid (10)')
    call finish_theta_grid
    current%level = init_level_list%theta_grid - 1
  end if
end subroutine theta_grid_subroutine

subroutine kt_grids_parameters_subroutine(going_up)
  use unit_tests, only: debug_message
  use kt_grids, only: init_kt_grids_parameters, finish_kt_grids_parameters
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- kt_grids_parameters (11)')
    call init_kt_grids_parameters
    current%level = init_level_list%kt_grids_parameters
  else
    call debug_message(1, 'gs2_init:init down reached init level -- kt_grids_parameters (11)')
    call finish_kt_grids_parameters
    current%level = init_level_list%kt_grids_parameters - 1
  end if
end subroutine kt_grids_parameters_subroutine

subroutine override_kt_grids_subroutine(going_up)
  use unit_tests, only: debug_message
  use kt_grids, only: ktso => set_overrides
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_kt_grids (12)')
    if (current%kt_ov%is_initialised()) call ktso(current%kt_ov)
    current%level = init_level_list%override_kt_grids
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_kt_grids (12)')

    current%level = init_level_list%override_kt_grids - 1
  end if
end subroutine override_kt_grids_subroutine

subroutine kt_grids_subroutine(going_up)
  use unit_tests, only: debug_message
  use kt_grids, only: init_kt_grids, finish_kt_grids
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- kt_grids (13)')
    call init_kt_grids
    current%level = init_level_list%kt_grids
  else
    call debug_message(1, 'gs2_init:init down reached init level -- kt_grids (13)')
    call finish_kt_grids
    current%level = init_level_list%kt_grids - 1
  end if
end subroutine kt_grids_subroutine

subroutine run_parameters_subroutine(going_up)
  use unit_tests, only: debug_message
  use run_parameters, only: init_run_parameters, finish_run_parameters
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- run_parameters (14)')
    call init_run_parameters
    current%level = init_level_list%run_parameters
  else
    call debug_message(1, 'gs2_init:init down reached init level -- run_parameters (14)')
    call finish_run_parameters
    current%level = init_level_list%run_parameters - 1
  end if
end subroutine run_parameters_subroutine

subroutine hyper_subroutine(going_up)
  use unit_tests, only: debug_message
  use hyper, only: init_hyper, finish_hyper
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- hyper (15)')
    call init_hyper
    current%level = init_level_list%hyper
  else
    call debug_message(1, 'gs2_init:init down reached init level -- hyper (15)')
    call finish_hyper
    current%level = init_level_list%hyper - 1
  end if
end subroutine hyper_subroutine

subroutine species_subroutine(going_up)
  use unit_tests, only: debug_message
  use species, only: init_species, finish_species
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- species (16)')
    call init_species
    current%level = init_level_list%species
  else
    call debug_message(1, 'gs2_init:init down reached init level -- species (16)')
    call finish_species
    current%level = init_level_list%species - 1
  end if
end subroutine species_subroutine

subroutine override_profiles_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: dfso=>set_overrides
  use species, only: sso=>set_overrides
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_profiles (17)')
    if (current%prof_ov%is_initialised()) call dfso(current%prof_ov)
    if (current%prof_ov%is_initialised()) call sso(current%prof_ov)
    current%level = init_level_list%override_profiles
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_profiles (17)')

    current%level = init_level_list%override_profiles - 1
  end if
end subroutine override_profiles_subroutine

subroutine override_timestep_subroutine(going_up)
  use unit_tests, only: debug_message
  use run_parameters, only: rso => set_overrides
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_timestep (18)')
    if (current%tstep_ov%is_initialised()) call rso(current%tstep_ov)
    current%level = init_level_list%override_timestep
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_timestep (18)')

    current%level = init_level_list%override_timestep - 1
  end if
end subroutine override_timestep_subroutine

subroutine le_grids_subroutine(going_up)
  use unit_tests, only: debug_message
  use le_grids, only: init_le_grids, finish_le_grids
  logical :: dummy1, dummy2
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- le_grids (19)')
    call init_le_grids(dummy1, dummy2)
    current%level = init_level_list%le_grids
  else
    call debug_message(1, 'gs2_init:init down reached init level -- le_grids (19)')
    call finish_le_grids
    current%level = init_level_list%le_grids - 1
  end if
end subroutine le_grids_subroutine

subroutine antenna_subroutine(going_up)
  use unit_tests, only: debug_message
  use antenna, only: init_antenna, finish_antenna
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- antenna (20)')
    call init_antenna
    current%level = init_level_list%antenna
  else
    call debug_message(1, 'gs2_init:init down reached init level -- antenna (20)')
    call finish_antenna
    current%level = init_level_list%antenna - 1
  end if
end subroutine antenna_subroutine

subroutine dist_fn_parameters_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_parameters, finish_dist_fn_parameters
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_parameters (21)')
    call init_dist_fn_parameters
    current%level = init_level_list%dist_fn_parameters
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_parameters (21)')
    call finish_dist_fn_parameters
    current%level = init_level_list%dist_fn_parameters - 1
  end if
end subroutine dist_fn_parameters_subroutine

subroutine dist_fn_layouts_subroutine(going_up)
  use unit_tests, only: debug_message
  use gs2_layouts, only: init_dist_fn_layouts, finish_dist_fn_layouts
  use kt_grids, only: naky, ntheta0
  use le_grids, only: nlambda, negrid
  use species, only: nspec
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_layouts (22)')
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
    current%level = init_level_list%dist_fn_layouts
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_layouts (22)')
    call finish_dist_fn_layouts
    current%level = init_level_list%dist_fn_layouts - 1
  end if
end subroutine dist_fn_layouts_subroutine

subroutine fields_level_1_subroutine(going_up)
  use unit_tests, only: debug_message
  use fields, only: init_fields_level_1, finish_fields_level_1
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- fields_level_1 (23)')
    call init_fields_level_1
    current%level = init_level_list%fields_level_1
  else
    call debug_message(1, 'gs2_init:init down reached init level -- fields_level_1 (23)')
    call finish_fields_level_1
    current%level = init_level_list%fields_level_1 - 1
  end if
end subroutine fields_level_1_subroutine

subroutine nonlinear_terms_subroutine(going_up)
  use unit_tests, only: debug_message
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- nonlinear_terms (24)')
    call init_nonlinear_terms
    current%level = init_level_list%nonlinear_terms
  else
    call debug_message(1, 'gs2_init:init down reached init level -- nonlinear_terms (24)')
    call finish_nonlinear_terms
    current%level = init_level_list%nonlinear_terms - 1
  end if
end subroutine nonlinear_terms_subroutine

subroutine dist_fn_arrays_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_arrays, finish_dist_fn_arrays
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_arrays (25)')
    call init_dist_fn_arrays
    current%level = init_level_list%dist_fn_arrays
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_arrays (25)')
    call finish_dist_fn_arrays
    current%level = init_level_list%dist_fn_arrays - 1
  end if
end subroutine dist_fn_arrays_subroutine

subroutine dist_fn_level_1_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_1, finish_dist_fn_level_1
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_level_1 (26)')
    call init_dist_fn_level_1
    current%level = init_level_list%dist_fn_level_1
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_level_1 (26)')
    call finish_dist_fn_level_1
    current%level = init_level_list%dist_fn_level_1 - 1
  end if
end subroutine dist_fn_level_1_subroutine

subroutine dist_fn_level_2_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_2, finish_dist_fn_level_2
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_level_2 (27)')
    call init_dist_fn_level_2
    current%level = init_level_list%dist_fn_level_2
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_level_2 (27)')
    call finish_dist_fn_level_2
    current%level = init_level_list%dist_fn_level_2 - 1
  end if
end subroutine dist_fn_level_2_subroutine

subroutine dist_fn_level_3_subroutine(going_up)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_3, finish_dist_fn_level_3
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- dist_fn_level_3 (28)')
    call init_dist_fn_level_3
    current%level = init_level_list%dist_fn_level_3
  else
    call debug_message(1, 'gs2_init:init down reached init level -- dist_fn_level_3 (28)')
    call finish_dist_fn_level_3
    current%level = init_level_list%dist_fn_level_3 - 1
  end if
end subroutine dist_fn_level_3_subroutine

subroutine collisions_subroutine(going_up)
  use unit_tests, only: debug_message
  use collisions, only: init_collisions, finish_collisions
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- collisions (29)')
    call init_collisions
    current%level = init_level_list%collisions
  else
    call debug_message(1, 'gs2_init:init down reached init level -- collisions (29)')
    call finish_collisions
    current%level = init_level_list%collisions - 1
  end if
end subroutine collisions_subroutine

subroutine fields_level_2_subroutine(going_up)
  use unit_tests, only: debug_message
  use fields, only: init_fields_level_2, finish_fields_level_2
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- fields_level_2 (30)')
    call init_fields_level_2
    current%level = init_level_list%fields_level_2
  else
    call debug_message(1, 'gs2_init:init down reached init level -- fields_level_2 (30)')
    call finish_fields_level_2
    current%level = init_level_list%fields_level_2 - 1
  end if
end subroutine fields_level_2_subroutine

subroutine override_initial_values_subroutine(going_up)
  use unit_tests, only: debug_message

  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- override_initial_values (31)')

    current%level = init_level_list%override_initial_values
  else
    call debug_message(1, 'gs2_init:init down reached init level -- override_initial_values (31)')

    current%level = init_level_list%override_initial_values - 1
  end if
end subroutine override_initial_values_subroutine

subroutine set_initial_values_subroutine(going_up)
  use unit_tests, only: debug_message

  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- set_initial_values (32)')
    call set_initial_field_and_dist_fn_values(current)
    current%level = init_level_list%set_initial_values
  else
    call debug_message(1, 'gs2_init:init down reached init level -- set_initial_values (32)')

    current%level = init_level_list%set_initial_values - 1
  end if
end subroutine set_initial_values_subroutine

subroutine full_subroutine(going_up)
  use unit_tests, only: debug_message

  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- full (33)')

    current%level = init_level_list%full
  else
    call debug_message(1, 'gs2_init:init down reached init level -- full (33)')

    current%level = init_level_list%full - 1
  end if
end subroutine full_subroutine