allocate( init_levels(33))
      init_levels(1) = init_level_type(name = 'basic', level = 1)
      init_levels(2) = init_level_type(name = 'gs2_layouts', level = 2)
      init_levels(3) = init_level_type(name = 'normalisations', level = 3)
      init_levels(4) = init_level_type(name = 'theta_grid_params', level = 4)
      init_levels(5) = init_level_type(name = 'fields_parameters', level = 5)
      init_levels(6) = init_level_type(name = 'gs2_save', level = 6)
      init_levels(7) = init_level_type(name = 'init_g', level = 7)
      init_levels(8) = init_level_type(name = 'override_optimisations', level = 8)
      init_levels(9) = init_level_type(name = 'override_miller_geometry', level = 9)
      init_levels(10) = init_level_type(name = 'theta_grid', level = 10)
      init_levels(11) = init_level_type(name = 'kt_grids_parameters', level = 11)
      init_levels(12) = init_level_type(name = 'override_kt_grids', level = 12)
      init_levels(13) = init_level_type(name = 'kt_grids', level = 13)
      init_levels(14) = init_level_type(name = 'run_parameters', level = 14)
      init_levels(15) = init_level_type(name = 'hyper', level = 15)
      init_levels(16) = init_level_type(name = 'species', level = 16)
      init_levels(17) = init_level_type(name = 'override_profiles', level = 17)
      init_levels(18) = init_level_type(name = 'le_grids', level = 18)
      init_levels(19) = init_level_type(name = 'antenna', level = 19)
      init_levels(20) = init_level_type(name = 'dist_fn_parameters', level = 20)
      init_levels(21) = init_level_type(name = 'dist_fn_layouts', level = 21)
      init_levels(22) = init_level_type(name = 'fields_level_1', level = 22)
      init_levels(23) = init_level_type(name = 'nonlinear_terms', level = 23)
      init_levels(24) = init_level_type(name = 'dist_fn_arrays', level = 24)
      init_levels(25) = init_level_type(name = 'dist_fn_level_1', level = 25)
      init_levels(26) = init_level_type(name = 'dist_fn_level_2', level = 26)
      init_levels(27) = init_level_type(name = 'override_timestep', level = 27)
      init_levels(28) = init_level_type(name = 'dist_fn_level_3', level = 28)
      init_levels(29) = init_level_type(name = 'collisions', level = 29)
      init_levels(30) = init_level_type(name = 'fields_level_2', level = 30)
      init_levels(31) = init_level_type(name = 'override_initial_values', level = 31)
      init_levels(32) = init_level_type(name = 'set_initial_values', level = 32)
      init_levels(33) = init_level_type(name = 'full', level = 33)
      init_levels(1)%change_level_specific => basic_subroutine
      init_levels(2)%change_level_specific => gs2_layouts_subroutine
      init_levels(3)%change_level_specific => normalisations_subroutine
      init_levels(4)%change_level_specific => theta_grid_params_subroutine
      init_levels(5)%change_level_specific => fields_parameters_subroutine
      init_levels(6)%change_level_specific => gs2_save_subroutine
      init_levels(7)%change_level_specific => init_g_subroutine
      init_levels(8)%change_level_specific => override_optimisations_subroutine
      init_levels(9)%change_level_specific => override_miller_geometry_subroutine
      init_levels(10)%change_level_specific => theta_grid_subroutine
      init_levels(11)%change_level_specific => kt_grids_parameters_subroutine
      init_levels(12)%change_level_specific => override_kt_grids_subroutine
      init_levels(13)%change_level_specific => kt_grids_subroutine
      init_levels(14)%change_level_specific => run_parameters_subroutine
      init_levels(15)%change_level_specific => hyper_subroutine
      init_levels(16)%change_level_specific => species_subroutine
      init_levels(17)%change_level_specific => override_profiles_subroutine
      init_levels(18)%change_level_specific => le_grids_subroutine
      init_levels(19)%change_level_specific => antenna_subroutine
      init_levels(20)%change_level_specific => dist_fn_parameters_subroutine
      init_levels(21)%change_level_specific => dist_fn_layouts_subroutine
      init_levels(22)%change_level_specific => fields_level_1_subroutine
      init_levels(23)%change_level_specific => nonlinear_terms_subroutine
      init_levels(24)%change_level_specific => dist_fn_arrays_subroutine
      init_levels(25)%change_level_specific => dist_fn_level_1_subroutine
      init_levels(26)%change_level_specific => dist_fn_level_2_subroutine
      init_levels(27)%change_level_specific => override_timestep_subroutine
      init_levels(28)%change_level_specific => dist_fn_level_3_subroutine
      init_levels(29)%change_level_specific => collisions_subroutine
      init_levels(30)%change_level_specific => fields_level_2_subroutine
      init_levels(31)%change_level_specific => override_initial_values_subroutine
      init_levels(32)%change_level_specific => set_initial_values_subroutine
      init_levels(33)%change_level_specific => full_subroutine