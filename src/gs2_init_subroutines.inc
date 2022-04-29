subroutine basic_subroutine(self, current, going_up)

  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

  else

  end if
end subroutine basic_subroutine

subroutine gs2_layouts_subroutine(self, current, going_up)
  use gs2_layouts, only: init_gs2_layouts, finish_gs2_layouts
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_gs2_layouts
  else
    call finish_gs2_layouts
  end if
end subroutine gs2_layouts_subroutine

subroutine normalisations_subroutine(self, current, going_up)
  use normalisations, only: init_normalisations, finish_normalisations
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_normalisations
  else
    call finish_normalisations
  end if
end subroutine normalisations_subroutine

subroutine theta_grid_params_subroutine(self, current, going_up)
  use theta_grid_params, only: init_theta_grid_params, finish_theta_grid_params
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_theta_grid_params
  else
    call finish_theta_grid_params
  end if
end subroutine theta_grid_params_subroutine

subroutine fields_parameters_subroutine(self, current, going_up)
  use fields, only: init_fields_parameters, finish_fields_parameters
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_parameters
  else
    call finish_fields_parameters
  end if
end subroutine fields_parameters_subroutine

subroutine gs2_save_subroutine(self, current, going_up)
  use gs2_save, only: init_gs2_save, finish_gs2_save
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_gs2_save
  else
    call finish_gs2_save
  end if
end subroutine gs2_save_subroutine

subroutine init_g_subroutine(self, current, going_up)
  use init_g, only: init_init_g, finish_init_g
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_init_g
  else
    call finish_init_g
  end if
end subroutine init_g_subroutine

subroutine override_optimisations_subroutine(self, current, going_up)
  use gs2_layouts, only: lso=>set_overrides
  use fields, only: fso=>set_overrides
  use dist_fn, only: dso=>set_overrides
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%opt_ov%is_initialised()) call lso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call fso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call dso(current%opt_ov)
  else

  end if
end subroutine override_optimisations_subroutine

subroutine override_miller_geometry_subroutine(self, current, going_up)
  use theta_grid_params, only: tgpso => set_overrides
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%mgeo_ov%is_initialised()) call tgpso(current%mgeo_ov)
  else

  end if
end subroutine override_miller_geometry_subroutine

subroutine theta_grid_subroutine(self, current, going_up)
  use theta_grid, only: init_theta_grid, finish_theta_grid
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_theta_grid
  else
    call finish_theta_grid
  end if
end subroutine theta_grid_subroutine

subroutine kt_grids_parameters_subroutine(self, current, going_up)
  use kt_grids, only: init_kt_grids_parameters, finish_kt_grids_parameters
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_kt_grids_parameters
  else
    call finish_kt_grids_parameters
  end if
end subroutine kt_grids_parameters_subroutine

subroutine override_kt_grids_subroutine(self, current, going_up)
  use kt_grids, only: ktso => set_overrides
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%kt_ov%is_initialised()) call ktso(current%kt_ov)
  else

  end if
end subroutine override_kt_grids_subroutine

subroutine kt_grids_subroutine(self, current, going_up)
  use kt_grids, only: init_kt_grids, finish_kt_grids
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_kt_grids
  else
    call finish_kt_grids
  end if
end subroutine kt_grids_subroutine

subroutine run_parameters_subroutine(self, current, going_up)
  use run_parameters, only: init_run_parameters, finish_run_parameters
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_run_parameters
  else
    call finish_run_parameters
  end if
end subroutine run_parameters_subroutine

subroutine hyper_subroutine(self, current, going_up)
  use hyper, only: init_hyper, finish_hyper
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_hyper
  else
    call finish_hyper
  end if
end subroutine hyper_subroutine

subroutine species_subroutine(self, current, going_up)
  use species, only: init_species, finish_species
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_species
  else
    call finish_species
  end if
end subroutine species_subroutine

subroutine override_profiles_subroutine(self, current, going_up)
  use dist_fn, only: dfso=>set_overrides
  use species, only: sso=>set_overrides
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%prof_ov%is_initialised()) call dfso(current%prof_ov)
    if (current%prof_ov%is_initialised()) call sso(current%prof_ov)
  else

  end if
end subroutine override_profiles_subroutine

subroutine le_grids_subroutine(self, current, going_up)
  use le_grids, only: init_le_grids, finish_le_grids
  logical :: dummy1, dummy2
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_le_grids(dummy1, dummy2)
  else
    call finish_le_grids
  end if
end subroutine le_grids_subroutine

subroutine antenna_subroutine(self, current, going_up)
  use antenna, only: init_antenna, finish_antenna
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_antenna
  else
    call finish_antenna
  end if
end subroutine antenna_subroutine

subroutine dist_fn_parameters_subroutine(self, current, going_up)
  use dist_fn, only: init_dist_fn_parameters, finish_dist_fn_parameters
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_parameters
  else
    call finish_dist_fn_parameters
  end if
end subroutine dist_fn_parameters_subroutine

subroutine dist_fn_layouts_subroutine(self, current, going_up)
  use gs2_layouts, only: init_dist_fn_layouts, finish_dist_fn_layouts
  use kt_grids, only: naky, ntheta0
  use le_grids, only: nlambda, negrid
  use species, only: nspec
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
  else
    call finish_dist_fn_layouts
  end if
end subroutine dist_fn_layouts_subroutine

subroutine fields_level_1_subroutine(self, current, going_up)
  use fields, only: init_fields_level_1, finish_fields_level_1
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_level_1
  else
    call finish_fields_level_1
  end if
end subroutine fields_level_1_subroutine

subroutine nonlinear_terms_subroutine(self, current, going_up)
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_nonlinear_terms
  else
    call finish_nonlinear_terms
  end if
end subroutine nonlinear_terms_subroutine

subroutine dist_fn_arrays_subroutine(self, current, going_up)
  use dist_fn, only: init_dist_fn_arrays, finish_dist_fn_arrays
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_arrays
  else
    call finish_dist_fn_arrays
  end if
end subroutine dist_fn_arrays_subroutine

subroutine dist_fn_level_1_subroutine(self, current, going_up)
  use dist_fn, only: init_dist_fn_level_1, finish_dist_fn_level_1
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_1
  else
    call finish_dist_fn_level_1
  end if
end subroutine dist_fn_level_1_subroutine

subroutine dist_fn_level_2_subroutine(self, current, going_up)
  use dist_fn, only: init_dist_fn_level_2, finish_dist_fn_level_2
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_2
  else
    call finish_dist_fn_level_2
  end if
end subroutine dist_fn_level_2_subroutine

subroutine override_timestep_subroutine(self, current, going_up)
  use run_parameters, only: rso => set_overrides
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%tstep_ov%is_initialised()) call rso(current%tstep_ov)
  else

  end if
end subroutine override_timestep_subroutine

subroutine dist_fn_level_3_subroutine(self, current, going_up)
  use dist_fn, only: init_dist_fn_level_3, finish_dist_fn_level_3
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_3
  else
    call finish_dist_fn_level_3
  end if
end subroutine dist_fn_level_3_subroutine

subroutine collisions_subroutine(self, current, going_up)
  use collisions, only: init_collisions, finish_collisions
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_collisions
  else
    call finish_collisions
  end if
end subroutine collisions_subroutine

subroutine fields_level_2_subroutine(self, current, going_up)
  use fields, only: init_fields_level_2, finish_fields_level_2
  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_level_2
  else
    call finish_fields_level_2
  end if
end subroutine fields_level_2_subroutine

subroutine override_initial_values_subroutine(self, current, going_up)

  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

  else

  end if
end subroutine override_initial_values_subroutine

subroutine set_initial_values_subroutine(self, current, going_up)

  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call set_initial_field_and_dist_fn_values(current)
  else

  end if
end subroutine set_initial_values_subroutine

subroutine full_subroutine(self, current, going_up)

  class(init_level_type), intent(in) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

  else

  end if
end subroutine full_subroutine