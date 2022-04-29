integer function basic_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message

  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function basic_function

integer function gs2_layouts_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use gs2_layouts, only: init_gs2_layouts, finish_gs2_layouts
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_gs2_layouts
    new_level = self%level
  else
    call finish_gs2_layouts
    new_level = self%level - 1
  end if
end function gs2_layouts_function

integer function normalisations_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use normalisations, only: init_normalisations, finish_normalisations
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_normalisations
    new_level = self%level
  else
    call finish_normalisations
    new_level = self%level - 1
  end if
end function normalisations_function

integer function theta_grid_params_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use theta_grid_params, only: init_theta_grid_params, finish_theta_grid_params
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_theta_grid_params
    new_level = self%level
  else
    call finish_theta_grid_params
    new_level = self%level - 1
  end if
end function theta_grid_params_function

integer function fields_parameters_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use fields, only: init_fields_parameters, finish_fields_parameters
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_parameters
    new_level = self%level
  else
    call finish_fields_parameters
    new_level = self%level - 1
  end if
end function fields_parameters_function

integer function gs2_save_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use gs2_save, only: init_gs2_save, finish_gs2_save
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_gs2_save
    new_level = self%level
  else
    call finish_gs2_save
    new_level = self%level - 1
  end if
end function gs2_save_function

integer function init_g_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use init_g, only: init_init_g, finish_init_g
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_init_g
    new_level = self%level
  else
    call finish_init_g
    new_level = self%level - 1
  end if
end function init_g_function

integer function override_optimisations_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use gs2_layouts, only: lso=>set_overrides
  use fields, only: fso=>set_overrides
  use dist_fn, only: dso=>set_overrides
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%opt_ov%is_initialised()) call lso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call fso(current%opt_ov)
    if (current%opt_ov%is_initialised()) call dso(current%opt_ov)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_optimisations_function

integer function override_miller_geometry_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use theta_grid_params, only: tgpso => set_overrides
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%mgeo_ov%is_initialised()) call tgpso(current%mgeo_ov)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_miller_geometry_function

integer function theta_grid_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use theta_grid, only: init_theta_grid, finish_theta_grid
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_theta_grid
    new_level = self%level
  else
    call finish_theta_grid
    new_level = self%level - 1
  end if
end function theta_grid_function

integer function kt_grids_parameters_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use kt_grids, only: init_kt_grids_parameters, finish_kt_grids_parameters
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_kt_grids_parameters
    new_level = self%level
  else
    call finish_kt_grids_parameters
    new_level = self%level - 1
  end if
end function kt_grids_parameters_function

integer function override_kt_grids_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use kt_grids, only: ktso => set_overrides
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%kt_ov%is_initialised()) call ktso(current%kt_ov)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_kt_grids_function

integer function kt_grids_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use kt_grids, only: init_kt_grids, finish_kt_grids
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_kt_grids
    new_level = self%level
  else
    call finish_kt_grids
    new_level = self%level - 1
  end if
end function kt_grids_function

integer function run_parameters_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use run_parameters, only: init_run_parameters, finish_run_parameters
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_run_parameters
    new_level = self%level
  else
    call finish_run_parameters
    new_level = self%level - 1
  end if
end function run_parameters_function

integer function hyper_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use hyper, only: init_hyper, finish_hyper
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_hyper
    new_level = self%level
  else
    call finish_hyper
    new_level = self%level - 1
  end if
end function hyper_function

integer function species_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use species, only: init_species, finish_species
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_species
    new_level = self%level
  else
    call finish_species
    new_level = self%level - 1
  end if
end function species_function

integer function override_profiles_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: dfso=>set_overrides
  use species, only: sso=>set_overrides
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%prof_ov%is_initialised()) call dfso(current%prof_ov)
    if (current%prof_ov%is_initialised()) call sso(current%prof_ov)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_profiles_function

integer function le_grids_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use le_grids, only: init_le_grids, finish_le_grids
  logical :: dummy1, dummy2
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_le_grids(dummy1, dummy2)
    new_level = self%level
  else
    call finish_le_grids
    new_level = self%level - 1
  end if
end function le_grids_function

integer function antenna_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use antenna, only: init_antenna, finish_antenna
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_antenna
    new_level = self%level
  else
    call finish_antenna
    new_level = self%level - 1
  end if
end function antenna_function

integer function dist_fn_parameters_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_parameters, finish_dist_fn_parameters
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_parameters
    new_level = self%level
  else
    call finish_dist_fn_parameters
    new_level = self%level - 1
  end if
end function dist_fn_parameters_function

integer function dist_fn_layouts_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use gs2_layouts, only: init_dist_fn_layouts, finish_dist_fn_layouts
  use kt_grids, only: naky, ntheta0
  use le_grids, only: nlambda, negrid
  use species, only: nspec
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
    new_level = self%level
  else
    call finish_dist_fn_layouts
    new_level = self%level - 1
  end if
end function dist_fn_layouts_function

integer function fields_level_1_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use fields, only: init_fields_level_1, finish_fields_level_1
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_level_1
    new_level = self%level
  else
    call finish_fields_level_1
    new_level = self%level - 1
  end if
end function fields_level_1_function

integer function nonlinear_terms_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_nonlinear_terms
    new_level = self%level
  else
    call finish_nonlinear_terms
    new_level = self%level - 1
  end if
end function nonlinear_terms_function

integer function dist_fn_arrays_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_arrays, finish_dist_fn_arrays
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_arrays
    new_level = self%level
  else
    call finish_dist_fn_arrays
    new_level = self%level - 1
  end if
end function dist_fn_arrays_function

integer function dist_fn_level_1_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_1, finish_dist_fn_level_1
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_1
    new_level = self%level
  else
    call finish_dist_fn_level_1
    new_level = self%level - 1
  end if
end function dist_fn_level_1_function

integer function dist_fn_level_2_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_2, finish_dist_fn_level_2
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_2
    new_level = self%level
  else
    call finish_dist_fn_level_2
    new_level = self%level - 1
  end if
end function dist_fn_level_2_function

integer function override_timestep_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use run_parameters, only: rso => set_overrides
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    if (current%tstep_ov%is_initialised()) call rso(current%tstep_ov)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_timestep_function

integer function dist_fn_level_3_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use dist_fn, only: init_dist_fn_level_3, finish_dist_fn_level_3
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_dist_fn_level_3
    new_level = self%level
  else
    call finish_dist_fn_level_3
    new_level = self%level - 1
  end if
end function dist_fn_level_3_function

integer function collisions_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use collisions, only: init_collisions, finish_collisions
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_collisions
    new_level = self%level
  else
    call finish_collisions
    new_level = self%level - 1
  end if
end function collisions_function

integer function fields_level_2_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message
  use fields, only: init_fields_level_2, finish_fields_level_2
  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call init_fields_level_2
    new_level = self%level
  else
    call finish_fields_level_2
    new_level = self%level - 1
  end if
end function fields_level_2_function

integer function override_initial_values_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message

  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function override_initial_values_function

integer function set_initial_values_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message

  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then
    call set_initial_field_and_dist_fn_values(current)
    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function set_initial_values_function

integer function full_function(self, current, going_up) result(new_level)
  use unit_tests, only: debug_message

  class(init_level_type), intent(in out) :: self
  type(init_type), intent(in) :: current
  logical, intent(in) :: going_up
  if (going_up) then

    new_level = self%level
  else

    new_level = self%level - 1
  end if
end function full_function