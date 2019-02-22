!!!!!!ADD_TEST_SUITE(test_diagnostics_fluxes_mod_suite)
!!!!!!ADD_TEST_SUITE(test_gs2_layouts_mod_suite)
!!!!!!ADD_TEST_SUITE(test_gs2_time_mod_suite)
!!!!!!
!!!!!!! working, but need tidying
!!!!!!ADD_TEST_SUITE(test_gs2_save_mod_suite)
!!!!!!ADD_TEST_SUITE(test_fields_mod_suite)
!!!!!!ADD_TEST_SUITE(test_normalisations_mod_suite)
!!!!!!ADD_TEST_SUITE(test_theta_grid_mod_suite)
!!!!!!ADD_TEST_SUITE(test_run_parameters_mod_suite)
!!!!!!ADD_TEST_SUITE(test_species_mod_suite)
!!!!!!ADD_TEST_SUITE(test_init_g_mod_suite)
!!!!!!
!! working, but need reformatting
ADD_TEST_SUITE(test_antenna_data_mod_suite)
!!!!!!
!!!!!!! working, but fail
!!!!!!!ADD_TEST_SUITE(test_theta_grid_params_mod_suite)
!!!!!!
!!!!!!! broken
!!!!!!!ADD_TEST_SUITE(test_gs2_layouts_mod_suite)
!!!!!!!ADD_TEST_SUITE(test_dist_fn_mod_suite)
!!!!!!!ADD_TEST_SUITE(test_antenna_mod_suite) ! test_antenna works before dist_fn_mod_suite, but not after. TODO check clean-up in dist fn
!!!!!!!ADD_TEST_SUITE(test_dist_fn_arrays_mod_suite)
!!!!!!
!!!!!!!Requires initialization to fields_level_2 - seems too high?
!!!!!!!ADD_TEST_SUITE(test_optimisation_config_mod_suite)
!!!!!!
!!!!!!! untested
!!!!!!!!!---!!!
!!!!!!!!!---!!!! Non-functional with shared memory
!!!!!!!!!---!!!
!!!!!!!!!---!!!! Non-functional tests
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_overrides_mod_suite)
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_gs2_init_mod_suite)
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_file_utils_suite)
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_mp_suite)
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_constants_suite)
!!!!!!!!!---!!!!ADD_TEST_SUITE(test_gs2_main_mod_suite)
