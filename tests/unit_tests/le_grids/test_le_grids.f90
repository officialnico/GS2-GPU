!> A program that runs unit tests on the le_grids module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under the MIT licence
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_le_grids
  use unit_tests
  use le_grids, only: energy, integrate_species
  use mp, only: init_mp, finish_mp, mp_comm, broadcast, proc0
  use file_utils, only: init_file_utils, input_unit_exist
  use species, only: init_species, nspec
  use constants, only: pi
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  use gs2_main, only: gs2_program_state_type
  use gs2_main, only: initialize_gs2, finalize_gs2, initialize_wall_clock_timer
  use gs2_init, only: init, init_level_list
  implicit none
  type(gs2_program_state_type) :: state
  real :: eps
  integer, dimension(:), allocatable :: sizes
  complex, dimension (:,:,:), allocatable :: integrate_species_results
  complex, dimension (:,:,:), allocatable :: g
  real, dimension(:), allocatable :: species_weights
  real :: ene, tolerance
  integer :: iglo, ie, in_file
  logical :: exist
  namelist /test_le_grids_parameters/ tolerance

  ! Set up depenencies
  call init_mp
  state%mp_comm = mp_comm
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call init(state%init, init_level_list%dist_fn_layouts)

  tolerance = 1.0e-16
  if(proc0) then
     in_file = input_unit_exist('test_le_grids_parameters', exist)
     if (exist) read(in_file, nml = test_le_grids_parameters)
  end if
  call broadcast(tolerance)
  eps = tolerance
  call announce_module_test('le_grids')
  call announce_test('init_le_grids')

  allocate(sizes(3))
  sizes(1) = naky
  sizes(2) = ntheta0
  sizes(3) = ntgrid

  allocate(integrate_species_results(-ntgrid:ntgrid,naky,ntheta0))
  allocate(species_weights(nspec))
  allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))

  species_weights = 0.0
  species_weights(1) = 1.0
  g = cmplx(1.0, 1.0)
  integrate_species_results = cmplx(1.0,1.0)
  call announce_test('integrate species ions only g = 1')
  call process_test(&
    test_integrate_species(&
      g, &
      species_weights, &
      sizes, &
      integrate_species_results, &
      eps),&
    'integrate species ions only g = 1')

  ! Integrate a polynomial
  call announce_test(&
    'integrate species ions only g = energy^4.0 * 3.0*energy^1.25')
  g = cmplx(0.0, 0.0)
  do iglo = g_lo%llim_proc,g_lo%ulim_proc
    ie = ie_idx(g_lo, iglo)
    ene = energy(ie,1)
    g(:,:,iglo) = cmplx(ene**4.0 + 3.0*ene**1.25,1.0)
  end do

  integrate_species_results = cmplx(64.5070177949108,1.0)
  call process_test(&
    test_integrate_species(&
      g, &
      species_weights, &
      sizes, &
      integrate_species_results, &
      eps),&
    'integrate species ions only g = energy^4.0 * 3.0*energy^1.25')

  deallocate(integrate_species_results)
  deallocate(species_weights)

  call close_module_test('le_grids')
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)
  call finish_mp

  contains
    function test_integrate_species(g, weights, sizes, result, err)
      use unit_tests, only: announce_check, process_check, agrees_with
      use theta_grid, only: ntgrid
      use kt_grids, only: naky, ntheta0
      use gs2_layouts, only: g_lo

      complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (inout) :: g
      real, dimension (:), intent (in out) :: weights
      integer, dimension (:), intent (in) :: sizes
      complex, dimension (:,:,:), allocatable :: total
      complex, dimension (-ntgrid:,:,:), intent(in) :: result
      real, intent(in) :: err
      logical :: test_integrate_species
      logical :: test_result   !! Overall test result, true if all pass
      logical :: check_result  !! Result of individual tests

      test_result = .true.

      call announce_check('Size of naky')
      check_result = agrees_with(naky, sizes(1))
      call process_check(test_result, check_result, 'Size of naky')
      call announce_check('Size of ntheta0')
      check_result = agrees_with(ntheta0, sizes(2))
      call process_check(test_result, check_result, 'Size of ntheta0')
      call announce_check('Size of ntgrid')
      check_result = agrees_with(ntgrid, sizes(3))
      call process_check(test_result, check_result, 'Size of ntgrid')

      allocate(total(-ntgrid:ntgrid,ntheta0,naky))

      total = cmplx(0.0,0.0)

      call integrate_species(g, weights, total)
      print*,'Got',total(:,1,1),'Expected',result(:,1,1)
      print*,'Error:',abs(result(:,1,1) - total(:,1,1))
      call announce_check('total')
      check_result = agrees_with(total(:,1,1), result(:,1,1), err)
      call process_check(test_result, check_result, 'total ')

      deallocate(total)

      test_integrate_species = test_result

    end function test_integrate_species
end program test_le_grids
