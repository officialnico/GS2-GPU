!> A test program to check that applying the linked boundary condition results
!> in consistent boundaries (i.e. values are in agreement at the duplicate theta
!> grid points).
program test_linked_bc
  use unit_tests
  use gs2_main, only: gs2_program_state_type
  use gs2_main, only: initialize_gs2, finalize_gs2, initialize_wall_clock_timer
  use mp, only: init_mp, finish_mp, mp_comm, mp_abort, proc0
  use gs2_init, only: init, init_level_list
  use dist_fn_arrays, only: gnew
  use dist_fn, only: enforce_linked_boundary_conditions
  use dist_fn, only: l_links, r_links
  use fields_arrays, only: phinew, bparnew
  use gs2_layouts, only: g_lo, it_idx, ik_idx, il_idx
  use kt_grids, only: naky, ntheta0, theta0
  use theta_grid, only: ntgrid, theta
  use gs2_diagnostics, only: save_restart_dist_fn
  implicit none
  real :: maxerr
  real, parameter :: tolerance = epsilon(maxerr) * 100.0
  logical, parameter :: debug = .false.
  type(gs2_program_state_type) :: state
  integer :: iglo, it, ik, il
  complex, dimension(:, :, :), allocatable :: target_function
  logical :: passes

  call init_mp
  state%mp_comm = mp_comm
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call init(state%init, init_level_list%dist_fn_level_3)
  call announce_module_test('linked_boundary_conditions')

  phinew = 0. ; bparnew = 0.
  gnew = 0.

  ! Set the continuous function that we expect to recover
  ! For now just (theta-theta0)**3, but we should probably
  ! set different functions in the real and imaginary parts
  ! and we may want to test different functions for different
  ! sigma. For now we just flip the sign for isig = 2
  allocate(target_function, mold = phinew)
  do ik = 1, naky
     do it = 1, ntheta0
        target_function(:, it, ik) = (theta - theta0(it, ik))**3
     end do
  end do

  ! Initialise the distribution function
  do iglo = g_lo%llim_proc, g_lo%ulim_proc
     it = it_idx(g_lo, iglo)
     ik = ik_idx(g_lo, iglo)
     il = il_idx(g_lo, iglo)
     ! Set to target function but offset so incoming boundary = 0 for all cells
     gnew(:, 1, iglo) = target_function(:, it, ik)
     if (l_links(it, ik) > 0) gnew(:, 1, iglo) = gnew(:, 1, iglo) - target_function(-ntgrid, it, ik)
     gnew(:, 2, iglo) = -target_function(:, it, ik)
     if (r_links(it, ik) > 0) gnew(:, 2, iglo) = gnew(:, 2, iglo) + target_function(ntgrid, it, ik)
  end do

  maxerr = check_error(gnew)
  if(proc0 .and. debug) print*,'Initial error : ',maxerr,maxval(abs(gnew)),maxval(abs(target_function)),ntheta0
  if (maxerr < epsilon(0.0)) call mp_abort(&
       "Function error pre boundary conditions already very small.", .true.)
  if (debug) call dump('.pre')
  call enforce_linked_boundary_conditions
  if (debug) call dump('.post')

  maxerr = check_error(gnew)

  passes = maxerr < tolerance

  if (.not. passes) then
     if (proc0) print*,'Test fails with error : ',maxerr
  end if

  call process_test(passes, 'linked_boundary_conditions')
  call close_module_test('linked_boundary_conditions')
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)
  call finish_mp

contains
  real function check_error(g_in) result(err)
    use mp, only: max_allreduce, sum_allreduce
    use dist_fn, only: itright, itleft
    use gs2_layouts, only: ie_idx, is_idx, idx
    use le_grids, only: il_is_trapped
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_in
    complex, dimension(:, :), allocatable :: right_boundaries, left_boundaries
    integer :: iglo, it, ik, il, iglo_right, iglo_left, ie, is
    err = 0.0
    ! This approach of gathering and checking all boundary values is
    ! not appropriate for inclusion in the main code as it could be
    ! prohibitively expensive. As this is a test case we expect the
    ! problem size to be small so this approach should be acceptable.
    ! In the main code we could use pass_right and pass_left to
    ! communicate just with the relevant processors.
    allocate(right_boundaries(g_lo%llim_world:g_lo%ulim_world, 2))
    allocate(left_boundaries(g_lo%llim_world:g_lo%ulim_world, 2))
    right_boundaries = 0
    left_boundaries = 0
    right_boundaries(g_lo%llim_proc:g_lo%ulim_proc, 1) = g_in(ntgrid, 1, :)
    right_boundaries(g_lo%llim_proc:g_lo%ulim_proc, 2) = g_in(ntgrid, 2, :)
    left_boundaries(g_lo%llim_proc:g_lo%ulim_proc, 1) = g_in(-ntgrid, 1, :)
    left_boundaries(g_lo%llim_proc:g_lo%ulim_proc, 2) = g_in(-ntgrid, 2, :)
    call sum_allreduce(right_boundaries)
    call sum_allreduce(left_boundaries)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       if (il_is_trapped(il)) cycle

       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)

       if (r_links(it, ik) > 0) then
          iglo_right = idx(g_lo, ik, itright(it, ik), il, ie, is)
          err = max(err, &
               abs(right_boundaries(iglo, 1) - left_boundaries(iglo_right, 1)))
       end if

       if (l_links(it, ik) > 0) then
          iglo_left = idx(g_lo, ik, itleft(it, ik), il, ie, is)
          err = max(err, &
               abs(left_boundaries(iglo, 2) - right_boundaries(iglo_left, 2)))
       end if

    end do
    call max_allreduce(err)
  end function check_error

  subroutine dump(ext)
    implicit none
    character(len=*), intent(in) :: ext
    call save_restart_dist_fn(.true., .false., .false., .false., 0.0, fileopt_base = ext)
  end subroutine dump

end program test_linked_bc
