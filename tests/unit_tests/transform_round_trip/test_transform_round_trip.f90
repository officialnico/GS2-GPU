program test_transform_round_trip
  use iso_fortran_env, only: output_unit
  use unit_tests
  use gs2_main, only: gs2_program_state_type
  use mp, only: init_mp, finish_mp
  use mp, only: proc0, max_allreduce
  use gs2_main, only: initialize_wall_clock_timer
  use gs2_main, only: initialize_gs2
  use gs2_main, only: initialize_equations
  use gs2_main, only: finalize_equations
  use gs2_main, only: finalize_gs2
  use dist_fn_arrays, only: gnew, g
  use gs2_layouts, only: yxf_lo, g_lo, il_idx, ik_idx
  use le_grids, only: forbid
  use gs2_transforms, only: transform2, inverse2
  implicit none
  type(gs2_program_state_type) :: state
  logical :: agrees
  real :: error = 0.0
  real, dimension(:, :), allocatable :: yxf
  real, parameter :: tolerance = epsilon(error)
  character(len=40) :: name
  integer :: i, iglo, il, ik

  call init_mp
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call initialize_equations(state)

  call announce_module_test('transform_round_trip')

  do iglo = g_lo%llim_proc, g_lo%ulim_proc
     il = il_idx(g_lo, iglo)
     where (.not. forbid(:, il))
        g(:, 1, iglo) = 1.0
        g(:, 2, iglo) = 1.0
     elsewhere
        g(:, 1, iglo) = 0.0
        g(:, 2, iglo) = 0.0
     end where
  end do

  gnew = g
  allocate(yxf(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
  call transform2(g, yxf)
  call inverse2(yxf, g)

  error = maxval(abs(gnew - g))
  call max_allreduce(error)
  if(proc0) print*,'Error:',error, maxval(abs(gnew)), maxval(abs(g))
  ! Initialise local variables
  agrees = error <= tolerance

  name = "Plain round trip"

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  call process_test(agrees,'transform_round_trip -- '//trim(name))

  call close_module_test('transform_round_trip')

  call finalize_equations(state)
  call finalize_gs2(state)
  call finish_mp

end program test_transform_round_trip
