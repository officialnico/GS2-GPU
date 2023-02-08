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
  use gs2_layouts, only: yxf_lo, g_lo, il_idx, ig_idx, xxf_lo
  use le_grids, only: forbid
  use gs2_transforms, only: transform2, inverse2
  use gs2_transforms, only: transform_x, transform_y, inverse_x, inverse_y
  implicit none
  type(gs2_program_state_type) :: state
  logical :: agrees, overall_result
  real :: error = 0.0
  real, dimension(:, :), allocatable :: yxf
  complex, dimension(:, :), allocatable :: xxf, xxf_original
  real, parameter :: tolerance = epsilon(error) * 100.0
  character(len=40) :: name
  integer :: i, iglo, il, ig

  call init_mp
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call initialize_equations(state)

  call announce_module_test('transform_round_trip')
  overall_result = .true.

  do iglo = g_lo%llim_proc, g_lo%ulim_proc
     il = il_idx(g_lo, iglo)
     where (.not. forbid(:, il))
        g(:, 1, iglo) = 10.0
        g(:, 2, iglo) = 10.0
     elsewhere
        g(:, 1, iglo) = 0.0
        g(:, 2, iglo) = 0.0
     end where
  end do
  gnew = g

  name = "Plain round trip"
  allocate(yxf(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
  call transform2(g, yxf)
  call inverse2(yxf, g)

  error = maxval(abs(gnew - g))
  call max_allreduce(error)
  if(proc0) print*,name,' Error:',error, maxval(abs(gnew)), maxval(abs(g))
  agrees = error <= tolerance
  overall_result = overall_result .and. agrees

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  name = "x-only, manual scaling"
  allocate(xxf(xxf_lo%nx, xxf_lo%llim_proc:xxf_lo%ulim_alloc))
  g = gnew
  call transform_x(g, xxf)
  call inverse_x(xxf, g)
  ! Currently the below scaling causes the test to fail, for some reason no scaling seems required
  g = g / xxf_lo%nx

  error = maxval(abs(gnew - g))
  call max_allreduce(error)
  if(proc0) print*,name,' Error:',error, maxval(abs(gnew)), maxval(abs(g)), xxf_lo%nx, yxf_lo%ny

  agrees = error <= tolerance
  overall_result = overall_result .and. agrees

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  name = "y-only, manual scaling"
  allocate(xxf_original(xxf_lo%nx, xxf_lo%llim_proc:xxf_lo%ulim_alloc))
  xxf = 0.0
  do iglo = xxf_lo%llim_proc, xxf_lo%ulim_proc
     il = il_idx(xxf_lo, iglo)
     ig = ig_idx(xxf_lo, iglo)
     if (.not. forbid(ig, il)) xxf(:, iglo) = 10.0
  end do
  xxf_original = xxf
  call transform_y(xxf, yxf)
  call inverse_y(yxf, xxf)
  xxf = xxf / yxf_lo%ny

  error = maxval(abs(xxf_original - xxf))
  call max_allreduce(error)
  if(proc0) print*,name,' Error:',error, maxval(abs(xxf_original)), maxval(abs(xxf)), xxf_lo%nx, yxf_lo%ny

  agrees = error <= tolerance
  overall_result = overall_result .and. agrees

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  name = "x+y, manual scaling"
  g = gnew
  call transform_x(g, xxf)
  call transform_y(xxf, yxf)
  call inverse_y(yxf, xxf)
  call inverse_x(xxf, g)
  g = g / ( yxf_lo%ny * xxf_lo%nx)

  error = maxval(abs(gnew - g))
  call max_allreduce(error)
  if(proc0) print*,name,' Error:',error, maxval(abs(gnew)), maxval(abs(g)), xxf_lo%nx, yxf_lo%ny

  agrees = error <= tolerance
  overall_result = overall_result .and. agrees

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with error = ", E12.5)') &
          trim(adjustl(name)), error
  end if

  call process_test(overall_result,'transform_round_trip')

  call close_module_test('transform_round_trip')

  call finalize_equations(state)
  call finalize_gs2(state)
  call finish_mp

end program test_transform_round_trip
