program test_field_equations_consistency
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
  use dist_fn, only: get_init_field, getfieldeq
  use dist_fn_arrays, only: gnew
  use fields_arrays, only: phi, apar, bpar
  use gs2_layouts, only: g_lo
  implicit none
  type(gs2_program_state_type) :: state
  complex, dimension(:, :, :), allocatable :: fq, fqa, fqp  
  real :: local_max, local_max_fq, local_max_fqa, local_max_fqp
  logical :: agrees
  real, parameter :: tolerance = epsilon(0.0)*1e2
  character(len=40) :: name
  
  call init_mp
  call initialize_wall_clock_timer
  call initialize_gs2(state)
  call initialize_equations(state)

  call announce_module_test('field_equations_consistency')

  ! Initialise local variables
  agrees = .true.

  phi = 0.0 ; apar = 0.0  ; bpar = 0.0
  allocate(fq, source = phi)
  allocate(fqa, source = apar)
  allocate(fqp, source = bpar)

  name = "initial condition"
  call get_init_field(phi, apar, bpar)
  call getfieldeq(phi, apar, bpar, fq, fqa, fqp)

  local_max_fq = maxval(abs(fq))
  local_max_fqa = maxval(abs(fqa))
  local_max_fqp = maxval(abs(fqp))
  
  call max_allreduce(local_max_fq)
  call max_allreduce(local_max_fqa)
  call max_allreduce(local_max_fqp)
  local_max = maxval([local_max_fq, local_max_fqa, local_max_fqp])  
  call max_allreduce(local_max)
  agrees = local_max <= tolerance

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with errors:")') trim(adjustl(name))
     write(output_unit, '("fq  :",E12.5)') local_max_fq
     write(output_unit, '("fqa :",E12.5)') local_max_fqa
     write(output_unit, '("fqp :",E12.5)') local_max_fqp
  end if
  
  call process_test(agrees,'field_equations_consistency -- '//trim(name))
  
  ! Now repeat but with slightly different gnew to allow
  ! for finite apar, for example.
  name = "modified g"
  phi = 0.0 ; apar = 0.0  ; bpar = 0.0
  fq = 0.0 ; fqa = 0.0 ; fqp = 0.0
  gnew(:, 2, :) = -gnew(:, 2, :)

  call get_init_field(phi, apar, bpar)
  call getfieldeq(phi, apar, bpar, fq, fqa, fqp)

  local_max_fq = maxval(abs(fq))
  local_max_fqa = maxval(abs(fqa))
  local_max_fqp = maxval(abs(fqp))
  
  call max_allreduce(local_max_fq)
  call max_allreduce(local_max_fqa)
  call max_allreduce(local_max_fqp)
  local_max = maxval([local_max_fq, local_max_fqa, local_max_fqp])
  call max_allreduce(local_max)
  agrees = local_max <= tolerance

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with errors:")') trim(adjustl(name))
     write(output_unit, '("fq  :",E12.5)') local_max_fq
     write(output_unit, '("fqa :",E12.5)') local_max_fqa
     write(output_unit, '("fqp :",E12.5)') local_max_fqp
  end if
  
  call process_test(agrees,'field_equations_consistency -- '//trim(name))  

  ! Now repeat but with slightly different gnew to allow
  ! for finite apar, phi and bpar, for example.
  name = "flipped g"
  phi = 0.0 ; apar = 0.0  ; bpar = 0.0
  fq = 0.0 ; fqa = 0.0 ; fqp = 0.0
  gnew(:, 2, :) = gnew(:, 2, g_lo%ulim_proc:g_lo%llim_proc:-1)

  call get_init_field(phi, apar, bpar)
  call getfieldeq(phi, apar, bpar, fq, fqa, fqp)

  local_max_fq = maxval(abs(fq))
  local_max_fqa = maxval(abs(fqa))
  local_max_fqp = maxval(abs(fqp))
  
  call max_allreduce(local_max_fq)
  call max_allreduce(local_max_fqa)
  call max_allreduce(local_max_fqp)
  local_max = maxval([local_max_fq, local_max_fqa, local_max_fqp])
  call max_allreduce(local_max)
  agrees = local_max <= tolerance

  if (.not. agrees .and. proc0) then
     write(output_unit, '("Test case ",A," failed with errors:")') trim(adjustl(name))
     write(output_unit, '("fq  :",E12.5)') local_max_fq
     write(output_unit, '("fqa :",E12.5)') local_max_fqa
     write(output_unit, '("fqp :",E12.5)') local_max_fqp
  end if
  
  call process_test(agrees,'field_equations_consistency -- '//trim(name))  
  
  call close_module_test('field_equations_consistency')
  
  call finalize_equations(state)
  call finalize_gs2(state)
  call finish_mp
end program test_field_equations_consistency

