!> A regression test for the old diagnostics module. This runs 
!> a linear cyclone test case and then checks the value of diffusivity
!> against the value produced by svn revision 4242.
!>
!> This is free software released under the MIT license
!>   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program test_gs2_diagnostics
  use gs2_main, only: run_gs2, finish_gs2, gs2_program_state_type
  use unit_tests
  use file_utils, only: run_name, error_unit
  use mp, only: init_mp, test_driver_flag, mp_abort, mp_comm, finish_mp
  implicit none

  type(gs2_program_state_type) :: state
  real :: eps

  ! We pick a loose tolerance here as we know that the test
  ! is not very well resolved.
  eps = 1.0e-2

  call init_mp()
  state%mp_comm = mp_comm

  test_driver_flag = .true.
  functional_test_flag = .true.

  call announce_module_test("gs2_diagnostics")
  call run_gs2(finalize=.false., state=state)
  call announce_test('diffusivity')

  if( trim(run_name) .eq. 'test_gs2_diagnostics' ) then
     call process_test(test_diffusivity(1.25, eps), 'diffusivity')
  else if( trim(run_name) .eq. 'test_gs2_diagnostics_nspec_2' ) then
     call process_test(test_diffusivity(2.53, eps), 'diffusivity')
  else
     call mp_abort('test_gs2_diagnostics calling invalid test name')
  end if

  call finish_gs2(state, quiet=.true.)
  call close_module_test("gs2_diagnostics")
  call finish_mp

contains

  !> Return true if calculated diffusivity is within tolerance
  function test_diffusivity(results, eps)
    use unit_tests, only: agrees_with
    use mp,only: proc0
    use gs2_diagnostics, only: diffusivity
    real, intent(in) :: results, eps
    logical :: test_diffusivity
    test_diffusivity = .true.
    if (proc0) then
       test_diffusivity = agrees_with(diffusivity(), results, eps)
       if (.not. test_diffusivity) write(*,*) "Expecting",results,"and got",diffusivity(),"with tolerance",eps
    end if
  end function test_diffusivity

end program test_gs2_diagnostics
