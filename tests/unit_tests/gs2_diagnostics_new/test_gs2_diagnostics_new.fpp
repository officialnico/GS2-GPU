#define CONCAT //

!> A program that tests the new diagnostics module. It  runs 
!! a  linear cyclone test case and then checks that the old and
!! new diagnostics give the same results
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)

module checks_mod
  use unit_tests
  public checks
  contains
    function checks()
      logical :: checks
      checks = .true.
    end function checks
end module checks_mod

program test_gs2_diagnostics_new
  !use functional_tests
  !use checks_mod
  !call test_gs2('Linear CBC (unit test) to test new diagnostics', checks)
    use gs2_main, only: run_gs2, finish_gs2, gs2_program_state_type
    use unit_tests
    use unit_tests, only: should_print
#ifdef NEW_DIAG
    use diagnostics_config, only: override_screen_printout_options
#endif
    use mp, only: init_mp, mp_comm, proc0, finish_mp
    use mp, only: broadcast
    use gs2_diagnostics, only: finish_gs2_diagnostics
    use gs2_diagnostics, only: pflux_avg, qflux_avg, heat_avg, vflux_avg
    use gs2_diagnostics, only: diffusivity
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
    use gs2_diagnostics_new, only: gnostics
#endif
    use run_parameters, only: use_old_diagnostics, nstep
    use species, only: nspec
    implicit none
    type(gs2_program_state_type) :: state
    integer :: n_vars
    integer :: i
    real :: eps
    real, dimension(:), allocatable :: pfluxav, qfluxav, heatav, vfluxav
    real :: diff
    !real :: vfluxav

    character(len=40), dimension(200) :: variables, new_variables, n_lines
    character(len=*), dimension(*), parameter :: file_names = [character(len=9):: &
         'heat',&
         'heat2', &
         'lpc', &
         'vres', &
         'phase', &
         'parity', &
         'eigenfunc' &
    ]

    eps = 1.0e-7
    if (precision(eps).lt. 11) eps = eps * 100.0

    ! n_lines - number of lines at the end of file "run_name.variables(i)" to
    !           compare as part of test
    variables(1) = 'lambda'
    n_lines(1) = '4'
    variables(2) = 'phi'
    n_lines(2) = '78'
    variables(3) = 'kx'
    n_lines(3) = '3'
    variables(4) = 'phi2_by_mode'
    n_lines(4) = '10'
    variables(5) = 't'
    n_lines(5) = '5'
    variables(6) = 'phi2_by_ky'
    n_lines(6) = '12'
    variables(7) = 'phi2_by_kx'
    n_lines(7) = '12'
    variables(8) = 'phi2_by_kx'
    n_lines(8) = '8'
    variables(9) = 'gds21'
    n_lines(9) = '4'
    variables(10) = 'bmag'
    n_lines(10) = '4'
    new_variables(1:10) = variables(1:10)

    variables(11) = 'es_heat_by_k'
    new_variables(11) = 'es_heat_flux_by_mode'
    n_lines(11) = '20'
    variables(12) = 'hflux_tot'
    new_variables(12) = 'heat_flux_tot'
    n_lines(12) = '4'
    variables(13) = 'vflux_tot'
    new_variables(13) = 'mom_flux_tot'
    n_lines(13) = '4'

    variables(14) = 'omega'
    new_variables(14) = 'omega'
    n_lines(14) = '4'
    variables(15) = 'omegaavg'
    new_variables(15) = 'omega_average'
    n_lines(15) = '4'
    
    variables(16) = 'ntot00'
    new_variables(16) = 'ntot_flxsurf_avg'
    n_lines(16) = '20'

    variables(17) = 'phi0'
    new_variables(17) = 'phi_igomega_by_mode'
    n_lines(17) = '15'

    variables(18) = 'apar_heat_by_k'
    new_variables(18) = 'apar_heat_flux_by_mode'
    n_lines(18) = '20'

    variables(19) = 'phase'
    new_variables(19) = 'phase'
    n_lines(19) = '20'

    variables(20) = 'phi_norm'
    new_variables(20) = 'phi_norm'
    n_lines(20) = '20'

    variables(21) = 'apar_norm'
    new_variables(21) = 'apar_norm'
    n_lines(21) = '20'

    variables(22) = 'bpar_norm'
    new_variables(22) = 'bpar_norm'
    n_lines(22) = '20'

    n_vars = 22

    call init_mp

  ! Here we switch on print_line and print_flux_line if we have 
  ! high verbosity... tests if they are working.

#ifdef NEW_DIAG  
  if (proc0) override_screen_printout_options = should_print(3)
  call broadcast(override_screen_printout_options)
#endif

  functional_test_flag = .true.
  state%mp_comm = mp_comm

  call announce_module_test("gs2_diagnostics_new")

  call run_gs2(finalize=.false., state=state)


  !call announce_test('results')
  !call process_test(test_function(), 'results')

#ifdef NEW_DIAG
    if (proc0) then
      if (use_old_diagnostics) then 
        open(120349, file='averages.dat', status="replace", action="write")
        write(120349, *) qflux_avg
        write(120349, *) pflux_avg
        write(120349, *) vflux_avg
        write(120349, *) heat_avg
        write(120349, *) diffusivity()
        close(120349)
      else
        if (nstep.eq.200) then
          allocate(qfluxav(nspec), pfluxav(nspec), heatav(nspec), vfluxav(nspec))
          open(120349, file='averages.dat')
          read(120349, *) qfluxav
          read(120349, *) pfluxav
          read(120349, *) vfluxav
          read(120349, *) heatav
          read(120349, *) diff
          close(120349)
          call announce_test("average heat flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_heat_flux_avg, qfluxav, eps), &
            "average heat flux")
          call announce_test("average momentum flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_momentum_flux_avg, vfluxav, eps), &
            "average momentum flux")
          call announce_test("average particle flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_particle_flux_avg, pfluxav, eps), &
            "average particle flux")
          call announce_test("diffusivity")
          call process_test( &
            agrees_with(gnostics%current_results%diffusivity, diff, eps), &
            "diffusivity")
        else
          call announce_test("Size of t array")
        end if
      end if
    end if
#endif

    call finish_gs2(state)

#ifdef NEW_DIAG
    if (proc0 .and. .not. use_old_diagnostics) then
      if (nstep==200) then 
        do i = 1,n_vars
          call announce_test("value of "//trim(new_variables(i)))
          call process_test(test_variable(trim(variables(i)), trim(new_variables(i))), &
            "value of "//trim(new_variables(i)))
        end do
      else if (gnostics%appending) then 
        do i = 1,n_vars
          ! omega and omega_average won't work because the 
          ! history is not stored in the restart file
          if (i.eq.14.or.i.eq.15) cycle
          call announce_test("value of "//trim(new_variables(i)))
          call process_test(&
            test_variable(trim(new_variables(i)), trim(new_variables(i))), &
            "value of "//trim(new_variables(i)))
        end do
      end if
      if (nstep==200) then 
        ! Note that heat and heat2 no longer pass because the comparison
        ! is no longer between files output from the same run but between
        ! files output from two successive runs. This causes changes greater
        ! than the precision of the heat and heat2 files. EGH
        do i = 3, size(file_names)
          call announce_test("content of "//trim(file_names(i)))
          call process_test(test_file(trim(file_names(i))), "content of "//trim(file_names(i)))
        end do
      end if
    end if
#endif

    call close_module_test("gs2_diagnostics_new")

    call finish_mp

contains
  
  function test_variable(var_name, new_var_name)
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif
    use unit_tests, only: should_print
    character(*), intent(in) :: var_name, new_var_name
    logical :: test_variable
    character(1000) ::  command 
    
#ifdef NEW_DIAG
    command = ''
    if (gnostics%append_old) then 
      ! Here we are comparing the new diagnostics 200 step run with 
      ! a 100-step run with an appended 100 step run.
      write(command,'("./compare.py ", A," ", A," ", A," ", A)') &
         "test_gs2_diagnostics_new.out.nc", var_name,&
         "test_gs2_diagnostics_new_append.out.nc", new_var_name
    else 
      write(command,'("./compare.py ", A," ", A," ", A," ", A)') &
         "old_diagnostics/test_gs2_diagnostics_new.out.nc", var_name,&
         "test_gs2_diagnostics_new.out.nc", new_var_name
     end if
    
    if (should_print(3)) write(*,*) trim(command)
    call execute_command_line(command)
    test_variable = .true.
    open(120349, file='tmpdata.dat')
    read(120349, '(L1)') test_variable
    close(120349)
#endif
  end function test_variable

  logical function test_file(file_name)
    use unit_tests, only: should_print
    character(*), intent(in) :: file_name
    character(1000) ::  command
    
    test_file=.true.
#ifdef NEW_DIAG
    ! bash command that compares the old and new diagnostics results for the last
    ! n_lines lines of a file, and writes T to tmpdata.dat if they agree
    write(command, &
         '("./compare_ascii.py old_diagnostics/test_gs2_diagnostics_new.", A," test_gs2_diagnostics_new.", A)') &
         file_name, file_name
    if (should_print(3)) write(*,*) trim(command)
    call execute_command_line(command)
    open(120349, file='tmpdata.dat')
    read(120349, '(L1)') test_file
    close(120349)
#endif
  end function test_file

end program test_gs2_diagnostics_new
