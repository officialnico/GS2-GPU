program box_to_ballooning_space
  use gs2_main, only: gs2_program_state_type, initialize_wall_clock_timer
  use mp, only: init_mp, finish_mp, mp_abort, proc0
  use gs2_main, only: initialize_gs2, finalize_gs2
  use gs2_init, only: init, init_level_list
  use kt_grids, only: box, naky, ntheta0, aky, akx, theta0
  use dist_fn, only: get_leftmost_it, n_supercells, n_cells, supercell_labels
  use theta_grid, only: nperiod
  use sorting, only: quicksort
  use config_collection, only: gs2_config_type
  use file_utils, only: open_output_file, close_output_file, run_name
  implicit none

  type :: supercell_properties
     real :: ky
     real, dimension(:), allocatable :: kxs, theta0s
     real :: kx_middle, theta0_middle
     integer :: ncell, nperiod
  end type supercell_properties

  type(gs2_program_state_type) :: state
  type(gs2_config_type) :: the_config
  type(supercell_properties), dimension(:), allocatable :: supercells
  integer :: ik, it, is, ic, supercell_index, nsupercells, ncell, unit, supercell_length
  integer, dimension(:), allocatable :: temporary_values, it_values, it_range
  real, dimension(:), allocatable :: kxs, theta0s
  character(len=64) :: extension
  logical :: shift_cells
  logical, parameter :: debug = .true.
  logical, parameter :: write_full_file = .true.

  call initialize_wall_clock_timer
  call init_mp(state%mp_comm)
  call initialize_gs2(state)

  ! Need to init to dist_fn_level_1 to ensure boundary
  ! conditions have been set up.
  call init(state%init, init_level_list%dist_fn_level_1)

  ! Exit if this isn't a box mode input
  if (.not. box) then
     call mp_abort('Cannot run box_to_ballooning_space on a non-box input.', .true.)
  end if

  allocate(it_range(ntheta0))
  it_range = [(it, it = 1, ntheta0)]

  nsupercells = sum(n_supercells)

  allocate(supercells(nsupercells))
  allocate(temporary_values(ntheta0))

  supercell_index = 0
  do ik = 1, naky
     temporary_values = supercell_labels(:, ik)
     do it = 1, ntheta0
        if (temporary_values(it) >= 0) then
           supercell_index = supercell_labels(it, ik)

           supercells(supercell_index)%ky = aky(ik)

           ncell = n_cells(it, ik)
           supercells(supercell_index)%ncell = ncell

           supercell_length =  (2 * nperiod - 1 ) * ncell

           ! If we have an even number of cells then we need to add
           ! one more to the supercell length as we can only have an
           ! odd number of 2 pi theta sections in ballooning space
           shift_cells = mod(ncell, 2) == 0
           if (shift_cells) supercell_length = supercell_length + 1

           ! What ballooning space nperiod gives the same supercell_length
           ! supercell_length = 2 * nperiod_bs - 1
           supercells(supercell_index)%nperiod = (supercell_length + 1) / 2

           allocate(kxs(ncell), theta0s(ncell), it_values(ncell))

           it_values = pack(it_range, temporary_values == temporary_values(it))

           do ic = 1, ncell
              kxs(ic) = akx(it_values(ic))
              theta0s(ic) = theta0(it_values(ic), ik)
           end do

           ! These are the raw kxs/theta0s which make up this supercell
           supercells(supercell_index)%kxs = kxs
           supercells(supercell_index)%theta0s = theta0s

           ! To identify the theta0/kx of the ballooning space run we need to
           ! get the value at the middle of the supercell. This is easy for
           ! supercells with an odd number of cells, but for even counts we
           ! will be adding an extra cell so need to take some more care. To assist
           ! with this we will sort the kxs/theta0s values and then index the sorted
           ! arrays at the middle location for odd ncell and rounding up for even ncell
           call quicksort(ncell, kxs)
           call quicksort(ncell, theta0s)
           if (shift_cells) then
              ic = 1 + ncell / 2
           else
              ic = nint( 0.5 + (ncell / 2.0))
           end if

           supercells(supercell_index)%kx_middle = kxs(ic)
           supercells(supercell_index)%theta0_middle = theta0s(ic)
           deallocate(kxs, theta0s, it_values)

           where (temporary_values == temporary_values(it))
              temporary_values = -1
           end where

        end if
     end do
  end do

  ! Debug output
  if (proc0 .and. debug) then
     do supercell_index = 1, nsupercells
        print*,'-------------------------------------------------------------'
        print*,'Supercell index:',supercell_index
        print*,'ky:',supercells(supercell_index)%ky
        print*,'akx/theta0:',supercells(supercell_index)%kx_middle, supercells(supercell_index)%theta0_middle
        print*,'kxs:',supercells(supercell_index)%kxs
        print*,'theta0s:',supercells(supercell_index)%theta0s
        print*,'nperiod:',supercells(supercell_index)%nperiod
        print*,'ncell:',supercells(supercell_index)%ncell
        print*,'-------------------------------------------------------------'
        print*,''
     end do
  end if

  ! Create ballooning space input files
  if (proc0) then
     ! Get the setup from the existing file
     call the_config%get_configs

     ! Override some general settings
     the_config%nonlinear_terms_config%nonlinear_mode = 'off'
     the_config%kt_grids_config%grid_option = 'single'

     do supercell_index = 1, nsupercells

        ! Setup this specific case
        the_config%kt_grids_single_config%aky = supercells(supercell_index)%ky
        the_config%kt_grids_single_config%akx = supercells(supercell_index)%kx_middle
        the_config%kt_grids_single_config%theta0 = supercells(supercell_index)%theta0_middle
        the_config%theta_grid_parameters_config%nperiod = supercells(supercell_index)%nperiod

        write(extension,'(".ballooning_supercell.",I0,".in")') supercell_index
        call open_output_file(unit, extension)

        if (write_full_file) then
           call the_config%write_to_unit(unit)
        else
           call the_config%nonlinear_terms_config%write(unit)
           call the_config%kt_grids_config%write(unit)
           call the_config%kt_grids_single_config%write(unit)
           call the_config%theta_grid_parameters_config%write(unit)
           write(unit,'(A,"!include ",A,".in")') new_line('a'),trim(run_name)
        end if

        call close_output_file(unit)

     end do
  end if

  ! Clean up
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)
  call finish_mp

end program box_to_ballooning_space
