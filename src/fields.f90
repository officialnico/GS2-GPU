!> FIXME : Add documentation
module fields
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use text_options, only: text_option

  implicit none

  private 

  public :: init_fields, finish_fields
  public :: read_parameters, wnml_fields, check_fields
  public :: advance, force_maxwell_reinit
  public :: reset_init, set_init_fields
  public :: fields_init_response, set_dump_and_read_response
  public :: dump_response_to_file
  public :: init_fields_parameters
  public :: finish_fields_parameters
  public :: set_overrides
  public :: init_fields_level_1, init_fields_level_2
  public :: finish_fields_level_1, finish_fields_level_2
  public :: fieldopt_switch, fieldopt_implicit, fieldopt_test, fieldopt_local, fieldopt_gf_local

  !> Made public for unit tests
  public :: fields_pre_init
  public :: remove_zonal_flows_switch
  !> Made public for replay
  public :: allocate_arrays

  public :: fields_config_type
  public :: set_fields_config
  
  !> FIXME : Add documentation
  interface fieldlineavgphi
     module procedure fieldlineavgphi_loc
     module procedure fieldlineavgphi_tot
  end interface

  ! knobs
  integer :: fieldopt_switch
  logical :: remove_zonal_flows_switch
  logical :: force_maxwell_reinit
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2, fieldopt_local = 3, fieldopt_gf_local = 4
  logical :: dump_response, read_response
  logical :: initialized = .false.
  logical :: parameters_read = .false.
  logical :: exist

  !EGH made fieldopts module level for overrides
  type (text_option), dimension (6), parameter :: fieldopts = &
       (/ text_option('default', fieldopt_implicit), &
          text_option('implicit', fieldopt_implicit), &
          text_option('test', fieldopt_test),&
          text_option('local', fieldopt_local),&
          text_option('gf_local', fieldopt_gf_local),&
          text_option('implicit_local', fieldopt_local)/)
  
  !> Used to represent the input configuration of fields
  type, extends(abstract_config_type) :: fields_config_type
     ! namelist : fields_knobs
     ! indexed : false
     !> Used with `field_option='local'`. If `true` and x/y are
     !> distributed then in time advance only update local part of
     !> field in operations like `phinew=phinew+phi` etc.
     logical :: do_smart_update = .false.
     !> Writes files containing the field response matrix after
     !> initialisation. This currently works for
     !> `field_option='implicit'` or `'local'`.  We write to netcdf
     !> files by default but fall back to fortran unformatted (binary)
     !> files (which are not really portable) in the absence of
     !> netcdf.
     logical :: dump_response = .false.
     !> Set to true to use an allreduce (on `mp_comm`) in field
     !> calculation (`field_option='local'` only) rather than a
     !> reduction on a sub-communicator followed by a global
     !> broadcast. Typically a little faster than default performance
     !> but may depend on MPI implementation.
     logical :: field_local_allreduce = .false.
     !> Set to `true`, along with [[field_local_allreduce]] and
     !> [[layouts_knobs:intspec_sub]] , to replace the allreduce used
     !> in the field calculation with an allreduce on a
     !> sub-communicator followed by a reduction on a "perpendicular"
     !> communicator.  Typically a bit faster than default and scales
     !> slightly more efficiently.  Note if this option is active only
     !> `proc0` has knowledge of the full field arrays. Other
     !> processors know the full field for any supercell (connected
     !> x-y domains) for which it has any of the xy indices local in
     !> the `g_lo` layout.
     logical :: field_local_allreduce_sub = .false.
     !> If `true` then use nonblocking collective operations in the
     !> [[fields_local]] field calculation. This may or may not
     !> improve performance.
     logical :: field_local_nonblocking_collectives = .false.
     !> Set to `true` when using `field_option='local'` to
     !> automatically tune and select the best performing minimum
     !> block size (in a single supercell) assigned to a single
     !> processor. This can improve performance, but is not guaranteed
     !> to.
     logical :: field_local_tuneminnrow = .false.
     !> The `field_option` variable controls which time-advance
     !> algorithm is used for the linear terms. Allowed values are:
     !>
     !> - 'implicit' Advance linear terms with Kotschenreuther's implicit algorithm.
     !> - 'default' Same as 'implicit'.
     !> - 'implicit_local' the same as 'local'.
     !> - 'local' Same implicit algorithm as 'implicit' but with
     !>   different data decomposition (typically much faster for flux
     !>   tube runs).
     !> - 'gf_local' Same as `'local'` but with `gf_lo` field
     !>   decomposition. To use this you also need to set
     !>   `gf_lo_integrate= .true.` in [[dist_fn_knobs]] and
     !>   `gf_local_fields = .true.` in [[layouts_knobs]].
     !> - 'test' Use for debugging.     
     character(len = 20) :: field_option = 'default'
     !> Set to true to use allgatherv to fetch parts of the field update
     !> vector calculated on other procs. When false uses a sum_allreduce
     !> instead. This doesn't rely on sub-communicators so should work for
     !> any layout and processor count.  Note: This only impacts
     !> `field_option='implicit'`
     logical :: field_subgath = .false.
     !> If `true` then recalculate the fields from the distribution
     !> function using [[get_init_field]] when restarting a
     !> simulation rather than using the values in the restart file.
     !>
     !> @todo Consider if this should this go into the reinit namelist.
     logical :: force_maxwell_reinit = .true. 
     !> Used with `field_option='local'` to set the minimum block size
     !> (in a single supercell) assigned to a single processor. Tuning
     !> this parameter changes the balance between work
     !> parallelisation and communication. As this value is lowered,
     !> more communication needs to be done, but more processors get
     !> assigned work. This may reduce the time spent in computation
     !> at the cost of time spent in communication. The optimal value
     !> is likely to depend upon the size of the problem and the
     !> number of processors being used. Furthermore it will affect
     !> intialisation and advance in different ways. Can be
     !> automatically tuned using [[fields_knobs:field_local_tuneminnrow]].
     integer :: minnrow = 16
     !> Reads files containing the field response matrix and uses to
     !> initialise GS2s response matrix rather than using the usual
     !> initialisation process.
     logical :: read_response = .false.
     !> Delete zonal flows at every timestep.
     logical :: remove_zonal_flows_switch = .false.
     !> Sets location in which to store/look for response dump files.
     !> We don't currently check that this location exists before
     !> attempting to use it, which could cause problems. The default
     !> is to save them in the working directory.
     character(len = 256) :: response_dir = ''
     !> Allows customisation of the base filename to be used for
     !> response files. If not set then we use `run_name` derived
     !> from the input file name.
     character(len = 256) :: response_file = ''
   contains
     procedure, public :: read => read_fields_config
     procedure, public :: write => write_fields_config
     procedure, public :: reset => reset_fields_config
     procedure, public :: broadcast => broadcast_fields_config
     procedure, public, nopass :: get_default_name => get_default_name_fields_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_fields_config
  end type fields_config_type
  
  type(fields_config_type) :: fields_config  
contains

  !> FIXME : Add documentation  
  subroutine check_fields(report_unit)
    use fields_local, only: do_smart_update, minNRow
    implicit none
    integer, intent(in) :: report_unit
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
    case (fieldopt_test)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will only be tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    case (fieldopt_local)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly with decomposition respecting g_lo layout.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
       write(report_unit, fmt="('Using a min block size of ',I0)") minNrow
       if(do_smart_update) write(report_unit, fmt="('Using optimised field update.')")
    case (fieldopt_gf_local)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly with decomposition respecting gf_lo layout.')")
    end select
  end subroutine check_fields

  !> FIXME : Add documentation  
  subroutine wnml_fields(unit)
    use fields_local, only: do_smart_update, minNRow
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return 
    write (unit, *)
    write (unit, fmt="(' &',a)") "fields_knobs"
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (unit, fmt="(' field_option = ',a)") '"implicit"'
    case (fieldopt_test)
       write (unit, fmt="(' field_option = ',a)") '"test"'
    case (fieldopt_local)
       write (unit, fmt="(' field_option = ',a)") '"local"'
       write (unit, fmt="(' minNrow = ',I0)") minNrow
       write (unit, fmt="(' do_smart_update = ',L1)") do_smart_update
    case (fieldopt_gf_local)
       write (unit, fmt="(' field_option = ',a)") '"local gf"'
    end select
    if(dump_response) write (unit, fmt="(' dump_response = ',L1)") dump_response
    if(read_response) write (unit, fmt="(' read_response = ',L1)") read_response
    write (unit, fmt="(' /')")
  end subroutine wnml_fields

  !> FIXME : Add documentation
  subroutine set_overrides(opt_ov)
    use overrides, only: optimisations_overrides_type
    use file_utils, only: error_unit
    use text_options, only: get_option_value
    use fields_local, only: minnrow
    use fields_implicit, only: field_subgath
    use fields_local, only: do_smart_update
    use fields_local, only: field_local_allreduce
    use fields_local, only: field_local_allreduce_sub
    type(optimisations_overrides_type), intent(in) :: opt_ov
    integer :: ierr
    if (opt_ov%override_field_option) then
       ierr = error_unit()
       call get_option_value &
            (opt_ov%field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in set_overrides",.true.)
    end if
    if (opt_ov%override_minnrow) minnrow = opt_ov%minnrow
    if (opt_ov%override_field_subgath) field_subgath = opt_ov%field_subgath
    if (opt_ov%override_do_smart_update) do_smart_update = & 
      opt_ov%do_smart_update
    if (opt_ov%override_field_local_allreduce) field_local_allreduce = & 
      opt_ov%field_local_allreduce
    if (opt_ov%override_field_local_allreduce_sub) field_local_allreduce_sub =& 
      opt_ov%field_local_allreduce_sub
  end subroutine set_overrides

  !> Calls all initialisations required for init_fields_implicit/local, 
  !! reads parameters and allocates field arrays
  subroutine fields_pre_init(fields_config_in)
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn, gf_lo_integrate
    use init_g, only: ginit, init_init_g
    use antenna, only: init_antenna
    use unit_tests, only: debug_message
    use kt_grids, only: naky, ntheta0
    use mp, only: nproc, proc0
    implicit none
    type(fields_config_type), intent(in), optional :: fields_config_in    
    integer, parameter :: verb=3
    
    call init_fields_parameters(fields_config_in)

    call debug_message(verb, "init_fields: init_theta_grid")
    call init_theta_grid
    
!CMR,30/3/2009:
! call init_init_g before init_run_parameters to read delt from restart file

    call debug_message(verb, "init_fields: init_init_g")
    call init_init_g
    call debug_message(verb, "init_fields: init_run_parameters")
    call init_run_parameters
    call debug_message(verb, "init_fields: init_dist_fn")
    call init_dist_fn

    if(nproc .lt. ntheta0*naky .and. fieldopt_switch .eq. fieldopt_gf_local) then
       fieldopt_switch = fieldopt_local
       gf_lo_integrate = .false.
       if(proc0) then
          write(*,*) 'gf local fields cannot be used as you are running less MPI processes than there are'
          write(*,*) 'naky*ntheta0 points.  Defaulting to local fields.  You need to use at least',naky*ntheta0
          write(*,*) 'MPI processes for this simulation'
       end if
    end if

    !call debug_message(verb, "init_fields: init_parameter_scan")
    !call init_parameter_scan
    call debug_message(verb, "init_fields: init_antenna")
    call init_antenna !Must come before allocate_arrays so we know if we need apar_ext
    !call debug_message(verb, "init_fields: read_parameters")
    !call read_parameters
    call debug_message(verb, "init_fields: allocate_arrays")
    call allocate_arrays
  end subroutine fields_pre_init

  !> FIXME : Add documentation  
  subroutine init_fields_parameters(fields_config_in)
    use unit_tests, only: debug_message
    implicit none
    type(fields_config_type), intent(in), optional :: fields_config_in    
    integer, parameter :: verb=3
    if (parameters_read) return
    call debug_message(verb, "init_fields: read_parameters")
    call read_parameters(fields_config_in)
    parameters_read = .true.
  end subroutine init_fields_parameters

  !> FIXME : Add documentation  
  subroutine finish_fields_parameters
    parameters_read = .false.
  end subroutine finish_fields_parameters

  !> FIXME : Add documentation  
  subroutine init_fields_level_1
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    call debug_message(verb, "init_fields: allocate_arrays")
    call allocate_arrays
  end subroutine init_fields_level_1

  !> FIXME : Add documentation  
  subroutine finish_fields_level_1
    call finish_fields
  end subroutine finish_fields_level_1

  !> FIXME : Add documentation  
  subroutine init_fields_level_2
    call init_fields
  end subroutine init_fields_level_2

  !> FIXME : Add documentation  
  subroutine finish_fields_level_2
    call reset_init
  end subroutine finish_fields_level_2

  !> FIXME : Add documentation  
  subroutine fields_init_response
    use fields_implicit, only: init_fields_implicit
    use fields_test, only: init_fields_test
    use fields_local, only: init_fields_local
    use fields_gf_local, only: init_fields_gf_local
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_implicit")
       call init_fields_implicit
    case (fieldopt_test)
       call debug_message(verb, "fields::fields_init_response init_fields_test")
       call init_fields_test
    case (fieldopt_local)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_local")
       call init_fields_local
    case (fieldopt_gf_local)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_gf_local")
       call init_fields_gf_local
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine fields_init_response

  !> FIXME : Add documentation  
  subroutine init_fields(fields_config_in)
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use init_g, only: ginit, init_init_g
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    use kt_grids, only: gridopt_switch, gridopt_box, kwork_filter
    implicit none
    type(fields_config_type), intent(in), optional :: fields_config_in    
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    call fields_pre_init(fields_config_in)

    call fields_init_response

    ! Turn on nonlinear terms.
    if (debug) write(6,*) "init_fields: nl_finish_init"
    call nl_finish_init

    !If running in flux tube disable evolution of ky=kx=0 mode
    if(gridopt_switch.eq.gridopt_box) kwork_filter(1,1)=.true.
  end subroutine init_fields

  !> Force the current response matrices to be written to file
  subroutine dump_response_to_file(suffix)
    use fields_implicit, only: dump_response_to_file_imp
    use fields_local, only: dump_response_to_file_local
    use fields_gf_local, only: dump_response_to_file_gf_local
    implicit none
    character(len=*), intent(in), optional :: suffix 
    !Note can pass optional straight through as long as also optional
    !in called routine (and not different routines combined in interface)
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call dump_response_to_file_imp(suffix)
    case (fieldopt_local)
       call dump_response_to_file_local(suffix)
    case (fieldopt_gf_local)
       call dump_response_to_file_gf_local(suffix)
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine dump_response_to_file

  !> FIXME : Add documentation
  subroutine set_init_fields
    use fields_implicit, only: init_allfields_implicit
    use fields_test, only: init_phi_test
    use mp, only: proc0
    use fields_local, only: init_allfields_local
    use fields_gf_local, only: init_allfields_gf_local
    use dist_fn, only: gf_lo_integrate
    use gs2_layouts, only: gf_local_fields
    use kt_grids, only: naky, ntheta0
    implicit none
    logical, parameter :: debug=.false.
    if(proc0.and.debug) write(6,*) "Syncing fields with g."
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_allfields_implicit"
       call init_allfields_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_phi_test"
       call init_phi_test
    case (fieldopt_local)
       if (debug) write(6,*) "init_fields: init_allfields_local"
       call init_allfields_local
    case (fieldopt_gf_local)
       if (debug) write(6,*) "init_fields: init_allfields_gf_local"
       if(.not. gf_lo_integrate .or. .not. gf_local_fields) then
          if(proc0) then
             write(*,*) 'gf_lo_integrate',gf_lo_integrate,'gf_local_fields',gf_local_fields
             write(*,*) 'gf local fields cannot be used by gf_lo_integrate'
             write(*,*) 'defaulting to local fields'
             write(*,*) 'if you want to use gf local fields then set gf_lo_integrate to true in dist_fn_knobs'
             write(*,*) 'and gf_local_fields to true in layouts_knobs and make sure you are running on more MPI'
             write(*,*) 'mpi processes than naky*ntheta0.  For this simulation that is: ',naky*ntheta0
          end if
          fieldopt_switch = fieldopt_local
          call init_allfields_local
       else
          call init_allfields_gf_local
       end if
    end select
  end subroutine set_init_fields

  !> FIXME : Add documentation  
  subroutine read_parameters(fields_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use fields_implicit, only: field_subgath
    use fields_local, only: minNRow
    use fields_local, only: do_smart_update, field_local_allreduce, field_local_allreduce_sub
    use fields_local, only: field_local_tuneminnrow,  field_local_nonblocking_collectives
    use fields_arrays, only: real_response_file => response_file
    use file_utils, only: run_name
    implicit none
    type(fields_config_type), intent(in), optional :: fields_config_in    
    character(20) :: field_option
    character(len=256) :: response_dir
    character(len=256) :: response_file
    integer :: ierr


    if (present(fields_config_in)) fields_config = fields_config_in

    call fields_config%init(name = 'fields_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    do_smart_update = fields_config%do_smart_update
    dump_response = fields_config%dump_response
    field_local_allreduce = fields_config%field_local_allreduce
    field_local_allreduce_sub = fields_config%field_local_allreduce_sub
    field_local_nonblocking_collectives = fields_config%field_local_nonblocking_collectives
    field_local_tuneminnrow = fields_config%field_local_tuneminnrow
    field_option = fields_config%field_option
    field_subgath = fields_config%field_subgath
    force_maxwell_reinit = fields_config%force_maxwell_reinit
    minnrow = fields_config%minnrow
    read_response = fields_config%read_response
    remove_zonal_flows_switch = fields_config%remove_zonal_flows_switch
    response_dir = fields_config%response_dir
    response_file = fields_config%response_file

    exist = fields_config%exist
    ierr = error_unit()
    call get_option_value &
         (field_option, fieldopts, fieldopt_switch, &
         ierr, "field_option in fields_knobs",.true.)

    if (trim(response_file) == '') then
       response_file = trim(run_name)
    end if

    if(trim(response_dir).eq.'')then
       write(real_response_file,'(A)') trim(response_file)
    else
       write(real_response_file,'(A,"/",A)') trim(response_dir),trim(response_file)
    endif

    !Set the solve type specific flags
    call set_dump_and_read_response(dump_response, read_response)
  end subroutine read_parameters

  !> FIXME : Add documentation  
  subroutine set_dump_and_read_response(dump_flag, read_flag)
    use fields_implicit, only: dump_response_imp => dump_response, read_response_imp=>read_response
    use fields_local, only: dump_response_loc => dump_response, read_response_loc=>read_response
    use fields_gf_local, only: dump_response_gf => dump_response, read_response_gf=>read_response
    implicit none
    logical, intent(in) :: dump_flag, read_flag
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       dump_response_imp=dump_flag
       read_response_imp=read_flag
    case (fieldopt_local)
       dump_response_loc=dump_flag
       read_response_loc=read_flag
    case (fieldopt_gf_local)
       dump_response_gf=dump_flag
       read_response_gf=read_flag
    case default
       !Silently ignore unsupported field types
    end select
  end subroutine set_dump_and_read_response

  !> FIXME : Add documentation  
  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use antenna, only: no_driver
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew, apar_ext
    use fields_arrays, only: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3

    if (.not. allocated(phi)) then
       call debug_message(verb, 'fields::allocate_arrays allocating')
       allocate (     phi (-ntgrid:ntgrid,ntheta0,naky))
       allocate (    apar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (   bpar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phinew (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( aparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparnew (-ntgrid:ntgrid,ntheta0,naky))
       if(fieldopt_switch .eq. fieldopt_gf_local) then
          !AJ It should be possible to reduce the size of these by only allocating them
          !AJ the extend of it and ik in gf_lo.  However, it would need to be done carefully 
          !AJ to ensure the proc0 has the full space if we are still reducing data to proc0 
          !AJ for diagnostics.
          !AJ With some careful thought it should also be possible to remove this all together and 
          !AJ simply use the arrays above, but that would need checking
          allocate (    gf_phi (-ntgrid:ntgrid,ntheta0,naky))
          allocate (   gf_apar (-ntgrid:ntgrid,ntheta0,naky))
          allocate (   gf_bpar (-ntgrid:ntgrid,ntheta0,naky))          
          allocate ( gf_phinew (-ntgrid:ntgrid,ntheta0,naky))
          allocate (gf_aparnew (-ntgrid:ntgrid,ntheta0,naky))
          allocate (gf_bparnew (-ntgrid:ntgrid,ntheta0,naky))          
       end if
    endif
    !AJ This shouldn't be necessary as it is done in the fields code....
    phi = 0.; phinew = 0.
    apar = 0.; aparnew = 0.
    bpar = 0.; bparnew = 0.
    if(fieldopt_switch .eq. fieldopt_gf_local) then
       gf_phi = 0.; gf_phinew = 0.
       gf_apar = 0.; gf_aparnew = 0.
       gf_bpar = 0.; gf_bparnew = 0.
    end if
    if(.not.allocated(apar_ext).and.(.not.no_driver))then
       allocate (apar_ext (-ntgrid:ntgrid,ntheta0,naky))
       apar_ext = 0.
    endif
  end subroutine allocate_arrays

  !> FIXME : Add documentation  
  subroutine advance (istep)
    use fields_implicit, only: advance_implicit
    use fields_test, only: advance_test
    use fields_local, only: advance_local
    use fields_gf_local, only: advance_gf_local

    implicit none
    integer, intent (in) :: istep

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep, remove_zonal_flows_switch)
    case (fieldopt_test)
       call advance_test (istep)
    case (fieldopt_local)
       call advance_local (istep, remove_zonal_flows_switch)
    case (fieldopt_gf_local)
       call advance_gf_local (istep, remove_zonal_flows_switch)
    end select
  end subroutine advance

  !> FIXME : Add documentation
  !>
  !> @note This routine has a potentially misleading name and isn't used anywhere
  subroutine kperp (ntgrid_output, akperp)
    use theta_grid, only: theta
    use kt_grids, only: naky, aky, ntheta0, kperp2
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use integration, only: trapezoidal_integration
    implicit none
    integer, intent (in) :: ntgrid_output
    real, dimension (:,:), intent (out) :: akperp
    real :: anorm
    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          anorm = trapezoidal_integration(theta, abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                         + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                         + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                      )
          if (anorm < 2.0*epsilon(0.0) .or. aky(ik) == 0.0) then
             akperp(it,ik) = 0.0
          else
             akperp(it,ik) &
                  = sqrt(trapezoidal_integration(theta, kperp2(-ntgrid_output:ntgrid_output,it,ik) &
                     *abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                          + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                          + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                          )/anorm)
          end if
       end do
    end do
  end subroutine kperp

  !> FIXME : Add documentation
  !>
  !> @note This routine isn't used anywhere
  subroutine fieldlineavgphi_loc (ntgrid_output, it, ik, phiavg)
    use theta_grid, only: ntgrid, jacob, theta
    use fields_arrays, only: phi
    use integration, only: trapezoidal_integration
    implicit none
    integer, intent (in) :: ntgrid_output, ik, it
    complex, intent (out) :: phiavg

    phiavg = &
         trapezoidal_integration(theta(-ntgrid_output:ntgrid_output), &
         phi(-ntgrid_output:ntgrid_output,it,ik) &
         *jacob(-ntgrid_output:ntgrid_output)) &
         /trapezoidal_integration(theta(-ntgrid_output:ntgrid_output), &
         jacob(-ntgrid_output:ntgrid_output))

  end subroutine fieldlineavgphi_loc

  !> FIXME : Add documentation
  !>
  !> @note This doesn't look like a useful routine and isn't used anywhere
  subroutine fieldlineavgphi_tot (phiavg)
    implicit none
    complex, dimension (:,:), intent (out) :: phiavg
    phiavg = 0.
    write(*,*) 'error in fields'
  end subroutine fieldlineavgphi_tot

  !> FIXME : Add documentation  
  subroutine reset_init
    use fields_implicit, only: fi_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use fields_local, only: fl_reset => reset_fields_local
    use fields_gf_local, only: flgf_reset => reset_fields_gf_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
    implicit none
    initialized  = .false.
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    if(fieldopt_switch .eq. fieldopt_gf_local) then
       gf_phi = 0.
       gf_apar = 0.
       gf_bpar = 0.
       gf_phinew = 0.
       gf_aparnew = 0.
       gf_bparnew = 0.
    end if
    !What about apar_ext?
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call fi_reset
    case (fieldopt_test)
       call ft_reset
    case (fieldopt_local)
       call fl_reset
    case (fieldopt_gf_local)
       call flgf_reset
    end select
  end subroutine reset_init

  !> FIXME : Add documentation  
  subroutine finish_fields

    use fields_implicit, only: implicit_reset => reset_init
    use fields_test, only: test_reset => reset_init
    use fields_local, only: finish_fields_local
    use fields_gf_local, only: finish_fields_gf_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: apar_ext, gf_phi, gf_apar, gf_bpar
    use fields_arrays, only: apar_ext, gf_phinew, gf_aparnew, gf_bparnew
    use unit_tests, only: debug_message

    implicit none
    integer, parameter :: verbosity = 3


    initialized  = .false.
!AJ Does these need zero'd if they are to be deallocated below?
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    call debug_message(verbosity, &
      'fields::finish_fields zeroed fields')
    
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       ! This line is no longer necessary
       ! because fields_implicit::reset_init is
       ! called by finish_fields_level_2
       !call implicit_reset
    case (fieldopt_test)
       ! This line is no longer necessary
       ! because fields_test::reset_init is
       ! called by finish_fields_level_2
       !call test_reset
    case (fieldopt_local)
       call finish_fields_local
    case (fieldopt_gf_local)
       call finish_fields_gf_local
    end select

    call debug_message(verbosity, &
      'fields::finish_fields called subroutines')

    if (allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew)
    if (allocated(gf_phi)) deallocate(gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew)
    if (allocated(apar_ext)) deallocate (apar_ext)
    call debug_message(verbosity, &
      'fields::finish_fields deallocated fields')

    call fields_config%reset()
  end subroutine finish_fields

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_fields_config(fields_config_in)
    use mp, only: mp_abort
    type(fields_config_type), intent(in), optional :: fields_config_in
    if (initialized) then
       call mp_abort("Trying to set fields_config when already initialized.", to_screen = .true.)
    end if
    if (present(fields_config_in)) then
       fields_config = fields_config_in
    end if
  end subroutine set_fields_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the fields_knobs namelist and populates the member variables
  subroutine read_fields_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(fields_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: field_option
    character(len = 256) :: response_dir, response_file
    integer :: minnrow
    logical :: do_smart_update, dump_response, field_local_allreduce, field_local_allreduce_sub, field_local_nonblocking_collectives, field_local_tuneminnrow, field_subgath, force_maxwell_reinit, read_response
    logical :: remove_zonal_flows_switch

    namelist /fields_knobs/ do_smart_update, dump_response, field_local_allreduce, field_local_allreduce_sub, field_local_nonblocking_collectives, field_local_tuneminnrow, field_option, field_subgath, force_maxwell_reinit, &
         minnrow, read_response, remove_zonal_flows_switch, response_dir, response_file

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    do_smart_update = self%do_smart_update
    dump_response = self%dump_response
    field_local_allreduce = self%field_local_allreduce
    field_local_allreduce_sub = self%field_local_allreduce_sub
    field_local_nonblocking_collectives = self%field_local_nonblocking_collectives
    field_local_tuneminnrow = self%field_local_tuneminnrow
    field_option = self%field_option
    field_subgath = self%field_subgath
    force_maxwell_reinit = self%force_maxwell_reinit
    minnrow = self%minnrow
    read_response = self%read_response
    remove_zonal_flows_switch = self%remove_zonal_flows_switch
    response_dir = self%response_dir
    response_file = self%response_file

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = fields_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = fields_knobs)
    endif

    ! Now copy from local variables into type members
    self%do_smart_update = do_smart_update
    self%dump_response = dump_response
    self%field_local_allreduce = field_local_allreduce
    self%field_local_allreduce_sub = field_local_allreduce_sub
    self%field_local_nonblocking_collectives = field_local_nonblocking_collectives
    self%field_local_tuneminnrow = field_local_tuneminnrow
    self%field_option = field_option
    self%field_subgath = field_subgath
    self%force_maxwell_reinit = force_maxwell_reinit
    self%minnrow = minnrow
    self%read_response = read_response
    self%remove_zonal_flows_switch = remove_zonal_flows_switch
    self%response_dir = response_dir
    self%response_file = response_file

    self%exist = exist
  end subroutine read_fields_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_fields_config(self, unit)
    implicit none
    class(fields_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("do_smart_update", self%do_smart_update, unit_internal)
    call self%write_key_val("dump_response", self%dump_response, unit_internal)
    call self%write_key_val("field_local_allreduce", self%field_local_allreduce, unit_internal)
    call self%write_key_val("field_local_allreduce_sub", self%field_local_allreduce_sub, unit_internal)
    call self%write_key_val("field_local_nonblocking_collectives", self%field_local_nonblocking_collectives, unit_internal)
    call self%write_key_val("field_local_tuneminnrow", self%field_local_tuneminnrow, unit_internal)
    call self%write_key_val("field_option", self%field_option, unit_internal)
    call self%write_key_val("field_subgath", self%field_subgath, unit_internal)
    call self%write_key_val("force_maxwell_reinit", self%force_maxwell_reinit, unit_internal)
    call self%write_key_val("minnrow", self%minnrow, unit_internal)
    call self%write_key_val("read_response", self%read_response, unit_internal)
    call self%write_key_val("remove_zonal_flows_switch", self%remove_zonal_flows_switch, unit_internal)
    call self%write_key_val("response_dir", self%response_dir, unit_internal)
    call self%write_key_val("response_file", self%response_file, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_fields_config

  !> Resets the config object to the initial empty state
  subroutine reset_fields_config(self)
    class(fields_config_type), intent(in out) :: self
    type(fields_config_type) :: empty
    select type (self)
    type is (fields_config_type)
       self = empty
    end select
  end subroutine reset_fields_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_fields_config(self)
    use mp, only: broadcast
    implicit none
    class(fields_config_type), intent(in out) :: self
    call broadcast(self%do_smart_update)
    call broadcast(self%dump_response)
    call broadcast(self%field_local_allreduce)
    call broadcast(self%field_local_allreduce_sub)
    call broadcast(self%field_local_nonblocking_collectives)
    call broadcast(self%field_local_tuneminnrow)
    call broadcast(self%field_option)
    call broadcast(self%field_subgath)
    call broadcast(self%force_maxwell_reinit)
    call broadcast(self%minnrow)
    call broadcast(self%read_response)
    call broadcast(self%remove_zonal_flows_switch)
    call broadcast(self%response_dir)
    call broadcast(self%response_file)

    call broadcast(self%exist)
  end subroutine broadcast_fields_config

  !> Gets the default name for this namelist
  function get_default_name_fields_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_fields_config
    get_default_name_fields_config = "fields_knobs"
  end function get_default_name_fields_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_fields_config()
    implicit none
    logical :: get_default_requires_index_fields_config
    get_default_requires_index_fields_config = .false.
  end function get_default_requires_index_fields_config
end module fields
