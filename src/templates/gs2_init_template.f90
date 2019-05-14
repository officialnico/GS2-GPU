! This file is used by gs2_init.rb as the template for gs2_init.f90.
! gs2_init.rb takes this template and replaces certain sections of it
! with automatically generated code. 
!
! Do not modify any lines containing the word TEMPLATE
! TEMPLATE_END_HEADER



!> This module is analogous to the init() function
!! in Linux-based operating systems: it initialises
!! gs2 to a certain init_level. At a given init level,
!! certain modules are initialised and certain are not.
!!
!! The gs2_init module is used by gs2_main to initialise modules. A typical
!! additional use case for this module is when it is desired
!! to override a given parameter (as in the override_* functions
!! in gs2_main). Gs2 must be taken down to the appropriate 
!! init_level, where all modules which contain any of those
!! parameters are uninitialized. The override is then set
!! and gs2 is brought back up to the highest init_level.
!! 
!! As in Linux, this module cannot be used until a certain
!! basic initialization has happened (think loading the kernel).
!! This basic initialization occurs in gs2_initialize in gs2_main,
!! and set the init_level to gs2_initialized.
!!
!! This is free software released under the MIT licence.
!! Written by:
!!            Edmund Highcock (edmundhighcock@users.sourceforge.net)
module gs2_init
  use abstract_config, only: abstract_config_type

  use overrides, only: miller_geometry_overrides_type
  use overrides, only: profiles_overrides_type
  use overrides, only: kt_grids_overrides_type
  use overrides, only: timestep_overrides_type
  use overrides, only: initial_values_overrides_type
  use overrides, only: optimisations_overrides_type
  public :: init_type
  !> A list of possible intialization levels.
  public :: init_level_list
  !> Bring gs2 to the target initialization
  !! level.
  public :: init
  !> Reads the gs2_init namelist
  public :: init_gs2_init
  !> Finalize the module
  public :: finish_gs2_init
  !> Save the current state of the fields and 
  !! distribution function, either to a file 
  !! or to a temporary array, depending on the
  !! value of in_memory. (NB if there is 
  !! sufficient memory to allocate a temporary
  !! array, in_memory will be overriden to 
  !! false).
  !public :: save_fields_and_dist_fn

  !> If initval_ov%in_memory is set to true (the default)
  !! write the current values of the dist fn into the 
  !! initval_ov object. If initval_ov%force_maxwell_reinit is
  !! set to false, also write the current field values into
  !! the object. 
  !!
  !! If initval_ov%in_memory is false, this will cause the current
  !! field and dist fn values to be written to file.
  !!
  !! You must have initialized the initval_ov
  !! object using init_initial_values_overrides before calling
  !! this function.
  !public :: set_initial_values_overrides_to_current_values

  public :: in_memory

  public :: init_config_type
  
  !> A type for labelling the different init
  !! levels available in gs2.
  type init_level_list_type
    !> The init_level reaches basic
    !! when initialize_gs2 has been called.
    ! TEMPLATE_GS2_INIT_LEVEL
    ! TEMPLATE_LEVEL_DECLARATIONS
  end type init_level_list_type

  type(init_level_list_type) :: init_level_list

  !> A type for storing the current initialization
  !! status, as well as all the overrides.
  type init_type
    !> The current init level
    integer :: level = 0
    !> Whether or not diagnostics have been initialized
    logical :: diagnostics_initialized = .false.
    !> An object for overriding all or selected
    !! Miller geometry parameters. You must call
    !! gs2_main::prepare_miller_geometry_overrides 
    !! before setting these overrides. See 
    !! documentation for the overrides::miller_geometry_overrides_type
    !! for more information.
    type(miller_geometry_overrides_type) :: mgeo_ov
    !> An object for overriding all or selected
    !! profile parameters such as species temperature, density, and gradients
    !! as well as the flow gradient and mach number. You must call
    !! gs2_main::prepare_profiles_overrides 
    !! before setting these overrides. See 
    !! documentation for the overrides::profiles_overrides_type
    !! for more information.
    type(kt_grids_overrides_type) :: kt_ov
    !> An object for overriding all or selected
    !! kt_grids parameters such as ny, nx, y0, jtwist etc.
    !! You must call
    !! gs2_main::prepare_kt_grids_overrides 
    !! before setting these overrides. See 
    !! documentation for the overrides::kt_grids_overrides_type
    !! for more information.
    type(profiles_overrides_type) :: prof_ov
    !> An object for overriding parameters connected
    !! to the timestep and cfl condition
    type(timestep_overrides_type) :: tstep_ov
    !> An object for overriding the initial values of 
    !! the fields and distribution function. You must call
    !! gs2_main::prepare_initial_values_overrides 
    !! before setting these overrides. This override
    !! is very complicated. See 
    !! documentation for the overrides::initial_values_overrides_type
    !! for more information.
    type(initial_values_overrides_type) :: initval_ov

    !> An object for overriding non physics parameters which
    !! may alter run time and efficiency. You must call
    !! gs2_main::prepare_optimisations_overrides 
    !! before setting these overrides. 
    type(optimisations_overrides_type) :: opt_ov


    !> A list of possible init levels
    !type(init_level_list_type) :: levels
  end type init_type

  private

  !complex, dimension(:,:,:), allocatable :: phi_tmp, apar_tmp, bpar_tmp
  !logical :: fields_and_dist_fn_saved = .false.
  logical :: in_memory = .false.

  
  !> Used to represent the input configuration of init
  type, extends(abstract_config_type) :: init_config_type
     ! namelist : init_knobs
     ! indexed : false    
     !> FIXME: Add documentation
     logical :: in_memory = .false.
   contains
     procedure, public :: read => read_init_config
     procedure, public :: write => write_init_config
     procedure, public :: reset => reset_init_config
     procedure, public :: broadcast => broadcast_init_config
  end type init_config_type

  interface init_config_type
     procedure :: init_config_constructor
  end interface init_config_type
  
  type(init_config_type) :: init_config  
  
contains
  !> Initialize gs2 to the level of target_level.
  !! The init_type current contains info
  !! about the current initialization level. At the end
  !! of the subroutine, current%level is set to target_level
  subroutine init(current, target_level)
    use fields, only: init_fields
    implicit none
    type(init_type), intent(inout) :: current
    integer, intent(in) :: target_level

    ! TEMPLATE_CHECK_GS2_INITIALISED

    if (current%level .eq. target_level) then
      return
    else
      if (up()) then 
        ! TEMPLATE_UP
      else if (down()) then
        ! TEMPLATE_DOWN
      end if
    end if
  contains
    function up()
      logical :: up
      up = (target_level>current%level) 
    end function up
    function down()
      logical :: down
      down = (target_level<current%level) 
    end function down
    ! TEMPLATE_SUBROUTINES


  end subroutine init

  subroutine init_gs2_init(current, init_config_in)
    implicit none
    type(init_type), intent(inout) :: current
    type(init_config_type), intent(in), optional :: init_config_in

    logical :: exist

    if (present(init_config_in)) init_config = init_config_in

    call init_config%init(name = 'init_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    in_memory = init_config%in_memory

    exist = init_config%exist
  end subroutine init_gs2_init

  subroutine finish_gs2_init(current)
    type(init_type), intent(inout) :: current
    call init_config%reset()
  end subroutine finish_gs2_init

  subroutine set_initial_field_and_dist_fn_values(current)
    use dist_fn_arrays, only: g, gnew
    use fields_arrays, only: phinew, aparnew, bparnew, phi, apar, bpar
    use fields, only: force_maxwell_reinit
    use fields, only: set_init_fields
    use file_utils, only: error_unit
    use init_g, only: ginit
    use init_g, only: ginitopt_restart_many
    use run_parameters, only: fphi, fapar, fbpar
    use mp, only: proc0
    use unit_tests, only: job_id
    implicit none
    type (init_type), intent(in) :: current
    logical :: restarted

    !write (*,*) 'set_init_field in_memory', in_memory, fields_and_dist_fn_saved

    if (.not. current%initval_ov%override) then
      ! This is the usual initial setup
      call ginit (restarted)
    else
      if (current%initval_ov%in_memory) then
        g = current%initval_ov%g
        gnew = g
      else
        call ginit(restarted, ginitopt_restart_many)
      end if
    end if

    if (current%initval_ov%force_maxwell_reinit)then
      call set_init_fields
      !if (.not. proc0) write (*,*) 'field value jjjj kkkk', phinew(-5, 3, 2), job_id
      !if (.not. proc0) write (*,*) 'field value sum jjjj kkkk', sum(real(conjg(phinew)*phinew)), job_id
    else

       if(current%initval_ov%override .and. current%initval_ov%in_memory) then
          if(fphi.gt.0) then
             phinew=current%initval_ov%phi
          else
             phinew = 0.
          endif
          if(fapar.gt.0) then
             aparnew=current%initval_ov%apar
          else
             aparnew = 0.
          endif
          if(fbpar.gt.0) then
             bparnew=current%initval_ov%bpar
          else
             bparnew = 0.
          endif
          phi = phinew
          apar = aparnew
          bpar = bparnew
       else
          ! No need to do anything: fields read from file.
       end if
    end if
  end subroutine set_initial_field_and_dist_fn_values

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the init_knobs namelist and populates the member variables
  subroutine read_init_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(init_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: in_memory

    namelist /init_knobs/ in_memory

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    in_memory = self%in_memory

    ! Now read in the main namelist
    if (self%requiresIndex) then
       call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
       read(in_file, nml = init_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%name, exist)
       if (exist) read(in_file, nml = init_knobs)
    endif

    ! Now copy from local variables into type members
    self%in_memory = in_memory

    self%exist = exist
  end subroutine read_init_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_init_config(self, unit)
    implicit none
    class(init_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("in_memory", self%in_memory, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_init_config

  !> Resets the config object to the initial empty state
  subroutine reset_init_config(self)
    class(init_config_type), intent(in out) :: self
    type(init_config_type) :: empty
    select type (self)
    type is (init_config_type)
       self = empty
    end select
  end subroutine reset_init_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_init_config(self)
    use mp, only: broadcast
    implicit none
    class(init_config_type), intent(in out) :: self
    call broadcast(self%in_memory)

    call broadcast(self%exist)
  end subroutine broadcast_init_config

  !> Sets the name and requiresIndex paremeters as appropriate
  function init_config_constructor()
    implicit none
    type(init_config_type) :: init_config_constructor
    init_config_constructor%name = "init_knobs"
    init_config_constructor%requiresIndex = .false.
  end function init_config_constructor
end module gs2_init

