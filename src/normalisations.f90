!> FIXME : Add documentation
module normalisations
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none
  
  private

  public :: norms
  public :: init_normalisations, finish_normalisations
  public :: check_normalisations, wnml_normalisations

  public :: normalisations_config_type
  public :: set_normalisations_config
  
  real, parameter :: default_value = -999.0
  
  !> Define a type to hold the set of normalisations, we
  !> could imagine extending this in the future to provide
  !> a normalise/denormalise routine which can be used to
  !> automatically convert quantities
  type norms_type
     private
     !Internal parameters
     real :: def_val = default_value !Default value, used to check if parameter has been set.

     !Note def_val is really a parameter, but not allowed to define params in derived type
     logical :: some_norms_set=.false.
     logical :: all_norms_set=.false.
     logical :: initialised = .false.

     !The following is used to help automate loops over normalisations
     integer :: nnorm
     character(len=6), dimension(:), allocatable, public :: names

     !The normalisations
     real :: mref !< Reference mass in atomic mass units
     real :: zref !< Reference charge in units of the elementary charge
     real :: nref !< Reference density in $m^3$
     real :: tref !< Reference temperature in $eV$
     real :: aref !< Reference length in $m$
     real :: vref !< Reference (thermal) velocity in $m/s$
     real :: bref !< Reference magnetic field in Tesla, $T$
     real :: rhoref !< Reference Larmor radius in $m$
   contains
     private
     procedure, public :: init => norms_init
     procedure, public :: finish => norms_finish
     procedure, public :: reset => norms_reset
     procedure :: read_parameters => norms_read_parameters
     procedure :: set_value => norms_set_value
     procedure, public :: get_value => norms_get_value
     procedure :: check_if_set => norms_check_if_set
     procedure :: set_logicals => norms_set_logicals
  end type norms_type

  type(norms_type) :: norms !The normalisation object, there should only be one of these
  logical :: exist !Has this namelist been specified?
  logical :: initialized = .false. !Have we setup the module yet?

  !> Used to represent the input configuration of normalisations
  type, extends(abstract_config_type) :: normalisations_config_type
     ! namelist : normalisations_knobs
     !> FIXME: Add documentation
     real :: aref = default_value
     !> FIXME: Add documentation
     real :: bref = default_value
     !> FIXME: Add documentation
     real :: mref = default_value
     !> FIXME: Add documentation
     real :: nref = default_value
     !> FIXME: Add documentation
     real :: rhoref = default_value
     !> FIXME: Add documentation
     real :: tref = default_value
     !> FIXME: Add documentation
     real :: vref = default_value
     !> FIXME: Add documentation
     real :: zref = default_value
   contains
     procedure, public :: read => read_normalisations_config
     procedure, public :: write => write_normalisations_config
     procedure, public :: reset => reset_normalisations_config
     procedure, public :: broadcast => broadcast_normalisations_config
     procedure, public, nopass :: get_default_name => get_default_name_normalisations_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_normalisations_config
  end type normalisations_config_type
  
  type(normalisations_config_type) :: normalisations_config  
contains
  !/////////////////////////////
  !// TYPE BOUND PROCEDURES
  !/////////////////////////////

  !> Sets the value of a particular normalisation, determined by the
  !> passed `val_name` string. Unknown `val_name` values result in a
  !> warning message on `proc0`.
  subroutine norms_set_value(self, val_name, val)
    use runtime_tests, only: verbosity
    use mp, only: proc0
    implicit none
    class(norms_type), intent(in out) :: self
    character(len=*), intent(in) :: val_name
    real, intent(in) :: val
    
    !Here we can have a message in very verbose runs
    if((verbosity().gt.5).and.proc0)then
       if(val.eq.self%def_val)then
          write(6,'("The ",A," normalisation has not been set. This may prevent conversion of some quantities.")') trim(val_name)
       else
          write(6,'("Setting the ",A," normalisation to ",F12.5)') trim(val_name), val
       endif
    endif

    !Should probably convert to lower case here, but until we
    !add some string utils to do this sort of stuff we'll rely
    !on developer doing the right thing.
    select case (trim(val_name))
    case("mref")
       self%mref=val
    case("zref")
       self%zref=val
    case("nref")
       self%nref=val
    case("tref")
       self%tref=val
    case("aref")
       self%aref=val
    case("vref")
       self%vref=val
    case("bref")
       self%bref=val
    case("rhoref")
       self%rhoref=val
    case default
       if(proc0) write(6,'("Warning : Attempt to set unknown normalisation ",A," --> Ignoring")')
    end select
  end subroutine norms_set_value

  !> Get the value of the normalisation associated with `val_name`.
  !> Unknown values of `val_name` result in a call to `mp_abort`.
  function norms_get_value(self, val_name)
    use mp, only: mp_abort
    implicit none
    class(norms_type), intent(in) :: self
    character(len=*), intent(in) :: val_name
    real :: norms_get_value

    !Should probably convert to lower case here, but until we
    !add some string utils to do this sort of stuff we'll rely
    !on developer doing the right thing.
    select case (trim(val_name))
    case("mref")
       norms_get_value=self%mref
    case("zref")
       norms_get_value=self%zref
    case("nref")
       norms_get_value=self%nref
    case("tref")
       norms_get_value=self%tref
    case("aref")
       norms_get_value=self%aref
    case("vref")
       norms_get_value=self%vref
    case("bref")
       norms_get_value=self%bref
    case("rhoref")
       norms_get_value=self%rhoref
    case default
      call mp_abort("Invalid normalisation requested")
    end select
  end function norms_get_value

  !> FIXME : Add documentation  
  subroutine norms_read_parameters(self, normalisations_config_in)
    use file_utils, only: input_unit_exist
    implicit none
    class(norms_type), intent(in out) :: self
    type(normalisations_config_type), intent(in), optional :: normalisations_config_in
    real :: mref,zref,nref,tref,aref,vref,bref,rhoref

    if (present(normalisations_config_in)) normalisations_config = normalisations_config_in

    call normalisations_config%init(name = 'normalisations_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    aref = normalisations_config%aref
    bref = normalisations_config%bref
    mref = normalisations_config%mref
    nref = normalisations_config%nref
    rhoref = normalisations_config%rhoref
    tref = normalisations_config%tref
    vref = normalisations_config%vref
    zref = normalisations_config%zref

    exist = normalisations_config%exist

    !Now copy parameters into holder type
    call self%set_value("mref",mref)
    call self%set_value("zref",zref)
    call self%set_value("nref",nref)
    call self%set_value("tref",tref)
    call self%set_value("aref",aref)
    call self%set_value("vref",vref)
    call self%set_value("bref",bref)
    call self%set_value("rhoref",rhoref)
  end subroutine norms_read_parameters

  !> Initialise the norms object
  subroutine norms_init(self, normalisations_config_in)
    implicit none
    class(norms_type), intent(in out) :: self
    type(normalisations_config_type), intent(in), optional :: normalisations_config_in

    !First setup allowed normalisations
    self%nnorm=8
    if(.not.allocated(self%names)) allocate(self%names(self%nnorm))
    self%names(1)="mref"
    self%names(2)="zref"
    self%names(3)="nref"
    self%names(4)="tref"
    self%names(5)="aref"
    self%names(6)="vref"
    self%names(7)="bref"
    self%names(8)="rhoref"

    if(self%initialised) return
    self%initialised = .true.

    call self%read_parameters(normalisations_config_in)
  end subroutine norms_init

  !> Reset the properties
  subroutine norms_reset(self)
    implicit none
    class(norms_type), intent(in out) :: self
    integer :: i
    
    !Loop over parameters and set them to def_val
    do i=1,len(self%names)
       call self%set_value(self%names(i),self%def_val)
    enddo

    !Set the logical vars
    self%initialised=.false.
    call self%set_logicals
  end subroutine norms_reset

  !> Reset and free memory
  subroutine norms_finish(self)
    implicit none
    class(norms_type), intent(in out) :: self
    call self%reset
    if(allocated(self%names)) deallocate(self%names)
  end subroutine norms_finish

  !> Decide if a given normalisation has been set
  function norms_check_if_set(self,val_name)
    implicit none
    class(norms_type), intent(in) :: self
    character(len=*), intent(in) :: val_name
    logical :: norms_check_if_set
    norms_check_if_set = self%get_value(val_name).ne.self%def_val
  end function norms_check_if_set

  !> Decide if all/some of the normalisations have been set
  subroutine norms_set_logicals(self)
    implicit none
    class(norms_type), intent(in out) :: self
    integer :: i
    logical :: some_set, all_set

    !Init internals
    some_set=.false.
    all_set=.true.

    !Loop over parameters and set them to def_val
    do i=1,len(self%names)
       all_set=all_set.and.self%check_if_set(self%names(i))
       some_set=some_set.or.self%check_if_set(self%names(i))
    enddo

    !Update object parameters
    self%some_norms_set=some_set
    self%all_norms_set=all_set
  end subroutine norms_set_logicals

  !/////////////////////////////
  !// MODULE LEVEL PROCDURES
  !/////////////////////////////
  
  !> FIXME : Add documentation
  subroutine check_normalisations(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    character(len=7) :: msg
    if(norms%all_norms_set)then
       msg = "All of"
    else
       if(norms%some_norms_set)then
          msg = "Some of"
       else
          msg = "None of"
       endif
    endif

    write(report_unit,fmt='(A," the normalisation parameters have been provided.")') trim(msg)
  end subroutine check_normalisations

  !> FIXME : Add documentation  
  subroutine wnml_normalisations(unit)
    implicit none
    integer, intent(in) :: unit
    integer :: i
    if(.not.exist) return

    write(unit,*)
    do i=1,len(norms%names)
       write(unit,fmt='(A," = ", F12.5)') trim(norms%names(i)),norms%get_value(norms%names(i))
    enddo
    write(unit,'("/")')
  end subroutine wnml_normalisations

  !> Read input file and populate the norms object
  subroutine init_normalisations(normalisations_config_in)
    implicit none
    type(normalisations_config_type), intent(in), optional :: normalisations_config_in
    if(initialized) return
    initialized = .true.
    call norms%init(normalisations_config_in)
  end subroutine init_normalisations

  !> Free memory etc. associated with normalisations
  subroutine finish_normalisations
    implicit none
    initialized = .false.
    call norms%finish
    call normalisations_config%reset()
  end subroutine finish_normalisations

  !> Reset normalisations object
  subroutine reset_normalisations
    implicit none
    call norms%reset    
  end subroutine reset_normalisations

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_normalisations_config(normalisations_config_in)
    use mp, only: mp_abort
    type(normalisations_config_type), intent(in), optional :: normalisations_config_in
    if (initialized) then
       call mp_abort("Trying to set normalisations_config when already initialized.", to_screen = .true.)
    end if
    if (present(normalisations_config_in)) then
       normalisations_config = normalisations_config_in
    end if
  end subroutine set_normalisations_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the normalisations_knobs namelist and populates the member variables
  subroutine read_normalisations_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(normalisations_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    real :: aref, bref, mref, nref, rhoref, tref, vref, zref

    namelist /normalisations_knobs/ aref, bref, mref, nref, rhoref, tref, vref, zref

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    aref = self%aref
    bref = self%bref
    mref = self%mref
    nref = self%nref
    rhoref = self%rhoref
    tref = self%tref
    vref = self%vref
    zref = self%zref

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = normalisations_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = normalisations_knobs)
    endif

    ! Now copy from local variables into type members
    self%aref = aref
    self%bref = bref
    self%mref = mref
    self%nref = nref
    self%rhoref = rhoref
    self%tref = tref
    self%vref = vref
    self%zref = zref

    self%exist = exist
  end subroutine read_normalisations_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_normalisations_config(self, unit)
    implicit none
    class(normalisations_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("aref", self%aref, unit_internal)
    call self%write_key_val("bref", self%bref, unit_internal)
    call self%write_key_val("mref", self%mref, unit_internal)
    call self%write_key_val("nref", self%nref, unit_internal)
    call self%write_key_val("rhoref", self%rhoref, unit_internal)
    call self%write_key_val("tref", self%tref, unit_internal)
    call self%write_key_val("vref", self%vref, unit_internal)
    call self%write_key_val("zref", self%zref, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_normalisations_config

  !> Resets the config object to the initial empty state
  subroutine reset_normalisations_config(self)
    class(normalisations_config_type), intent(in out) :: self
    type(normalisations_config_type) :: empty
    select type (self)
    type is (normalisations_config_type)
       self = empty
    end select
  end subroutine reset_normalisations_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_normalisations_config(self)
    use mp, only: broadcast
    implicit none
    class(normalisations_config_type), intent(in out) :: self
    call broadcast(self%aref)
    call broadcast(self%bref)
    call broadcast(self%mref)
    call broadcast(self%nref)
    call broadcast(self%rhoref)
    call broadcast(self%tref)
    call broadcast(self%vref)
    call broadcast(self%zref)

    call broadcast(self%exist)
  end subroutine broadcast_normalisations_config

  !> Gets the default name for this namelist
  function get_default_name_normalisations_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_normalisations_config
    get_default_name_normalisations_config = "normalisations_knobs"
  end function get_default_name_normalisations_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_normalisations_config()
    implicit none
    logical :: get_default_requires_index_normalisations_config
    get_default_requires_index_normalisations_config = .false.
  end function get_default_requires_index_normalisations_config
end module normalisations
