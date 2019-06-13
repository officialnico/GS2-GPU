!> Define abstract type to use as a base type for derived types used to represent
!! the input configuration of a module (actual or logical). In other words it can
!! be used to represent a namelist state.
module abstract_config
  implicit none
  private

  !Types
  public :: abstract_config_type

  integer, parameter :: MAX_NAME_LEN = 64

  !> Define abstract type to represent the config type,
  !! allows common code to be moved here
  type, abstract :: abstract_config_type
     logical, private :: initialised = .false. !< Have we initialised the instance?
     logical :: exist = .false. !< Does the related namelist exist in the target input file?
     logical :: requires_index = .false. !< Is this a numbered namelist with an index
     integer :: index = 0 !< Used to hold the specific index of numbered namelists
     character(len = MAX_NAME_LEN) :: name="Undefined name." !< The name of the namelist that we represent.
   contains
     procedure :: is_initialised => is_initialised_generic
     procedure :: init => init_generic
     procedure :: setup => setup_generic
     procedure :: write_namelist_header
     procedure :: get_name => get_name_generic
     procedure, nopass :: write_namelist_footer
     procedure, private, nopass :: write_key_val_string
     procedure, private, nopass :: write_key_val_real
     procedure, private :: write_key_val_real_array     
     procedure, private, nopass :: write_key_val_complex
     procedure, private :: write_key_val_complex_array     
     procedure, private, nopass :: write_key_val_integer
     procedure, private :: write_key_val_integer_array
     procedure, private, nopass :: write_key_val_logical
     generic :: write_key_val => write_key_val_string, write_key_val_real, write_key_val_complex, write_key_val_integer, write_key_val_logical, write_key_val_real_array, write_key_val_complex_array, write_key_val_integer_array
     procedure(read_interface), deferred :: read
     procedure(write_interface), deferred :: write
     procedure(reset_interface), deferred :: reset
     procedure(broadcast_interface), deferred :: broadcast
  end type abstract_config_type

  interface
     subroutine read_interface(self)
       import
       class(abstract_config_type), intent(in out) :: self
     end subroutine read_interface

     subroutine write_interface(self, unit)
       import
       class(abstract_config_type), intent(in) :: self
       integer, intent(in), optional :: unit
     end subroutine write_interface

     subroutine reset_interface(self)
       import
       class(abstract_config_type), intent(in out) :: self
     end subroutine reset_interface

     subroutine broadcast_interface(self)
       import
       class(abstract_config_type), intent(in out) :: self
     end subroutine broadcast_interface
  end interface
contains

  !> Is this instance initialised?
  function is_initialised_generic(self)
    class(abstract_config_type), intent(in) :: self
    logical :: is_initialised_generic
    is_initialised_generic = self%initialised
  end function is_initialised_generic

  !> Fully initialise the config object
  subroutine init_generic(self, name, requires_index, index)
    implicit none
    class(abstract_config_type), intent(inout) :: self
    character(len = *), intent(in) :: name
    logical, intent(in) :: requires_index
    integer, intent(in), optional :: index
    if (self%is_initialised()) return
    call self%setup(name, requires_index, index)
    call self%read()
    call self%broadcast()
  end subroutine init_generic
  
  !> Do some standard setup/checking
  subroutine setup_generic(self, name, requires_index, index)
    implicit none
    class(abstract_config_type), intent(inout) :: self
    character(len = *), intent(in), optional :: name
    logical, intent(in), optional :: requires_index
    integer, intent(in), optional :: index
    if(present(name)) self%name = name
    if(present(requires_index)) self%requires_index = requires_index

    !Set a default index number if required and not set (could abort)
    if(present(index)) then
       self%index = index
    else
       if(self%requires_index) self%index = 1
    endif

    self%initialised = .true.
  end subroutine setup_generic

  !> Returns the namelist name. Not very useful at the moment
  !> but may want to do more interesting things in the future
  function get_name_generic(self)
    implicit none
    class(abstract_config_type), intent(in) :: self
    character(len = MAX_NAME_LEN) :: get_name_generic
    get_name_generic = self%name
  end function get_name_generic
  
  !> Write the namelist header for this instance
  subroutine write_namelist_header(self, unit)
    implicit none
    class(abstract_config_type), intent(in) :: self
    integer, intent(in) :: unit
    !Decide if we should include the index or not
    if(self%requires_index) then
       write(unit,'("&",A,"_",I0)') trim(self%get_name()), self%index
    else
       write(unit,'("&",A)') trim(self%get_name())
    endif
  end subroutine write_namelist_header

  !> Write the namelist footer
  subroutine write_namelist_footer(unit)
    implicit none
    integer, intent(in) :: unit
    write(unit,'("/")')
    write(unit,'()')
  end subroutine write_namelist_footer

  !> Writes a {key,val} pair where the value is of type character
  subroutine write_key_val_string(key, val, unit)
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    integer, intent(in) :: unit
    write(unit,'("  ",A20," = "," ",A)') key, '"'//trim(val)//'"'
  end subroutine write_key_val_string

  !> Writes a {key,val} pair where the value is of type real
  subroutine write_key_val_real(key, val, unit)
    character(len=*), intent(in) :: key
    real, intent(in) :: val
    integer, intent(in) :: unit
    write(unit,'("  ",A20," = ",e18.11)') key, val
  end subroutine write_key_val_real

  !> Writes a {key,val} pair where the value is of type real array
  subroutine write_key_val_real_array(self, key, val, unit)
    class(abstract_config_type), intent(in) :: self
    character(len = *), intent(in) :: key
    real, dimension(:), intent(in) :: val
    integer, intent(in) :: unit
    integer :: i
    character(len = 12) :: subscript
    do i = 1, size(val)
       write(subscript,'("(",I0,")")') i
       call self%write_key_val(trim(key)//subscript, val(i), unit)
    end do
  end subroutine write_key_val_real_array
  
  !> Writes a {key,val} pair where the value is of type complex
  subroutine write_key_val_complex(key, val, unit)
    character(len=*), intent(in) :: key
    complex, intent(in) :: val
    integer, intent(in) :: unit
    write(unit,'("  ",A20," = (",e18.11,", ",e18.11,")")') key, real(val), aimag(val)
  end subroutine write_key_val_complex

  !> Writes a {key,val} pair where the value is of type complex array
  subroutine write_key_val_complex_array(self, key, val, unit)
    class(abstract_config_type), intent(in) :: self
    character(len = *), intent(in) :: key
    complex, dimension(:), intent(in) :: val
    integer, intent(in) :: unit
    integer :: i
    character(len = 12) :: subscript
    do i = 1, size(val)
       write(subscript,'("(",I0,")")') i
       call self%write_key_val(trim(key)//subscript, val(i), unit)
    end do
  end subroutine write_key_val_complex_array
  
  !> Writes a {key,val} pair where the value is of type integer
  subroutine write_key_val_integer(key, val, unit)
    character(len=*), intent(in) :: key
    integer, intent(in) :: val
    integer, intent(in) :: unit
    if (val.lt.0) then
       write(unit,'("  ",A20," = ",I0)') key, val
    else
       write(unit,'("  ",A20," = "," ",I0)') key, val
    endif
  end subroutine write_key_val_integer

  !> Writes a {key,val} pair where the value is of type integer array
  subroutine write_key_val_integer_array(self, key, val, unit)
    class(abstract_config_type), intent(in) :: self
    character(len = *), intent(in) :: key
    integer, dimension(:), intent(in) :: val
    integer, intent(in) :: unit
    integer :: i
    character(len = 12) :: subscript
    do i = 1, size(val)
       write(subscript,'("(",I0,")")') i
       call self%write_key_val(trim(key)//subscript, val(i), unit)
    end do
  end subroutine write_key_val_integer_array
  
  !> Writes a {key,val} pair where the value is of type logical
  subroutine write_key_val_logical(key, val, unit)
    character(len=*), intent(in) :: key
    logical, intent(in) :: val
    integer, intent(in) :: unit
    write(unit,'("  ",A20," = "," ",L1)') key, val
  end subroutine write_key_val_logical
end module abstract_config
