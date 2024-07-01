
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the species_parameters namelist and populates the member variables
subroutine read_species_element_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(species_element_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: f0type, type
  real :: bess_fac, dens, dens0, fprim, mass, nu_h, temp, tpar0, tperp0, tprim, u0, uprim, uprim2, vcprim, vcrit, vnewk, z

  namelist /species_parameters/ bess_fac, dens, dens0, f0type, fprim, mass, nu_h, temp, tpar0, tperp0, tprim, type, u0, uprim, uprim2, vcprim, vcrit, vnewk, z

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  bess_fac = self%bess_fac
  dens = self%dens
  dens0 = self%dens0
  f0type = self%f0type
  fprim = self%fprim
  mass = self%mass
  nu_h = self%nu_h
  temp = self%temp
  tpar0 = self%tpar0
  tperp0 = self%tperp0
  tprim = self%tprim
  type = self%type
  u0 = self%u0
  uprim = self%uprim
  uprim2 = self%uprim2
  vcprim = self%vcprim
  vcrit = self%vcrit
  vnewk = self%vnewk
  z = self%z

  ! Now read in the main namelist
  call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index, exist)
  if (exist) read(in_file, nml = species_parameters)
  close(unit = in_file)

  ! Now copy from local variables into type members
  self%bess_fac = bess_fac
  self%dens = dens
  self%dens0 = dens0
  self%f0type = f0type
  self%fprim = fprim
  self%mass = mass
  self%nu_h = nu_h
  self%temp = temp
  self%tpar0 = tpar0
  self%tperp0 = tperp0
  self%tprim = tprim
  self%type = type
  self%u0 = u0
  self%uprim = uprim
  self%uprim2 = uprim2
  self%vcprim = vcprim
  self%vcrit = vcrit
  self%vnewk = vnewk
  self%z = z

  self%exist = exist
end subroutine read_species_element_config

!> Writes out a namelist representing the current state of the config object
subroutine write_species_element_config(self, unit)
  implicit none
  class(species_element_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("bess_fac", self%bess_fac, unit_internal)
  call self%write_key_val("dens", self%dens, unit_internal)
  call self%write_key_val("dens0", self%dens0, unit_internal)
  call self%write_key_val("f0type", self%f0type, unit_internal)
  call self%write_key_val("fprim", self%fprim, unit_internal)
  call self%write_key_val("mass", self%mass, unit_internal)
  call self%write_key_val("nu_h", self%nu_h, unit_internal)
  call self%write_key_val("temp", self%temp, unit_internal)
  call self%write_key_val("tpar0", self%tpar0, unit_internal)
  call self%write_key_val("tperp0", self%tperp0, unit_internal)
  call self%write_key_val("tprim", self%tprim, unit_internal)
  call self%write_key_val("type", self%type, unit_internal)
  call self%write_key_val("u0", self%u0, unit_internal)
  call self%write_key_val("uprim", self%uprim, unit_internal)
  call self%write_key_val("uprim2", self%uprim2, unit_internal)
  call self%write_key_val("vcprim", self%vcprim, unit_internal)
  call self%write_key_val("vcrit", self%vcrit, unit_internal)
  call self%write_key_val("vnewk", self%vnewk, unit_internal)
  call self%write_key_val("z", self%z, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_species_element_config

!> Resets the config object to the initial empty state
subroutine reset_species_element_config(self)
  class(species_element_config_type), intent(in out) :: self
  type(species_element_config_type) :: empty
  select type (self)
  type is (species_element_config_type)
     self = empty
  end select
end subroutine reset_species_element_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_species_element_config(self)
  use mp, only: broadcast
  implicit none
  class(species_element_config_type), intent(in out) :: self
  call broadcast(self%bess_fac)
  call broadcast(self%dens)
  call broadcast(self%dens0)
  call broadcast(self%f0type)
  call broadcast(self%fprim)
  call broadcast(self%mass)
  call broadcast(self%nu_h)
  call broadcast(self%temp)
  call broadcast(self%tpar0)
  call broadcast(self%tperp0)
  call broadcast(self%tprim)
  call broadcast(self%type)
  call broadcast(self%u0)
  call broadcast(self%uprim)
  call broadcast(self%uprim2)
  call broadcast(self%vcprim)
  call broadcast(self%vcrit)
  call broadcast(self%vnewk)
  call broadcast(self%z)

  call broadcast(self%exist)
end subroutine broadcast_species_element_config

!> Gets the default name for this namelist
function get_default_name_species_element_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_species_element_config
  get_default_name_species_element_config = "species_parameters"
end function get_default_name_species_element_config

!> Gets the default requires index for this namelist
function get_default_requires_index_species_element_config()
  implicit none
  logical :: get_default_requires_index_species_element_config
  get_default_requires_index_species_element_config = .true.
end function get_default_requires_index_species_element_config

!> Get the array of module level config instances. If it isn't allocated,
!> then return a zero-length array
pure function get_species_element_config()
  type(species_element_config_type), allocatable, dimension(:) :: get_species_element_config
  if (allocated(species_element_config)) then
    get_species_element_config = species_element_config
  else
    allocate(get_species_element_config(0))
  end if
end function get_species_element_config
