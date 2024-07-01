
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the ballstab_knobs namelist and populates the member variables
subroutine read_ballstab_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(ballstab_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: n_beta, n_shat
  logical :: make_salpha
  real :: beta_div, beta_mul, diff, shat_max, shat_min, theta0

  namelist /ballstab_knobs/ beta_div, beta_mul, diff, make_salpha, n_beta, n_shat, shat_max, shat_min, theta0

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  beta_div = self%beta_div
  beta_mul = self%beta_mul
  diff = self%diff
  make_salpha = self%make_salpha
  n_beta = self%n_beta
  n_shat = self%n_shat
  shat_max = self%shat_max
  shat_min = self%shat_min
  theta0 = self%theta0

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = ballstab_knobs)

  ! Now copy from local variables into type members
  self%beta_div = beta_div
  self%beta_mul = beta_mul
  self%diff = diff
  self%make_salpha = make_salpha
  self%n_beta = n_beta
  self%n_shat = n_shat
  self%shat_max = shat_max
  self%shat_min = shat_min
  self%theta0 = theta0

  self%exist = exist
end subroutine read_ballstab_config

!> Writes out a namelist representing the current state of the config object
subroutine write_ballstab_config(self, unit)
  implicit none
  class(ballstab_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("beta_div", self%beta_div, unit_internal)
  call self%write_key_val("beta_mul", self%beta_mul, unit_internal)
  call self%write_key_val("diff", self%diff, unit_internal)
  call self%write_key_val("make_salpha", self%make_salpha, unit_internal)
  call self%write_key_val("n_beta", self%n_beta, unit_internal)
  call self%write_key_val("n_shat", self%n_shat, unit_internal)
  call self%write_key_val("shat_max", self%shat_max, unit_internal)
  call self%write_key_val("shat_min", self%shat_min, unit_internal)
  call self%write_key_val("theta0", self%theta0, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_ballstab_config

!> Resets the config object to the initial empty state
subroutine reset_ballstab_config(self)
  class(ballstab_config_type), intent(in out) :: self
  type(ballstab_config_type) :: empty
  select type (self)
  type is (ballstab_config_type)
     self = empty
  end select
end subroutine reset_ballstab_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_ballstab_config(self)
  use mp, only: broadcast
  implicit none
  class(ballstab_config_type), intent(in out) :: self
  call broadcast(self%beta_div)
  call broadcast(self%beta_mul)
  call broadcast(self%diff)
  call broadcast(self%make_salpha)
  call broadcast(self%n_beta)
  call broadcast(self%n_shat)
  call broadcast(self%shat_max)
  call broadcast(self%shat_min)
  call broadcast(self%theta0)

  call broadcast(self%exist)
end subroutine broadcast_ballstab_config

!> Gets the default name for this namelist
function get_default_name_ballstab_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_ballstab_config
  get_default_name_ballstab_config = "ballstab_knobs"
end function get_default_name_ballstab_config

!> Gets the default requires index for this namelist
function get_default_requires_index_ballstab_config()
  implicit none
  logical :: get_default_requires_index_ballstab_config
  get_default_requires_index_ballstab_config = .false.
end function get_default_requires_index_ballstab_config

!> Get the module level config instance
pure function get_ballstab_config()
  type(ballstab_config_type) :: get_ballstab_config
  get_ballstab_config = ballstab_config
end function get_ballstab_config
