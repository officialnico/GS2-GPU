
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the eigval_knobs namelist and populates the member variables
subroutine read_eigval_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(eigval_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: extraction_option, solver_option, transform_option, which_option
  integer :: max_iter, n_eig, nadv
  logical :: analyse_ddt_operator, save_restarts, use_ginit
  real :: targ_im, targ_re, tolerance

  namelist /eigval_knobs/ analyse_ddt_operator, extraction_option, max_iter, n_eig, nadv, save_restarts, solver_option, targ_im, targ_re, tolerance, transform_option, use_ginit, which_option

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  analyse_ddt_operator = self%analyse_ddt_operator
  extraction_option = self%extraction_option
  max_iter = self%max_iter
  n_eig = self%n_eig
  nadv = self%nadv
  save_restarts = self%save_restarts
  solver_option = self%solver_option
  targ_im = self%targ_im
  targ_re = self%targ_re
  tolerance = self%tolerance
  transform_option = self%transform_option
  use_ginit = self%use_ginit
  which_option = self%which_option

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = eigval_knobs)

  ! Now copy from local variables into type members
  self%analyse_ddt_operator = analyse_ddt_operator
  self%extraction_option = extraction_option
  self%max_iter = max_iter
  self%n_eig = n_eig
  self%nadv = nadv
  self%save_restarts = save_restarts
  self%solver_option = solver_option
  self%targ_im = targ_im
  self%targ_re = targ_re
  self%tolerance = tolerance
  self%transform_option = transform_option
  self%use_ginit = use_ginit
  self%which_option = which_option

  self%exist = exist
end subroutine read_eigval_config

!> Writes out a namelist representing the current state of the config object
subroutine write_eigval_config(self, unit)
  implicit none
  class(eigval_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("analyse_ddt_operator", self%analyse_ddt_operator, unit_internal)
  call self%write_key_val("extraction_option", self%extraction_option, unit_internal)
  call self%write_key_val("max_iter", self%max_iter, unit_internal)
  call self%write_key_val("n_eig", self%n_eig, unit_internal)
  call self%write_key_val("nadv", self%nadv, unit_internal)
  call self%write_key_val("save_restarts", self%save_restarts, unit_internal)
  call self%write_key_val("solver_option", self%solver_option, unit_internal)
  call self%write_key_val("targ_im", self%targ_im, unit_internal)
  call self%write_key_val("targ_re", self%targ_re, unit_internal)
  call self%write_key_val("tolerance", self%tolerance, unit_internal)
  call self%write_key_val("transform_option", self%transform_option, unit_internal)
  call self%write_key_val("use_ginit", self%use_ginit, unit_internal)
  call self%write_key_val("which_option", self%which_option, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_eigval_config

!> Resets the config object to the initial empty state
subroutine reset_eigval_config(self)
  class(eigval_config_type), intent(in out) :: self
  type(eigval_config_type) :: empty
  select type (self)
  type is (eigval_config_type)
     self = empty
  end select
end subroutine reset_eigval_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_eigval_config(self)
  use mp, only: broadcast
  implicit none
  class(eigval_config_type), intent(in out) :: self
  call broadcast(self%analyse_ddt_operator)
  call broadcast(self%extraction_option)
  call broadcast(self%max_iter)
  call broadcast(self%n_eig)
  call broadcast(self%nadv)
  call broadcast(self%save_restarts)
  call broadcast(self%solver_option)
  call broadcast(self%targ_im)
  call broadcast(self%targ_re)
  call broadcast(self%tolerance)
  call broadcast(self%transform_option)
  call broadcast(self%use_ginit)
  call broadcast(self%which_option)

  call broadcast(self%exist)
end subroutine broadcast_eigval_config

!> Gets the default name for this namelist
function get_default_name_eigval_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_eigval_config
  get_default_name_eigval_config = "eigval_knobs"
end function get_default_name_eigval_config

!> Gets the default requires index for this namelist
function get_default_requires_index_eigval_config()
  implicit none
  logical :: get_default_requires_index_eigval_config
  get_default_requires_index_eigval_config = .false.
end function get_default_requires_index_eigval_config

!> Get the module level config instance
pure function get_eigval_config()
  type(eigval_config_type) :: get_eigval_config
  get_eigval_config = eigval_config
end function get_eigval_config
