!> Types with bound routines for representing/holding the "input configuration"
!! set from the contents of the provided input file. This allows us to override
!! what's in the input file by changing values in instances of these types.
!!
!! @todo Allow a specific input file to be specified rather than relying on file_utils choice
!!
!! @todo Testing
module config_types
  use abstract_config, only: abstract_config_type
  implicit none
  private
  
!contains

end module config_types
