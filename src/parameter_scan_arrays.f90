!> Data module holding variables associated with the parameter
!> scan functionality. Held separately from parameter_scan
!> module to avoid complicating dependencies.
module parameter_scan_arrays
  implicit none
  
  private

  public :: current_scan_parameter_value, write_scan_parameter
  public :: hflux_tot, momflux_tot, phi2_tot, nout

  real :: current_scan_parameter_value
  logical :: write_scan_parameter
  real, dimension(:), allocatable :: hflux_tot, momflux_tot, phi2_tot
  integer :: nout
end module parameter_scan_arrays
