!> Data module holding variables associated with the parameter
!> scan functionality. Held separately from parameter_scan
!> module to avoid complicating dependencies.
module parameter_scan_arrays
  implicit none
  
  private

  public :: current_scan_parameter_value, write_scan_parameter
  public :: scan_parameter_switch, scan_spec
  public :: scan_parameter_tprim, scan_parameter_g_exb
  public :: hflux_tot, momflux_tot, phi2_tot
  public :: growth_rate, nout

  real :: current_scan_parameter_value
  logical :: write_scan_parameter
  integer :: scan_parameter_switch
  integer, parameter :: scan_parameter_tprim = 1, scan_parameter_g_exb = 2
  integer :: scan_spec, nout
  integer :: nout
  real, dimension(:), allocatable :: hflux_tot, momflux_tot, phi2_tot
  real, dimension(:), allocatable :: growth_rate 
end module parameter_scan_arrays
