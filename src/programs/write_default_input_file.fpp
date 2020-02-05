program write_default_input_file
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: get_unused_unit, open_output_file, close_output_file

  use config_collection, only: config_collection_type
  implicit none

  type(config_collection_type) :: the_config

  logical :: list
  integer :: output_unit
  integer :: i

  call init_mp

  if (proc0) then
     call get_unused_unit(output_unit)
     open(unit = output_unit, action = "write", file = "defaults.in")
     call the_config%write_to_unit(output_unit)
     call close_output_file(output_unit)
  end if

  call finish_mp
end program write_default_input_file
