!! Write some standard metadata to a netCDF file
!!
!! Standalone module to avoid circular dependencies
module gs2_metadata
  implicit none
contains
  !> Add some standard metadata as global attributes to the file:
  !>
  !> - title: description of file
  !> - software_name: GS2
  !> - software_version: <current version>
  !> - netcdf_version: version of netCDF used
  !> - date_created: date/time of this simulation
  !> - id: UUID of this simulation
  !> - comment: human-friendly description of run
  !> - Conventions: URL of GS2 website
  !>
  !> Uses [[standard_header_type]] to get the UUID and date-time consistent
  !> across the simulation
  subroutine create_metadata(file_id, title, header, comments)
#ifdef NETCDF
    use netcdf, only : NF90_GLOBAL, nf90_inq_libvers
    use netcdf_utils, only: ensure_netcdf_var_exists, netcdf_error, ensure_netcdf_dim_exists, ensure_netcdf_att_exists
    use git_version_mod, only : get_git_version
#endif
    use standard_header, only : standard_header_type

    !> NetCDF ID of the file
    integer, intent(in) :: file_id
    !> Brief description of the file contents
    character(len=*), intent(in) :: title
    !> Pre-generated header with UUID and date information
    type(standard_header_type), optional, intent(in) :: header
    !> Comments on this file, for example a human-friendly description
    character(len=*), optional, intent(in) :: comments
#ifdef NETCDF
    type(standard_header_type) :: local_header

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "title", trim(title))
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "software_name", "GS2")
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "software_version", trim(get_git_version()))
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "netcdf_version", trim(nf90_inq_libvers()))
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "date_created", trim(local_header%date_time))
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "id", trim(local_header%run_uuid))
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "Conventions", "https://gyrokinetics.gitlab.io/gs2/")
    call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "references", "https://doi.org/10.5281/zenodo.2551066")
    if (present(comments)) then
      call ensure_netcdf_att_exists(file_id, NF90_GLOBAL, "comment", trim(comments))
    end if
#endif
  end subroutine create_metadata
end module gs2_metadata
