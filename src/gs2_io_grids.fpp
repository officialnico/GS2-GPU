module gs2_io_grids

# ifdef NETCDF
  use netcdf_utils, only: kind_nf
  use neasyf, only: neasyf_open, neasyf_close
  use neasyf, only: neasyf_read, neasyf_write
  use neasyf, only: neasyf_dim
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: NF90_NOERR
# else
  use mp, only: mp_abort
# endif

  implicit none

  private

  public :: nc_get_grid_sizes, nc_write_grid_sizes
  public :: nc_get_lambda_grid_size
  public :: nc_get_grid_scalars, nc_write_grid_scalars
  public :: nc_get_grids, nc_write_grids
  public :: nc_get_lambda_grid, nc_write_lambda_grid
  public :: nc_grid_file_new

contains

  subroutine nc_grid_file_new(filename)
    use file_utils, only: get_unused_unit
    character (len=2000), intent(in) :: filename
    logical :: is_exist
    integer :: unit
    inquire(file=filename, exist=is_exist)
    if (is_exist) then
       write(6,*) 'WARNING: existing file will be over-written', trim(filename)
       open(unit, file=filename)
       close(unit, status='delete')
    end if
  end subroutine nc_grid_file_new

  subroutine nc_get_grid_sizes(filename, ntheta, ntgrid, nperiod)
    character (len=2000), intent(in) :: filename
    integer, intent(out), optional :: ntheta, ntgrid, nperiod

# ifdef NETCDF
    integer (kind_nf) :: ncid

    ncid = neasyf_open(trim(filename), "r")

    if (present(ntgrid)) call neasyf_read(ncid, "ntgrid", ntgrid)
    if (present(nperiod)) call neasyf_read(ncid, "nperiod", nperiod)
    if (present(ntheta)) call neasyf_read(ncid, "ntheta", ntheta)

    call neasyf_close(ncid)

# else
    call mp_abort("nc_get_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_get_grid_sizes

  subroutine nc_write_grid_sizes(filename, ntheta, ntgrid, nperiod)
    character (len=2000), intent(in) :: filename
    integer, intent(in), optional :: ntheta, ntgrid, nperiod
# ifdef NETCDF
    integer (kind_nf) :: ncid
    character (len=2) :: action
    logical :: is_exist

    ! always append mode if file exists
    action="w"
    inquire(file=trim(filename), exist=is_exist)
    if (is_exist) then
       action="rw"
    end if

    ncid = neasyf_open(trim(filename), action)

    if (present(ntgrid)) call neasyf_write(ncid, "ntgrid", ntgrid)
    if (present(nperiod)) call neasyf_write(ncid, "nperiod", nperiod)
    if (present(ntheta)) call neasyf_write(ncid, "ntheta", ntheta)

    call neasyf_close(ncid)

# else
    call mp_abort("nc_write_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_write_grid_sizes

  subroutine nc_get_lambda_grid_size(filename, nlambda)
    character (len=2000), intent(in) :: filename
    integer, intent(out) :: nlambda
# ifdef NETCDF
    integer (kind_nf) :: ncid, dimid
    integer :: status
    ncid = neasyf_open(trim(filename), "r")
    status = nf90_inq_dimid(ncid, "nlambda", dimid)
    if (status == NF90_NOERR) then
       status = nf90_inquire_dimension(ncid, dimid, len=nlambda)
    end if
    call neasyf_close(ncid)
# else
    call mp_abort("nc_get_lambda_grid_size does not work without NetCDF")
# endif
  end subroutine nc_get_lambda_grid_size

  subroutine nc_get_grid_scalars(filename, shat, drhodpsi, kxfac, qval, rmaj)
    character (len=2000), intent(in) :: filename
    real, intent(out), optional :: shat, drhodpsi, kxfac, qval, rmaj
# ifdef NETCDF
    integer (kind_nf) :: ncid
    ncid = neasyf_open(trim(filename), "r")
    call neasyf_read(ncid, "drhodpsi", drhodpsi)
    call neasyf_read(ncid, "rmaj", rmaj)
    call neasyf_read(ncid, "shat", shat)
    call neasyf_read(ncid, "kxfac", kxfac)
    call neasyf_read(ncid, "q", qval)
    call neasyf_close(ncid)
# else
    call mp_abort("nc_get_grid_scalars does not work without NetCDF")
# endif
  end subroutine nc_get_grid_scalars
  
  subroutine nc_write_grid_scalars(filename, shat, drhodpsi, kxfac, qval, rmaj)
    character (len=2000), intent(in) :: filename
    real, intent(in), optional :: shat, drhodpsi, kxfac, qval, rmaj
# ifdef NETCDF
    integer (kind_nf) :: ncid
    character (len=2) :: action
    logical :: is_exist

    ! always append mode if file exists
    action="w"
    inquire(file=trim(filename), exist=is_exist)
    if (is_exist) then
       action="rw"
    end if

    ncid = neasyf_open(trim(filename), action)
    call neasyf_write(ncid, "drhodpsi", drhodpsi)
    call neasyf_write(ncid, "rmaj", rmaj)
    call neasyf_write(ncid, "shat", shat)
    call neasyf_write(ncid, "kxfac", kxfac)
    call neasyf_write(ncid, "q", qval)
    call neasyf_close(ncid)
# else
    call mp_abort("nc_write_grid_scalars does not work without NetCDF")
# endif
  end subroutine nc_write_grid_scalars

  subroutine nc_get_grids(filename, ntgrid, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, theta, &
       cdrift, cdrift0, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       Bpol, jacob)
    character (len=2000), intent(in) :: filename
    integer, intent(in) :: ntgrid
    real, dimension(-ntgrid:ntgrid), intent(out) :: bmag, gradpar, &
         gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho, theta
    real, dimension(-ntgrid:ntgrid), intent(out), optional :: &
         cdrift, cdrift0, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
         Bpol, jacob
# ifdef NETCDF
    integer (kind_nf) :: ncid, varid
    integer :: status
    integer :: no_end_point_int
    logical :: no_end_point

    ncid = neasyf_open(trim(filename), "r")

    no_end_point = .false.
    status = nf90_inq_varid(ncid, "no_end_point", varid)
    if (status == NF90_NOERR) then
       call neasyf_read(ncid, "no_end_point", no_end_point_int)
       if (no_end_point_int /= 0) no_end_point = .true.
    end if

    call read_and_pad(ncid, gbdrift, "gbdrift", no_end_point)
    call read_and_pad(ncid, gbdrift0, "gbdrift0", no_end_point)
    call read_and_pad(ncid, cvdrift, "cvdrift", no_end_point)
    call read_and_pad(ncid, cvdrift0, "cvdrift0", no_end_point)
    call read_and_pad(ncid, gradpar, "gradpar", no_end_point)
    call read_and_pad(ncid, bmag, "bmag", no_end_point)
    call read_and_pad(ncid, gds2, "gds2", no_end_point)
    call read_and_pad(ncid, gds21, "gds21", no_end_point)
    call read_and_pad(ncid, gds22, "gds22", no_end_point)
    call read_and_pad(ncid, theta, "theta", no_end_point)
    call read_and_pad(ncid, grho, "gho", no_end_point)

    call read_and_pad(ncid, cdrift, "cdrift", no_end_point)
    call read_and_pad(ncid, cdrift0, "cdrift0", no_end_point)
    call read_and_pad(ncid, Bpol, "bpol", no_end_point)
    call read_and_pad(ncid, jacob, "jacob", no_end_point)

    call read_and_pad(ncid, Rplot, "Rplot", no_end_point)
    call read_and_pad(ncid, Rprime, "Rpime", no_end_point)
    call read_and_pad(ncid, Zplot, "Zplot", no_end_point)
    call read_and_pad(ncid, Zprime, "Zpime", no_end_point)
    call read_and_pad(ncid, aplot, "aplot", no_end_point)
    call read_and_pad(ncid, aprime, "apime", no_end_point)

    call neasyf_close(ncid)

# else
    call mp_abort("nc_get_grids does not work without NetCDF")
# endif
  end subroutine nc_get_grids

  subroutine nc_write_grids(filename, ntgrid, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, theta, &
       cdrift, cdrift0, &
       Rplot, Rprime, Zplot, Zprime, aplot, aprime, &
       Bpol, jacob, no_end_point_in)
    character (len=2000), intent(in) :: filename
    integer, intent(in) :: ntgrid
    real, dimension(-ntgrid:ntgrid), intent(in) :: bmag, gradpar, &
         gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho, theta
    real, dimension(-ntgrid:ntgrid), intent(in), optional :: &
         cdrift, cdrift0, &
         Rplot, Rprime, Zplot, Zprime, aplot, aprime, &
         Bpol, jacob
    logical, intent(in), optional :: no_end_point_in
# ifdef NETCDF
    integer (kind_nf) :: ncid, nt_dim
    integer :: theta_end_point
    logical :: no_end_point
    integer :: no_end_point_int

    character (len=2) :: action
    logical :: is_exist

    ! always append mode if file exists
    action="w"
    inquire(file=trim(filename), exist=is_exist)
    if (is_exist) then
       action="rw"
    end if

    ncid = neasyf_open(trim(filename), action)

    no_end_point = .false.
    no_end_point_int = 0
    if (present(no_end_point_in)) then
       no_end_point = no_end_point_in
       if (no_end_point) then
          no_end_point_int = 1
       end if
    end if
    call neasyf_write(ncid, "no_end_point", no_end_point_int)

    if (no_end_point) then
       theta_end_point = ntgrid-1
    else
       theta_end_point = ntgrid
    end if
    call neasyf_dim(ncid, "nt", dim_size=theta_end_point+ntgrid+1, dimid=nt_dim)
    call neasyf_write(ncid, "gbdrift", gbdrift(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "gbdrift0", gbdrift0(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "cvdrift", cvdrift(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "cvdrift0", cvdrift0(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "gradpar", gradpar(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "bmag", bmag(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "gds2", gds2(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "gds21", gds21(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "gds22", gds22(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "theta", theta(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    call neasyf_write(ncid, "grho", grho(-ntgrid:theta_end_point), dim_ids=[nt_dim])

    if (present(Bpol)) then
       call neasyf_write(ncid, "bpol", Bpol(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if

    if (present(jacob)) then
       call neasyf_write(ncid, "jacob", jacob(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if

    if (present(Rplot)) then
       call neasyf_write(ncid, "Rplot", Rplot(-ntgrid:theta_end_point), dim_ids=[nt_dim])
   end if
    if (present(Rprime)) then
       call neasyf_write(ncid, "Rprime", Rprime(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if
    if (present(Zplot)) then
       call neasyf_write(ncid, "Zplot", Zplot(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if
    if (present(Zprime)) then
       call neasyf_write(ncid, "Zprime", Zprime(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if
    if (present(aplot)) then
       call neasyf_write(ncid, "aplot", aplot(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if
    if (present(aprime)) then
       call neasyf_write(ncid, "aprime", aprime(-ntgrid:theta_end_point), dim_ids=[nt_dim])
    end if

    call neasyf_close(ncid)

# else
    call mp_abort("nc_write_grids does not work without NetCDF")
# endif
  end subroutine nc_write_grids

  subroutine nc_get_lambda_grid(filename, nlambda, lambda)
    character (len=2000), intent(in) :: filename
    integer, intent(in) :: nlambda
    real, dimension(nlambda), intent(out) :: lambda
# ifdef NETCDF
    integer (kind_nf) :: ncid
    ncid = neasyf_open(trim(filename), "r")
    call neasyf_read(ncid, "lambda", lambda)
    call neasyf_close(ncid)
# else
    call mp_abort("nc_get_lambda_grid does not work without NetCDF")
# endif
  end subroutine nc_get_lambda_grid

  subroutine nc_write_lambda_grid(filename, nlambda, lambda)
    character (len=2000), intent(in) :: filename
    integer, intent(in) :: nlambda
    real, dimension(nlambda), intent(in) :: lambda
# ifdef NETCDF
    integer (kind_nf) :: ncid, nl_dim
    character (len=2) :: action
    logical :: is_exist

    ! always append mode if file exists
    action="w"
    inquire(file=trim(filename), exist=is_exist)
    if (is_exist) then
       action="rw"
    end if

    ncid = neasyf_open(trim(filename), action)
    call neasyf_dim(ncid, "nlambda", dim_size=nlambda, dimid=nl_dim)
    call neasyf_write(ncid, "lambda", lambda, dim_ids=[nl_dim])
    call neasyf_close(ncid)
# else
    call mp_abort("nc_write_lambda_grid does not work without NetCDF")
# endif
  end subroutine nc_write_lambda_grid

# ifdef NETCDF
  subroutine read_and_pad(ncid, array, name, no_end_point)
    implicit none
    integer (kind_nf), intent(in) :: ncid
    real, dimension(:), intent(in out) :: array
    character (len=*), intent(in) :: name
    logical, intent(in) :: no_end_point
    integer (kind_nf) :: varid
    integer :: end_point
    integer :: status
    end_point = size(array)
    ! array corresponding to the dummy argument array in the upstream subroutine
    ! may be non-present optional argument, whose size will be 1 /= 2*ntgrid+1.!py
    if (end_point == 1) return
    if (no_end_point) then
       end_point = end_point -1
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if (status == NF90_NOERR) then
       call neasyf_read(ncid, name, array(:end_point))
    end if
    if (no_end_point) call ntgrid_padding(array)
  end subroutine read_and_pad
# endif

  subroutine ntgrid_padding(array)
    ! RN>
    ! If grid file does not have the end point data, they are generated
    ! using information on the other end by assuming (anti-)symmetry
    ! w.r.t theta=0.
    ! If symmetric, array(ntgrid) = array(-ntgrid)
    ! If anti-symmetric, array(ntgrid) = -array(-ntgrid)
    !
    ! If this deos not apply, do not use no_end_point = T grid file.
    ! <RN
    real, dimension (:), intent (in out) :: array
    integer :: s
    integer :: l
    l = size(array)
    s = sign(1.,(array(l-1)-array(l-2))*(array(2)-array(3)))
    array(l) = array(l-1) + s * (array(1)-array(2))
  end subroutine ntgrid_padding

end module gs2_io_grids
