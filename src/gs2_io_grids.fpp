!> A module for writing geometry to netcdf.
module gs2_io_grids

# ifdef NETCDF
  use netcdf_utils, only: kind_nf
  use neasyf, only: neasyf_open, neasyf_close
  use neasyf, only: neasyf_read, neasyf_write
  use neasyf, only: neasyf_dim
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: NF90_NOERR
# endif
  use mp, only: proc0, mp_abort

  implicit none

  private

  public :: nc_get_grid_sizes, nc_write_grid_sizes
  public :: nc_get_lambda_grid_size
  public :: nc_get_grid_scalars, nc_write_grid_scalars
  public :: nc_get_grids, nc_write_grids
  public :: nc_get_lambda_grid, nc_write_lambda_grid
  public :: nc_grid_file_open, nc_grid_file_close

# ifdef NETCDF
  integer (kind_nf) :: ncid
# endif

contains

  !> Remove the output netCDF file if already exists
  subroutine nc_grid_file_new(filename)
    use file_utils, only: get_unused_unit
    character (len=*), intent(in) :: filename
    logical :: is_exist
    integer :: unit
    if (.not.proc0) return
    inquire(file=filename, exist=is_exist)
    if (is_exist) then
       write(6,*) 'WARNING: existing file will be over-written ', trim(filename)
       open(unit, file=filename)
       close(unit, status='delete')
    end if
  end subroutine nc_grid_file_new

  !> Open the netCDF file
  !> If the file already exists, remove it.
  subroutine nc_grid_file_open(filename, mode)
    use file_utils, only: get_unused_unit
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: mode
    logical :: is_exist
    integer :: unit
# ifdef NETCDF
    if (.not.proc0) return
    if (mode == "w") then
       inquire(file=filename, exist=is_exist)
       if (is_exist) then
          write(6,*) 'WARNING: existing file will be over-written ', trim(filename)
          open(unit, file=filename)
          close(unit, status='delete')
       end if
    end if
    ncid = neasyf_open(trim(filename), mode)
# else
    call mp_abort("nc_grid_file_open does not work without NetCDF")
# endif
  end subroutine nc_grid_file_open

  !> Close the netCDF file
  subroutine nc_grid_file_close()
# ifdef NETCDF
    if (.not.proc0) return
    call neasyf_close(ncid)
# else
    call mp_abort("nc_grid_file_close does not work without NetCDF")
# endif
  end subroutine nc_grid_file_close

  !> Read sizes of geometry data from the input netCDF file
  subroutine nc_get_grid_sizes(ntheta, ntgrid, nperiod)
    integer, intent(out), optional :: ntheta, ntgrid, nperiod
# ifdef NETCDF
    if (.not.proc0) return
    if (present(ntgrid)) call neasyf_read(ncid, "ntgrid", ntgrid)
    if (present(nperiod)) call neasyf_read(ncid, "nperiod", nperiod)
    if (present(ntheta)) call neasyf_read(ncid, "ntheta", ntheta)
# else
    call mp_abort("nc_get_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_get_grid_sizes

  !> Writes sizes of geometry data to the output netCDF file
  subroutine nc_write_grid_sizes(ntheta, ntgrid, nperiod)
    integer, intent(in), optional :: ntheta, ntgrid, nperiod
# ifdef NETCDF
    if (.not.proc0) return
    if (present(ntgrid)) call neasyf_write(ncid, "ntgrid", ntgrid)
    if (present(nperiod)) call neasyf_write(ncid, "nperiod", nperiod)
    if (present(ntheta)) call neasyf_write(ncid, "ntheta", ntheta)
# else
    call mp_abort("nc_write_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_write_grid_sizes

  !> Read sizes of lambda grid data from the input netCDF file
  subroutine nc_get_lambda_grid_size(nlambda)
    integer, intent(out) :: nlambda
# ifdef NETCDF
    integer (kind_nf) :: dimid
    integer :: status
    if (.not.proc0) return
    status = nf90_inq_dimid(ncid, "nlambda", dimid)
    if (status == NF90_NOERR) then
       status = nf90_inquire_dimension(ncid, dimid, len=nlambda)
    else
       call mp_abort("failed to get nlambda fron netCDF grid file")
    end if
# else
    call mp_abort("nc_get_lambda_grid_size does not work without NetCDF")
# endif
  end subroutine nc_get_lambda_grid_size

  !> Read scalar geometry data from the input netCDF file
  subroutine nc_get_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj, B_T, aminor, grhoavg, surfarea)
    real, intent(out), optional :: shat, drhodpsi, kxfac, qval, rmaj, B_T, aminor, grhoavg, surfarea
# ifdef NETCDF
    if (.not.proc0) return
    if (present(drhodpsi)) call neasyf_read(ncid, "drhodpsi", drhodpsi)
    if (present(rmaj)) call neasyf_read(ncid, "rmaj", rmaj)
    if (present(shat)) call neasyf_read(ncid, "shat", shat)
    if (present(kxfac)) call neasyf_read(ncid, "kxfac", kxfac)
    if (present(qval)) call neasyf_read(ncid, "q", qval)
    if (present(B_T)) call neasyf_read(ncid, "B_T", B_T)
    if (present(aminor)) call neasyf_read(ncid, "aminor", aminor)
    if (present(grhoavg)) call neasyf_read(ncid, "grhoavg", grhoavg)
    if (present(surfarea)) call neasyf_read(ncid, "surfarea", surfarea)
# else
    call mp_abort("nc_get_grid_scalars does not work without NetCDF")
# endif
  end subroutine nc_get_grid_scalars
  
  !> Write scalar geometry data to the output netCDF file
  subroutine nc_write_grid_scalars(shat, drhodpsi, kxfac, qval, rmaj, B_T, aminor, grhoavg, surfarea)
    real, intent(in), optional :: shat, drhodpsi, kxfac, qval, rmaj, B_T, aminor, grhoavg, surfarea
# ifdef NETCDF
    if (.not.proc0) return
    if (present(drhodpsi)) call neasyf_write(ncid, "drhodpsi", drhodpsi)
    if (present(rmaj)) call neasyf_write(ncid, "rmaj", rmaj)
    if (present(shat)) call neasyf_write(ncid, "shat", shat)
    if (present(kxfac)) call neasyf_write(ncid, "kxfac", kxfac)
    if (present(qval)) call neasyf_write(ncid, "q", qval)
    if (present(B_T)) call neasyf_write(ncid, "B_T", B_T)
    if (present(aminor)) call neasyf_write(ncid, "aminor", aminor)
    if (present(grhoavg)) call neasyf_write(ncid, "grhoavg", grhoavg)
    if (present(surfarea)) call neasyf_write(ncid, "surfarea", surfarea)
# else
    call mp_abort("nc_write_grid_scalars does not work without NetCDF")
# endif
  end subroutine nc_write_grid_scalars

  !> Read geometry data from the input netCDF file
  subroutine nc_get_grids(ntgrid, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, theta, &
       cdrift, cdrift0, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       Bpol, jacob)
    integer, intent(in) :: ntgrid
    real, dimension(-ntgrid:ntgrid), intent(out) :: bmag, gradpar, &
         gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho, theta
    real, dimension(-ntgrid:ntgrid), intent(out), optional :: &
         cdrift, cdrift0, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
         Bpol, jacob
# ifdef NETCDF
    integer (kind_nf) :: varid
    integer :: status
    integer :: no_end_point_int
    logical :: no_end_point

    if (.not.proc0) return

    !> The no_end_point flag in the input netCDF file tells if the grid data
    !> have the duplicate end point data at ntgrid.
    !> The default no_end_point behavior:
    !>  The grid file generated by this module should always store the no_end_point flag.
    !>  The grid file came from outside (GX) may not have this flag.
    !>  (GX grid file does not have the end point data.)
    no_end_point = .true.
    status = nf90_inq_varid(ncid, "no_end_point", varid)
    if (status == NF90_NOERR) then
       call neasyf_read(ncid, "no_end_point", no_end_point_int)
       if (no_end_point_int /= 0) no_end_point = .true.
    end if

    !> read_and_pad handles the end point data
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
    call read_and_pad(ncid, grho, "grho", no_end_point)

    if (present(cdrift)) call read_and_pad(ncid, cdrift, "cdrift", no_end_point)
    if (present(cdrift0)) call read_and_pad(ncid, cdrift0, "cdrift0", no_end_point)
    if (present(Bpol)) call read_and_pad(ncid, Bpol, "bpol", no_end_point)
    if (present(jacob)) call read_and_pad(ncid, jacob, "jacob", no_end_point)

    if (present(Rplot)) call read_and_pad(ncid, Rplot, "Rplot", no_end_point)
    if (present(Rprime)) call read_and_pad(ncid, Rprime, "Rprime", no_end_point)
    if (present(Zplot)) call read_and_pad(ncid, Zplot, "Zplot", no_end_point)
    if (present(Zprime)) call read_and_pad(ncid, Zprime, "Zprime", no_end_point)
    if (present(aplot)) call read_and_pad(ncid, aplot, "aplot", no_end_point)
    if (present(aprime)) call read_and_pad(ncid, aprime, "aprime", no_end_point)
# else
    call mp_abort("nc_get_grids does not work without NetCDF")
# endif
  end subroutine nc_get_grids

  !> Write geometry data to the output netCDF file
  subroutine nc_write_grids(ntgrid, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, theta, &
       cdrift, cdrift0, &
       Rplot, Rprime, Zplot, Zprime, aplot, aprime, &
       Bpol, jacob, no_end_point_in)
    use optionals, only: get_option_with_default
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
    integer (kind_nf) :: nt_dim
    integer :: theta_end_point
    logical :: no_end_point
    integer :: no_end_point_int

    if (.not.proc0) return

    !> Get the no_end_point flag from the input
    !> and store it in the output netCDF file.
    no_end_point = get_option_with_default(no_end_point_in, .false.)
    no_end_point_int = 0
    if (no_end_point) no_end_point_int = 1
    call neasyf_write(ncid, "no_end_point", no_end_point_int)

    !> If no_end_point=.true., data at ntgrid are not stored.
    !> Data at ntgrid may be generated from the data at -ntgrid with some symmetry assumption.
    !> This option is meant to provide compatibility with GX.
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

    if (present(Bpol)) call neasyf_write(ncid, "bpol", Bpol(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])

    if (present(jacob)) call neasyf_write(ncid, "jacob", jacob(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])

    if (present(Rplot)) call neasyf_write(ncid, "Rplot", Rplot(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(Rprime)) call neasyf_write(ncid, "Rprime", Rprime(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(Zplot)) call neasyf_write(ncid, "Zplot", Zplot(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(Zprime)) call neasyf_write(ncid, "Zprime", Zprime(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(aplot)) call neasyf_write(ncid, "aplot", aplot(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(aprime)) call neasyf_write(ncid, "aprime", aprime(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(cdrift)) call neasyf_write(ncid, "cdrift", cdrift(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
    if (present(cdrift0)) call neasyf_write(ncid, "cdrift0", cdrift0(-ntgrid:theta_end_point), &
         dim_ids=[nt_dim])
# else
    call mp_abort("nc_write_grids does not work without NetCDF")
# endif
  end subroutine nc_write_grids

  !> Read lambda grid data from the input netCDF file
  subroutine nc_get_lambda_grid(nlambda, lambda)
    integer, intent(in) :: nlambda
    real, dimension(nlambda), intent(out) :: lambda
# ifdef NETCDF
    if (.not.proc0) return
    call neasyf_read(ncid, "lambda", lambda)
# else
    call mp_abort("nc_get_lambda_grid does not work without NetCDF")
# endif
  end subroutine nc_get_lambda_grid

  !> Write lambda grid data to the output netCDF file
  subroutine nc_write_lambda_grid(nlambda, lambda)
    integer, intent(in) :: nlambda
    real, dimension(nlambda), intent(in) :: lambda
# ifdef NETCDF
    integer (kind_nf) :: nl_dim
    if (.not.proc0) return
    call neasyf_dim(ncid, "nlambda", dim_size=nlambda, dimid=nl_dim)
    call neasyf_write(ncid, "lambda", lambda, dim_ids=[nl_dim])
# else
    call mp_abort("nc_write_lambda_grid does not work without NetCDF")
# endif
  end subroutine nc_write_lambda_grid

# ifdef NETCDF
  !> Helper subroutine to appropriately get geometry data depending on the no_end_point flag.
  !> Read geometry data from the input netCDF file with an appropriate size
  !> and pad ntgrid data if no_end_point=.true.
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
    !> The array corresponding to the dummy argument array in the upstream subroutine
    !> may be non-present optional argument, whose size will be 1 /= 2*ntgrid+1.
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

  !> Generate data at ntgrid point using some symmetry assumption
  subroutine ntgrid_padding(array)
    !> If grid file does not have the end point data, they are generated
    !> using information on the other end by assuming (anti-)symmetry
    !> w.r.t theta=0.
    !> If symmetric, array(ntgrid) = array(-ntgrid)
    !> If anti-symmetric, array(ntgrid) = -array(-ntgrid)
    !>
    !> If this deos not apply, do not use no_end_point = T grid file.
    real, dimension (:), intent (in out) :: array
    integer :: s
    integer :: l
    l = size(array)
    s = int(sign(1.,(array(l-1)-array(l-2))*(array(2)-array(3))))
    array(l) = array(l-1) + s * (array(1)-array(2))
  end subroutine ntgrid_padding

end module gs2_io_grids
