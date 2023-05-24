module gs2_io_grids

# ifdef NETCDF
  use netcdf_utils, only: kind_nf
  use neasyf, only: neasyf_open, neasyf_close
  use neasyf, only: neasyf_read, neasyf_write
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: NF90_NOERR
  use mp, only: mp_abort
# endif

  implicit none

  private

  public :: nc_get_grid_sizes, nc_get_grid_scalars, nc_get_grids

contains

  subroutine nc_get_grid_sizes(filename, ntheta, ntgrid, nperiod, nlambda)
    character (len=2000), intent(in) :: filename
    integer, intent(out), optional :: ntheta
    integer, intent(out), optional :: ntgrid
    integer, intent(out), optional :: nperiod
    integer, intent(out), optional :: nlambda

# ifdef NETCDF
    integer (kind_nf) :: ncid, dimid
    integer :: status

    ncid = neasyf_open(trim(filename), "r")

    if (present(nlambda)) then
       status = nf90_inq_dimid(ncid, "nlambda", dimid)
       if (status == NF90_NOERR) then
          status = nf90_inquire_dimension(ncid, dimid, len=nlambda)
       end if
    end if

    if (present(ntgrid)) call neasyf_read(ncid, "ntgrid", ntgrid)
    if (present(nperiod)) call neasyf_read(ncid, "nperiod", nperiod)
    if (present(ntheta)) call neasyf_read(ncid, "ntheta", ntheta)

    call neasyf_close(ncid)

# else
    call mp_abort("nc_get_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_get_grid_sizes

  subroutine nc_get_grid_scalars(filename, shat, drhodpsi, kxfac, qval, rmaj)
    character (len=2000), intent(in) :: filename
    real, intent(out), optional :: shat, drhodpsi, kxfac, qval, rmaj
# ifdef NETCDF
    integer (kind_nf) :: ncid, varid
    integer :: status
    ncid = neasyf_open(trim(filename), "r")
    call neasyf_read(ncid, "drhodpsi", drhodpsi)
    call neasyf_read(ncid, "rmaj", rmaj)
    call neasyf_read(ncid, "shat", shat)
    call neasyf_read(ncid, "kxfac", kxfac)
    call neasyf_read(ncid, "q", qval)
    call neasyf_close(ncid)
# else
    call mp_abort("nc_get_grid_sizes does not work without NetCDF")
# endif
  end subroutine nc_get_grid_scalars
  
  subroutine nc_get_grids(filename, nlambda, ntgrid, lambda, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gds2, gds21, gds22, grho, theta, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       Bpol, jacob)
    character (len=2000), intent(in) :: filename
    integer, intent(in) :: nlambda, ntgrid
    real, dimension(nlambda), intent(out) :: lambda
    real, dimension(-ntgrid:ntgrid), intent(out) :: bmag, gradpar, &
         gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gds2, gds21, gds22, grho, theta
    real, dimension(-ntgrid:ntgrid), intent(out), optional :: &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
         Bpol, jacob
# ifdef NETCDF
    integer (kind_nf) :: ncid, varid
    integer :: status
    integer :: ulim
    integer :: no_end_point_int
    logical :: no_end_point

    ncid = neasyf_open(trim(filename), "r")

    no_end_point = .false.
    status = nf90_inq_varid(ncid, "no_end_point", varid)
    if (status == NF90_NOERR) then
       call neasyf_read(ncid, "no_end_point", no_end_point_int)
       if (no_end_point_int /= 0) no_end_point = .true.
    end if
    
    call neasyf_read(ncid, "lambda", lambda)

    ulim = ntgrid
    if (no_end_point) ulim = ntgrid-1
    
    call neasyf_read(ncid, "gbdrift", gbdrift(-ntgrid:ulim))
    call neasyf_read(ncid, "gbdrift0", gbdrift0(-ntgrid:ulim))
    call neasyf_read(ncid, "cvdrift", cvdrift(-ntgrid:ulim))
    call neasyf_read(ncid, "cvdrift0", cvdrift0(-ntgrid:ulim))
    call neasyf_read(ncid, "gradpar", gradpar(-ntgrid:ulim))
    call neasyf_read(ncid, "bmag", bmag(-ntgrid:ulim))
    call neasyf_read(ncid, "gds2", gds2(-ntgrid:ulim))
    call neasyf_read(ncid, "gds21", gds21(-ntgrid:ulim))
    call neasyf_read(ncid, "gds22", gds22(-ntgrid:ulim))
    call neasyf_read(ncid, "theta", theta(-ntgrid:ulim))
    call neasyf_read(ncid, "grho", grho(-ntgrid:ulim))

    if (no_end_point) then
       call ntgrid_padding(gbdrift)
       call ntgrid_padding(gbdrift0)
       call ntgrid_padding(cvdrift)
       call ntgrid_padding(cvdrift0)
       call ntgrid_padding(gradpar)
       call ntgrid_padding(bmag)
       call ntgrid_padding(gds2)
       call ntgrid_padding(gds21)
       call ntgrid_padding(gds22)
       call ntgrid_padding(theta)
       call ntgrid_padding(grho)
    end if

    if (present(Bpol)) then
       status = nf90_inq_varid(ncid, "bpol", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "bpol", Bpol(-ntgrid:ulim))
       end if
       if (no_end_point) call ntgrid_padding(Bpol)
    end if

    if (present(jacob)) then
       status = nf90_inq_varid(ncid, "jacob", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "jacob", jacob(-ntgrid:ulim))
       end if
       if (no_end_point) call ntgrid_padding(jacob)
    end if

    if (present(Rplot)) then
       status = nf90_inq_varid(ncid, "Rplot", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "Rplot", Rplot(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(Rplot)
   end if
    if (present(Rprime)) then
       status = nf90_inq_varid(ncid, "Rprime", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "Rprime", Rprime(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(Rprime)
    end if
    if (present(Zplot)) then
       status = nf90_inq_varid(ncid, "Zplot", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "Zplot", Zplot(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(Zplot)
    end if
    if (present(Zprime)) then
       status = nf90_inq_varid(ncid, "Zprime", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "Zprime", Zprime(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(Zprime)
    end if

    if (present(aplot)) then
       status = nf90_inq_varid(ncid, "aplot", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "aplot", aplot(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(aplot)
    end if
    if (present(aprime)) then
       status = nf90_inq_varid(ncid, "aprime", varid)
       if (status == NF90_NOERR) then
          call neasyf_read(ncid, "aprime", aprime(-ntgrid:ulim))
       end if
        if (no_end_point) call ntgrid_padding(aprime)
    end if

    call neasyf_close(ncid)

# else
    call mp_abort("nc_get_grids does not work without NetCDF")
# endif
  end subroutine nc_get_grids

  subroutine ntgrid_padding(array)
    real, dimension (:), intent (in out) :: array
    integer :: s
    integer :: l
    l = size(array)
    s = sign(1.,(array(l-1)-array(l-2))*(array(2)-array(3)))
    array(l) = array(l-1) + s * (array(1)-array(2))
  end subroutine ntgrid_padding

end module gs2_io_grids
