!> netCDF writer for the GS2 GPU export program `gs2_export`.
!>
!> File/group handling, variable attributes and the canonicalisers that expand GS2's
!> distributed layouts (g_lo, le_lo) into the layout-independent canonical index order of
!> gpu/spec/09_export_infra.md §0.2 and §2.4:
!>
!>  * every array is written in Fortran dimension order (first index fastest), so NumPy shows
!>    the dimensions reversed;
!>  * complex values get a leading Fortran dimension `ri` (1 = real, 2 = imaginary);
!>  * logicals are written as int8 0/1 with the attribute fortran_type = "logical";
!>  * integers keep their Fortran values (1-based indices stay 1-based, -1 sentinels stay -1).
!>
!> Only proc0 touches netCDF. The exporter is restricted to a single MPI rank (hard gate in
!> gs2_export), so the g_lo/le_lo arrays of rank 0 hold every point and no gather is needed.
!> All routines are no-ops on other ranks.
module gs2_export_writer
  use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real64
  implicit none

  private

  public :: ex_create, ex_close, ex_group, ex_sync, ex_add_dim, ex_has_dim, ex_dim_len
  public :: ex_att, ex_absent, ex_scalar, ex_vec
  public :: ex_g7c, ex_g7r, ex_g6r, ex_g7m1c, ex_glc, ex_g5i, ex_g5l, ex_g5sl
  public :: ex_le6r, ex_le6c
  public :: ex_compression
  public :: ex_fail

  !> Deflate level for new variables (0 = none; shuffle+deflate are lossless)
  integer :: ex_compression = 0

  interface ex_att
     module procedure ex_att_text
     module procedure ex_att_int
     module procedure ex_att_real
  end interface ex_att

  interface ex_scalar
     module procedure ex_scalar_r
     module procedure ex_scalar_i
     module procedure ex_scalar_l
  end interface ex_scalar

  interface ex_vec
     module procedure ex_vec_r1, ex_vec_r2, ex_vec_r3, ex_vec_r4
     module procedure ex_vec_i1, ex_vec_i2, ex_vec_i3
     module procedure ex_vec_l1, ex_vec_l2, ex_vec_l3
     module procedure ex_vec_c1, ex_vec_c2, ex_vec_c3
  end interface ex_vec

contains

  !> Abort with a message (collective abort through mp_abort)
  subroutine ex_fail(message)
    use mp, only: mp_abort
    implicit none
    character(len=*), intent(in) :: message
    call mp_abort('gs2_export: '//trim(message), .true.)
  end subroutine ex_fail

  !> Check a netCDF status
  subroutine check(status, what)
    use netcdf, only: nf90_noerr, nf90_strerror
    implicit none
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) call ex_fail('netCDF error ('//trim(what)//'): '//trim(nf90_strerror(status)))
  end subroutine check

  !> Create a NETCDF4 file (proc0 only; returns -1 elsewhere) with the GS2 standard metadata
  integer function ex_create(fname, title, header) result(ncid)
    use netcdf, only: nf90_create, NF90_CLOBBER, NF90_NETCDF4
    use mp, only: proc0
    use gs2_metadata, only: create_metadata
    use standard_header, only: standard_header_type
    implicit none
    character(len=*), intent(in) :: fname, title
    type(standard_header_type), intent(in), optional :: header
    ncid = -1
    if (.not. proc0) return
    call check(nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid), 'create '//trim(fname))
    call create_metadata(ncid, title, header)
  end function ex_create

  !> Close a file
  subroutine ex_close(ncid)
    use netcdf, only: nf90_close
    use mp, only: proc0
    implicit none
    integer, intent(in) :: ncid
    if (.not. proc0) return
    if (ncid < 0) return
    call check(nf90_close(ncid), 'close')
  end subroutine ex_close

  !> Flush a file to disk
  subroutine ex_sync(ncid)
    use netcdf, only: nf90_sync
    use mp, only: proc0
    implicit none
    integer, intent(in) :: ncid
    if (.not. proc0) return
    if (ncid < 0) return
    call check(nf90_sync(ncid), 'sync')
  end subroutine ex_sync

  !> Return the id of the group `path` (slash-separated, relative to `parent`), creating
  !> the groups that do not exist yet
  integer function ex_group(parent, path) result(gid)
    use netcdf, only: nf90_inq_ncid, nf90_def_grp, nf90_noerr
    use mp, only: proc0
    implicit none
    integer, intent(in) :: parent
    character(len=*), intent(in) :: path
    integer :: i0, i1, pid, status
    character(len=len(path)) :: p
    gid = -1
    if (.not. proc0) return
    p = adjustl(path)
    pid = parent
    i0 = 1
    if (p(1:1) == '/') i0 = 2
    do while (i0 <= len_trim(p))
       i1 = index(p(i0:), '/')
       if (i1 == 0) then
          i1 = len_trim(p) + 1
       else
          i1 = i0 + i1 - 1
       end if
       if (i1 > i0) then
          status = nf90_inq_ncid(pid, p(i0:i1-1), gid)
          if (status /= nf90_noerr) call check(nf90_def_grp(pid, p(i0:i1-1), gid), 'def_grp '//p(i0:i1-1))
          pid = gid
       end if
       i0 = i1 + 1
    end do
    gid = pid
  end function ex_group

  !> Define a dimension in group `gid` (normally the root group), optionally with a
  !> coordinate variable of the same name
  subroutine ex_add_dim(gid, name, n, rvalues, ivalues)
    use netcdf, only: nf90_def_dim, nf90_def_var, nf90_put_var, NF90_DOUBLE, NF90_INT
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid, n
    character(len=*), intent(in) :: name
    real, dimension(:), intent(in), optional :: rvalues
    integer, dimension(:), intent(in), optional :: ivalues
    integer :: dimid, varid
    if (.not. proc0) return
    call check(nf90_def_dim(gid, name, n, dimid), 'def_dim '//name)
    if (present(rvalues)) then
       call check(nf90_def_var(gid, name, NF90_DOUBLE, [dimid], varid), 'def_var '//name)
       call check(nf90_put_var(gid, varid, real(rvalues, real64)), 'put_var '//name)
    else if (present(ivalues)) then
       call check(nf90_def_var(gid, name, NF90_INT, [dimid], varid), 'def_var '//name)
       call check(nf90_put_var(gid, varid, int(ivalues, int32)), 'put_var '//name)
    end if
  end subroutine ex_add_dim

  !> Does dimension `name` exist (searching `gid` and its ancestors)?
  logical function ex_has_dim(gid, name)
    use netcdf, only: nf90_inq_dimid, nf90_noerr
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    integer :: dimid
    ex_has_dim = .false.
    if (.not. proc0) return
    ex_has_dim = nf90_inq_dimid(gid, trim(name), dimid) == nf90_noerr
  end function ex_has_dim

  !> Length of dimension `name`
  integer function ex_dim_len(gid, name)
    use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    integer :: dimid
    ex_dim_len = -1
    if (.not. proc0) return
    call check(nf90_inq_dimid(gid, trim(name), dimid), 'inq_dimid '//trim(name))
    call check(nf90_inquire_dimension(gid, dimid, len=ex_dim_len), 'inquire_dimension '//trim(name))
  end function ex_dim_len

  !> Split a space separated list of dimension names
  subroutine split_dims(dims, names, n)
    implicit none
    character(len=*), intent(in) :: dims
    character(len=64), dimension(:), intent(out) :: names
    integer, intent(out) :: n
    integer :: i, i0
    logical :: inword
    n = 0 ; inword = .false. ; i0 = 1
    do i = 1, len(dims) + 1
       if (i <= len(dims)) then
          if (dims(i:i) /= ' ') then
             if (.not. inword) then
                inword = .true. ; i0 = i
             end if
             cycle
          end if
       end if
       if (inword) then
          n = n + 1
          if (n > size(names)) call ex_fail('too many dimensions in "'//trim(dims)//'"')
          names(n) = dims(i0:i-1)
          inword = .false.
       end if
    end do
  end subroutine split_dims

  !> Look up a variable or define it (dims: space separated, Fortran order)
  integer function def_var(gid, name, xtype, dims) result(varid)
    use netcdf, only: nf90_inq_varid, nf90_inq_dimid, nf90_def_var, nf90_def_var_deflate
    use netcdf, only: nf90_noerr
    implicit none
    integer, intent(in) :: gid, xtype
    character(len=*), intent(in) :: name, dims
    character(len=64), dimension(16) :: names
    integer, dimension(16) :: dimids
    integer :: n, i
    if (nf90_inq_varid(gid, name, varid) == nf90_noerr) return
    call split_dims(dims, names, n)
    do i = 1, n
       call check(nf90_inq_dimid(gid, trim(names(i)), dimids(i)), &
            'inq_dimid '//trim(names(i))//' for '//trim(name))
    end do
    if (n == 0) then
       call check(nf90_def_var(gid, name, xtype, varid), 'def_var '//trim(name))
    else
       call check(nf90_def_var(gid, name, xtype, dimids(1:n), varid), 'def_var '//trim(name))
       if (ex_compression > 0) call check(nf90_def_var_deflate(gid, varid, 1, 1, ex_compression), &
            'def_var_deflate '//trim(name))
    end if
  end function def_var

  !> The Fortran lower bounds implied by the canonical dimension names
  function lbounds_of(dims) result(lb)
    use theta_grid, only: ntgrid
    implicit none
    character(len=*), intent(in) :: dims
    character(len=256) :: lb
    character(len=64), dimension(16) :: names
    character(len=16) :: tmp
    integer :: n, i
    call split_dims(dims, names, n)
    lb = ''
    do i = 1, n
       select case (trim(names(i)))
       case ('ri')
          cycle
       case ('theta', 'theta_m1')
          write (tmp, '(I0)') -ntgrid
       case ('glo')
          tmp = '0'
       case default
          tmp = '1'
       end select
       if (len_trim(lb) == 0) then
          lb = trim(tmp)
       else
          lb = trim(lb)//','//trim(tmp)
       end if
    end do
  end function lbounds_of

  !> Write the standard variable attributes (gpu/spec/09_export_infra.md §6.1.2)
  subroutine set_meta(gid, varid, ftype, tag, dims, fname, src, owner, dtdep, ibase, decl)
    use netcdf, only: nf90_put_att
    use gs2_time, only: code_dt
    implicit none
    integer, intent(in) :: gid, varid
    character(len=*), intent(in) :: ftype, tag, dims
    character(len=*), intent(in), optional :: fname, src, owner, decl
    logical, intent(in), optional :: dtdep
    integer, intent(in), optional :: ibase
    character(len=256) :: lb
    call check(nf90_put_att(gid, varid, 'fortran_type', ftype), 'put_att fortran_type')
    call check(nf90_put_att(gid, varid, 'canonical', tag), 'put_att canonical')
    lb = lbounds_of(dims)
    if (len_trim(lb) > 0) call check(nf90_put_att(gid, varid, 'fortran_lbounds', trim(lb)), 'put_att lbounds')
    if (present(fname)) call check(nf90_put_att(gid, varid, 'fortran_name', fname), 'put_att fortran_name')
    if (present(decl)) call check(nf90_put_att(gid, varid, 'fortran_decl', decl), 'put_att fortran_decl')
    if (present(src)) call check(nf90_put_att(gid, varid, 'source', src), 'put_att source')
    if (present(owner)) call check(nf90_put_att(gid, varid, 'owner_spec', owner), 'put_att owner_spec')
    if (present(ibase)) call check(nf90_put_att(gid, varid, 'index_base', int(ibase, int32)), 'put_att index_base')
    if (present(dtdep)) then
       if (dtdep) then
          call check(nf90_put_att(gid, varid, 'dt_dependent', 1_int32), 'put_att dt_dependent')
          call check(nf90_put_att(gid, varid, 'code_dt', real(code_dt, real64)), 'put_att code_dt')
       else
          call check(nf90_put_att(gid, varid, 'dt_dependent', 0_int32), 'put_att dt_dependent')
       end if
    end if
  end subroutine set_meta

  !-------------------------------------------------------------------------- attributes

  !> Text attribute of a variable, or of the group/file if varname == ''
  subroutine ex_att_text(gid, varname, name, value)
    use netcdf, only: nf90_put_att, nf90_inq_varid, NF90_GLOBAL
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: varname, name, value
    integer :: varid
    if (.not. proc0) return
    if (len_trim(varname) == 0) then
       varid = NF90_GLOBAL
    else
       call check(nf90_inq_varid(gid, varname, varid), 'inq_varid '//varname)
    end if
    if (len(value) == 0) then
       call check(nf90_put_att(gid, varid, name, ' '), 'put_att '//name)
    else
       call check(nf90_put_att(gid, varid, name, value), 'put_att '//name)
    end if
  end subroutine ex_att_text

  !> Integer attribute
  subroutine ex_att_int(gid, varname, name, value)
    use netcdf, only: nf90_put_att, nf90_inq_varid, NF90_GLOBAL
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: varname, name
    integer, intent(in) :: value
    integer :: varid
    if (.not. proc0) return
    if (len_trim(varname) == 0) then
       varid = NF90_GLOBAL
    else
       call check(nf90_inq_varid(gid, varname, varid), 'inq_varid '//varname)
    end if
    call check(nf90_put_att(gid, varid, name, int(value, int32)), 'put_att '//name)
  end subroutine ex_att_int

  !> Real attribute
  subroutine ex_att_real(gid, varname, name, value)
    use netcdf, only: nf90_put_att, nf90_inq_varid, NF90_GLOBAL
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: varname, name
    real, intent(in) :: value
    integer :: varid
    if (.not. proc0) return
    if (len_trim(varname) == 0) then
       varid = NF90_GLOBAL
    else
       call check(nf90_inq_varid(gid, varname, varid), 'inq_varid '//varname)
    end if
    call check(nf90_put_att(gid, varid, name, real(value, real64)), 'put_att '//name)
  end subroutine ex_att_real

  !> Append `name` to the space separated group attribute `absent` (arrays not allocated in
  !> GS2, or deliberately not exported, gpu/spec/09_export_infra.md §6.1.2)
  subroutine ex_absent(gid, name)
    use netcdf, only: nf90_inquire_attribute, nf90_get_att, nf90_put_att, NF90_GLOBAL, nf90_noerr
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: old
    integer :: alen
    if (.not. proc0) return
    if (nf90_inquire_attribute(gid, NF90_GLOBAL, 'absent', len=alen) == nf90_noerr) then
       allocate(character(len=alen) :: old)
       call check(nf90_get_att(gid, NF90_GLOBAL, 'absent', old), 'get_att absent')
       if (len_trim(old) > 0) then
          call check(nf90_put_att(gid, NF90_GLOBAL, 'absent', trim(old)//' '//trim(name)), 'put_att absent')
       else
          call check(nf90_put_att(gid, NF90_GLOBAL, 'absent', trim(name)), 'put_att absent')
       end if
    else
       call check(nf90_put_att(gid, NF90_GLOBAL, 'absent', trim(name)), 'put_att absent')
    end if
  end subroutine ex_absent

  !-------------------------------------------------------------------------- scalars

  !> 0-D real variable
  subroutine ex_scalar_r(gid, name, value, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    real, intent(in) :: value
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    integer :: varid
    if (.not. proc0) return
    varid = def_var(gid, name, NF90_DOUBLE, '')
    call check(nf90_put_var(gid, varid, real(value, real64)), 'put_var '//name)
    call set_meta(gid, varid, 'real', 'scalar', '', fname, src, owner, dtdep)
  end subroutine ex_scalar_r

  !> 0-D integer variable
  subroutine ex_scalar_i(gid, name, value, fname, src, owner, ibase)
    use netcdf, only: nf90_put_var, NF90_INT
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    character(len=*), intent(in), optional :: fname, src, owner
    integer, intent(in), optional :: ibase
    integer :: varid
    if (.not. proc0) return
    varid = def_var(gid, name, NF90_INT, '')
    call check(nf90_put_var(gid, varid, int(value, int32)), 'put_var '//name)
    call set_meta(gid, varid, 'integer', 'scalar', '', fname, src, owner, ibase=ibase)
  end subroutine ex_scalar_i

  !> 0-D logical variable (int8 0/1)
  subroutine ex_scalar_l(gid, name, value, fname, src, owner)
    use netcdf, only: nf90_put_var, NF90_BYTE
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    character(len=*), intent(in), optional :: fname, src, owner
    integer :: varid
    if (.not. proc0) return
    varid = def_var(gid, name, NF90_BYTE, '')
    call check(nf90_put_var(gid, varid, merge(1_int8, 0_int8, value)), 'put_var '//name)
    call set_meta(gid, varid, 'logical', 'scalar', '', fname, src, owner)
  end subroutine ex_scalar_l

  !-------------------------------------------------------------------------- flat writers

  !> Write a flat real buffer as variable `name` with dims (Fortran order)
  subroutine put_flat_r(gid, name, buf, dims, ftype, tag, fname, src, owner, dtdep, decl)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims, ftype, tag
    real(real64), dimension(:), intent(in) :: buf
    character(len=*), intent(in), optional :: fname, src, owner, decl
    logical, intent(in), optional :: dtdep
    character(len=64), dimension(16) :: names
    integer, dimension(16) :: counts
    integer :: varid, n, i
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    call split_dims(dims, names, n)
    do i = 1, n
       counts(i) = ex_dim_len(gid, names(i))
    end do
    if (product(int(counts(1:n), int64)) /= size(buf, kind=int64)) &
         call ex_fail('size mismatch writing '//trim(name)//' with dims "'//trim(dims)//'"')
    call check(nf90_put_var(gid, varid, buf, start=[(1, i=1, n)], count=counts(1:n)), 'put_var '//trim(name))
    call set_meta(gid, varid, ftype, tag, dims, fname, src, owner, dtdep, decl=decl)
  end subroutine put_flat_r

  !> Write a flat int32 buffer
  subroutine put_flat_i(gid, name, buf, dims, ftype, tag, fname, src, owner, ibase, decl)
    use netcdf, only: nf90_put_var, NF90_INT
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims, ftype, tag
    integer(int32), dimension(:), intent(in) :: buf
    character(len=*), intent(in), optional :: fname, src, owner, decl
    integer, intent(in), optional :: ibase
    character(len=64), dimension(16) :: names
    integer, dimension(16) :: counts
    integer :: varid, n, i
    varid = def_var(gid, name, NF90_INT, dims)
    call split_dims(dims, names, n)
    do i = 1, n
       counts(i) = ex_dim_len(gid, names(i))
    end do
    if (product(int(counts(1:n), int64)) /= size(buf, kind=int64)) &
         call ex_fail('size mismatch writing '//trim(name)//' with dims "'//trim(dims)//'"')
    call check(nf90_put_var(gid, varid, buf, start=[(1, i=1, n)], count=counts(1:n)), 'put_var '//trim(name))
    call set_meta(gid, varid, ftype, tag, dims, fname, src, owner, ibase=ibase, decl=decl)
  end subroutine put_flat_i

  !> Write a flat int8 buffer (logicals)
  subroutine put_flat_b(gid, name, buf, dims, tag, fname, src, owner, decl)
    use netcdf, only: nf90_put_var, NF90_BYTE
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims, tag
    integer(int8), dimension(:), intent(in) :: buf
    character(len=*), intent(in), optional :: fname, src, owner, decl
    character(len=64), dimension(16) :: names
    integer, dimension(16) :: counts
    integer :: varid, n, i
    varid = def_var(gid, name, NF90_BYTE, dims)
    call split_dims(dims, names, n)
    do i = 1, n
       counts(i) = ex_dim_len(gid, names(i))
    end do
    if (product(int(counts(1:n), int64)) /= size(buf, kind=int64)) &
         call ex_fail('size mismatch writing '//trim(name)//' with dims "'//trim(dims)//'"')
    call check(nf90_put_var(gid, varid, buf, start=[(1, i=1, n)], count=counts(1:n)), 'put_var '//trim(name))
    call set_meta(gid, varid, 'logical', tag, dims, fname, src, owner, decl=decl)
  end subroutine put_flat_b

  !> Complex flat buffer -> (ri, ...) real buffer
  function c_to_ri(c) result(r)
    implicit none
    complex, dimension(:), intent(in) :: c
    real(real64), dimension(2*size(c)) :: r
    integer :: i
    do i = 1, size(c)
       r(2*i-1) = real(c(i), real64)
       r(2*i) = aimag(c(i))
    end do
  end function c_to_ri

  !-------------------------------------------------------------------------- V arrays

  !> Real array, rank 1 (dims: Fortran order, space separated)
  subroutine ex_vec_r1(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, real(arr, real64), dims, 'real', tagv(tag), fname, src, owner, dtdep)
  end subroutine ex_vec_r1

  !> Real array, rank 2
  subroutine ex_vec_r2(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, real(reshape(arr, [size(arr)]), real64), dims, 'real', tagv(tag), &
         fname, src, owner, dtdep)
  end subroutine ex_vec_r2

  !> Real array, rank 3
  subroutine ex_vec_r3(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, real(reshape(arr, [size(arr)]), real64), dims, 'real', tagv(tag), &
         fname, src, owner, dtdep)
  end subroutine ex_vec_r3

  !> Real array, rank 4
  subroutine ex_vec_r4(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    real, dimension(:,:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, real(reshape(arr, [size(arr)]), real64), dims, 'real', tagv(tag), &
         fname, src, owner, dtdep)
  end subroutine ex_vec_r4

  !> Integer array, rank 1
  subroutine ex_vec_i1(gid, name, arr, dims, fname, src, owner, ibase, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    integer, intent(in), optional :: ibase
    if (.not. proc0) return
    call put_flat_i(gid, name, int(arr, int32), dims, 'integer', tagv(tag), fname, src, owner, ibase)
  end subroutine ex_vec_i1

  !> Integer array, rank 2
  subroutine ex_vec_i2(gid, name, arr, dims, fname, src, owner, ibase, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    integer, intent(in), optional :: ibase
    if (.not. proc0) return
    call put_flat_i(gid, name, int(reshape(arr, [size(arr)]), int32), dims, 'integer', tagv(tag), &
         fname, src, owner, ibase)
  end subroutine ex_vec_i2

  !> Integer array, rank 3
  subroutine ex_vec_i3(gid, name, arr, dims, fname, src, owner, ibase, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    integer, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    integer, intent(in), optional :: ibase
    if (.not. proc0) return
    call put_flat_i(gid, name, int(reshape(arr, [size(arr)]), int32), dims, 'integer', tagv(tag), &
         fname, src, owner, ibase)
  end subroutine ex_vec_i3

  !> Logical array, rank 1 (int8)
  subroutine ex_vec_l1(gid, name, arr, dims, fname, src, owner, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    logical, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    if (.not. proc0) return
    call put_flat_b(gid, name, merge(1_int8, 0_int8, arr), dims, tagv(tag), fname, src, owner)
  end subroutine ex_vec_l1

  !> Logical array, rank 2 (int8)
  subroutine ex_vec_l2(gid, name, arr, dims, fname, src, owner, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    logical, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    if (.not. proc0) return
    call put_flat_b(gid, name, merge(1_int8, 0_int8, reshape(arr, [size(arr)])), dims, tagv(tag), &
         fname, src, owner)
  end subroutine ex_vec_l2

  !> Logical array, rank 3 (int8)
  subroutine ex_vec_l3(gid, name, arr, dims, fname, src, owner, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    logical, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    if (.not. proc0) return
    call put_flat_b(gid, name, merge(1_int8, 0_int8, reshape(arr, [size(arr)])), dims, tagv(tag), &
         fname, src, owner)
  end subroutine ex_vec_l3

  !> Complex array, rank 1: dims must start with 'ri'
  subroutine ex_vec_c1(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    complex, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, c_to_ri(arr), dims, 'complex', tagv(tag), fname, src, owner, dtdep)
  end subroutine ex_vec_c1

  !> Complex array, rank 2: dims must start with 'ri'
  subroutine ex_vec_c2(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    complex, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, c_to_ri(reshape(arr, [size(arr)])), dims, 'complex', tagv(tag), &
         fname, src, owner, dtdep)
  end subroutine ex_vec_c2

  !> Complex array, rank 3: dims must start with 'ri'
  subroutine ex_vec_c3(gid, name, arr, dims, fname, src, owner, dtdep, tag)
    use mp, only: proc0
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, dims
    complex, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner, tag
    logical, intent(in), optional :: dtdep
    if (.not. proc0) return
    call put_flat_r(gid, name, c_to_ri(reshape(arr, [size(arr)])), dims, 'complex', tagv(tag), &
         fname, src, owner, dtdep)
  end subroutine ex_vec_c3

  !> Default canonical tag of V arrays
  function tagv(tag)
    implicit none
    character(len=*), intent(in), optional :: tag
    character(len=:), allocatable :: tagv
    if (present(tag)) then
       tagv = tag
    else
       tagv = 'V'
    end if
  end function tagv

  !-------------------------------------------------------------------------- g_lo canonicalisers

  !> Check the shape of a g_lo-distributed argument
  subroutine check_glo_arg(name, n3, n1, n1_expected, n2, n2_expected)
    use gs2_layouts, only: g_lo
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: n3, n1, n1_expected, n2, n2_expected
    if (n1 /= n1_expected .or. n2 /= n2_expected .or. n3 < g_lo%ulim_proc - g_lo%llim_proc + 1) &
         call ex_fail('unexpected shape of g_lo array '//trim(name))
  end subroutine check_glo_arg

  !> G7 complex: arr(1:nth, 1:2, 1:nloc) -> (ri, theta, sign, kx, ky, lambda, egrid, species),
  !> written per species slab (gpu/spec/09_export_infra.md §2.4)
  subroutine ex_g7c(gid, name, arr, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    complex, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    character(len=*), parameter :: dims = 'ri theta sign kx ky lambda egrid species'
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base, k
    integer :: nth, is, ie, il, ik, it, isgn, ig, iloc, varid
    if (.not. proc0) return
    nth = 2*ntgrid + 1
    call check_glo_arg(name, size(arr, 3), size(arr, 1), nth, size(arr, 2), 2)
    nslab = 2_int64*nth*2*ntheta0*naky*nlambda*negrid
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky
                do it = 1, ntheta0
                   iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
                   base = 2_int64*nth*2*((it-1) + int(ntheta0, int64)*((ik-1) + &
                        int(naky, int64)*((il-1) + int(nlambda, int64)*(ie-1))))
                   do isgn = 1, 2
                      do ig = 1, nth
                         k = 1 + base + 2_int64*((ig-1) + nth*(isgn-1))
                         slab(k) = real(arr(ig, isgn, iloc), real64)
                         slab(k+1) = aimag(arr(ig, isgn, iloc))
                      end do
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, 1, 1, is], &
            count=[2, nth, 2, ntheta0, naky, nlambda, negrid, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'complex', 'G7', dims, fname, src, owner, dtdep, &
         decl='complex(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)')
  end subroutine ex_g7c

  !> G7 real: arr(1:nth, 1:2, 1:nloc) -> (theta, sign, kx, ky, lambda, egrid, species)
  subroutine ex_g7r(gid, name, arr, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    real, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    character(len=*), parameter :: dims = 'theta sign kx ky lambda egrid species'
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base
    integer :: nth, is, ie, il, ik, it, isgn, iloc, varid
    if (.not. proc0) return
    nth = 2*ntgrid + 1
    call check_glo_arg(name, size(arr, 3), size(arr, 1), nth, size(arr, 2), 2)
    nslab = int(nth, int64)*2*ntheta0*naky*nlambda*negrid
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky
                do it = 1, ntheta0
                   iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
                   base = int(nth, int64)*2*((it-1) + int(ntheta0, int64)*((ik-1) + &
                        int(naky, int64)*((il-1) + int(nlambda, int64)*(ie-1))))
                   do isgn = 1, 2
                      slab(base+1+nth*(isgn-1):base+nth*isgn) = real(arr(:, isgn, iloc), real64)
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, 1, is], &
            count=[nth, 2, ntheta0, naky, nlambda, negrid, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'real', 'G7', dims, fname, src, owner, dtdep, &
         decl='real(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)')
  end subroutine ex_g7r

  !> G6 real: arr(1:nth, 1:nloc) -> (theta, kx, ky, lambda, egrid, species)
  subroutine ex_g6r(gid, name, arr, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    real, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    character(len=*), parameter :: dims = 'theta kx ky lambda egrid species'
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base
    integer :: nth, is, ie, il, ik, it, iloc, varid
    if (.not. proc0) return
    nth = 2*ntgrid + 1
    call check_glo_arg(name, size(arr, 2), size(arr, 1), nth, 1, 1)
    nslab = int(nth, int64)*ntheta0*naky*nlambda*negrid
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky
                do it = 1, ntheta0
                   iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
                   base = int(nth, int64)*((it-1) + int(ntheta0, int64)*((ik-1) + &
                        int(naky, int64)*((il-1) + int(nlambda, int64)*(ie-1))))
                   slab(base+1:base+nth) = real(arr(:, iloc), real64)
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, is], &
            count=[nth, ntheta0, naky, nlambda, negrid, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'real', 'G6', dims, fname, src, owner, dtdep, &
         decl='real(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc)')
  end subroutine ex_g6r

  !> G7m1 complex: arr(1:2*ntgrid, 1:2, 1:nloc) -> (ri, theta_m1, sign, kx, ky, lambda, egrid, species)
  subroutine ex_g7m1c(gid, name, arr, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    complex, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    character(len=*), parameter :: dims = 'ri theta_m1 sign kx ky lambda egrid species'
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base, k
    integer :: nt, is, ie, il, ik, it, isgn, ig, iloc, varid
    if (.not. proc0) return
    nt = 2*ntgrid
    call check_glo_arg(name, size(arr, 3), size(arr, 1), nt, size(arr, 2), 2)
    nslab = 2_int64*nt*2*ntheta0*naky*nlambda*negrid
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky
                do it = 1, ntheta0
                   iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
                   base = 2_int64*nt*2*((it-1) + int(ntheta0, int64)*((ik-1) + &
                        int(naky, int64)*((il-1) + int(nlambda, int64)*(ie-1))))
                   do isgn = 1, 2
                      do ig = 1, nt
                         k = 1 + base + 2_int64*((ig-1) + nt*(isgn-1))
                         slab(k) = real(arr(ig, isgn, iloc), real64)
                         slab(k+1) = aimag(arr(ig, isgn, iloc))
                      end do
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, 1, 1, is], &
            count=[2, nt, 2, ntheta0, naky, nlambda, negrid, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'complex', 'G7m1', dims, fname, src, owner, dtdep, &
         decl='complex(-ntgrid:ntgrid-1,2,g_lo%llim_proc:g_lo%ulim_alloc)')
  end subroutine ex_g7m1c

  !> GL complex: arr(1:nl, 1:2, 1:nloc) -> (ri, nlinks_max, sign, kx, ky, lambda, egrid, species)
  subroutine ex_glc(gid, name, arr, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    complex, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    character(len=*), parameter :: dims = 'ri nlinks_max sign kx ky lambda egrid species'
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base, k
    integer :: nl, is, ie, il, ik, it, isgn, i, iloc, varid
    if (.not. proc0) return
    nl = size(arr, 1)
    call check_glo_arg(name, size(arr, 3), nl, nl, size(arr, 2), 2)
    nslab = 2_int64*nl*2*ntheta0*naky*nlambda*negrid
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky
                do it = 1, ntheta0
                   iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
                   base = 2_int64*nl*2*((it-1) + int(ntheta0, int64)*((ik-1) + &
                        int(naky, int64)*((il-1) + int(nlambda, int64)*(ie-1))))
                   do isgn = 1, 2
                      do i = 1, nl
                         k = 1 + base + 2_int64*((i-1) + nl*(isgn-1))
                         slab(k) = real(arr(i, isgn, iloc), real64)
                         slab(k+1) = aimag(arr(i, isgn, iloc))
                      end do
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, 1, 1, is], &
            count=[2, nl, 2, ntheta0, naky, nlambda, negrid, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'complex', 'GL', dims, fname, src, owner, dtdep, &
         decl='complex(n_links_max,2,g_lo%llim_proc:g_lo%ulim_alloc)')
  end subroutine ex_glc

  !> Index (it,ik,il,ie,is) -> position in a G5 canonical buffer (1-based)
  integer(int64) function g5_pos(it, ik, il, ie, is)
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    implicit none
    integer, intent(in) :: it, ik, il, ie, is
    g5_pos = 1 + (it-1) + int(ntheta0, int64)*((ik-1) + int(naky, int64)*((il-1) + &
         int(nlambda, int64)*((ie-1) + int(negrid, int64)*(is-1))))
  end function g5_pos

  !> G5 integer: arr(1:nloc) -> (kx, ky, lambda, egrid, species)
  subroutine ex_g5i(gid, name, arr, fname, src, owner, ibase)
    use mp, only: proc0
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    integer, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer, intent(in), optional :: ibase
    integer(int32), dimension(:), allocatable :: buf
    integer :: is, ie, il, ik, it
    if (.not. proc0) return
    allocate (buf(int(ntheta0, int64)*naky*nlambda*negrid*nspec))
    do is = 1, nspec ; do ie = 1, negrid ; do il = 1, nlambda ; do ik = 1, naky ; do it = 1, ntheta0
       buf(g5_pos(it, ik, il, ie, is)) = int(arr(idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1), int32)
    end do ; end do ; end do ; end do ; end do
    call put_flat_i(gid, name, buf, 'kx ky lambda egrid species', 'integer', 'G5', fname, src, owner, ibase)
  end subroutine ex_g5i

  !> G5 logical: arr(1:nloc) -> int8 (kx, ky, lambda, egrid, species)
  subroutine ex_g5l(gid, name, arr, fname, src, owner)
    use mp, only: proc0
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    logical, dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer(int8), dimension(:), allocatable :: buf
    integer :: is, ie, il, ik, it
    if (.not. proc0) return
    allocate (buf(int(ntheta0, int64)*naky*nlambda*negrid*nspec))
    do is = 1, nspec ; do ie = 1, negrid ; do il = 1, nlambda ; do ik = 1, naky ; do it = 1, ntheta0
       buf(g5_pos(it, ik, il, ie, is)) = merge(1_int8, 0_int8, arr(idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1))
    end do ; end do ; end do ; end do ; end do
    call put_flat_b(gid, name, buf, 'kx ky lambda egrid species', 'G5', fname, src, owner)
  end subroutine ex_g5l

  !> G5s logical: arr(1:2, 1:nloc) -> int8 (sign, kx, ky, lambda, egrid, species)
  subroutine ex_g5sl(gid, name, arr, fname, src, owner)
    use mp, only: proc0
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: g_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name
    logical, dimension(:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    integer(int8), dimension(:), allocatable :: buf
    integer(int64) :: k
    integer :: is, ie, il, ik, it, iloc
    if (.not. proc0) return
    allocate (buf(2_int64*ntheta0*naky*nlambda*negrid*nspec))
    do is = 1, nspec ; do ie = 1, negrid ; do il = 1, nlambda ; do ik = 1, naky ; do it = 1, ntheta0
       iloc = idx(g_lo, ik, it, il, ie, is) - g_lo%llim_proc + 1
       k = 2*(g5_pos(it, ik, il, ie, is) - 1)
       buf(k+1) = merge(1_int8, 0_int8, arr(1, iloc))
       buf(k+2) = merge(1_int8, 0_int8, arr(2, iloc))
    end do ; end do ; end do ; end do ; end do
    call put_flat_b(gid, name, buf, 'sign kx ky lambda egrid species', 'G5s', fname, src, owner)
  end subroutine ex_g5sl

  !-------------------------------------------------------------------------- le_lo canonicalisers

  !> LE6 real: arr(1:NX, 1:NE, 1:nloc) -> (xdim, edim, theta, kx, ky, species), per species slab
  !> (gpu/spec/06_collisions.md §2.4)
  subroutine ex_le6r(gid, name, arr, xdim, edim, fname, src, owner, dtdep)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: le_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, xdim, edim
    real, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    logical, intent(in), optional :: dtdep
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base, nxe
    integer :: nth, nx, ne, is, ik, it, ig, iloc, varid, ix, ie
    character(len=128) :: dims
    if (.not. proc0) return
    nth = 2*ntgrid + 1
    nx = size(arr, 1) ; ne = size(arr, 2)
    if (nx /= ex_dim_len(gid, xdim) .or. ne /= ex_dim_len(gid, edim) .or. &
         size(arr, 3) < le_lo%ulim_proc - le_lo%llim_proc + 1) &
         call ex_fail('unexpected shape of le_lo array '//trim(name))
    dims = trim(xdim)//' '//trim(edim)//' theta kx ky species'
    nxe = int(nx, int64)*ne
    nslab = nxe*nth*ntheta0*naky
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                iloc = idx(le_lo, ig, ik, it, is) - le_lo%llim_proc + 1
                base = nxe*((ig+ntgrid) + int(nth, int64)*((it-1) + int(ntheta0, int64)*(ik-1)))
                do ie = 1, ne
                   do ix = 1, nx
                      slab(base + ix + int(nx, int64)*(ie-1)) = real(arr(ix, ie, iloc), real64)
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, is], &
            count=[nx, ne, nth, ntheta0, naky, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'real', 'LE6', dims, fname, src, owner, dtdep, &
         decl='real(nx,ne,le_lo%llim_proc:le_lo%ulim_alloc)')
  end subroutine ex_le6r

  !> LE6 complex: arr(1:NX, 1:NE, 1:nloc) -> (ri, xdim, edim, theta, kx, ky, species)
  subroutine ex_le6c(gid, name, arr, xdim, edim, fname, src, owner)
    use netcdf, only: nf90_put_var, NF90_DOUBLE
    use mp, only: proc0
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: le_lo, idx
    implicit none
    integer, intent(in) :: gid
    character(len=*), intent(in) :: name, xdim, edim
    complex, dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in), optional :: fname, src, owner
    real(real64), dimension(:), allocatable :: slab
    integer(int64) :: nslab, base, nxe, k
    integer :: nth, nx, ne, is, ik, it, ig, iloc, varid, ix, ie
    character(len=128) :: dims
    if (.not. proc0) return
    nth = 2*ntgrid + 1
    nx = size(arr, 1) ; ne = size(arr, 2)
    if (nx /= ex_dim_len(gid, xdim) .or. ne /= ex_dim_len(gid, edim) .or. &
         size(arr, 3) < le_lo%ulim_proc - le_lo%llim_proc + 1) &
         call ex_fail('unexpected shape of le_lo array '//trim(name))
    dims = 'ri '//trim(xdim)//' '//trim(edim)//' theta kx ky species'
    nxe = int(nx, int64)*ne
    nslab = 2*nxe*nth*ntheta0*naky
    allocate (slab(nslab))
    varid = def_var(gid, name, NF90_DOUBLE, dims)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                iloc = idx(le_lo, ig, ik, it, is) - le_lo%llim_proc + 1
                base = nxe*((ig+ntgrid) + int(nth, int64)*((it-1) + int(ntheta0, int64)*(ik-1)))
                do ie = 1, ne
                   do ix = 1, nx
                      k = 2*(base + ix + int(nx, int64)*(ie-1)) - 1
                      slab(k) = real(arr(ix, ie, iloc), real64)
                      slab(k+1) = aimag(arr(ix, ie, iloc))
                   end do
                end do
             end do
          end do
       end do
       call check(nf90_put_var(gid, varid, slab, start=[1, 1, 1, 1, 1, 1, is], &
            count=[2, nx, ne, nth, ntheta0, naky, 1]), 'put_var '//trim(name))
    end do
    call set_meta(gid, varid, 'complex', 'LE6', dims, fname, src, owner, &
         decl='complex(nx,ne,le_lo%llim_proc:le_lo%ulim_alloc)')
  end subroutine ex_le6c
end module gs2_export_writer
