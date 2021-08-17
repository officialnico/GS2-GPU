
!> Wrapper routines for reading variables from netcdf file.
module simpledataio_read
  use iso_fortran_env, only: real32, real64
  use simpledataio_write, only: real_imaginary_dimension_name
  implicit none
  private
  public :: read_variable_with_offset, read_variable


  interface read_variable_with_offset
    module procedure read_variable_with_offset_real_real32_0
    module procedure read_variable_with_offset_real_real32_1
    module procedure read_variable_with_offset_real_real32_2
    module procedure read_variable_with_offset_real_real32_3
    module procedure read_variable_with_offset_real_real32_4
    module procedure read_variable_with_offset_real_real32_5
    module procedure read_variable_with_offset_real_real32_6
    module procedure read_variable_with_offset_real_real64_0
    module procedure read_variable_with_offset_real_real64_1
    module procedure read_variable_with_offset_real_real64_2
    module procedure read_variable_with_offset_real_real64_3
    module procedure read_variable_with_offset_real_real64_4
    module procedure read_variable_with_offset_real_real64_5
    module procedure read_variable_with_offset_real_real64_6
    module procedure read_variable_with_offset_integer_0
    module procedure read_variable_with_offset_integer_1
    module procedure read_variable_with_offset_integer_2
    module procedure read_variable_with_offset_integer_3
    module procedure read_variable_with_offset_integer_4
    module procedure read_variable_with_offset_integer_5
    module procedure read_variable_with_offset_integer_6
    module procedure read_variable_with_offset_character_0
    module procedure read_variable_with_offset_character_1
    module procedure read_variable_with_offset_character_2
    module procedure read_variable_with_offset_character_3
    module procedure read_variable_with_offset_character_4
    module procedure read_variable_with_offset_character_5
    module procedure read_variable_with_offset_character_6
    module procedure read_variable_with_offset_complex_real32_0
    module procedure read_variable_with_offset_complex_real32_1
    module procedure read_variable_with_offset_complex_real32_2
    module procedure read_variable_with_offset_complex_real32_3
    module procedure read_variable_with_offset_complex_real32_4
    module procedure read_variable_with_offset_complex_real32_5
    module procedure read_variable_with_offset_complex_real32_6
    module procedure read_variable_with_offset_complex_real64_0
    module procedure read_variable_with_offset_complex_real64_1
    module procedure read_variable_with_offset_complex_real64_2
    module procedure read_variable_with_offset_complex_real64_3
    module procedure read_variable_with_offset_complex_real64_4
    module procedure read_variable_with_offset_complex_real64_5
    module procedure read_variable_with_offset_complex_real64_6
  end interface read_variable_with_offset

  interface read_variable
    module procedure read_variable_real_real32_0
    module procedure read_variable_real_real32_1
    module procedure read_variable_real_real32_2
    module procedure read_variable_real_real32_3
    module procedure read_variable_real_real32_4
    module procedure read_variable_real_real32_5
    module procedure read_variable_real_real32_6
    module procedure read_variable_real_real64_0
    module procedure read_variable_real_real64_1
    module procedure read_variable_real_real64_2
    module procedure read_variable_real_real64_3
    module procedure read_variable_real_real64_4
    module procedure read_variable_real_real64_5
    module procedure read_variable_real_real64_6
    module procedure read_variable_integer_0
    module procedure read_variable_integer_1
    module procedure read_variable_integer_2
    module procedure read_variable_integer_3
    module procedure read_variable_integer_4
    module procedure read_variable_integer_5
    module procedure read_variable_integer_6
    module procedure read_variable_character_0
    module procedure read_variable_character_1
    module procedure read_variable_character_2
    module procedure read_variable_character_3
    module procedure read_variable_character_4
    module procedure read_variable_character_5
    module procedure read_variable_character_6
    module procedure read_variable_complex_real32_0
    module procedure read_variable_complex_real32_1
    module procedure read_variable_complex_real32_2
    module procedure read_variable_complex_real32_3
    module procedure read_variable_complex_real32_4
    module procedure read_variable_complex_real32_5
    module procedure read_variable_complex_real32_6
    module procedure read_variable_complex_real64_0
    module procedure read_variable_complex_real64_1
    module procedure read_variable_complex_real64_2
    module procedure read_variable_complex_real64_3
    module procedure read_variable_complex_real64_4
    module procedure read_variable_complex_real64_5
    module procedure read_variable_complex_real64_6
  end interface read_variable
contains

  subroutine read_variable_with_offset_real_real32_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind=real32), dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_0

  subroutine read_variable_with_offset_real_real32_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_1

  subroutine read_variable_with_offset_real_real32_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_2

  subroutine read_variable_with_offset_real_real32_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_3

  subroutine read_variable_with_offset_real_real32_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_4

  subroutine read_variable_with_offset_real_real32_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_5

  subroutine read_variable_with_offset_real_real32_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real32_6

  subroutine read_variable_with_offset_real_real64_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind=real64), dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_0

  subroutine read_variable_with_offset_real_real64_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_1

  subroutine read_variable_with_offset_real_real64_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_2

  subroutine read_variable_with_offset_real_real64_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_3

  subroutine read_variable_with_offset_real_real64_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_4

  subroutine read_variable_with_offset_real_real64_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_5

  subroutine read_variable_with_offset_real_real64_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_real_real64_6

  subroutine read_variable_with_offset_integer_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    integer, dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_0

  subroutine read_variable_with_offset_integer_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_1

  subroutine read_variable_with_offset_integer_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_2

  subroutine read_variable_with_offset_integer_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_3

  subroutine read_variable_with_offset_integer_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_4

  subroutine read_variable_with_offset_integer_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_5

  subroutine read_variable_with_offset_integer_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_integer_6

  subroutine read_variable_with_offset_character_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    character, dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_0

  subroutine read_variable_with_offset_character_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_1

  subroutine read_variable_with_offset_character_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_2

  subroutine read_variable_with_offset_character_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_3

  subroutine read_variable_with_offset_character_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_4

  subroutine read_variable_with_offset_character_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_5

  subroutine read_variable_with_offset_character_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+0
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val(offsets(1):, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_character_6

  subroutine read_variable_with_offset_complex_real32_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1), temp_storage(2))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_0

  subroutine read_variable_with_offset_complex_real32_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :), temp_storage(2, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_1

  subroutine read_variable_with_offset_complex_real32_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :), temp_storage(2, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_2

  subroutine read_variable_with_offset_complex_real32_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :), temp_storage(2, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_3

  subroutine read_variable_with_offset_complex_real32_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :), temp_storage(2, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_4

  subroutine read_variable_with_offset_complex_real32_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :), temp_storage(2, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_5

  subroutine read_variable_with_offset_complex_real32_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):, offsets(7):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :, :), temp_storage(2, :, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real32_6

  subroutine read_variable_with_offset_complex_real64_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+0+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1), temp_storage(2))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_0

  subroutine read_variable_with_offset_complex_real64_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+1+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :), temp_storage(2, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_1

  subroutine read_variable_with_offset_complex_real64_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+2+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :), temp_storage(2, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_2

  subroutine read_variable_with_offset_complex_real64_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+3+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :), temp_storage(2, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_3

  subroutine read_variable_with_offset_complex_real64_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+4+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :), temp_storage(2, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_4

  subroutine read_variable_with_offset_complex_real64_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+5+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :), temp_storage(2, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_5

  subroutine read_variable_with_offset_complex_real64_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+6+1
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage(:, offsets(2):, offsets(3):, offsets(4):, offsets(5):, offsets(6):, offsets(7):), start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :, :), temp_storage(2, :, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_with_offset_complex_real64_6

  subroutine read_variable_real_real32_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind=real32), dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_0

  subroutine read_variable_real_real32_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_1

  subroutine read_variable_real_real32_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_2

  subroutine read_variable_real_real32_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_3

  subroutine read_variable_real_real32_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_4

  subroutine read_variable_real_real32_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_5

  subroutine read_variable_real_real32_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real32), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real32_6

  subroutine read_variable_real_real64_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind=real64), dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_0

  subroutine read_variable_real_real64_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_1

  subroutine read_variable_real_real64_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_2

  subroutine read_variable_real_real64_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_3

  subroutine read_variable_real_real64_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_4

  subroutine read_variable_real_real64_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_5

  subroutine read_variable_real_real64_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    real(kind=real64), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_real_real64_6

  subroutine read_variable_integer_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    integer, dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_0

  subroutine read_variable_integer_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_1

  subroutine read_variable_integer_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_2

  subroutine read_variable_integer_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_3

  subroutine read_variable_integer_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_4

  subroutine read_variable_integer_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_5

  subroutine read_variable_integer_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    integer, dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_integer_6

  subroutine read_variable_character_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    character, dimension(1) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = temp_storage(1)

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_0

  subroutine read_variable_character_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_1

  subroutine read_variable_character_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_2

  subroutine read_variable_character_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_3

  subroutine read_variable_character_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_4

  subroutine read_variable_character_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_5

  subroutine read_variable_character_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    character, dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         val, start=int(starts), count=int(counts))

    

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_character_6

  subroutine read_variable_complex_real32_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1), temp_storage(2))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_0

  subroutine read_variable_complex_real32_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :), temp_storage(2, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_1

  subroutine read_variable_complex_real32_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :), temp_storage(2, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_2

  subroutine read_variable_complex_real32_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :), temp_storage(2, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_3

  subroutine read_variable_complex_real32_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :), temp_storage(2, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_4

  subroutine read_variable_complex_real32_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :), temp_storage(2, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_5

  subroutine read_variable_complex_real32_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real32), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :, :), temp_storage(2, :, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real32_6

  subroutine read_variable_complex_real64_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1), temp_storage(2))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_0

  subroutine read_variable_complex_real64_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :), temp_storage(2, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_1

  subroutine read_variable_complex_real64_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :), temp_storage(2, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_2

  subroutine read_variable_complex_real64_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :), temp_storage(2, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_3

  subroutine read_variable_complex_real64_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :), temp_storage(2, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_4

  subroutine read_variable_complex_real64_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :), temp_storage(2, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_5

  subroutine read_variable_complex_real64_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    complex(kind=real64), dimension(:, :, :, :, :, :), intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    real(kind = kind(val)), dimension(2, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: temp_storage
  
#ifdef FORTRAN_NETCDF
    n2 = number_of_dimensions(sfile, variable_name)
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         temp_storage, start=int(starts), count=int(counts))

    val = cmplx( temp_storage(1, :, :, :, :, :, :), temp_storage(2, :, :, :, :, :, :))

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
  end subroutine read_variable_complex_real64_6
end module simpledataio_read
