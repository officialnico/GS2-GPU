!> Provide some small utilty functions for trying to
!> monitor memory use.
module memory_usage
  use iso_fortran_env, only: int64

  implicit none

  private

  public :: print_memory_usage
  public :: number_of_bits_in_use, peak_number_of_bits_in_use
  public :: bits_to_human_readable

contains

  !> Writes current estimate of memory use and peak memory use
  !> to provided unit.
  subroutine print_memory_usage(show_peak, show_current, unit)
    use iso_fortran_env, only: output_unit
    use optionals, only: get_option_with_default
    implicit none
    logical, intent(in), optional :: show_peak, show_current
    integer, intent(in), optional :: unit
    logical :: peak, current
    integer :: output
    integer(int64) :: bits, peak_bits
    ! Handle optionals
    peak = get_option_with_default(show_peak, .true.)
    current = get_option_with_default(show_current, .true.)
    output = get_option_with_default(unit, output_unit)

    ! Leave early if possible
    if (.not. (peak.or.current)) return

    write(output, '(35("-"))')

    if (current) then
       bits = number_of_bits_in_use()
       if ( bits >= 0 ) then
          write(output, '(" Current memory usage  : ", A)') &
               trim(adjustl(bits_to_human_readable(bits)))
       else
          write(output, '(" Cannot estimate memory usage.")')
       end if
    end if

    if (peak) then
       peak_bits = peak_number_of_bits_in_use()
       if ( peak_bits >= 0 ) then
          write(output, '(" Peak memory usage     : ", A)') &
               trim(adjustl(bits_to_human_readable(peak_bits)))
       else
          write(output, '(" Cannot estimate peak memory usage.")')
       end if
    end if

    write(output, '(35("-"))')
  end subroutine print_memory_usage

  !> Internal routine which can be used for debugging. Tries to dump the contents
  !> of /proc/self/status to output_unit
  subroutine print_status()
    use iso_fortran_env, only: output_unit
    implicit none
    character(len=:), allocatable :: string
    integer :: unit, i, err

    open(action = "read", file = '/proc/self/status', newunit = unit, iostat = err)
    if (err /= 0) then
       return
    end if

    allocate(character(len=60)::string)
    do while(err == 0)
       read(unit, '(A)', iostat = err) string
       if (err == 0 ) write(output_unit, '(A)') trim(string)
    end do
    close(unit)
  end subroutine print_status

  !> Try to read the number of kilobytes in use currently
  !> from /proc/self/status and return number of bits in use.
  !> Will return negative if error encountered.
  integer(int64) function number_of_bits_in_use() result(in_use)
    implicit none
    ! We multiply by 1024 * 8 to convert to bits from kilobytes
    in_use = get_status_entry(status_entry = 18) * 1024 * 8
  end function number_of_bits_in_use

  !> Try to read the peak number of kilobytes currently
  !> from /proc/self/status and return number of bits in use.
  !> Will return negative if error encountered.
  integer(int64) function peak_number_of_bits_in_use() result(in_use)
    implicit none
    ! We multiply by 1024 * 8 to convert to bits from kilobytes
    in_use = get_status_entry(status_entry = 17) * 1024 * 8
  end function peak_number_of_bits_in_use

  !> Try to read an integer from /proc/self/status. Returns negative
  !> if not possible to open file.
  integer(int64) function get_status_entry(status_entry) result(val)
    implicit none
    integer, intent(in) :: status_entry
    character(len=:), allocatable :: string
    integer :: unit, i, err
    open(action = "read", file = '/proc/self/status', newunit = unit, iostat = err)
    if (err /= 0) then
       val = -1
       return
    end if

    allocate(character(len=60)::string)
    do i = 1, status_entry
       read(unit, '(A)') string
    end do
    close(unit)

    i = scan(string, ":")
    string = string(i+1:)
    i = scan(string, "k")
    read(string(:i-1), *) val
  end function get_status_entry

  !> Converts the passed number of bits into a slightly more human
  !> readable form.
  pure function bits_to_human_readable(bits) result(human_form)
    use iso_fortran_env, only: real64
    implicit none
    integer(kind = int64), intent(in) :: bits
    character(len = 16) :: human_form
    ! Note the following are really KiB, MiB, GiB etc. as measured in units of 1024
    ! but we drop this here to avoid complexity
    character(len = 2), dimension(*), parameter :: extensions = [" B", "KB", "MB", "GB", "TB", "PB", "EB"]
    real(kind = real64) :: bytes, factor
    integer :: extension_index

    ! Special handling for exactly 0
    if (bits .eq. 0) then
       human_form = "0.00 B"
       return
    end if

    ! We round anything less than a byte up to 1 byte
    bytes = max(bits/8.0_real64,1.0_real64)

    ! Find which range we're in
    extension_index = floor(log(bytes)/log(1024.0_real64))
    ! Limit range to known extensions
    extension_index = min(extension_index, size(extensions))

    ! Put bytes into right range
    factor = 1024.0_real64**(extension_index)
    bytes = bytes / factor

    ! Form the message
    write(human_form,'(F7.2," ",A)') bytes, trim(adjustl(extensions(extension_index+1)))

    ! Remove leading and trailing whitespace
    write(human_form,'(A)') trim(adjustl(trim(human_form)))
  end function bits_to_human_readable
end module memory_usage
