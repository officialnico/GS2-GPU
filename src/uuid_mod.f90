!> Create a random (version 4) Universally Unique Identifier (UUID)
!>
!> A UUID looks like "4042E716-2556-4715-90F0-C6518463B4E5" and is a
!> random 128-bit number. This is sufficiently random that we can
!> except zero collisions in any relevant timeframe, and so is useful
!> for uniquely identifing simulations, for example.
!>
!> This module contains its own implementation of the MT19937
!> psuedo-random number generator (PRNG) so that it doesn't intefere
!> with the PRNG state used in the rest of the program.
module uuid_mod
  use, intrinsic :: iso_c_binding, only: c_signed_char
  use, intrinsic :: iso_fortran_env, only: int32, int8
  implicit none

  private

  public :: generate_uuid
  public :: uuid_len

  !> Length of the UUID character
  integer, parameter :: uuid_len = 36

  integer, parameter :: N_mt19937 = 624
  integer, parameter :: M_mt19937 = 397
  integer, parameter :: MATA = -1727483681 ! constant vector a

  integer, save :: mt(0:N_mt19937 - 1) ! the array for the state vector
  integer, save :: mti = N_mt19937 + 1 ! mti==N+1 means mt[N] is not initialized
  integer, save :: mag01(0:1) = [0, MATA] ! mag01(x) = x * MATA for x=0,1

contains

  !> Convert a 1-byte integer into hexadecimal
  elemental function c_signed_char_to_hex(char) result(hex)
    integer(c_signed_char), intent(in) :: char
    character(len=2) :: hex
    ! Note this almost certainly writes the result in uppercase
    write (hex, '(z2.2)') char
  end function c_signed_char_to_hex

  !> Convert an array of 1-byte integers into UUID format:
  !> "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
  function format_uuid(uuid_in)
    character(len=uuid_len) :: format_uuid
    character(2), dimension(0:15) :: temp1
    character(32) :: temp
    integer(c_signed_char), dimension(16), intent(in) :: uuid_in
    integer :: i

    temp1 = c_signed_char_to_hex(uuid_in)

    do i = 0, 15
      temp(2*i + 1:2*i + 2) = temp1(i)
    end do

    format_uuid = temp(1:8)//"-"//temp(9:12)//"-" &
         //temp(13:16)//"-"//temp(17:20)//"-"//temp(21:)
  end function format_uuid

  !> Generate a version 4 UUID using either the libuuid library, or
  !> our own random number generator wrapper
  function generate_uuid()
    use, intrinsic :: iso_fortran_env, only : int16
    character(len=uuid_len) :: generate_uuid
    integer(c_signed_char), dimension(16) :: uuid_version_4
    integer :: i

    ! Always re-start the PRNG
    call init_uuid_mt19937

    ! Generate random numbers between 0 and 255
    do i = 1, size(uuid_version_4)
      uuid_version_4(i) = get_random_number()
    end do

    ! Version must be 0100xxxx (i.e. 0x4X)
    uuid_version_4(7) = iand(uuid_version_4(7), 79_c_signed_char)
    uuid_version_4(7) = ior(uuid_version_4(7), 64_c_signed_char)

    ! Variant must be 10xxxxxx (i.e 0x8X - 0xBX)
    uuid_version_4(9) = iand(uuid_version_4(9), -65_c_signed_char)
    uuid_version_4(9) = ior(uuid_version_4(9), -127_c_signed_char)

    ! Convert to a lovely string
    generate_uuid = format_uuid(uuid_version_4)
  end function generate_uuid

  !> Initialise the mt19937 generator according to the generator Line
  !> 25 of Table 1 in KNUTH 1981, The Art of Computer Programming
  !> Vol. 2 (2nd Ed.), pp102. Uses the current system time as the seed
  !>
  !> Uses a method for getting a nice seed for the PRNG taken from
  !> [gfortran
  !> documentation](https://gcc.gnu.org/onlinedocs/gcc-6.4.0/gfortran/RANDOM_005fSEED.html)
  subroutine init_uuid_mt19937
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer :: seed

    integer(int64) :: time
    integer :: time_values(8)
    ! Convert from milliseconds to larger units
    integer, parameter :: second = 1000
    integer, parameter :: minute = 60*second
    integer, parameter :: hour = 60*minute
    integer, parameter :: day = 24*hour
    integer, parameter :: month = 31*day
    integer(int64), parameter :: year = 365_int64*day

    integer :: un, istat

    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      ! Fallback to using current system time
      call system_clock(time)
      if (time == 0) then
        ! It's plausible the system_clock isn't POSIX/epoch time so
        ! let's try a different method to get the time
        call date_and_time(values=time_values)
        time = (time_values(1) - 1970)*year &
             + time_values(2)*month &
             + time_values(3)*day &
             + time_values(5)*hour &
             + time_values(6)*minute &
             + time_values(7)*second &
             + time_values(8)
      end if

      ! Ideally here we'd also XOR the PID, which would help when
      ! launching multiple processes at the same time. Note that
      ! `standard_header::simulation_run_uuid` should be used to get a
      ! consistent UUID across MPI ranks, so we don't care about
      ! getting different UUIDs on different ranks

      ! One step of linear congruential generator
      time = mod(time, 4294967296_int64)
      time = mod(time*279470273_int64, 4294967291_int64)
      seed = int(mod(int64, int(huge(0), int64)), kind(seed))
    end if

    mt(0) = iand(seed, -1)
    do mti = 1, N_mt19937 - 1
      mt(mti) = iand(69069*mt(mti - 1), -1)
    end do

  end subroutine init_uuid_mt19937

  !> Generate a random number between 0 and 255 (inclusive)
  function get_random_number()
    implicit none
    integer(c_signed_char) :: get_random_number

    ! Period parameters
    integer, parameter :: LMASK = 2147483647 ! least significant r bits
    integer, parameter :: UMASK = -LMASK - 1 ! most significant w-r bits

    ! Tempering parameters
    integer, parameter :: TMASKB = -1658038656
    integer, parameter :: TMASKC = -272236544

    real, parameter :: pow = 2.**32
    real, parameter :: div = 1./pow       ! divided by 2**32 [0,1)-real-interval
    real :: random_temp
    integer :: y, kk

    if (mti >= N_mt19937) then ! generate N words at one time

      do kk = 0, N_mt19937 - M_mt19937 - 1
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + M_mt19937), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N_mt19937 - M_mt19937, N_mt19937 - 2
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + (M_mt19937 - N_mt19937)), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(mt(N_mt19937 - 1), UMASK), iand(mt(0), LMASK))
      mt(N_mt19937 - 1) = ieor(ieor(mt(M_mt19937 - 1), ishft(y, -1)), mag01(iand(y, 1)))
      mti = 0
    end if

    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    random_temp = real(y)
    if (random_temp < 0.) random_temp = random_temp + pow
    random_temp = random_temp*div

    get_random_number = int(random_temp*huge(get_random_number), kind=c_signed_char)

  end function get_random_number

end module uuid_mod
