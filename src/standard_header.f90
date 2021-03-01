!> Contains functions for creating a standard header for files and output
module standard_header_mod
  implicit none

  ! Explicit interfaces for the functions in git_hash.fpp
  interface
    function get_git_hash(length_in)
      integer, optional, intent(in) :: length_in
      integer :: length
      character(:), allocatable :: get_git_hash
    end function get_git_hash

    function get_git_state()
      character(:), allocatable :: get_git_state
    end function get_git_state
  end interface

contains

  !> Return the UUID for this run
  !>
  !> Generated when first called, guaranteed not to change during the
  !> simulation. Communicates when first called to broadcast the UUID
  !> to all processors
  function simulation_run_uuid()
    use mp, only : proc0, broadcast
    use uuid_mod, only : generate_uuid, uuid_len
    character(len=uuid_len) :: simulation_run_uuid
    character(len=uuid_len), save :: run_uuid
    logical, save :: first_run = .true.

    if (first_run) then
      if (proc0) then
        run_uuid = generate_uuid()
      end if
      call broadcast(run_uuid)
      first_run = .false.
    end if

    simulation_run_uuid = run_uuid
  end function simulation_run_uuid

  !> Return the current date and time in ISO8601 format:
  !>     YYYY-MM-DDThh:mm:ss.ssssZhh:mm
  function date_iso8601()
    character(:), allocatable :: date_iso8601
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    call date_and_time(date, time, zone)

    date_iso8601 = date(1:4) // "-" // date(5:6) // "-" // date (7:8) &
         // "T" // time(1:2) // ":" // time(3:4) // ":" // time(5:10) &
         // "Z" // zone(1:3) // ":" // zone(4:5)
  end function date_iso8601

  !> Return a multiline string with a standard header of the form:
  !>
  !>     Created by GS2 at 2021-02-02T15:01:26.370Z+00:00
  !>     Run UUID: <unique ID for this run>
  !>     GK_SYSTEM: local
  !>     Compiler: gfortran
  !>     <optional file description>
  !>
  !> If \p comment_character is passed, it is prepended to the start
  !> of each line. Note that it is used as-is, including any whitespace.
  function standard_header(comment_character, file_description) result(header)
    use runtime_tests, only: runtime_info
    character(:), allocatable :: header
    !> Optional character(s) to start each line. Whitespace is not
    !> stripped or added. Default is empty string
    character(*), optional, intent(in) :: comment_character
    !> Optional file description. Treated as a single line, that is,
    !> new lines in this string will not begin with
    !> [[comment_character]]
    character(*), optional, intent(in) :: file_description

    ! Indicates a dirty git repo, may be empty
    character(:), allocatable :: git_dirty
    ! Actual comment character to use, may be empty
    character(:), allocatable :: comment_char
    ! Note that description_line will contain the comment character if
    ! present
    character(:), allocatable :: description_line
    ! A literal new line `\n`, because Fortran
    character(*), parameter :: newline = new_line('a')

    if (present(comment_character)) then
      comment_char = comment_character
    else
      comment_char = ""
    end if

    if (present(file_description)) then
      ! If we add more lines to the header, we should add a newline here
      description_line = comment_char // file_description
    else
      description_line = ""
    end if

    header = comment_char // "Created by GS2 at " // date_iso8601() // newline &
         // comment_char // "Run UUID: " // simulation_run_uuid() // newline &
         // comment_char // "GK_SYSTEM: " // trim(runtime_info%get_gk_system()) // newline &
         // comment_char // "Compiler: " // trim(runtime_info%get_compiler_name()) // newline &
         // comment_char // "Commit hash: " // get_git_hash() // get_git_state() // newline &
         // description_line

  end function standard_header
end module standard_header_mod
