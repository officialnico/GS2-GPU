!> Contains functions for creating a standard header for files and output
module standard_header_mod
  implicit none
contains
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
  !>     Build ID: local.gfortran.b3393a8-dirty
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

    if (runtime_info%get_git_modified()) then
      git_dirty = "-dirty"
    else
      git_dirty = ""
    end if

    header = comment_char // "Created by GS2 at " // date_iso8601() // newline &
         // comment_char // "Build ID: " // trim(runtime_info%build_identifier()) // git_dirty // newline &
         // description_line

  end function standard_header
end module standard_header_mod
