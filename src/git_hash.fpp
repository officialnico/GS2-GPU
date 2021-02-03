! This file is deliberately not in a module, so please don't add one! We want to
! recompile this file when the git hash or state changes, but we don't want to
! recompile the entire project when that happens. If this were in a module, then
! anything that `use`d it would pick up a dependency and therefore need
! rebuilding on new commits. Instead, we just rebuild this file and link the
! object in at the link stage -- so the executables get relinked, but that's it.
! standard_header.f90 has explicit interfaces for these functions, so please
! `use` that module in order to call these function.

#ifndef GIT_HASH
#define GIT_HASH "unknown"
#endif

#ifndef GIT_STATE
#define GIT_STATE "unknown"
#endif


!> Returns the git hash
function get_git_hash(length_in)
  implicit none
  !> How many characters of the hash to return, default is 7
  integer, optional, intent(in) :: length_in
  integer :: length
  character(:), allocatable :: get_git_hash

  length = 7
  if (present(length_in)) then
    if (length_in <= 40) then
      length = length_in
    end if
  end if

  allocate(character(length)::get_git_hash)
  get_git_hash = GIT_HASH(1:length)
end function get_git_hash

!> Return "-dirty" if the repository has modifications to tracked files,
!> or the empty string otherwise
function get_git_state()
  implicit none
  character(:), allocatable :: get_git_state

  if (GIT_STATE == "clean") then
    get_git_state = ""
  else
    get_git_state = "-dirty"
  endif
end function get_git_state
