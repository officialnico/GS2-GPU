! This file is deliberately not in a module, so please don't add one! We want to
! recompile this file when the git hash or state changes, but we don't want to
! recompile the entire project when that happens. If this were in a module, then
! anything that `use`d it would pick up a dependency and therefore need
! rebuilding on new commits. Instead, we just rebuild this file and link the
! object in at the link stage -- so the executables get relinked, but that's it.
! `git_version_mod.f90` has explicit interfaces for these functions, so please
! `use` that module in order to call these functions.

#ifndef GIT_HASH
#define GIT_HASH "unknown"
#endif

#ifndef GIT_STATE
#define GIT_STATE "unknown"
#endif

#ifndef GIT_VERSION
#define GIT_VERSION "unknown"
#endif

!> Returns the git version from `git describe`
!>
!> This looks like: `{tag}-g{hash}[-dirty]`
function get_git_version()
  implicit none
  character(:), allocatable :: get_git_version
  integer, parameter :: max_length = 40
  integer :: length

  length = min(max_length, len(GIT_VERSION))
  allocate(character(length)::get_git_version)
  get_git_version = GIT_VERSION(1:length)
  get_git_version = trim(get_git_version)
end function get_git_version


!> Returns the git hash of the current commit
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
