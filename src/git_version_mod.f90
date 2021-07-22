!> Explicit interfaces for the functions in git_hash.fpp
module git_version_mod
  implicit none
  interface
    !> Returns the git version from `git describe`
    !>
    !> This looks like: `{tag}-g{hash}[-dirty]`
    function get_git_version()
      character(:), allocatable :: get_git_version
    end function get_git_version

    !> Returns the git hash of the current commit
    function get_git_hash(length_in)
      integer, optional, intent(in) :: length_in
      integer :: length
      character(:), allocatable :: get_git_hash
    end function get_git_hash

    !> Return "-dirty" if the repository has modifications to tracked
    !> files, or the empty string otherwise
    function get_git_state()
      character(:), allocatable :: get_git_state
    end function get_git_state
  end interface
end module git_version_mod
