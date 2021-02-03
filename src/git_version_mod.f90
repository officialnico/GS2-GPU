module git_version_mod
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
end module git_version_mod
