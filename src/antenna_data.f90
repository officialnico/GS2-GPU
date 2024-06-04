!> FIXME : Add documentation
module antenna_data

! quick fix for now; put amplitudes here and thus be 
! certain that there are no weird module dependencies
  implicit none

  private

  public :: init_antenna_data, finish_antenna_data
  public :: ant_on, a_ant, b_ant
  protected :: ant_on

  complex, dimension(:), allocatable :: a_ant, b_ant
  logical :: ant_on = .false., initialized = .false.

contains
  !> FIXME : Add documentation
  subroutine init_antenna_data (nk_stir)
    integer, intent (in) :: nk_stir

! do not reallocate this array on this processor...

    if (initialized) return
    initialized = .true.

    ant_on = nk_stir > 0
    if (.not. ant_on) return
    allocate (a_ant(nk_stir), b_ant(nk_stir))
  end subroutine init_antenna_data

  !> FIXME : Add documentation
  subroutine finish_antenna_data
    implicit none
    ant_on = .false.
    !CMR,17/7/2009: add allocated test to avoid severe runtime library error
    if (allocated(a_ant) .and. allocated(b_ant)) deallocate (a_ant, b_ant)
    initialized = .false.
  end subroutine finish_antenna_data
end module antenna_data
