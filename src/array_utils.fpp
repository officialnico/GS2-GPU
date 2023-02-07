module array_utils
  implicit none

  private

  public :: copy

  !> Interface to helper subroutines for copying one array into another.
  interface copy
     module procedure copy_complex_array_3
  end interface copy

contains

  !> Copy 3D complex arrays
  subroutine copy_complex_array_3(array_in, array_out)
    implicit none
    complex, dimension(:, :, :), intent(in) :: array_in
    complex, dimension(:, :, :), intent(out) :: array_out
    integer :: counter, start, end
    ! In a debug build we might want to add a check here that the arrays have
    ! compatible shape.

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(dynamic)
    do counter = start, end
       array_out(:, :, counter) = array_in(:, :, counter)
    end do
    !$OMP END PARALLEL DO
  end subroutine copy_complex_array_3

end module array_utils
