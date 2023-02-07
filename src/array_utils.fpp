module array_utils
  implicit none

  private

  public :: copy, gs2_maxval

  !> Interface to helper subroutines for copying one array into another.
  interface copy
     module procedure copy_complex_array_3
  end interface copy

  !> Interface to helper functions for finding the maximum in an array
  interface gs2_maxval
     module procedure gs2_maxval_real_array_2
     module procedure gs2_maxval_real_array_3
  end interface gs2_maxval

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

  !> Find maximum in 2D array
  real function gs2_maxval_real_array_2(array_in) result(the_max)
    implicit none
    real, dimension(:, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(dynamic)
    do counter = start, end
       the_max = max(maxval(array_in(:, counter)), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_maxval_real_array_2

  !> Find maximum in 3D array
  real function gs2_maxval_real_array_3(array_in) result(the_max)
    implicit none
    real, dimension(:, :, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(dynamic)
    do counter = start, end
       the_max = max(maxval(array_in(:, :, counter)), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_maxval_real_array_3

end module array_utils
