module array_utils
  implicit none

  private

  public :: copy, gs2_maxval, gs2_abs, gs2_max_abs
  public :: zero_array

  !> Interface to helper subroutines for copying one array into another.
  interface copy
     module procedure copy_complex_array_3
  end interface copy

  !> Interface to helper functions for finding the maximum in an array
  interface gs2_maxval
     module procedure gs2_maxval_real_array_2
     module procedure gs2_maxval_real_array_3
  end interface gs2_maxval

  !> Interface to helper functions for finding magnitude of an array
  interface gs2_abs
     module procedure gs2_abs_real_array_2
     module procedure gs2_abs_real_array_3
     module procedure gs2_abs_complex_array_2
     module procedure gs2_abs_complex_array_3
  end interface gs2_abs

  !> Interface to helper functions for finding maxval of magnitude of an array
  interface gs2_max_abs
     module procedure gs2_max_abs_real_array_2
     module procedure gs2_max_abs_real_array_3
     module procedure gs2_max_abs_complex_array_2
     module procedure gs2_max_abs_complex_array_3
  end interface gs2_max_abs

  !> Interface to helper functions to zero an array
  interface zero_array
     module procedure zero_array_real_array_2
     module procedure zero_array_real_array_3
     module procedure zero_array_complex_array_2
     module procedure zero_array_complex_array_3
     module procedure zero_array_complex_array_4
  end interface zero_array

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
    !$OMP SCHEDULE(static)
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
    !$OMP SCHEDULE(static)
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
    !$OMP SCHEDULE(static)
    do counter = start, end
       the_max = max(maxval(array_in(:, :, counter)), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_maxval_real_array_3

  !> Return absolute of 2D array
  function gs2_abs_real_array_2(array_in) result(array_out)
    implicit none
    real, dimension(:, :), intent(in) :: array_in
    real, dimension(:, :), allocatable :: array_out
    integer :: counter, start, end

    allocate(array_out, mold = array_in)
    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_out(:, counter) = abs(array_in(:, counter))
    end do
    !$OMP END PARALLEL DO
  end function gs2_abs_real_array_2

  !> Return absolute of 3D array
  function gs2_abs_real_array_3(array_in) result(array_out)
    implicit none
    real, dimension(:, :, :), intent(in) :: array_in
    real, dimension(:, :, :), allocatable :: array_out
    integer :: counter, start, end

    allocate(array_out, mold = array_in)
    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_out(:, :, counter) = abs(array_in(:, :, counter))
    end do
    !$OMP END PARALLEL DO
  end function gs2_abs_real_array_3

  !> Return absolute of 2D array
  function gs2_abs_complex_array_2(array_in) result(array_out)
    implicit none
    complex, dimension(:, :), intent(in) :: array_in
    real, dimension(:, :), allocatable :: array_out
    integer :: counter, start, end

    allocate(array_out(size(array_in, dim=1),size(array_in, dim=2)))
    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_out(:, counter) = abs(array_in(:, counter))
    end do
    !$OMP END PARALLEL DO
  end function gs2_abs_complex_array_2

  !> Return absolute of 3D array
  function gs2_abs_complex_array_3(array_in) result(array_out)
    implicit none
    complex, dimension(:, :, :), intent(in) :: array_in
    real, dimension(:, :, :), allocatable :: array_out
    integer :: counter, start, end

    allocate(array_out(size(array_in, dim=1),size(array_in, dim=2),size(array_in, dim=3)))
    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_out(:, :, counter) = abs(array_in(:, :, counter))
    end do
    !$OMP END PARALLEL DO
  end function gs2_abs_complex_array_3

  !> Find maximum of abs 2D array
  real function gs2_max_abs_real_array_2(array_in) result(the_max)
    implicit none
    real, dimension(:, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       the_max = max(maxval(abs(array_in(:, counter))), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_max_abs_real_array_2

  !> Find maximum of abs 3D array
  real function gs2_max_abs_real_array_3(array_in) result(the_max)
    implicit none
    real, dimension(:, :, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       the_max = max(maxval(abs(array_in(:, :, counter))), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_max_abs_real_array_3

  !> Find maximum of abs 2D array
  real function gs2_max_abs_complex_array_2(array_in) result(the_max)
    implicit none
    complex, dimension(:, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       the_max = max(maxval(abs(array_in(:, counter))), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_max_abs_complex_array_2

  !> Find maximum of abs 3D array
  real function gs2_max_abs_complex_array_3(array_in) result(the_max)
    implicit none
    complex, dimension(:, :, :), intent(in) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)
    the_max = 0.0
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP REDUCTION(max: the_max) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       the_max = max(maxval(abs(array_in(:, :, counter))), the_max)
    end do
    !$OMP END PARALLEL DO
  end function gs2_max_abs_complex_array_3

  !> Zero out a 2D array
  subroutine zero_array_real_array_2(array_in)
    implicit none
    real, dimension(:, :), intent(in out) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_in(:, counter) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine zero_array_real_array_2

  !> Zero out a 3D array
  subroutine zero_array_real_array_3(array_in)
    implicit none
    real, dimension(:, :, :), intent(in out) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_in(:, :, counter) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine zero_array_real_array_3

  !> Zero out a 2D array
  subroutine zero_array_complex_array_2(array_in)
    implicit none
    complex, dimension(:, :), intent(in out) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 2)
    end = ubound(array_in, dim = 2)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_in(:, counter) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine zero_array_complex_array_2

  !> Zero out a 3D array
  subroutine zero_array_complex_array_3(array_in)
    implicit none
    complex, dimension(:, :, :), intent(in out) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 3)
    end = ubound(array_in, dim = 3)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_in(:, :, counter) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine zero_array_complex_array_3

  !> Zero out a 4D array
  subroutine zero_array_complex_array_4(array_in)
    implicit none
    complex, dimension(:, :, :, :), intent(in out) :: array_in
    integer :: counter, start, end

    start = lbound(array_in, dim = 4)
    end = ubound(array_in, dim = 4)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(static)
    do counter = start, end
       array_in(:, :, :, counter) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine zero_array_complex_array_4

end module array_utils
