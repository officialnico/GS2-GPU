!> Times the matrix inversion code for a particular method for various sizes.
!> Can also check the correctness.
program mat_inv_benchmark
  use mp, only: init_mp, finish_mp, barrier, proc0, iproc
  use mat_inv, only: invert_serial, matrix_inversion_method_type
  use mat_inv, only: mat_inv_serial_lapack
  use mat_inv, only: mat_inv_serial_gauss_jordan, mat_inv_serial_block
  use ran, only: ranf
  use job_manage, only: timer_local
  use iso_fortran_env, only: output_unit
  implicit none
  complex, dimension(:, :), allocatable :: array, original_array
  integer, dimension(:), allocatable :: sizes
  integer :: icase, irep, nrep, current_size, ncases, j, k
  real, dimension(:), allocatable :: current_times
  real, dimension(:), allocatable :: average_times, max_times, min_times
  real :: time_a, time_b
  integer :: output
  type(matrix_inversion_method_type), parameter :: the_method = mat_inv_serial_gauss_jordan

  call init_mp

  output = output_unit

  if (proc0) then
     sizes =  [2,4,8,16,32,64,128,256,257,400,512,1024,2048]
     ncases = size(sizes)
     allocate(average_times(ncases))
     allocate(max_times(ncases))
     allocate(min_times(ncases))

     do icase = 1, ncases
        current_size = sizes(icase)
        nrep = max(1,sizes(ncases)/current_size)

        allocate( array(current_size, current_size))

        time_a = timer_local()
        do j = 1, current_size
           do k = 1, current_size
              array(j, k) = cmplx(ranf(), ranf())
           end do

           ! Try to ensure array is diagonally dominant
           array(j, j) = array(j, j) + 1.0
        end do

        original_array = array

        time_b = timer_local() - time_a
        write(output, '("Took ",F12.5,"s to construct with size ",I0," for ",I0," repeats")') time_b, current_size, nrep

        allocate(current_times(nrep))
        current_times = 0

        do irep = 1, nrep
           time_a = timer_local()

           call invert_serial(array, current_size, the_method)

           time_b = timer_local()
           current_times(irep) = time_b - time_a

           ! Test the inversion, if it fails abort.
           if(irep==1) then
              if(.not.test(original_array,array, print_error = .true.)) then
                 exit
              end if
           end if

        end do

        average_times(icase) = sum(current_times) / size(current_times)
        max_times(icase) = maxval(current_times)
        min_times(icase) = minval(current_times)
        deallocate(array)
        deallocate(current_times)
     end do

     write(output,'()')
     write(output,'("---------------------------------------------------")')
     write(output,'()')
     write(output,'("Size : Average, max, min")')

     do icase = 1, ncases
        write(output,'(I0," : ",3(F12.5,"s "))') sizes(icase), average_times(icase), max_times(icase), min_times(icase)
     end do

  end if


  call finish_mp

contains

  logical function test(A, B, print_error)
    implicit none
    complex, dimension(:,:), intent(in) :: A, B
    logical, intent(in), optional :: print_error
    complex, dimension(:,:), allocatable :: C
    real, dimension(:,:), allocatable :: identity
    real, parameter :: tolerance = 1.0e-8
    real :: maxerr
    integer :: i, N
    integer, dimension(2) :: error_position
    logical :: print_error_internal

    C = matmul(A,B)
    N = size(A, dim = 1)
    allocate( identity(N, N))
    identity = 0
    do i = 1, N
       identity(i, i) = 1.0
    end do

    maxerr = maxval(abs(C - identity))

    test = maxerr <= tolerance

    print_error_internal = .false.
    if (present(print_error)) print_error_internal = print_error

    if (.not. test .and. print_error) then
       print*,"ERROR: Inversion failed with maxerr = ",maxerr," compared to tolerance ",tolerance
       error_position = maxloc(abs(C - identity))
       print*,"Location of first maximum error is given at (",error_position(1),",",error_position(2),")"
    end if

  end function test

end program mat_inv_benchmark
