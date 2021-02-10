program mat_inv_benchmark
  use mp, only: init_mp, finish_mp, barrier, proc0, iproc
  use mat_inv, only: inverse_gj
  !use mat_inv, only: inverse_lapack
  use ran, only: ranf
  use job_manage, only: timer_local
  use iso_fortran_env, only: output_unit
  complex, dimension(:, :), allocatable :: array
  integer, dimension(:), allocatable :: sizes
  integer :: icase, nrep, current_size, ncases
  real, dimension(:), allocatable :: current_times
  real, dimension(:), allocatable :: average_times, max_times, min_times
  real :: time_a, time_b
  integer :: output

  call init_mp

  output = output_unit

  if (proc0) then
     sizes =  [2,4,8,16,32,64,256,512,1024,2048,4096,8192]
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
              array(j,k) = cmplx(ranf(), ranf())
           end do
        end do
        time_b = timer_local() - time_a
        write(output, '("Took ",F12.5,"s to construct with size ",I0," for ",I0," repeats")') time_b, current_size, nrep

        allocate(current_times(nrep))
        current_times = 0

        do irep = 1, nrep
           time_a = timer_local()

           call inverse_gj(array, current_size)
           !call inverse_lapack(array, current_size)

           time_b = timer_local()
           current_times(irep) = time_b - time_a
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
end program mat_inv_benchmark
