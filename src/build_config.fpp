!> Compile-time parameters describing the built configuration
module build_config
  implicit none

#ifdef MPI
  logical, parameter :: gs2_has_mpi = .true.
#else
  logical, parameter :: gs2_has_mpi = .false.
#endif

#ifdef MPI3
  logical, parameter :: gs2_has_mpi3 = .true.
#else
  logical, parameter :: gs2_has_mpi3 = .false.
#endif

#ifdef FFT
  logical, parameter :: gs2_has_fft = .true.
#else
  logical, parameter :: gs2_has_fft = .false.
#endif

#ifdef ISO_C_BINDING
  logical, parameter :: gs2_has_iso_c_binding = .true.
#else
  logical, parameter :: gs2_has_iso_c_binding = .false.
#endif

#ifdef NEW_DIAG
  logical, parameter :: gs2_has_new_diag = .true.
#else
  logical, parameter :: gs2_has_new_diag = .false.
#endif

#ifdef GENQUAD
  logical, parameter :: gs2_has_genquad = .true.
#else
  logical, parameter :: gs2_has_genquad = .false.
#endif

#ifdef NETCDF
  logical, parameter :: gs2_has_netcdf = .true.
#else
  logical, parameter :: gs2_has_netcdf = .false.
#endif

#ifdef LAPACK
  logical, parameter :: gs2_has_lapack = .true.
#else
  logical, parameter :: gs2_has_lapack = .false.
#endif

#ifdef LOWFLOW
  logical, parameter :: gs2_has_lowflow = .true.
#else
  logical, parameter :: gs2_has_lowflow = .false.
#endif

#ifdef WITH_EIG
  logical, parameter :: gs2_has_with_eig = .true.
#else
  logical, parameter :: gs2_has_with_eig = .false.
#endif

#ifdef HDF
  logical, parameter :: gs2_has_hdf5 = .true.
#else
  logical, parameter :: gs2_has_hdf5 = .false.
#endif

#ifdef F200X_INTRINSICS
  logical, parameter :: gs2_has_f200x_intrinsics = .true.
#else
  logical, parameter :: gs2_has_f200x_intrinsics = .false.
#endif

#ifdef SINGLE_PRECISION
  logical, parameter :: gs2_is_single_precision = .true.
#else
  logical, parameter :: gs2_is_single_precision = .false.
#endif

#ifdef QUAD
  logical, parameter :: gs2_has_quad = .true.
#else
  logical, parameter :: gs2_has_quad = .false.
#endif

contains

  function formatted_build_config() result(output)
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
    character(len=:), allocatable :: output

    character(len=*), parameter :: nl = new_line('a')
    integer, parameter :: max_len = 2000

    allocate(character(len=max_len)::output)

    write(output, '(a, a)') "Built with ", compiler_version()
    write(output, '(a, a)') trim(output) // nl // "Build flags: ", compiler_options()
    write(output, '(a)') trim(output) // nl // nl // "Current build configuration: "
    write(output, '(a, l)') trim(output) // nl // "MPI: ", gs2_has_mpi
    write(output, '(a, l)') trim(output) // nl // "MPI3: ", gs2_has_mpi3
    write(output, '(a, l)') trim(output) // nl // "FFTW: ", gs2_has_fft
    write(output, '(a, l)') trim(output) // nl // "iso_c_binding: ", gs2_has_iso_c_binding
    write(output, '(a, l)') trim(output) // nl // "New diagnostics: ", gs2_has_new_diag
    write(output, '(a, l)') trim(output) // nl // "GENQUAD: ", gs2_has_genquad
    write(output, '(a, l)') trim(output) // nl // "NetCDF: ", gs2_has_netcdf
    write(output, '(a, l)') trim(output) // nl // "LAPACK: ", gs2_has_lapack
    write(output, '(a, l)') trim(output) // nl // "Lowflow: ", gs2_has_lowflow
    write(output, '(a, l)') trim(output) // nl // "Eigensolver: ", gs2_has_with_eig
    write(output, '(a, l)') trim(output) // nl // "HDF5: ", gs2_has_hdf5
    write(output, '(a, l)') trim(output) // nl // "F200X intrinsics: ", gs2_has_f200x_intrinsics
    write(output, '(a, l)') trim(output) // nl // "Single precision: ", gs2_is_single_precision
    write(output, '(a, l)') trim(output) // nl // "Quad precision: ", gs2_has_quad

    output = trim(output)
  end function formatted_build_config

end module build_config
