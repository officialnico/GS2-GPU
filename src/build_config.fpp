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

end module build_config
