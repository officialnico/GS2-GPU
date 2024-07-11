#include "unused_dummy.inc"
!> A module which writes out the electrostatic drive of zonal flows.
!> The code is based on a Python3 post-processing tool from Stephen Biggs-Fox
!> which can be found at [gitlab.com/stephen-biggs-fox/zonal-transfer-functions]
!>
!> The code calculates the Reynolds stress drive through a bispectra-like
!> transfer function by directly performing three-wave-coupling calculations.
!> The couplings take place between a "target mode", a "source modes", and
!> a "mediator mode", the latter being determined by a selection rule.
!>
!> The total electrostatic zonal flow drive is resolved over the poloidal angle
!> in the output "kinetic_energy_transfer_theta". Additional insight into the
!> poloidal distributions of contributions from different source and target scales
!> are provided by "kinetic_energy_transfer_theta_kxsource",
!> "kinetic_energy_transfer_theta_kxtarget", and "kinetic_energy_transfer_theta_kysource".
module diagnostics_kinetic_energy_transfer

  implicit none

  private

  public :: init_diagnostics_kinetic_energy_transfer
  public :: finish_diagnostics_kinetic_energy_transfer
  public :: write_kinetic_energy_transfer
  public :: calculate_kinetic_energy_transfer

  real, dimension(:), allocatable :: kxnew_array, kynew_array
  logical, dimension(:,:,:), allocatable :: valid_mediator_lookup
  real, dimension(:,:,:), allocatable :: z_factor
  integer :: ikx0, iky0, nkynew

  real, dimension(:), allocatable :: T_result_theta
  real, dimension(:,:), allocatable :: T_result_theta_kxsource
  real, dimension(:,:), allocatable :: T_result_theta_kxtarget
  real, dimension(:,:), allocatable :: T_result_theta_kysource

  interface fftshift
     module procedure fftshift1D_real
     module procedure fftshift2D_complex
  end interface fftshift

  !> If true then distribute theta grid points over processors. This requires
  !> some additional global communication so at scale (and with small ntheta)
  !> it might be more efficient to do all the calculation on proc0 and avoid
  !> the associated communcations.
  logical, parameter :: distribute_work = .true.
contains

  subroutine init_diagnostics_kinetic_energy_transfer
    use mp, only: iproc
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky, akx

    nkynew = 2 * naky - 1
    ikx0 = (ntheta0 + 1) / 2
    iky0 = (nkynew + 1) / 2

    ! allocate regardless of whether this proc does any work because
    ! sum_reduce needs this array to exist on all procs and a selected
    ! sum_reduce function does not seem to exist
    allocate (T_result_theta(-ntgrid:ntgrid))
    allocate (T_result_theta_kxsource(ntheta0, -ntgrid:ntgrid))
    allocate (T_result_theta_kxtarget(ntheta0, -ntgrid:ntgrid))
    allocate (T_result_theta_kysource(nkynew, -ntgrid:ntgrid))

    ! Some allocations only needed on procs with work
    if (is_responsible_for_any_theta_grid_point(iproc, ntgrid)) then

       ! allocate variables needed on each proc involved in calculation
       allocate (kxnew_array(ntheta0)) ; kxnew_array = akx
       allocate (kynew_array(nkynew)) ; kynew_array = 0.
       allocate (valid_mediator_lookup(nkynew, ntheta0, ntheta0)) ; valid_mediator_lookup = .false.
       allocate (z_factor(nkynew, ntheta0, ntheta0)) ; z_factor = 0.

       call fftshift(kxnew_array)
       call extend_ky(kynew_array)
       call get_z_factor(kxnew_array, kynew_array, z_factor)
       call get_mediator_helpers(valid_mediator_lookup)
    end if
  end subroutine init_diagnostics_kinetic_energy_transfer

  subroutine define_kinetic_energy_transfer_dims(file_id)
#ifdef NETCDF
    use neasyf, only: neasyf_dim
    use gs2_io, only: kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
#endif
    integer, intent(in) :: file_id

#ifdef NETCDF
    call neasyf_dim(file_id, trim(kxs_shifted_dim), values=kxnew_array, &
         long_name="Values of kx-source (shifted), used for kinetic energy transfer diagnostic", units="1/rho_r")
    call neasyf_dim(file_id, trim(kxt_shifted_dim), values=kxnew_array, &
         long_name="Values of kx-target (shifted), used for entropy transfer", units="1/rho_r")
    call neasyf_dim(file_id, trim(kys_extended_dim), values=kynew_array, &
         long_name="Values of ky-source extended (shifted), used for kinetic energy transfer diagnostic", units="1/rho_r")
#else
    UNUSED_DUMMY(file_id)
#endif
  end subroutine define_kinetic_energy_transfer_dims

  subroutine finish_diagnostics_kinetic_energy_transfer
    if (allocated(T_result_theta)) deallocate(T_result_theta)
    if (allocated(T_result_theta_kxsource)) deallocate(T_result_theta_kxsource)
    if (allocated(T_result_theta_kxtarget)) deallocate(T_result_theta_kxtarget)
    if (allocated(T_result_theta_kysource)) deallocate(T_result_theta_kysource)
    if (allocated(kynew_array)) deallocate(kynew_array)
    if (allocated(kxnew_array)) deallocate(kxnew_array)
    if (allocated(valid_mediator_lookup)) deallocate(valid_mediator_lookup)
    if (allocated(z_factor)) deallocate(z_factor)
  end subroutine finish_diagnostics_kinetic_energy_transfer

  elemental logical function is_responsible_for_any_theta_grid_point(iproc, ntgrid) result(any_responsible)
    integer, intent(in) :: iproc, ntgrid
    integer :: ntheta
    if (distribute_work) then
       ntheta = 2 * ntgrid + 1
       any_responsible = iproc < ntheta
    else
       any_responsible = iproc == 0
    end if
  end function is_responsible_for_any_theta_grid_point

  elemental logical function is_responsible_for_this_theta_grid_index(iproc, theta_index, nproc, ntgrid) result(responsible)
    integer, intent(in) :: iproc, theta_index, nproc, ntgrid
    integer :: ntheta, theta_index_shifted

    ! If not distributing work then only iproc = 0 is responsible
    if (.not. distribute_work) then
       responsible = iproc == 0
       return
    end if

    ntheta = 2 * ntgrid + 1
    ! theta_index is in [-ntgrid,...,0,...,ntgrid]
    ! theta_index_shifted is in [0, 2*ntgrid]
    ! iproc is in [0, 1, ... , nproc-1]
    theta_index_shifted = theta_index + ntgrid

    ! theta_index_shifted = [0, nproc, 2*nproc, ...] are dealt with by
    ! proc0. [1, nproc+1, 2*nproc+1, ...] by proc1 and so on so have
    ! to take mod when nproc < ntheta
    ! In opposite case the first ntheta procs get one point each.
    if (nproc < ntheta) theta_index_shifted = mod(theta_index_shifted, nproc)
    responsible = iproc == theta_index_shifted
  end function is_responsible_for_this_theta_grid_index

  !> extend ky to negative values and return in "shifted" format.
  !> due to reality condition only ky >= 0 is used in GS2.
  subroutine extend_ky(kynew)
    use kt_grids, only: naky, aky
    real, dimension(:), intent(out) :: kynew ! shape (2*naky-1)
    integer :: i
    kynew(naky:) = aky
    do i = 1, naky-1
       kynew(i) = -aky(naky+1-i)
    end do
  end subroutine extend_ky

  !> calculate:
  !> -(b_hat x k_t . k_s)(b_hat x k_t) . (b_hat x k_s) =
  !> -(k^x_t * k^y_s - k^y_t * k^x_s) * (k^x_t * k^x_s + k^y_t * k^y_s)
  !> where the second term, (b_hat x k_t) . (b_hat x k_s),
  !> is unique for kinetic energy transfer and does not appear for
  !> internal energy transfer or entropy transfer
  subroutine get_z_factor(kx_array, ky_array, result_array)
    use kt_grids, only: ntheta0
    real, dimension(:), intent(in) :: kx_array, ky_array
    real, dimension(:,:,:), intent(out) :: result_array
    integer :: ikys, ikyt, ikxs, ikxt

    ! Pick iky0 as the target ky index here
    ikyt = iky0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ikxt, ikxs, ikys) &
    !$OMP SHARED(ntheta0, nkynew, result_array, kx_array, ky_array, ikyt) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do ikxt = 1, ntheta0
       do ikxs = 1, ntheta0
          do ikys = 1, nkynew
             result_array(ikys, ikxs, ikxt) = &
                  - ( kx_array(ikxt) * ky_array(ikys) - ky_array(ikyt) * kx_array(ikxs) ) &
                  * ( kx_array(ikxt) * kx_array(ikxs) + ky_array(ikyt) * ky_array(ikys) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine get_z_factor

  !> extend a 2D data array to negative ky-values
  !> with the original "unshifted" order
  subroutine get_full(field_arr, result_array)
    use kt_grids, only: ntheta0, naky
    complex, dimension(:,:), intent(in) :: field_arr ! shape (nkx,nky)
    complex, dimension(:,:), intent(out) :: result_array ! shape (2*nky-1,nkx)
    complex, dimension(naky-1, ntheta0) :: temp
    complex, dimension(naky, ntheta0) :: field_arr_transpose

    ! phinew comes in with shape (nkx, nky) but we want
    ! to move to general shape of (kys, kxs, kxt) which means
    ! ky dimensions before the kx dimensions
    field_arr_transpose = transpose(field_arr)

    ! copy in the values for ky >= 0 without change
    result_array(1:naky,:) = field_arr_transpose

    ! reverse order along all kx's and all non-zero ky's.
    ! Also take complex conjugate for reality condition.
    temp = conjg(field_arr_transpose(naky:2:-1,ntheta0:1:-1))

    ! shift the entire negative ky-subarray one index along kx.
    ! After this step the reality condition is satisfied.
    ! This step is required due to the asymmetry which is introduced
    ! by the "unshiftedness" of the array.
    result_array(naky+1:,:) = cshift(temp, shift = -1, dim = 2)
  end subroutine get_full

  !> similar to numpy.fft.fftshift in Python
  !> shift data such that the element for kx=0 or ky=0 is in the center of the array
  subroutine fftshift1D_real(array)
    real, dimension(:), intent(inout) :: array
    integer :: shift_no
    shift_no = -(size(array) - 1) / 2
    array = cshift(array, shift = shift_no, dim = 1)
  end subroutine fftshift1D_real

  !> similar to numpy.fft.fftshift in Python
  !> shift data such that the element for kx=0 and ky=0 is in the center of the matrix
  subroutine fftshift2D_complex(array)
    complex, dimension(:,:), intent(inout) :: array
    integer :: shift_no, dim

    dim = 1
    shift_no = -(size(array, dim = dim) - 1) / 2
    array = cshift(array, shift = shift_no, dim = dim)

    dim = 2
    shift_no = -(size(array, dim = dim) - 1) / 2
    array = cshift(array, shift = shift_no, dim = dim)

  end subroutine fftshift2D_complex

  !> return the mediator index for either kx or ky
  elemental integer function get_ikm(ikt, iks, ik0) result(ikm)
    integer, intent(in) :: ikt, iks, ik0
    ikm = ikt - iks + ik0
  end function get_ikm

  !> return bolean whether a source-mode and target-mode combination yields a valid mediator
  elemental logical function valid_mediator(ikxt, ikxs, ikx0, ikyt, ikys, iky0, nkx, nky)
    integer, intent(in) :: ikxt, ikxs, ikx0, ikyt, ikys, iky0, nkx, nky
    integer :: ikxm, ikym
    ikxm = get_ikm(ikxt, ikxs, ikx0)
    ikym = get_ikm(ikyt, ikys, iky0)
    valid_mediator = (ikxm > 0 .and. ikxm <= nkx) .and. (ikym > 0 .and. ikym <= nky)
  end function valid_mediator

  !> create a boolean lookup for the valid mediator attribute of a
  !> pair of target wavevector and source wavevector
  subroutine get_mediator_helpers(valid_mediator_lookup)
    use kt_grids, only: ntheta0
    logical, dimension(:,:,:), intent(out) :: valid_mediator_lookup ! shape (nky,nkx,nkx)
    integer :: ikys, ikyt, ikxs, ikxt

    ! Pick iky0 as the target ky index here
    ikyt = iky0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ikxt, ikxs, ikys) &
    !$OMP SHARED(ntheta0, nkynew, valid_mediator_lookup, ikx0, ikyt, iky0) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do ikxt = 1, ntheta0
       do ikxs = 1, ntheta0
          do ikys = 1, nkynew
             valid_mediator_lookup(ikys, ikxs, ikxt) = valid_mediator(ikxt, ikxs, ikx0, ikyt, ikys, iky0, ntheta0, nkynew)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine get_mediator_helpers

  !> calculate the mediator mode for all possible pairs of source and target wavevectors
  subroutine get_phi_mediator(phi_source, result_array)
    use kt_grids, only: ntheta0
    complex, dimension(:, :), intent(in) :: phi_source
    complex, dimension(:, :, :), intent(out) :: result_array
    integer :: ikys, ikyt, ikxs, ikxt, ikxm, ikym

    ! Pick iky0 as the target ky index here
    ikyt = iky0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ikxt, ikxs, ikys, ikym, ikxm) &
    !$OMP SHARED(ntheta0, nkynew, result_array, valid_mediator_lookup, &
    !$OMP phi_source, iky0, ikx0, ikyt) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do ikxt = 1, ntheta0
       do ikxs = 1, ntheta0
          do ikys = 1, nkynew
             ikym = get_ikm(ikyt, ikys, iky0)
             ikxm = get_ikm(ikxt, ikxs, ikx0)

             if (valid_mediator_lookup(ikys, ikxs, ikxt)) then
                result_array(ikys, ikxs, ikxt) = phi_source(ikym, ikxm)
             else
                result_array(ikys, ikxs, ikxt) = 0.0
             endif
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine get_phi_mediator

  !> compute the three-wave couplings
  subroutine real_calculate_kinetic_energy_transfer(B_factor, z_factor, &
       phi_target, phi_mediator, phi_source, result_array)
    use kt_grids, only: ntheta0
    real, intent(in) :: B_factor
    real, dimension(:, :, :), intent(in) :: z_factor
    complex, dimension(:), intent(in) :: phi_target
    complex, dimension(:, :, :), intent(in) :: phi_mediator
    complex, dimension(:, :), intent(in) :: phi_source
    real, dimension(:,:,:), intent(out) :: result_array ! shape (nky,nkx,nkx)
    integer :: i, j, k

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i, j, k) &
    !$OMP SHARED(ntheta0, nkynew, result_array, B_factor, z_factor, phi_target, &
    !$OMP phi_mediator, phi_source) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do k = 1, ntheta0
       do j = 1, ntheta0
          do i = 1, nkynew
             result_array(i, j, k) = 2.0 &
                  * B_factor &
                  * z_factor(i, j, k) &
                  * real( &
                  phi_target(k) &
                  * phi_mediator(i, j, k) &
                  * phi_source(i, j) &
                  )
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine real_calculate_kinetic_energy_transfer

  subroutine write_kinetic_energy_transfer(file_id, nout)
#ifdef NETCDF
    use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
    use neasyf, only: neasyf_write
#endif
    integer, intent(in) :: file_id, nout
#ifdef NETCDF
    call define_kinetic_energy_transfer_dims(file_id)
    call neasyf_write(file_id, "kinetic_energy_transfer_theta", T_result_theta, &
         dim_names=[theta_dim, time_dim], start=starts(2, nout), &
         long_name="kinetic energy transfer as a function of theta", &
         units="GS2 units (to be specified)")

    call neasyf_write(file_id, "kinetic_energy_transfer_theta_kxsource", T_result_theta_kxsource, &
         dim_names=[kxs_shifted_dim, theta_dim, time_dim], start=starts(3, nout), &
         long_name="kinetic energy transfer as a function of theta and kx-source", &
         units="GS2 units (to be specified)")

    call neasyf_write(file_id, "kinetic_energy_transfer_theta_kxtarget", T_result_theta_kxtarget, &
         dim_names=[kxt_shifted_dim, theta_dim, time_dim], start=starts(3, nout), &
         long_name="kinetic energy transfer as a function of theta and kx-target", &
         units="GS2 units (to be specified)")

    call neasyf_write(file_id, "kinetic_energy_transfer_theta_kysource", T_result_theta_kysource, &
         dim_names=[kys_extended_dim, theta_dim, time_dim], start=starts(3, nout), &
         long_name="kinetic energy transfer as a function of theta and ky-source", &
         units="GS2 units (to be specified)")
#else
    UNUSED_DUMMY(file_id); UNUSED_DUMMY(nout)
#endif

  end subroutine write_kinetic_energy_transfer

  subroutine calculate_kinetic_energy_transfer_this_theta(theta_index, T_result_3D)
    use fields_arrays, only: phinew
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: bmag
    integer, intent(in) :: theta_index
    real, dimension(:, :, :), intent(out) :: T_result_3D
    real :: B_factor
    complex, dimension(:), allocatable :: phi_target
    complex, dimension(:,:), allocatable :: phi_source
    complex, dimension(:,:,:), allocatable :: phi_mediator

    ! bmag is indexed like theta [-ntgrid,...,0,...,ntgrid]
    B_factor = 1. / bmag(theta_index)**3

    ! phi_source has shape (kys, kxs) because we want to keep the
    ! general order (kys, kxs, kxt) in this file
    allocate (phi_source(2*naky-1, ntheta0)) ; phi_source = 0.
    allocate (phi_target(ntheta0)) ; phi_target = 0.
    allocate (phi_mediator(2*naky-1, ntheta0, ntheta0)) ; phi_mediator = 0.

    call get_full(phinew(theta_index, :, :), phi_source)
    call fftshift(phi_source)
    call get_phi_mediator(phi_source, phi_mediator)
    phi_target = conjg(phi_source(iky0, :))
    call real_calculate_kinetic_energy_transfer(B_factor, z_factor, &
         phi_target, phi_mediator, phi_source, T_result_3D)
  end subroutine calculate_kinetic_energy_transfer_this_theta

  subroutine calculate_kinetic_energy_transfer
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use mp, only: iproc, nproc, sum_reduce
    integer :: theta_index
    real, dimension(:, :, :), allocatable :: T_result_3D

    allocate (T_result_3D(nkynew, ntheta0, ntheta0))

    T_result_theta = 0.
    T_result_theta_kxsource = 0.
    T_result_theta_kxtarget = 0.
    T_result_theta_kysource = 0.

    do theta_index = -ntgrid, ntgrid
       if (.not. is_responsible_for_this_theta_grid_index(&
            iproc, theta_index, nproc, ntgrid)) cycle
       call calculate_kinetic_energy_transfer_this_theta(theta_index, T_result_3D)

       ! Should these sums really use volume average or similar?
       ! Should we just write T_result_3D for each theta? Might get large!

       ! sum over (kxs, kys, kxt). shape (kys, kxs, kxt) -> scalar
       T_result_theta(theta_index) = sum(T_result_3D)
       ! sum over (kxt, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxs)
       T_result_theta_kxsource(:,theta_index) = sum(sum(T_result_3D, dim=1), dim=2)
       ! sum over (kxs, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxt)
       T_result_theta_kxtarget(:,theta_index) = sum(sum(T_result_3D, dim=1), dim=1)
       ! sum over (kxs, kxt) for kyt = 0. shape (kys, kxs, kxt) -> (kys)
       T_result_theta_kysource(:,theta_index) = sum(sum(T_result_3D, dim=2), dim=2)
    end do

    ! only proc0 needs to know the result
    if (distribute_work) then
       call sum_reduce (T_result_theta, 0)
       call sum_reduce (T_result_theta_kxsource, 0)
       call sum_reduce (T_result_theta_kxtarget, 0)
       call sum_reduce (T_result_theta_kysource, 0)
    end if

  end subroutine calculate_kinetic_energy_transfer
end module diagnostics_kinetic_energy_transfer
