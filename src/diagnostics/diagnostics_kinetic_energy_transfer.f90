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
    real, dimension(:,:,:), allocatable :: T_result_3D
    complex, dimension(:,:), allocatable :: phi_source
    complex, dimension(:), allocatable :: phi_target
    complex, dimension(:,:,:), allocatable :: phi_mediator
    integer :: ikx0, iky0, nkynew
    logical :: any_responsibility

    real, dimension(:), allocatable :: T_result_theta_kxsource_temp
    real, dimension(:), allocatable :: T_result_theta_kxtarget_temp
    real, dimension(:), allocatable :: T_result_theta_kysource_temp

    real, dimension(:), allocatable :: T_result_theta
    real, dimension(:,:), allocatable :: T_result_theta_kxsource
    real, dimension(:,:), allocatable :: T_result_theta_kxtarget
    real, dimension(:,:), allocatable :: T_result_theta_kysource

contains
    subroutine init_diagnostics_kinetic_energy_transfer(gnostics)
        use diagnostics_config, only: diagnostics_type
        use mp, only: iproc, nproc
        use theta_grid, only: ntgrid
        use kt_grids, only: ntheta0, naky, aky, akx

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer) return
        
        gnostics%current_results%kinetic_energy_transfer_theta = 0.0
        gnostics%current_results%kinetic_energy_transfer_theta_kxsource = 0.0
        gnostics%current_results%kinetic_energy_transfer_theta_kxtarget = 0.0
        gnostics%current_results%kinetic_energy_transfer_theta_kysource = 0.0

        ! check if this proc every does any work for this diagnostic
        any_responsibility = is_responsible_for_any_theta_grid_point(iproc, nproc, ntgrid)

        ! allocate regardless of whether this proc does any work because sum_reduce needs
        ! this array to exist on all procs and a selected sum_reduce function does not seem to exist
        allocate (T_result_theta(-ntgrid:ntgrid)) ; T_result_theta = 0.
        allocate (T_result_theta_kxsource(ntheta0,-ntgrid:ntgrid)) ; T_result_theta_kxsource = 0.
        allocate (T_result_theta_kxtarget(ntheta0,-ntgrid:ntgrid)) ; T_result_theta_kxtarget = 0.
        allocate (T_result_theta_kysource(2*naky-1,-ntgrid:ntgrid)) ; T_result_theta_kysource = 0.

        if (any_responsibility) then
            ! allocate variables needed on each proc involved in calculation
            allocate (kxnew_array(ntheta0)) ; kxnew_array = 0.
            allocate (kynew_array(2*naky-1)) ; kynew_array = 0.
            allocate (valid_mediator_lookup(2*naky-1,ntheta0,ntheta0)) ; valid_mediator_lookup = .false.
            allocate (z_factor(2*naky-1,ntheta0,ntheta0)) ; z_factor = 0.

            allocate (T_result_3D(2*naky-1,ntheta0,ntheta0)) ; T_result_3D = 0.
            allocate (phi_source(2*naky-1,ntheta0)) ; phi_source = 0.
            allocate (phi_target(ntheta0)) ; phi_target = 0.
            allocate (phi_mediator(2*naky-1,ntheta0,ntheta0)) ; phi_mediator = 0.

            allocate (T_result_theta_kxsource_temp(ntheta0)) ; T_result_theta_kxsource_temp = 0.
            allocate (T_result_theta_kxtarget_temp(ntheta0)) ; T_result_theta_kxtarget_temp = 0.
            allocate (T_result_theta_kysource_temp(2*naky-1)) ; T_result_theta_kysource_temp = 0.

            nkynew = 2*naky-1
            ikx0 = (ntheta0+1)/2
            iky0 = (nkynew+1)/2
            kxnew_array = akx
            call fftshift1D_real(kxnew_array)
            call extend_ky(kynew_array)
            call z_hat_factor(kxnew_array, kynew_array, z_factor)
            call get_mediator_helpers(valid_mediator_lookup)
        endif

    end subroutine init_diagnostics_kinetic_energy_transfer

    subroutine finish_diagnostics_kinetic_energy_transfer(gnostics)
        use diagnostics_config, only: diagnostics_type
        use mp, only: iproc, nproc
        use theta_grid, only: ntgrid

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer) return

        deallocate(T_result_theta)
        deallocate(T_result_theta_kxsource)
        deallocate(T_result_theta_kxtarget)
        deallocate(T_result_theta_kysource)

        if (any_responsibility) then
            deallocate(kynew_array)
            deallocate(kxnew_array)
            deallocate(valid_mediator_lookup)
            deallocate(z_factor)
            deallocate(T_result_3D)
            deallocate(phi_source)
            deallocate(phi_target)
            deallocate(phi_mediator)
            deallocate(T_result_theta_kxsource_temp)
            deallocate(T_result_theta_kxtarget_temp)
            deallocate(T_result_theta_kysource_temp)
        endif

    end subroutine finish_diagnostics_kinetic_energy_transfer

    elemental logical function is_responsible_for_any_theta_grid_point(iproc, nproc, ntgrid) result(any_responsible)
        integer, intent(in) :: iproc, nproc, ntgrid
        integer :: ntheta

        ntheta = 2 * ntgrid + 1 ! itheta is in [-ntgrid,...,0,...,ntgrid]
        any_responsible = .true.

        if (iproc >= ntheta) then ! these procs can relax. iproc is in [0, 1, ... , nproc-1]
            any_responsible = .false.
        end if

    end function is_responsible_for_any_theta_grid_point

    elemental logical function is_responsible_for_this_theta_grid_index(iproc, theta_index, nproc, ntgrid) result(responsible)
        integer, intent(in) :: iproc, theta_index, nproc, ntgrid
        integer :: ntheta, theta_index_shifted, theta_index_modulo

        ntheta = 2 * ntgrid + 1 
        responsible = .false.

        ! theta_index is in [-ntgrid,...,0,...,ntgrid]
        ! iproc is in [0, 1, ... , nproc-1]
        if (nproc >= ntheta) then ! each proc deals with one or no theta points
            theta_index_shifted = theta_index + ntgrid ! obtains min of 0 and max of 2*ntgrid
            if (iproc .eq. theta_index_shifted) then
                responsible = .true.
            endif
        else ! each proc deals with at least one theta index
            theta_index_shifted = theta_index + ntgrid ! obtains min of 0 and max of 2*ntgrid
            theta_index_modulo = MOD(theta_index_shifted, nproc) ! theta_index_shifted = [0, nproc, 2*nproc, ...] are dealt with by proc0. [1, nproc+1, 2*nproc+1, ...] by proc1 and so on.
            if (iproc .eq. theta_index_modulo) then
                responsible = .true.
            endif
        endif

    end function is_responsible_for_this_theta_grid_index

    !> extend ky to negative values.
    !> due to reality condition only ky >= 0 is used in GS2.
    subroutine extend_ky(kynew)
        use kt_grids, only: naky, aky

        real, dimension(:), intent(out) :: kynew ! shape (2*naky-1)
        integer :: i

        kynew(naky:) = aky
    
        do i=1, naky-1
            kynew(i) = -aky(naky+1-i)
        end do

    end subroutine extend_ky

    !> calculate:
    !> -(b_hat x k_t . k_s)(b_hat x k_t) . (b_hat x k_s) = -(k^x_t * k^y_s - k^y_t * k^x_s) * (k^x_t * k^x_s + k^y_t * k^y_s)
    !> where the second term
    !> (b_hat x k_t) . (b_hat x k_s)
    !> is unique for kinetic energy transfer and does not appear for internal energy transfer or entropy transfer
    subroutine z_hat_factor(kx_array, ky_array, result_array)
        use kt_grids, only: ntheta0

        ! integer, intent(in) :: nkynew, iky0
        real, dimension(:), intent(in) :: kx_array, ky_array
        real, dimension(:,:,:), intent(out) :: result_array
        integer :: ikys, ikyt, ikxs, ikxt

        ikyt = iky0
    
        do ikxt=1, ntheta0
            do ikxs=1, ntheta0
                do ikys=1, nkynew
                    result_array(ikys,ikxs,ikxt) = -1.0 * ( kx_array(ikxt)*ky_array(ikys) - ky_array(ikyt)*kx_array(ikxs) ) * ( kx_array(ikxt)*kx_array(ikxs) + ky_array(ikyt)*ky_array(ikys) )
                end do
            end do
        end do

    end subroutine z_hat_factor

    !> extend a 2D data array to negative ky-values
    !> with the original "unshifted" order
    subroutine get_full(field_arr, result_array)
        use kt_grids, only: ntheta0, naky

        complex, dimension(:,:), intent(in) :: field_arr ! shape (nkx,nky)
        complex, dimension(:,:), intent(out) :: result_array ! shape (2*nky-1,nkx)
        complex, dimension(naky-1,ntheta0) :: temp
        complex, dimension(naky,ntheta0) :: field_arr_transpose

        ! phinew comes in with shape (nkx, nky) but we want
        ! to move to general shape of (kys, kxs, kxt) which means 
        ! ky dimensions before the kx dimensions
        field_arr_transpose = TRANSPOSE(field_arr)

        ! copy in the values for ky >= 0 without change
        result_array(1:naky,:) = field_arr_transpose

        ! reverse order along all kx's and all non-zero ky's.
        ! Also take complex conjugate for reality condition.
        temp = CONJG(field_arr_transpose(naky:2:-1,ntheta0:1:-1))

        ! shift the entire negative ky-subarray one index along kx.
        ! After this step the reality condition is satisfied.
        ! This step is required due to the asymmetry which is introduced
        ! by the "unshiftedness" of the array.
        result_array(naky+1:,:) = cshift(temp, SHIFT=-1,DIM=2)

    end subroutine get_full

    !> similar to numpy.fft.fftshift in Python
    !> shift data such that the element for kx=0 or ky=0 is in the center of the array
    subroutine fftshift1D_real(array)
        real, dimension(:), intent(inout) :: array
        integer :: shift_no

        shift_no = -(size(array)-1)/2
        array = cshift(array, SHIFT=shift_no, DIM=1)

    end subroutine fftshift1D_real

    !> similar to numpy.fft.fftshift in Python
    !> shift data such that the element for kx=0 and ky=0 is in the center of the matrix
    subroutine fftshift2D_complex(array)
        complex, dimension(:,:), intent(inout) :: array
        integer :: shift_no, dim

        dim = 1
        shift_no = -(size(array, DIM=dim)-1)/2
        array = cshift(array, SHIFT=shift_no, DIM=dim)

        dim = 2
        shift_no = -(size(array, DIM=dim)-1)/2
        array = cshift(array, SHIFT=shift_no, DIM=dim)

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
        integer :: ikys, ikyt, ikxs, ikxt, ikxm, ikym

        ikyt = iky0

        do ikxt = 1, ntheta0
            do ikxs = 1, ntheta0
                do ikys = 1, nkynew
                    valid_mediator_lookup(ikys, ikxs, ikxt) = valid_mediator(ikxt, ikxs, ikx0, ikyt, ikys, iky0, ntheta0, nkynew)
                end do
            end do
        end do

    end subroutine get_mediator_helpers

    !> calculate the mediator mode for all possible pairs of source and target wavevectors
    subroutine get_mediator(phi, result_array, valid_mediator_lookup)
        use kt_grids, only: ntheta0

        complex, dimension(:,:), intent(in) :: phi
        logical, dimension(:,:,:), intent(in) :: valid_mediator_lookup
        complex, dimension(:,:,:), intent(out) :: result_array
        integer :: ikys, ikyt, ikxs, ikxt, ikxm, ikym
        logical :: valid_mediator

        ikyt = iky0
    
        do ikxt = 1, ntheta0
            do ikxs = 1, ntheta0
                do ikys = 1, nkynew
                    ikym = get_ikm(ikyt, ikys, iky0)
                    ikxm = get_ikm(ikxt, ikxs, ikx0)

                    valid_mediator = valid_mediator_lookup(ikys, ikxs, ikxt)

                    if (valid_mediator) then
                        result_array(ikys, ikxs, ikxt) = phi(ikym,ikxm)
                    else
                        result_array(ikys, ikxs, ikxt) = 0.0
                    endif
                end do
            end do
        end do

    end subroutine get_mediator

    !> compute the three-wave couplings
    subroutine compute_kinetic_energy_transfer(phi_source, phi_target, phi_mediator, z_hat_factor, B_factor, result_array)
        use kt_grids, only: ntheta0

        complex, dimension(:,:), intent(in) :: phi_source ! shape (nky,nkx)
        complex, dimension(:), intent(in) :: phi_target ! shape (nkx)
        complex, dimension(:,:,:), intent(in) :: phi_mediator ! shape (nky,nkx,nkx)
        real, dimension(:,:,:), intent(in) :: z_hat_factor ! shape (nky,nkx,nkx)
        real, intent(in) :: B_factor
        real, dimension(:,:,:), intent(out) :: result_array ! shape (nky,nkx,nkx)
        integer :: i, j, k
    
        do k=1,ntheta0
            do j=1,ntheta0
                do i=1,nkynew
                    result_array(i,j,k) = 2.0 &
                                        * B_factor &
                                        * z_hat_factor(i,j,k) & 
                                        * REAL( &
                                            phi_target(k) &
                                            * phi_mediator(i,j,k) &
                                            * phi_source(i,j) &
                                        )
                end do
            end do
        end do

    end subroutine compute_kinetic_energy_transfer

    subroutine write_kinetic_energy_transfer(gnostics)
        use diagnostics_config, only: diagnostics_type
        use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
        use neasyf, only: neasyf_write
        
        type(diagnostics_type), intent(in out) :: gnostics

        if (.not. gnostics%writing) return
        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta", gnostics%current_results%kinetic_energy_transfer_theta, &
        dim_names=[theta_dim, time_dim], start=starts(2, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta", &
        units="GS2 units (to be specified)")

        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kxsource", gnostics%current_results%kinetic_energy_transfer_theta_kxsource, &
        dim_names=[kxs_shifted_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and kx-source", &
        units="GS2 units (to be specified)")

        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kxtarget", gnostics%current_results%kinetic_energy_transfer_theta_kxtarget, &
        dim_names=[kxt_shifted_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and kx-target", &
        units="GS2 units (to be specified)")

        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kysource", gnostics%current_results%kinetic_energy_transfer_theta_kysource, &
        dim_names=[kys_extended_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and ky-source", &
        units="GS2 units (to be specified)")

    end subroutine write_kinetic_energy_transfer

    subroutine calculate_kinetic_energy_transfer_this_theta(theta_index, T_result_theta_temp, T_result_theta_kxsource_temp, T_result_theta_kxtarget_temp, T_result_theta_kysource_temp)
        use fields_arrays, only: phinew
        use kt_grids, only: ntheta0, naky
        use theta_grid, only: ntgrid, bmag

        integer, intent(in) :: theta_index
        real, intent(out) :: T_result_theta_temp
        real, dimension(:), intent(out) :: T_result_theta_kxsource_temp, T_result_theta_kxtarget_temp, T_result_theta_kysource_temp

        complex, dimension(ntheta0,naky) :: phi_temp
        real :: B_factor
        
        B_factor = 1. / bmag(theta_index)**3 ! bmag is indexed like theta [-ntgrid,...,0,...,ntgrid]
        phi_temp = phinew(theta_index,:,:) ! shape (kx, ky)
        call get_full(phi_temp, phi_source) ! phi_source has shape (kys, kxs) because we want to keep the general order (kys, kxs, kxt) in this file
        call fftshift2D_complex(phi_source)
        call get_mediator(phi_source, phi_mediator, valid_mediator_lookup)
        phi_target = CONJG(phi_source(iky0,:))
        call compute_kinetic_energy_transfer(phi_source, phi_target, phi_mediator, z_factor, B_factor, T_result_3D)

        T_result_theta_temp = SUM(SUM(SUM(T_result_3D, DIM=1), DIM=1), DIM=1) ! sum over (kxs, kys, kxt). shape (kys, kxs, kxt) -> ()
        T_result_theta_kxsource_temp = SUM(SUM(T_result_3D, DIM=1), DIM=2) ! sum over (kxt, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxs)
        T_result_theta_kxtarget_temp = SUM(SUM(T_result_3D, DIM=1), DIM=1) ! sum over (kxs, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxt)
        T_result_theta_kysource_temp = SUM(SUM(T_result_3D, DIM=2), DIM=2) ! sum over (kxs, kxt) for kyt = 0. shape (kys, kxs, kxt) -> (kys)

    end subroutine calculate_kinetic_energy_transfer_this_theta

    subroutine calculate_kinetic_energy_transfer(gnostics)

        use diagnostics_config, only: diagnostics_type
        use theta_grid, only: ntgrid
        use mp, only: proc0, iproc, nproc, sum_reduce
        
        type(diagnostics_type), intent(in out) :: gnostics

        integer :: timestep = 0
        integer :: theta_index
        real :: T_result_theta_temp
        logical :: responsible

        T_result_theta = 0.
        T_result_theta_kxsource = 0.
        T_result_theta_kxtarget = 0.
        T_result_theta_kysource = 0.

        if (any_responsibility) then
            do theta_index = -ntgrid, ntgrid
                responsible = is_responsible_for_this_theta_grid_index(iproc, theta_index, nproc, ntgrid)
                if (responsible) then
                    call calculate_kinetic_energy_transfer_this_theta(theta_index, T_result_theta_temp, T_result_theta_kxsource_temp, T_result_theta_kxtarget_temp, T_result_theta_kysource_temp)
                    T_result_theta(theta_index) = T_result_theta_temp
                    T_result_theta_kxsource(:,theta_index) = T_result_theta_kxsource_temp
                    T_result_theta_kxtarget(:,theta_index) = T_result_theta_kxtarget_temp
                    T_result_theta_kysource(:,theta_index) = T_result_theta_kysource_temp
                endif
            end do
        endif

        ! only proc0 needs to know the result
        call sum_reduce (T_result_theta, 0)
        call sum_reduce (T_result_theta_kxsource, 0)
        call sum_reduce (T_result_theta_kxtarget, 0)
        call sum_reduce (T_result_theta_kysource, 0)

        if (proc0) then
            gnostics%current_results%kinetic_energy_transfer_theta = T_result_theta
            gnostics%current_results%kinetic_energy_transfer_theta_kxsource  = T_result_theta_kxsource
            gnostics%current_results%kinetic_energy_transfer_theta_kxtarget  = T_result_theta_kxtarget 
            gnostics%current_results%kinetic_energy_transfer_theta_kysource = T_result_theta_kysource
        endif

    end subroutine calculate_kinetic_energy_transfer

end module diagnostics_kinetic_energy_transfer