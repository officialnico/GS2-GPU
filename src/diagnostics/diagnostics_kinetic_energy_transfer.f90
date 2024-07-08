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

    public :: init_diagnostics_kinetic_energy_transfer_theta, init_diagnostics_kinetic_energy_transfer_theta_kxsource, init_diagnostics_kinetic_energy_transfer_theta_kxtarget, init_diagnostics_kinetic_energy_transfer_theta_kysource
    public :: write_kinetic_energy_transfer_theta
    public :: write_kinetic_energy_transfer_theta_kxsource
    public :: write_kinetic_energy_transfer_theta_kxtarget
    public :: write_kinetic_energy_transfer_theta_kysource
    public :: calculate_kinetic_energy_transfer

contains
    subroutine init_diagnostics_kinetic_energy_transfer_theta(gnostics)
        use diagnostics_config, only: diagnostics_type

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer_theta) return
        
        gnostics%current_results%kinetic_energy_transfer_theta = 0.0

    end subroutine init_diagnostics_kinetic_energy_transfer_theta

    subroutine init_diagnostics_kinetic_energy_transfer_theta_kxsource(gnostics)
        use diagnostics_config, only: diagnostics_type

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer_theta_kxsource) return
        
        gnostics%current_results%kinetic_energy_transfer_theta_kxsource = 0.0

    end subroutine init_diagnostics_kinetic_energy_transfer_theta_kxsource

    subroutine init_diagnostics_kinetic_energy_transfer_theta_kxtarget(gnostics)
        use diagnostics_config, only: diagnostics_type

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer_theta_kxtarget) return
        
        gnostics%current_results%kinetic_energy_transfer_theta_kxtarget = 0.0

    end subroutine init_diagnostics_kinetic_energy_transfer_theta_kxtarget

    subroutine init_diagnostics_kinetic_energy_transfer_theta_kysource(gnostics)
        use diagnostics_config, only: diagnostics_type

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_kinetic_energy_transfer_theta_kysource) return
        
        gnostics%current_results%kinetic_energy_transfer_theta_kysource = 0.0

    end subroutine init_diagnostics_kinetic_energy_transfer_theta_kysource

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
    subroutine z_hat_factor(kx_array, ky_array, result_array, nkynew, iky0)
        use kt_grids, only: ntheta0

        integer, intent(in) :: nkynew, iky0
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
    subroutine get_full(field_arr, result_array, nkx, nky)
        integer, intent(in) :: nkx, nky
        complex, dimension(:,:), intent(in) :: field_arr ! shape (nky,nkx)
        complex, dimension(:,:), intent(out) :: result_array ! shape (2*nky-1,nkx)
        complex, dimension(nky-1,nkx) :: temp

        ! copy in the values for ky >= 0 without change
        result_array(1:nky,:) = field_arr

        ! reverse order along all kx's and all non-zero ky's.
        ! Also take complex conjugate for reality condition.
        temp = CONJG(field_arr(nky:2:-1,nkx:1:-1))

        ! shift the entire negative ky-subarray one index along kx.
        ! After this step the reality condition is satisfied.
        ! This step is required due to the asymmetry which is introduced
        ! by the "unshiftedness" of the array.
        result_array(nky+1:,:) = cshift(temp, SHIFT=-1,DIM=2)

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
    elemental integer function ikm(ikt, iks, ik0)
        integer, intent(in) :: ikt, iks, ik0

        ikm = ikt - iks + ik0

    end function ikm

    !> return bolean whether a source-mode and target-mode combination yields a valid mediator
    elemental logical function valid_mediator(ikxt, ikxs, ikx0, ikyt, ikys, iky0, nkx, nky)
        integer, intent(in) :: ikxt, ikxs, ikx0, ikyt, ikys, iky0, nkx, nky
        integer :: ikxm, ikym

        ikxm = ikm(ikxt, ikxs, ikx0)
        ikym = ikm(ikyt, ikys, iky0)
        valid_mediator = (ikxm > 0 .and. ikxm <= nkx) .and. (ikym > 0 .and. ikym <= nky)

    end function valid_mediator

    !> create a boolean lookup for the valid mediator attribute of a
    !> pair of target wavevector and source wavevector
    subroutine get_mediator_helpers(ikx0, iky0, nkx, nky, valid_mediator_lookup)
        integer, intent(in) :: ikx0, iky0, nkx, nky
        logical, dimension(:,:,:), intent(out) :: valid_mediator_lookup ! shape (nky,nkx,nkx)
        integer :: ikys, ikyt, ikxs, ikxt, ikxm, ikym

        ikyt = iky0

        do ikxt = 1, nkx
            do ikxs = 1, nkx
                do ikys = 1, nky
                    valid_mediator_lookup(ikys, ikxs, ikxt) = valid_mediator(ikxt, ikxs, ikx0, ikyt, ikys, iky0, nkx, nky)
                end do
            end do
        end do

    end subroutine get_mediator_helpers

    !> calculate the mediator mode for all possible pairs of source and target wavevectors
    subroutine get_mediator(phi, ikx0, iky0, nkx, nky, result_array, valid_mediator_lookup)
        integer, intent(in) :: ikx0, iky0, nkx, nky
        complex, dimension(:,:), intent(in) :: phi
        logical, dimension(:,:,:), intent(in) :: valid_mediator_lookup
        complex, dimension(:,:,:), intent(out) :: result_array
        integer :: ikys, ikyt, ikxs, ikxt, ikxm, ikym
        logical :: valid_mediator

        ikyt = iky0
    
        do ikxt = 1, nkx
            do ikxs = 1, nkx
                do ikys = 1, nky
                    ikym = ikm(ikyt, ikys, iky0)
                    ikxm = ikm(ikxt, ikxs, ikx0)

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
    subroutine compute_kinetic_energy_transfer(phi_source, phi_target, phi_mediator, z_hat_factor, B_factor, nkx, nky, result_array)
        integer, intent(in) :: nkx, nky
        complex, dimension(:,:), intent(in) :: phi_source ! shape (nky,nkx)
        complex, dimension(:), intent(in) :: phi_target ! shape (nkx)
        complex, dimension(:,:,:), intent(in) :: phi_mediator ! shape (nky,nkx,nkx)
        real, dimension(:,:,:), intent(in) :: z_hat_factor ! shape (nky,nkx,nkx)
        real, intent(in) :: B_factor
        real, dimension(:,:,:), intent(out) :: result_array ! shape (nky,nkx,nkx)

        complex, dimension(1, 1, nkx) :: phi_target_new
        complex, dimension(nky, nkx, 1) :: phi_source_new
        integer :: i, j, k

        phi_target_new = reshape(phi_target,(/1,1,nkx/))
        phi_source_new = reshape(phi_source,(/nky,nkx,1/))
    
        do k=1,nkx
            do j=1,nkx
                do i=1,nky
                    result_array(i,j,k) = 2.0 &
                                        * B_factor &
                                        * z_hat_factor(i,j,k) & 
                                        * REAL( &
                                            phi_target_new(1,1,k) &
                                            * phi_mediator(i,j,k) &
                                            * phi_source_new(i,j,1) &
                                        )
                end do
            end do
        end do

    end subroutine compute_kinetic_energy_transfer

    subroutine write_kinetic_energy_transfer_theta(gnostics)
        use diagnostics_config, only: diagnostics_type
        use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
        use neasyf, only: neasyf_write
        
        type(diagnostics_type), intent(in out) :: gnostics

        if (.not. gnostics%writing) return
        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta", gnostics%current_results%kinetic_energy_transfer_theta, &
        dim_names=[theta_dim, time_dim], start=starts(2, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta", &
        units="GS2 units (to be specified)")

    end subroutine write_kinetic_energy_transfer_theta

    subroutine write_kinetic_energy_transfer_theta_kxsource(gnostics)
        use diagnostics_config, only: diagnostics_type
        use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
        use neasyf, only: neasyf_write
        
        type(diagnostics_type), intent(in out) :: gnostics

        if (.not. gnostics%writing) return
        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kxsource", gnostics%current_results%kinetic_energy_transfer_theta_kxsource, &
        dim_names=[kxs_shifted_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and kx-source", &
        units="GS2 units (to be specified)")

    end subroutine write_kinetic_energy_transfer_theta_kxsource

    subroutine write_kinetic_energy_transfer_theta_kxtarget(gnostics)
        use diagnostics_config, only: diagnostics_type
        use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
        use neasyf, only: neasyf_write
        
        type(diagnostics_type), intent(in out) :: gnostics

        if (.not. gnostics%writing) return
        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kxtarget", gnostics%current_results%kinetic_energy_transfer_theta_kxtarget, &
        dim_names=[kxt_shifted_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and kx-target", &
        units="GS2 units (to be specified)")

    end subroutine write_kinetic_energy_transfer_theta_kxtarget

    subroutine write_kinetic_energy_transfer_theta_kysource(gnostics)
        use diagnostics_config, only: diagnostics_type
        use gs2_io, only: starts, time_dim, theta_dim, kxt_shifted_dim, kxs_shifted_dim, kys_extended_dim
        use neasyf, only: neasyf_write
        
        type(diagnostics_type), intent(in out) :: gnostics

        if (.not. gnostics%writing) return
        call neasyf_write(gnostics%file_id, "kinetic_energy_transfer_theta_kysource", gnostics%current_results%kinetic_energy_transfer_theta_kysource, &
        dim_names=[kys_extended_dim, theta_dim, time_dim], start=starts(3, gnostics%nout), &
        long_name="kinetic energy transfer as a function of theta and ky-source", &
        units="GS2 units (to be specified)")

    end subroutine write_kinetic_energy_transfer_theta_kysource

    subroutine calculate_kinetic_energy_transfer(gnostics)

        use diagnostics_config, only: diagnostics_type
        use kt_grids, only: ntheta0, naky, aky, akx
        use fields_arrays, only: phinew
        use theta_grid, only: ntgrid, bmag, theta
        use mp, only: proc0, iproc, nproc, barrier, send, receive, broadcast
        use species, only: nspec, spec
        
        type(diagnostics_type), intent(in out) :: gnostics

        ! declare 0D arrays
        integer :: timestep = 0
        integer :: ikx0, iky0, nkynew
        integer :: theta_index
        integer :: ntheta, n_batches, ibatch, itheta_temp, nproc_counter, nproc_max_temp
        real :: T_result_theta_temp, T_result_theta_temp_temp
        real :: B_factor

        ! prepare arrays with dimension >= 1
        complex, dimension(naky,ntheta0) :: phi_2D
        real, dimension(:,:,:), allocatable :: T_result_3D
        real, dimension(:), allocatable :: T_result_theta
        real, dimension(:,:), allocatable :: T_result_theta_kxsource
        real, dimension(:,:), allocatable :: T_result_theta_kxtarget
        real, dimension(:,:), allocatable :: T_result_theta_kysource
        real, dimension(:), allocatable :: T_result_theta_kxsource_temp
        real, dimension(:), allocatable :: T_result_theta_kxtarget_temp
        real, dimension(:), allocatable :: T_result_theta_kysource_temp
        real, dimension(:), allocatable :: T_result_theta_kxsource_temp_temp
        real, dimension(:), allocatable :: T_result_theta_kxtarget_temp_temp
        real, dimension(:), allocatable :: T_result_theta_kysource_temp_temp
        complex, dimension(:,:), allocatable :: phi_source
        complex, dimension(:), allocatable :: phi_target
        complex, dimension(:,:,:), allocatable :: phi_mediator
        real, dimension(:), allocatable :: kxnew_array, kynew_array 
        real, dimension(:,:,:), allocatable :: z_factor
        logical, dimension(:,:,:), allocatable :: valid_mediator_lookup

        ntheta = 2 * ntgrid + 1 ! itheta is in [-ntgrid,...,0,...,ntgrid]

        ! allocate some variables for the
        ! data collection which happens on proc0
        if (proc0) then
            allocate (T_result_theta(-ntgrid:ntgrid)) ; T_result_theta = 0.
            allocate (T_result_theta_kxsource(ntheta0,-ntgrid:ntgrid)) ; T_result_theta_kxsource = 0.
            allocate (T_result_theta_kxtarget(ntheta0,-ntgrid:ntgrid)) ; T_result_theta_kxtarget = 0.
            allocate (T_result_theta_kysource(2*naky-1,-ntgrid:ntgrid)) ; T_result_theta_kysource = 0.
            allocate (T_result_theta_kxsource_temp_temp(ntheta0)) ; T_result_theta_kxsource_temp_temp = 0.
            allocate (T_result_theta_kxtarget_temp_temp(ntheta0)) ; T_result_theta_kxtarget_temp_temp = 0.
            allocate (T_result_theta_kysource_temp_temp(2*naky-1)) ; T_result_theta_kysource_temp_temp = 0.
        endif

        ! more procs than parallel grid points (expected). 
        ! Some procs deal with one grid point, some with none.
        if (nproc .ge. ntheta) then

            if (proc0 .AND. (timestep .eq. 0) .AND. (nproc > ntheta)) write(*,*) "nproc > ntheta detected"
            if (proc0 .AND. (timestep .eq. 0) .AND. (nproc .eq. ntheta)) write(*,*) "nproc = ntheta detected"

            ! all other procs can relax. iproc is in [0, 1, ... , nproc-1]
            if (iproc < ntheta) then 
                theta_index = iproc - ntgrid ! obtains min of -ntgrid, and max of +ntgrid
                ! ***************************
                ! **** calculation start ****
                ! ***************************
                ! allocate variables needed on each proc involved in calculation
                allocate (kxnew_array(ntheta0)) ; kxnew_array = 0.
                allocate (kynew_array(2*naky-1)) ; kynew_array = 0.
                allocate (valid_mediator_lookup(2*naky-1,ntheta0,ntheta0)) ; valid_mediator_lookup = .false.
                allocate (z_factor(2*naky-1,ntheta0,ntheta0)) ; z_factor = 0.

                ! Do common calculations.
                nkynew = 2*naky-1
                ikx0 = (ntheta0+1)/2
                iky0 = (nkynew+1)/2
                kxnew_array = akx
                call fftshift1D_real(kxnew_array)
                call extend_ky(kynew_array)
                call z_hat_factor(kxnew_array, kynew_array, z_factor, nkynew, iky0)
                call get_mediator_helpers(ikx0, iky0, ntheta0, nkynew, valid_mediator_lookup)
                B_factor = 1. / bmag(theta_index)**3 ! bmag is indexed like theta [-ntgrid,...,0,...,ntgrid]

                allocate (T_result_3D(2*naky-1,ntheta0,ntheta0)) ; T_result_3D = 0.
                allocate (T_result_theta_kxsource_temp(ntheta0)) ; T_result_theta_kxsource_temp = 0.
                allocate (T_result_theta_kxtarget_temp(ntheta0)) ; T_result_theta_kxtarget_temp = 0.
                allocate (T_result_theta_kysource_temp(2*naky-1)) ; T_result_theta_kysource_temp = 0.
                allocate (phi_source(2*naky-1,ntheta0)) ; phi_source = 0.
                allocate (phi_target(ntheta0)) ; phi_target = 0.
                allocate (phi_mediator(2*naky-1,ntheta0,ntheta0)) ; phi_mediator = 0.

                phi_2D(:,:) = TRANSPOSE(phinew(theta_index,:,:)) ! shape (kx, ky) -> (ky, kx)
                call get_full(phi_2D, phi_source, ntheta0,  naky)
                call fftshift2D_complex(phi_source)
                call get_mediator(phi_source, ikx0, iky0, ntheta0, nkynew, phi_mediator, valid_mediator_lookup)
                phi_target = CONJG(phi_source(iky0,:))
                call compute_kinetic_energy_transfer(phi_source, phi_target, phi_mediator, z_factor, B_factor, ntheta0, nkynew, T_result_3D)

                T_result_theta_temp = SUM(SUM(SUM(T_result_3D, DIM=1), DIM=1), DIM=1) ! sum over (kxs, kys, kxt). shape (kys, kxs, kxt) -> ()
                T_result_theta_kxsource_temp = SUM(SUM(T_result_3D, DIM=1), DIM=2) ! sum over (kxt, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxs)
                T_result_theta_kxtarget_temp = SUM(SUM(T_result_3D, DIM=1), DIM=1) ! sum over (kxs, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxt)
                T_result_theta_kysource_temp = SUM(SUM(T_result_3D, DIM=2), DIM=2) ! sum over (kxs, kxt) for kyt = 0. shape (kys, kxs, kxt) -> (kys)
                ! ***************************
                ! ***** calculation end *****
                ! ***************************

                ! ***************************************
                ! *** Collect results on proc0 start ****
                ! ***************************************
                if (.not. proc0) then
                    ! send data to proc0
                    call send(theta_index, 0) 
                    call send(T_result_theta_temp, 0)
                    call send(T_result_theta_kxsource_temp, 0)
                    call send(T_result_theta_kxtarget_temp, 0)
                    call send(T_result_theta_kysource_temp, 0)
                else
                    ! save values from proc0
                    T_result_theta(theta_index) = T_result_theta_temp
                    T_result_theta_kxsource(:,theta_index) = T_result_theta_kxsource_temp
                    T_result_theta_kxtarget(:,theta_index) = T_result_theta_kxtarget_temp
                    T_result_theta_kysource(:,theta_index) = T_result_theta_kysource_temp

                    ! collect values from all other procs
                    ! procs start counting from 0, skip proc0 because that is the one we collect on
                    do nproc_counter = 1, ntheta-1 
                        call receive(itheta_temp, nproc_counter)
                        call receive(T_result_theta_temp_temp, nproc_counter)
                        call receive(T_result_theta_kxsource_temp_temp, nproc_counter)
                        call receive(T_result_theta_kxtarget_temp_temp, nproc_counter)
                        call receive(T_result_theta_kysource_temp_temp, nproc_counter)
                        T_result_theta(itheta_temp) = T_result_theta_temp_temp
                        T_result_theta_kxsource(:,itheta_temp) = T_result_theta_kxsource_temp_temp
                        T_result_theta_kxtarget(:,itheta_temp) = T_result_theta_kxtarget_temp_temp
                        T_result_theta_kysource(:,itheta_temp) = T_result_theta_kysource_temp_temp
                    end do
                endif
                ! *************************************
                ! *** Collect results on proc0 end ****
                ! *************************************
            endif

        ! more parallel grid points than procs. 
        ! Some or all procs deal with more than one grid point. 
        ! Some procs will deal with one more point than others if ntheta % nproc != 0.
        else
            ! allocate here once before cycling through batches
            ! allocate variables needed on each proc involved in calculation
            allocate (kxnew_array(ntheta0)) ; kxnew_array = 0.
            allocate (kynew_array(2*naky-1)) ; kynew_array = 0.
            allocate (valid_mediator_lookup(2*naky-1,ntheta0,ntheta0)) ; valid_mediator_lookup = .false.
            allocate (z_factor(2*naky-1,ntheta0,ntheta0)) ; z_factor = 0.

            ! Do common calculations.
            nkynew = 2*naky-1
            ikx0 = (ntheta0+1)/2
            iky0 = (nkynew+1)/2
            kxnew_array = akx
            call fftshift1D_real(kxnew_array)
            call extend_ky(kynew_array)
            call z_hat_factor(kxnew_array, kynew_array, z_factor, nkynew, iky0)
            call get_mediator_helpers(ikx0, iky0, ntheta0, nkynew, valid_mediator_lookup)

            allocate (T_result_3D(2*naky-1,ntheta0,ntheta0)) ; T_result_3D = 0.
            allocate (T_result_theta_kxsource_temp(ntheta0)) ; T_result_theta_kxsource_temp = 0.
            allocate (T_result_theta_kxtarget_temp(ntheta0)) ; T_result_theta_kxtarget_temp = 0.
            allocate (T_result_theta_kysource_temp(2*naky-1)) ; T_result_theta_kysource_temp = 0.
            allocate (phi_source(2*naky-1,ntheta0)) ; phi_source = 0.
            allocate (phi_target(ntheta0)) ; phi_target = 0.
            allocate (phi_mediator(2*naky-1,ntheta0,ntheta0)) ; phi_mediator = 0.

            n_batches = CEILING(REAL(ntheta)/REAL(nproc)) ! Round up to the nearest integer

            if (proc0 .AND. (timestep .eq. 0) .AND. (MOD(ntheta, nproc) .eq. 0.)) write(*,*) "nproc < ntheta with ntheta % nproc = 0 detected"
            if (proc0 .AND. (timestep .eq. 0) .AND. (MOD(ntheta, nproc) .ne. 0.)) write(*,*) "nproc < ntheta with ntheta % nproc != 0 detected"

            do ibatch = 0, n_batches-1
                theta_index = (iproc - ntgrid) + ibatch * nproc
                
                ! this could only be false for last batch 
                ! where not all procs are needed anymore
                if (theta_index .le. ntgrid) then
                    ! ***************************
                    ! **** calculation start ****
                    ! ***************************
                    B_factor = 1. / bmag(theta_index)**3 ! bmag is indexed like theta [-ntgrid,...,0,...,ntgrid]
                    phi_2D(:,:) = TRANSPOSE(phinew(theta_index,:,:)) ! shape (kx, ky) -> (ky, kx)                    
                    call get_full(phi_2D, phi_source, ntheta0,  naky)
                    call fftshift2D_complex(phi_source)
                    call get_mediator(phi_source, ikx0, iky0, ntheta0, nkynew, phi_mediator, valid_mediator_lookup)
                    phi_target = CONJG(phi_source(iky0,:))
                    call compute_kinetic_energy_transfer(phi_source, phi_target, phi_mediator, z_factor, B_factor, ntheta0, nkynew, T_result_3D)

                    T_result_theta_temp = SUM(SUM(SUM(T_result_3D, DIM=1), DIM=1), DIM=1) ! sum over (kxs, kys, kxt). shape (kys, kxs, kxt) -> ()
                    T_result_theta_kxsource_temp = SUM(SUM(T_result_3D, DIM=1), DIM=2) ! sum over (kxt, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxs)
                    T_result_theta_kxtarget_temp = SUM(SUM(T_result_3D, DIM=1), DIM=1) ! sum over (kxs, kys) for kyt = 0. shape (kys, kxs, kxt) -> (kxt)
                    T_result_theta_kysource_temp = SUM(SUM(T_result_3D, DIM=2), DIM=2) ! sum over (kxs, kxt) for kyt = 0. shape (kys, kxs, kxt) -> (kys)
                    ! ***************************
                    ! ***** calculation end *****
                    ! ***************************

                    ! ***************************************
                    ! *** Collect results on proc0 start ****
                    ! ***************************************
                    if (.not. proc0) then
                        ! send data to proc0
                        call send(theta_index, 0) 
                        call send(T_result_theta_temp, 0)
                        call send(T_result_theta_kxsource_temp, 0)
                        call send(T_result_theta_kxtarget_temp, 0)
                        call send(T_result_theta_kysource_temp, 0)
                    else
                        ! save values from proc0
                        T_result_theta(theta_index) = T_result_theta_temp
                        T_result_theta_kxsource(:,theta_index) = T_result_theta_kxsource_temp
                        T_result_theta_kxtarget(:,theta_index) = T_result_theta_kxtarget_temp
                        T_result_theta_kysource(:,theta_index) = T_result_theta_kysource_temp

                        ! collect values from all other procs
                        ! catch case where ntheta a multiple of nproc
                        ! in which case nproc_max_temp = nproc for all batches
                        ! OR if it is not the last batch then also nproc_max_temp = nproc
                        if ((MOD(ntheta, nproc) .eq. 0.) .OR. (ibatch < n_batches-1)) then
                            nproc_max_temp = nproc
                        else
                            nproc_max_temp = MOD(ntheta, nproc)
                        endif
                        ! procs start counting from 0, skip proc0 because that is the one we collect on
                        do nproc_counter = 1, nproc_max_temp-1 
                            call receive(itheta_temp, nproc_counter)
                            call receive(T_result_theta_temp_temp, nproc_counter)
                            call receive(T_result_theta_kxsource_temp_temp, nproc_counter)
                            call receive(T_result_theta_kxtarget_temp_temp, nproc_counter)
                            call receive(T_result_theta_kysource_temp_temp, nproc_counter)
                            T_result_theta(itheta_temp) = T_result_theta_temp_temp
                            T_result_theta_kxsource(:,itheta_temp) = T_result_theta_kxsource_temp_temp
                            T_result_theta_kxtarget(:,itheta_temp) = T_result_theta_kxtarget_temp_temp
                            T_result_theta_kysource(:,itheta_temp) = T_result_theta_kysource_temp_temp
                        end do
                    endif
                    ! *************************************
                    ! *** Collect results on proc0 end ****
                    ! *************************************
                end if
            end do
        end if

        ! write results to diagnostic attribute on proc0 and deallocate proc0-specific variables
        if (proc0) then
            gnostics%current_results%kinetic_energy_transfer_theta = T_result_theta
            gnostics%current_results%kinetic_energy_transfer_theta_kxsource  = T_result_theta_kxsource
            gnostics%current_results%kinetic_energy_transfer_theta_kxtarget  = T_result_theta_kxtarget 
            gnostics%current_results%kinetic_energy_transfer_theta_kysource = T_result_theta_kysource
            deallocate(T_result_theta)
            deallocate(T_result_theta_kxsource)
            deallocate(T_result_theta_kxtarget)
            deallocate(T_result_theta_kysource)
            deallocate(T_result_theta_kxsource_temp_temp)
            deallocate(T_result_theta_kxtarget_temp_temp)
            deallocate(T_result_theta_kysource_temp_temp)
        endif

        ! deallocate variabels used on all procs which did work
        if (iproc < ntheta) then
            deallocate(z_factor)
            deallocate(valid_mediator_lookup)
            deallocate(kxnew_array)
            deallocate(kynew_array)
            deallocate(T_result_3D)
            deallocate(T_result_theta_kxsource_temp)
            deallocate(T_result_theta_kxtarget_temp)
            deallocate(T_result_theta_kysource_temp)
            deallocate(phi_source)
            deallocate(phi_target)
            deallocate(phi_mediator)
        endif

        timestep = timestep + 1 ! advance timestep counter on all procs

    end subroutine calculate_kinetic_energy_transfer

end module diagnostics_kinetic_energy_transfer