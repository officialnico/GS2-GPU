module diagnostics_entropy_transfer

    implicit none

    private

    public :: init_diagnostics_entropy_transfer
    public :: calculate_entropy_transfer, write_entropy_transfer

contains
    subroutine init_diagnostics_entropy_transfer(gnostics) !complete for test

        use diagnostics_config, only: diagnostics_type

        implicit none

        type(diagnostics_type), intent(inout) :: gnostics

        if(.not. gnostics%write_entropy_transfer) return
        
        !> initialize entropy_transfer to zero
        gnostics%current_results%entropy_transfer = 0.0

    end subroutine init_diagnostics_entropy_transfer

    ! insert all actual modules needed for calculation
    subroutine get_kxky_values(grid_array)
        use kt_grids, only: ntheta0, naky, aky, akx
        use gs2_layouts, only: g_lo, ik_idx, it_idx

        implicit none

        real, dimension(:,:), intent(out) :: grid_array
        integer :: iglo, it, ik
        real :: kx_value, ky_value

        do iglo = g_lo%llim_proc, g_lo%ulim_proc
                it = it_idx(g_lo,iglo)
                ik = ik_idx(g_lo,iglo)
                kx_value = akx(it)
                ky_value = aky(ik)
                grid_array(it,ik) = ky_value
        end do

    end subroutine get_kxky_values


    subroutine write_entropy_transfer(gnostics) !complete for test
        use diagnostics_create_and_write, only: create_and_write_variable
        use diagnostics_dimensions, only: dim_string
        use diagnostics_config, only: diagnostics_type
    
        implicit none
        
        type(diagnostics_type), intent(in out) :: gnostics

        call create_and_write_variable(gnostics, gnostics%rtype, "entropy_transfer", &
         dim_string([gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         "real valued test output, as a function of kx, ky and time", &
         "no units", gnostics%current_results%entropy_transfer)

    end subroutine write_entropy_transfer

    subroutine calculate_entropy_transfer(gnostics)
        use diagnostics_config, only: diagnostics_type
        use kt_grids, only: ntheta0, naky
    
        implicit none
        
        type(diagnostics_type), intent(in out) :: gnostics

        real, dimension (:,:), allocatable ::  entropy_transfer, grid_array

        allocate (entropy_transfer(ntheta0,naky)) ; entropy_transfer = 0.
        allocate (grid_array(ntheta0,naky)) ; grid_array = 0.

        entropy_transfer = 1.0 !fake calculation

        call get_kxky_values(grid_array)

        gnostics%current_results%entropy_transfer = grid_array

        deallocate (entropy_transfer,grid_array)

    end subroutine calculate_entropy_transfer

end module diagnostics_entropy_transfer
