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

        complex, dimension (:,:), allocatable ::  entropy_transfer
        allocate (entropy_transfer(ntheta0,naky)) ; entropy_transfer = 0.

        entropy_transfer = 1.0 !fake calculation

        gnostics%current_results%entropy_transfer = entropy_transfer

        deallocate (entropy_transfer)

    end subroutine calculate_entropy_transfer

end module diagnostics_entropy_transfer
