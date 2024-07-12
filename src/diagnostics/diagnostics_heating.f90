!> Subroutines for writing out the turbulent heating...
!! the turbulent energy dissipated in collisions and hyperviscosity
module diagnostics_heating
  use gs2_heating, only: heating_diagnostics
  implicit none
  private

  public :: init_diagnostics_heating, finish_diagnostics_heating
  public :: calculate_heating, write_heating

  !> Set averages in gnostics\%current_results to 0.
  public :: reset_averages_and_counters

  type (heating_diagnostics) :: h
  type (heating_diagnostics), dimension(:), allocatable :: h_hist
  type (heating_diagnostics), dimension(:,:), allocatable :: hk
  type (heating_diagnostics), dimension(:,:,:), allocatable :: hk_hist
contains

  !> FIXME : Add documentation  
  subroutine finish_diagnostics_heating(gnostics)
    use gs2_heating, only: del_htype
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (in) :: gnostics

    if (gnostics%write_heating) then
       call del_htype (h)
       call del_htype (h_hist)
       call del_htype (hk_hist)
       call del_htype (hk)
    endif
    
    if (allocated(h_hist)) deallocate (h_hist, hk_hist, hk)
  end subroutine finish_diagnostics_heating

  !> FIXME : Add documentation  
  subroutine init_diagnostics_heating(gnostics)
    use gs2_heating, only: init_htype
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (in) :: gnostics
    integer :: navg
    navg = gnostics%navg
    
    ! allocate heating diagnostic data structures
    if (gnostics%write_heating) then
       allocate (h_hist(0:navg-1))
       call init_htype (h_hist,  nspec)
       
       allocate (hk_hist(ntheta0,naky,0:navg-1))
       call init_htype (hk_hist, nspec)
       
       call init_htype (h,  nspec)

       allocate (hk(ntheta0, naky))
       call init_htype (hk, nspec)
    else
       allocate (h_hist(0))
       allocate (hk(1,1))
       allocate (hk_hist(1,1,0))
    end if
  end subroutine init_diagnostics_heating

  !> FIXME : Add documentation
  subroutine reset_averages_and_counters(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (inout) :: gnostics
    gnostics%current_results%species_heating_avg = 0.0
  end subroutine reset_averages_and_counters

  !> FIXME : Add documentation  
  subroutine calculate_heating (gnostics)
    use diagnostics_config, only: diagnostics_type
    use gs2_diagnostics, only: heating
    implicit none
    type(diagnostics_type), intent (inout) :: gnostics
    call heating(gnostics%istep, gnostics%navg, h, hk)
  end subroutine calculate_heating

  !> FIXME : Add documentation  
  subroutine write_heating(gnostics)
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use gs2_diagnostics, only: do_write_heating
    use gs2_io, only: nc_write_heating
    implicit none
    type(diagnostics_type), intent (inout) :: gnostics

    if (proc0) call nc_write_heating(gnostics%file_id, gnostics%nout, h, hk)
    
    gnostics%current_results%species_heating = h%imp_colls
    gnostics%current_results%species_heating_avg = gnostics%current_results%species_heating_avg + &
         h%imp_colls*(gnostics%user_time - gnostics%user_time_old)
    
    if (gnostics%write_ascii .and. proc0) then
      call do_write_heating(gnostics%user_time, gnostics%ascii_files%heat, gnostics%ascii_files%heat2, h)
    end if
  end subroutine write_heating
end module diagnostics_heating
