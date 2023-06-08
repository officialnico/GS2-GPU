!> A module for printing summary information to stdout during
!! the simulation
module diagnostics_printout
  implicit none

  private

  public :: print_line
  public :: write_line
  public :: print_flux_line
  public :: write_flux_line
  public :: write_phi_fluxes_to_unit, write_fluxes_to_unit

  !> Format string for writing fluxes to screen/text
  character(*), parameter :: flux_format = &
       "('t= ', e17.10, ' <', a, '**2>= ', e13.6, ' ', a, ':', *(1x, e13.6))"

contains

  !> Print out growth rates to screen.
  !! Doesn't quite conform to the principle of having all
  !! useful results stored in gnostics... takes stuff
  !! straight from diagnostics_omega
  subroutine print_line(gnostics)
    use file_utils, only: stdout_unit
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    call output_line(gnostics, stdout_unit)
  end subroutine print_line

  !> Print growth rates to run_name.new.out
  subroutine write_line(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    if (.not. gnostics%write_ascii) return
    if (gnostics%create) return
    call output_line(gnostics, gnostics%ascii_files%out)
  end subroutine write_line

  !> Write growth rates to given unit
  subroutine output_line(gnostics, outunit)
    use mp, only: proc0, sum_reduce
    use kt_grids, only: naky, ntheta0, aky, akx, theta0
    use gs2_time, only: woutunits
    use volume_averages, only: average_theta
    use diagnostics_omega, only: omegahist, omega_average
    use diagnostics_config, only: diagnostics_type

    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    !> Pre-existing opened unit
    integer, intent(in) :: outunit
    real, dimension (ntheta0, naky) :: phitot
    integer :: ik, it
    phitot=0.
    !call average_theta(phinew, phinew, phitot, gnostics%distributed)
    call phinorm (phitot)
!WARNING: Phitot here has not been set
    if (gnostics%distributed) call sum_reduce(phitot, 0)

    if(.not.(proc0 .and. gnostics%writing)) return
    write (unit=outunit, fmt="('time=',1pe12.4)") gnostics%user_time
    do ik = 1, naky
       do it = 1, ntheta0
          write (unit=outunit, fmt="('t=',1pe11.4,' ky=',1pe9.2, ' kx=',1pe9.2, &
               & ' om= ',e9.2,1x,e9.2,' omav= ',e9.2,1x,e9.2, &
               & ' phtot= ',e9.2,' theta0= ',1pe9.2)") &
               gnostics%user_time, aky(ik), akx(it), &
               real( omegahist(mod(gnostics%istep,gnostics%navg),it,ik)*woutunits(ik)), &
               aimag(omegahist(mod(gnostics%istep,gnostics%navg),it,ik)*woutunits(ik)), &
               real( omega_average(it,ik)*woutunits(ik)), &
               aimag(omega_average(it,ik)*woutunits(ik)), &
               phitot(it,ik), theta0(it,ik)
       end do
    end do
    write (outunit,*)

  contains
      !> Compute normalised flux-surface averaged field
      subroutine phinorm (phitot)
        use fields_arrays, only: phinew, aparnew, bparnew
        use theta_grid, only: theta
        use kt_grids, only: naky, ntheta0
        use constants, only: pi
        use integration, only: trapezoidal_integration
        implicit none
        real, dimension (:,:), intent (out) :: phitot
        integer :: ik, it
        do ik = 1, naky
           do it = 1, ntheta0
              phitot(it,ik) =0.5/pi &
               *(trapezoidal_integration(theta, (abs(phinew(:,it,ik))**2 + abs(aparnew(:,it,ik))**2 &
                      + abs(bparnew(:,it,ik))**2)))
           end do
        end do
      end subroutine phinorm
  end subroutine output_line

  !> Print instantaneous heat fluxes to run_name.new.out
  subroutine write_flux_line(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    if (.not. gnostics%ascii_files%write_to_out) return
    if (gnostics%create) return
    call output_flux_line(gnostics, gnostics%ascii_files%out)
  end subroutine write_flux_line

  !> Print instantaneous heat fluxes to screen
  subroutine print_flux_line(gnostics)
    use file_utils, only: stdout_unit
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    call output_flux_line(gnostics, stdout_unit)
  end subroutine print_flux_line

  !> Write instantaneous heat fluxes to given unit
  !>
  !> Used for new diagnostics only
  subroutine output_flux_line(gnostics, outunit)
    use run_parameters, only: has_phi, has_apar, has_bpar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    integer, intent(in) :: outunit
    if (has_phi) then
      call write_phi_fluxes_to_unit(outunit, gnostics%user_time, &
           gnostics%current_results%phi2, &
           gnostics%current_results%species_es_heat_flux, &
           gnostics%current_results%species_energy_exchange)
    end if
    if (has_apar) then
      call write_fluxes_to_unit(outunit, gnostics%user_time, &
           gnostics%current_results%apar2, "apar", &
           gnostics%current_results%species_apar_heat_flux)
    end if
    if (has_bpar) then
      call write_fluxes_to_unit(outunit, gnostics%user_time, &
           gnostics%current_results%bpar2, "bpar", &
           gnostics%current_results%species_bpar_heat_flux)
    end if
  end subroutine output_flux_line

  !> Write instantaneous fluxes from electrostatic potential to screen
  !> or file
  subroutine write_phi_fluxes_to_unit(unit, time, phi2, &
       heat_fluxes, energy_exchange, part_fluxes, mom_fluxes, &
       lfmom_fluxes, vflux1_avg, use_lowflow)
    use optionals, only: get_option_with_default
    !> An existing, open unit
    integer, intent(in) :: unit
    !> Current simulation time
    real, intent(in) :: time
    !> Averaged amplitude squared of electrostatic potential
    real, intent(in) :: phi2
    !> Heat flux due to electrostatic potential
    real, dimension(:), intent(in) :: heat_fluxes
    !> Energy exchange -- see [[dist_fn:eexchange]]
    real, dimension(:), intent(in) :: energy_exchange
    !> Particle flux due to electrostatic potential
    real, dimension(:), intent(in), optional :: part_fluxes
    !> Momentum flux due to electrostatic potential
    real, dimension(:), intent(in), optional :: mom_fluxes
    !> "Low flow" momentum flux due to electrostatic potential
    real, dimension(:), intent(in), optional :: lfmom_fluxes
    !> FIXME: Add documentation (some other low-flow flux term)
    real, dimension(:), intent(in), optional :: vflux1_avg
    !> If true, write the low-flow momentum fluxes
    logical, intent(in), optional :: use_lowflow

    call write_fluxes_to_unit(unit, time, phi2, "phi", heat_fluxes, part_fluxes)

    if (present(mom_fluxes)) then
      write (unit=unit, fmt=flux_format) time, "phi", phi2, "mom fluxes", mom_fluxes
    end if

    write (unit=unit, fmt=flux_format) time, "phi", phi2, "energy exchange", energy_exchange

    ! Bail if we don't want to write the low-flow momentum fluxes
    if (.not. get_option_with_default(use_lowflow, .false.)) return

    if (present(use_lowflow) .and. present(lfmom_fluxes) .and. present(vflux1_avg)) then
      ! This one writes two separate fluxes on the same line, so can't
      ! reuse format string. Intel also requires unlimited-format '*'
      ! to be last item, so we split over two statements
      write (unit=unit, fmt="('t= ',e17.10,' <phi**2>= ',e13.6, &
           & ' lfmom fluxes:', *(1x,e13.6))", advance="no") time, phi2, lfmom_fluxes
      write (unit=unit, fmt="(' lfvflx1:', *(1x,e13.6))") vflux1_avg
    end if
  end subroutine write_phi_fluxes_to_unit

  !> Write instantaneous heat and particle fluxes to screen/file
  !>
  !> Used directly for apar/bpar, but phi should use [[write_phi_fluxes_to_unit]]
  subroutine write_fluxes_to_unit(unit, time, field, field_name, heat_fluxes, part_fluxes)
    !> An existing, open unit
    integer, intent(in) :: unit
    !> Current simulation time
    real, intent(in) :: time
    !> Average amplitude squared of the parallel magnetic potential
    real, intent(in) :: field
    !> Name of the field
    character(len=*), intent(in) :: field_name
    !> Heat flux due to A//
    real, dimension(:), intent(in) :: heat_fluxes
    !> Particle flux due to A//
    real, dimension(:), intent(in), optional :: part_fluxes

    ! Spelling of "fluxes" for different fields
    character(len=6) :: flux_name
    select case (field_name)
      case("apar"); flux_name = "mluxes"
      case("bpar"); flux_name = "bluxes"
      case default; flux_name = "fluxes"
    end select

    write (unit=unit, fmt=flux_format) time, field_name, field, "heat " // flux_name, heat_fluxes

    if (present(part_fluxes)) then
      write (unit=unit, fmt=flux_format) time, field_name, field, "part " // flux_name, part_fluxes
    end if
  end subroutine write_fluxes_to_unit
end module diagnostics_printout
