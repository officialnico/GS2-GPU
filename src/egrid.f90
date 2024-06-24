!> FIXME : Add documentation
module egrid
  ! By Tomo Tatsuno, Aug 2005
  ! Improved accuracy and speed and maximum number of energy grid points

  implicit none

  private

  public :: setvgrid, setvgrid_genquad

contains

  !> FIXME : Add documentation
  subroutine setvgrid (vcut_in, negrid, epts, wgts, nesub_in)
    use constants, only: pi => dpi
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, f0_sdanalytic, f0_maxwellian, calculate_f0_arrays, f0_values

    implicit none
    integer, intent (in) :: negrid
    real, intent (in) :: vcut_in
    integer, intent (in) :: nesub_in
    real, dimension(:,:), intent (out) :: epts, wgts
    integer:: nesub, is
    real:: vcut

    do is = 1,nspec
       select case (spec(is)%f0type)
       case(f0_maxwellian)
          vcut = vcut_in
          nesub = nesub_in
       case(f0_sdanalytic)
          vcut = 1.0
          nesub = negrid
       case default
          vcut = vcut_in
          nesub = nesub_in
       end select

       ! get grid points in v up to vcut (epts is not E yet)
       call get_legendre_grids_from_cheb (0., vcut, epts(:nesub,is), wgts(:nesub,is))

       ! change from v to E
       epts(:nesub,is) = epts(:nesub,is)**2

       ! absorb exponential and volume element in weights
       !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
       !> GW: F0 now gets multiplied later
       wgts(:nesub, is) = wgts(:nesub, is) * epts(:nesub,is) * pi

       if (negrid > nesub) then

          ! get grid points in y = E - vcut**2 (epts not E yet)
          call get_laguerre_grids (epts(nesub+1:,is), wgts(nesub+1:,is))

          ! change from y to E
          epts(nesub+1:,is) = epts(nesub+1:,is) + vcut**2

          ! Do not absorb exponential in weights, since we'll be multiplying by f0 later
          wgts(nesub+1:, is) = wgts(nesub+1:, is) * exp(epts(nesub+1:, is) - vcut**2) &
               * pi * 0.5 * sqrt(epts(nesub+1:, is))
       end if
    end do

    call calculate_f0_arrays(epts)
    wgts = wgts * f0_values
  end subroutine setvgrid

  !> FIXME : Add documentation
  subroutine setvgrid_genquad (negrid, epts, wgts)
    use genquad, only: get_quadrature_rule
    use constants, only: pi => dpi
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, calculate_f0_arrays, eval_f0, set_current_f0_species, f0_maxwellian
    implicit none

    integer, intent (in) :: negrid
    real, dimension(:,:), intent (out) :: epts, wgts
    integer :: is
    logical :: upper_bound_is_inf

    do is = 1, nspec

       ! Call the function that sets the species number so that eval_f0
       ! takes only one argument
       call set_current_f0_species(is)

       ! If Maxwellian, use an semi-infinite domain
       ! Otherwise, use a bounded domain with vmax = vref_s.
       upper_bound_is_inf = spec(is)%f0type == f0_maxwellian
       call get_quadrature_rule(eval_f0, negrid, 0.0, 1.0, &
            epts(:, is), wgts(:, is), .false., upper_bound_is_inf)

       ! change from v to E
       epts(:, is) = epts(:, is)**2

       ! absorb exponential and volume element in weights
       !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
       ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
       ! See eq. 4.12 of M. Barnes's thesis
       wgts(:, is) = wgts(:, is) * pi * epts(:, is)

    end do

    call calculate_f0_arrays(epts)
    !Note we don't do wgts = wgts * f0_values unlike setvgrid
  end subroutine setvgrid_genquad

end module egrid
