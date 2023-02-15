module test_le_grids_integrate_utils
  implicit none
  private


  interface mean
     module procedure mean_4c
     module procedure mean_r
     module procedure mean_2r
  end interface mean

  public :: mean
contains

  complex function mean_4c(array) result(mean)
    implicit none
    complex, dimension(:,:,:,:), intent(in) :: array
    mean = sum(array)/size(array)
  end function mean_4c

  real function mean_r(array) result(mean)
    implicit none
    real, dimension(:), intent(in) :: array
    mean = sum(array)/size(array)
  end function mean_r

  real function mean_2r(array) result(mean)
    implicit none
    real, dimension(:,:), intent(in) :: array
    mean = sum(array)/size(array)
  end function mean_2r

end module test_le_grids_integrate_utils

!> Possible integration tests -- check check at all theta and mean
!> 1. Int(1) = 1.0 | Integral of normalised Maxwellian, normalised to integral of Maxwellian
!> 2, Int(0) = 0.0 | Trivial sanity check
!> 3. Int(forbid? 1.0 : 0.0) = 0.0 | Check that no contribution from the forbidden region
!> 4. Int(forbid? 0.0 : 1.0) = 1.0 | Same as 1 plus no contribution from the forbidden region
!> 5. Int(a) = a | Generalise 1 & 2
!> 6. Int(x^n J_0(x)) = .... for n>0 | Just vperp dependence, should be expressible in terms of other Bessel functions
!> 7. Int(Poly(vperp2))
!> 8. Int(Poly(lambda))
!> 9. Int(Poly(vpar))?
!>
!> In GS2 the integral over energy of a quantity f is given by
!> Integral(f*v^2*exp(-v^2) dv, {v,0,infinity})/Integral(v^2*exp(-v^2) dv, {v,0,infinity})
!> i.e.
!> Integral(f*v^2*exp(-v^2) dv, {v,0,infinity})*4/sqrt(pi)
!>
!> The integral over lambda is of a quantity f is given by
!> Integral(f sqrt(lambda*Bmag) d alpha, {alpha, 0, pi})
!> where alpha is the pitch angle, alpha = arcsin(sqrt(lambda*bmag))
!>
!> We use quadrature rules such that we expect accurate integral results
!> with relatively small grids. Our integrals are just sums of the quantity
!> weighted with the velocity space weights, wl and w, (lambda and energy).
!> We expect the sum of the energy weights, w, to be 0.25, the sum of lambda
!> weights to be 2.0 such that the integral of 1 is exactly 1. This is given by
!> Sum_{sigma,lambda,energy}(w*wl*sigma) at fixed theta. This allows us to
!> estimate the error in the integral weights (very roughly).
!>
!> Out of the default test cases the nesuper version has the most
!> accurate results. The majority of tests then pass at ~1e-8. This
!> level does not reduce with increasing ngauss or nesuper/nesub
!> indicating we are perhaps limited by the trapped pitch
!> angles. Switching to new_trap_int = T reduces the error by 2-4
!> orders of magnitude at fixed resolution. (The same reduction can be
!> achieved with ntheta ~= 132 and new_trap_int = F) The error then
!> reduceds further with increasing ngauss/negrid.  Setting nesuper
!> =12, nesub = 24, ngauss = 10 leads to errors in the range
!> 1e-15. Switching to radau_gauss_grid = F drops some of the tests
!> further to ~4e-16 (close to the minimum possible > 0). Changing
!> ntheta doesn't significantly improve most of these tests with
!> new_trap_int active and can actually increase error.  The
!> integration is a weighted sum over all velocity space. This sum
!> will have associated round off error which is likely setting a
!> lower limit on the error, and this limit may be expected to
!> increase as the number of velocity space points increases. We could
!> perhaps consider using compensated summation instead of direct
!> summation in order to further minimise this. In practice, however,
!> we are probably more than happy with an error in the range 1e-12.
!>
!> The mod vpa test case is an interesting example as this is the only
!> (non-failing) test which doesn't pass when turning on new_trap_int
!> (with radau_gauss_grid = F).
program test_le_grids_integrate
  use test_le_grids_integrate_utils, only: mean
  use unit_tests
  use le_grids, only: integrate_moment, energy, al, forbid, wl, w
  use dist_fn, only: init_dist_fn
  use dist_fn_arrays, only: gnew
  use file_utils, only: init_file_utils, input_unit_exist
  use mp, only: init_mp, finish_mp, proc0, broadcast, barrier
  use species, only: nspec
  use kt_grids, only: naky, ntheta0
  use theta_grid, only: ntgrid, bmag
  use gs2_layouts, only: g_lo, ie_idx, il_idx, is_idx
  use constants, only: pi
  use iso_fortran_env, only: output_unit

  implicit none

  logical :: dummy
  complex, dimension (:,:,:,:), allocatable :: integrate_results
  real, dimension(:, :,:), allocatable :: test_sum
  integer :: ig, il, ie, is, iglo
  real :: tolerance = 1.0e-8
  logical :: agrees
  integer :: in_file
  logical :: exist
  logical :: reduce_via_mean = .false.
  logical :: check_all_values = .false.
  logical :: run_failing_tests = .false.
  real, parameter :: expected_lambda = 2.0, expected_energy = 0.25
  real, dimension(:), allocatable :: lambda_error, energy_error
  real, dimension(:, :), allocatable :: mixed_error, mixed_value
  real :: relevant_value, tmp
  character(len=22) :: error_name
  namelist /le_grids_integrate/ tolerance, check_all_values, run_failing_tests

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy)

  call init_dist_fn
  call announce_module_test('le_grids_integrate')
  allocate(integrate_results(-ntgrid:ntgrid, ntheta0, naky, nspec))

  if(proc0) then
     in_file = input_unit_exist('le_grids_integrate', exist)
     if (exist) read(in_file, nml = le_grids_integrate)
  end if
  call broadcast(tolerance)
  call broadcast(check_all_values)
  call broadcast(run_failing_tests)

  !First get estimates of the error on the integration weights.
  if (proc0) then
     allocate(lambda_error(-ntgrid:ntgrid))
     lambda_error = 0.0
     do ig = -ntgrid, ntgrid
        lambda_error(ig) = expected_lambda - sum(wl(ig, :))
     end do
     !Introduce factor 2 here to account for sigma doubling
     lambda_error = 2*lambda_error

     allocate(energy_error(nspec))
     energy_error = 0.0
     do is = 1, nspec
        energy_error(is) = expected_energy - sum(w(:,is))
     end do

     allocate(mixed_error(-ntgrid:ntgrid, nspec))
     allocate(mixed_value(size(wl(0,:)), size(w(:,1))))
     mixed_error = 0.0
     do is = 1, nspec
        do ig = -ntgrid, ntgrid
           mixed_value = 0.0
           do ie = 1, size(w(:,1))
              mixed_value(:, ie) = wl(ig,:)*w(ie, is)
           end do
           mixed_error(ig, is) = expected_lambda*expected_energy - sum(mixed_value)
        end do
     end do

     !Introduce factor 2 here to account for sigma doubling
     mixed_error = 2 * mixed_error

     write(output_unit,'("Lambda weights errors, max/mean/mean(|error|) :",3(E14.6E3," "))') maxval(abs(lambda_error)),mean(lambda_error), mean(abs(lambda_error))
     write(output_unit,'("Energy weights errors, max/mean/mean(|error|) :",3(E14.6E3," "))') maxval(abs(energy_error)),mean(energy_error), mean(abs(energy_error))

     write(output_unit,'("Mixed weights errors, max/mean/mean(|error|)  :",3(E14.6E3," "))') maxval(abs(mixed_error)),mean(mixed_error), mean(abs(mixed_error))
     write(output_unit,'("")')

     ! Crudely we probably expect integrals to be no more accurate than mean(mixed_error)
     ! for theta = 0 ror max(abs(error)) for checking all so lets print a warning about
     ! that here.

     if (check_all_values) then
        relevant_value = maxval(abs(mixed_error))
        error_name = 'Max(abs(mixed_error))'
     else
        relevant_value = abs(mean(mixed_error))
        error_name = 'Abs(mean(mixed_error))'
     end if

     if (relevant_value >= tolerance * 2) then
        write(output_unit, '("Warning : ",A," >= tolerance (",E14.6E3,") so expect tests may fail.")') error_name, tolerance
     else if (relevant_value >= tolerance / 2) then
        write(output_unit, '("Warning : ",A," similar to tolerance (",E14.6E3,") so tests may fail.")') error_name, tolerance
     else
        write(output_unit, '(A," < tolerance (",E14.6E3,") so expect tests should pass.")') error_name, tolerance
     end if
  end if
  call barrier

  ! Default to passing
  agrees = .true.

  agrees = integrate_zero_test() .and. agrees
  agrees = integrate_zero_forbid_test() .and. agrees
  agrees = integrate_one_test() .and. agrees
  agrees = integrate_one_forbid_test() .and. agrees
  agrees = integrate_one_signed_test() .and. agrees
  agrees = integrate_sqrt_e_test() .and. agrees
  agrees = integrate_inverse_sqrt_e_test() .and. agrees
  agrees = integrate_e_test() .and. agrees
  agrees = integrate_pitch_angle_test() .and. agrees
  agrees = integrate_one_minus_sin_pitch_angle_squared_test() .and. agrees
  !The following test is surprisingly hard to get converged to
  !a reasonable degree. It does not appear to converge unless
  !new_trap_int = T and even then it only converges within 1e-6 if
  !ngauss >= 30. We might need to explore why this is the case.
  if (run_failing_tests) agrees = integrate_pitch_angle_squared_test() .and. agrees
  agrees = integrate_pitch_angle_v_test() .and. agrees
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vperp_test() .and. agrees
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vperp2_test() .and. agrees
  agrees = integrate_vpa_test() .and. agrees
  agrees = integrate_mod_vpa_test() .and. agrees
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vpa2_test() .and. agrees
  !The following tests the mod cell centred vpa integral. We might expect this
  !to agree with the mod vpa value however this currently seems to fail with
  !a relatively large disagreement. The disagreement reduces as ntheta increases
  !but still fails to pass even for excessively high resolution. This may indicate
  !a potential issue.
  if (run_failing_tests) agrees = integrate_mod_vpac_test() .and. agrees
  !The following test fails with large errors
  if (run_failing_tests) agrees = integrate_vperp_inverse_test() .and. agrees
  agrees = integrate_polynomial_test() .and. agrees

  call process_test(agrees,'le_grids_integrate')

  call close_module_test('le_grids_integrate')
  call finish_mp

contains
  complex function reduce_to_single_value(array) result(val)
    implicit none
    complex, dimension(-ntgrid:,:,:,:), intent(in) :: array
    if(reduce_via_mean) then
       val = mean(array)
    else
       val = array(0,1,1,1)
    end if
  end function reduce_to_single_value

  logical function integrate_and_test(name, expected_value) result(agrees)
    implicit none
    character(len=*), intent(in) :: name
    real, intent(in) :: expected_value
    real :: integral
    logical :: dummy_flag
    call announce_test('le_grids_integrate :: '//name)
    integrate_results = 0.0
    call integrate_moment(gnew, integrate_results, all = .true.)

    if (check_all_values) then
       agrees = all(abs(integrate_results-expected_value) <= tolerance)
       ! Hacky representation of reduced integral result
       integral = maxval(abs(integrate_results-expected_value)) + expected_value
       if(.not. agrees .and. proc0) print*,"Failure : ",name," : Got ",integral," but expected ",expected_value,"(disagreement ~ ",maxval(abs(integrate_results-expected_value)),"). Note the quoted obtained value is not accurate here."

    else
       integral = reduce_to_single_value(integrate_results)
       agrees = abs(integral - expected_value) <= tolerance
       if(.not. agrees .and. proc0) print*,"Failure : ",name," : Got ",integral," but expected ",expected_value,"(disagreement ~ ",abs(integral-expected_value),")"
    end if

    dummy_flag = .true.
    call process_check(dummy_flag, agrees, "TEST CASE "//name)
  end function integrate_and_test

  logical function integrate_zero_test() result(agrees)
    implicit none
    ! Integral[0 * v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 0.0
    character(len=*), parameter :: name = '0'
    gnew = 0.0
    agrees = integrate_and_test(name, expected_value)
  end function integrate_zero_test

  logical function integrate_one_test() result(agrees)
    implicit none
    ! Integral[1 * v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 1.0
    character(len=*), parameter :: name = '1'
    gnew = 1.0
    agrees = integrate_and_test(name, expected_value)
  end function integrate_one_test

  logical function integrate_one_signed_test() result(agrees)
    implicit none
    ! Sum_sigma{ sigma *Integral[1 * v^2 exp[-v^2],{v,0,infty}]}/
    ! Sum_sigma{Integral[v^2 exp[-v^2],{v,0,infty}]}
    real, parameter :: expected_value = 0.0
    character(len=*), parameter :: name = '1 signed'
    gnew(:,1,:) = 1.0
    gnew(:,2,:) = -1.0
    agrees = integrate_and_test(name, expected_value)
  end function integrate_one_signed_test

  logical function integrate_sqrt_e_test() result(agrees)
    implicit none
    ! Integral[v v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 2.0/sqrt(pi)
    character(len=*), parameter :: name = 'sqrt(energy)'
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       gnew(:,:,iglo) = sqrt(energy(ie, is))
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_sqrt_e_test

  logical function integrate_inverse_sqrt_e_test() result(agrees)
    implicit none
    ! Integral[v^2 exp[-v^2] / v,{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 2.0/sqrt(pi)
    character(len=*), parameter :: name = 'inverse sqrt(energy)'
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       gnew(:,:,iglo) = 1.0/sqrt(energy(ie, is))
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_inverse_sqrt_e_test

  logical function integrate_e_test() result(agrees)
    implicit none
    ! Integral[v^2 v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 3.0/2.0
    character(len=*), parameter :: name = 'energy'
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       gnew(:,:,iglo) = energy(ie,is)
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_e_test

  logical function integrate_one_forbid_test() result(agrees)
    implicit none
    ! Integral[1 * v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 1.0
    character(len=*), parameter :: name = '1 only in allowed'
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          gnew(ig,:,iglo) = 1.0
       end do
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_one_forbid_test

  logical function integrate_zero_forbid_test() result(agrees)
    implicit none
    ! Integral[0 * v^2 exp[-v^2],{v,0,infty}]/Integral[v^2 exp[-v^2],{v,0,infty}]
    real, parameter :: expected_value = 0.0
    character(len=*), parameter :: name = '0 only in allowed'
    gnew = 1.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          gnew(ig,:,iglo) = 0.0
       end do
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_zero_forbid_test

  logical function integrate_pitch_angle_test() result(agrees)
    implicit none
    ! Integral[alpha*sin[alpha],{alpha,0,pi}]/Integral[sin[alpha],{alpha,0,pi}]
    real, parameter :: expected_value = pi / 2.0
    character(len=*), parameter :: name = 'pitch angle'
    real :: pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          pitch_angle = asin(sqrt(al(il)*bmag(ig)))
          gnew(ig,1,iglo) = pitch_angle
          gnew(ig,2,iglo) = pi-pitch_angle
       end do
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_pitch_angle_test

  logical function integrate_one_minus_sin_pitch_angle_squared_test() result(agrees)
    implicit none
    ! Integral[sqrt(1-sin^2[alpha])sin[alpha],{alpha,0,pi}]/Integral[sin[alpha],{alpha,0,pi}]
    real, parameter :: expected_value = 1.0/2.0
    character(len=*), parameter :: name = 'one minus sin squared pitch angle'
    real :: sin_squared_pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          sin_squared_pitch_angle = al(il)*bmag(ig)
          gnew(ig,1,iglo) = sqrt(1 - sin_squared_pitch_angle)
          gnew(ig,2,iglo) = sqrt(1 - sin_squared_pitch_angle)
       end do
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_one_minus_sin_pitch_angle_squared_test

  logical function integrate_pitch_angle_squared_test() result(agrees)
    implicit none
    ! Integral[alpha^2*sin[alpha],{alpha,0,pi}]/Integral[sin[alpha],{alpha,0,pi}]
    real, parameter :: expected_value = (pi**2 - 4.0)/2.0
    character(len=*), parameter :: name = 'pitch angle squared'
    real :: pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          pitch_angle = asin(sqrt(al(il)*bmag(ig)))
          gnew(ig,1,iglo) = pitch_angle**2
          gnew(ig,2,iglo) = (pi-pitch_angle)**2
       end do
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_pitch_angle_squared_test

  logical function integrate_pitch_angle_v_test() result(agrees)
    implicit none
    ! Integral[alpha*sin[alpha]*v^3*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = sqrt(pi)
    character(len=*), parameter :: name = 'pitch angle times v'
    real :: pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          pitch_angle = asin(sqrt(al(il)*bmag(ig)))
          gnew(ig,1,iglo) = pitch_angle
          gnew(ig,2,iglo) = (pi-pitch_angle)
       end do
       gnew(:, :, iglo) = gnew(:, :, iglo)*sqrt(energy(ie,is))
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_pitch_angle_v_test

  logical function integrate_vperp_test() result(agrees)
    implicit none
    ! vperp = sqrt(energy) * sin(alpha) = sqrt(energy*lambda*Bmag)
    ! Integral[vperp*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = sqrt(pi)/2.0
    character(len=*), parameter :: name = 'vperp'
    real :: sin_pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          sin_pitch_angle = sqrt(al(il)*bmag(ig))
          gnew(ig,1,iglo) = sin_pitch_angle
          gnew(ig,2,iglo) = sin_pitch_angle
       end do
       gnew(:, :, iglo) = gnew(:, :, iglo)*sqrt(energy(ie,is))
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_vperp_test

  logical function integrate_vperp_inverse_test() result(agrees)
    implicit none
    ! vperp = sqrt(energy) * sin(alpha) = sqrt(energy*lambda*Bmag)
    ! Integral[sin[alpha]*v^2*exp[-v^2]/sqrt(vperp2),{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    ! Which becomes
    ! Integral[v*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    ! Note that our integrand contains 1/sin(alpha) and our
    ! pitch angle domain ranges from 0 to pi. Whilst we don't have
    ! pitch angle points exactly on the boundaries (i.e. at 0 and pi)
    ! they can approach the boundaries and hence 1/sin(alpha) can
    ! become large. This results in this test case showing relatively large
    ! errors and slow convergence (first order only).
    real, parameter :: expected_value = sqrt(pi)
    character(len=*), parameter :: name = 'vperp inverse'
    real :: sin_pitch_angle
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          if(forbid(ig, il)) cycle
          sin_pitch_angle = sqrt(al(il)*bmag(ig))
          gnew(ig,1,iglo) = 1.0/sin_pitch_angle
          gnew(ig,2,iglo) = 1.0/sin_pitch_angle
       end do
       gnew(:, :, iglo) = gnew(:, :, iglo)/sqrt(energy(ie,is))
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_vperp_inverse_test

  logical function integrate_vperp2_test() result(agrees)
    use dist_fn_arrays, only: vperp2
    implicit none
    ! vperp = sqrt(energy) * sin(alpha) = sqrt(energy*lambda*Bmag)
    ! vperp2 = energy*lambda*Bmag
    ! Integral[vperp^2*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = 1.0
    character(len=*), parameter :: name = 'vperp^2'
    gnew = 0.0d0
    gnew(:,1,:) = vperp2
    gnew(:,2,:) = vperp2
    agrees = integrate_and_test(name, expected_value)
  end function integrate_vperp2_test

  logical function integrate_vpa_test() result(agrees)
    use dist_fn_arrays, only: vpa
    implicit none
    ! vpa = +/- sqrt(energy*(1-lambda*bmag)) = v * sqrt(1-sin^2(alpha))
    ! Integral[vpa*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = 0.0
    character(len=*), parameter :: name = 'vpa'
    gnew = vpa
    agrees = integrate_and_test(name, expected_value)
  end function integrate_vpa_test

  logical function integrate_mod_vpac_test() result(agrees)
    use dist_fn_arrays, only: vpac
    implicit none
    ! vpa = +/- sqrt(energy*(1-lambda*bmag)) = v * sqrt(1-sin^2(alpha))
    ! vpac = cell centred vpa
    ! Integral[vpa*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = 1.0/sqrt(pi)
    character(len=*), parameter :: name = 'mod vpac'
    gnew = abs(vpac)
    ! As vpac is cell centred the ntgrid location does not get a value
    ! so here we copy in the -ntgrid value to ensure that if we check
    ! all locations it doesn't automatically fail as one theta point
    ! was not set non-zero.
    gnew(ntgrid,:,:) = gnew(-ntgrid,:,:)
    agrees = integrate_and_test(name, expected_value)
  end function integrate_mod_vpac_test

  logical function integrate_mod_vpa_test() result(agrees)
    use dist_fn_arrays, only: vpa
    implicit none
    ! vpa = +/- sqrt(energy*(1-lambda*bmag)) = v * sqrt(1-sin^2(alpha))
    ! Integral[vpa*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = 1.0/sqrt(pi)
    character(len=*), parameter :: name = 'mod vpa'
    gnew = abs(vpa)
    agrees = integrate_and_test(name, expected_value)
  end function integrate_mod_vpa_test

  logical function integrate_vpa2_test() result(agrees)
    use dist_fn_arrays, only: vpa
    implicit none
    ! vpa^2 = (energy*(1-lambda*bmag)) = v^2 * (1-sin^2(alpha))
    ! Integral[vpa*sin[alpha]*v^2*exp[-v^2],{alpha,0,pi},{v,0,infty}]/
    ! Integral[sin[alpha]*v^2*exp(-v^2),{alpha,0,pi}{v,0,infty}]
    real, parameter :: expected_value = 1.0/2.0
    character(len=*), parameter :: name = 'vpa2'
    gnew = vpa*vpa
    agrees = integrate_and_test(name, expected_value)
  end function integrate_vpa2_test

  logical function integrate_polynomial_test() result(agrees)
    use le_grids, only: nesub, vcut
    implicit none
    ! Sum_sigma{ sigma *Integral[sqrt[pi]*v^N * (exp[v^2]/v^2) * v^2 exp[-v^2],{v,0,infty}]}/
    ! Sum_sigma{Integral[v^2 exp[-v^2],{v,0,infty}]}
    ! This is a slightly odd test which attempts to check the standard
    ! Legendre quadrature used in the lower energy domain
    ! (i.e. between 0 and Ecut=vcut^2).  We expect the Gauss-Legendre
    ! to exactly capture integrals of polynomials of degree 2*nesub-1
    ! or lower. This region is formulated in terms of velocity, hence
    ! the use of sqrt(energy). The exp(E)*sqrt(pi)/E factor comes from
    ! cancelling the contribution to the weights which comes from the
    ! fact our normal integrals require us to absorb the inverse of
    ! these factors (see eqn 28/29 of Numata et al 2010, the AstroGK
    ! paper).
    real :: expected_value
    character(len=*), parameter :: name = 'exp poly'
    integer, parameter :: poly_order = 1
    expected_value = 2 * expected_lambda * (vcut**(poly_order+1))/(poly_order+1)
    gnew = 0.0d0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       if(ie > nesub) cycle
       gnew(:,:,iglo) = sqrt(energy(ie,is))**poly_order*&
            exp(energy(ie, is))*sqrt(pi)/energy(ie,is)
    end do
    agrees = integrate_and_test(name, expected_value)
  end function integrate_polynomial_test

end program test_le_grids_integrate
