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
program test_le_grids_integrate
  use unit_tests
  use le_grids, only: integrate_moment, energy, al, forbid
  use dist_fn, only: init_dist_fn
  use dist_fn_arrays, only: gnew
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, input_unit_exist
  use species, only: nspec
  use kt_grids, only: naky, ntheta0
  use theta_grid, only: ntgrid, bmag
  use gs2_layouts, only: g_lo, ie_idx, il_idx, is_idx
  use constants, only: pi

  implicit none
  logical :: dummy
  complex, dimension (:,:,:,:), allocatable :: integrate_results
  real, dimension(:, :,:), allocatable :: test_sum
  integer :: ig, il, ie, is, iglo
  real :: tolerance = 1.0e-8
  logical :: agrees
  integer :: in_file
  logical :: exist
  logical :: global_pass
  logical :: reduce_via_mean = .false.
  logical :: check_all_values = .false.
  logical :: run_failing_tests = .false.
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

  global_pass = .true.

  agrees = integrate_zero_test()
  agrees = integrate_zero_forbid_test()
  agrees = integrate_one_test()
  agrees = integrate_one_forbid_test()
  agrees = integrate_one_signed_test()
  agrees = integrate_sqrt_e_test()
  agrees = integrate_e_test()
  agrees = integrate_pitch_angle_test()
  !The following test is surprisingly hard to get converged to
  !a reasonable degree. It does not appear to converge unless
  !new_trap_int = T and even then it only converges within 1e-6 if
  !ngauss >= 30. We might need to explore why this is the case.
  if (run_failing_tests) agrees = integrate_pitch_angle_squared_test()
  agrees = integrate_pitch_angle_v_test()
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vperp_test()
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vperp2_test()
  agrees = integrate_vpa_test()
  !The following test passes with default settings but has been seen to fail
  !at higher resolution in certain cases.
  agrees = integrate_mod_vpa_test()
  !The following test doesn't currently converge at reasonable
  !resolution.
  if (run_failing_tests) agrees = integrate_vpa2_test()
  !The following tests the mod cell centred vpa integral. We might expect this
  !to agree with the mod vpa value however this currently seems to fail with
  !a relatively large disagreement. The disagreement reduces as ntheta increases
  !but still fails to pass even for excessively high resolution. This may indicate
  !a potential issue.
  if (run_failing_tests) agrees = integrate_mod_vpac_test()

  call process_test(global_pass,'le_grids_integrate')

  call close_module_test('le_grids_integrate')
  call finish_mp

contains
  complex function mean(array)
    implicit none
    complex, dimension(:,:,:,:), intent(in) :: array
    mean = sum(array)/size(array)
  end function mean

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

    call process_check(global_pass,agrees,"TEST CASE "//name)
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

end program test_le_grids_integrate
