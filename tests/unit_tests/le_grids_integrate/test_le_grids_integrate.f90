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
  namelist /le_grids_integrate/ tolerance

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
  global_pass = .true.

  agrees = integrate_zero_test()
  agrees = integrate_zero_forbid_test()
  agrees = integrate_one_test()
  agrees = integrate_one_forbid_test()
  agrees = integrate_sqrt_e_test()
  agrees = integrate_e_test()
  agrees = integrate_pitch_angle_test()

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
    integral = reduce_to_single_value(integrate_results)
    agrees = abs(integral - expected_value) <= tolerance
    if(.not. agrees .and. proc0) print*,"Failure : ",name," : Got ",integral," but expected ",expected_value
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

  logical function integrate_pitch_angle_test() result(agrees)
    implicit none
    ! Integral[alpha*sin[alpha],{alpha,0,pi}]/Integral[sin[alpha],{alpha,0,pi}]
    real, parameter :: expected_value = pi / 2.0
    character(len=*), parameter :: name = 'pitch angle'
    real :: pitch_angle
    gnew = 1.0d0
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

end program test_le_grids_integrate
