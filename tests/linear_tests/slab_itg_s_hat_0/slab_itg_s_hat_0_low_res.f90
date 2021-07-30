!> A test program to run a linear slab itg benchmark
!! at low resolution and check the growth rate

module checks_mod
  use functional_tests
  public checks
  contains
    function checks()
      logical :: checks
      ! This test case is relatively low res so whilst we test to 0.1%
      ! we do anticipate the golden answer varying as code numerics are
      ! updated/changed etc.
      checks =  check_growth_rate((/0.1618/), 0.001)
    end function checks
end module checks_mod

program slab_itg_s_hat_0_low_res
  use functional_tests
  use checks_mod


  call test_gs2('Linear slab ITG with zero magnetic shear', checks)


end program slab_itg_s_hat_0_low_res
