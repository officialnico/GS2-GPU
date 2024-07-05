!> FIXME : Add documentation
module fields_test
  implicit none

  private

  public :: init_fields_test
  public :: advance_test
  public :: init_phi_test
  public :: reset_init

  logical :: initialized = .false.

contains
  !> FIXME : Add documentation
  subroutine init_fields_test
    implicit none
    if (initialized) return
    initialized = .true.
  end subroutine init_fields_test

  !> FIXME : Add documentation  
  subroutine init_phi_test
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use array_utils, only: zero_array
    implicit none
    call zero_array(phi) ; call zero_array(phinew)
    call zero_array(apar) ; call zero_array(aparnew)
    call zero_array(bpar) ; call zero_array(bparnew)
  end subroutine init_phi_test

  !> FIXME : Add documentation
  subroutine advance_test (istep)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: timeadv
    use dist_fn_arrays, only: g, gnew
    use array_utils, only: copy
    implicit none
    integer, intent (in) :: istep
    call copy(gnew, g)
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep)
  end subroutine advance_test

  !> FIXME : Add documentation
  subroutine reset_init
    initialized = .false.
  end subroutine reset_init

end module fields_test
