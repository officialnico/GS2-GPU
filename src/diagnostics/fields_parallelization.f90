#include "unused_dummy.inc"
!> FIXME : Add documentation
module fields_parallelization
   implicit none
   private
   public :: field_k_local
contains
  !> FIXME : Add documentation
  !>
  !> @note It's not clear that this function is correct in general as it is
  !> assigning {it,ik} field ownership in a round-robin sense but this might
  !> not align with fields_local (i.e. g_lo) ownership)
  elemental logical function field_k_local(it, ik)
    use mp, only: iproc, nproc
    implicit none
    integer, intent(in) :: ik, it
    UNUSED_DUMMY(ik)
    ! This is temporary while the fields are being parallelised
    field_k_local = (mod(it-1, nproc) == iproc)
  end function field_k_local

end module fields_parallelization
