!> Provides storage for the EM potentials
!! and some related data
module fields_arrays
  use constants, only: run_name_size
  implicit none

  private

  public :: phi, apar, bpar, phinew, aparnew, bparnew
  public :: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
  public :: apar_ext, time_field, response_file
  public :: get_specific_response_file_name

  !Main fields
  complex, dimension (:,:,:), allocatable :: phi,    apar,    bpar
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
  complex, dimension (:,:,:), allocatable :: gf_phi, gf_apar, gf_bpar
  complex, dimension (:,:,:), allocatable :: gf_phinew, gf_aparnew, gf_bparnew
  !For antenna
  complex, dimension (:,:,:), allocatable :: apar_ext
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  !Timing data
  real :: time_field(2)=0.
  
  !For response data
  character(run_name_size) :: response_file = 'NOT_SET'

contains

  !> This function returns the response file name when given
  !> characteristic data.
  function get_specific_response_file_name(ik, is, suffix) result(file_name)
    implicit none
    integer, intent(in) :: ik, is
    character(len = *), optional, intent(in) :: suffix
    character(len = 256) :: file_name

    character(len = 64) :: suffix_local
    character(len = 64), parameter :: suffix_default='.response'

    !Set file suffix
    suffix_local=suffix_default
    if(present(suffix)) suffix_local=suffix

    write(file_name,'(A,"_ik_",I0,"_is_",I0,A)') trim(response_file),ik,is,trim(suffix_local)
  end function get_specific_response_file_name

end module fields_arrays
