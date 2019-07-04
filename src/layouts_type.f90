!> FIXME : Add documentation
module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90
  use mp, only: comm_type
  
  implicit none

  private
  public :: ikit, ikitprocs
  public :: g_layout_type, lz_layout_type, e_layout_type
  public :: le_layout_type, p_layout_type, gf_layout_type
  public :: deallocate_layout

  !> FIXME : Add documentation  
  type :: ikit
     sequence
     integer :: ik, it
  end type ikit

  !> FIXME : Add documentation  
  type :: ikitprocs
     sequence
     logical :: mine
     integer :: num_procs
     type(comm_type) :: comm
     integer,dimension(:),allocatable :: proc_list
     integer,dimension(:),allocatable :: sub_proc_list
  end type ikitprocs

  !> FIXME : Add documentation  
  type :: g_layout_type
     sequence
     integer :: iproc
     integer :: naky, ntheta0, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ik_ord, it_ord, il_ord, ie_ord, is_ord
     integer :: ik_comp, it_comp, il_comp, ie_comp, is_comp
     integer,dimension(5) :: compound_count, dim_size
!    For replacing loops over iglo with 5 separate loops
!    This is not currently implemented, but involves replacing "do iglo=...." with
!    EQUIVALENCE(dind(1),d1),(dind(2),d2).... !Note this is just a short hand for storing currently d values in an array   
!    ik=>dind(g_lo%ik_ord) !This is a pointer assignment, so that whatever layout we use ik points to the correct d value
!    do d5=d5_min,d5max ; do d4=d4_min,d4_max ...
     integer :: d1_min,d1_max,d2_min,d2_max,d3_min,d3_max
     integer :: d4_min,d4_max,d5_min,d5_max
!    Used in deciding colour for sub-communicator creation.
!    Also used to select slices of array to take part in sub-comm global reductions 
!    (in order to reduce message size by avoiding regions of arrays outside this procs domain)
     integer :: ik_min, ik_max, it_min, it_max, il_min, il_max, ie_min, ie_max, is_min, is_max
     integer :: xyblock_comm !Sub comms for xy blocks (i.e. l-e-s integrals)
     integer :: xysblock_comm !Sub comms for xys blocks (i.e. l-e integrals)
     integer :: lesblock_comm !Sub comms for les blocks (i.e. x-y integrals)
     integer,dimension(:,:),allocatable :: les_kxky_range !Used in integrate species, holds start and stop indices for flattened kxky dimension
     type(ikitprocs),dimension(:,:),allocatable :: ikit_procs_assignment !Used to record processes that have been assigned specific it ik points in g
     integer :: ikitrange
     integer :: max_ikit_nprocs ! Used to record the largest proc_list size for all the ikitprocs points.  This is useful for setting the size of non-blocking request arrays for communications
     type(ikit),dimension(:),allocatable :: local_ikit_points
     logical :: x_before_y !Information about layout order
     logical :: x_local, y_local, l_local, e_local, s_local !Used to record if various dimensions are entirely local
! shm data
     integer :: ppn, nnd, bpn, ndbklsize, llim_node, ulim_node
     integer, allocatable :: proc_map(:) ! llimit on all ranks
  end type g_layout_type

  !> FIXME : Add documentation  
  type :: lz_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, isgn_ord, ik_ord, it_ord, ie_ord, is_ord
     integer :: ig_comp, isgn_comp, ik_comp, it_comp, ie_comp, is_comp
     integer,dimension(5) :: compound_count
     !> Array that holds the size of each dimension in the order that they
     !> appear in the layout. This allows the removal of slow
     !> 'select case(layout)' from indexing routines.
     integer,dimension(5) :: dim_size

  end type lz_layout_type

  !> FIXME : Add documentation
  type :: e_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nlambda, nspec, nsign, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, isgn_ord, ik_ord, it_ord, il_ord, is_ord
     integer :: ig_comp, isgn_comp, ik_comp, it_comp, il_comp, is_comp
     integer,dimension(6) :: compound_count
     !> Array that holds the size of each dimension in the order that they
     !> appear in the layout. This allows the removal of slow
     !> 'select case(layout)' from indexing routines.
     integer,dimension(6) :: dim_size
  end type e_layout_type

  !> FIXME : Add documentation
  type :: le_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, ik_ord, it_ord, is_ord
     integer :: ig_comp, ik_comp, it_comp, is_comp
     integer,dimension(4) :: compound_count
     integer :: ik_min, it_min, ig_min, is_min 
     integer :: ik_max, it_max, ig_max, is_max 
     logical :: x_local, y_local, t_local, s_local !Used to record if various dimensions are entirely local
     !> Array that holds the size of each dimension in the order that they
     !> appear in the layout. This allows the removal of slow
     !> 'select case(layout)' from indexing routines.
     integer,dimension(4) :: dim_size
  end type le_layout_type

  !> FIXME : Add documentation  
  type :: gf_layout_type
     sequence
     integer :: iproc, nproc
     integer :: ntgrid, naky, ntheta0, nspec, ntgridtotal, npsec, nlambda, negrid
     integer :: xypoints
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer :: smallblocksize, largeblocksize, largeblocklimit, divisionblock
     integer :: largeregionlimit, smallgapsize, largegapsize
  end type gf_layout_type

  !> FIXME : Add documentation
  type :: p_layout_type
     sequence
     integer :: iproc
     integer :: naky, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ik_ord, il_ord, ie_ord, is_ord
     integer :: ik_comp, il_comp, ie_comp, is_comp
     integer,dimension(4) :: compound_count
  end type p_layout_type

  interface deallocate_layout
     module procedure deallocate_g_layout
  end interface

  contains
    !> Deallocate g_layout's arrays
    subroutine deallocate_g_layout(g_lo)
      implicit none
      type(g_layout_type), intent(inout) :: g_lo

      if (allocated(g_lo%les_kxky_range)) deallocate(g_lo%les_kxky_range)
      ! Note, we do not deallocate ikit_procs_assignment()%proc_list or ikit_procs_assignment()%sub_proc_list as in Fortran these 
      ! should automatically be deallocated when ikit_procs_assignment is deallocated.
      if(allocated(g_lo%ikit_procs_assignment)) deallocate(g_lo%ikit_procs_assignment)
      if(allocated(g_lo%local_ikit_points)) deallocate(g_lo%local_ikit_points)

      if(allocated(g_lo%proc_map)) deallocate(g_lo%proc_map)

    end subroutine deallocate_g_layout

end module layouts_type
