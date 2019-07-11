! Modifications for unbalanced decomposition functionality:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project

!> FIXME : Add documentation
module gs2_layouts
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  use layouts_type, only: ikit, ikitprocs
  use layouts_type, only: g_layout_type, lz_layout_type, e_layout_type
  use layouts_type, only: le_layout_type, p_layout_type, gf_layout_type
  use constants, only: run_name_size

  implicit none

  private

  public :: factors
  public :: pe_layout, layout, local_field_solve
  
  public :: init_dist_fn_layouts, init_gs2_layouts
  public :: finish_dist_fn_layouts
  public :: finish_gs2_layouts
  public :: wnml_gs2_layouts
  public :: init_parity_layouts
  public :: g_lo
  public :: p_lo

  public :: finish_layouts

  public :: set_overrides

  public :: finish_fields_layouts, finish_jfields_layouts

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

  public :: init_lambda_layouts
  public :: lz_lo

  public :: init_energy_layouts
  public :: e_lo

  public :: init_le_layouts
  public :: le_lo

  public :: init_gf_layouts
  public :: gf_lo

  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: calculate_unbalanced_x, calculate_unbalanced_y, calculate_idle_processes
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: gidx2xxfidx, xxfidx2yxfidx

  public :: init_accel_transform_layouts
  public :: accel_lo, accel_layout_type, dealiasing
  public :: accelx_lo, accelx_layout_type

  public :: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx, if_idx, isign_idx
  public :: im_idx, in_idx, ij_idx, ifield_idx
  public :: idx, proc_id, idx_local

  public :: opt_local_copy
  public :: opt_redist_init !Do we use optimised redist init routines?
  public :: intmom_sub !Do we use sub-communicators in velocity space?
  public :: intspec_sub !Do we use sub-communicators in integrate_species?

  public :: fft_use_wisdom, fft_wisdom_file, fft_measure_plan

  public :: simple_gf_decomposition, gf_local_fields

  public :: layouts_config_type
  public :: set_gs2_layouts_config

  logical :: initialized = .false.
  logical :: initialized_x_transform = .false.
  logical :: initialized_y_transform = .false.
  logical :: initialized_layouts = .false.
  logical :: initialized_dist_fn_layouts = .false.
  logical :: initialized_fields_layouts = .false.
  logical :: initialized_jfields_layouts = .false.
  logical :: initialized_le_layouts = .false.
  logical :: initialized_gf_layouts = .false.
  logical :: initialized_energy_layouts = .false.
  logical :: initialized_lambda_layouts = .false.
  logical :: initialized_parity_layouts = .false.
  logical :: initialized_accel_transform_layouts = .false.

  logical :: opt_local_copy
  logical :: opt_redist_init
  logical :: intmom_sub
  logical :: intspec_sub
  logical :: local_field_solve, accel_lxyes, lambda_local, unbalanced_xxf, unbalanced_yxf
  real :: max_unbalanced_xxf, max_unbalanced_yxf
  character (len=5) :: layout
  character(run_name_size) :: fft_wisdom_file
  logical :: fft_use_wisdom, fft_measure_plan
  logical :: exist
  logical :: simple_gf_decomposition
  logical :: gf_local_fields

  type (g_layout_type) :: g_lo
  type (le_layout_type) :: le_lo
  !> layout for fields_gf_local and gf_lo integration
  type (gf_layout_type) :: gf_lo
  type (p_layout_type) :: p_lo

  !> FIXME : Add documentation  
  type :: f_layout_type
     integer :: iproc
     integer :: nidx, nfield, ntgrid, nindex, naky, ntheta0, M_class, N_class, i_class
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer, dimension(:,:), pointer :: ik => null ()
     integer, dimension(:,:), pointer :: it => null ()
  end type f_layout_type

  type (f_layout_type), dimension(:), allocatable :: f_lo

  !> FIXME : Add documentation
  type :: jf_layout_type
     integer :: iproc
     integer :: nindex, naky, ntheta0, ntgrid, nfield
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer :: ig_min,ig_max,ifield_min,ifield_max,it_min,it_max,ik_min,ik_max
  end type jf_layout_type

  type (jf_layout_type) :: jf_lo

  integer, dimension(:), allocatable :: ij, mj
  integer, dimension(:,:), allocatable :: dj

  type (lz_layout_type) :: lz_lo

  type (e_layout_type) :: e_lo

  !> FIXME : Add documentation
  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ntheta0, nx, nadd, negrid, nlambda, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     !> Order of iX in the distributed index. e.g. ig_ord = 1 => theta is the fastest moving distributed index
     integer :: ig_ord, isgn_ord, ik_ord, il_ord, ie_ord, is_ord
     !> `it_comp` is the product of the sizes of dimensions which are more
     !> local than kx in the g_lo layout.
     integer :: ig_comp, isgn_comp, ik_comp, il_comp, ie_comp, is_comp
     !> Work array that holds iX_comp in the order that they appear in the
     !> layout.
     integer,dimension(6) :: compound_count
     !> Array that holds the size of each dimension in the order that they
     !> appear in the layout.
     integer,dimension(6) :: dim_size
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  !> FIXME : Add documentation  
  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, negrid, nlambda, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     !> Order of iX in the distributed index. e.g. ig_ord = 1 => theta is the fastest moving distributed index
     integer :: ig_ord, isgn_ord, it_ord, il_ord, ie_ord, is_ord
     !> `it_comp` is the product of the sizes of dimensions which are more
     !> local than kx in the g_lo layout.
     integer :: ig_comp, isgn_comp, it_comp, il_comp, ie_comp, is_comp
     !> Work array that holds iX_comp in the order that they appear in the
     !> layout.
     integer,dimension(6) :: compound_count
     !> Array that holds the size of each dimension in the order that they
     !> appear in the layout.
     integer,dimension(6) :: dim_size
  end type yxf_layout_type

  type (yxf_layout_type) :: yxf_lo

  !> FIXME : Add documentation  
  type :: accel_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ndky, ny, ntheta0
     integer :: nx, nxny, nxnky, negrid, nlambda, nspec, nia
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
#ifdef SHMEM
     integer :: ppn, nnd, bpn, llim_node, ulim_node !shm data
#endif
  end type accel_layout_type

  type (accel_layout_type) :: accel_lo

  !> FIXME : Add documentation  
  type :: accelx_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, nxny, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
#ifdef SHMEM
     integer :: ppn, nnd, bpn, llim_node, ulim_node ! shm data
     integer, allocatable :: proc_map(:)
#endif
  end type accelx_layout_type

  type (accelx_layout_type) :: accelx_lo
  
  !> Used to represent the input configuration of layouts
  type, extends(abstract_config_type) :: layouts_config_type
     ! namelist : layouts_knobs
     ! indexed : false     
     !> False=>use heuristics rather than measure in fft plan create
     logical :: fft_measure_plan = .false.
     !> Try to load and save wisdom about fftw plans to fft_wisdom_file
     logical :: fft_use_wisdom = .true. 
     !> Location of fftw wisdom... if left as default, this is set to
     !> run_name//'.fftw_wisdom', unless overriden by the environment
     !> variable GK_FFTW3_WISDOM. If set to anything other default,
     !> overrides GK_FFTW3_WISDOM
     character(len = run_name_size) :: fft_wisdom_file = 'default'
     !> FIXME: Add documentation
     logical :: gf_local_fields = .false.
     !> FIXME: Add documentation
     logical :: intmom_sub = .false.
     !> FIXME: Add documentation
     logical :: intspec_sub = .false.
     !> FIXME: Add documentation
     character(len = 5) :: layout = 'lxyes'
     !> FIXME: Add documentation
     logical :: local_field_solve = .false.
     !> FIXME: Add documentation
     real :: max_unbalanced_xxf = 0.0
     !> FIXME: Add documentation
     logical :: opt_local_copy = .false. 
     !> True=>Use optimised routines to init redist objects
     logical :: opt_redist_init = .false.
     !> True=>Use nonblocking redistributes
     logical :: opt_redist_nbk = .true.
     !> True=>Use persistent communications in redistributes
     logical :: opt_redist_persist = .false.
     !> True=>Start comms before doing local copy
     logical :: opt_redist_persist_overlap = .false.
     !> FIXME: Add documentation
     logical :: simple_gf_decomposition = .true.
     !> FIXME: Add documentation
     logical :: unbalanced_xxf = .false.
     !> FIXME: Add documentation
     logical :: unbalanced_yxf = .false.
   contains
     procedure, public :: read => read_layouts_config
     procedure, public :: write => write_layouts_config
     procedure, public :: reset => reset_layouts_config
     procedure, public :: broadcast => broadcast_layouts_config
     procedure, public, nopass :: get_default_name => get_default_name_layouts_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_layouts_config     
  end type layouts_config_type
  
  type(layouts_config_type) :: layouts_config
  
  interface ij_idx
     module procedure ij_idx_f
     module procedure ij_idx_jf
  end interface

  interface im_idx
     module procedure im_idx_f
  end interface

  interface in_idx
     module procedure in_idx_f
  end interface

  interface ifield_idx
     module procedure ifield_idx_f
     module procedure ifield_idx_jf
  end interface

  interface if_idx
     module procedure if_idx_f
     module procedure if_idx_jf
  end interface

  interface ig_idx
     module procedure ig_idx_lz
     module procedure ig_idx_e
     module procedure ig_idx_le
     module procedure ig_idx_xxf
     module procedure ig_idx_yxf
     module procedure ig_idx_f
     module procedure ig_idx_jf
  end interface

  interface ik_idx
     module procedure ik_idx_g
     module procedure ik_idx_jf
     module procedure ik_idx_lz
     module procedure ik_idx_e
     module procedure ik_idx_le
     module procedure ik_idx_gf
     module procedure ik_idx_xxf
     module procedure ik_idx_accel
     module procedure ik_idx_accelx
     module procedure ik_idx_parity
     module procedure ik_idx_f
  end interface

  interface it_idx
     module procedure it_idx_g
     module procedure it_idx_jf
     module procedure it_idx_lz
     module procedure it_idx_e
     module procedure it_idx_le
     module procedure it_idx_gf
     module procedure it_idx_yxf
     module procedure it_idx_accel
     module procedure it_idx_accelx
     module procedure it_idx_f
  end interface

  interface il_idx
     module procedure il_idx_g
     module procedure il_idx_e
     module procedure il_idx_xxf
     module procedure il_idx_yxf
     module procedure il_idx_accelx
     module procedure il_idx_parity
  end interface

  interface ie_idx
     module procedure ie_idx_g
     module procedure ie_idx_lz
     module procedure ie_idx_xxf
     module procedure ie_idx_yxf
     module procedure ie_idx_accelx
     module procedure ie_idx_parity
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_lz
     module procedure is_idx_e
     module procedure is_idx_le
     module procedure is_idx_xxf
     module procedure is_idx_yxf
     module procedure is_idx_accelx
     module procedure is_idx_parity
  end interface

  interface isign_idx
     module procedure isign_idx_xxf
     module procedure isign_idx_yxf
     module procedure isign_idx_e
  end interface

  interface proc_id
     module procedure proc_id_g
     module procedure proc_id_f
     module procedure proc_id_jf
     module procedure proc_id_lz
     module procedure proc_id_e
     module procedure proc_id_le
     module procedure proc_id_gf
     module procedure proc_id_xxf
     module procedure proc_id_yxf
     module procedure proc_id_accelx
     module procedure proc_id_parity
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_f
     module procedure idx_jf
     module procedure idx_lz
     module procedure idx_e
     module procedure idx_le
     module procedure idx_gf
     module procedure idx_xxf
     module procedure idx_yxf
     module procedure idx_accelx
     module procedure idx_parity
  end interface

  interface idx_local
     module procedure idx_local_g,      ig_local_g
     module procedure idx_local_f,      ig_local_f
     module procedure idx_local_jf,     ig_local_jf
     module procedure idx_local_lz,     ig_local_lz
     module procedure idx_local_e,      ig_local_e
     module procedure idx_local_le,      ig_local_le
     module procedure idx_local_gf,      ig_local_gf
     module procedure idx_local_xxf,    ig_local_xxf
     module procedure idx_local_yxf,    ig_local_yxf
     module procedure idx_local_accelx, ig_local_accelx
     module procedure idx_local_parity, ig_local_parity
  end interface

  interface set_dimension_order
     module procedure set_dimension_order_g
     module procedure set_dimension_order_le
     module procedure set_dimension_order_lz
     module procedure set_dimension_order_e
     module procedure set_dimension_order_xxf
     module procedure set_dimension_order_yxf
  end interface set_dimension_order
  
contains 
  !> FIXME : Add documentation 
  subroutine wnml_gs2_layouts(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "layouts_knobs"
    write (unit, fmt="(' layout = ',a)") '"'//trim(layout)//'"'
    write (unit, fmt="(' local_field_solve = ',L1)") local_field_solve
    write (unit, fmt="(' /')")
  end subroutine wnml_gs2_layouts

  !> FIXME : Add documentation   
  subroutine init_gs2_layouts(layouts_config_in)
    use mp, only: proc0
    implicit none
    type(layouts_config_type), intent(in), optional :: layouts_config_in
    
    if (initialized_layouts) return
    initialized_layouts = .true.

    if (initialized) return
    initialized = .true.
    call read_parameters(layouts_config_in)
  end subroutine init_gs2_layouts

  !> FIXME : Add documentation 
  subroutine finish_gs2_layouts
    call finish_layouts
    initialized = .false.
    initialized_layouts = .false.
  end subroutine finish_gs2_layouts

  !> FIXME : Add documentation   
  subroutine read_parameters(layouts_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist, error_unit
    use redistribute, only: opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap
    use mp, only: mp_abort, proc0, broadcast
    implicit none
    type(layouts_config_type), intent(in), optional :: layouts_config_in    
    integer :: in_file

    if (present(layouts_config_in)) layouts_config = layouts_config_in

    call layouts_config%init(name = 'layouts_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    fft_measure_plan = layouts_config%fft_measure_plan
    fft_use_wisdom = layouts_config%fft_use_wisdom
    fft_wisdom_file = layouts_config%fft_wisdom_file
    gf_local_fields = layouts_config%gf_local_fields
    intmom_sub = layouts_config%intmom_sub
    intspec_sub = layouts_config%intspec_sub
    layout = layouts_config%layout
    local_field_solve = layouts_config%local_field_solve
    max_unbalanced_xxf = layouts_config%max_unbalanced_xxf
    opt_local_copy = layouts_config%opt_local_copy
    opt_redist_init = layouts_config%opt_redist_init
    opt_redist_nbk = layouts_config%opt_redist_nbk
    opt_redist_persist = layouts_config%opt_redist_persist
    opt_redist_persist_overlap = layouts_config%opt_redist_persist_overlap
    simple_gf_decomposition = layouts_config%simple_gf_decomposition
    unbalanced_xxf = layouts_config%unbalanced_xxf
    unbalanced_yxf = layouts_config%unbalanced_yxf

    exist = layouts_config%exist
    
    if (layout.ne.'yxels' .and. layout.ne.'yxles' .and. layout.ne.'lexys'&
    .and. layout.ne.'lxyes' .and. layout.ne.'lyxes' .and. layout.ne.'xyles') &
    then
       write(6,*) "gs2_layouts: read_parameters finds illegal layout=",layout," =>stop"
       call mp_abort("gs2_layouts: read_parameters finds illegal layout=")
    endif

    ! max_unbalanced_xxf and max_unbalanced_yxf have a maximum range of 0.0 - 1.0
    ! so check here whether the user has specified the range correctly and adjust if 
    ! they have not.
    if(max_unbalanced_xxf .gt. 1.0) then
       max_unbalanced_xxf = 1.0
       write(*,*) 'max_unbalanced_xxf too large, setting to 1.0'
    else if(max_unbalanced_xxf .lt. 0.0) then
       max_unbalanced_xxf = 0.0
       write(*,*) 'max_unbalanced_xxf too small, setting to 0.0'
    end if
    if(max_unbalanced_yxf .gt. 1.0) then
       max_unbalanced_yxf = 1.0
       write(*,*) 'max_unbalanced_yxf too large, setting to 1.0'
    else if(max_unbalanced_yxf .lt. 0.0) then
       max_unbalanced_yxf = 0.0
       write(*,*) 'max_unbalanced_yxf too small, setting to 0.0'
    end if

    call ensure_redist_persist_consistency(opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap)

    ! For now only proc0 does this as we can't promise that all procs know the run_name
    if(proc0) call get_wisdom_file(fft_wisdom_file)
    call broadcast(fft_wisdom_file)

  end subroutine read_parameters

  !> FIXME : Add documentation 
  subroutine set_overrides(opt_ov)
    use overrides, only: optimisations_overrides_type
    use redistribute, only: opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap
    type(optimisations_overrides_type), intent(in) :: opt_ov
    if (opt_ov%override_layout) layout = opt_ov%layout
    if (opt_ov%override_opt_redist_nbk) &
      opt_redist_nbk = opt_ov%opt_redist_nbk
    if (opt_ov%override_opt_redist_persist) &
      opt_redist_persist = opt_ov%opt_redist_persist
    if (opt_ov%override_intmom_sub) &
      intmom_sub = opt_ov%intmom_sub
    if (opt_ov%override_intspec_sub) &
      intspec_sub = opt_ov%intspec_sub

    call ensure_redist_persist_consistency(opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap)

  end subroutine set_overrides

  !> Disable opt_redist_persist and opt_redist_persist_overlap if dependent settings not set
  subroutine ensure_redist_persist_consistency(opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap)

    implicit none
    logical, intent(in) :: opt_redist_nbk
    logical, intent(inout) :: opt_redist_persist
    logical, intent(inout) :: opt_redist_persist_overlap

    opt_redist_persist = opt_redist_persist .and. opt_redist_nbk
    opt_redist_persist_overlap = opt_redist_persist_overlap .and. opt_redist_persist
    
  end subroutine ensure_redist_persist_consistency

  !> Set the fft_wisdom_file based on the value set in the namelist (if any)
  !! and the value in the GK_FFTW3_WISDOM environment variable
  !!
  !! @note This doesn't seem like the correct place for this routine. Can we move into
  !! the fft module?
  subroutine get_wisdom_file(wisdom_file)
    use file_utils, only: run_name
    use constants, only: run_name_size
    character(len=*), intent(inout) :: wisdom_file
    character(run_name_size) :: env_wisdom_file
    call get_environment_variable("GK_FFTW3_WISDOM", env_wisdom_file)
    !read (env_wisdom_file,'(I10)') verbosity

    if (wisdom_file .eq. 'default') then
      if (trim(env_wisdom_file) .eq. '') then
        wisdom_file = trim(run_name)//'.fftw_wisdom'
      else
        wisdom_file = env_wisdom_file
      end if
    end if
  end subroutine get_wisdom_file

  !> FIXME : Add documentation   
  subroutine check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)
    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    integer, intent (out) :: nblock
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs
    integer :: i

    if (.not. layout == 'lxyes') then
       accel_lxyes = .false.
       return
    end if

    accel_lxyes = .true.

    allocate (facs(max(nspec,negrid,naky,ntheta0,nlambda)/2+1,5))
    call factors (nspec,   nsfacs, facs(:,1))
    call factors (negrid,  nefacs, facs(:,2))
    call factors (naky,    nyfacs, facs(:,3))
    call factors (ntheta0, nxfacs, facs(:,4))
    call factors (nlambda, nlfacs, facs(:,5))
    
!    do i=1,nsfacs-1
!       if (nproc == facs(i,1)) then
!          nblock = facs(i,1)
!          goto 100
!       end if
!    end do
!    
!    do i=1,nefacs-1
!       if (nproc == facs(i,2)*nspec) then
!          nblock = facs(i,2)*nspec
!          goto 100
!       end if
!    end do
    
    nblock = negrid*nspec
    do i=1,nyfacs-1
       if (nproc == facs(i,3)*negrid*nspec) goto 100
    end do
    
    do i=1,nxfacs-1
       if (nproc == facs(i,4)*naky*negrid*nspec) goto 100
    end do
    
    do i=1,nlfacs-1
       if (nproc == facs(i,5)*ntheta0*naky*negrid*nspec) goto 100
    end do
    
    accel_lxyes = .false.

100 deallocate (facs)

  end subroutine check_accel

  !> FIXME : Add documentation 
  subroutine check_llocal (ntheta0, naky, nlambda, negrid, nspec)
    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs
    integer :: i

    lambda_local = .true.

    allocate (facs(max(nspec,negrid,naky,ntheta0,nlambda)/2+1,5))

    select case (layout)

    case ('lxyes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (naky,    nyfacs, facs(:,3))
       call factors (ntheta0, nxfacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nyfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
       do i=1,nxfacs-1
          if (nproc == facs(i,4)*naky*negrid*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*ntheta0*naky*negrid*nspec) goto 100
       end do
          
    case ('lyxes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (ntheta0, nxfacs, facs(:,3))
       call factors (naky,    nyfacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nxfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
       do i=1,nyfacs-1
          if (nproc == facs(i,4)*ntheta0*negrid*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*naky*ntheta0*negrid*nspec) goto 100
       end do
          
    case ('lexys')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (naky,    nyfacs, facs(:,2))
       call factors (ntheta0, nxfacs, facs(:,3))
       call factors (negrid,  nefacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nyfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nxfacs-1
          if (nproc == facs(i,3)*naky*nspec) goto 100
       end do
          
       do i=1,nefacs-1
          if (nproc == facs(i,4)*ntheta0*naky*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*negrid*ntheta0*naky*nspec) goto 100
       end do

    end select

    lambda_local = .false.    
                       
100 deallocate (facs)


  end subroutine check_llocal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation 
  subroutine init_dist_fn_layouts &
       (naky, ntheta0, nlambda, negrid, nspec)
    use mp, only: iproc, nproc, proc0, comm_type, mp_abort
    use mp, only: split, free_comm, mp_comm, nproc_comm, rank_comm
    use mp, only: sum_allreduce_sub, sum_allreduce, sum_reduce, max_reduce
    use file_utils, only: error_unit
#ifdef SHMEM
    use shm_mpi3, only: shm_info
    use mp, only: allgather
#endif
    implicit none
    integer, intent (in) :: naky, ntheta0, nlambda, negrid, nspec
    integer :: iglo,ik,it,il,ie,is,col,mycol
    integer :: ik_min,ik_max,it_min,it_max,il_min,il_max,ie_min,ie_max,is_min,is_max,ip
    integer :: nproc_subcomm, rank_subcomm
    integer :: nproc_tmp, idim, tmp, i
    real :: tmp_r
#ifdef SHMEM
    integer nxy, les, xyb, nd_les, llim_nd, left_lim_nd, ind, sfft
#endif
    integer,dimension(5) :: nproc_dim
    integer,dimension(:),allocatable :: les_kxky_range
    logical,dimension(5) :: dim_divides,dim_local
    logical,dimension(:),allocatable :: it_local, ik_local, il_local, ie_local, is_local
    integer,dimension(nproc) :: tmp_proc_list
    integer :: colour
    type(comm_type) :: tmpcomm

    if (initialized_dist_fn_layouts) return
    initialized_dist_fn_layouts = .true.
   
    g_lo%iproc = iproc
    g_lo%naky = naky
    g_lo%ntheta0 = ntheta0
    g_lo%nlambda = nlambda
    g_lo%negrid = negrid
    g_lo%nspec = nspec
    g_lo%llim_world = 0
    g_lo%ulim_world = naky*ntheta0*negrid*nlambda*nspec - 1
      
#ifndef SHMEM
    g_lo%blocksize = g_lo%ulim_world/nproc + 1
    g_lo%llim_proc = g_lo%blocksize*iproc
    g_lo%ulim_proc = min(g_lo%ulim_world, g_lo%llim_proc + g_lo%blocksize - 1)
    g_lo%ulim_alloc = max(g_lo%llim_proc, g_lo%ulim_proc)

#else

    nxy=naky*ntheta0
    les = nlambda*negrid*nspec

    !shm 
    g_lo%ppn = shm_info%size
    g_lo%nnd = nproc / g_lo%ppn
    ! needs a test and to see if every node has the same number of processors
    !if ( mod(nproc, g_lo%ppn) /= 0) g_lo%nnd =g_lo%nnd + 1
    ind = iproc / g_lo%ppn
    g_lo%bpn = les / g_lo%nnd
    if (ind < mod(les,  g_lo%nnd )) g_lo%bpn = g_lo%bpn + 1

    llim_nd = ind * (les / g_lo%nnd) + min(ind, mod(les, g_lo%nnd))
    g_lo%llim_node = llim_nd * nxy
    g_lo%ulim_node =  g_lo%llim_node + g_lo%bpn * nxy -1

    nd_les = g_lo%bpn
    xyb = (nxy * nd_les) / g_lo%ppn
    g_lo%blocksize =  xyb
    if (shm_info%id < mod(nxy * nd_les, g_lo%ppn)) g_lo%blocksize = g_lo%blocksize + 1
    !g_lo%llim_proc = g_lo%blocksize*iproc
    
    g_lo%llim_proc = g_lo%llim_node &
         + shm_info%id * xyb + min(shm_info%id, mod(nxy * nd_les, g_lo%ppn))
    g_lo%ulim_proc = min(g_lo%ulim_world, g_lo%llim_proc + g_lo%blocksize - 1)
    g_lo%ulim_alloc = max(g_lo%llim_proc, g_lo%ulim_proc)

    allocate(g_lo%proc_map(0:nproc-1))
    call allgather((/ g_lo%llim_proc /), 1, g_lo%proc_map, 1)  

    ! report on shm details
    if (proc0) then 
       write(*,'(/,a)') "g_lo shared memory parameters (rank 0)" 
       write(*,*) "ranks per node :", g_lo%ppn
       write(*,*) "total number of blocks", les, "; and per node :", g_lo%bpn
       write(*,*) "block size inside node :", g_lo%blocksize
    endif
    
#endif

    ! Calculate constants used in the index lookup routines
    ! (rather than calculating them at each call as done previously)
    ! This allows the removal of the select case statements in the index routines
    select case (layout)
    case ('yxels')
       call set_dimension_order(g_lo,ik_ord=1,it_ord=2,ie_ord=3,il_ord=4,is_ord=5)
    case ('yxles')
       call set_dimension_order(g_lo,ik_ord=1,it_ord=2,il_ord=3,ie_ord=4,is_ord=5)
    case ('lexys')
       call set_dimension_order(g_lo,il_ord=1,ie_ord=2,it_ord=3,ik_ord=4,is_ord=5)
    case ('lxyes')
       call set_dimension_order(g_lo,il_ord=1,it_ord=2,ik_ord=3,ie_ord=4,is_ord=5)
    case ('lyxes')
       call set_dimension_order(g_lo,il_ord=1,ik_ord=2,it_ord=3,ie_ord=4,is_ord=5)
    case ('xyles')
       call set_dimension_order(g_lo,it_ord=1,ik_ord=2,il_ord=3,ie_ord=4,is_ord=5)
    end select

    g_lo%dim_size(g_lo%it_ord)=ntheta0
    g_lo%dim_size(g_lo%ik_ord)=naky
    g_lo%dim_size(g_lo%il_ord)=nlambda
    g_lo%dim_size(g_lo%ie_ord)=negrid
    g_lo%dim_size(g_lo%is_ord)=nspec
    
    g_lo%compound_count(1)=1
    g_lo%compound_count(2)=g_lo%compound_count(1)*g_lo%dim_size(1)
    g_lo%compound_count(3)=g_lo%compound_count(2)*g_lo%dim_size(2)
    g_lo%compound_count(4)=g_lo%compound_count(3)*g_lo%dim_size(3)
    g_lo%compound_count(5)=g_lo%compound_count(4)*g_lo%dim_size(4)

    g_lo%ik_comp=g_lo%compound_count(g_lo%ik_ord)
    g_lo%it_comp=g_lo%compound_count(g_lo%it_ord)
    g_lo%il_comp=g_lo%compound_count(g_lo%il_ord)
    g_lo%ie_comp=g_lo%compound_count(g_lo%ie_ord)
    g_lo%is_comp=g_lo%compound_count(g_lo%is_ord)

    ! Work out min and max of the five compound indices and the locality
    ! These are currently only used to slice the arrays used in the velocity space integrals
    ! but could well be useful in a number of places.
    ! NOTE: These are typically only useful if the dimensions are split nicely (i.e. we're using
    ! a sweet spot) otherwise our local ik section may look like, [naky-1,naky,1,2], so whilst the 
    ! min and max values (1,naky) are technically correct, we cannot assume we cover the range 1-naky
    ik_max=0
    ik_min=naky+1
    it_max=0
    it_min=ntheta0+1
    il_max=0
    il_min=nlambda+1
    ie_max=0
    ie_min=negrid+1
    is_max=0
    is_min=nspec+1
    allocate(ik_local(naky),it_local(ntheta0),il_local(nlambda),ie_local(negrid),is_local(nspec))
    ik_local=.false.
    it_local=.false.
    il_local=.false.
    ie_local=.false.
    is_local=.false.
    g_lo%ikitrange = 0
    ! g_lo%ikit_procs_assignment is used to record which ik,it points I have from the range of 
    ! ik,it points.
    allocate(g_lo%ikit_procs_assignment(ntheta0,naky))
    g_lo%ikit_procs_assignment(:,:)%mine = .false.
    g_lo%ikit_procs_assignment(:,:)%num_procs = 0
    do iglo=g_lo%llim_proc,g_lo%ulim_proc
       ik=ik_idx(g_lo,iglo)
       it=it_idx(g_lo,iglo)
       il=il_idx(g_lo,iglo)
       ie=ie_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)
       ik_max=MAX(ik,ik_max)
       it_max=MAX(it,it_max)
       il_max=MAX(il,il_max)
       ie_max=MAX(ie,ie_max)
       is_max=MAX(is,is_max)
       ik_min=MIN(ik,ik_min)
       it_min=MIN(it,it_min)
       il_min=MIN(il,il_min)
       ie_min=MIN(ie,ie_min)
       is_min=MIN(is,is_min)
       ik_local(ik)=.true.
       it_local(it)=.true.
       il_local(il)=.true.
       ie_local(ie)=.true.
       is_local(is)=.true.
       if(.not. g_lo%ikit_procs_assignment(it,ik)%mine) then
          g_lo%ikit_procs_assignment(it,ik)%mine = .true.
          g_lo%ikitrange = g_lo%ikitrange + 1
       end if
    enddo

    ! If the number of naky*ntheta0 points is less than the number of MPI
    ! processes then disable gf_local_fields
    if(nproc .lt. naky*ntheta0 .and. gf_local_fields) then
       if(proc0) then
          write(*,*) 'Disabling gf_local_fields in layouts because not enough MPI processes are being used'
          write(*,*) 'To use gf_local_fields with this simulation ensure you use more than ',naky*ntheta0
          write(*,*) 'MPI processes.  You will also need to set the field_option to gf_local in the fields'
          write(*,*) 'properties list, and gf_lo_integrate to true in dist_fn properties'
       end if
       gf_local_fields = .false.
    end if                                                                                                                                                                                                 

    ! Record the ik,it points this process have in a contigious array.
    allocate(g_lo%local_ikit_points(g_lo%ikitrange))
    tmp = 1
    do ik=1,naky
       do it=1,ntheta0
          colour = 0
          if(g_lo%ikit_procs_assignment(it,ik)%mine) then
             g_lo%local_ikit_points(tmp)%ik = ik 
             g_lo%local_ikit_points(tmp)%it = it 
             tmp = tmp + 1
             colour = 1
          end if
          if(gf_local_fields) then
             call split(colour,tmpcomm)
             if(g_lo%ikit_procs_assignment(it,ik)%mine) then
                g_lo%ikit_procs_assignment(it,ik)%comm = tmpcomm
             else
                call free_comm(tmpcomm)
             end if
          end if
       end do
    end do


    if(gf_local_fields) then
       g_lo%max_ikit_nprocs = 0
       do ik=1,naky
          do it=1,ntheta0
             tmp = 0
             tmp_proc_list = 0
             if(g_lo%ikit_procs_assignment(it,ik)%mine) then
                tmp = 1
                ! +1 is required here to ensure that proc0 of the sub comm is also counted
                tmp_proc_list(iproc+1) = g_lo%ikit_procs_assignment(it,ik)%comm%iproc+1
             end if
             call sum_allreduce(tmp)
             call sum_allreduce(tmp_proc_list)
             ! proc_list is the array process ranks (of the processes that own this ik,it point) in the communicator used for g_lo
             allocate(g_lo%ikit_procs_assignment(it,ik)%proc_list(tmp))
             ! sub_proc_list is the array of process ranks (of the processes that own this ik,it point) in the sub-communicator used for the fields calculations
             allocate(g_lo%ikit_procs_assignment(it,ik)%sub_proc_list(tmp))
             g_lo%ikit_procs_assignment(it,ik)%num_procs = tmp
             g_lo%max_ikit_nprocs = max(tmp, g_lo%max_ikit_nprocs)
             tmp = 1
             ! Iterate through all the processes and record if they are part of this ik,it point in g_lo
             do i=1,nproc
                if(tmp_proc_list(i) .gt. 0) then
                   g_lo%ikit_procs_assignment(it,ik)%proc_list(tmp) = i-1
                   g_lo%ikit_procs_assignment(it,ik)%sub_proc_list(tmp) = tmp_proc_list(i)-1
                   tmp = tmp + 1
                end if
             end do
          end do
       end do
       
       tmp = g_lo%ikitrange
       call max_reduce(tmp,0)
#ifdef DEBUG
       if(proc0) then
          write(*,*) 'maximum number of ik,it points a given process owns in g_lo',tmp
       end if
#endif
       tmp = g_lo%ikitrange
       call sum_reduce(tmp, 0)
#ifdef DEBUG
       if(proc0) then
          write(*,*) 'average number of ik,ik points',tmp/nproc
       end if
#endif
    end if
    !Store dimension locality
    g_lo%y_local=all(ik_local)
    g_lo%x_local=all(it_local)
    g_lo%l_local=all(il_local)
    g_lo%e_local=all(ie_local)
    g_lo%s_local=all(is_local)
    
    !//Deallocate locality data. Note this may actually be useful in some
    !places, so we may want to think about attaching it to the g_lo object
    deallocate(ik_local,it_local,il_local,ie_local,is_local)

    !Now store min/max in appropriate dimension vars
    !This data isn't currently used but is useful to allow us to
    !replace loops over iglo with nested loops over the 5 generic dimensions
    !which are certain to be in optimised memory order.
    !Here is an example of this.
    !
    !ORIGINAL ROUTINE
    !subroutine Test
    ! use gs2_layouts, only: g_lo, ik_idx,it_idx
    ! use theta_grids, only: ntgrid
    ! use kt_grids, only: naky, ntheta0
    ! implicit none
    ! integer :: iglo, it,ik
    ! real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lim_alloc)::g
    ! real, dimension(-ntgrid:ntgrid,naky,ntheta0) :: total
    ! do iglo=g_lo%llim_proc,g_lim_alloc
    !     it=it_idx(g_lo,iglo)
    !     kt=ik_idx(g_lo,iglo)
    !     total(:,ik,it)=total(:,ik,it)+g(:,1,iglo)+g(:,2,iglo)
    ! enddo
    !end subroutine
    !
    !NEW ROUTINE
    !subroutine Test
    ! use gs2_layouts, only: g_lo, ik_idx,it_idx
    ! use theta_grids, only: ntgrid
    ! use kt_grids, only: naky, ntheta0
    ! implicit none
    ! integer :: iglo, d1,d2,d3,d4,d5
    ! integer, pointer :: it,ik
    ! integer, dimension(5),target :: dind
    ! real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lim_alloc)::g
    ! real, dimension(-ntgrid:ntgrid,naky,ntheta0) :: total
    ! EQUIVALENCE(dind(1),d1),(dind(2),d2),(dind(3),d3),(dind(4),d4),(dind(5),d5) !<--This just allows us to access d1..d5 using array notation
    ! !Setup ik and it to point to appropriate d1..d5
    ! ik=>dind(g_lo%ik_ord)
    ! it=>dind(g_lo%it_ord)
    ! iglo=g_lo%llim_proc
    ! do d5=g_lo%d5_min,g_lo%d5_max
    !    do d4=.... ; do d3=.... ; do d2=.... ; do d1=....
    !          total(:,ik,it)=total(:,ik,it)+g(:,1,iglo)+g(:,2,iglo)
    !          iglo=iglo+1
    !    enddo;enddo;enddo;enddo
    ! enddo
    !end subroutine
    select case (g_lo%ik_ord)
    case (1)
       g_lo%d1_min=ik_min
       g_lo%d1_max=ik_max
    case (2)
       g_lo%d2_min=ik_min
       g_lo%d2_max=ik_max
    case (3)
       g_lo%d3_min=ik_min
       g_lo%d3_max=ik_max
    case (4)
       g_lo%d4_min=ik_min
       g_lo%d4_max=ik_max
    case (5)
       g_lo%d5_min=ik_min
       g_lo%d5_max=ik_max
    end select

    select case (g_lo%it_ord)
    case (1)
       g_lo%d1_min=it_min
       g_lo%d1_max=it_max
    case (2)
       g_lo%d2_min=it_min
       g_lo%d2_max=it_max
    case (3)
       g_lo%d3_min=it_min
       g_lo%d3_max=it_max
    case (4)
       g_lo%d4_min=it_min
       g_lo%d4_max=it_max
    case (5)
       g_lo%d5_min=it_min
       g_lo%d5_max=it_max
    end select
       
    select case (g_lo%il_ord)
    case (1)
       g_lo%d1_min=il_min
       g_lo%d1_max=il_max
    case (2)
       g_lo%d2_min=il_min
       g_lo%d2_max=il_max
    case (3)
       g_lo%d3_min=il_min
       g_lo%d3_max=il_max
    case (4)
       g_lo%d4_min=il_min
       g_lo%d4_max=il_max
    case (5)
       g_lo%d5_min=il_min
       g_lo%d5_max=il_max
    end select

    select case (g_lo%ie_ord)
    case (1)
       g_lo%d1_min=ie_min
       g_lo%d1_max=ie_max
    case (2)
       g_lo%d2_min=ie_min
       g_lo%d2_max=ie_max
    case (3)
       g_lo%d3_min=ie_min
       g_lo%d3_max=ie_max
    case (4)
       g_lo%d4_min=ie_min
       g_lo%d4_max=ie_max
    case (5)
       g_lo%d5_min=ie_min
       g_lo%d5_max=ie_max
    end select

    select case (g_lo%is_ord)
    case (1)
       g_lo%d1_min=is_min
       g_lo%d1_max=is_max
    case (2)
       g_lo%d2_min=is_min
       g_lo%d2_max=is_max
    case (3)
       g_lo%d3_min=is_min
       g_lo%d3_max=is_max
    case (4)
       g_lo%d4_min=is_min
       g_lo%d4_max=is_max
    case (5)
       g_lo%d5_min=is_min
       g_lo%d5_max=is_max
    end select

    g_lo%ik_min=ik_min
    g_lo%ik_max=ik_max
    g_lo%it_min=it_min
    g_lo%it_max=it_max
    g_lo%il_min=il_min
    g_lo%il_max=il_max
    g_lo%ie_min=ie_min
    g_lo%ie_max=ie_max
    g_lo%is_min=is_min
    g_lo%is_max=is_max
    
    !Work out if x comes before y in layout
    do iglo=1,LEN_TRIM(layout)
       if(layout(iglo:iglo).eq.'x') then
          g_lo%x_before_y=.true.
          exit
       elseif (layout(iglo:iglo).eq.'y') then
          g_lo%x_before_y=.false.
          exit
       endif
    enddo

    !Now work out how the compound dimension splits amongst processors
    nproc_tmp=nproc
    !/Initialise
    nproc_dim=-1
    dim_divides=.false.
    dim_local=.false.
    do idim=5,1,-1

       !If we are down to 1 proc then everything else is local
       !so say dim_divides and exit loop
       !If nproc_tmp is 0 then this indicates the last dimension
       !had more entries than procs available
       if(nproc_tmp.le.1) then
          dim_divides(idim:1:-1)=.true.
          dim_local(idim:1:-1)=.true.
          exit
       endif
       
       !Does the current dimension divide evenly
       !amongst the remaining processors
       tmp_r=(1.0*nproc_tmp)/g_lo%dim_size(idim)

       !Check if the current dimension is larger than the number of processors
       !if so then we can still divide nicely but can't rely on division being integer
       !Invert division to check!
       if(int(tmp_r).eq.0) then
          tmp_r=g_lo%dim_size(idim)/(1.0*nproc_tmp)
       endif

       !Store logical to indicate if things divide perfectly
       dim_divides(idim)=((tmp_r-int(tmp_r)).eq.0)

       !If this dimension doesn't divide perfectly then we exit
       !loop and leave other dims as "un-initialised"
       if(.not.dim_divides(idim)) then
          nproc_dim(idim)=int(tmp_r-int(tmp_r))
          exit
       endif

       !If we do split nicely then store how many procs there are
       !to deal with the subsequent dims
       !If this is zero then it indicates a dimension that is not
       !fully distributed.
       nproc_dim(idim)=nproc_tmp/g_lo%dim_size(idim)

       !Update the available processor count
       nproc_tmp=nproc_dim(idim)
    enddo

    !Now print some stuff for debugging
!     if (proc0) then
!        write(6,'("Idim | Divides | Nproc_per_block | Size | Local")')
!        do idim=5,1,-1
!           write(6,'(I0,6(" "),L,9(" "),I0,16(" "),I0,8(" "),L)') idim,dim_divides(idim),nproc_dim(idim),g_lo%dim_size(idim),dim_local(idim)
!        enddo
!     endif

    !if (.not. (layout .eq. 'xyles' .or. layout .eq. 'yxles')) then
      !intmom_sub = .false.
      !intspec_sub = .false.
    !end if

    !Now use the dimension splitting to work out if the various subcommunicators are
    !allowed.
    if(.not.(dim_divides(g_lo%is_ord).and.dim_divides(g_lo%it_ord).and.dim_divides(g_lo%ik_ord))) then
       if(intmom_sub.and.proc0) write(error_unit(),'("Disabling intmom_sub -- is,it,ik split nicely ? ",3(L1," "))') dim_divides(g_lo%is_ord),dim_divides(g_lo%it_ord),dim_divides(g_lo%ik_ord)
       intmom_sub=.false.
    endif
    if(.not.(dim_divides(g_lo%it_ord).and.dim_divides(g_lo%ik_ord))) then
       if(intspec_sub.and.proc0) write(error_unit(),'("Disabling intspec_sub -- it,ik split nicely ? ",2(L1," "))') dim_divides(g_lo%it_ord),dim_divides(g_lo%ik_ord)
       intspec_sub=.false.
    endif
    !Note that because we currently need to gather amongst LES blocks at the end of
    !integrate_species we need to make sure that LES blocks are also sensible
    !This could be removed if we don't want to gather in integrate_species
    if(.not.(dim_divides(g_lo%il_ord).and.dim_divides(g_lo%ie_ord).and.dim_divides(g_lo%is_ord))) then
       if(intspec_sub.and.proc0) write(error_unit(),'("Disabling intspec_sub -- il,ie,is split nicely ? ",3(L1," "))') dim_divides(g_lo%il_ord),dim_divides(g_lo%ie_ord),dim_divides(g_lo%is_ord)
       intspec_sub=.false.
    endif

    !Update various flags based on locality
    !/Don't use subcommunicator and  gather for integrate_species if x and y are entirely local
    if((g_lo%x_local.and.g_lo%y_local)) then
       if(proc0.and.intspec_sub) write(error_unit(),'("Disabling intspec_sub as x and y are local")')
       intspec_sub=.false.
    endif

    !We could(should) use a logical allreduction here to ensure that ALL processors have the same flag settings
    !(if not then we'll get a program hang in the following splits).
    !The logic above **should** ensure everyone has the right answer anyway.

    !Now work out which xy block we live in for velocity space comms
    !Set sub communicators to mp_comm if we want to keep collective comms global.
    !NOTE: We might want to consider using a communicator type so that we can attach
    !some useful information to the communicators such as number of procs in comm,
    !this procs rank in comm, if this comm is equivalent to mp_comm (i.e. if split
    !doesn't actually split the global comm) etc.
    !This could help reduce the amount of mp(i) calls placed throughout the code
    !as things like the number of procs and the local rank are useful when using the sub
    !comms but don't currently have a logical place to be stored and hence tend to be
    !(re)calculated when needed.
    if(intspec_sub) then
       !This is for unique xy blocks
       col=1
       mycol=0
       do it=1,ntheta0
          do ik=1,naky
             !This is inefficient but who cares
             if(it_min.eq.it.and.ik_min.eq.ik) then
                mycol=col
                exit
             endif
             col=col+1
          enddo
       enddo

       !Now split the global communicator
       call split(mycol,g_lo%xyblock_comm)
       !This is for unique les blocks
       col=1
       mycol=0
       do is=1,nspec
          do il=1,nlambda
             do ie=1,negrid
                !This is inefficient but who cares
                if(il_min.eq.il.and.ie_min.eq.ie.and.is_min.eq.is) then
                   mycol=col
                   exit
                endif
                col=col+1
             enddo
          enddo
       enddo

       !Now split
       call split(mycol,g_lo%lesblock_comm)
    else
       g_lo%xyblock_comm=mp_comm
       g_lo%lesblock_comm=mp_comm
    endif

    if(intmom_sub)then
       !NOTE: A simple test that all blocks are equal would be:
       !tmp=it_max*ik_max
       !call max_allreduce(tmp,g_lo%xyblock_comm)
       !if(tmp.ne.(it_max*ik_max)) then
       !   print'(a)',"ERROR: In subcomm g_lo%xyblock_comm"
       !   print'(a,I0,a)',"Proc with global rank ",iproc," doesn't have the same it*ik extent as other procs in its subcomm"
       !   call mp_abort("Sub-comm creation error -- Consider changing the number of procs/number of xy gridpoints or setting intmom_sub=.false.")
       !endif
       !Similar tests can be done for all subcomms created here

       !Now work out which xys block we live in for velocity space comms
       !This is for unique xys blocks
       col=1
       mycol=0
       do it=1,ntheta0
          do ik=1,naky
             do is=1,nspec
                !This is inefficient but who cares
                if(it_min.eq.it.and.ik_min.eq.ik.and.is_min.eq.is) then
                   mycol=col
                   exit
                endif
                col=col+1
             enddo
          enddo
       enddo

       !Now split the global communicator
       call split(mycol,g_lo%xysblock_comm)
    else
       g_lo%xysblock_comm=mp_comm
    endif

    !Now allocate and fill various arrays
    ! EGH changed the line below. We assume that if this function
    ! has been called then g_lo%les_kxky_range has not been allocated.
    !if(intspec_sub.and.(.not.allocated(g_lo%les_kxky_range))) then
    if (allocated(g_lo%les_kxky_range)) then
       call mp_abort("g_lo%les_kxky_range should never be allocated&
       &at this point", .true.)
    end if
    
    if(intspec_sub) then
       !->Kx-Ky range for procs in lesblock_comm sub comm
       !Get number of procs in sub-comm
       call nproc_comm(g_lo%lesblock_comm,nproc_subcomm)
       !Get rank of proc in sub-comm
       call rank_comm(g_lo%lesblock_comm,rank_subcomm)
       !Allocate array | This will store the start and stop indices 
       !for flattened kxky indices in velocity space independent (e.g. field)
       !variables
       allocate(g_lo%les_kxky_range(2,0:nproc_subcomm-1))
       allocate(les_kxky_range(0:nproc_subcomm*2-1))
       !Initialise as we do sum_allreduce later
       !can remove if we use a gather
       les_kxky_range=0
       
       !Now calculate our personal kxky range
       it_min=it_idx(g_lo,g_lo%llim_proc)
       it_max=it_idx(g_lo,g_lo%ulim_alloc)
       ik_min=ik_idx(g_lo,g_lo%llim_proc)
       ik_max=ik_idx(g_lo,g_lo%ulim_alloc)
       
       if(g_lo%x_before_y) then
          les_kxky_range(rank_subcomm*2)=it_min-1+ntheta0*(ik_min-1)
          tmp=les_kxky_range(rank_subcomm*2)
          do iglo=g_lo%llim_proc,g_lo%ulim_alloc
             tmp=MAX(tmp,it_idx(g_lo,iglo)-1+ntheta0*(ik_idx(g_lo,iglo)-1))
          enddo
          les_kxky_range(rank_subcomm*2+1)=tmp
       else
          les_kxky_range(rank_subcomm*2)=ik_min-1+naky*(it_min-1)
          tmp=les_kxky_range(rank_subcomm*2)
          do iglo=g_lo%llim_proc,g_lo%ulim_alloc
             tmp=MAX(tmp,ik_idx(g_lo,iglo)-1+naky*(it_idx(g_lo,iglo)-1))
          enddo
          les_kxky_range(rank_subcomm*2+1)=tmp
       endif

       !Now gather values and store in g_lo, sum_allreduce would also work
       !    call mpi_allgather(les_kxky_range(rank_subcomm*2:rank_subcomm*2+1),2,MPI_INTEGER,&
       !         les_kxky_range,2,MPI_INTEGER,g_lo%lesblock_comm)
       !Should really use gather but this is easier
       call sum_allreduce_sub(les_kxky_range,g_lo%lesblock_comm)
       do ip=0,nproc_subcomm-1
          g_lo%les_kxky_range(:,ip)=les_kxky_range(ip*2:ip*2+1)
       enddo

       !Free memory
       deallocate(les_kxky_range)
    endif

  end subroutine init_dist_fn_layouts

  !> Return the species index, given g-layout and a global index
  elemental function is_idx_g (lo, i)
    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_g=1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
  end function is_idx_g

  !> Return the lambda index, given g-layout and a global index
  elemental function il_idx_g (lo, i)
    implicit none
    integer :: il_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_g=1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
  end function il_idx_g

  !> Return the energy index, given g-layout and a global index
  elemental function ie_idx_g (lo, i)
    implicit none
    integer :: ie_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_g=1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
  end function ie_idx_g

  !> Return the kx index, given g-layout and a global index
  elemental function it_idx_g (lo, i)
    implicit none
    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_g=1+mod((i-lo%llim_world)/lo%it_comp,lo%ntheta0)
  end function it_idx_g

  !> Return the ky index, given g-layout and a global index
  elemental function ik_idx_g (lo, i)
    implicit none
    integer :: ik_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_g=1+mod((i-lo%llim_world)/lo%ik_comp,lo%naky)
  end function ik_idx_g

  !> Return the processor id that has the point i in the g_lo layout
  elemental function proc_id_g (lo, i)
#ifdef SHMEM
    use mp, only : nproc
#endif
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

#ifdef SHMEM
    integer k

    proc_id_g = nproc - 1
    do k =1, nproc-1
       if ( i < lo%proc_map(k)) then 
          proc_id_g = k - 1
          exit
       endif
    enddo
#else
    proc_id_g = i/lo%blocksize
#endif

  end function proc_id_g

  !> Return the global index in the g_lo layout given the dimensional indices
  elemental function idx_g (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    integer, dimension(5) :: i
    i(lo%ik_ord) = ik - 1
    i(lo%it_ord) = it - 1
    i(lo%il_ord) = il - 1
    i(lo%ie_ord) = ie - 1
    i(lo%is_ord) = is - 1
    idx_g = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*(i(4) + lo%dim_size(4)*i(5))))

  end function idx_g

  !> Return whether a point is local to a processor given the dimensional indices
  elemental function idx_local_g (lo, ik, it, il, ie, is)
    implicit none
    logical :: idx_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    idx_local_g = idx_local(lo, idx(lo, ik, it, il, ie, is))
  end function idx_local_g

  !> Return whether a point is local to a processor given the global index
  elemental function ig_local_g (lo, ig)
    implicit none
    logical :: ig_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

#ifndef SHMEM
    ig_local_g = lo%iproc == proc_id(lo, ig)
#else
    ig_local_g = .false.
    if ( ig >= lo%llim_proc .and. ig <= lo%ulim_proc ) &
         ig_local_g = .true.
#endif
  end function ig_local_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Field layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! M_class(i) = number of i_th supercell
! N_class(i) = size of i_th supercell
! i_class = number of classes of (different sized) supercells.  

  !> FIXME : Add documentation   
  subroutine init_fields_layouts (nfield, nindex, naky, ntheta0, &
       M_class, N_class, i_class)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nfield, nindex, naky, ntheta0
    integer, dimension(:), intent (in) :: M_class, N_class
    integer, intent (in) :: i_class
    integer :: i, utmp, btmp

    if (initialized_fields_layouts) return
    initialized_fields_layouts = .true.
    
    allocate (f_lo(i_class))
    do i = 1, i_class
       allocate (f_lo(i)%ik(M_class(i), N_class(i)))
       allocate (f_lo(i)%it(M_class(i), N_class(i)))
    end do

    do i = 1, i_class
       f_lo(i)%iproc = iproc
       f_lo(i)%nfield = nfield
       f_lo(i)%nidx = nindex                 ! cell size
       f_lo(i)%ntgrid = (nindex/nfield-1)/2
       f_lo(i)%nindex = nindex*N_class(i)    ! supercell size
       f_lo(i)%naky = naky
       f_lo(i)%ntheta0 = ntheta0
       f_lo(i)%M_class = M_class(i)
       f_lo(i)%N_class = N_class(i)
       f_lo(i)%i_class = i_class
       f_lo(i)%llim_world = 0
       f_lo(i)%ulim_world = f_lo(i)%nindex*M_class(i) - 1
       if (local_field_solve) then    ! guarantee local operations for matrix inversion
          utmp = M_class(i) - 1
          btmp = utmp/nproc + 1
          f_lo(i)%blocksize = f_lo(i)%nindex*btmp
       else
          f_lo(i)%blocksize = f_lo(i)%ulim_world/nproc + 1
       end if
       f_lo(i)%llim_proc = f_lo(i)%blocksize*iproc
       f_lo(i)%ulim_proc = min(f_lo(i)%ulim_world, f_lo(i)%llim_proc + f_lo(i)%blocksize - 1)
       f_lo(i)%ulim_alloc = max(f_lo(i)%llim_proc, f_lo(i)%ulim_proc)
    end do

  end subroutine init_fields_layouts

  !> FIXME : Add documentation   
  subroutine finish_fields_layouts
    implicit none

    integer :: i,f_size

    if(allocated(f_lo))then
       f_size = size(f_lo)
       do i=1,f_size
          if(associated(f_lo(i)%ik)) deallocate (f_lo(i)%ik)
          if(associated(f_lo(i)%it)) deallocate (f_lo(i)%it)
       enddo
       deallocate(f_lo)    
    endif
    initialized_fields_layouts = .false.
  end subroutine finish_fields_layouts

  !> FIXME : Add documentation
  !!
  !! Finish routine so that gs2 can be restarted with a different
  !! number of processors by trinity
  subroutine finish_layouts
    implicit none
    call finish_fields_layouts
    call finish_jfields_layouts

    initialized = .false.
    initialized_layouts = .false.
    initialized_le_layouts = .false.
    initialized_gf_layouts = .false.
    initialized_energy_layouts = .false.
    initialized_lambda_layouts = .false.
    initialized_parity_layouts = .false.
    initialized_accel_transform_layouts = .false.
    initialized_x_transform = .false.
    initialized_y_transform = .false.

    call layouts_config%reset()
  end subroutine finish_layouts

  !> FIXME : Add documentation
  subroutine finish_dist_fn_layouts

    use layouts_type, only: deallocate_layout

    call deallocate_layout(g_lo)
    initialized_dist_fn_layouts = .false.

  end subroutine finish_dist_fn_layouts

  !> FIXME : Add documentation
  function ik_idx_f (lo, i)
    implicit none
    integer :: ik_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_f = lo%ik(im_idx(lo,i),in_idx(lo,i))
  end function ik_idx_f

  !> FIXME : Add documentation
  function it_idx_f (lo, i)
    implicit none
    integer :: it_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_f = lo%it(im_idx(lo,i),in_idx(lo,i))
  end function it_idx_f

  !> FIXME : Add documentation  
  function im_idx_f (lo, i)
    implicit none
    integer :: im_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    im_idx_f = 1 + mod((i - lo%llim_world)/lo%nindex, lo%M_class)
  end function im_idx_f

  !> FIXME : Add documentation  
  function if_idx_f (lo, i)
    implicit none
    integer :: if_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_f = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_f

  !> FIXME : Add documentation  
  function idx_f (lo, if, im)
    implicit none
    integer :: idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, im
    idx_f = if-1 + lo%nindex*(im-1)
  end function idx_f

  !> FIXME : Add documentation  
  function proc_id_f (lo, i)
    implicit none
    integer :: proc_id_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_f = i/lo%blocksize
  end function proc_id_f

  !> FIXME : Add documentation  
  function idx_local_f (lo, if, im)
    implicit none
    logical :: idx_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, im

    idx_local_f = idx_local(lo, idx(lo, if, im))
  end function idx_local_f

  !> FIXME : Add documentation  
  function ig_local_f (lo, ig)
    implicit none
    logical :: ig_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_f = lo%iproc == proc_id(lo, ig)
  end function ig_local_f

  !> FIXME : Add documentation  
  function in_idx_f (lo, i)
    implicit none
    integer :: in_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    in_idx_f = 1 + mod((i - lo%llim_world)/lo%nidx, lo%N_class)
  end function in_idx_f

  !> FIXME : Add documentation  
  function ig_idx_f (lo, i)
    implicit none
    integer :: ig_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_f = -lo%ntgrid + mod((i - lo%llim_world)/lo%nfield, (2*lo%ntgrid+1))
  end function ig_idx_f

  !> FIXME : Add documentation  
  function ij_idx_f (lo, ig, if, n)
    implicit none
    integer :: ij_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, if, n
    ij_idx_f = ig+lo%ntgrid + (2*lo%ntgrid+1)*(if-1 + lo%nfield*(n-1))
  end function ij_idx_f

  !> FIXME : Add documentation  
  function ifield_idx_f (lo, i)
    implicit none
    integer :: ifield_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ifield_idx_f = 1 + mod((i - lo%llim_world), lo%nfield)
  end function ifield_idx_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fast field layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_jfields_layouts (nfield, nindex, naky, ntheta0, i_class)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nfield, nindex, naky, ntheta0, i_class
    integer :: jlo
    if (initialized_jfields_layouts) return
    initialized_jfields_layouts = .true.
    
    jf_lo%iproc = iproc
    jf_lo%nindex = nindex
    jf_lo%ntgrid = (nindex/nfield-1)/2
    jf_lo%nfield = nfield
    jf_lo%naky = naky
    jf_lo%ntheta0 = ntheta0
    jf_lo%llim_world = 0
    jf_lo%ulim_world = nindex*ntheta0*naky - 1
    jf_lo%blocksize = jf_lo%ulim_world/nproc + 1
    jf_lo%llim_proc = jf_lo%blocksize*iproc
    jf_lo%ulim_proc = min(jf_lo%ulim_world, jf_lo%llim_proc + jf_lo%blocksize - 1)
    jf_lo%ulim_alloc = max(jf_lo%llim_proc, jf_lo%ulim_proc)

    !Now work out the range of dimensions stored on this proc !<DD>
    jf_lo%ig_max=-jf_lo%ntgrid-1
    jf_lo%ig_min=jf_lo%ntgrid+1
    jf_lo%ifield_max=0
    jf_lo%ifield_min=nfield+1
    jf_lo%ik_max=0
    jf_lo%ik_min=naky+1
    jf_lo%it_max=0
    jf_lo%it_min=0
    do jlo=jf_lo%llim_proc,jf_lo%ulim_alloc
       jf_lo%ig_max=MAX(jf_lo%ig_max,ig_idx(jf_lo,jlo))
       jf_lo%ig_min=MIN(jf_lo%ig_min,ig_idx(jf_lo,jlo))
       jf_lo%ifield_max=MAX(jf_lo%ifield_max,ifield_idx(jf_lo,jlo))
       jf_lo%ifield_min=MIN(jf_lo%ifield_min,ifield_idx(jf_lo,jlo))
       jf_lo%ik_max=MAX(jf_lo%ik_max,ik_idx(jf_lo,jlo))
       jf_lo%ik_min=MIN(jf_lo%ik_min,ik_idx(jf_lo,jlo))
       jf_lo%it_max=MAX(jf_lo%it_max,it_idx(jf_lo,jlo))
       jf_lo%it_min=MIN(jf_lo%it_min,it_idx(jf_lo,jlo))
    enddo

    allocate (ij(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (mj(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (dj(i_class,jf_lo%llim_proc:jf_lo%ulim_alloc))
    ij = 1  ; mj = 1;  dj = 0
    
  end subroutine init_jfields_layouts

  !> FIXME : Add documentation  
  subroutine finish_jfields_layouts
    implicit none
    if(allocated(ij)) deallocate(ij)
    if(allocated(mj)) deallocate(mj)
    if(allocated(dj)) deallocate(dj)
    initialized_jfields_layouts = .false.
  end subroutine finish_jfields_layouts

  !> FIXME : Add documentation  
  function ik_idx_jf (lo, i)
    implicit none
    integer :: ik_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_jf = 1 + mod((i - lo%llim_world)/lo%nindex/lo%ntheta0, lo%naky)
  end function ik_idx_jf

  !> FIXME : Add documentation  
  function it_idx_jf (lo, i)
    implicit none
    integer :: it_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_jf = 1 + mod((i - lo%llim_world)/lo%nindex, lo%ntheta0)
  end function it_idx_jf

  !> FIXME : Add documentation
  !!
  !! @note Here if is the compound theta*field index
  function if_idx_jf (lo, i)
    implicit none
    integer :: if_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_jf = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_jf

  !> FIXME : Add documentation
  function ig_idx_jf (lo, i)
    implicit none
    integer :: ig_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_jf = -lo%ntgrid + mod(i - lo%llim_world, lo%ntgrid*2+1)
  end function ig_idx_jf

  !> FIXME : Add documentation
  function ifield_idx_jf (lo, i)
    implicit none
    integer :: ifield_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ifield_idx_jf = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid+1), lo%nfield)
  end function ifield_idx_jf

  !> FIXME : Add documentation  
  function idx_jf (lo, if, ik, it)
    implicit none
    integer :: idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it
    idx_jf = if-1 + lo%nindex*(it-1 + lo%ntheta0*(ik-1))
  end function idx_jf

  !> FIXME : Add documentation  
  function ij_idx_jf (lo, ig, if, ik, it)
    implicit none
    integer :: ij_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, if, ik, it
    ij_idx_jf = ig+lo%ntgrid + (2*lo%ntgrid+1)*(if-1+lo%nfield*(it-1+lo%ntheta0*(ik-1)))
  end function ij_idx_jf

  !> FIXME : Add documentation  
  function proc_id_jf (lo, i)
    implicit none
    integer :: proc_id_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_jf = i/lo%blocksize
  end function proc_id_jf

  !> FIXME : Add documentation  
  function idx_local_jf (lo, if, ik, it)
    implicit none
    logical :: idx_local_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it

    idx_local_jf = idx_local(lo, idx(lo, if, ik, it))
  end function idx_local_jf

  !> FIXME : Add documentation  
  function ig_local_jf (lo, ig)
    implicit none
    logical :: ig_local_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_jf = lo%iproc == proc_id(lo, ig)
  end function ig_local_jf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lambda-Energy layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_le_layouts (ntgrid, naky, ntheta0, nspec)
    use mp, only: iproc, nproc
    use file_utils, only: error_unit
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nspec
    integer :: ik_min, it_min, ig_min, is_min
    integer :: ik_max, it_max, ig_max, is_max
    integer :: ik, it, ig, is, ile

    if (initialized_le_layouts) return
    initialized_le_layouts = .true.
    
    le_lo%iproc = iproc
    le_lo%ntgrid = ntgrid
    le_lo%ntgridtotal=(2*ntgrid+1)
    le_lo%ntheta0 = ntheta0
    le_lo%naky = naky
    le_lo%nspec = nspec
    le_lo%llim_world = 0
    le_lo%ulim_world = (2*ntgrid+1) * naky * ntheta0 * nspec - 1
    le_lo%blocksize = le_lo%ulim_world / nproc + 1
    le_lo%llim_proc = le_lo%blocksize * iproc
    le_lo%ulim_proc &
         = min(le_lo%ulim_world, le_lo%llim_proc + le_lo%blocksize - 1)
    le_lo%ulim_alloc = max(le_lo%llim_proc, le_lo%ulim_proc)

    ! See g_lo init
    select case (layout)
    case ('yxels')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=3,ik_ord=2,is_ord=4)
    case ('yxles')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=3,ik_ord=2,is_ord=4)
    case ('lexys')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=2,ik_ord=3,is_ord=4)
    case ('lxyes')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=2,ik_ord=3,is_ord=4)
    case ('lyxes')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=3,ik_ord=2,is_ord=4)
    case ('xyles')
       call set_dimension_order(le_lo,ig_ord=1,it_ord=2,ik_ord=3,is_ord=4)
    end select
    
    le_lo%dim_size(le_lo%ig_ord)=le_lo%ntgridtotal
    le_lo%dim_size(le_lo%it_ord)=ntheta0
    le_lo%dim_size(le_lo%ik_ord)=naky
    le_lo%dim_size(le_lo%is_ord)=nspec

    le_lo%compound_count(1)=1
    le_lo%compound_count(2)=le_lo%compound_count(1)*le_lo%dim_size(1)
    le_lo%compound_count(3)=le_lo%compound_count(2)*le_lo%dim_size(2)
    le_lo%compound_count(4)=le_lo%compound_count(3)*le_lo%dim_size(3)

    le_lo%ig_comp=le_lo%compound_count(le_lo%ig_ord)
    le_lo%ik_comp=le_lo%compound_count(le_lo%ik_ord)
    le_lo%it_comp=le_lo%compound_count(le_lo%it_ord)
    le_lo%is_comp=le_lo%compound_count(le_lo%is_ord)
    
    !Work out min and max of the five compound indices
    !These are currently only used to slice the arrays used in the velocity space integrals
    !but could well be useful in a number of places.
    ik_max=0
    ik_min=naky+1
    it_max=0
    it_min=ntheta0+1
    ig_max=-ntgrid-1
    ig_min=ntgrid+1
    is_max=0
    is_min=nspec+1
    do ile=le_lo%llim_proc,le_lo%ulim_proc
       ik=ik_idx(le_lo,ile)
       it=it_idx(le_lo,ile)
       ig=ig_idx(le_lo,ile)
       is=is_idx(le_lo,ile)
       ik_max=MAX(ik,ik_max)
       it_max=MAX(it,it_max)
       ig_max=MAX(ig,ig_max)
       is_max=MAX(is,is_max)
       ik_min=MIN(ik,ik_min)
       it_min=MIN(it,it_min)
       ig_min=MIN(ig,ig_min)
       is_min=MIN(is,is_min)
    enddo

    !Store dimension locality
    le_lo%x_local=(it_min.eq.1).and.(it_max.eq.ntheta0)
    le_lo%y_local=(ik_min.eq.1).and.(ik_max.eq.naky)
    le_lo%t_local=(ig_min.eq.-ntgrid).and.(ig_max.eq.ntgrid)
    le_lo%s_local=(is_min.eq.1).and.(is_max.eq.nspec)

  end subroutine init_le_layouts

  !> FIXME : Add documentation  
  elemental function is_idx_le (lo, i)
    implicit none
    integer :: is_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_le = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_le

  !> FIXME : Add documentation    
  elemental function it_idx_le (lo, i)
    implicit none
    integer :: it_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_le = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
  end function it_idx_le

  !> FIXME : Add documentation    
  elemental function ik_idx_le (lo, i)
    implicit none
    integer :: ik_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_le = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
  end function ik_idx_le

  !> FIXME : Add documentation  
  elemental function ig_idx_le (lo, i)
    implicit none
    integer :: ig_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_le = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_le

  !> FIXME : Add documentation    
  elemental function idx_le (lo, ig, ik, it, is)
    implicit none
    integer :: idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, is
    integer, dimension(4) :: i
    i(lo%ik_ord) = ik - 1
    i(lo%it_ord) = it - 1
    i(lo%is_ord) = is - 1
    i(lo%ig_ord) = ig + lo%ntgrid
    idx_le = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*i(4)))
  end function idx_le

  !> FIXME : Add documentation  
  elemental function proc_id_le (lo, i)
    implicit none
    integer :: proc_id_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_le = i / lo%blocksize

  end function proc_id_le

  !> FIXME : Add documentation    
  elemental function idx_local_le (lo, ig, ik, it, is)
    implicit none
    logical :: idx_local_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, is
    idx_local_le = idx_local(lo, idx(lo, ig, ik, it, is))
  end function idx_local_le

  !> FIXME : Add documentation  
  elemental function ig_local_le (lo, ig)
    implicit none
    logical :: ig_local_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig
    ig_local_le = lo%iproc == proc_id(lo, ig)
  end function ig_local_le

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fields Local layout:  This is used to calculate velocity space integration
! and other fields related calculations on a restricted process count for 
! large core count simulations.  It takes naky and ntheta and distributed 
! that amongst theta so that a single naky,ntheta point is only owned by 
! a single process.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_gf_layouts (ntgrid, naky, ntheta0, negrid, nlambda, nspec)
    use mp, only: iproc, nproc, mp_abort
    use file_utils, only: error_unit
    use layouts_type, only: ikit
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, negrid, nlambda, nspec

    if (initialized_gf_layouts) return
    initialized_gf_layouts = .true.
    
    gf_lo%iproc = iproc
    gf_lo%nproc = nproc
    gf_lo%ntgrid = ntgrid
    gf_lo%ntgridtotal=(2*ntgrid+1)
    gf_lo%negrid = negrid
    gf_lo%nlambda = nlambda
    gf_lo%nspec = nspec
    gf_lo%ntheta0 = ntheta0
    gf_lo%naky = naky
    gf_lo%divisionblock = gf_lo%naky * gf_lo%ntheta0
    gf_lo%llim_world = 1
    gf_lo%blocksize = 1
    gf_lo%ulim_world = gf_lo%divisionblock * gf_lo%blocksize
! This logic has been setup to deal with all the possible cases when constructing the decomposition of
! gf_lo over processes, namely; having more gf_lo points than processes (i.e. each process has multiple gf_lo points),
! having the same number of gf_lo points and processes, and having less gf_lo points than processes (i.e. some processes 
! have 1 gf_lo point and some have none).  We have two ways of dealing with the final scenario (less points than processes), 
! the first (simple_gf_decomposition) just assigns the points to the first n processes (where n is the number of points) and leaves the 
! rest of the points empty.  The second way is to try and evenly distribute the points throughout the process space.
! As a consequence of dealing with both main scenarios (i.e. less or more points than processes) we use some of the
! gf_lo variables in different ways depending on the scenario presented.  More specifically, the two variables that 
! are used different in the different scenarios are largeregionlimit and largeblocklimit.  In the scenario where 
! we have more processes than points then largeregionlimit == process rank that is the last that has a large gap between
! it and the next point, and largeblocklimit == gf_lo point where that large gaps size stop; whereas in the scenario 
! where we have more blocks than processes, largeregionlimit == gf_lo point where large blocks stop and 
! largeblocklimit == process number which is the last to have large blocks.  It is not very clean to have to 
! confusion, so TODO: Consider refactoring this to make the meaning of the variables described above consistent, or
! have different variables for the different scenarios. 
! TODO: Currently the decompositions do not respect naky boundaries.  Because if the linked nature of supercells 
! in field calculations (i.e. a supercell is a single naky point for one or more ntheta points) if a process has more than 
! one gf_lo point it should make sure those points are still all within a single supercell (i.e. different ntheta for 
! the same naky).  This would make the supercell calculations more efficient in this case (processes would not be 
! involved in multiple sets of collective communications in the fields calcs).
    if((gf_lo%divisionblock / gf_lo%nproc) .eq. 0) then
! x == ntheta0 y == naky
! Here we have less xy points than processes so we will have one process assigned the xy point (the cell) and 
! that will be the cell head.  The other processes involved with that y point won't be assigned anything.  We aim 
! to choose cell head from the group of processes that are already have that x y data
! We are currently using a very simple (and stupid) way of mapping x y points to processes, just taking the 
! first N processes (where N is naky*ntheta0).
       if(simple_gf_decomposition) then
          gf_lo%largeblocklimit = 0
          gf_lo%largeblocksize = 0
          gf_lo%largeregionlimit = 0
          gf_lo%smallblocksize = gf_lo%blocksize
          gf_lo%smallgapsize = 0
          gf_lo%largegapsize = 0
          if(gf_lo%divisionblock/(gf_lo%iproc + 1) .ne. 0) then
             gf_lo%llim_proc = gf_lo%smallblocksize * iproc + 1
             gf_lo%ulim_proc = gf_lo%llim_proc + gf_lo%smallblocksize - 1
             gf_lo%xypoints = 1
          else
             gf_lo%llim_proc = gf_lo%llim_world
             gf_lo%ulim_proc = 0
             gf_lo%xypoints = 0
          end if
! This variant if the situation where we have less blocks than processes attempts to spread the blocks out through 
! the range of processes we have.  As such there will be gaps in between each assigned block where processes have 
! no blocks assigned to them.  It is designed to spread out the blocks as evenly as possible, and can cope with the 
! situation where the number of processes does not exactly divide by the number of blocks.  In that scenario there 
! are bigger gaps between processes below a certain limit, and then smaller gaps above that limit (the largeblocklimit).
! In this context largeregionlimit is the process id of the top process in the region that has larger blocks.
       else
          gf_lo%smallblocksize = 1
          gf_lo%largeblocksize = 1
          gf_lo%smallgapsize = gf_lo%nproc/gf_lo%divisionblock
          gf_lo%largegapsize = gf_lo%smallgapsize + 1
          gf_lo%largeblocklimit = mod(gf_lo%nproc,gf_lo%divisionblock)
          gf_lo%largeregionlimit = gf_lo%largeblocklimit * gf_lo%largegapsize
          if(gf_lo%iproc .le. gf_lo%largeregionlimit) then
            if(mod(gf_lo%iproc,gf_lo%largegapsize) .eq. 0) then 
               gf_lo%xypoints = gf_lo%largeblocksize
               gf_lo%llim_proc = (gf_lo%iproc / gf_lo%largegapsize) + 1
               gf_lo%ulim_proc = gf_lo%llim_proc
            else
               gf_lo%llim_proc = gf_lo%llim_world
               gf_lo%ulim_proc = 0
               gf_lo%xypoints = 0
            end if
         else
            if(mod((gf_lo%iproc - (gf_lo%largeregionlimit)),gf_lo%smallgapsize) .eq. 0) then 
               gf_lo%xypoints = gf_lo%smallblocksize
               gf_lo%llim_proc = gf_lo%largeblocklimit + ((gf_lo%iproc - (gf_lo%largeregionlimit))/gf_lo%smallgapsize) + 1
               gf_lo%ulim_proc = gf_lo%llim_proc
            else
               gf_lo%llim_proc = gf_lo%llim_world
               gf_lo%ulim_proc = 0
               gf_lo%xypoints = 0
            end if
         end if
      end if
    else if((gf_lo%divisionblock / gf_lo%nproc) .ge. 1 .and. mod(gf_lo%divisionblock, gf_lo%nproc) .eq. 0) then
! Here we have exactly the same number of xy points as processes so each get one, or the number of points divides 
! exactly between processes so each gets exactly the same amount. 
       gf_lo%largeblocklimit = 0
       gf_lo%smallblocksize = ((gf_lo%divisionblock / gf_lo%nproc)) * gf_lo%blocksize
       gf_lo%largeblocksize = 0
       gf_lo%largeregionlimit = 0
       gf_lo%llim_proc = gf_lo%smallblocksize * iproc + 1
       gf_lo%ulim_proc = gf_lo%llim_proc + gf_lo%smallblocksize - 1
       gf_lo%xypoints = gf_lo%smallblocksize
       gf_lo%smallgapsize = 0
       gf_lo%largegapsize = 0
    else
! Here we have more xy points than processes so some or all will have multiple xy points.  In this scenario
! the large region limit is the top of the igf domain where the large blocks finish and we go into small blocks,
! and the large block limit is the highest process number where the large block region finishes
! AJ TODO: Make this respect supercell boundaries so each process only has a single y point but multiple x points 
! rather than multiple y points.  This would benefit fields calculations (see discussions in previous comment).
       gf_lo%largeblocklimit = mod(gf_lo%divisionblock,gf_lo%nproc)
       gf_lo%smallblocksize = (gf_lo%divisionblock/gf_lo%nproc) * gf_lo%blocksize
       gf_lo%largeblocksize = ((gf_lo%divisionblock/gf_lo%nproc)+1) * gf_lo%blocksize
       gf_lo%largeregionlimit = gf_lo%largeblocksize * gf_lo%largeblocklimit
       gf_lo%smallgapsize = 0
       gf_lo%largegapsize = 0
       if(gf_lo%iproc .lt. gf_lo%largeblocklimit) then
          gf_lo%xypoints = gf_lo%largeblocksize
          gf_lo%llim_proc = gf_lo%largeblocksize * iproc + 1
          gf_lo%ulim_proc = gf_lo%llim_proc + gf_lo%largeblocksize - 1
       else
          gf_lo%llim_proc = (gf_lo%largeblocklimit * gf_lo%largeblocksize) + (gf_lo%iproc - gf_lo%largeblocklimit) * gf_lo%smallblocksize + 1
          gf_lo%ulim_proc = gf_lo%llim_proc + gf_lo%smallblocksize - 1
          gf_lo%xypoints = gf_lo%smallblocksize
       end if
    end if
        
    gf_lo%ulim_alloc = gf_lo%ulim_proc

  end subroutine init_gf_layouts

  !> FIXME : Add documentation    
  elemental function it_idx_gf (lo, i)
    implicit none
    integer :: it_idx_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_gf = 1 + mod((i - lo%llim_world), lo%ntheta0)
  end function it_idx_gf

  !> FIXME : Add documentation    
  elemental function ik_idx_gf (lo, i)
    implicit none
    integer :: ik_idx_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_gf = 1 + ((i - lo%llim_world)/lo%ntheta0)
  end function ik_idx_gf

  !> FIXME : Add documentation    
  elemental function idx_gf (lo, ik, it)
    implicit none
    integer :: idx_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it
    idx_gf = (ik-1)*lo%ntheta0 + it
  end function idx_gf

  !> FIXME : Add documentation  
  elemental function proc_id_gf (lo, i)
    implicit none
    integer :: proc_id_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

! For this scenario we have more blocks than processes so each process has multiple blocks, but the number of blocks
! does not exactly divide by the number of processes, so some have larger blocks than others.
    if(lo%largeblocksize .gt. 1) then
       if(i .gt. lo%largeregionlimit) then
           proc_id_gf = (((i - (lo%largeblocksize * lo%largeblocklimit) - lo%llim_world)/lo%smallblocksize) + lo%largeblocklimit)
       else
           proc_id_gf = (i - lo%llim_world) / lo%largeblocksize
       end if
! For this scenario we have less blocks than processes so we try to spread the blocks between processes.
    else if(lo%largeblocksize .eq. 1) then
       if(i .gt. lo%largeblocklimit) then
          proc_id_gf = lo%largeregionlimit + ((i - lo%largeblocklimit - lo%llim_world)*lo%smallgapsize)
       else
          proc_id_gf = (i - lo%llim_world) * lo%largegapsize
       end if
! For this scenario we either have exactly the same number of blocks as processes, or the number of blocks divides exactly
! by the number of processes so each process has the same size blocks of the domain.
! It is also used by the "simple_gf_decomposition" functionality for the scenario where we have less blocks than processes and we 
! are not trying to distributed those blocks through the whole process space, just assigning them to the first n processes 
! (where n is the number of xypoints we have in the gf_lo domain).
    else
       proc_id_gf = (i - lo%llim_world) / lo%smallblocksize
    end if

  end function proc_id_gf

  !> FIXME : Add documentation    
  elemental function idx_local_gf (lo, ik, it)
    implicit none
    logical :: idx_local_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it
    idx_local_gf = idx_local(lo, idx(lo, ik, it))
  end function idx_local_gf

  !> FIXME : Add documentation  
  elemental function ig_local_gf (lo, ig)
    implicit none
    logical :: ig_local_gf
    type (gf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig
    ig_local_gf = lo%iproc == proc_id(lo, ig)
  end function ig_local_gf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation
  subroutine init_energy_layouts &
       (ntgrid, naky, ntheta0, nlambda, nspec)
    use mp, only: iproc, nproc
    use file_utils, only: error_unit
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, nspec

    if (initialized_energy_layouts) return
    initialized_energy_layouts = .true.
    
    e_lo%iproc = iproc
    e_lo%ntgrid = ntgrid
    e_lo%nsign = 2
    e_lo%ntheta0 = ntheta0
    e_lo%naky = naky
    e_lo%nspec = nspec
    e_lo%nlambda = nlambda
    e_lo%ntgridtotal=2*ntgrid+1
    e_lo%llim_world = 0
    e_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*nlambda*nspec*2 - 1
    e_lo%blocksize = e_lo%ulim_world/nproc + 1
    e_lo%llim_proc = e_lo%blocksize*iproc
    e_lo%ulim_proc &
         = min(e_lo%ulim_world, e_lo%llim_proc + e_lo%blocksize - 1)
    e_lo%ulim_alloc = max(e_lo%llim_proc, e_lo%ulim_proc)

    !See g_lo init
    select case (layout)
    case ('yxels')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=4,ik_ord=3,il_ord=5,is_ord=6)
    case ('yxles')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=4,ik_ord=3,il_ord=5,is_ord=6)
    case ('lexys')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=4,ik_ord=5,il_ord=3,is_ord=6)
    case ('lxyes')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=4,ik_ord=5,il_ord=3,is_ord=6)
    case ('lyxes')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=5,ik_ord=4,il_ord=3,is_ord=6)
    case ('xyles')
       call set_dimension_order(e_lo,ig_ord=1,isgn_ord=2,it_ord=3,ik_ord=4,il_ord=5,is_ord=6)
    end select

    e_lo%dim_size(e_lo%ig_ord)=e_lo%ntgridtotal
    e_lo%dim_size(e_lo%isgn_ord)=2
    e_lo%dim_size(e_lo%it_ord)=ntheta0
    e_lo%dim_size(e_lo%ik_ord)=naky
    e_lo%dim_size(e_lo%il_ord)=nlambda
    e_lo%dim_size(e_lo%is_ord)=nspec
    
    e_lo%compound_count(1)=1
    e_lo%compound_count(2)=e_lo%compound_count(1)*e_lo%dim_size(1)
    e_lo%compound_count(3)=e_lo%compound_count(2)*e_lo%dim_size(2)
    e_lo%compound_count(4)=e_lo%compound_count(3)*e_lo%dim_size(3)
    e_lo%compound_count(5)=e_lo%compound_count(4)*e_lo%dim_size(4)
    e_lo%compound_count(6)=e_lo%compound_count(5)*e_lo%dim_size(5)    

    e_lo%ig_comp=e_lo%compound_count(e_lo%ig_ord)
    e_lo%isgn_comp=e_lo%compound_count(e_lo%isgn_ord)
    e_lo%ik_comp=e_lo%compound_count(e_lo%ik_ord)
    e_lo%it_comp=e_lo%compound_count(e_lo%it_ord)
    e_lo%il_comp=e_lo%compound_count(e_lo%il_ord)
    e_lo%is_comp=e_lo%compound_count(e_lo%is_ord)

  end subroutine init_energy_layouts

  !> FIXME : Add documentation    
  elemental function is_idx_e (lo, i)
    implicit none
    integer :: is_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_e = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_e

  !> FIXME : Add documentation    
  elemental function il_idx_e (lo, i)
    implicit none
    integer :: il_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_e = 1 + mod((i - lo%llim_world)/lo%il_comp, lo%nlambda)
  end function il_idx_e

  !> FIXME : Add documentation
  elemental function isign_idx_e (lo, i)
    implicit none
    integer :: isign_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    isign_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1),lo%nsign)
  end function isign_idx_e

  !> FIXME : Add documentation    
  elemental function it_idx_e (lo, i)
    implicit none
    integer :: it_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_e = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
  end function it_idx_e

  !> FIXME : Add documentation    
  elemental function ik_idx_e (lo, i)
    implicit none
    integer :: ik_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_e = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
  end function ik_idx_e

  !> FIXME : Add documentation
  elemental function ig_idx_e (lo, i)
    implicit none
    integer :: ig_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_e = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_e

  !> FIXME : Add documentation    
  elemental function idx_e (lo, ig, isign, ik, it, il, is)
    implicit none
    integer :: idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, it, il, is
    integer, dimension(6) :: i
    i(lo%ik_ord) = ik - 1
    i(lo%it_ord) = it - 1
    i(lo%is_ord) = is - 1
    i(lo%il_ord) = il - 1
    i(lo%isgn_ord) = isign - 1
    i(lo%ig_ord) = ig + lo%ntgrid
    idx_e = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*(i(4) + lo%dim_size(4)*(i(5) + lo%dim_size(5)*i(6)))))
  end function idx_e

  !> FIXME : Add documentation  
  elemental function proc_id_e (lo, i)
    implicit none
    integer :: proc_id_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_e = i/lo%blocksize

  end function proc_id_e

  !> FIXME : Add documentation
  elemental function idx_local_e (lo, ig, isign, ik, it, il, is)
    implicit none
    logical :: idx_local_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, it, il, is
    idx_local_e = idx_local(lo, idx(lo, ig, isign, ik, it, il, is))
  end function idx_local_e

  !> FIXME : Add documentation  
  elemental function ig_local_e (lo, ig)
    implicit none
    logical :: ig_local_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig
    ig_local_e = lo%iproc == proc_id(lo, ig)
  end function ig_local_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lambda layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_lambda_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    use mp, only: iproc, nproc
    use file_utils, only: error_unit
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2
    integer :: ngroup, nprocset

    if (initialized_lambda_layouts) return
    initialized_lambda_layouts = .true.

    lz_lo%iproc = iproc
    lz_lo%ntgrid = ntgrid
    lz_lo%naky = naky
    lz_lo%ntheta0 = ntheta0
    lz_lo%negrid = negrid
    lz_lo%nspec = nspec
    lz_lo%ntgridtotal=2*ntgrid+1
    lz_lo%ng2 = ng2
    lz_lo%llim_world = 0
    lz_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*negrid*nspec - 1

    lambda_local = .false.

    if (layout == 'lxyes' .or. layout == 'lexys' .or. layout == 'lyxes') &
         call check_llocal (ntheta0, naky, nlambda, negrid, nspec)

    if (lambda_local) then

       lz_lo%groupblocksize = (naky*ntheta0*negrid*nspec-1)/nproc + 1
       
       ngroup = min (nproc, naky*ntheta0*negrid*nspec)
       lz_lo%ngroup = ngroup
       
       nprocset = nproc / lz_lo%ngroup
       lz_lo%nprocset = nprocset

       lz_lo%iset   = mod (iproc, nprocset)
       lz_lo%igroup = mod (iproc/nprocset, ngroup)

       lz_lo%llim_group = 0
       lz_lo%ulim_group = (2*ntgrid+1)*lz_lo%groupblocksize - 1
       lz_lo%gsize      = (2*ntgrid+1)*lz_lo%groupblocksize 
       
       lz_lo%nset = lz_lo%ulim_group/lz_lo%nprocset + 1
       
       lz_lo%llim_proc = lz_lo%igroup*lz_lo%gsize + lz_lo%iset*lz_lo%nset
       lz_lo%ulim_proc = min(lz_lo%ulim_group+lz_lo%igroup*lz_lo%gsize, &
            lz_lo%llim_proc + lz_lo%nset - 1)
       lz_lo%ulim_alloc = max(lz_lo%llim_proc, lz_lo%ulim_proc)
    else

       lz_lo%blocksize = lz_lo%ulim_world/nproc + 1
       lz_lo%llim_proc = lz_lo%blocksize*iproc
       lz_lo%ulim_proc &
            = min(lz_lo%ulim_world, lz_lo%llim_proc + lz_lo%blocksize - 1)
       lz_lo%ulim_alloc = max(lz_lo%llim_proc, lz_lo%ulim_proc)
    end if

    !See g_lo init
    select case (layout)
    case ('yxels')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=3,ik_ord=2,ie_ord=4,is_ord=5)
    case ('yxles')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=3,ik_ord=2,ie_ord=4,is_ord=5)
    case ('lexys')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=3,ik_ord=4,ie_ord=2,is_ord=5)
    case ('lxyes')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=2,ik_ord=3,ie_ord=4,is_ord=5)
    case ('lyxes')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=3,ik_ord=2,ie_ord=4,is_ord=5)
    case ('xyles')
       call set_dimension_order(lz_lo,ig_ord=1,it_ord=2,ik_ord=3,ie_ord=4,is_ord=5)
    end select

    lz_lo%dim_size(lz_lo%ig_ord)=lz_lo%ntgridtotal
    lz_lo%dim_size(lz_lo%it_ord)=ntheta0
    lz_lo%dim_size(lz_lo%ik_ord)=naky
    lz_lo%dim_size(lz_lo%ie_ord)=negrid
    lz_lo%dim_size(lz_lo%is_ord)=nspec
    
    lz_lo%compound_count(1)=1
    lz_lo%compound_count(2)=lz_lo%compound_count(1)*lz_lo%dim_size(1)
    lz_lo%compound_count(3)=lz_lo%compound_count(2)*lz_lo%dim_size(2)
    lz_lo%compound_count(4)=lz_lo%compound_count(3)*lz_lo%dim_size(3)
    lz_lo%compound_count(5)=lz_lo%compound_count(4)*lz_lo%dim_size(4)
    
    lz_lo%ig_comp=lz_lo%compound_count(lz_lo%ig_ord)
    lz_lo%ik_comp=lz_lo%compound_count(lz_lo%ik_ord)
    lz_lo%it_comp=lz_lo%compound_count(lz_lo%it_ord)
    lz_lo%ie_comp=lz_lo%compound_count(lz_lo%ie_ord)
    lz_lo%is_comp=lz_lo%compound_count(lz_lo%is_ord)
  end subroutine init_lambda_layouts

  !> FIXME : Add documentation    
  elemental function is_idx_lz (lo, i)
    implicit none
    integer :: is_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_lz = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_lz

  !> FIXME : Add documentation    
  elemental function ie_idx_lz (lo, i)
    implicit none
    integer :: ie_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_lz = 1 + mod((i - lo%llim_world)/lo%ie_comp, lo%negrid)
  end function ie_idx_lz

  !> FIXME : Add documentation    
  elemental function it_idx_lz (lo, i)
    implicit none
    integer :: it_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_lz = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
  end function it_idx_lz

  !> FIXME : Add documentation    
  elemental function ik_idx_lz (lo, i)
    implicit none
    integer :: ik_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_lz = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
  end function ik_idx_lz

  !> FIXME : Add documentation  
  elemental function ig_idx_lz (lo, i)
    implicit none
    integer :: ig_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_lz = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_lz

  !> FIXME : Add documentation    
  elemental function idx_lz (lo, ig, ik, it, ie, is)
    implicit none
    integer :: idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is
    integer, dimension(5) :: i
    i(lo%ik_ord) = ik - 1
    i(lo%it_ord) = it - 1
    i(lo%ie_ord) = ie - 1
    i(lo%is_ord) = is - 1
    i(lo%ig_ord) = ig + lo%ntgrid
    idx_lz = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*(i(4) + lo%dim_size(4)*i(5))))
  end function idx_lz

  !> FIXME : Add documentation  
  elemental function proc_id_lz (lo, i)
    implicit none
    integer :: proc_id_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    if (lambda_local) then
       proc_id_lz = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       proc_id_lz = i/lo%blocksize
    end if

  end function proc_id_lz

  !> FIXME : Add documentation    
  elemental function idx_local_lz (lo, ig, ik, it, ie, is)
    implicit none
    logical :: idx_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is
    idx_local_lz = idx_local(lo, idx(lo, ig, ik, it, ie, is))
  end function idx_local_lz

  !> FIXME : Add documentation  
  elemental function ig_local_lz (lo, ig)
    implicit none
    logical :: ig_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig
    ig_local_lz = lo%iproc == proc_id(lo, ig)
  end function ig_local_lz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! layouts for parity diagnostic 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_parity_layouts &
       (naky, nlambda, negrid, nspec)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: naky, nlambda, negrid, nspec

    if (initialized_parity_layouts) return
    initialized_parity_layouts = .true.
    
    p_lo%iproc = iproc
    p_lo%naky = naky
    p_lo%nlambda = nlambda
    p_lo%negrid = negrid
    p_lo%nspec = nspec
    p_lo%llim_world = 0
    p_lo%ulim_world = naky*negrid*nlambda*nspec - 1
    p_lo%blocksize = p_lo%ulim_world/nproc + 1
    p_lo%llim_proc = p_lo%blocksize*iproc
    p_lo%ulim_proc = min(p_lo%ulim_world, p_lo%llim_proc + p_lo%blocksize - 1)
    p_lo%ulim_alloc = max(p_lo%llim_proc, p_lo%ulim_proc)

    !See g_lo init
    select case (layout)
    case ('yxels')
       p_lo%ie_ord=2
       p_lo%il_ord=3
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*negrid
       p_lo%compound_count(4)=p_lo%compound_count(3)*nlambda
    case ('yxles')
       p_lo%ie_ord=3
       p_lo%il_ord=2
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*nlambda
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('lexys')
       p_lo%ie_ord=2
       p_lo%il_ord=1
       p_lo%ik_ord=3
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*negrid
       p_lo%compound_count(4)=p_lo%compound_count(3)*naky
    case ('lxyes')
       p_lo%ie_ord=3
       p_lo%il_ord=1
       p_lo%ik_ord=2
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*naky
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('lyxes')
       p_lo%ie_ord=3
       p_lo%il_ord=1
       p_lo%ik_ord=2
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*naky
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('xyles')
       p_lo%ie_ord=3
       p_lo%il_ord=2
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*nlambda
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    end select
    p_lo%ie_comp=p_lo%compound_count(p_lo%ie_ord)
    p_lo%ik_comp=p_lo%compound_count(p_lo%ik_ord)
    p_lo%il_comp=p_lo%compound_count(p_lo%il_ord)
    p_lo%is_comp=p_lo%compound_count(p_lo%is_ord)

  end subroutine init_parity_layouts

  !> FIXME : Add documentation  
  elemental function idx_parity (lo, ik, il, ie, is)
    implicit none
    integer :: idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, il, ie, is
    select case (layout)
    case ('yxels')
       idx_parity = ik-1 + lo%naky*(ie-1 + lo%negrid*(il-1 + lo%nlambda*(is-1)))
    case ('yxles')
       idx_parity = ik-1 + lo%naky*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))
    case ('lexys')
       idx_parity = il-1 + lo%nlambda*(ie-1 + lo%negrid*(ik-1 + lo%naky*(is-1)))
    case ('lxyes')
       idx_parity = il-1 + lo%nlambda*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
    case ('lyxes')
       idx_parity = il-1 + lo%nlambda*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
    case ('xyles')
       idx_parity = ik-1 + lo%naky*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))
    end select
  end function idx_parity

  !> FIXME : Add documentation  
  elemental function idx_local_parity (lo, ik, il, ie, is)
    implicit none
    logical :: idx_local_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, il, ie, is

    idx_local_parity = idx_local(lo, idx(lo, ik, il, ie, is))
  end function idx_local_parity

  !> FIXME : Add documentation  
  elemental function ig_local_parity (lo, ip)
    implicit none
    logical :: ig_local_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ip

    ig_local_parity = lo%iproc == proc_id(lo, ip)
  end function ig_local_parity

  !> FIXME : Add documentation  
  elemental function proc_id_parity (lo, i)
    implicit none
    integer :: proc_id_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_parity = i/lo%blocksize
  end function proc_id_parity

  !> FIXME : Add documentation  
  elemental function ik_idx_parity (lo, i)
    implicit none
    integer :: ik_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_parity = 1+mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
  end function ik_idx_parity

  !> FIXME : Add documentation  
  elemental function is_idx_parity (lo, i)
    
    implicit none

    integer :: is_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_parity = 1+mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_parity

  !> FIXME : Add documentation  
  elemental function il_idx_parity (lo, i)

    implicit none

    integer :: il_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_parity = 1+mod((i - lo%llim_world)/lo%il_comp, lo%nlambda)
  end function il_idx_parity

  !> FIXME : Add documentation  
  elemental function ie_idx_parity (lo, i)

    implicit none

    integer :: ie_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_parity = 1+mod((i - lo%llim_world)/lo%ie_comp, lo%negrid)
  end function ie_idx_parity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> FIXME : Add documentation  
  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use mp, only: iproc, nproc, proc0, barrier
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    integer :: nprocset, ngroup, nblock, ntgridtotal, nsign
    real :: unbalanced_amount

    if (initialized_x_transform) return
    initialized_x_transform = .true.

    xxf_lo%iproc = iproc
    xxf_lo%ntgrid = ntgrid
    ntgridtotal=(2*ntgrid+1)
    xxf_lo%ntgridtotal = ntgridtotal
    nsign=2
    xxf_lo%nsign = nsign
    xxf_lo%naky = naky
    xxf_lo%ntheta0 = ntheta0
    if (nx > ntheta0) then
       xxf_lo%nx = nx
    else
       xxf_lo%nx = (3*ntheta0+1)/2
    end if
    xxf_lo%nadd = xxf_lo%nx - ntheta0
    xxf_lo%nlambda = nlambda
    xxf_lo%negrid = negrid
    xxf_lo%nspec = nspec
    xxf_lo%llim_world = 0
    xxf_lo%ulim_world = naky*(2*ntgrid+1)*2*nlambda*negrid*nspec - 1

    !See g_lo init
    select case (layout)
    case ('yxels')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,ie_ord=4,il_ord=5,is_ord=6)
    case ('yxles')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lexys')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lxyes')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lyxes')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('xyles')
       call set_dimension_order(xxf_lo,ik_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    end select

    xxf_lo%dim_size(xxf_lo%ik_ord)=naky
    xxf_lo%dim_size(xxf_lo%ig_ord)=ntgridtotal
    xxf_lo%dim_size(xxf_lo%isgn_ord)=nsign
    xxf_lo%dim_size(xxf_lo%il_ord)=nlambda
    xxf_lo%dim_size(xxf_lo%ie_ord)=negrid
    xxf_lo%dim_size(xxf_lo%is_ord)=nspec
    
    xxf_lo%compound_count(1)=1
    xxf_lo%compound_count(2)=xxf_lo%compound_count(1)*xxf_lo%dim_size(1)
    xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*xxf_lo%dim_size(2)
    xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*xxf_lo%dim_size(3)
    xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*xxf_lo%dim_size(4)
    xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*xxf_lo%dim_size(5)    

    xxf_lo%ig_comp=xxf_lo%compound_count(xxf_lo%ig_ord)
    xxf_lo%isgn_comp=xxf_lo%compound_count(xxf_lo%isgn_ord)
    xxf_lo%ik_comp=xxf_lo%compound_count(xxf_lo%ik_ord)
    xxf_lo%ie_comp=xxf_lo%compound_count(xxf_lo%ie_ord)
    xxf_lo%il_comp=xxf_lo%compound_count(xxf_lo%il_ord)
    xxf_lo%is_comp=xxf_lo%compound_count(xxf_lo%is_ord)

    call check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)
    if (accel_lxyes) then  

       xxf_lo%groupblocksize = (nblock-1)/nproc + 1

       ngroup = min (nproc, nblock)
       xxf_lo%ngroup = ngroup
       
       nprocset = nproc / xxf_lo%ngroup
       xxf_lo%nprocset = nprocset

       xxf_lo%iset   = mod (iproc, nprocset)
       xxf_lo%igroup = mod (iproc/nprocset, ngroup)

       xxf_lo%llim_group = 0
       xxf_lo%ulim_group = naky*(2*ntgrid+1)*2*nlambda*xxf_lo%groupblocksize - 1
       xxf_lo%gsize      = naky*(2*ntgrid+1)*2*nlambda*xxf_lo%groupblocksize 
       
       xxf_lo%nset = xxf_lo%ulim_group/xxf_lo%nprocset + 1
       
       xxf_lo%llim_proc = xxf_lo%igroup*xxf_lo%gsize + xxf_lo%iset*xxf_lo%nset
       xxf_lo%ulim_proc = min(xxf_lo%ulim_group+xxf_lo%igroup*xxf_lo%gsize, &
            xxf_lo%llim_proc + xxf_lo%nset - 1)
       xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)

    else
       ! AJ November 2011
       ! unbalanced_xxf is a variable initialised in the input file
       ! which if set to .true. will enable the code below to
       ! investigate whether an unbalanced decomposition can be
       ! constructed.  Whether the unbalanced decomposition is used
       ! or not is also dependent on the variable max_unbalanced_xxf
       ! which is also set in the input file and defines the maximum
       ! amount of imbalance permitted.
       ! The code that constructs the unbalance decomposition works
       ! through the different indicies that are used in the xxf_lo
       ! to construct the full xxf_lo data space, namely:
       ! naky*(2*ntgrid+1)*nsign*nlambda*negrid*nspec
       ! Note: We precalculate 2*ntgrid+1 as ntgridtotal.
       if (unbalanced_xxf) then
          call calculate_unbalanced_x(nproc, iproc, unbalanced_amount)
       end if

       ! If we are not using the unbalanced code, either because the
       ! user has chosen not to do this in the code or because the
       ! calculated imbalance was too large then use the original
       ! decomposition.
       if (.not. unbalanced_xxf) then

          !write (*,*) 'XXF_LO%ulim_world', xxf_lo%ulim_world, nproc
          xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
          xxf_lo%llim_proc = xxf_lo%blocksize*iproc
          xxf_lo%ulim_proc &
               = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
          xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)

       else
             
          if(proc0) then
             write(*, fmt="('Using unbalanced decomposition for xxf. ',&
              & 'Unbalanced fraction',F6.2)") unbalanced_amount
          end if

       end if

    end if

  end subroutine init_x_transform_layouts

  !> This subroutine (calculate_unbalanced_x) wraps up all the 
  !! functionality required to calculate the xxf blocksizes and 
  !! work out if and unbalanced decomposition should be used, and if so 
  !! the size of the unbalanced blocks.
  !!
  !! This routine takes the number of processes to be used and the rank 
  !! of the calling process and returns the amount of unbalance in the 
  !! calculate decomposition (in a real variable with a range between 
  !! 0.0 and 1.0 where 1.0 is a 100% difference between the small and 
  !! large block size and 0.0 is no difference between the small and 
  !! large block sizes).
  !!
  !! This routine can be called from routines within the gs2_layouts
  !! module (specifically init_x_transform_layouts) or from outside 
  !! this module (specifically it is currently called from ingen to 
  !! provide suggestions for unbalanced decompositions and process 
  !! counts).  However, this routine relies on init_x_transform_layouts 
  !! having been called to initialise some of the values it used 
  !! (those in xxf_lo%) which is why there is a check whether the 
  !! initialize_x_transform variable has been set to true.  If not this 
  !! routine will not run and will notify the user of the problem.
  !!
  !! AJ, November 2011: New code from DCSE project  
  subroutine calculate_unbalanced_x(nproc, iproc, unbalanced_amount)
    implicit none
    integer, intent (in) :: nproc, iproc
    real, intent (out) :: unbalanced_amount
    integer :: level_proc_num, i, j, k, m

    ! Ensure that the xxf_lo% data has been properly initialized as this 
    ! routine relies on some data from that data structure.  If it has not 
    ! then abort this routine.
    if(.not. initialized_x_transform) then
       write(*,*) 'X Transform data structures not initialized so calculate_unbalanced_x subroutine will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_unbalanced_x'
       return
    end if
    
    ! The xxf_lo%blocksize is initialised to the value that
    ! is used if this unbalanced code is not used (the
    ! original code).
    xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
    xxf_lo%small_block_balance_factor = 1
    xxf_lo%large_block_balance_factor = 1                                                        
    
    ! Initialise the number of processors/cores used in the
    ! decomposition as the total number available.
    level_proc_num = nproc
    
    ! The first factor to be considered is nspec (number
    ! of species).
    k = xxf_lo%nspec
    
    ! Calculate the factors of nspec and the number of
    ! processors/cores we have left from the decomposition.
    ! Further details of i, m, and j are provided in the
    ! comments for the subroutine calculate_block_breakdown
    ! but briefly described i is the remainder left from k/level_proc_num,
    ! m is the remainder of level_proc_num/k and j is the integer
    ! result of k/level_proc_num.
    call calculate_block_breakdown(k, i, m, j, level_proc_num)
    
    ! If j = 0 this means that level_proc_num is larger than k
    ! the factor we are looking to divide so we have reached the
    ! point where we can't divide anymore and can now build the
    ! decomposition.
    if(j .eq. 0) then
       
       ! If we have got to this point then k is larger than the
       ! number of processes we have so we can try and split it
       ! up over the processes.
       ! If m = 0 then we know that k exactly divides into the
       ! number of processes we have at this point.  This means
       ! we can eliminate this factor from the decomposition
       ! and reduce the number of processes we have available
       ! by that factor.   Therefore we divide the number
       ! of processes we have by k and set the factor we need
       ! to decompose to 1.
       ! If m does not equal zero then we keep the current number
       ! of processes and pass the factor (k) on to the next
       ! level of the decomposition.
       if(m .eq. 0) then
          level_proc_num = level_proc_num/k
          k = 1
       end if
       
       ! For the xxf decomposition the next factor to be used
       ! depends on the layout chosen.  If it is the yxels
       ! decompsoition then the next factor is nlambda, otherwise
       ! it is negrid.
       ! Multiple the previous factor by the next factor to be
       ! be used.  If the previous factor divided exactly by the
       ! processes available then k will be 1 so we will only be
       ! considering this factor.
       select case(layout)
       case ('yxels')
          k = xxf_lo%nlambda * k
       case default
          k = xxf_lo%negrid * k
       end select
       
       ! The next code follows the pattern described above
       ! moving through all the other factor, i.e.:
       ! for yxels: negrid*nsign*ntgridtotal*naky
       ! and for the other layouts: nlambda*nsign*ntgridtotal*naky
       ! until it gets to a point where j = 0 (i.e. level_proc_num
       ! is greater than the factors being decomposed), or we
       ! have run out of factors to be decomposed.
       ! Once either of those points is reached then the
       ! code moves on to calculate the blocksizes in the
       ! else clauses of this code below.
       ! This part of the code is commented in the first occurence
       ! of "if(i .ne. 0) then" below.
       call calculate_block_breakdown(k, i, m, j, level_proc_num)
       
       if(j .eq. 0) then
          
          if(m .eq. 0) then
             level_proc_num = level_proc_num/k
             k = 1
          end if
          
          select case(layout)
          case ('yxels')
             k = xxf_lo%negrid * k
          case default
             k = xxf_lo%nlambda * k
          end select
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = xxf_lo%nsign * k
             call calculate_block_breakdown(k, i, m, j, level_proc_num)
             
             if(j .eq. 0) then
                
                if(m .eq. 0) then
                   level_proc_num = level_proc_num/k
                   k = 1
                end if
                
                k = xxf_lo%ntgridtotal * k
                call calculate_block_breakdown(k, i, m, j, level_proc_num)
                
                if(j .eq. 0) then
                   
                   if(m .eq. 0) then
                      level_proc_num = level_proc_num/k
                      k = 1
                   end if
                   
                   k = xxf_lo%naky * k
                   call calculate_block_breakdown(k, i, m, j, level_proc_num)
                   
                   ! At this point we have run out of factors to
                   ! decompose.  Now we check whether i is not
                   ! equal to 0.  If i is equal to 0 (remember that
                   ! i is the remainder left when k/level_proc_num)
                   ! then this is not an unbalanced decomposition so
                   ! we do not have to calculate the unbalanced decompostion.
                   ! If i is not equal to 0 then k cannot be exactly divided
                   ! by level_proc_num.  This means we cannot evenly divide
                   ! the decomposition over the processes we have available so
                   ! we need to find a way to create a different decomposition.
                   if(i .ne. 0) then
                      
                      ! Calculate the least unbalanced split of k over level_proc_num and also
                      ! how what factor of processes are assigned to each block size.
                      call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
                      ! Calculate the actual block sizes using the factors calculated above and the remaining data space to be
                      ! decomposed.  In this instance we are at the lowest level of the decomposition so there is nothing left to decompose
                      ! so the remaining data space is 1.  For other levels of the decomposition this 1 is replaced by the part
                      ! of the data space that has not been split up yet.
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                           1, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)     
                      ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                      ! naky for this level of the decomposition.
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                           xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple,  xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)                         
                   ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                   ! naky,ntgridtotal for this level of the decomposition.
                   call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                        xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                        xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign for this level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                     xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
             
             select case(layout)
             case('yxels')
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign,negrid for this layout and level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%negrid*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, &
                     xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
             case default
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign,nlambda for these layouts and this level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%nlambda*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, &
                     xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
             end select
             
          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
          ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
          ! naky,ntgridtotal,nsign,nlambda,negrid for this level of the decomposition.
          call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
               xxf_lo%negrid*xxf_lo%nlambda*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, &
               xxf_lo%large_block_size, xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
          
       end if
       
    end if
    
    
    ! If the decomposition has calculated that the large block and small blocks are equal then
    ! we don't have an unbalanced decomposition so set the unbalanced_amount to 0.
    ! large_block_balance_factor and small_block_balance_factor are initialised to 1 at
    ! the start of this decomposition code so if they haven't been changed we already
    ! have a balanced decomposition.
    if(xxf_lo%large_block_balance_factor .eq. 1 .and. xxf_lo%small_block_balance_factor .eq. 1) then
       
       unbalanced_amount = 0
       
    else 
       
       ! If there is an unbalanced decomposition work out the percentage of
       ! difference between the two blocks.  This is used to ensure that
       ! we don't create decompositions that have significant differences
       ! between the two block sizes which would impact the amount of
       ! computation the different groups of processes have to perform.
       unbalanced_amount = real(xxf_lo%large_block_size)/real(xxf_lo%small_block_size)
       unbalanced_amount = unbalanced_amount - 1
       
    end if
    
    ! If we calculate that there is not an unbalanced decomposition or that the
    ! amount of unbalance is larger than a integer the user sets in the input file
    ! called: max_unbalanced_xxf ; then do not use the unbalanced decomposition.
    if (unbalanced_amount .gt. max_unbalanced_xxf .or. unbalanced_amount .eq. 0) then
       
       unbalanced_xxf = .false.       
       
    end if
    
    
  end subroutine calculate_unbalanced_x
  
  !> This subroutine (calculate_block_breakdown) is used in the code
  !! to create unbalanced xxf and yxf layouts to address MPI
  !! communication overheads in the redistribution of data for the
  !! FFTs used in the nonlinear calculations to move from Fourier
  !! space to real space and back again.
  !!
  !! The routine calculates the factor k divided by a provided number
  !! of processors/cores and returns it in integer j.  As j is an
  !! integer the number returned will be more than zero if k is larger
  !! than level_proc_num and zero if k is smaller than level_proc_num.
  !!
  !! The integer i that is returned is the remained when k is divided
  !! by level_proc_num.  This is used in the unbalanced decomposition
  !! code to work out whether the factor divides exactly by the number
  !! of procs  (in level_proc_num).  If it does then it is not an
  !! unbalanced decomposition (so the original decomposition can be
  !! used.
  !!
  !! The integer m which is returned by this subroutine is the remainder
  !! of level_proc_num divided by k.  This is used in the unbalanced
  !! decomposition code to work out if the factor k exactly divides
  !! the number of processors/cores in level_proc_num.  If it does
  !! (i.e. m is zero) then the decomposition code can perform this
  !! division and move on to the next factor.  If it does not then
  !! the code keeps the factor k and multiplies it by the next factor
  !! and checks that one (and so on).
  !!
  !! AJ, November 2011: New code from DCSE project  
  subroutine calculate_block_breakdown (k, i, m, j, level_proc_num)
    implicit none
   
    integer, intent(in) :: k, level_proc_num
    integer, intent(out) :: i, m, j

    i = mod(k, level_proc_num)
    m = mod(level_proc_num, k)
    j = k/level_proc_num

  end subroutine calculate_block_breakdown


  !> This subroutine (calculate_unbalanced_decomposition) is used in the code
  !! to create unbalanced xxf and yxf layouts to address MPI
  !! communication overheads in the redistribution of data for the
  !! FFTs used in the nonlinear calculations to move from Fourier
  !! space to real space and back again.
  !!
  !! k is an integer input which is the factor that needs to be decomposed
  !! across the number of processes (the input integer localprocs)
  !! that we have at this point in the decomposition.
  !!
  !! When this subroutine is called k is larger than
  !! localprocs (this is enforced through the if and else statements
  !! in the code that calls this routine).  It is also only called when
  !! localprocs does not completely/equally divide k (so k/localprocs
  !! will not given a round number).
  !!
  !! The subroutine returns four integers; i, j, numlarge, and numsmall
  !!
  !! j is calculated as the integer division of k/localprocs.  As this
  !! is an integer division the result is rounded down (remembering
  !! that k/localprocs will not be an equal division).  j is used by the
  !! small block factor (the relative size of the small block).
  !!
  !! i is the large block factor.  As j is the rounded down division,
  !! i = j + 1.  This gives two different block sizes (when these are
  !! later calculated).
  !!
  !! numlarge is used to work out how many processes will use the large
  !! block size in the decomposition.
  !!
  !! numsmall is used to work out how many processes will use the small
  !! block size in the decomposition.
  !!
  !! AJ, November 2011: New code from DCSE project
  subroutine calculate_unbalanced_decomposition(k, j, i, numsmall, numlarge, localprocs)
    implicit none

    integer, intent(in) :: k, localprocs
    integer, intent(out) :: i, j, numlarge, numsmall

    ! To calculate the number of processes to give the large block size
    ! we work out the remainder left from dividing the input factor
    ! but the processes we have left.  This remainder is proportionate
    ! to the difference in the two block sizes.  For example, if k = 62
    ! and localprocs = 3 the code will work out i to be 21 and j to be
    ! 20, although k/localprocs is actually equal to 20 and 2/3.
    ! The code will also work out numlarge to be 2 and numsmall to
    ! be 1, so of the total processes in the simulatino 1/3 will have
    ! the small block and 2/3 the large block.  This maps to the 2/3
    ! in the actual division as opposed to the integer division (i.e
    ! the factor we have removed and replaced with small and large
    ! sizes, so replacing 20 2/3 by 1/3 x 20 and 2/3 * 21.
    numlarge = mod(k,localprocs)
    j  = k/localprocs
    i = j + 1  
    numsmall = localprocs - numlarge

  end subroutine calculate_unbalanced_decomposition

  !> This subroutine (calculate_block_size) is used in the code
  !! to create unbalanced xxf and yxf layouts to address MPI
  !! communication overheads in the redistribution of data for the
  !! FFTs used in the nonlinear calculations to move from Fourier
  !! space to real space and back again.
  !!
  !! Input iproc is this processes MPI id.
  !!
  !! Input numsmall is the number of processes that will use the small
  !! block size (calculated in subroutine
  !! calculate_unbalanced_decomposition).
  !!
  !! Input numlarge is the number of processes that will use the large
  !! block size (calculated in suboutine
  !! calcuate_unbalanced_decomposition).
  !!
  !! Input smalldecomp is the small_block_balance_factor of the
  !! decomposition (calculated in subroutine
  !! calculate_unbalanced_decomposition where it is referred to as j.
  !!
  !! Input largedecomp is the large_block_balance_factor of the
  !! decomposition (calculated in subroutine
  !! calculate_unbalanced_decomposition where it is referred to as i.
  !!
  !! Input sizeblock is the factor that needs to be split over the
  !! cores/processors we have (i.e. everything that has not yet been
  !! decomposed).
  !!
  !! Output blocksize is the blocksize assigned to this actual process
  !! (to be used throughout the code of the decompositions).
  !!
  !! Output smallblocksize is the small block size.
  !!
  !! Output largeblocksize is the large block size.
  !!
  !! Output block_multiple is the size of the small block and large block
  !! added together.  This is used to calculate which process owns a particular
  !! data point at later points in the code (particular in the subroutine 
  !! proc_id).
  !! 
  !! Output llim is the lower limit of the data block for this process.
  !!
  !! Output ulim is the upper limit of the data block for this process.
  !!
  !! Output ulim_alloc is the upper allocation limit for this processes block.
  !!
  !! AJ, November 2011: New code from DCSE project
  subroutine calculate_block_size(iproc, numsmall, numlarge, smalldecomp, largedecomp, sizeblock, blocksize, smallblocksize, largeblocksize, block_multiple, llim, ulim, ulim_alloc)
    implicit none

    integer, intent(in) :: iproc, numsmall, numlarge, smalldecomp, largedecomp, sizeblock
    integer, intent(out) :: blocksize, smallblocksize, largeblocksize, block_multiple, llim, ulim, ulim_alloc
    integer :: modproc, procfactors

    ! Small blocksize is simply calculated by taking the factors left to be
    ! distributed and multiplying it by the small_block_balance_factor.
    smallblocksize = sizeblock * smalldecomp
    ! Likewise large blocksize is the factors multiplied by the
    ! large_block_balance_factor.
    largeblocksize = sizeblock * largedecomp

    ! The block multiple is the chunks that decomposition blocks are arranged
    ! in.  For instance if we have a decomposition with 3 blocks, one small 
    ! of size 640 elements and two large of 672 elements then block_multiple
    ! will be equal to 1*640+2*672.  This is then used in the proc_id code 
    ! to isolate which block an element id belongs to by factorising the id 
    ! with this multiple to isolate the id to a small set of points within 
    ! a single chunk (i.e. set of small and large blocks).  So, for instance, 
    ! if the id 15646 is presented to proc_id then you can work out that it is 
    ! in a large block 
    ! (by doing 15646 - ((15646/block_multiple)*block_multiple) in integer 
    ! arithmetic) which would give 1758, which is in the 3 block of that chunk 
    ! of blocks.  With this information it is possible to work out the
    ! particular proc that owns that chunk.  For more information the
    ! subroutines proc_id_xxf and proc_id_yxf
    
    block_multiple = (smallblocksize * numsmall) + (largeblocksize * numlarge)

    ! This is also used in the proc_id functionality as well as working out 
    ! whether this process has a small or large block.
    procfactors = numsmall + numlarge

    ! This is used to calculate whether this process has a small or large block.
    modproc = mod(iproc, procfactors)

    ! Set this processes blocksize depending if it is to use the smallblock
    ! or large block.
    if(modproc .lt. numsmall) then
       blocksize = smallblocksize
    else
       blocksize = largeblocksize
    end if

    ! Calculate the lower limit to the block this process owns
    llim = ((iproc / procfactors) * block_multiple)
    if(modproc .ne. 0) then
       if(modproc .lt. numsmall) then
          llim = llim + smallblocksize * (modproc)
       else
          llim = llim + (smallblocksize * numsmall) + (largeblocksize  * (modproc - numsmall))
       end if
    end if

    ! The upper block limit is the lower limit plus the blocksize
    ulim = llim + blocksize  - 1
    ! The allocation upper limit is the upper block limit unless this process 
    ! has no elements of this index (in which situation the ulim will be less
    ! than or equal to the llim.  This ensures that ulim_alloc = llim for zero 
    ! sized blocks.
    ulim_alloc = max(llim, ulim)

  end subroutine calculate_block_size

  !> This subroutine (calculate_idle_processes) is used to calculate the 
  !! difference between the number of processes used in the xxf and yxf
  !! data layouts.  This is important as it can affect the amount of 
  !! communication that the code has to undertake when moving between 
  !! linear and non-linear calculations.  
  !!
  !! This routine is used by ingen when it is suggesting optimal process
  !! counts for users to flag up when suggested process counts will 
  !! results in there being a significant difference in the processes 
  !! used in the two layouts, and therefore a significant communication 
  !! overhead moving between the two layouts.
  !!
  !! AJ, November 2011: New code from DCSE project
  subroutine calculate_idle_processes(nprocs, idle_percentage)
    implicit none

    integer, intent(in) :: nprocs
    real, intent(out) :: idle_percentage
    integer :: xxf_blocksize, yxf_blocksize
    real :: xxf_usedprocs, xxf_idleprocs
    real :: yxf_usedprocs, yxf_idleprocs
    real :: delta_idle_procs

    ! Ensure that the xxf_lo% and yxf_lo%data has been properly initialized as this 
    ! routine relies on some data from those data structures.  If it has not 
    ! then abort this routine.
    if(.not. initialized_x_transform .and. .not. initialized_y_transform) then
       write(*,*) 'X and/or Y transform data structures not initialized so calculate_idle_processes will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_idle_processes'
       return
    end if

    if(nprocs .lt. 1) then
       write(*,*) 'nprocs value in calculate_idle_processes subroutine is less than 1 which is incorrect.'
       write(*,*) 'calculate_idle_processes aborting.'
       return
    end if

    ! Calculate the standard xxf_blocksize
    xxf_blocksize = xxf_lo%ulim_world/nprocs + 1
    ! Use the blocksize calculated above to calculate how many processes the
    ! xxf space maps to using this block size
    xxf_usedprocs = (xxf_lo%ulim_world+1)/real(xxf_blocksize)  
    ! Now work out how many processes do not have any xxf data space assigned
    ! to them.  This is calculated using real arthimetic so it will also 
    ! include partial data spaces (so for instance it will calculate where 
    ! a process only has half a block assigned to it).
    xxf_idleprocs = nprocs - xxf_usedprocs
 
    ! Calculate the standard yxf_blocksize
    yxf_blocksize = yxf_lo%ulim_world/nprocs + 1
    ! Use the blocksize calculated above to calculate how many processes the
    ! yxf space maps to using this block size    
    yxf_usedprocs = (yxf_lo%ulim_world+1)/real(yxf_blocksize)
    ! Now work out how many processes do not have any yxf data space assigned
    ! to them.  This is calculated using real arthimetic so it will also 
    ! include partial data spaces (so for instance it will calculate where 
    ! a process only has half a block assigned to it).
    yxf_idleprocs = nprocs - yxf_usedprocs

    ! Calculate the difference between the idle processes in the yxf and xxf 
    ! decompositions.  A high delta_idle_procs will cause high communication
    ! costs in the transform routines.
    delta_idle_procs = abs(yxf_idleprocs - xxf_idleprocs)
 
    ! Roughly calculate the percentage of data to be transferred in the
    ! transform between the xxf and yxf data spaces using the delta_idle_procs
    ! variable calculated above.
    if ( delta_idle_procs .le. 1 ) then
       idle_percentage = 0.5d0 * delta_idle_procs
    else
       idle_percentage = (1.0d0 - 1.0d0/(2.0d0 * delta_idle_procs))
    end if

  end subroutine calculate_idle_processes

  !> Return the species index, given xxf-layout and a global index
  elemental function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_xxf = 1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
  end function is_idx_xxf

  !> Return the energy index, given xxf-layout and a global index
  elemental function ie_idx_xxf (lo, i)
    implicit none
    integer :: ie_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_xxf=1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
  end function ie_idx_xxf

  !> Return the lambda index, given xxf-layout and a global index
  elemental function il_idx_xxf (lo, i)
    implicit none
    integer :: il_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_xxf = 1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
  end function il_idx_xxf

  !> Return the sign index, given xxf-layout and a global index
  elemental function isign_idx_xxf (lo, i)
    implicit none
    integer :: isign_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_xxf = 1+mod((i-lo%llim_world)/lo%isgn_comp,lo%nsign)
  end function isign_idx_xxf

  !> Return the theta index, given xxf-layout and a global index
  elemental function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_xxf = -lo%ntgrid+mod((i-lo%llim_world)/lo%ig_comp,lo%ntgridtotal)
  end function ig_idx_xxf

  !> Return the ky index, given xxf-layout and a global index
  elemental function ik_idx_xxf (lo, i)
    implicit none
    integer :: ik_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_xxf = 1+mod((i-lo%llim_world)/lo%ik_comp,lo%naky)
  end function ik_idx_xxf

  !> Return the global index in the xxf_lo layout given the dimensional indices
  elemental function idx_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    integer :: idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is
    integer, dimension(6) :: i

    i(lo%ik_ord) = ik - 1
    i(lo%il_ord) = il - 1
    i(lo%ie_ord) = ie - 1
    i(lo%is_ord) = is - 1
    i(lo%ig_ord) = ig + lo%ntgrid
    i(lo%isgn_ord) = isign - 1
    idx_xxf = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*(i(4) + lo%dim_size(4)*(i(5) + lo%dim_size(5)*i(6)))))
  end function idx_xxf

  !> Return the processor id that has the point i in the xxf_lo layout
  elemental function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    integer :: block_offset, offset_block_number, j, k, tempi
    
    if (accel_lxyes) then
       proc_id_xxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       !AJ This code has been added to deal with the unbalanced decomposition functionality.
       !AJ If an unbalanced xxf decomposition is being used then the proc_id cannot
       !AJ use a simple lo%blocksize as there will be two separate block
       !AJ sizes used so we have to work out which block this i lives in and 
       !AJ therefore which process it belongs to.
       if (unbalanced_xxf) then
          !AJ block_offset works out how many groups of blocks this i point is
          !AJ after (the small and large blocks used in the decomposition are
          !AJ grouped together).
          block_offset = (i / lo%block_multiple)          
          !AJ j represents how many blocks are in each group of blocks
          j = lo%num_small + lo%num_large
          !AJ offset_block_number is the number of blocks up to the start of the 
          !AJ group of blocks we are considering.
          offset_block_number = block_offset * j
          !AJ tempi represents where this index is inside the group of blocks
          !AJ this i point sits.
          tempi = i - (block_offset * lo%block_multiple)
          !AJ Work through each block in the group of blocks and see if this i 
          !AJ is within that block.  If it is set the proc_id_xxf as this block 
          !AJ owner.          
           do k=1,j
             if(k .le. lo%num_small) then
                !AJ TODO: The if-else construct used below could potentially be rationalised 
                !AJ TODO: to a more efficient formula where:
                !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%small_block_size
                !AJ TODO: Although a method for selecting if tempi is in the small or large 
                !AJ TODO: blocks would have to be provided.
                if(tempi .lt. lo%small_block_size) then
                   !AJ (k -1) is the number of blocks that we have already considered 
                   !AJ within this group of blocks we have selected. 
                   proc_id_xxf =  offset_block_number + (k - 1)
                   exit 
                else
                   !AJ If the index is not in this block then reduce tempi by a small
                   !AJ block size and move on to the next block.
                   tempi = tempi - lo%small_block_size
                end if
             else
                !AJ TODO: The if-else construct used below could potentially be rationalised 
                !AJ TODO: to a more efficient formula where:
                !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%large_block_size
                !AJ TODO: Although a method for selecting if tempi is in the small or large 
                !AJ TODO: blocks would have to be provided.
                if(tempi .lt. lo%large_block_size) then
                   proc_id_xxf = offset_block_number + (k - 1)
                   exit 
                else
                   !AJ If the index is not in this block then reduce tempi by a large
                   !AJ block size and move on to the next block.
                   tempi = tempi - lo%large_block_size
                end if
             end if
          end do
       else
          !AJ This code is called if the unbalanced decomposition is not being used.
          proc_id_xxf = i/lo%blocksize
       end if
    end if
  end function proc_id_xxf

  !> Return whether a point is local to a processor given the dimensional indices
  elemental function idx_local_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    logical :: idx_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is
    idx_local_xxf = idx_local (lo, idx(lo, ig, isign, ik, il, ie, is))
  end function idx_local_xxf

  !> Return whether a point is local to a processor given the global index
  elemental function ig_local_xxf (lo, i)
    implicit none
    logical ig_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_xxf = lo%iproc == proc_id(lo, i)
  end function ig_local_xxf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Y-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_y_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use mp, only: iproc, nproc, proc0
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny, ngroup, nprocset, nblock, ntgridtotal, nsign
    real :: unbalanced_amount

    if (initialized_y_transform) return
    initialized_y_transform = .true.

    if (nx > ntheta0) then
       nnx = nx
    else
       nnx = (3*ntheta0+1)/2
    end if
    if (ny > naky) then
       nny = ny
    else
       nny = 3*naky
    end if

    yxf_lo%iproc = iproc
    yxf_lo%ntgrid = ntgrid
    ntgridtotal=(2*ntgrid+1)
    yxf_lo%ntgridtotal = ntgridtotal
    nsign=2
    yxf_lo%nsign = nsign
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nlambda = nlambda
    yxf_lo%negrid = negrid
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*2*nlambda*negrid*nspec - 1

    ! Calculate compound index divisors
    select case (layout)
    case ('yxels')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,ie_ord=4,il_ord=5,is_ord=6)
    case ('yxles')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lexys')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lxyes')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('lyxes')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    case ('xyles')
       call set_dimension_order(yxf_lo,it_ord=1,ig_ord=2,isgn_ord=3,il_ord=4,ie_ord=5,is_ord=6)
    end select

    yxf_lo%dim_size(yxf_lo%it_ord)=nx
    yxf_lo%dim_size(yxf_lo%ig_ord)=ntgridtotal
    yxf_lo%dim_size(yxf_lo%isgn_ord)=nsign
    yxf_lo%dim_size(yxf_lo%il_ord)=nlambda
    yxf_lo%dim_size(yxf_lo%ie_ord)=negrid
    yxf_lo%dim_size(yxf_lo%is_ord)=nspec
    
    yxf_lo%compound_count(1)=1
    yxf_lo%compound_count(2)=yxf_lo%compound_count(1)*yxf_lo%dim_size(1)
    yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*yxf_lo%dim_size(2)
    yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*yxf_lo%dim_size(3)
    yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*yxf_lo%dim_size(4)
    yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*yxf_lo%dim_size(5)    

    yxf_lo%ig_comp=yxf_lo%compound_count(yxf_lo%ig_ord)
    yxf_lo%isgn_comp=yxf_lo%compound_count(yxf_lo%isgn_ord)
    yxf_lo%it_comp=yxf_lo%compound_count(yxf_lo%it_ord)
    yxf_lo%ie_comp=yxf_lo%compound_count(yxf_lo%ie_ord)
    yxf_lo%il_comp=yxf_lo%compound_count(yxf_lo%il_ord)
    yxf_lo%is_comp=yxf_lo%compound_count(yxf_lo%is_ord)

    call check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)

    if (accel_lxyes) then  

       yxf_lo%groupblocksize = (nblock-1)/nproc + 1

       ngroup = min (nproc, nblock)
       yxf_lo%ngroup = ngroup
       
       nprocset = nproc / yxf_lo%ngroup
       yxf_lo%nprocset = nprocset

       yxf_lo%iset   = mod (iproc, nprocset)
       yxf_lo%igroup = mod (iproc/nprocset, ngroup)

       yxf_lo%llim_group = 0
       yxf_lo%ulim_group = nnx*(2*ntgrid+1)*2*nlambda - 1
       yxf_lo%gsize      = nnx*(2*ntgrid+1)*2*nlambda*yxf_lo%groupblocksize
       
       yxf_lo%nset = yxf_lo%ulim_group/yxf_lo%nprocset + 1
       
       yxf_lo%llim_proc = yxf_lo%igroup*yxf_lo%gsize + yxf_lo%iset*yxf_lo%nset
       yxf_lo%ulim_proc = min(yxf_lo%ulim_group+yxf_lo%igroup*yxf_lo%gsize, &
            yxf_lo%llim_proc + yxf_lo%nset - 1)
       yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)

    else
       ! AJ November 2011
       ! unbalanced_yxf is a variable initialised in the input file
       ! which if set to .true. will enable the code below to
       ! investigate whether an unbalanced decomposition can be
       ! constructed.  Whether the unbalanced decomposition is used
       ! or not is also dependent on the variable max_unbalanced_yxf
       ! which is also set in the input file and defines the maximum
       ! amount of imbalance permitted.
       ! The code that constructs the unbalance decomposition works
       ! through the different indicies that are used in the yxf_lo
       ! to construct the full xxf_lo data space, namely:
       ! nxx*(2*ntgrid+1)*nsign*nlambda*negrid*nspec
       ! Note: We precalculate 2*ntgrid+1 as ntgridtotal.
       ! This functionality is the same as the functionality in
       ! init_x_transform_layouts which is extensively commented.
       ! Please see init_x_transform_layouts for more details on the 
       ! functionality below.
       if (unbalanced_yxf) then

          call calculate_unbalanced_y(nproc, iproc, unbalanced_amount)
          
       end if
       
       if (.not. unbalanced_yxf) then
          
          yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
          yxf_lo%llim_proc = yxf_lo%blocksize*iproc
          yxf_lo%ulim_proc &
               = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
          yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
       
       else
             
          if(proc0) then
             write(*, fmt="('Using unbalanced decomposition for yxf. ',&
             & 'Unbalanced fraction',F6.2)") unbalanced_amount
          end if

       end if

    end if


  end subroutine init_y_transform_layouts

  !> FIXME : Add documentation
  subroutine calculate_unbalanced_y(nproc, iproc,  unbalanced_amount)
    implicit none
    integer, intent (in) :: nproc, iproc
    real, intent (out) :: unbalanced_amount
    integer :: level_proc_num, i, j, k, m

    if(.not. initialized_y_transform) then
       write(*,*) 'Y Transform data structures not initialized so calculate_unbalanced_y subroutine will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_unbalanced_y'
       return
    end if

    yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
    yxf_lo%small_block_balance_factor = 1
    yxf_lo%large_block_balance_factor = 1                                                        
    
    level_proc_num = nproc

    k = yxf_lo%nspec

    call calculate_block_breakdown(k, i, m, j, level_proc_num)
    
    if(j .eq. 0) then
       
       if(m .eq. 0) then
          level_proc_num = level_proc_num/k
          k = 1
       end if
       
       select case(layout)
       case ('yxels')
          k = yxf_lo%nlambda * k
       case default
          k = yxf_lo%negrid * k
       end select
       call calculate_block_breakdown(k, i, m, j, level_proc_num)
       
       if(j .eq. 0) then
          
          if(m .eq. 0) then
             level_proc_num = level_proc_num/k
             k = 1
          end if
          
          select case(layout)
          case ('yxels')
             k = yxf_lo%negrid * k
          case default
             k = yxf_lo%nlambda * k
          end select
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = yxf_lo%nsign * k
             call calculate_block_breakdown(k, i, m, j, level_proc_num)
             
             if(j .eq. 0) then
                
                if(m .eq. 0) then
                   level_proc_num = level_proc_num/k
                   k = 1
                end if
                
                k = yxf_lo%ntgridtotal * k
                call calculate_block_breakdown(k, i, m, j, level_proc_num)
                
                if(j .eq. 0) then
                   
                   if(m .eq. 0) then
                      level_proc_num = level_proc_num/k
                      k = 1
                   end if
                   
                   k = yxf_lo%nx * k
                   call calculate_block_breakdown(k, i, m, j, level_proc_num)
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                           1, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)     
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                           yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)                         
                   call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                        yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                        yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                     yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
             
             select case(layout)
             case('yxels')
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%negrid*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, &
                     yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
             case default
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%nlambda*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, &
                     yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
             end select
             
          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
          call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
               yxf_lo%negrid*yxf_lo%nlambda*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, &
               yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
          
       end if
       
    end if
    
    
    if(yxf_lo%large_block_balance_factor .eq. 1 .and. yxf_lo%small_block_balance_factor .eq. 1) then
       
       unbalanced_amount = 0
       
    else 
       
       unbalanced_amount = real(yxf_lo%large_block_size)/real(yxf_lo%small_block_size)
       unbalanced_amount = unbalanced_amount - 1
              
    end if
    
    if (unbalanced_amount .gt. max_unbalanced_yxf .or. unbalanced_amount .eq. 0) then
       
       unbalanced_yxf = .false.
       
    end if
    

  end subroutine calculate_unbalanced_y

  !> Return the species index, given yxf-layout and a global index
  elemental function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_yxf = 1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
  end function is_idx_yxf

  !> Return the energy index, given yxf-layout and a global index
  elemental function ie_idx_yxf (lo, i)
    implicit none
    integer :: ie_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_yxf = 1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
  end function ie_idx_yxf

  !> Return the lambda index, given yxf-layout and a global index
  elemental function il_idx_yxf (lo, i)
    implicit none
    integer :: il_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_yxf = 1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
  end function il_idx_yxf

  !> Return the sign index, given yxf-layout and a global index
  elemental function isign_idx_yxf (lo, i)
    implicit none
    integer :: isign_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_yxf = 1+mod((i-lo%llim_world)/lo%isgn_comp,lo%nsign)
  end function isign_idx_yxf
 
  !> Return the theta index, given yxf-layout and a global index of a point on this proc
  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_yxf = -lo%ntgrid+mod((i-lo%llim_world)/lo%ig_comp,lo%ntgridtotal)
  end function ig_idx_yxf

  !> Return the kx index, given yxf-layout and a global index
  elemental function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_yxf = 1+mod((i-lo%llim_world)/lo%it_comp,lo%nx)
  end function it_idx_yxf

  !> FIXME : Add documentation  
  elemental function idx_accelx (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is

    select case (layout)
    case ('yxels')
       idx_accelx = ik-1 + lo%ny*(it-1 + lo%nx*(ie-1 + &
            lo%negrid*(il-1 + lo%nlambda*(is-1))))
    case ('yxles')
       idx_accelx = ik-1 + lo%ny*(it-1 + lo%nx*(il-1 + &
            lo%nlambda*(ie-1 + lo%negrid*(is-1))))
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       idx_accelx = it-1 + lo%nx*(ik-1 + lo%ny*(il-1 + &
!            lo%nlambda*(ie-1 + lo%negrid*(is-1))))
    end select
  end function idx_accelx

  !> Return the global index in the yxf_lo layout given the dimensional indices
  elemental function idx_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    integer :: idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is
    integer, dimension(6) :: i

    i(lo%it_ord) = it - 1
    i(lo%il_ord) = il - 1
    i(lo%ie_ord) = ie - 1
    i(lo%is_ord) = is - 1
    i(lo%ig_ord) = ig + lo%ntgrid
    i(lo%isgn_ord) = isign - 1
    idx_yxf = i(1) + lo%dim_size(1)*(i(2) + lo%dim_size(2)*(i(3) + lo%dim_size(3)*(i(4) + lo%dim_size(4)*(i(5) + lo%dim_size(5)*i(6)))))
  end function idx_yxf

  elemental function proc_id_accelx (lo, i)
#ifdef SHMEM
    use mp, only : nproc
#endif
    implicit none
    integer :: proc_id_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

#ifdef SHMEM
    integer k 

    proc_id_accelx = nproc - 1
    do k =1, nproc-1
       if ( i < lo%proc_map(k)) then 
          proc_id_accelx = k - 1
          exit
       endif
    enddo
#else
    proc_id_accelx = i/lo%blocksize
#endif

  end function proc_id_accelx

  !> Return the processor id that has the point i in the yxf_lo layout
  elemental function proc_id_yxf (lo, i)
    implicit none
    integer :: proc_id_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    integer :: block_offset, offset_block_number, j, k, tempi

    if (accel_lxyes) then
       proc_id_yxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       !AJ This code has been added to deal with the unbalanced decomposition functionality.
       !AJ If an unbalanced yxf decomposition is being used then the proc_id cannot
       !AJ use a simple lo%blocksize as there will be two separate block
       !AJ sizes used so we have to work out which block this i lives in and 
       !AJ therefore which process it belongs to.
       if (unbalanced_yxf) then
          !AJ block_offset works out how many groups of blocks this i point is
          !AJ after (the small and large blocks used in the decomposition are
          !AJ grouped together).
          block_offset = (i / lo%block_multiple)
          !AJ j represents how many blocks are in each group of blocks
          j = lo%num_small + lo%num_large
          !AJ offset_block_number is the number of blocks up to the start of the 
          !AJ group of blocks we are considering.
          offset_block_number = block_offset * j
          !AJ tempi represents where this index is inside the group of blocks
          !AJ this i point sits.
          tempi = i - (block_offset * lo%block_multiple)
          !AJ Work through each block in the group of blocks and see if this i 
          !AJ is within that block.  If it is set the proc_id_xxf as this block 
          !AJ owner.          
          do k=1,j
             if(k .le. lo%num_small) then
                !AJ TODO: The if-else construct used below could potentially be rationalised 
                !AJ TODO: to a more efficient formula where:
                !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%small_block_size
                !AJ TODO: Although a method for selecting if tempi is in the small or large 
                !AJ TODO: blocks would have to be provided.
                if(tempi .lt. lo%small_block_size) then
                   !AJ (k -1) is the number of blocks that we have already considered 
                   !AJ within this group of blocks we have selected. 
                   proc_id_yxf = offset_block_number + (k - 1)
                   exit 
                else
                   !AJ If the index is not in this block then reduce tempi by a small
                   !AJ block size and move on to the next block.
                   tempi = tempi - lo%small_block_size
                end if
             else
                !AJ TODO: The if-else construct used below could potentially be rationalised 
                !AJ TODO: to a more efficient formula where:
                !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%large_block_size
                !AJ TODO: Although a method for selecting if tempi is in the small or large 
                !AJ TODO: blocks would have to be provided.
                if(tempi .lt. lo%large_block_size) then
                   proc_id_yxf = (block_offset * j) + (k - 1)
                   exit 
                else
                   !AJ If the index is not in this block then reduce tempi by a large
                   !AJ block size and move on to the next block.
                   tempi = tempi - lo%large_block_size
                end if
             end if
          end do
       else
          proc_id_yxf = i/lo%blocksize
       end if
    end if
  end function proc_id_yxf

  !> FIXME : Add documentation  
  elemental function idx_local_accelx (lo, ik, it, il, ie, is)
    implicit none
    logical :: idx_local_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    idx_local_accelx = idx_local (lo, idx(lo, ik, it, il, ie, is))
  end function idx_local_accelx

  !> FIXME : Add documentation  
  elemental function ig_local_accelx (lo, i)
    implicit none
    logical ig_local_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i
#ifndef SHMEM
    ig_local_accelx = lo%iproc == proc_id(lo, i)
#else
    ig_local_accelx = .false.
    if ( i >= lo%llim_proc .and. i <= lo%ulim_proc) &
         ig_local_accelx = .true.
#endif
  end function ig_local_accelx

  !> Return whether a point is local to a processor given the dimensional indices
  elemental function idx_local_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    logical :: idx_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is
    idx_local_yxf = idx_local (lo, idx(lo, ig, isign, it, il, ie, is))
  end function idx_local_yxf

  !> Return whether a point is local to a processor given the global index
  elemental function ig_local_yxf (lo, i)
    implicit none
    logical ig_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_yxf = lo%iproc == proc_id(lo, i)
  end function ig_local_yxf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Accelerated FFT layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation  
  subroutine init_accel_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use mp, only: iproc, nproc
#ifdef SHMEM
    use mp, only: proc0, allgather
    use shm_mpi3, only : shm_info
#endif
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny
#ifdef SHMEM
    integer :: nxy, les, xyb, llim_nd, nd_les, ind, left_lim_nd, sfft
#endif

    if (initialized_accel_transform_layouts) return
    initialized_accel_transform_layouts = .true.

    if (nx > ntheta0) then
       nnx = nx
    else
       nnx = (3*ntheta0+1)/2
    end if
    if (ny > naky) then
       nny = ny
    else
       nny = 3*naky
    end if

    accelx_lo%iproc = iproc
    accelx_lo%ntgrid = ntgrid
    accelx_lo%nsign = 2
    accelx_lo%naky = naky
    accelx_lo%ny = nny
    accelx_lo%ntheta0 = ntheta0
    accelx_lo%nx = nnx
    accelx_lo%nxny = nnx*nny
    accelx_lo%nlambda = nlambda
    accelx_lo%negrid = negrid
    accelx_lo%nspec = nspec
    accelx_lo%llim_world = 0
    accelx_lo%ulim_world = nnx*nny*nlambda*negrid*nspec - 1
#ifndef SHMEM
    accelx_lo%blocksize = accelx_lo%ulim_world/nproc + 1
    accelx_lo%llim_proc = accelx_lo%blocksize*iproc
    accelx_lo%ulim_proc &
         = min(accelx_lo%ulim_world, accelx_lo%llim_proc + accelx_lo%blocksize - 1)
    accelx_lo%ulim_alloc = max(accelx_lo%llim_proc, accelx_lo%ulim_proc)

    accel_lo%iproc = iproc
    accel_lo%ntgrid = ntgrid
    accel_lo%nsign = 2
    accel_lo%naky = naky
    accel_lo%ndky = nny/2+1  ! Note: this is for dealiased k space quantities
    accel_lo%ny = nny
    accel_lo%ntheta0 = ntheta0
    accel_lo%nx = nnx
    accel_lo%nxny = nnx*nny
    accel_lo%nxnky = nnx*(nny/2+1)
    accel_lo%nlambda = nlambda
    accel_lo%negrid = negrid
    accel_lo%nspec = nspec
    accel_lo%nia = negrid*nlambda*nspec/nproc
    accel_lo%llim_world = 0
    accel_lo%ulim_world = nnx*(nny/2+1)*nlambda*negrid*nspec - 1
    accel_lo%blocksize = accel_lo%ulim_world/nproc + 1
    accel_lo%llim_proc = accel_lo%blocksize*iproc
    accel_lo%ulim_proc &
         = min(accel_lo%ulim_world, accel_lo%llim_proc + accel_lo%blocksize - 1)
    accel_lo%ulim_alloc = max(accel_lo%llim_proc, accel_lo%ulim_proc)

#else

    nxy=nnx*nny
    les = nlambda*negrid*nspec

    accelx_lo%ppn = shm_info%size
    accelx_lo%nnd = nproc / accelx_lo%ppn
    ! needs a test and to see if every node has the same number of processors
    !if ( mod(nproc,accelx_lo%ppn) /= 0) accelx_lo%nnd = accelx_lo%nnd + 1
    ind = iproc / accelx_lo%ppn

    accelx_lo%bpn = les / accelx_lo%nnd
    if (ind < mod(les,  accelx_lo%nnd )) accelx_lo%bpn = accelx_lo%bpn + 1
    llim_nd = ind * (les / accelx_lo%nnd) + min(ind, mod(les, accelx_lo%nnd))
    accelx_lo%llim_node = llim_nd * nxy
    accelx_lo%ulim_node = accelx_lo%llim_node + accelx_lo%bpn * nxy -1
    nd_les = accelx_lo%bpn

    xyb = (nxy *nd_les) / accelx_lo%ppn
    accelx_lo%blocksize = xyb
    if ( shm_info%id < mod(nxy * nd_les, accelx_lo%ppn) ) accelx_lo%blocksize = accelx_lo%blocksize + 1
    
    !accelx_lo%llim_proc = accelx_lo%blocksize*iproc
    accelx_lo%llim_proc = accelx_lo%llim_node &
         + shm_info%id * xyb +  min(shm_info%id, mod(nxy * nd_les, accelx_lo%ppn))
    accelx_lo%ulim_proc &
         = min(accelx_lo%ulim_world, accelx_lo%llim_proc + accelx_lo%blocksize - 1)
    accelx_lo%ulim_alloc = max(accelx_lo%llim_proc, accelx_lo%ulim_proc)
    
    allocate(accelx_lo%proc_map(0:nproc-1))
    call allgather((/accelx_lo%llim_proc/), 1, accelx_lo%proc_map, 1)

    if (proc0) then 
       write(*,'(/,a)') "accelx_lo shared memory parammeters (rank 0)" 
       write(*,*) "ranks per node :", accelx_lo%ppn
       write(*,*) "total number of blocks: ", les, "; and per node :", accelx_lo%bpn
       write(*,*) "block size per rank :", accelx_lo%blocksize
    endif

    accel_lo%iproc = iproc
    accel_lo%ntgrid = ntgrid
    accel_lo%nsign = 2
    accel_lo%naky = naky
    accel_lo%ndky = nny/2+1  ! Note: this is for dealiased k space quantities
    accel_lo%ny = nny
    accel_lo%ntheta0 = ntheta0
    accel_lo%nx = nnx
    accel_lo%nxny = nnx*nny
    accel_lo%nxnky = nnx*(nny/2+1)
    accel_lo%nlambda = nlambda
    accel_lo%negrid = negrid
    accel_lo%nspec = nspec
    !accel_lo%nia = negrid*nlambda*nspec/nproc
    accel_lo%llim_world = 0
    accel_lo%ulim_world = nnx*(nny/2+1)*nlambda*negrid*nspec - 1
    !accel_lo%blocksize = accel_lo%ulim_world/nproc + 1
    nxy = nnx*(nny/2+1)

    accel_lo%ppn = shm_info%size
    accel_lo%nnd = nproc / accel_lo%ppn
    ! needs a test and to see if every node has the same number of processors
    !if ( mod(nproc, accel_lo%ppn) /= 0) accel_lo%nnd = accel_lo%nnd + 1
    ind = iproc/accel_lo%ppn

    accel_lo%bpn = les / accel_lo%nnd
    if (ind < mod(les,  accel_lo%nnd )) accel_lo%bpn = accel_lo%bpn + 1
    llim_nd = ind * (les / accel_lo%nnd) + min(ind, mod(les, accel_lo%nnd))
    accel_lo%llim_node = llim_nd * nxy
    accel_lo%ulim_node = accel_lo%llim_node + accel_lo%bpn * nxy -1 
    nd_les = accel_lo%bpn

    xyb = (nxy * nd_les) / accel_lo%ppn
    accel_lo%blocksize = xyb
    if ( shm_info%id < mod(nxy * nd_les, accel_lo%ppn) ) accel_lo%blocksize = accel_lo%blocksize + 1

    accel_lo%nia = accel_lo%bpn
    !accel_lo%llim_proc = accel_lo%blocksize*iproc
    accel_lo%llim_proc =  accel_lo%llim_node &
         + shm_info%id * xyb +  min(shm_info%id, mod(nxy * nd_les, accel_lo%ppn))
    accel_lo%ulim_proc &
         = min(accel_lo%ulim_world, accel_lo%llim_proc + accel_lo%blocksize - 1)
    accel_lo%ulim_alloc = max(accel_lo%llim_proc, accel_lo%ulim_proc)

   if (proc0) then 
       write(*,'(/,a)') "accel_lo shared memory parammeters (rank 0) " 
       write(*,*) "ranks per node :", accel_lo%ppn
       write(*,*) "total number of blocks:", les, "; per node :", accel_lo%bpn
       write(*,*) "block size per rank :", accel_lo%blocksize
    endif

#endif
  end subroutine init_accel_transform_layouts

  !> FIXME : Add documentation
  elemental function is_idx_accelx (lo, i)
    implicit none
    integer :: is_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%negrid/lo%nlambda, lo%nspec)
  end function is_idx_accelx

  !> FIXME : Add documentation  
  elemental function ie_idx_accelx (lo, i)
    implicit none
    integer :: ie_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%negrid)
    case ('yxles')
       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%nlambda, lo%negrid)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%nlambda, lo%negrid)
    end select
  end function ie_idx_accelx

  !> FIXME : Add documentation  
  elemental function il_idx_accelx (lo, i)
    implicit none
    integer :: il_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%negrid, lo%nlambda)
    case ('yxles')
       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%nlambda)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!   case ('xyles')
!       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%nlambda)
    end select
  end function il_idx_accelx

  !> FIXME : Add documentation  
  elemental function it_idx_accelx (lo, i)
    implicit none
    integer :: it_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
!    select case (layout)
!    case ('yxels')
!       it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
!    case ('yxles')
!       it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       it_idx_accelx = 1 + mod((i - lo%llim_world), lo%nx)
!    end select
  end function it_idx_accelx

  !> FIXME : Add documentation  
  elemental function ik_idx_accelx (lo, i)
    implicit none
    integer :: ik_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
!    select case (layout)
!    case ('yxels')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
!    case ('yxles')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nx, lo%ny)
!    end select
  end function ik_idx_accelx

  !> FIXME : Add documentation  
  elemental function it_idx_accel (lo, i)
    implicit none
    integer :: it_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
!    select case (layout)
!    case ('yxels')
!       it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
!    case ('yxles')
!       it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       it_idx_accel = 1 + mod((i - lo%llim_world), lo%nx)
!    end select
  end function it_idx_accel

  !> FIXME : Add documentation
  elemental function ik_idx_accel (lo, i)
    implicit none
    integer :: ik_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
!    select case (layout)
!    case ('yxels')
!       ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
!    case ('yxles')
!       ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ik_idx_accel = 1 + mod((i - lo%llim_world)/lo%nx, lo%ndky)
!    end select
  end function ik_idx_accel

  !> FIXME : Add documentation  
  elemental function dealiasing (lo, ia)
    implicit none
    logical dealiasing
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: ia
    
    dealiasing = .true.

    if (it_idx(accel_lo, ia) > lo%ntheta0/2+1 &
  .and. it_idx(accel_lo, ia) <= lo%nx-lo%ntheta0/2) return
    if (ik_idx(accel_lo, ia) > lo%naky) return

    dealiasing = .false.
  end function dealiasing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transformation subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation    
  elemental subroutine gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
    implicit none
    integer, intent (in) :: ig, isign, iglo
    type (g_layout_type), intent (in) :: g_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: it, ixxf

    it = it_idx(g_lo,iglo)
    if (it > (xxf_lo%ntheta0+1)/2) then
       it = it - xxf_lo%ntheta0 + xxf_lo%nx
    end if

    ixxf = idx(xxf_lo, ig, isign, ik_idx(g_lo,iglo), il_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2xxfidx
  
  !> FIXME : Add documentation  
  elemental subroutine xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (yxf_layout_type), intent (in) :: yxf_lo
    integer, intent (out) :: ik, iyxf

    ik = ik_idx(xxf_lo,ixxf)
    iyxf = idx(yxf_lo, ig_idx(xxf_lo,ixxf), isign_idx(xxf_lo,ixxf), &
         it, il_idx(xxf_lo,ixxf), ie_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2yxfidx

  !> FIXME : Add documentation  
  subroutine pe_layout (char)
    implicit none
    character (1), intent (out) :: char

    select case (layout)
    case ('yxels')
       char = 'v'
    case ('yxles')
       char = 'v'
    case ('lexys')
       char = 'x'
    case ('lxyes')
       char = 'm'    ! mixed
    case ('lyxes')
       char = 'b'    ! big, since this layout makes sense mainly for big runs?
    case ('xyles')
!CMR: set char="?"
! char only acts if set to "v", to force use of accel layout.
! accel layouts NOT implemented for "xyles", which assume throughout that 
! y is fastest index, including in calls to 2D FFTs
       char = '?'    
    end select

  end subroutine pe_layout

  !> Get the factors of an integer n
  subroutine factors (n, j, div)
    !> The number to factorize
    integer, intent (in) :: n
    !> The number of factors
    integer, intent (out) :: j
    !> Array containing the j factors of n
    integer, dimension (:), intent (out) :: div
    integer :: i, imax

    ! find: i = lowest factor of n
    ! and therefore imax=n/i is the HIGHEST factor of n
    do i=2,n
       if (mod(n,i)==0) exit
    end do
    imax = n/i
    j=1
    ! loop over all possible factors of n, and return in div(j)
    do i=1,imax
       if (mod(n,i)==0) then
          div(j) = i
          j=j+1
       end if
    end do
    div(j) = n
  end subroutine factors

  !> Sets the order of x y l e s in g_lo
  subroutine set_dimension_order_g(lo, it_ord, ik_ord, il_ord, ie_ord, is_ord)
    implicit none
    type (g_layout_type), intent(in out) :: lo
    integer, intent(in) :: ie_ord, il_ord, it_ord, ik_ord, is_ord
    integer, dimension(5) :: check

    ! Check input is the integers 1 to 5
    check = [it_ord, ik_ord, il_ord, ie_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (g_lo) is not unique integers 1 to 5')

    lo%ie_ord=ie_ord
    lo%il_ord=il_ord
    lo%it_ord=it_ord
    lo%ik_ord=ik_ord
    lo%is_ord=is_ord
  end subroutine set_dimension_order_g

  !> Sets the order of x y th s in le_lo
  subroutine set_dimension_order_le(lo, ig_ord, it_ord, ik_ord, is_ord)
    implicit none
    type (le_layout_type), intent(in out) :: lo
    integer, intent(in) :: ig_ord, it_ord, ik_ord, is_ord
    integer, dimension(4) :: check

    ! Check input is the integers 1 to N
    check = [ig_ord, it_ord, ik_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (le_lo) is not unique integers 1 to 4')

    lo%ig_ord=ig_ord
    lo%it_ord=it_ord
    lo%ik_ord=ik_ord
    lo%is_ord=is_ord
  end subroutine set_dimension_order_le

  !> Sets the order of x y e th s in lz_lo
  subroutine set_dimension_order_lz(lo, ig_ord, it_ord, ik_ord, ie_ord, is_ord)
    implicit none
    type (lz_layout_type), intent(in out) :: lo
    integer, intent(in) :: ig_ord, it_ord, ik_ord, ie_ord, is_ord
    integer, dimension(5) :: check

    ! Check input is the integers 1 to N
    check = [ig_ord, it_ord, ik_ord, ie_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (lz_lo) is not unique integers 1 to 5')

    lo%ig_ord=ig_ord
    lo%it_ord=it_ord
    lo%ik_ord=ik_ord
    lo%ie_ord=ie_ord    
    lo%is_ord=is_ord
  end subroutine set_dimension_order_lz

  !> Sets the order of th sgn x y l s in e_lo
  subroutine set_dimension_order_e(lo, ig_ord, isgn_ord, it_ord, ik_ord, il_ord, is_ord)
    implicit none
    type (e_layout_type), intent(in out) :: lo
    integer, intent(in) :: ig_ord, isgn_ord, it_ord, ik_ord, il_ord, is_ord
    integer, dimension(6) :: check

    ! Check input is the integers 1 to N
    check = [ig_ord, isgn_ord, it_ord, ik_ord, il_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (e_lo) is not unique integers 1 to 5')

    lo%ig_ord=ig_ord
    lo%isgn_ord=isgn_ord    
    lo%it_ord=it_ord
    lo%ik_ord=ik_ord
    lo%il_ord=il_ord    
    lo%is_ord=is_ord
  end subroutine set_dimension_order_e
  
  !> Sets the order of y th sgn l e s in xxf_lo
  !> Note: usually the order is y, theta, sign, then "les" in the order they appear in
  !> the g_lo layout
  subroutine set_dimension_order_xxf(lo,ik_ord,ig_ord,isgn_ord,il_ord,ie_ord,is_ord)
    implicit none
    type (xxf_layout_type), intent(in out) :: lo
    integer, intent(in) :: ie_ord, il_ord, ig_ord, ik_ord, is_ord, isgn_ord
    integer, dimension(6) :: check

    ! Check input is the integers 1 to N
    check = [ik_ord, ig_ord, isgn_ord, il_ord, ie_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (xxf_lo) is not unique integers 1 to 6')

    xxf_lo%ik_ord=ik_ord
    xxf_lo%ig_ord=ig_ord
    xxf_lo%isgn_ord=isgn_ord
    xxf_lo%il_ord=il_ord
    xxf_lo%ie_ord=ie_ord
    xxf_lo%is_ord=is_ord
  end subroutine set_dimension_order_xxf

  !> Sets the order of x th sgn l e s in xxf_lo
  !> Note: usually the order is x, theta, sign, then "les" in the order they appear in
  !> the g_lo layout
  subroutine set_dimension_order_yxf(lo,it_ord,ig_ord,isgn_ord,il_ord,ie_ord,is_ord)
    implicit none
    type (yxf_layout_type), intent(in out) :: lo
    integer, intent(in) :: ie_ord, il_ord, ig_ord, it_ord, is_ord, isgn_ord
    integer, dimension(6) :: check

    ! Check input is the integers 1 to N
    check = [it_ord, ig_ord, isgn_ord, il_ord, ie_ord, is_ord]
    call check_unique_integers(check,'Input to set_dimension_order (yxf_lo) is not unique integers 1 to 6')

    yxf_lo%it_ord=it_ord
    yxf_lo%ig_ord=ig_ord
    yxf_lo%isgn_ord=isgn_ord
    yxf_lo%il_ord=il_ord
    yxf_lo%ie_ord=ie_ord
    yxf_lo%is_ord=is_ord
  end subroutine set_dimension_order_yxf

  !> Checks that the given array contains unique integers 1 to (array size) (in any order),
  !> aborts if not.
  subroutine check_unique_integers(array,message)

    use mp, only: mp_abort
    use sorting, only: quicksort

    implicit none
    !> Dimension of input 
    integer, dimension(:), intent(inout) :: array
    !> Error message
    character(*), intent(in) :: message
    integer, dimension(:), allocatable :: check
    integer :: n
    integer :: i

    n = size(array)
    check = [(i, i = 1, n)]

    call quicksort(n,array)
    if (any(check.ne.array)) then
      call mp_abort(message)
    end if

  end subroutine check_unique_integers

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_gs2_layouts_config(gs2_layouts_config_in)
    use mp, only: mp_abort
    type(layouts_config_type), intent(in), optional :: gs2_layouts_config_in
    if (initialized) then
       call mp_abort("Trying to set gs2_layouts_config when already initialized.", to_screen = .true.)
    end if
    if (present(gs2_layouts_config_in)) then
       layouts_config = gs2_layouts_config_in
    end if
  end subroutine set_gs2_layouts_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the layouts_knobs namelist and populates the member variables
  subroutine read_layouts_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(layouts_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: fft_measure_plan, fft_use_wisdom, gf_local_fields, intmom_sub, intspec_sub, local_field_solve, opt_local_copy, opt_redist_init, opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap, simple_gf_decomposition, unbalanced_xxf, unbalanced_yxf
    real :: max_unbalanced_xxf
    character(len = 5) :: layout
    character(len = run_name_size) :: fft_wisdom_file

    namelist /layouts_knobs/ fft_measure_plan, fft_use_wisdom, fft_wisdom_file, gf_local_fields, intmom_sub, intspec_sub, layout, local_field_solve, max_unbalanced_xxf, opt_local_copy, opt_redist_init, opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap, simple_gf_decomposition, unbalanced_xxf, unbalanced_yxf

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    fft_measure_plan = self%fft_measure_plan
    fft_use_wisdom = self%fft_use_wisdom
    fft_wisdom_file = self%fft_wisdom_file
    gf_local_fields = self%gf_local_fields
    intmom_sub = self%intmom_sub
    intspec_sub = self%intspec_sub
    layout = self%layout
    local_field_solve = self%local_field_solve
    max_unbalanced_xxf = self%max_unbalanced_xxf
    opt_local_copy = self%opt_local_copy
    opt_redist_init = self%opt_redist_init
    opt_redist_nbk = self%opt_redist_nbk
    opt_redist_persist = self%opt_redist_persist
    opt_redist_persist_overlap = self%opt_redist_persist_overlap
    simple_gf_decomposition = self%simple_gf_decomposition
    unbalanced_xxf = self%unbalanced_xxf
    unbalanced_yxf = self%unbalanced_yxf

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = layouts_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = layouts_knobs)
    endif

    ! Now copy from local variables into type members
    self%fft_measure_plan = fft_measure_plan
    self%fft_use_wisdom = fft_use_wisdom
    self%fft_wisdom_file = fft_wisdom_file
    self%gf_local_fields = gf_local_fields
    self%intmom_sub = intmom_sub
    self%intspec_sub = intspec_sub
    self%layout = layout
    self%local_field_solve = local_field_solve
    self%max_unbalanced_xxf = max_unbalanced_xxf
    self%opt_local_copy = opt_local_copy
    self%opt_redist_init = opt_redist_init
    self%opt_redist_nbk = opt_redist_nbk
    self%opt_redist_persist = opt_redist_persist
    self%opt_redist_persist_overlap = opt_redist_persist_overlap
    self%simple_gf_decomposition = simple_gf_decomposition
    self%unbalanced_xxf = unbalanced_xxf
    self%unbalanced_yxf = unbalanced_yxf

    self%exist = exist
  end subroutine read_layouts_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_layouts_config(self, unit)
    implicit none
    class(layouts_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("fft_measure_plan", self%fft_measure_plan, unit_internal)
    call self%write_key_val("fft_use_wisdom", self%fft_use_wisdom, unit_internal)
    call self%write_key_val("fft_wisdom_file", self%fft_wisdom_file, unit_internal)
    call self%write_key_val("gf_local_fields", self%gf_local_fields, unit_internal)
    call self%write_key_val("intmom_sub", self%intmom_sub, unit_internal)
    call self%write_key_val("intspec_sub", self%intspec_sub, unit_internal)
    call self%write_key_val("layout", self%layout, unit_internal)
    call self%write_key_val("local_field_solve", self%local_field_solve, unit_internal)
    call self%write_key_val("max_unbalanced_xxf", self%max_unbalanced_xxf, unit_internal)
    call self%write_key_val("opt_local_copy", self%opt_local_copy, unit_internal)
    call self%write_key_val("opt_redist_init", self%opt_redist_init, unit_internal)
    call self%write_key_val("opt_redist_nbk", self%opt_redist_nbk, unit_internal)
    call self%write_key_val("opt_redist_persist", self%opt_redist_persist, unit_internal)
    call self%write_key_val("opt_redist_persist_overlap", self%opt_redist_persist_overlap, unit_internal)
    call self%write_key_val("simple_gf_decomposition", self%simple_gf_decomposition, unit_internal)
    call self%write_key_val("unbalanced_xxf", self%unbalanced_xxf, unit_internal)
    call self%write_key_val("unbalanced_yxf", self%unbalanced_yxf, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_layouts_config

  !> Resets the config object to the initial empty state
  subroutine reset_layouts_config(self)
    class(layouts_config_type), intent(in out) :: self
    type(layouts_config_type) :: empty
    select type (self)
    type is (layouts_config_type)
       self = empty
    end select
  end subroutine reset_layouts_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_layouts_config(self)
    use mp, only: broadcast
    implicit none
    class(layouts_config_type), intent(in out) :: self
    call broadcast(self%fft_measure_plan)
    call broadcast(self%fft_use_wisdom)
    call broadcast(self%fft_wisdom_file)
    call broadcast(self%gf_local_fields)
    call broadcast(self%intmom_sub)
    call broadcast(self%intspec_sub)
    call broadcast(self%layout)
    call broadcast(self%local_field_solve)
    call broadcast(self%max_unbalanced_xxf)
    call broadcast(self%opt_local_copy)
    call broadcast(self%opt_redist_init)
    call broadcast(self%opt_redist_nbk)
    call broadcast(self%opt_redist_persist)
    call broadcast(self%opt_redist_persist_overlap)
    call broadcast(self%simple_gf_decomposition)
    call broadcast(self%unbalanced_xxf)
    call broadcast(self%unbalanced_yxf)

    call broadcast(self%exist)
  end subroutine broadcast_layouts_config

  !> Gets the default name for this namelist
  function get_default_name_layouts_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_layouts_config
    get_default_name_layouts_config = "layouts_knobs"
  end function get_default_name_layouts_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_layouts_config()
    implicit none
    logical :: get_default_requires_index_layouts_config
    get_default_requires_index_layouts_config = .false.
  end function get_default_requires_index_layouts_config
end module gs2_layouts
