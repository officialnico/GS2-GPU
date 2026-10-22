!> Stage probes for the GS2 GPU export program `gs2_export`.
!>
!> The physics modules call the probe routines of this module at fixed points of the time
!> step (see gpu/spec/09_export_infra.md §5). The probes forward to procedure pointers which
!> only `gs2_export` installs, so this module depends on nothing in GS2 and can be used from
!> every level of the code without dependency cycles.
!>
!> Every call site has the form `if (probe_active) call probe_...(...)`, so when the probes
!> are not active (always the case in `bin/gs2`) the cost is one logical test per call site
!> and no data are touched.
module gs2_export_probe
  implicit none

  private

  public :: probe_active, probe_in_advance, probe_istep, probe_pass
  public :: probe_point, probe_g, probe_le, probe_real, probe_c2
  public :: probe_point_hook, probe_g_hook, probe_le_hook, probe_real_hook, probe_c2_hook

  !> Master switch, only ever set by gs2_export.
  logical :: probe_active = .false.
  !> True between the start and the end of fields::advance (set by the probes in fields.f90).
  logical :: probe_in_advance = .false.
  !> istep of the current advance, and which timeadv call (1 or 2; 0 outside timeadv).
  integer :: probe_istep = -1, probe_pass = 0

  abstract interface
     !> A probe without data: the hook dumps module state chosen by the tag.
     subroutine point_iface(tag, istep)
       implicit none
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
     end subroutine point_iface

     !> A g-layout array arr(1:2*ntgrid+1, 1:2, 1:nloc), nloc = g_lo%ulim_alloc-g_lo%llim_proc+1
     subroutine g_iface(tag, istep, arr)
       implicit none
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       complex, dimension(:,:,:), intent(in) :: arr
     end subroutine g_iface

     !> An le-layout array arr(1:nxi+1, 1:negrid+1, 1:nloc), nloc = le_lo%ulim_alloc-le_lo%llim_proc+1
     subroutine le_iface(tag, istep, arr)
       implicit none
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       complex, dimension(:,:,:), intent(in) :: arr
     end subroutine le_iface

     !> A small real vector (e.g. nonlinear_terms::max_vel_components)
     subroutine real_iface(tag, istep, vec)
       implicit none
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       real, dimension(:), intent(in) :: vec
     end subroutine real_iface

     !> Response build: am(1:R_c, 1:local columns) of class ic before inversion
     subroutine c2_iface(tag, ic, arr)
       implicit none
       character(len=*), intent(in) :: tag
       integer, intent(in) :: ic
       complex, dimension(:,:), intent(in) :: arr
     end subroutine c2_iface
  end interface

  procedure(point_iface), pointer :: probe_point_hook => null()
  procedure(g_iface), pointer :: probe_g_hook => null()
  procedure(le_iface), pointer :: probe_le_hook => null()
  procedure(real_iface), pointer :: probe_real_hook => null()
  procedure(c2_iface), pointer :: probe_c2_hook => null()

contains

  !> Probe without data
  subroutine probe_point(tag, istep)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in), optional :: istep
    if (.not. probe_active) return
    if (.not. associated(probe_point_hook)) return
    if (present(istep)) then
       call probe_point_hook(tag, istep)
    else
       call probe_point_hook(tag, probe_istep)
    end if
  end subroutine probe_point

  !> Probe passing a g-layout array (a dummy argument or local of the caller)
  subroutine probe_g(tag, arr, istep)
    implicit none
    character(len=*), intent(in) :: tag
    complex, dimension(:,:,:), intent(in) :: arr
    integer, intent(in), optional :: istep
    if (.not. probe_active) return
    if (.not. associated(probe_g_hook)) return
    if (present(istep)) then
       call probe_g_hook(tag, istep, arr)
    else
       call probe_g_hook(tag, probe_istep, arr)
    end if
  end subroutine probe_g

  !> Probe passing an le-layout array (a dummy argument or local of the caller)
  subroutine probe_le(tag, arr, istep)
    implicit none
    character(len=*), intent(in) :: tag
    complex, dimension(:,:,:), intent(in) :: arr
    integer, intent(in), optional :: istep
    if (.not. probe_active) return
    if (.not. associated(probe_le_hook)) return
    if (present(istep)) then
       call probe_le_hook(tag, istep, arr)
    else
       call probe_le_hook(tag, probe_istep, arr)
    end if
  end subroutine probe_le

  !> Probe passing a real vector
  subroutine probe_real(tag, vec, istep)
    implicit none
    character(len=*), intent(in) :: tag
    real, dimension(:), intent(in) :: vec
    integer, intent(in), optional :: istep
    if (.not. probe_active) return
    if (.not. associated(probe_real_hook)) return
    if (present(istep)) then
       call probe_real_hook(tag, istep, vec)
    else
       call probe_real_hook(tag, probe_istep, vec)
    end if
  end subroutine probe_real

  !> Probe passing a complex matrix of class ic (response build, probe P-R1)
  subroutine probe_c2(tag, ic, arr)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: ic
    complex, dimension(:,:), intent(in) :: arr
    if (.not. probe_active) return
    if (.not. associated(probe_c2_hook)) return
    call probe_c2_hook(tag, ic, arr)
  end subroutine probe_c2
end module gs2_export_probe
