!> Set the perpendicular box size and resolution for linear or nonlinear runs.
module kt_grids_box
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: read_parameters_box
  public :: finish_parameters_box
  public :: box_set_overrides
  public :: check_kt_grids_box, wnml_kt_grids_box
  public :: x0, y0, jtwist !RN> Caution: these are not broadcasted!

  public :: kt_grids_box_config_type
  public :: set_kt_grids_box_config
  public :: get_kt_grids_box_config

  real :: ly, y0, x0, rtwist, rhostar_box
  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: n0
  integer :: jtwist
  integer :: naky, ntheta0, nx, ny
  logical :: parameters_read = .false.
  logical :: initialized = .false.

  !> Used to represent the input configuration of kt_grids_box
  type, extends(abstract_config_type) :: kt_grids_box_config_type
     ! namelist : kt_grids_box_parameters
     ! indexed : false
     !> For finite magnetic shear determines the box size in the x
     !> direction according to \(L_x = L_y
     !> \textrm{jtwist}/2\pi\hat{s}\).  This also affects the number
     !> of connections at each `ky` when linked boundary conditions
     !> are selected in the [[dist_fn_knobs]] namelist.
     !>
     !> This gets a smart default Initialised to `jtwist =
     !> max(int(2*pi*shat+0.5,1))` so that \(L_x\approx L_y\). If the
     !> magnetic shear is less than around 0.16 then `jtwist` will
     !> default to the minimum allowed value of 1.
     integer :: jtwist = 1
     !> Sets the box size in the y direction. If set to 0 (the default) then
     !> we set `ly=2*pi*y0`.
     real :: ly = 0.0
     !> If set greater than zero (the default) then this sets the
     !> toroidal mode number of the first non-zero `ky` by overriding
     !> the value given for `y0` through `y0 =
     !> 1.0/(n0*rhostar_box*drhodpsi)` where `drhodpsi` is determined
     !> during geometry setup.
     integer :: n0 = 0
     !> The actual number of ky modes. For nonlinear runs it is
     !> generally recommended to use `ny` instead. If set to 0 (the
     !> default) this will be automatically set to `1 + (ny-1)/3`.  If
     !> both `ny` and `naky` are set then GS2 will check that `ny`
     !> is sufficiently high to ensure de-aliasing. It can be larger
     !> than this minimum value. Setting both values can be useful to
     !> allow values to be selected which are performant for the FFTs
     !> and provide a good range of sweetspots.
     integer :: naky = 0
     !> The actual number of theta0 modes. For nonlinear runs it is
     !> generally recommended to use `nx` instead. If set to 0 (the
     !> default) this will be automatically set to `1 + 2*(nx-1)/3`.
     !> If both `nx` and `ntheta0` are set then GS2 will check that
     !> `nx` is sufficiently high to ensure de-aliasing. It can be
     !> larger than this minimum value. Setting both values can be
     !> useful to allow values to be selected which are performant for
     !> the FFTs and provide a good range of sweetspots.
     integer :: ntheta0 = 0
     !> The number of kx points in inputs to the fft routines, and
     !> hence the number of radial points used in real space
     !> calculations. This differs from the actual number of kx points
     !> simulated in the rest of the code due to the need to prevent
     !> aliasing. The number of kx modes actually simulated
     !> (`ntheta0`) is, by default, equal to `1 + 2*(nx - 1)/3`.
     integer :: nx = 0
     !> The number of ky points in inputs to the fft routines, and
     !> hence the number of binormal points used in real space
     !> calculations. This differs from the actual number of ky points
     !> simulated in the rest of the code due to the need to prevent
     !> aliasing. The number of ky modes actually simulated
     !> (`naky`) is, by default, equal to `1 + (ny - 1)/3`.
     integer :: ny = 0
     !> The rhostar (`rhoref/Lref`) to use. Only used if `n0` also set
     !> greater than zero.  If `rhostar_box` and `n0` are greater than
     !> zero then `y0=1.0/(n0*rhostar_box*drhodpsi)`, which
     !> effectively sets the minimum non-zero `ky` used in the
     !> simulation.
     real :: rhostar_box = 0.0
     !> Expert usage only -- more documentation required.
     !>
     !> Used to control the kx spacing in simulations with effectively
     !> zero shear (\(\< 10^{-5}\)) where linked boundaries are not
     !> appropriate so periodic boundaries are used. Also only used if
     !> `x0` has been set to zero (the default). If `rtwist` is set to
     !> 0.0 (the default) then it is set to the value of
     !> `jtwist`. Effectively ends up setting the box size in the x
     !> direction, as \(L_x = L_y \textrm{rtwist}\) if `rtwist > 0 `
     !> and \(L_x = L_y / \textrm{rtwist}\) if `rtwist < 0`. See
     !> [[kt_grids_box_parameters:x0]] for an alternative.
     real :: rtwist = 0.0
     !> Controls the box length in the x direction (measured in the
     !> reference Larmour radius) if the magnetic shear is small (\(\<
     !> 10^{-5}\)). The box size in the x direction is given by \(L_x
     !> = 2\pi x_0\). See [[kt_grids_box_parameters:rtwist]] for an
     !> alternative.
     real :: x0 = 0.
     !> Controls the box length in the y direction (measured in the
     !> reference Larmour radius). The box size in the y direction is
     !> given by \(L_y = 2\pi y_0\). Note if `y0` is set negative
     !> then, it is replaced with `-1/y0` and Effectively sets the
     !> minimum wavelength captured by the box.
     real :: y0 = 2.0
   contains
     procedure, public :: read => read_kt_grids_box_config
     procedure, public :: write => write_kt_grids_box_config
     procedure, public :: reset => reset_kt_grids_box_config
     procedure, public :: broadcast => broadcast_kt_grids_box_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_box_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_box_config
  end type kt_grids_box_config_type

  type(kt_grids_box_config_type) :: kt_grids_box_config

contains

  !> FIXME : Add documentation
  subroutine read_parameters_box(kt_grids_box_config_in)
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: init_theta_grid, shat
    use constants, only: pi
    implicit none
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in
    logical :: exist

    if (parameters_read) return
    parameters_read = .true.

    call init_theta_grid

    if (present(kt_grids_box_config_in)) kt_grids_box_config = kt_grids_box_config_in

    ! Set the smart defaults
    if (.not.kt_grids_box_config%is_initialised()) then
       kt_grids_box_config%jtwist = max(int(2.0*pi*shat + 0.5),1)  ! default jtwist -- MAB
    end if

    call kt_grids_box_config%init(name = 'kt_grids_box_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    jtwist = kt_grids_box_config%jtwist
    ly = kt_grids_box_config%ly
    n0 = kt_grids_box_config%n0
    naky = kt_grids_box_config%naky
    ntheta0 = kt_grids_box_config%ntheta0
    nx = kt_grids_box_config%nx
    ny = kt_grids_box_config%ny
    rhostar_box = kt_grids_box_config%rhostar_box
    rtwist = kt_grids_box_config%rtwist
    x0 = kt_grids_box_config%x0
    y0 = kt_grids_box_config%y0

    exist = kt_grids_box_config%exist
  end subroutine read_parameters_box

  !> FIXME : Add documentation
  subroutine finish_parameters_box
    implicit none
    parameters_read = .false.
    initialized = .false.
    call kt_grids_box_config%reset()
  end subroutine finish_parameters_box

  !> FIXME : Add documentation
  subroutine init_kt_grids_box(kt_grids_box_config_in)
!CMR, 14/10/2013:
! New namelist variables: n0, rhostar_box.
! If n0 and rhostar_box defined, set ky(1) using toroidal mode number.

    use file_utils, only: error_unit
    use theta_grid, only: drhodpsi
    use constants, only: twopi
    use mp, only: mp_abort, proc0
    implicit none
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in

    if (initialized) return
    initialized = .true.

    call read_parameters_box(kt_grids_box_config_in)

    if (ny==0 .and. naky==0) call mp_abort("ERROR: ny==0 .and. naky==0", .true.)
    if (nx==0 .and. ntheta0==0) call mp_abort("ERROR: nx==0 .and. ntheta0==0", .true.)

    if (rhostar_box > 0.0 .and. n0 > 0) y0 = 1.0 / (n0 * rhostar_box * drhodpsi)

    if (y0 < 0) y0 = -1 / y0
    if (ly == 0.) ly = twopi * y0
    if (naky == 0) naky = (ny - 1) / 3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx - 1) / 3) + 1
    if (rtwist == 0.) rtwist = real(jtwist)

    if (mod(ntheta0, 2) /= 1) then
       call mp_abort("ERROR: ntheta0 must be an odd number in box mode", .true.)
    end if

    ! Now we make sure that we set ny and nx for given
    ! choices of naky and ntheta0. If e.g. both ny and naky
    ! are set and they are not consistent with each other
    ! raise an error.

    if (ny == 0) then
       ! If ny hasn't been set the determine it from naky
       ny = (naky - 1) * 3 + 1
       if (proc0) write (error_unit(), '("INFO: ny (",I0,") set from naky (",I0,").")') ny, naky
    else
       ! If both naky and ny are set then check that the resulting padding is at least
       ! as much as required for dealiasing
       if (naky < (ny - 1) / 3 + 1) then ! Excess padding
          if (proc0) then
             write (error_unit(), '("INFO: Both ny (",I0,") and naky (",I0,") have been set by the user.")') ny, naky
             write (error_unit(), '("      these values lead to excess padding (",I0,".) for the FFTs.")') ny - ((naky - 1)*  3 + 1)
             write (error_unit(), '("      This may be desirable if the values chosen increase the number of parallelization sweetspots,")')
             write (error_unit(), '("      but this also makes the resolution in ny larger than is strictly necessary.")')
          end if

       else if (naky > (ny - 1) / 3 + 1) then ! Insufficient padding
          if (proc0) then
             write (error_unit(), '("ERROR: ny (",I0,") and naky (",I0,") have been set by the user.")') ny, naky
             write (error_unit(), '("       but these values lead to insufficient padding (",I0,") for the FFTS.")') ny - ((naky - 1)*  3 + 1)
             write (error_unit(), '("       The zero-padding must satisfy the one thirds rule to avoid aliasing, which requires naky <= (ny-1)/3 + 1.")')
             write (error_unit(), '("  Please do one of the following:")')
             write (error_unit(), '("       1. Increase ny or decrease naky to ensure naky <= (ny-1)/3 + 1")')
             write (error_unit(), '("       2. Set just one of ny and naky, GS2 will then set the other appropriately.")')
          end if

          call mp_abort("ERROR: naky and ny both set resulting in insufficient padding. See error file for more details.", .true.)
       end if
    end if

    if (nx == 0) then
       ! If nx hasn't been set the determine it from ntheta0
       nx = ((ntheta0 - 1) /  2) * 3 + 1
       if (proc0) write (error_unit(), '("INFO: nx (",I0,") set from ntheta0 (",I0,").")') nx, ntheta0
    else
       ! If both ntheta0 and nx are set then check that the resulting padding is at least
       ! as much as required for dealiasing
       if (ntheta0 < 2*((nx - 1) / 3) + 1) then ! Excess padding
          if (proc0) then
             write (error_unit(), '("INFO: Both nx (",I0,") and ntheta0 (",I0,") have been set by the user.")') nx, ntheta0
             write (error_unit(), '("      these values lead to excess padding (",I0,".) for the FFTs.")') nx - (((ntheta0 - 1) /  2) * 3 + 1)
             write (error_unit(), '("      This may be desirable if the values chosen increase the number of parallelization sweetspots,")')
             write (error_unit(), '("      but this also makes the resolution in ny larger than is strictly necessary.")')
          end if

       else if (ntheta0 > 2*((nx - 1) / 3) + 1) then ! Insufficient padding
          if (proc0) then
             write (error_unit(), '("ERROR: nx (",I0,") and ntheta0 (",I0,") have been set by the user.")') nx, ntheta0
             write (error_unit(), '("       but these values lead to insufficient padding (",I0,") for the FFTS.")') nx - (((ntheta0 - 1) /  2) * 3 + 1)
             write (error_unit(), '("       The zero-padding must satisfy the two thirds rule to avoid aliasing, which requires ntheta0 <= 2*((nx-1)/3)+1.")')
             write (error_unit(), '("  Please do one of the following:")')
             write (error_unit(), '("       1. Increase nx or decrease ntheta0 to ensure ntheta0 <= 2*((nx-1)/3)+1")')
             write (error_unit(), '("       2. Set just one of nx and ntheta0, GS2 will then set the other appropriately.")')
          end if

          call mp_abort("ERROR: ntheta0 and nx both set resulting in insufficient padding. See error file for more details.", .true.)
       end if
    end if

    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  !> Write namelist for kt_grids_box
  subroutine wnml_kt_grids_box (unit)
    implicit none
    integer, intent(in) :: unit

    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_box_parameters"
    write (unit, fmt="(' nx = ',i4)") nx_private
    write (unit, fmt="(' ntheta0 = ',i4)") ntheta0_private
    write (unit, fmt="(' ny = ',i4)") ny_private
    write (unit, fmt="(' naky = ',i4)") naky_private
    write (unit, fmt="(' Ly = ',e17.10)") ly
    if (rtwist /= 0.) then
       write (unit, fmt="(' rtwist = ',e17.10)") rtwist
    else
       write (unit, fmt="(' jtwist = ',i4)") jtwist
    end if
    write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_box

  !> Get the various grid sizes
  subroutine box_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    !> Number of \(k_y \rho\) modes. See [[kt_grids_box_parameters::naky]]
    integer, intent (out) :: naky
    !> Number of \(\theta_0\) points. See [[kt_grids_box_parameters::ntheta0]]
    integer, intent (out) :: ntheta0
    !> Number of x points in real space. See [[kt_grids_box_parameters::nx]]
    integer, intent (out) :: nx
    !> Number of y points in real space. See [[kt_grids_box_parameters::ny]]
    integer, intent (out) :: ny
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
  end subroutine box_get_sizes

  !> FIXME : Add documentation
  subroutine box_set_overrides(grids_ov)
    use overrides, only: kt_grids_overrides_type
    type(kt_grids_overrides_type), intent(in) :: grids_ov
    if (grids_ov%override_naky) naky = grids_ov%naky
    if (grids_ov%override_ny) ny = grids_ov%ny
    if (grids_ov%override_ntheta0) ntheta0 = grids_ov%ntheta0
    if (grids_ov%override_nx) nx = grids_ov%nx
    if (grids_ov%override_x0) x0 = grids_ov%x0
    if (grids_ov%override_y0) y0 = grids_ov%y0
  end subroutine box_set_overrides

  !> Calculate the grid of wavenumbers for box mode
  subroutine box_get_grids (aky, theta0, akx, ikx)
    use theta_grid, only: is_effectively_zero_shear, nperiod, shat
    use constants, only: pi
    use theta_grid, only: nperiod
    implicit none
    !> The \(k_x \rho\) grid
    real, dimension (:), intent (out) :: akx
    !> The \(k_y \rho\) grid
    real, dimension (:), intent (out) :: aky
    !> The \(\theta_0(k_x, k_y)\) grid
    real, dimension (:,:), intent (out) :: theta0
    !> Discrete kx wavenumber grid indices
    integer, dimension (:), intent (out) :: ikx

    real :: dkx, dky, ratio
    integer :: i, naky, ntheta0

    naky = size(aky)
    ntheta0 = size(akx)

    dky = 2.0*pi/ly

    if(is_effectively_zero_shear()) then   ! non-quantized b/c assumed to be periodic instead linked boundary conditions

       if (x0 == 0.) then

          if (rtwist > 0) then
             ratio = rtwist
          else
             ratio = 1. / abs(rtwist)
          end if

          dkx = dky / ratio

       else

          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if

       end if

    else
       if (jtwist /= 0) then
          dkx = dky * 2.0*pi*(2*nperiod-1)*abs(shat)/real(jtwist)
       else
          dkx = dky
       end if
    endif

    do i = 1, naky
       aky(i) = real(i-1)*dky
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

    if (shat /= 0.) then
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = akx(i)/(aky(2:)*shat)
       end do
    else
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = - akx(i)/aky(2:)   ! meaningless, so be careful
       end do
    end if

  end subroutine box_get_grids

  !> FIXME : Add documentation
  subroutine check_kt_grids_box(report_unit)
    use theta_grid, only: shat, is_effectively_zero_shear
    use constants, only: pi
    implicit none
    integer, intent(in) :: report_unit
    real :: lx
    integer :: nx, ny, naky, ntheta0
    integer :: naky_calc, ntheta0_calc

    naky=naky_private
    ntheta0=ntheta0_private
    nx=nx_private
    ny=ny_private

    if (y0 /= 2.) then
       if (abs(2.0*pi*y0 - ly) > 1.0e-7) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You cannot specify both ly and y0.')")
          write (report_unit, fmt="(' ly=',e12.4,'  2.0*pi*y0=',2e12.4)") ly
          write (report_unit, fmt="('THIS IS AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end if
    end if

    write (report_unit, *)
    write (report_unit, fmt="('A rectangular simulation domain has been selected.')")
    if (rhostar_box .gt. 0.0 .and. n0 .gt. 0) write (report_unit, fmt="('The flux-tube size corresponds to toroidal mode number, n0=',i8/T44,' at rhostar_box=',1pe12.4)") n0,rhostar_box
    write (report_unit, *)
    write (report_unit, fmt="('The domain is ',f10.4,' rho in the y direction.')") ly


    if (is_effectively_zero_shear()) then
       if (x0 == 0.) then
          if (rtwist > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  ly*rtwist
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -ly/rtwist
          end if
       else
          if (x0 > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  2.*pi*x0
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -2.*pi/x0
          end if
       end if
    else
       lx = ly * jtwist / (2.*pi*abs(shat))
       write (report_unit, fmt="('At theta=0, the domain is ',f10.4,' rho in the x direction.')") lx
       write (report_unit,*) ly, rtwist, jtwist, pi, shat
    end if

    write (report_unit, *)
    write (report_unit, fmt="('The nonlinear terms will be evaluated on a grid with ',&
         & i4,' points in x and ',i4,' points in y.')") nx, ny
    write (report_unit, *)
    naky_calc = (ny-1)/3+1
    ntheta0_calc = 2*((nx-1)/3)+1
    write (report_unit, fmt="('After de-aliasing, there will be ',i4,'  ky >= 0 modes and ',i4,' kx modes.')") naky_calc, ntheta0_calc
    write (report_unit, fmt="('The modes with ky < 0 are determined by the reality condition.')")

    if ( naky_calc < naky ) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('ERROR : The requested number of ky >= 0 modes is ',I0,' which exceeds that available after de-aliasing.')") naky
       write (report_unit, fmt="('THIS IS AN ERROR.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    else if ( naky_calc > naky ) then
       write (report_unit, *)
       write (report_unit, fmt="('The requested number of ky >= 0 modes is ',I0,' which is less than available after de-aliasing.')") naky
       write (report_unit, fmt="('This corresponds to extra padding in the y grid of ',I0,' points.')") ny - ((naky - 1)*  3 + 1)
       write (report_unit, *)
    end if

    if ( ntheta0_calc < ntheta0 ) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('ERROR : The requested number of kx modes is ',I0,' which exceeds that available after de-aliasing.')") ntheta0
       write (report_unit, fmt="('THIS IS AN ERROR.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    else if ( ntheta0_calc > ntheta0 ) then
       write (report_unit, *)
       write (report_unit, fmt="('The requested number of kx modes is ',I0,' which is less than available after de-aliasing.')") ntheta0
       write (report_unit, fmt="('This corresponds to extra padding in the x grid of ',I0,' points.')") nx - (((ntheta0 - 1) /  2) * 3 + 1)
       write (report_unit, *)
    end if

  end subroutine check_kt_grids_box

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_box_config(kt_grids_box_config_in)
    use mp, only: mp_abort
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_box_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_box_config_in)) then
       kt_grids_box_config = kt_grids_box_config_in
    end if
  end subroutine set_kt_grids_box_config

#include "kt_grids_box_auto_gen.inc"
end module kt_grids_box
