!> FIXME : Add documentation
!!
!! @note nexp has been changed in a rush.  Only known to be correct for nexp=2
module hyper
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  implicit none

  private

  public :: init_hyper, finish_hyper, hyper_diff, D_res
  public :: read_parameters, wnml_hyper, check_hyper
  public :: D_v, D_eta, nexp
  public :: hypervisc_filter

  public :: hyper_config_type
  public :: set_hyper_config
  public :: get_hyper_config
  
! D_v, D_eta are hyper coefficients, normalized correctly 
! i.e., either by unity or by 1/k_perp**2*nexp

  real :: D_v, D_eta
  real :: D_hypervisc, D_hyperres, omega_osc, D_hyper
  real :: akx4_max, aky4_max, aky_max, akperp4_max

  integer :: hyper_option_switch, nexp
  integer, parameter :: hyper_option_none = 1, &
       hyper_option_visc = 2, &
       hyper_option_res  = 3, &
       hyper_option_both = 4, &
       hyper_option_simple3D = 5 ! MRH
       
  character(9) :: hyper_option
  logical :: const_amp, include_kpar, isotropic_shear, damp_zonal_only
  logical :: hyper_on = .false.
  logical :: gridnorm

  real, dimension (:,:), allocatable :: D_res
  ! (it, ik)

  real, dimension (:,:,:), allocatable :: hypervisc_filter
  ! (-ntgrid:ntgrid, ntheta0, naky)

  ! ! MRH Coefficients for simple3D hyperviscous model
  real :: kperp2_max, kperp2_max_zonal
  ! the maximum k_perp ^ 2 including geometric coefficients for non-zonal and zonal modes
  real :: D_hyper3D, P_hyper3D, ky_cut, kx_cut
  ! Hyperviscous coefficient, hyperviscous power, ky cut-off, kx cut-off
  logical :: isotropic_model
  ! Controls if zonal modes are treated identically to non-zonal modes

  logical :: initialized = .false.


  !> Used to represent the input configuration of hyper
  type, extends(abstract_config_type) :: hyper_config_type
     ! namelist : hyper_knobs
     ! indexed : false
     !> Determines whether hyperviscosity includes time dependent amplitude
     !> factor when calculating damping rate. Recommend `true` for
     !> linear runs and `false` for nolinear runs, since amplutide of
     !> turbulence grows linearly with time in linear run.
     logical :: const_amp = .false.
     !> If `hyper_option = 'both'` is used then this sets both the
     !> hyperresistivity and hyperviscosity damping coefficients. Can
     !> override the individual coefficients with
     !> [[hyper_knobs:D_hyperres]] and [[hyper_knobs:D_hypervisc]].
     real :: d_hyper = -10.0
     !> Used with the simple3D hyperviscosity model of the form
     !> D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D 
     !> and ky_cut, kx_cut set max |kperp|
     real :: d_hyper3d = -10.
     !> Sets hyperresistivity parameter multiplying damping term.
     real :: d_hyperres = -10.0
     !> Sets hyperviscosity parameter multiplying damping term. See
     !> [E. Belli (2006)
     !> thesis](https://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=50BFD54A8F8D72FC225D025FDEDFFFA5?doi=10.1.1.706.9568&rep=rep1&type=pdf)
     !> for more information.
     real :: d_hypervisc = -10.0
     !> If `true` then hyperdissipation only applied to the zonal
     !> mode.
     logical :: damp_zonal_only = .false.
     !> If `true` (default) then set wavenumber parameters entering
     !> the models based on the maximum `ky` and `kx` included in the
     !> current simulation. If `false` then these values are set to 1.
     logical :: gridnorm = .true.
     !> Selects the type of hyper terms included. Should be one of
     !>
     !> - 'default' -- no hyper terms included.
     !> - 'none' -- the same as default.
     !> - 'visc_only' -- only hyperviscosity included.
     !> - 'res_only' -- only hyperresistivity included.
     !> - 'both' -- both viscosity and resistivity included.
     !> - 'simple3D' -- simple hyperviscous dissipation rate of the
     !>    form D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D
     !>    described in [“Multiscale turbulence in magnetic
     !>    confinement fusion devices”, M. Hardman, DPhil Thesis,
     !>    appendix
     !>    B.4](https://ora.ox.ac.uk/objects/uuid:233a22cb-3c8b-4fe0-a689-4a37d8fe0314)
     !>    note that the dissipation in this version is applied to g, not h, as
     !>    in the reference. 
     character(len = 9) :: hyper_option = 'default'
     !> Not used.
     !>
     !> @todo Remove this variable.
     logical :: include_kpar = .false.
     !> if true damp zonal and drift waves with same dissipation formula
     logical :: isotropic_model = .true.
     !> If `true` then use isotropic shear model.
     logical :: isotropic_shear = .true.
     !> Used with the simple3D hyperviscosity model of the form
     !> D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D 
     !> and ky_cut, kx_cut set max |kperp|
     real :: kx_cut = -10.
     !> Used with the simple3D hyperviscosity model of the form
     !> D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D 
     !> and ky_cut, kx_cut set max |kperp|
     real :: ky_cut = -10.
     !> Sets the power to which \(k_\bot^2\) is raised in the dissipation filter.
     integer :: nexp = 2
     !> Sets a parameter in the anisotropic shearing rate calculation.
     real :: omega_osc = 0.4
     !> Used with the simple3D hyperviscosity model of the form
     !> D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D 
     !> and ky_cut, kx_cut set max |kperp|
     real :: p_hyper3d = 4.
   contains
     procedure, public :: read => read_hyper_config
     procedure, public :: write => write_hyper_config
     procedure, public :: reset => reset_hyper_config
     procedure, public :: broadcast => broadcast_hyper_config
     procedure, public, nopass :: get_default_name => get_default_name_hyper_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_hyper_config
  end type hyper_config_type

  type(hyper_config_type) :: hyper_config
  
contains

  !> FIXME : Add documentation
  subroutine check_hyper(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    select case (hyper_option_switch)
    case (hyper_option_none)
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
               &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  &
               &D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_visc)
       hyper_on = .true.
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
               &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but &
               &D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       endif

    case (hyper_option_res)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write(report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_both)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       endif
       if (D_hypervisc < 0. .and. D_hyperres < 0.) hyper_on = .false.

    end select

    if (hyper_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('------------------------------------------------------------')")
       write (report_unit, *) 

       select case (hyper_option_switch)

       case (hyper_option_visc)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperviscosity included without hyperresistivity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hypervisc
          else
             write (report_unit, fmt="('The damping coefficent is ',e11.4)") D_hypervisc
             write (report_unit, fmt="('The damping rate is proportional to the RMS amplitude of the turbulence.')")
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperviscosity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('The hyperviscosity is anisotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for drift-type calculations.')")
             write (report_unit, fmt="('omega_osc = ',e11.4)") omega_osc
          end if

       case (hyper_option_res)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity included without hyperviscosity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperresistivity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if

       case (hyper_option_both)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity and hyperviscosity included.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The damping is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
       end select
    end if
  end subroutine check_hyper

  !> FIXME : Add documentation
  subroutine wnml_hyper(unit)
    implicit none
    integer, intent(in) :: unit          
    if (.not. hyper_on) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "hyper_knobs"

    select case (hyper_option_switch)

    case (hyper_option_visc) 
       write (unit, fmt="(' hyper_option = ',a)") '"visc_only"'
       write (unit, fmt="(' D_hypervisc = ',e17.10)") D_hypervisc

    case (hyper_option_res) 
       write (unit, fmt="(' hyper_option = ',a)") '"res_only"'
       write (unit, fmt="(' D_hyperres = ',e17.10)") D_hyperres

    case (hyper_option_both) 
       write (unit, fmt="(' hyper_option = ',a)") '"both"'
       if (D_hyperres == D_hypervisc) then
          write (unit, fmt="(' D_hyper = ',e17.10)") D_hyper
       else
          write (unit, fmt="(' D_hypervisc = ',e17.10)") D_hypervisc
          write (unit, fmt="(' D_hyperres = ',e17.10)") D_hyperres
       end if
    end select

    select case (hyper_option_switch)
    
    case (hyper_option_visc, hyper_option_res, hyper_option_both)
!    write (unit, fmt="(' include_kpar = ',L1)") include_kpar
    
    write (unit, fmt="(' const_amp = ',L1)") const_amp
    write (unit, fmt="(' isotropic_shear = ',L1)") isotropic_shear
    if (.not. isotropic_shear) &
         write (unit, fmt="(' omega_osc = ',e17.10)") omega_osc

    write (unit, fmt="(' gridnorm = ',L1)") gridnorm
    write (unit, fmt="(' /')")
    
    case (hyper_option_simple3D)
    write (unit, fmt="(' hyper_option = ',a)") '"simple3D"'
    write (unit, fmt="(' D_hyper3D = ',e17.10)") D_hyper3D
    write (unit, fmt="(' P_hyper3D = ',e17.10)") P_hyper3D
    if (ky_cut > 0.) &
         write (unit, fmt="(' ky_cut = ',e17.10)") ky_cut
    if (kx_cut > 0.) &
         write (unit, fmt="(' kx_cut = ',e17.10)") kx_cut
    write (unit, fmt="(' isotropic_model = ',L1)") isotropic_model
    write (unit, fmt="(' /')")
        
    end select
  end subroutine wnml_hyper

  !> FIXME : Add documentation
  subroutine init_hyper(hyper_config_in)
    use kt_grids, only: ntheta0, naky, akx, aky
    use kt_grids, only: kperp2 ! MRH
    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat ! MRH
    use gs2_time, only: code_dt
    use gs2_layouts, only: init_gs2_layouts
    implicit none
    type(hyper_config_type), intent(in), optional :: hyper_config_in        
    integer :: ik, it

    if (initialized) return
    initialized = .true.
    
    call init_gs2_layouts
    call read_parameters(hyper_config_in)
    call allocate_arrays

    select case (hyper_option_switch)
        case(hyper_option_both, hyper_option_res, hyper_option_visc)
    
        ! Define variables used in hyperviscosity and hyperresistivity models

            if (gridnorm) then
               akx4_max    = akx(ntheta0/2 + 1) ** (2*nexp)
               aky_max     = aky(naky)
               aky4_max     = aky(naky) ** (2*nexp)
               akperp4_max = ( akx(ntheta0/2 + 1) ** 2  +  aky(naky) ** 2) ** (nexp)
            else
               akx4_max = 1.
               aky_max = 1.
               aky4_max = 1.
               akperp4_max = 1.
            end if

        ! Get D_res set up if needed

            if (D_hyperres > 0.) then
               do ik = 1, naky
                  do it = 1, ntheta0
                     D_res(it, ik) = D_hyperres*code_dt &
                          * (aky(ik)**2 + akx(it)**2)**nexp/akperp4_max
                  end do
               end do
               D_eta = D_hyperres/akperp4_max
            else
               D_res = 0.
               D_eta = 0.
            end if

            if (D_hypervisc > 0.) then
               D_v = D_hypervisc/akperp4_max
            else
               D_v = 0.
            end if
            
        case(hyper_option_simple3D)  ! MRH
            ! Make the default cut-offs the grid scale
            if(.not. ky_cut > 0.) ky_cut = aky(naky)
            if(.not. kx_cut > 0.) kx_cut = akx((ntheta0+1)/2)
            
            kperp2_max = maxval(ky_cut*ky_cut*gds2(:) + &
                              2.0*ky_cut*kx_cut*gds21(:)/shat + &
                              kx_cut*kx_cut*gds22/(shat*shat))
            if (isotropic_model) then 
                kperp2_max_zonal = kperp2_max
            else 
                kperp2_max_zonal = maxval(kx_cut*kx_cut*gds22/(shat*shat))
            endif
            ! Calculate the filter
            ik =1  ! Zonal modes
            do it = 1, ntheta0
                 hypervisc_filter(:,it,ik) = exp(- ( D_hyper3D * code_dt * &
               ( kperp2(:,it,ik)/kperp2_max_zonal )**(P_hyper3D/2.0)) ) 
            end do
              
            do ik = 2, naky ! Non-zonal modes
              do it = 1, ntheta0
                 hypervisc_filter(:,it,ik) = exp(- ( D_hyper3D * code_dt * &
               ( kperp2(:,it,ik)/kperp2_max )**(P_hyper3D/2.0)) ) 
              end do
            end do

                    
    end select

  end subroutine init_hyper

  !> FIXME : Add documentation
  subroutine read_parameters(hyper_config_in)
    use file_utils, only: input_unit, error_unit
    use text_options, only: text_option, get_option_value
    use mp, only: proc0
    implicit none
    type(hyper_config_type), intent(in), optional :: hyper_config_in    
    type (text_option), dimension(6), parameter :: hyperopts = &
         (/ text_option('default', hyper_option_none), &
            text_option('none', hyper_option_none), &
            text_option('visc_only', hyper_option_visc), &
            text_option('res_only', hyper_option_res), &
            text_option('both', hyper_option_both), &            
            text_option('simple3D', hyper_option_simple3D) /)
    integer :: ierr
    logical :: exist

    if (present(hyper_config_in)) hyper_config = hyper_config_in

    call hyper_config%init(name = 'hyper_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    const_amp = hyper_config%const_amp
    d_hyper = hyper_config%d_hyper
    d_hyper3d = hyper_config%d_hyper3d
    d_hyperres = hyper_config%d_hyperres
    d_hypervisc = hyper_config%d_hypervisc
    damp_zonal_only = hyper_config%damp_zonal_only
    gridnorm = hyper_config%gridnorm
    hyper_option = hyper_config%hyper_option
    include_kpar = hyper_config%include_kpar
    isotropic_model = hyper_config%isotropic_model
    isotropic_shear = hyper_config%isotropic_shear
    kx_cut = hyper_config%kx_cut
    ky_cut = hyper_config%ky_cut
    nexp = hyper_config%nexp
    omega_osc = hyper_config%omega_osc
    p_hyper3d = hyper_config%p_hyper3d

    exist = hyper_config%exist
    
    ierr = error_unit()
    
    call get_option_value &
         (hyper_option, hyperopts, hyper_option_switch, &
         ierr, "hyper_option in hyper_knobs",.true.)
    
    if (.not. isotropic_shear .and. nexp /=2) then
       if (proc0) write (ierr, *) 'Forcing nexp = 2.  Higher values not implemented for anisotropic shear model.'
       nexp = 2
    end if

    
    select case (hyper_option_switch)

    case (hyper_option_none)
       if (D_hyperres > 0.) then
          if (proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses no hyperresistivity.  D_hyperres ignored.'
          D_hyperres = -10.
       end if
       if (D_hypervisc > 0.) then
          if (proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses no hyperviscosity.  D_hypervisc ignored.'
          D_hypervisc = -10.
       endif

    case (hyper_option_visc)
       hyper_on = .true.
       if (proc0) write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
            'D_hyperres.  Recommend: Set them equal.'
       if (D_hyperres > 0.) then
          if (proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses no hyperresistivity.  D_hyperres ignored.'
          D_hyperres = -10.
       end if
       if (D_hypervisc < 0.) then
          if (proc0) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperviscosity but D_hypervisc < 0 is illegal.'
             write(ierr, *) 'No hyperviscosity used.'
          end if
          hyper_on = .false.
       endif

    case (hyper_option_res)
       hyper_on = .true.
       if (proc0) write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
            'D_hyperres.  Recommend: Set them equal.'
       if (D_hyperres < 0.) then
          if (proc0) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperresistivity but  D_hyperres < 0 is illegal.'
             write(ierr, *) 'No hyperresistivity used.'
          end if
          hyper_on = .false.
       end if
       if (D_hypervisc > 0.) then
          if (proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses no hyperviscosity.  D_hypervisc ignored.'
          D_hypervisc = -10.
       endif

    case (hyper_option_both)
       hyper_on = .true.

       if (D_hyper < 0.) then
          if (D_hyperres /= D_hyperres) then
             if (proc0) write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
                  'D_hyperres.  Recommend: Set them equal.'
          end if
       else
          if (proc0) write (ierr, *) 'WARNING: Setting D_hypervisc = D_hyperres, each to value of D_hyper'
          D_hyperres  = D_hyper
          D_hypervisc = D_hyper
       end if

       if (D_hyperres < 0.) then
          if (proc0) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperresistivity but  D_hyperres < 0 is illegal.'
             write(ierr, *) 'No hyperresistivity used.'
          end if
       end if
       if (D_hypervisc < 0.) then
          if (proc0) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperviscosity but D_hypervisc < 0 is illegal.'
             write(ierr, *) 'No hyperviscosity used.'
          end if
       endif
       if (D_hypervisc < 0. .and. D_hyperres < 0.) hyper_on = .false.

    case (hyper_option_simple3D) ! MRH
       hyper_on = .true.

       if (P_hyper3D < 4.) then
          hyper_on = .false.
          if(proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses a simple hyperviscous filter but  P_hyper3D < 4. is illegal.'
          if(proc0) write(ierr, *) 'No hyperviscosity used.'
       endif

       if (D_hyper3D < 0.) then
          hyper_on = .false.
          if(proc0) write(ierr, *) 'hyper_option = ',hyper_option, &
               ' chooses a simple hyperviscous filter but  D_hyper3D < 0 is illegal.'
          if(proc0) write(ierr, *) 'No hyperviscosity used.'
       endif
    end select
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    implicit none

    if (.not. allocated(D_res)) then
       allocate (D_res(ntheta0, naky)) 
    endif
    if (.not. allocated(hypervisc_filter)) then
       allocate (hypervisc_filter(-ntgrid:ntgrid,ntheta0,naky)) ; hypervisc_filter = 1.0
    end if
    D_res = 0.
    
  end subroutine allocate_arrays

  !> FIXME : Add documentation
  subroutine hyper_diff (g0, phi)

    use gs2_layouts, only: ik_idx, it_idx, is_idx
    use theta_grid, only: ntgrid
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo
    use kt_grids, only: aky, akx, naky, ntheta0

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi

    real, dimension (-ntgrid:ntgrid) :: shear_rate_nz, shear_rate_z, shear_rate_z_nz

    integer :: iglo, ik, it
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Includes models by Belli and Hammett
! to calculate the x-y avged shearing rate
! S(theta)^2 = <|grad_perp|^4 |phi|^2> 
!            =  sum_over_k(kperp^4 * |phi|^2)
!               (times crazy fac factor due to FFT conventions.)
!
! and to implement this anisotropically in k_perp, taking into 
! account properties of zonal flows.
!
! Begun December 2001
!
! Number conservation added April 2006
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not. hyper_on) return
    
    select case (hyper_option_switch)
        case(hyper_option_both, hyper_option_res, hyper_option_visc)
            
            if (D_hypervisc < 0.) return

            if(isotropic_shear) then
               call iso_shear
            else
               call aniso_shear
            end if
            
        case(hyper_option_simple3D) ! MRH
            
            call simple3Dfilter
            
    end select    

  contains
    !> FIXME : Add documentation
    subroutine simple3Dfilter
    
        implicit none
        
        do iglo = g_lo%llim_proc, g_lo%ulim_alloc
         ik = ik_idx(g_lo, iglo)
         it = it_idx(g_lo, iglo)
                  
         g0(:,1,iglo) = g0(:,1,iglo) * hypervisc_filter(:,it,ik)
         g0(:,2,iglo) = g0(:,2,iglo) * hypervisc_filter(:,it,ik)
        end do
      
    end subroutine simple3Dfilter

    !> FIXME : Add documentation
    subroutine aniso_shear

      real, dimension(naky) :: fac
      
! Do the Belli-Hammett anisotropic calculation 
! which accounts for some zonal/non-zonal differences
       
      fac = 0.5
      fac(1) = 1.0
         
! shearing rate due to non-zonal modes (on nonzonal modes)
      shear_rate_nz = 0.
      do ik = 2, naky
         do it = 1, ntheta0
            shear_rate_nz(:) = shear_rate_nz(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
         end do
      end do
      shear_rate_nz = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_nz) ** 0.5 )
       
! shearing rate due to zonal modes (on nonzonal modes)
      shear_rate_z = 0.
      do ik = 1, 1
         do it = 1, ntheta0
            shear_rate_z(:) = shear_rate_z(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
         end do
      end do
! shear_rate_z = shear_rate_z ** 0.5
      shear_rate_z = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_z) ** 0.5 )
      
! shearing rate due to nonzonal modes (on zonal modes)
      shear_rate_z_nz = 0.
      do ik = 2, naky
         do it = 1, ntheta0
            shear_rate_z_nz(:) = shear_rate_z_nz(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) *  aky(ik)**4 * fac(ik)
         end do
      end do
! shear_rate_z_nz = shear_rate_z_nz ** 0.5
      shear_rate_z_nz = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_z_nz) ** 0.5 )
       
! end of anisotropic shearing rate calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo, iglo)
         it = it_idx(g_lo, iglo)
         
         if(aky(ik) == 0.) then
            hypervisc_filter(:,it,ik) = exp(- (D_hypervisc * code_dt &
                 * ( shear_rate_z_nz(:) * akx(it) ** 4 / akx4_max )))
         else
            hypervisc_filter(:,it,ik) = exp(- ( D_hypervisc * code_dt & 
                 * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max & 
                 + shear_rate_z(:) * akx(it) ** 4/ akx4_max * aky(ik) / aky_max )))
         endif
         
         g0(:,1,iglo) = g0(:,1,iglo) * hypervisc_filter(:,it,ik)
         g0(:,2,iglo) = g0(:,2,iglo) * hypervisc_filter(:,it,ik)
      end do
    
    end subroutine aniso_shear

    !> FIXME : Add documentation
    subroutine iso_shear

      real, dimension(naky) :: fac
      
      if (const_amp) then
         shear_rate_nz = 1.
      else
         fac = 0.5
         fac(1) = 1.0
         shear_rate_nz = 0.
         do ik = 1, naky
            do it = 1, ntheta0
               shear_rate_nz(:) = shear_rate_nz(:) &
                    + real(conjg(phi(:,it,ik))*phi(:,it,ik)) &
                    * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
            end do
         end do
         shear_rate_nz = shear_rate_nz**0.5
      end if

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo, iglo)
         it = it_idx(g_lo, iglo)
         if (damp_zonal_only .and. .not. aky(ik)==epsilon(0.0)) cycle
         hypervisc_filter(:,it,ik) = exp(- ( D_hypervisc * code_dt &
              * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max)))
         
         g0(:,1,iglo) = g0(:,1,iglo) * hypervisc_filter(:,it,ik)
         g0(:,2,iglo) = g0(:,2,iglo) * hypervisc_filter(:,it,ik)
      end do
      
    end subroutine iso_shear

  end subroutine hyper_diff

  !> FIXME : Add documentation
  subroutine finish_hyper

    implicit none

    hyper_on = .false.
    initialized = .false.
    if (allocated(D_res)) deallocate (D_res)
    if (allocated(hypervisc_filter)) deallocate (hypervisc_filter)

    call hyper_config%reset()
  end subroutine finish_hyper

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_hyper_config(hyper_config_in)
    use mp, only: mp_abort
    type(hyper_config_type), intent(in), optional :: hyper_config_in
    if (initialized) then
       call mp_abort("Trying to set hyper_config when already initialized.", to_screen = .true.)
    end if
    if (present(hyper_config_in)) then
       hyper_config = hyper_config_in
    end if
  end subroutine set_hyper_config

    !> Get the module level config instance
  function get_hyper_config()
    type(hyper_config_type) :: get_hyper_config
    get_hyper_config = hyper_config
  end function get_hyper_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the hyper_knobs namelist and populates the member variables
  subroutine read_hyper_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(hyper_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 9) :: hyper_option
    integer :: nexp
    logical :: const_amp, damp_zonal_only, gridnorm, include_kpar, isotropic_model, isotropic_shear
    real :: d_hyper, d_hyper3d, d_hyperres, d_hypervisc, kx_cut, ky_cut, omega_osc, p_hyper3d

    namelist /hyper_knobs/ const_amp, d_hyper, d_hyper3d, d_hyperres, d_hypervisc, damp_zonal_only, gridnorm, hyper_option, include_kpar, isotropic_model, isotropic_shear, kx_cut, ky_cut, nexp, omega_osc, p_hyper3d

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    const_amp = self%const_amp
    d_hyper = self%d_hyper
    d_hyper3d = self%d_hyper3d
    d_hyperres = self%d_hyperres
    d_hypervisc = self%d_hypervisc
    damp_zonal_only = self%damp_zonal_only
    gridnorm = self%gridnorm
    hyper_option = self%hyper_option
    include_kpar = self%include_kpar
    isotropic_model = self%isotropic_model
    isotropic_shear = self%isotropic_shear
    kx_cut = self%kx_cut
    ky_cut = self%ky_cut
    nexp = self%nexp
    omega_osc = self%omega_osc
    p_hyper3d = self%p_hyper3d

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = hyper_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = hyper_knobs)
    endif

    ! Now copy from local variables into type members
    self%const_amp = const_amp
    self%d_hyper = d_hyper
    self%d_hyper3d = d_hyper3d
    self%d_hyperres = d_hyperres
    self%d_hypervisc = d_hypervisc
    self%damp_zonal_only = damp_zonal_only
    self%gridnorm = gridnorm
    self%hyper_option = hyper_option
    self%include_kpar = include_kpar
    self%isotropic_model = isotropic_model
    self%isotropic_shear = isotropic_shear
    self%kx_cut = kx_cut
    self%ky_cut = ky_cut
    self%nexp = nexp
    self%omega_osc = omega_osc
    self%p_hyper3d = p_hyper3d

    self%exist = exist
  end subroutine read_hyper_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_hyper_config(self, unit)
    implicit none
    class(hyper_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("const_amp", self%const_amp, unit_internal)
    call self%write_key_val("d_hyper", self%d_hyper, unit_internal)
    call self%write_key_val("d_hyper3d", self%d_hyper3d, unit_internal)
    call self%write_key_val("d_hyperres", self%d_hyperres, unit_internal)
    call self%write_key_val("d_hypervisc", self%d_hypervisc, unit_internal)
    call self%write_key_val("damp_zonal_only", self%damp_zonal_only, unit_internal)
    call self%write_key_val("gridnorm", self%gridnorm, unit_internal)
    call self%write_key_val("hyper_option", self%hyper_option, unit_internal)
    call self%write_key_val("include_kpar", self%include_kpar, unit_internal)
    call self%write_key_val("isotropic_model", self%isotropic_model, unit_internal)
    call self%write_key_val("isotropic_shear", self%isotropic_shear, unit_internal)
    call self%write_key_val("kx_cut", self%kx_cut, unit_internal)
    call self%write_key_val("ky_cut", self%ky_cut, unit_internal)
    call self%write_key_val("nexp", self%nexp, unit_internal)
    call self%write_key_val("omega_osc", self%omega_osc, unit_internal)
    call self%write_key_val("p_hyper3d", self%p_hyper3d, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_hyper_config

  !> Resets the config object to the initial empty state
  subroutine reset_hyper_config(self)
    class(hyper_config_type), intent(in out) :: self
    type(hyper_config_type) :: empty
    select type (self)
    type is (hyper_config_type)
       self = empty
    end select
  end subroutine reset_hyper_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_hyper_config(self)
    use mp, only: broadcast
    implicit none
    class(hyper_config_type), intent(in out) :: self
    call broadcast(self%const_amp)
    call broadcast(self%d_hyper)
    call broadcast(self%d_hyper3d)
    call broadcast(self%d_hyperres)
    call broadcast(self%d_hypervisc)
    call broadcast(self%damp_zonal_only)
    call broadcast(self%gridnorm)
    call broadcast(self%hyper_option)
    call broadcast(self%include_kpar)
    call broadcast(self%isotropic_shear)
    call broadcast(self%kx_cut)
    call broadcast(self%ky_cut)
    call broadcast(self%nexp)
    call broadcast(self%omega_osc)
    call broadcast(self%p_hyper3d)
    call broadcast(self%exist)
  end subroutine broadcast_hyper_config

  !> Gets the default name for this namelist
  function get_default_name_hyper_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_hyper_config
    get_default_name_hyper_config = "hyper_knobs"
  end function get_default_name_hyper_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_hyper_config()
    implicit none
    logical :: get_default_requires_index_hyper_config
    get_default_requires_index_hyper_config = .false.
  end function get_default_requires_index_hyper_config
end module hyper
