!> FIXME : Add documentation
module theta_grid_eik
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use geometry, only: eikcoefs_output_type
  use geo_utils, only: EQFILE_LENGTH
  implicit none

  private

  public :: init_theta_grid_eik, finish_theta_grid_eik, check_theta_grid_eik, wnml_theta_grid_eik
  public :: eik_get_sizes, eik_get_grids

  public :: theta_grid_eik_config_type
  public :: set_theta_grid_eik_config
  public :: get_theta_grid_eik_config

  logical :: exist, initialized = .false.
  type(eikcoefs_output_type) :: eikcoefs_results
  integer :: ntheta_geometry

  !> Used to represent the input configuration of theta_grid
  type, extends(abstract_config_type) :: theta_grid_eik_config_type
     ! namelist : theta_grid_eik_knobs
     ! indexed : false
     !> Used in calculation of `dp_new = -alpha_input/qval**2/rmaj*drhodpsi`
     !> @note This gets a "smart" default.
     real :: alpha_input = 0.0
     !> The gradient of the pressure. Strictly speaking this parameter
     !> is not \(\frac{\partial \beta}{\partial \rho}\) but \(\beta
     !> \frac{1}{p}\frac{\partial p}{\partial \rho}\): in other words,
     !> the gradient of the magnetic field is ignored. Used only if
     !> `bishop` = 4 or 9.
     !> @note This gets a "smart" default.
     real :: beta_prime_input = 0.0
     !> Use Bishop relations to generate metric coefficients.
     !>
     !> -  0: Use high-aspect ratio coefficients (only for debugging)
     !> -  1: Use actual equilibrium values of shat, p' recommended
     !> -  2: Use numerical equilibrium + s_hat_input and p_prime_input
     !> -  3: Use numerical equilibrium + s_hat_input and inv_Lp_input
     !> -  4: Use numerical equilibrium + s_hat_input and beta_prime_input
     !> -  5: Use numerical equilibrium + s_hat_input and alpha_input
     !> -  6: Use numerical equilibrium + beta_prime_input
     !> -  7: Use numerical equilibrium and multiply pressure gradient by dp_mult
     !> -  8: Use numerical equilibrium + s_hat_input and multiply pressure gradient by dp_mult
     !> -  9: Use numerical equilibrium + s_hat_input and beta_prime_input
     !> -  Otherwise: Use magnetic shear and pressure gradient as set elsewhere.
     !>
     integer :: bishop = 5
     !> Use equilbrium data from the CHEASE file ogyropsi.dat
     logical :: chs_eq = .false.
     !> Step size for radial derivatives, \(\Delta r_{\psi N}\). Should be
     !> "small enough", typically 0.001.
     real :: delrho = 1e-3
     !> Vacuum magnetic dipole geometry
     logical :: dfit_eq = .false.
     !> Used to scale the pressure gradient, only if bishop = 7 or 8.
     real :: dp_mult = 1.0
     !> Use EFIT equilibrium (EFIT, codes with eqdsk format)
     logical :: efit_eq = .false.
     !> Name of file with numerical equilibrium data (if required)
     character(len = EQFILE_LENGTH) :: eqfile = "default_unset_value"
     !> Name of file with numerical equilibrium normalization data (if required)
     !> currently, only used for dipole equilibrium (deq)
     character(len = EQFILE_LENGTH) :: eqnormfile = "default_unset_value"
     !> Change field-line coordinate. Recommended value: F
     !> @note We recommend `.false.` but default to `.true.`. We should consider
     !> changing the default.
     logical :: equal_arc = .true.
     !> If true then forces up-down symmetry in some geometrical quantities
     logical :: force_sym = .false.
     !> Use Toq-style NetCDF equilibrium (TOQ)
     logical :: gen_eq = .false.
     !> Read Colin Roach's GS2D equilibrium file
     logical :: gs2d_eq = .false.
     !> Unknown equilibrium file. You probably don't want this.
     !> FIXME: Add documentation
     logical :: idfit_eq = .false.
     !> Deprecated -- redundant information.
     integer :: iflux = -1
     !> Used with [[theta_grid_eik_knobs:bishop]] == 3: controls pressure length
     !> scale by multiplying \(p \frac{d \rho}{d \psi}\)
     real :: invlp_input = 0.
     !> Choose definition of flux surface coordinate
     !>
     !> -  1: rho == sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)
     !> -  2: rho == midplane diameter/LCFS diameter - Recommended
     !> -  3: rho == poloidal flux/poloidal flux of LCFS
     !> -  4: rho == rho_mid (vacuum ring dipole, `dfit_eq = T`, only)
     !>
     !> NB For consistency fprim, tprim, shat, uprim, g_exb etc
     !>  *must be computed using same radial variable*,
     !>  i.e. *depend on choice of irho*!
     integer :: irho = 2
     !> Deprecated, see force_sym instead
     integer :: isym = -1
     !> Deprecated, see use_large_aspect instead
     integer :: itor = -1
     !> If `.true.` use Miller-style local equilibrium else use other
     !> numerical equilibrium types
     logical :: local_eq = .true.
     !> The number of theta grid points to use in eikcoefs calls.
     !> Currently may not have an effect for all equilibrium types.
     !> If not set then defaults to [[theta_grid_parameters:ntheta]]
     integer :: ntheta_geometry = -1
     !> Use Menard-style NetCDF equilibrium (JSOLVER)
     logical :: ppl_eq = .false.
     !> Used to overrides s_hat prescribed by the numerical
     !> equilibrium, but _only_ if bishop=2,3,4,5,8, or 9.
     !> @note This gets a "smart" default.
     real :: s_hat_input = 0.0
     !> Use PPL NetCDF equilibrium (psipqgrz equilibrium from TRANSP/TRXPL)
     logical :: transp_eq = .false.
     !> If true use large aspect ratio expansions in eik geometry to get ~s-alpha
     logical :: use_large_aspect = .false.
     !> Write a little extra about geometry to the screen.
     logical :: writelots = .false.
   contains
     procedure, public :: read => read_theta_grid_eik_config
     procedure, public :: write => write_theta_grid_eik_config
     procedure, public :: reset => reset_theta_grid_eik_config
     procedure, public :: broadcast => broadcast_theta_grid_eik_config
     procedure, public, nopass :: get_default_name => get_default_name_theta_grid_eik_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_theta_grid_eik_config
  end type theta_grid_eik_config_type

  type(theta_grid_eik_config_type) :: theta_grid_eik_config

contains
  !> FIXME : Add documentation
  subroutine wnml_theta_grid_eik(unit)
    use geometry, only: surf, alpha_input, beta_prime_input, invLp_input
    use geometry, only: dp_mult, bishop, irho, force_sym, use_large_aspect
    use geometry, only: gen_eq, efit_eq, ppl_eq, local_eq, dfit_eq, idfit_eq
    use geometry, only: gs2d_eq, chs_eq, transp_eq, writelots, equal_arc, eqfile, eqnormfile
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_eik_knobs"
    write (unit, fmt="(' irho =  ',i2)") irho
    write (unit, fmt="(' ppl_eq =   ',L1)") ppl_eq
    write (unit, fmt="(' efit_eq =  ',L1)") efit_eq
    write (unit, fmt="(' gen_eq =   ',L1)") gen_eq
    write (unit, fmt="(' dfit_eq =  ',L1)") dfit_eq
    write (unit, fmt="(' idfit_eq = ',L1)") idfit_eq
    write (unit, fmt="(' local_eq =  ',L1)") local_eq
    write (unit, fmt="(' transp_eq =  ',L1)") transp_eq
    write (unit, fmt="(' gs2d_eq =  ',L1)") gs2d_eq
    write (unit, fmt="(' chs_eq =  ',L1)") chs_eq
    write (unit, fmt="(' equal_arc =  ',L1)") equal_arc
    write (unit, fmt="(' bishop =  ',i2)") bishop
    write (unit, fmt="(' s_hat_input =  ',e13.6)") surf%shat
    write (unit, fmt="(' alpha_input =  ',e13.6)") alpha_input
    write (unit, fmt="(' invLp_input =  ',e13.6)") invLp_input
    write (unit, fmt="(' beta_prime_input =  ',e13.6)") beta_prime_input
    write (unit, fmt="(' dp_mult =  ',e13.6)") dp_mult
    write (unit, fmt="(' delrho =  ',e13.6)") surf%dr
    write (unit, fmt="(' force_sym =  ',L1)") force_sym
    write (unit, fmt="(' use_large_aspect =  ',L1)") use_large_aspect
    write (unit, fmt="(' writelots =  ',L1)") writelots
    write (unit, fmt="(' eqfile = ',a)") '"'//trim(eqfile)//'"'
    write (unit, fmt="(' eqnormfile = ',a)") '"'//trim(eqnormfile)//'"'
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_eik

  !> FIXME : Add documentation
  subroutine check_theta_grid_eik(report_unit, dbetadrho)
    use theta_grid_params, only: eps
    use geometry, only: surf, dp_mult, alpha_input, beta_prime_input, invLp_input
    use geometry, only: bishop, irho, eqfile, eqnormfile
    use geometry, only: idfit_eq, gen_eq, efit_eq, ppl_eq, local_eq, dfit_eq, chs_eq
    use geo_utils, only: geo_type_miller, geo_type_global, &
         geo_type_generalized_elongation, geo_type_fourier_series, &
         geo_type_miller_extended_harmonic
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: dbetadrho
    real :: beta_prime_new, shat, rhoc
    beta_prime_new = eikcoefs_results%dbetadrho
    shat = eikcoefs_results%shat
    rhoc = eikcoefs_results%rhoc
    call checklogic_theta_grid_eik(report_unit)
    write (report_unit, *)
    if (local_eq) then
       write (report_unit, fmt="('The following local equilibrium model has been selected:')")
       select case (surf%geoType)
          case (geo_type_miller)
             write (report_unit, fmt="('Miller')")
          case (geo_type_global)
             write (report_unit, fmt="('Global')")
          case (geo_type_generalized_elongation)
             write (report_unit, fmt="('GeneralizedEllipticity')")
          case (geo_type_fourier_series)
             write (report_unit, fmt="('FourierSeries')")
          case (geo_type_miller_extended_harmonic)
             write (report_unit, fmt="('MillerExtendedHarmonic')")
          case default
             write (report_unit, fmt="('ERROR: invalid analytic geometry specification')")
       end select
       if (surf%rmaj == 1.0) then
          write (report_unit, &
               & fmt="('Scale lengths are normalized to the major radius, R')")
       else
          write (report_unit, fmt="('The aspect ratio R/a = ',f7.4)") surf%rmaj
       end if
       if (surf%rmaj /= surf%r_geo) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('R_geo is not equal to Rmaj.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
       end if
       if (irho /= 2) then
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You have selected irho = ',i2)") irho
          write (report_unit, fmt="('For local equilibria, irho=2 is required.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
       end if
       write (report_unit, *)
       write (report_unit, fmt="('The safety factor q =      ',f7.4)") surf%q
       eps = rhoc / surf%R_geo
       write (report_unit, fmt="('and epsilon == r/R =       ',f7.4)") eps
       write (report_unit, *)
       if (eps > epsilon(0.0)) then
          write (report_unit, fmt="('Trapped particles are included.')")
       else
          write (report_unit, fmt="('Trapped particles are neglected.')")
       end if
       write (report_unit, *)
       write (report_unit, fmt="('B_poloidal is determined by:')")
       if (surf%geoType == geo_type_global) then
          write (report_unit, *)
          write (report_unit, fmt="('    minor radius of shaped surface: aSurf = ',f7.4)") surf%aSurf
       elseif (surf%geoType == geo_type_generalized_elongation .or. surf%geoType == geo_type_fourier_series) then
          write (report_unit, *)
          write (report_unit, fmt="('    first shaping mode number: mMode =  ',i2)") surf%mMode
          write (report_unit, fmt="('    second shaping mode number: nMode = ',i2)") surf%nMode
       end if
       write (report_unit, *)
       write (report_unit, fmt="('    deltam =       ',f7.4)") surf%delm
       write (report_unit, fmt="('    deltan =        ',f7.4)") surf%deln
       if (surf%geoType /= geo_type_global) then
          write (report_unit, *)
          write (report_unit, fmt="('  & gradient: d deltam /d rho =   ',f7.4)") surf%delmp
          write (report_unit, fmt="('  & gradient: d deltan /d rho = ',f7.4)") surf%delnp
       end if
       write (report_unit, *)
       write (report_unit, fmt="('    thetam =       ',f7.4)") surf%thm
       write (report_unit, fmt="('    thetan =        ',f7.4)") surf%thn
       write (report_unit, *)
       write (report_unit, fmt="('    shift =        ',f7.4)") surf%sHorz
       write (report_unit, fmt="('    shiftVert =    ',f7.4)") surf%sVert

       write (report_unit, *)
       write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") shat
       write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')")
       if (abs(shat) <= 1.e-5) then
          write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
       end if
       select case (bishop)
       case (3)
          write (report_unit, fmt="('The normalized inverse pressure gradient scale length = ',f8.4)") invLp_input
       case (4)
          write (report_unit, fmt="('The beta gradient d beta / d rho = ',f8.4)") beta_prime_input
          if (beta_prime_input > epsilon(0.0)) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('beta_prime > 0.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
          if (abs(beta_prime_input - dbetadrho) > 1.e-2*max(abs(beta_prime_input),abs(dbetadrho))) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('beta_prime_input is not consistent with beta and Lp.')")
             write (report_unit, fmt="('beta_prime_input=',f6.4,' dbeta_drho (from species)=',f6.4)") beta_prime_input, dbetadrho
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
       case (5)
          write (report_unit, fmt="('The alpha parameter (R beta_prime q**2) = ',f8.4)") alpha_input
          !              write (*,*) alpha_input, dbetadrho, qinp, Rmaj
          if (abs(alpha_input + dbetadrho*surf%q**2*surf%rmaj) > 1.e-2) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('alpha is not consistent with beta, q, and Lp.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
       case default
          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You have selected bishop = ',i2)") bishop
          write (report_unit, fmt="('For local equilibria, bishop = 4 is recommended.')")
          if (bishop == 1) then
             write (report_unit, fmt="('For d beta / d rho = 0, bishop = 1 is ok.')")
             write (report_unit, fmt="('Otherwise, ')")
          end if
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)
       end select
    end if
    if (.not. local_eq) then
       if (gen_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Equilibrium information obtained from NetCDF file:')")
          write (report_unit, fmt="(a)") trim(eqfile)
       end if
       if (ppl_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Equilibrium information obtained from NetCDF file:')")
          write (report_unit, fmt="(a)") trim(eqfile)
       end if
       if (dfit_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Dipole equilibrium information obtained from file:')")
          write (report_unit, fmt="(a)") trim(eqfile)
          write (report_unit, fmt="(a)") trim(eqnormfile)
       end if
       if (idfit_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Dipole equilibrium information obtained from file:')")
          write (report_unit, fmt="(a)") trim(eqfile)
       end if
       if (efit_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Equilibrium information obtained from eqdsk:')")
          write (report_unit, fmt="(a)") trim(eqfile)
       end if
       if (chs_eq) then
          write (report_unit, *)
          write (report_unit, fmt="('Equilibrium information obtained from chease:')")
          write (report_unit, fmt="(a)") trim(eqfile)
       end if
       select case (bishop)
       case (1)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=1, so dp/drho and s_hat will be found from the equilibrium file.')")
          write (report_unit, *)
       case (3)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=3.')")
          write (report_unit, *)
          write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") surf%shat
          write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')")
          if (abs(shat) <= 1.e-5) then
             write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
          end if
          write (report_unit, fmt="('The normalized inverse pressure gradient scale length = ',f8.4)") invLp_input
       case (4)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=4.')")
          write (report_unit, *)
          write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") surf%shat
          write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')")
          if (abs(shat) <= 1.e-5) then
             write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
          end if
          write (report_unit, fmt="('The beta gradient d beta / d rho = ',f8.4)") beta_prime_input
          if (beta_prime_input > epsilon(0.0)) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('beta_prime > 0.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
          if (abs(beta_prime_input - dbetadrho) > 1.e-2*abs(dbetadrho)) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('beta_prime_input is not consistent with beta and Lp.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
       case (5)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=5.')")
          write (report_unit, *)
          write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") surf%shat
          write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')")
          if (abs(shat) <= 1.e-5) then
             write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
          end if
          write (report_unit, fmt="('The alpha parameter (R beta_prime q**2) = ',f8.4)") alpha_input
          write (*,*) alpha_input, dbetadrho, surf%q, surf%rmaj
          if (abs(alpha_input + dbetadrho*surf%q**2*surf%rmaj) > 1.e-2) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('alpha is not consistent with beta, q, and Lp.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if
       case (6)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=6.')")
          write (report_unit, *)
          write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") surf%shat
          write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')")
          if (abs(shat) <= 1.e-5) then
             write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
          end if
          write (report_unit, fmt="('The value of dp/drho will be found from the equilibrium file.')")
       case (7)
          write (report_unit, *)
          write (report_unit, fmt="('You have set bishop=7.')")
          write (report_unit, fmt="('The value of s_hat will be found from the equilibrium file.')")
          write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") shat
          write (report_unit, fmt="('The value of dp/drho found from the equilibrium file will be multiplied by',f10.4)") dp_mult
          write (report_unit, fmt="('to give beta gradient d beta / d rho = ',f8.4)") beta_prime_new

          if (abs(beta_prime_new - dbetadrho) > 1.e-2*abs(dbetadrho)) then
             write (report_unit, *)
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('beta_prime_new is not consistent with beta and Lp.')")
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *)
          end if

       case default

          write (report_unit, *)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You have selected a value for bishop that is not recommended.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *)

       end select
    end if
  end subroutine check_theta_grid_eik

  !> FIXME : Add documentation
  subroutine checklogic_theta_grid_eik(report_unit)
    use geometry, only: geq=>gen_eq, eeq=>efit_eq, peq=>ppl_eq
    use geometry, only: leq=>local_eq, deq=>dfit_eq, ceq=>chs_eq
    integer, intent (in) :: report_unit

    if(geq .and. deq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing gen_eq = .true. AND dfit_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

    if(geq .and. eeq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing gen_eq = .true. AND efit_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

     if(geq .and. ceq) then
        write (report_unit, *)
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing gen_eq = .true. AND chs_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *)
     endif

    if(geq .and. peq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing gen_eq = .true. AND ppl_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

    if(geq .and. leq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing gen_eq = .true. AND local_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

    if(eeq .and. deq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing efit_eq = .true. AND dfit_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

    if(eeq .and. leq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing efit_eq = .true. AND local_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

     if(eeq .and. ceq) then
        write (report_unit, *)
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing efit_eq = .true. AND chs_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *)
     endif

    if(eeq .and. peq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing efit_eq = .true. AND ppl_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

    if(deq .and. leq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing dfit_eq = .true. AND local_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

     if(deq .and. ceq) then
        write (report_unit, *)
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing dfit_eq = .true. AND chs_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *)
     endif

    if(deq .and. peq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing dfit_eq = .true. AND ppl_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif

     if(ceq .and. peq) then
        write (report_unit, *)
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing chs_eq = .true. AND ppl_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *)
     endif

     if(ceq .and. leq) then
        write (report_unit, *)
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing chs_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *)
     endif

    if(peq .and. leq) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write(report_unit,fmt="('Choosing ppl_eq = .true. AND local_eq = .true. is not permitted.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    endif
  end subroutine checklogic_theta_grid_eik

  !> FIXME : Add documentation
  subroutine init_theta_grid_eik(theta_grid_eik_config_in)
    use geometry, only: eikcoefs, surf, use_large_aspect, verb_geo => verb
    use geometry, only: run_eikcoefs_and_resample
    use unit_tests, only: job_id
    use theta_grid_params, only: init_theta_grid_params, ntheta, nperiod
    use runtime_tests, only: verbosity
    use unit_tests, only: debug_message
    use mp, only: proc0
    implicit none
    type(theta_grid_eik_config_type), intent(in), optional :: theta_grid_eik_config_in
    real :: rhoc_save
    integer, parameter :: verb = 3
    character(4) :: ntheta_char
    verb_geo = verbosity()

    if (initialized) return
    initialized = .true.
    write(ntheta_char, "(I4)") ntheta
    call debug_message(verb, "init_theta_grid_eik: call init_theta_grid_params, ntheta="//ntheta_char)
    ! After this call, would think you have ntheta from input file
    ! stored in theta_grid_params data structure.
    ! but when running from numerical equilibrium, this is not right
    ! Instead, get it stored via the eikcoefs call below.
    call init_theta_grid_params

    write(ntheta_char, "(I4)") ntheta
    call debug_message(verb, "init_theta_grid_eik: call read_parameters, ntheta="//ntheta_char)
    call read_parameters(theta_grid_eik_config_in)

    rhoc_save = surf%r
    if (use_large_aspect) surf%r = 1.5 * surf%dr

    if(proc0) then
       if (ntheta_geometry == -1) then
          write(ntheta_char, "(I4)") ntheta
          call debug_message(verb, "init_theta_grid_eik: call eikcoefs, ntheta="//ntheta_char)
          call eikcoefs(ntheta, nperiod, eikcoefs_results, job_id)
       else
          write(ntheta_char, "(I4)") ntheta_geometry
          call debug_message(verb, "init_theta_grid_eik: call run_eikcoefs_and_resample, ntheta_geometry="//ntheta_char)
          call run_eikcoefs_and_resample (ntheta_geometry, ntheta, nperiod, eikcoefs_results, job_id)
       end if
    end if

    write(ntheta_char, "(I4)") ntheta
    call debug_message(verb, "init_theta_grid_eik: done, ntheta="//ntheta_char)

    surf%r = rhoc_save
  end subroutine init_theta_grid_eik

  subroutine finish_theta_grid_eik
    use geometry, only: finish_geometry

    initialized = .false.
    call finish_geometry
    call theta_grid_eik_config%reset()
  end subroutine finish_theta_grid_eik

  !> FIXME : Add documentation
  subroutine eik_get_sizes (nthetaout, nperiodout, nbsetout)
    use theta_grid_params, only: ntheta, nperiod
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound
  end subroutine eik_get_sizes

  !> FIXME : Add documentation
  subroutine eik_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag,&
       gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gds2, gds21, gds22, gds23, gds24, gds24_noq, &
       grho, Rplot, Zplot, Rprime, Zprime, aplot, aprime, shat, drhodpsi,&
       kxfac, qval, gb_to_cv, Bpol, surfarea, dvdrhon, rhoc)
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    use unit_tests, only: debug_message
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    real, intent (out) :: shat, drhodpsi, kxfac, qval, surfarea, dvdrhon, rhoc
    logical, intent (in) :: gb_to_cv
    integer, parameter :: verb=3
    character(4) :: ntgrid_char
    write(ntgrid_char, "(I4)") ntgrid
    call debug_message(verb, 'eik_get_grids: ntgrid= '//ntgrid_char)
    theta(-ntgrid:ntgrid) = eikcoefs_results%theta(-ntgrid:ntgrid)
    gradpar(-ntgrid:ntgrid) = eikcoefs_results%gradpar(-ntgrid:ntgrid)
    bmag(-ntgrid:ntgrid) = eikcoefs_results%bmag(-ntgrid:ntgrid)
    cvdrift(-ntgrid:ntgrid) = eikcoefs_results%cvdrift(-ntgrid:ntgrid)
    cvdrift0(-ntgrid:ntgrid) = eikcoefs_results%cvdrift0(-ntgrid:ntgrid)
    gbdrift(-ntgrid:ntgrid) = eikcoefs_results%gbdrift(-ntgrid:ntgrid)
    gbdrift0(-ntgrid:ntgrid) = eikcoefs_results%gbdrift0(-ntgrid:ntgrid)
    cdrift(-ntgrid:ntgrid) = eikcoefs_results%cdrift(-ntgrid:ntgrid)
    cdrift0(-ntgrid:ntgrid) = eikcoefs_results%cdrift0(-ntgrid:ntgrid)
    gds2(-ntgrid:ntgrid) = eikcoefs_results%gds2(-ntgrid:ntgrid)
    gds21(-ntgrid:ntgrid) = eikcoefs_results%gds21(-ntgrid:ntgrid)
    gds22(-ntgrid:ntgrid) = eikcoefs_results%gds22(-ntgrid:ntgrid)
    gds23(-ntgrid:ntgrid) = eikcoefs_results%gds23(-ntgrid:ntgrid)
    gds24(-ntgrid:ntgrid) = eikcoefs_results%gds24(-ntgrid:ntgrid)
    gds24_noq(-ntgrid:ntgrid) = eikcoefs_results%gds24_noq(-ntgrid:ntgrid)
    grho(-ntgrid:ntgrid) = eikcoefs_results%grho(-ntgrid:ntgrid)
    Rplot(-ntgrid:ntgrid) = eikcoefs_results%Rplot(-ntgrid:ntgrid)
    Zplot(-ntgrid:ntgrid) = eikcoefs_results%Zplot(-ntgrid:ntgrid)
    aplot(-ntgrid:ntgrid) = eikcoefs_results%aplot(-ntgrid:ntgrid)
    Rprime(-ntgrid:ntgrid) = eikcoefs_results%Rprime(-ntgrid:ntgrid)
    Zprime(-ntgrid:ntgrid) = eikcoefs_results%Zprime(-ntgrid:ntgrid)
    aprime(-ntgrid:ntgrid) = eikcoefs_results%aprime(-ntgrid:ntgrid)
    Bpol(-ntgrid:ntgrid) = eikcoefs_results%Bpol(-ntgrid:ntgrid)

    if (gb_to_cv) then
       gbdrift(-ntgrid:ntgrid) = cvdrift(-ntgrid:ntgrid)
       gbdrift0(-ntgrid:ntgrid) = cvdrift0(-ntgrid:ntgrid)
    end if

    call debug_message(verb, 'eik_get_grids: call theta_grid_gridgen_init')
    call theta_grid_gridgen_init
    call debug_message(verb, 'eik_get_grids: call gridgen_get_grids')
    call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
         theta, bset, bmag, &
         gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gds2, gds21, gds22, gds23, gds24, gds24_noq, &
         grho, Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    shat = eikcoefs_results%shat
    drhodpsi = eikcoefs_results%drhodpsi
    kxfac = eikcoefs_results%kxfac
    qval = eikcoefs_results%qsf
    surfarea = eikcoefs_results%surfarea
    dvdrhon = eikcoefs_results%dvdrhon
    rhoc = eikcoefs_results%rhoc
    call debug_message(verb, 'eik_get_grids: end')
  end subroutine eik_get_grids

  !> FIXME : Add documentation
  subroutine read_parameters(theta_grid_eik_config_in)
    use mp, only: proc0
    use file_utils, only: input_unit, input_unit_exist
    use geometry, only: surf, irho, eqnormfile
    use geometry, only: ppl_eq, gen_eq, efit_eq, eqfile, local_eq, dfit_eq, gs2d_eq
    use geometry, only: equal_arc, transp_eq, idfit_eq, chs_eq, bishop
    use geometry, only: alpha_input, invLp_input, beta_prime_input, dp_mult
    use geometry, only: writelots, force_sym, use_large_aspect
    use theta_grid_params, only: get_parameters_as_surf
    use theta_grid_params, only: shat_in => shat
    use theta_grid_params, only: alpmhd_in => alpmhd
    use theta_grid_params, only: betaprim_in => betaprim
    implicit none
    type(theta_grid_eik_config_type), intent(in), optional :: theta_grid_eik_config_in
    integer :: isym, iflux, itor
    ! Copy some values from theta_grid_params into this module (actually into geometry)
    surf = get_parameters_as_surf()
    isym = -1 ; iflux = -1 ; itor = -1
    if (present(theta_grid_eik_config_in)) theta_grid_eik_config = theta_grid_eik_config_in

    ! Smart defaults
    if (.not.theta_grid_eik_config%is_initialised()) then
       theta_grid_eik_config%alpha_input = alpmhd_in
       theta_grid_eik_config%s_hat_input = shat_in
       theta_grid_eik_config%beta_prime_input = betaprim_in
    end if

    call theta_grid_eik_config%init(name = 'theta_grid_eik_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    alpha_input = theta_grid_eik_config%alpha_input
    beta_prime_input = theta_grid_eik_config%beta_prime_input
    bishop = theta_grid_eik_config%bishop
    chs_eq = theta_grid_eik_config%chs_eq
    surf%dr = theta_grid_eik_config%delrho
    dfit_eq = theta_grid_eik_config%dfit_eq
    dp_mult = theta_grid_eik_config%dp_mult
    efit_eq = theta_grid_eik_config%efit_eq
    eqfile = theta_grid_eik_config%eqfile
    eqnormfile = theta_grid_eik_config%eqnormfile
    equal_arc = theta_grid_eik_config%equal_arc
    force_sym = theta_grid_eik_config%force_sym
    gen_eq = theta_grid_eik_config%gen_eq
    gs2d_eq = theta_grid_eik_config%gs2d_eq
    idfit_eq = theta_grid_eik_config%idfit_eq
    iflux = theta_grid_eik_config%iflux
    invlp_input = theta_grid_eik_config%invlp_input
    irho = theta_grid_eik_config%irho
    isym = theta_grid_eik_config%isym
    itor = theta_grid_eik_config%itor
    local_eq = theta_grid_eik_config%local_eq
    ntheta_geometry = theta_grid_eik_config%ntheta_geometry
    ppl_eq = theta_grid_eik_config%ppl_eq
    surf%shat = theta_grid_eik_config%s_hat_input
    transp_eq = theta_grid_eik_config%transp_eq
    use_large_aspect = theta_grid_eik_config%use_large_aspect
    writelots = theta_grid_eik_config%writelots

    if (proc0) then
       if (isym /= -1) write(6,*) 'Warning: isym input is deprecated -- set force_sym instead'
       if (iflux /= -1) write(6,*) 'Warning: iflux input is deprecated'
       if (itor /= -1) write(6,*) 'Warning: itor input is deprecated -- set use_large_aspect instead'
    end if

    exist = theta_grid_eik_config%exist
  end subroutine read_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_theta_grid_eik_config(theta_grid_eik_config_in)
    use mp, only: mp_abort
    type(theta_grid_eik_config_type), intent(in), optional :: theta_grid_eik_config_in
    if (initialized) then
       call mp_abort("Trying to set theta_grid_eik_config when already initialized.", to_screen = .true.)
    end if
    if (present(theta_grid_eik_config_in)) then
       theta_grid_eik_config = theta_grid_eik_config_in
    end if
  end subroutine set_theta_grid_eik_config

  !> Get the module level config instance
  function get_theta_grid_eik_config()
    type(theta_grid_eik_config_type) :: get_theta_grid_eik_config
    get_theta_grid_eik_config = theta_grid_eik_config
  end function get_theta_grid_eik_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------

  !> Reads in the theta_grid_eik_knobs namelist and populates the member variables
  subroutine read_theta_grid_eik_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(theta_grid_eik_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = EQFILE_LENGTH) :: eqfile, eqnormfile
    integer :: bishop, iflux, irho, isym, itor, ntheta_geometry
    logical :: chs_eq, dfit_eq, efit_eq, equal_arc, force_sym, gen_eq, gs2d_eq, idfit_eq, local_eq, ppl_eq, transp_eq, use_large_aspect, writelots
    real :: alpha_input, beta_prime_input, delrho, dp_mult, invlp_input, s_hat_input

    namelist /theta_grid_eik_knobs/ alpha_input, beta_prime_input, bishop, chs_eq, delrho, dfit_eq, dp_mult, efit_eq, eqfile, eqnormfile, equal_arc, force_sym, gen_eq, gs2d_eq, idfit_eq, iflux, invlp_input, irho, isym, itor, local_eq, &
         ntheta_geometry, ppl_eq, s_hat_input, transp_eq, use_large_aspect, writelots

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    alpha_input = self%alpha_input
    beta_prime_input = self%beta_prime_input
    bishop = self%bishop
    chs_eq = self%chs_eq
    delrho = self%delrho
    dfit_eq = self%dfit_eq
    dp_mult = self%dp_mult
    efit_eq = self%efit_eq
    eqfile = self%eqfile
    eqnormfile = self%eqnormfile
    equal_arc = self%equal_arc
    force_sym = self%force_sym
    gen_eq = self%gen_eq
    gs2d_eq = self%gs2d_eq
    idfit_eq = self%idfit_eq
    iflux = self%iflux
    invlp_input = self%invlp_input
    irho = self%irho
    isym = self%isym
    itor = self%itor
    local_eq = self%local_eq
    ntheta_geometry = self%ntheta_geometry
    ppl_eq = self%ppl_eq
    s_hat_input = self%s_hat_input
    transp_eq = self%transp_eq
    use_large_aspect = self%use_large_aspect
    writelots = self%writelots

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = theta_grid_eik_knobs)

    ! Now copy from local variables into type members
    self%alpha_input = alpha_input
    self%beta_prime_input = beta_prime_input
    self%bishop = bishop
    self%chs_eq = chs_eq
    self%delrho = delrho
    self%dfit_eq = dfit_eq
    self%dp_mult = dp_mult
    self%efit_eq = efit_eq
    self%eqfile = eqfile
    self%eqnormfile = eqnormfile
    self%equal_arc = equal_arc
    self%force_sym = force_sym
    self%gen_eq = gen_eq
    self%gs2d_eq = gs2d_eq
    self%idfit_eq = idfit_eq
    self%iflux = iflux
    self%invlp_input = invlp_input
    self%irho = irho
    self%isym = isym
    self%itor = itor
    self%local_eq = local_eq
    self%ntheta_geometry = ntheta_geometry
    self%ppl_eq = ppl_eq
    self%s_hat_input = s_hat_input
    self%transp_eq = transp_eq
    self%use_large_aspect = use_large_aspect
    self%writelots = writelots

    self%exist = exist
  end subroutine read_theta_grid_eik_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_theta_grid_eik_config(self, unit)
    implicit none
    class(theta_grid_eik_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("alpha_input", self%alpha_input, unit_internal)
    call self%write_key_val("beta_prime_input", self%beta_prime_input, unit_internal)
    call self%write_key_val("bishop", self%bishop, unit_internal)
    call self%write_key_val("chs_eq", self%chs_eq, unit_internal)
    call self%write_key_val("delrho", self%delrho, unit_internal)
    call self%write_key_val("dfit_eq", self%dfit_eq, unit_internal)
    call self%write_key_val("dp_mult", self%dp_mult, unit_internal)
    call self%write_key_val("efit_eq", self%efit_eq, unit_internal)
    call self%write_key_val("eqfile", self%eqfile, unit_internal)
    call self%write_key_val("eqnormfile", self%eqnormfile, unit_internal)
    call self%write_key_val("equal_arc", self%equal_arc, unit_internal)
    call self%write_key_val("force_sym", self%force_sym, unit_internal)
    call self%write_key_val("gen_eq", self%gen_eq, unit_internal)
    call self%write_key_val("gs2d_eq", self%gs2d_eq, unit_internal)
    call self%write_key_val("idfit_eq", self%idfit_eq, unit_internal)
    call self%write_key_val("iflux", self%iflux, unit_internal)
    call self%write_key_val("invlp_input", self%invlp_input, unit_internal)
    call self%write_key_val("irho", self%irho, unit_internal)
    call self%write_key_val("isym", self%isym, unit_internal)
    call self%write_key_val("itor", self%itor, unit_internal)
    call self%write_key_val("local_eq", self%local_eq, unit_internal)
    call self%write_key_val("ntheta_geometry", self%ntheta_geometry, unit_internal)
    call self%write_key_val("ppl_eq", self%ppl_eq, unit_internal)
    call self%write_key_val("s_hat_input", self%s_hat_input, unit_internal)
    call self%write_key_val("transp_eq", self%transp_eq, unit_internal)
    call self%write_key_val("use_large_aspect", self%use_large_aspect, unit_internal)
    call self%write_key_val("writelots", self%writelots, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_theta_grid_eik_config

  !> Resets the config object to the initial empty state
  subroutine reset_theta_grid_eik_config(self)
    class(theta_grid_eik_config_type), intent(in out) :: self
    type(theta_grid_eik_config_type) :: empty
    select type (self)
    type is (theta_grid_eik_config_type)
       self = empty
    end select
  end subroutine reset_theta_grid_eik_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_theta_grid_eik_config(self)
    use mp, only: broadcast
    implicit none
    class(theta_grid_eik_config_type), intent(in out) :: self
    call broadcast(self%alpha_input)
    call broadcast(self%beta_prime_input)
    call broadcast(self%bishop)
    call broadcast(self%chs_eq)
    call broadcast(self%delrho)
    call broadcast(self%dfit_eq)
    call broadcast(self%dp_mult)
    call broadcast(self%efit_eq)
    call broadcast(self%eqfile)
    call broadcast(self%eqnormfile)
    call broadcast(self%equal_arc)
    call broadcast(self%force_sym)
    call broadcast(self%gen_eq)
    call broadcast(self%gs2d_eq)
    call broadcast(self%idfit_eq)
    call broadcast(self%iflux)
    call broadcast(self%invlp_input)
    call broadcast(self%irho)
    call broadcast(self%isym)
    call broadcast(self%itor)
    call broadcast(self%local_eq)
    call broadcast(self%ntheta_geometry)
    call broadcast(self%ppl_eq)
    call broadcast(self%s_hat_input)
    call broadcast(self%transp_eq)
    call broadcast(self%use_large_aspect)
    call broadcast(self%writelots)

    call broadcast(self%exist)
  end subroutine broadcast_theta_grid_eik_config

  !> Gets the default name for this namelist
  function get_default_name_theta_grid_eik_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_eik_config
    get_default_name_theta_grid_eik_config = "theta_grid_eik_knobs"
  end function get_default_name_theta_grid_eik_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_theta_grid_eik_config()
    implicit none
    logical :: get_default_requires_index_theta_grid_eik_config
    get_default_requires_index_theta_grid_eik_config = .false.
  end function get_default_requires_index_theta_grid_eik_config
end module theta_grid_eik
