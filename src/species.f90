!> FIXME : Add documentation  
module species
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private
  public :: init_species, finish_species

  public :: write_trinity_parameters
  public :: wnml_species, check_species
  public :: nspec, specie, spec
  public :: ion_species, electron_species, tracer_species
  public :: f0_maxwellian, f0_tabulated, f0_sdanalytic
  public :: has_ion_species, has_electron_species, has_nonmaxw_species
  public :: calculate_slowingdown_parameters
  public :: ions, electrons, impurity
  public :: set_overrides
  public :: determine_species_order
  public :: calculate_f0_arrays, set_current_f0_species, eval_f0
  public :: nonmaxw_corr, f0_values, dlnf0drho

  public :: species_config_type
  public :: species_element_config_type
  public :: set_species_config
  public :: get_species_config
  public :: get_species_element_config

  !> FIXME : Add documentation 
  type :: specie
     real :: z
     real :: mass
     real :: dens, dens0, u0
     real :: tpar0,tperp0
     real :: temp
     real :: tprim
     real :: fprim
     real :: uprim, uprim2
     real :: vnewk, nustar
     real :: nu, nu_h  ! nu will be the preferred collisionality parameter moving forward
                       ! as it will have a consistent v_t scaling
                       ! nu_h controls a hyperviscous term embedded in the collision operator
     real :: stm, zstm, tz, smz, zt
     real :: bess_fac !Artificial factor multiplying the Bessel function argument
     real :: vcrit
     real :: vcprim
     integer :: type
     integer :: f0type
  end type specie

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: tracer_species = 3 ! for test particle diffusion studies

  integer, parameter :: f0_maxwellian = 1, &           !< Regular Maxwellian species
                        f0_tabulated = 2, &            !< f0 determined from table in external input file
                        f0_sdanalytic = 3              !< Analytic Gaffey slowing-down distribution

  ! FIXME
  !> Various quantities needed in the calculation of the slowing-down distribution in case electrons are
  !! adiabatic
  real:: me, ZI_fac

  !> Defines a working species so that eval_f0 only needs to take one argument (needed for genquad)
  integer:: is_global = -1

  integer :: nspec
  type (specie), dimension (:), allocatable :: spec

  integer :: ions, electrons, impurity
  real, dimension (:,:), allocatable :: nonmaxw_corr, f0_values, dlnf0drho

  logical :: initialized = .false.
  logical :: exist

  !> Used to represent the input configuration of species
  type, extends(abstract_config_type) :: species_config_type
     ! namelist : species_knobs
     ! indexed : false
     !> Specifies the electron mass used in calculating the slowing
     !> down parameters in simulations without kinetic electrons.
     real :: me = -1.0
     !> Number of kinetic species to evolve.
     integer :: nspec = 2
     !> Specify the ion charge factor used in calculating the slowing
     !> down parameters. If not specified that will attempt to
     !> calculate from the kinetic ion species.
     real :: zi_fac = -1.0
   contains
     procedure, public :: read => read_species_config
     procedure, public :: write => write_species_config
     procedure, public :: reset => reset_species_config
     procedure, public :: broadcast => broadcast_species_config
     procedure, public, nopass :: get_default_name => get_default_name_species_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_species_config
  end type species_config_type
  
  type(species_config_type) :: species_config

  !> Used to represent the input configuration of species
  type, extends(abstract_config_type) :: species_element_config_type
     ! namelist : species_parameters
     ! indexed : true
     !> Multiplies the argument of the Bessel function for this
     !> species. Useful for exploring the effect of the Bessel
     !> functions.
     real :: bess_fac = 1.0
     !> Set the normalised density for this species
     real :: dens = 1.0
     !> Parameter used for a few specific initial condition types.
     !>
     !> @todo Consider moving this parameter to [[init_g_knobs]].
     real :: dens0 = 1.0
     !> Select which type of background distribution should be assumed
     !> for this species. Must be one of:
     !>
     !> - 'maxwellian' : Standard maxwellian
     !> - 'tabulated' : A general tabulated form
     !> - 'sdanalytic' : A analytic slowing down form.
     !>
     !> See [G. Wilkie
     !> thesis](https://drum.lib.umd.edu/handle/1903/17302).
     character(len = 20) :: f0type = "maxwellian"
     !> Normalised inverse density gradient: \(-\frac{1}{n}\frac{dn}{d
     !> \rho_N}\) (Note here \(\rho_N\) is the normalised radius
     !> \(\frac{\rho}{L_\textrm{ref}}\).
     real :: fprim = 2.2
     !> Normalised mass of this species
     real :: mass = 1.0
     !> Not used.
     !>
     !> @todo remove
     real :: nu = -1.0
     !> Set the species hyper collision frequency used in the pitch
     !> angle scattering collision operator if
     !> [[collision_knobs::hypermult]] is `true`.
     real :: nu_h = 0.0
     !> Not used.
     !>
     !> @todo remove
     real :: nustar = -1.0
     !> Set the normalised species temperature
     real :: temp = 1.0
     !> Parameter used for a few specific initial condition types.
     !>
     !> @todo Consider moving this parameter to [[init_g_knobs]].
     real :: tpar0 = 0.
     !> Parameter used for a few specific initial condition types.
     !>
     !> @todo Consider moving this parameter to [[init_g_knobs]].
     real :: tperp0 = 0.
     !> Normalised inverse density gradient: \(-\frac{1}{T}\frac{dT}{d
     !> \rho_N}\) (Note here \(\rho_N\) is the normalised radius
     !> \(\frac{\rho}{L_\textrm{ref}}\).
     real :: tprim = 6.9
     !> Sets the characterisation of the species. Should be one of
     !>
     !> - 'ion' : Thermal ion species
     !> - 'default' :  Same as 'ion'
     !> - 'electron' : Thermal electron species
     !> - 'e' : Same as 'electron'
     !> - 'trace' : Tracer species
     !>
     !> This primarily controls the treatment in the collision
     !> operator, the detection of adiabatic species. Tracer species
     !> do not contribute to velocity space integrals.
     character(len = 20) :: type = "default"
     !> Parameter used for a few specific initial condition types.
     !>
     !> @todo Consider moving this parameter to [[init_g_knobs]].
     real :: u0 = 1.0
     !> Controls an energy independent part of a source term.
     !>
     !> @note Documentation needs improving, looks like it's an early
     !> attempt at adding in the parallel flow shear drive term
     real :: uprim = 0.0
     !> Controls an energy dependent part of a source term.
     !>
     !> @note Documentation needs improving, looks like it's an early
     !> attempt at adding in the parallel flow shear drive term
     real :: uprim2 = 0.0
     !> Part of the `sdanalytic` slowing down model. See [G. Wilkie
     !> thesis](https://drum.lib.umd.edu/handle/1903/17302).
     real :: vcprim = -999.9
     !> Part of the `sdanalytic` slowing down model. See [G. Wilkie
     !> thesis](https://drum.lib.umd.edu/handle/1903/17302).
     real :: vcrit = -1.0
     !> Sets the normalised species-species collisionality
     !> frequency. For species s, `vnewk = nu_ss Lref/vref` where Lref
     !> is the reference length, `vref` is the reference thermal speed
     !> (not the thermal speed of species s!), \(\nu_{ss}=\sqrt{2}\pi
     !> n_s Z_s^4 e^4 \log(\lambda)/\sqrt{m_s}T^{3/2}_s\) is a
     !> dimensional collision frequency, e is the proton charge,
     !> \(\log(\lambda)\) is the Coloumb logarithm, and
     !> \(Z_s,T_s,n_s,m_s\) are the (charge, temperature, density,
     !> mass) of species s. (The dimensional collision frequency given
     !> here is in Gaussian units. For SI units, include an extra
     !> factor \((4\pi\epsilon_0)^2\) in the denominator of the
     !> definition of \(\nu_{ss}\).
     real :: vnewk = 0.0
     !> Normalised species charge
     real :: z = 1.0
   contains
     procedure, public :: read => read_species_element_config
     procedure, public :: write => write_species_element_config
     procedure, public :: reset => reset_species_element_config
     procedure, public :: broadcast => broadcast_species_element_config
     procedure, public, nopass :: get_default_name => get_default_name_species_element_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_species_element_config
  end type species_element_config_type
  
  type(species_element_config_type), dimension(:), allocatable :: species_element_config

contains
  !> FIXME : Add documentation   
  subroutine check_species(report_unit,beta,tite,alne,dbetadrho_spec)
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: beta, tite
    real, intent(out) :: alne, dbetadrho_spec
    integer :: is
    real :: aln, alp, charge, ee, ne, ptot, zeff_calc
    write (report_unit, fmt="('Number of species: ',i3)") nspec
    zeff_calc = 0.
    charge = 0.
    aln = 0.
    alne = 0.
    alp = 0.
    ptot = 0.
    do is=1, nspec
       write (report_unit, *) 
       write (report_unit, fmt="('  Species ',i3)") is
       if (spec(is)%type == 1) write (report_unit, fmt="('    Type:             Ion')")
       if (spec(is)%type == 2) write (report_unit, fmt="('    Type:             Electron')")
       if (spec(is)%type == 3) write (report_unit, fmt="('    Type:             Slowing-down')")
       write (report_unit, fmt="('    Charge:         ',f7.3)") spec(is)%z
       write (report_unit, fmt="('    Mass:             ',es11.4)") spec(is)%mass
       write (report_unit, fmt="('    Density:        ',f7.3)") spec(is)%dens
       write (report_unit, fmt="('    Temperature:    ',f7.3)") spec(is)%temp
       write (report_unit, fmt="('    Collisionality:   ',es11.4)") spec(is)%vnewk
       write (report_unit, fmt="('    Normalized Inverse Gradient Scale Lengths:')")
       write (report_unit, fmt="('      Temperature:  ',f7.3)") spec(is)%tprim
       write (report_unit, fmt="('      Density:      ',f7.3)") spec(is)%fprim
       write (report_unit, fmt="('      Parallel v:   ',f7.3)") spec(is)%uprim
       if (spec(is)%bess_fac.ne.1.0) &
            write (report_unit, fmt="('      Bessel function arg multiplied by:   ',f7.3)") spec(is)%bess_fac
       !        write (report_unit, fmt="('    Ignore this:')")
       !        write (report_unit, fmt="('    D_0: ',es10.4)") spec(is)%dens0
       if (spec(is)%type /= 2) then
          zeff_calc = zeff_calc + spec(is)%dens*spec(is)%z**2
          charge = charge + spec(is)%dens*spec(is)%z
          aln = aln + spec(is)%dens*spec(is)%z*spec(is)%fprim
       else
          alne = alne + spec(is)%dens*spec(is)%z*spec(is)%fprim
          ne = spec(is)%dens
          ee = spec(is)%z
       end if
       alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
       ptot = ptot + spec(is)%dens * spec(is)%temp
    end do
    
    if (.not. has_electron_species(spec)) then
       ptot = ptot + 1./tite   ! electron contribution to pressure
       alp = alp + aln/tite    ! assuming charge neutrality, electron contribution to alp
    end if
    
    alp = alp / ptot
    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    
    write (report_unit, fmt="('Calculated Z_eff: ',f7.3)") zeff_calc
    
    if (has_electron_species(spec)) then
       if (abs(charge+ne*ee) > 1.e-2) then
          if (charge+ne*ee < 0.) then
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('You are neglecting an ion species.')")
             write (report_unit, fmt="('This species has a charge fraction of ',f7.3)") abs(charge+ne*ee)
             write (report_unit, &
                  & fmt="('and a normalized inverse density gradient scale length of ',f7.3)") &
                  (aln+alne)/(charge+ne*ee)
             write (report_unit, fmt="('################# WARNING #######################')")
          else
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('There is an excess ion charge fraction of ',f7.3)") abs(charge+ne*ee)
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
             write (report_unit, fmt="('################# WARNING #######################')")
          end if
       else
          if (abs(aln+alne) > 1.e-2) then
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('The density gradients are inconsistent'/' a/lni =',e12.4,' but alne =',e12.4)") aln, alne
             write (report_unit, fmt="('################# WARNING #######################')")
          end if
       end if
    else
       if (charge > 1.01) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('There is an excess ion charge fraction of ',f7.3)") charge-1.
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
       end if
       if (charge < 0.99) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You are neglecting an ion species.')")
          write (report_unit, fmt="('This species has a charge fraction of ',f7.3)") abs(charge-1.)
          write (report_unit, fmt="('################# WARNING #######################')")
       end if
    end if
    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    
    write (report_unit, *) 
    write (report_unit, fmt="('GS2 beta parameter = ',f9.4)") beta
    write (report_unit, fmt="('Total beta = ',f9.4)") beta*ptot
    write (report_unit, *) 
    write (report_unit, fmt="('The total normalized inverse pressure gradient scale length is ',f10.4)") alp
    dbetadrho_spec = -beta*ptot*alp
    write (report_unit, fmt="('corresponding to d beta / d rho = ',f10.4)") dbetadrho_spec
  end subroutine check_species

  !> FIXME : Add documentation   
  subroutine wnml_species(unit)
    implicit none
    integer, intent(in) :: unit
    integer :: i
    character (100) :: line
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "species_knobs"
    write (unit, fmt="(' nspec = ',i2)") nspec
    write (unit, fmt="(' /')")

    do i=1,nspec
       write (unit, *)
       write (line, *) i
       write (unit, fmt="(' &',a)") &
            & trim("species_parameters_"//trim(adjustl(line)))
       write (unit, fmt="(' z = ',e13.6)") spec(i)%z
       write (unit, fmt="(' mass = ',e13.6)") spec(i)%mass
       write (unit, fmt="(' dens = ',e13.6)") spec(i)%dens
       write (unit, fmt="(' temp = ',e13.6)") spec(i)%temp
       write (unit, fmt="(' tprim = ',e13.6)") spec(i)%tprim
       write (unit, fmt="(' fprim = ',e13.6)") spec(i)%fprim
       write (unit, fmt="(' uprim = ',e13.6)") spec(i)%uprim
       write (unit, fmt="(' vcrit = ',e13.6)") spec(i)%vcrit
       write (unit, fmt="(' vcprim = ',e13.6)") spec(i)%vcprim
       if (spec(i)%uprim2 /= 0.) write (unit, fmt="(' uprim2 = ',e13.6)") spec(i)%uprim2
       write (unit, fmt="(' vnewk = ',e13.6)") spec(i)%vnewk
       if (spec(i)%type == ion_species) &
            write (unit, fmt="(a)") ' type = "ion" /'
       if (spec(i)%type == electron_species) &
            write (unit, fmt="(a)") ' type = "electron"  /'
       if (spec(i)%type == tracer_species) &
            write (unit, fmt="(a)") ' type = "trace"  /'
       if (spec(i)%f0type == f0_maxwellian) &
            write (unit, fmt="(a)") ' f0type = "maxwellian" /'
       if (spec(i)%f0type == f0_sdanalytic) &
            write (unit, fmt="(a)") ' f0type = "sdanalytic"  /'
       if (spec(i)%f0type == f0_tabulated) &
            write (unit, fmt="(a)") ' f0type = "tabulated"  /'
       write (unit, fmt="(' dens0 = ',e13.6)") spec(i)%dens0
       if(spec(i)%bess_fac.ne.1.0) &
            write (unit, fmt="(' bess_fac = ',e13.6)") spec(i)%bess_fac
    end do
  end subroutine wnml_species

  !> FIXME : Add documentation   
  subroutine init_species(species_config_in, species_elements_config_in)
    implicit none
    type(species_config_type), intent(in), optional :: species_config_in
    type(species_element_config_type), intent(in), dimension(:), optional :: species_elements_config_in

    if (initialized) return
    initialized = .true.

    call read_parameters(species_config_in, species_elements_config_in)

  end subroutine init_species

  !> FIXME : Add documentation   
  subroutine read_parameters(species_config_in, species_elements_config_in)
    use file_utils, only: input_unit, error_unit, get_indexed_namelist_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast, mp_abort
    implicit none
    type(species_config_type), intent(in), optional :: species_config_in
    type(species_element_config_type), intent(in), dimension(:), optional :: species_elements_config_in
    real :: z, mass, dens, dens0, u0, temp, tprim, fprim, uprim, uprim2, vnewk, nustar, nu, nu_h
    real :: tperp0, tpar0, bess_fac, vcrit, vcprim
    real, dimension(:), allocatable :: tmp_bcast
    character(20) :: type, f0type
    integer :: is
    integer :: ierr

    type (text_option), dimension (5), parameter :: typeopts = &
         (/ text_option('default', ion_species), &
            text_option('ion', ion_species), &
            text_option('electron', electron_species), &
            text_option('e', electron_species), &
            text_option('trace', tracer_species) /)

    type (text_option), dimension (3), parameter :: f0_opts = &
         (/ text_option('maxwellian', f0_maxwellian), &
            text_option('tabulated', f0_tabulated), &
            text_option('sdanalytic', f0_sdanalytic) /)

    if (present(species_config_in)) species_config = species_config_in

    call species_config%init(name = 'species_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    me = species_config%me
    nspec = species_config%nspec
    zi_fac = species_config%zi_fac

    exist = species_config%exist
    
    if (proc0) then
       if (nspec < 1) then
          ierr = error_unit()
          write (unit=ierr, &
               fmt="('Invalid nspec in species_knobs: ', i5)") nspec
          call mp_abort('Invalid nspec in species_knobs')
       end if
    end if

    allocate (spec(nspec))

    if (present(species_elements_config_in)) species_element_config = species_elements_config_in

    if (.not.allocated(species_element_config)) allocate(species_element_config(nspec))

    if (size(species_element_config) .ne. nspec) then
       if (proc0) print*,"inconsistent number of config elements"
    endif
       
    do is = 1, nspec
       call species_element_config(is)%init(name = 'species_parameters', requires_index = .true., index = is)

       ! Copy out internal values into module level parameters
       bess_fac = species_element_config(is)%bess_fac
       dens = species_element_config(is)%dens
       dens0 = species_element_config(is)%dens0
       f0type = species_element_config(is)%f0type
       fprim = species_element_config(is)%fprim
       mass = species_element_config(is)%mass
       nu = species_element_config(is)%nu
       nu_h = species_element_config(is)%nu_h
       nustar = species_element_config(is)%nustar
       temp = species_element_config(is)%temp
       tpar0 = species_element_config(is)%tpar0
       tperp0 = species_element_config(is)%tperp0
       tprim = species_element_config(is)%tprim
       type = species_element_config(is)%type
       u0 = species_element_config(is)%u0
       uprim = species_element_config(is)%uprim
       uprim2 = species_element_config(is)%uprim2
       vcprim = species_element_config(is)%vcprim
       vcrit = species_element_config(is)%vcrit
       vnewk = species_element_config(is)%vnewk
       z = species_element_config(is)%z

       spec(is)%z = z
       spec(is)%mass = mass
       spec(is)%dens = dens
       spec(is)%dens0 = dens0
       spec(is)%u0 = u0
       spec(is)%tperp0 = tperp0
       spec(is)%tpar0 = tpar0
       spec(is)%temp = temp
       spec(is)%tprim = tprim
       spec(is)%fprim = fprim
       spec(is)%uprim = uprim
       spec(is)%uprim2 = uprim2
       spec(is)%vnewk = vnewk
       spec(is)%nustar = nustar
       spec(is)%nu = nu
       spec(is)%nu_h = nu_h

       spec(is)%stm = sqrt(temp/mass)
       spec(is)%zstm = z/sqrt(temp*mass)
       spec(is)%tz = temp/z
       spec(is)%zt = z/temp
       spec(is)%smz = abs(sqrt(temp*mass)/z)

       spec(is)%bess_fac=bess_fac
       spec(is)%vcrit=vcrit
       spec(is)%vcprim=vcprim

       ierr = error_unit()
       call get_option_value (type, typeopts, spec(is)%type, ierr, "type in species_parameters_x",.true.)
       call get_option_value (f0type, f0_opts, spec(is)%f0type, ierr, "f0type in species_parameters_x",.true.)
    end do

    !We use a temporary array here to allow us to broadcast all species data at once
    !May be possible to avoid this by simply calling broadcast(spec%z) etc.
    allocate(tmp_bcast(nspec))
    tmp_bcast=spec%z ; call broadcast(tmp_bcast) ; spec%z=tmp_bcast
    tmp_bcast=spec%mass ; call broadcast(tmp_bcast) ; spec%mass=tmp_bcast
    tmp_bcast=spec%dens ; call broadcast(tmp_bcast) ; spec%dens=tmp_bcast
    tmp_bcast=spec%dens0 ; call broadcast(tmp_bcast) ; spec%dens0=tmp_bcast
    tmp_bcast=spec%u0 ; call broadcast(tmp_bcast) ; spec%u0=tmp_bcast
    tmp_bcast=spec%tperp0 ; call broadcast(tmp_bcast) ; spec%tperp0=tmp_bcast
    tmp_bcast=spec%tpar0 ; call broadcast(tmp_bcast) ; spec%tpar0=tmp_bcast
    tmp_bcast=spec%temp ; call broadcast(tmp_bcast) ; spec%temp=tmp_bcast
    tmp_bcast=spec%tprim ; call broadcast(tmp_bcast) ; spec%tprim=tmp_bcast
    tmp_bcast=spec%fprim ; call broadcast(tmp_bcast) ; spec%fprim=tmp_bcast
    tmp_bcast=spec%uprim ; call broadcast(tmp_bcast) ; spec%uprim=tmp_bcast
    tmp_bcast=spec%uprim2 ; call broadcast(tmp_bcast) ; spec%uprim2=tmp_bcast
    tmp_bcast=spec%vnewk ; call broadcast(tmp_bcast) ; spec%vnewk=tmp_bcast
    tmp_bcast=spec%nu ; call broadcast(tmp_bcast) ; spec%nu=tmp_bcast
    tmp_bcast=spec%nu_h ; call broadcast(tmp_bcast) ; spec%nu_h=tmp_bcast
    tmp_bcast=spec%nustar ; call broadcast(tmp_bcast) ; spec%nustar=tmp_bcast
    tmp_bcast=spec%stm ; call broadcast(tmp_bcast) ; spec%stm=tmp_bcast
    tmp_bcast=spec%zstm ; call broadcast(tmp_bcast) ; spec%zstm=tmp_bcast
    tmp_bcast=spec%tz ; call broadcast(tmp_bcast) ; spec%tz=tmp_bcast
    tmp_bcast=spec%zt ; call broadcast(tmp_bcast) ; spec%zt=tmp_bcast
    tmp_bcast=spec%smz ; call broadcast(tmp_bcast) ; spec%smz=tmp_bcast
    tmp_bcast=spec%type ; call broadcast(tmp_bcast) ; spec%type=int(tmp_bcast)
    tmp_bcast=spec%f0type ; call broadcast(tmp_bcast) ; spec%f0type=int(tmp_bcast)
    tmp_bcast=spec%bess_fac ; call broadcast(tmp_bcast) ; spec%bess_fac=tmp_bcast
    tmp_bcast=spec%vcrit ; call broadcast(tmp_bcast) ; spec%vcrit=tmp_bcast
    tmp_bcast=spec%vcprim ; call broadcast(tmp_bcast) ; spec%vcprim=tmp_bcast
    deallocate(tmp_bcast)

  end subroutine read_parameters

  !> FIXME : Add documentation   
  pure function has_ion_species (spec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    logical :: has_ion_species
    has_ion_species = any(spec%type == ion_species)
  end function has_ion_species

  !> FIXME : Add documentation   
  subroutine calculate_slowingdown_parameters(alpha_index,vc,vcprim)
    use constants, only: pi
    use mp, only: proc0
    implicit none
    integer, intent(in):: alpha_index
    real, intent(inout):: vc, vcprim
    integer:: electron_spec, main_ion_species, is,js
    real:: vinj, Te_prim, ZI_fac_local, ne_prim, ni_prim, ne, vte
    
    electron_spec = -1
    main_ion_species = -1

    do is = 1,nspec
       if (spec(is)%type .eq. electron_species)                            electron_spec = is
       if (main_ion_species < 1 .and. spec(is)%type .eq. ion_species)      main_ion_species = is
    end do

    if (electron_spec .GT. 0) then
       me = spec(electron_spec)%mass
       ne = spec(electron_spec)%dens
       ne_prim = spec(electron_spec)%fprim
       Te_prim = spec(electron_spec)%tprim
       vte = sqrt(spec(electron_spec)%temp/spec(electron_spec)%mass)
    else
       if (me .LT. 0.0) then
          if (proc0) write(*,*) "WARNING: vte not specified in species_knobs. Assuming reference species is deuterium, so me = 0.0002723085"
          me = 0.0002723085
       end if             

       ne = sum(spec(:)%z*spec(:)%dens)
       ne_prim = sum(spec(:)%z*spec(:)%dens*spec(:)%fprim)/ne
       Te_prim = spec(main_ion_species)%tprim
       vte = sqrt(spec(main_ion_species)%temp/spec(main_ion_species)%mass)
       if (proc0) then
          write(*,*) "Since electrons are adiabatic, we need to improvise electron parameters to "
          write(*,*) "caluclate alpha species properties. Imposed by global quasineutrality:"
          write(*,*) "  ne = ", ne
          write(*,*) "  vte = ", vte
          write(*,*) "  ne_prim = ", ne_prim
          write(*,*) "  Te_prim = ", Te_prim
       end if
    end if

    ni_prim = 0.0
    ! If ZI_fac is not specified, calculate it from existing ion species:
    if (ZI_fac .LT. 0.0) then
       ZI_fac_local = 0.0
       do js = 1,nspec
          if (spec(js)%type .eq. ion_species) then 
             ZI_fac_local = ZI_fac_local   + spec(alpha_index)%mass*spec(js)%dens*spec(js)%z**2/(spec(js)%mass*ne)
             ni_prim = ni_prim + spec(alpha_index)%mass*spec(js)%dens*spec(js)%z**2*spec(js)%fprim/(spec(js)%mass*ne)
          end if
       end do
    else
       ZI_fac_local = ZI_fac
    end if
    ni_prim = ni_prim / ZI_fac_local

    ! Calculate vc/vinj
    vinj = sqrt(spec(alpha_index)%temp/spec(alpha_index)%mass)
    vc = (0.25*3.0*sqrt(pi)*ZI_fac_local*me/spec(alpha_index)%mass)**(1.0/3.0)*vte/vinj
    vcprim = ni_prim - ne_prim + 1.5*Te_prim

  end subroutine calculate_slowingdown_parameters

  !> FIXME : Add documentation   
  subroutine calculate_f0_arrays(egrid)
     implicit none
     real, dimension (:,:), intent (in) :: egrid
     integer:: negrid,is,ie
     real:: v, f0, df0dE, f0prim
    
     negrid = size(egrid(:,1))

     allocate(nonmaxw_corr(negrid,nspec))
     allocate(f0_values(negrid,nspec))
     allocate(dlnf0drho(negrid,nspec))

     do is = 1,nspec
        is_global = is
        do ie = 1,negrid
           v = sqrt(egrid(ie,is))

           select case (spec(is)%f0type)
           case (f0_maxwellian)
              f0_values(ie,is) = eval_f0(v)
              nonmaxw_corr(ie,is) = 1.0
              dlnf0drho(ie,is) = -spec(is)%fprim - (egrid(ie,is) - 1.5)*spec(is)%tprim
           case (f0_sdanalytic)
              call eval_f0_sdanalytic(is,v,f0,df0dE,f0prim)
              f0_values(ie,is) = f0
              nonmaxw_corr(ie,is) = -spec(is)%temp*df0dE
              dlnf0drho(ie,is) = f0prim
           end select
        end do
        if (spec(is)%f0type .EQ. f0_tabulated) then
           call calculate_f0_arrays_tabulated(is,egrid)
        end if
        
     end do

  end subroutine calculate_f0_arrays

  !> FIXME : Add documentation   
  subroutine set_current_f0_species(is)
     implicit none
     integer, intent(in):: is
     is_global = is
  end subroutine set_current_f0_species


  !> FIXME : Add documentation
  !!
  !! @note Relies upon is_global to determine species. This function only takes one argument so that
  !! it can be passed as an argument to genquad
  function eval_f0(v)
     use constants, only: pi
     use mp, only: mp_abort
     implicit none
     real, intent(in):: v
     real:: eval_f0
     real::f0temp, dummy1,dummy2

     if (is_global .LT. 0)  call mp_abort("is_global not specified when calling eval_f0")

     select case(spec(is_global)%f0type)
     case (f0_maxwellian)
        f0temp = exp(-v**2)/pi**1.5        
     case (f0_sdanalytic)
        call eval_f0_sdanalytic(is_global,v,f0temp,dummy1,dummy2)
     case (f0_tabulated)
        call mp_abort("ERROR: eval_f0 cannot yet handle tabulated option.")
     end select

     eval_f0 = f0temp

  end function eval_f0

  !> FIXME : Add documentation   
  subroutine eval_f0_sdanalytic(is,v,f0,df0dE_out,f0prim)
    use constants, only: pi
    use mp, only: mp_abort
    use spfunc, only: erf => erf_ext
    implicit none
    integer,intent(in):: is
    real,intent(in):: v
    real,intent(inout):: f0,df0dE_out, f0prim
    real:: A, df0dv,vcprim, vcrit

    if (spec(is)%vcrit .LE. 0.0) then
       call calculate_slowingdown_parameters(is,vcrit,vcprim)
       if (abs(spec(is)%vcprim + 999.0) .LT. 1.0e-3)  vcprim = spec(is)%vcprim
    else
       !> If vcrit is specified, there is no way of knowing its radial derivative, so set dvc/drho = 0
       vcrit = spec(is)%vcrit
       if (abs(spec(is)%vcprim + 999.0) .LT. 1.0e-3)   vcprim = 0.0
    end if


    A = (3.0/(4.0*pi))/log( (vcrit**3 + 1.0)/(vcrit**3))

    f0 = A/(vcrit**3 + v**3)
    df0dv = -A*3.0*v**2/(vcrit**3 + v**3)**2
    df0dE_out =  (0.5/v)*df0dv / ( f0 * spec(is)%temp )
    f0prim = -spec(is)%fprim - (-(vcrit**3/(vcrit**3 + v**3)) + (1.0/( log(1.0 + vcrit**(-3)) * (1.0 + vcrit**(3)))))*vcprim

    if (v .GT. 1.0) then
       f0 = 0.0
       f0prim = 0.0
       df0dE_out = 1.0
    end if

!    write(*,*) v, vcrit,vcprim,f0prim

  end subroutine eval_f0_sdanalytic

  !> This subroutine calculates f0 on the grid from an external
  !! input file. The grid on the input file can differ from that 
  !! of gs2. A cubic spline is used to interpolate between the two.
  !! The user can either specify df0/dE or it can be estimated internally.
  !! The radial derivative of f0 as a function of energy should also be given.
  subroutine calculate_f0_arrays_tabulated(is,egrid)
    !
    ! The egrid is already given by get_legendre_grids. 
    !
    ! The input file can take two forms, as controlled by the
    ! control parameter mode, the first integer read from the file
    !  mode = 1: E points in first columns, F0(E) in second. To calculate
    !            generalized temperature, a cubic spline is used to 
    !            interpolate and the slope is taken from that.
    !       = 2: First column is E, second is F0(E), the third is (1/F0)*dF0/dE
    !       = 3: E, F0(E), -(1/F0)*dF0/drho 
    !       = 4: E, F0(E),dF0/dE, -(1/F0)*dF0/drho 
    use mp, only: broadcast
    use constants, only: pi
    use file_utils, only: run_name
    use splines, only: fitp_curvd, fitp_curv1, fitp_curv2
    implicit none
    integer, intent(in) :: is
    real, dimension(:,:),intent(in):: egrid
    integer:: f0in_unit = 21, mode, ie, ierr, numdat, negrid
    logical:: acceptable_spline
    real:: moment0, sigma, sigma_fac, df0dE
    real, dimension(:), allocatable:: f0_values_dat, df0dE_dat, egrid_dat, &
                                      yp, temp, df0drho_dat
    
    ! Open file and read column option
    open(unit=f0in_unit,file=trim(run_name)//'.f0in',status='old',action='read')
    read(f0in_unit,*) mode
    read(f0in_unit,*) numdat

    allocate(f0_values_dat(numdat))
    allocate(df0dE_dat(numdat))
    allocate(df0drho_dat(numdat))
    allocate(egrid_dat(numdat))
    allocate(yp(numdat))
    allocate(temp(numdat))

    negrid = size(egrid(:,1))
    
    sigma = 1.0
    sigma_fac = 2.0
    acceptable_spline = .false.

    if (mode .EQ. 1) then
       ! Read f0 values
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie)
       end do
       close(f0in_unit)

       ! Perform cubic spline to get F0 and its slope at grid points
       ! If F0 ends up being negative or zero anywhere on the grid, double the tension
       ! factor and try again

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat,yp,sigma)

             ! Calculate d/dE lnF0 to get generalised temperature
             df0dE =  fitp_curvd(egrid(ie,is),numdat,egrid_dat, f0_values_dat,yp,sigma)  / &
                                (f0_values(ie,is) * spec(is)%temp)
             nonmaxw_corr(ie,is) = -spec(is)%temp*df0dE

          end do

          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if
       end do

       dlnf0drho(:,is) = - spec(is)%fprim

    else if (mode .EQ. 2) then
       ! Read both f0 and df0/dE
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid

             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                f0_values_dat,yp,sigma)
          end do

          ! Generate spline parameters for df0/dE
          call fitp_curv1(numdat,egrid_dat,df0dE_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             df0dE     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0dE_dat,yp,sigma)  / spec(is)%temp
             nonmaxw_corr(ie,is) = -spec(is)%temp*df0dE
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

       dlnf0drho(:,is) = - spec(is)%fprim

    else if (mode .EQ. 3) then
       ! Read f0 and df0/drho
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0drho_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                 f0_values_dat,yp,sigma)
             ! Calculate d/dE lnF0 to get generalised temperature
             df0dE = fitp_curvd(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat,yp,sigma) / (f0_values(ie,is) * spec(is)%temp)
             nonmaxw_corr(ie,is) = -spec(is)%temp*df0dE
          end do

          ! Generate spline parameters for df0/drho
          call fitp_curv1(numdat,egrid_dat,df0drho_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0prim at grid points
             dlnf0drho(ie,is)    = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0drho_dat,yp,sigma)
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

    else if (mode .EQ. 4) then
       ! Read f0, df0/dE, df0/drho
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie), df0drho_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid

             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                f0_values_dat,yp,sigma)
          end do

          ! Generate spline parameters for df0/dE
          call fitp_curv1(numdat,egrid_dat,df0dE_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             df0dE     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0dE_dat,yp,sigma)  / spec(is)%temp
             nonmaxw_corr(ie,is) = -spec(is)%temp*df0dE
          end do
 
          ! Generate spline parameters for df0/drho
          call fitp_curv1(numdat,egrid_dat,df0drho_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             dlnf0drho(ie,is)     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                 df0drho_dat,yp,sigma) 
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

       dlnf0drho(:,is) = - spec(is)%fprim

    else
       write(*,*) "ERROR. First line in f0 input file should be mode" 
       stop 1
    end if

    ! Now calculate moments of F0 from trapezoidal rule.
    ! Should probably do this by quadrature, 
    ! but le_grids is not defined yet

    moment0 = 0.0
    do ie = 1,numdat-1
       moment0  = moment0 + 0.5*(egrid_dat(ie+1)-egrid_dat(ie)) * &
                           ( sqrt(egrid_dat(ie))*f0_values_dat(ie) + &
                             sqrt(egrid_dat(ie+1))*f0_values_dat(ie+1) )
    end do
    moment0 = moment0*4.0*pi

    
    deallocate(egrid_dat,f0_values_dat,df0dE_dat,df0drho_dat,yp,temp)

  end subroutine calculate_f0_arrays_tabulated

  !> FIXME : Add documentation   
  pure function has_electron_species (spec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    logical :: has_electron_species
    has_electron_species = any(spec%type == electron_species)
  end function has_electron_species

  !> FIXME : Add documentation   
  pure function has_nonmaxw_species (spec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    logical :: has_nonmaxw_species
    has_nonmaxw_species = .NOT. all(spec%f0type == f0_maxwellian)
  end function has_nonmaxw_species

  !> FIXME : Add documentation   
  subroutine finish_species
    implicit none
    if(allocated(spec)) deallocate (spec)
    !    if(allocated(f0_values)) deallocate (f0_values,nonmaxw_corr,dlnf0drho)
    call species_config%reset()
    if (allocated(species_element_config)) deallocate(species_element_config)
    initialized = .false.
  end subroutine finish_species

  !> This routine is used only for Trinity
  subroutine determine_species_order
    use mp, only: broadcast, proc0, mp_abort

    implicit none

    integer is
     if (nspec == 1) then
        ions = 1
        electrons = 0
        impurity = 0
     else
        ! if 2 or more species in GS2 calculation, figure out which is main ion
        ! and which is electron via mass (main ion mass assumed to be one)
        do is = 1, nspec
           if (abs(spec(is)%mass-1.0) <= epsilon(0.0)) then
              ions = is
           else if (spec(is)%mass < 0.3) then
              ! for electrons, assuming electrons are at least a factor of 3 less massive
              ! than main ion and other ions are no less than 30% the mass of the main ion
              electrons = is
           else if (spec(is)%mass > 1.0 + epsilon(0.0)) then
              impurity = is
           else
              if (proc0) write (*,*) &
                   "Error: TRINITY requires the main ions to have mass 1", &
                   "and the secondary ions to be impurities (mass > 1)"
              call mp_abort("Error: TRINITY requires the main ions to have mass 1 and the secondary ions to be impurities (mass > 1)")
           end if
        end do
     end if
  end subroutine determine_species_order

  !> FIXME : Add documentation   
  subroutine set_overrides(prof_ov)
    use overrides, only: profiles_overrides_type
    use unit_tests, only: debug_message
    type(profiles_overrides_type), intent(in) :: prof_ov
    integer :: is
    integer, parameter :: verbosity=3
    character(len=800) :: message
    call debug_message(verbosity, 'species::set_overrides starting')
    do is = 1,nspec
      message = ''
      write(message,*) 'species::set_overrides setting species ', is
      call debug_message(verbosity, trim(message))
      if (prof_ov%override_tprim(is)) spec(is)%tprim = prof_ov%tprim(is)
      if (prof_ov%override_fprim(is)) spec(is)%fprim = prof_ov%fprim(is)
      if (prof_ov%override_temp(is)) spec(is)%temp = prof_ov%temp(is)
      if (prof_ov%override_dens(is)) spec(is)%dens = prof_ov%dens(is)
      if (prof_ov%override_vnewk(is)) spec(is)%vnewk = prof_ov%vnewk(is)
    end do
    call debug_message(verbosity,&
    'species::set_overrides calling calculate_and_broadcast_species_properties')
    call calculate_and_broadcast_species_properties
    call debug_message(verbosity, 'species::set_overrides finishing')
  end subroutine set_overrides

  !> FIXME : Add documentation 
  subroutine calculate_and_broadcast_species_properties
    use mp, only: proc0, broadcast
    use kt_grids, only: gryfx
    integer :: is
    if (proc0) then 
      do is = 1, nspec
         spec(is)%stm = sqrt(spec(is)%temp/spec(is)%mass)
         spec(is)%zstm = spec(is)%z/sqrt(spec(is)%temp*spec(is)%mass)
         spec(is)%tz = spec(is)%temp/spec(is)%z
         spec(is)%zt = spec(is)%z/spec(is)%temp
         spec(is)%smz = abs(sqrt(spec(is)%temp*spec(is)%mass)/spec(is)%z)
         if (gryfx()) then
            !WARNING THIS NORMALISATION IS A GUESS!!!
            !FEEL FREE TO QUESTION IT AND WORRY ABOUT IT
           spec(is)%vnewk = 0.01 ! spec(is)%vnewk / sqrt(2.0)
         end if
      end do
    end if
    call broadcast (nspec)

    do is = 1, nspec
       call broadcast (spec(is)%dens)
       call broadcast (spec(is)%temp)
       call broadcast (spec(is)%fprim)
       call broadcast (spec(is)%tprim)
       call broadcast (spec(is)%vnewk)
       call broadcast (spec(is)%stm)
       call broadcast (spec(is)%zstm)
       call broadcast (spec(is)%tz)
       call broadcast (spec(is)%zt)
       call broadcast (spec(is)%smz)
    end do
  end subroutine calculate_and_broadcast_species_properties

  !> FIXME : Add documentation   
  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    integer :: is
    do is = 1,nspec
      if (is<10) then
        write(trinpars_unit, "(A20,I1)") '&species_parameters_', is
      else
        write(trinpars_unit, "(A20,I2)") '&species_parameters_', is
      end if
      write (trinpars_unit, *) ' temp = ', spec(is)%temp
      write (trinpars_unit, *) ' dens = ', spec(is)%dens
      write (trinpars_unit, *) ' tprim = ', spec(is)%tprim
      write (trinpars_unit, *) ' fprim = ', spec(is)%fprim
      write (trinpars_unit, *) ' vnewk = ', spec(is)%vnewk
      write (trinpars_unit, "(A1)") '/'
    end do
  end subroutine write_trinity_parameters

  !> Set the module level config types
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_species_config(species_config_in, species_element_config_in)
    use mp, only: mp_abort
    type(species_config_type), intent(in), optional :: species_config_in
    type(species_element_config_type), intent(in), dimension(:), optional :: species_element_config_in

    if (initialized) then
       call mp_abort("Trying to set species_config when already initialized.", to_screen = .true.)
    end if
    if (present(species_config_in)) then
       species_config = species_config_in
    end if
    if (present(species_element_config_in)) then
       species_element_config = species_element_config_in
    end if    
  end subroutine set_species_config   

  !> Get the module level config instance
  function get_species_config()
    type(species_config_type) :: get_species_config
    get_species_config = species_config
  end function get_species_config

  !> Get the module level config instance
  function get_species_element_config()
    type(species_element_config_type), allocatable :: get_species_element_config(:)
    if (allocated(species_element_config)) then
      get_species_element_config = species_element_config
    else
      allocate(get_species_element_config(0))
    end if
  end function get_species_element_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the species_knobs namelist and populates the member variables
  subroutine read_species_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(species_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: nspec
    real :: me, zi_fac

    namelist /species_knobs/ me, nspec, zi_fac

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    me = self%me
    nspec = self%nspec
    zi_fac = self%zi_fac

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = species_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = species_knobs)
    endif

    ! Now copy from local variables into type members
    self%me = me
    self%nspec = nspec
    self%zi_fac = zi_fac

    self%exist = exist
  end subroutine read_species_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_species_config(self, unit)
    implicit none
    class(species_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("me", self%me, unit_internal)
    call self%write_key_val("nspec", self%nspec, unit_internal)
    call self%write_key_val("zi_fac", self%zi_fac, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_species_config

  !> Resets the config object to the initial empty state
  subroutine reset_species_config(self)
    class(species_config_type), intent(in out) :: self
    type(species_config_type) :: empty
    select type (self)
    type is (species_config_type)
       self = empty
    end select
  end subroutine reset_species_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_species_config(self)
    use mp, only: broadcast
    implicit none
    class(species_config_type), intent(in out) :: self
    call broadcast(self%me)
    call broadcast(self%nspec)
    call broadcast(self%zi_fac)

    call broadcast(self%exist)
  end subroutine broadcast_species_config
    
  !> Gets the default name for this namelist
  function get_default_name_species_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_species_config
    get_default_name_species_config = "species_knobs"
  end function get_default_name_species_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_species_config()
    implicit none
    logical :: get_default_requires_index_species_config
    get_default_requires_index_species_config = .false.
  end function get_default_requires_index_species_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the species_parameters namelist and populates the member variables
  subroutine read_species_element_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(species_element_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    character(len = 20) :: f0type, type
    real :: bess_fac, dens, dens0, fprim, mass, nu, nu_h, nustar, temp, tpar0, tperp0, tprim, u0, uprim, uprim2, vcprim, vcrit, vnewk, z

    namelist /species_parameters/ bess_fac, dens, dens0, f0type, fprim, mass, nu, nu_h, nustar, temp, tpar0, tperp0, tprim, type, u0, uprim, uprim2, vcprim, vcrit, vnewk, z

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    bess_fac = self%bess_fac
    dens = self%dens
    dens0 = self%dens0
    f0type = self%f0type
    fprim = self%fprim
    mass = self%mass
    nu = self%nu
    nu_h = self%nu_h
    nustar = self%nustar
    temp = self%temp
    tpar0 = self%tpar0
    tperp0 = self%tperp0
    tprim = self%tprim
    type = self%type
    u0 = self%u0
    uprim = self%uprim
    uprim2 = self%uprim2
    vcprim = self%vcprim
    vcrit = self%vcrit
    vnewk = self%vnewk
    z = self%z

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = species_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = species_parameters)
    endif

    ! Now copy from local variables into type members
    self%bess_fac = bess_fac
    self%dens = dens
    self%dens0 = dens0
    self%f0type = f0type
    self%fprim = fprim
    self%mass = mass
    self%nu = nu
    self%nu_h = nu_h
    self%nustar = nustar
    self%temp = temp
    self%tpar0 = tpar0
    self%tperp0 = tperp0
    self%tprim = tprim
    self%type = type
    self%u0 = u0
    self%uprim = uprim
    self%uprim2 = uprim2
    self%vcprim = vcprim
    self%vcrit = vcrit
    self%vnewk = vnewk
    self%z = z

    self%exist = exist
  end subroutine read_species_element_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_species_element_config(self, unit)
    implicit none
    class(species_element_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("bess_fac", self%bess_fac, unit_internal)
    call self%write_key_val("dens", self%dens, unit_internal)
    call self%write_key_val("dens0", self%dens0, unit_internal)
    call self%write_key_val("f0type", self%f0type, unit_internal)
    call self%write_key_val("fprim", self%fprim, unit_internal)
    call self%write_key_val("mass", self%mass, unit_internal)
    call self%write_key_val("nu", self%nu, unit_internal)
    call self%write_key_val("nu_h", self%nu_h, unit_internal)
    call self%write_key_val("nustar", self%nustar, unit_internal)
    call self%write_key_val("temp", self%temp, unit_internal)
    call self%write_key_val("tpar0", self%tpar0, unit_internal)
    call self%write_key_val("tperp0", self%tperp0, unit_internal)
    call self%write_key_val("tprim", self%tprim, unit_internal)
    call self%write_key_val("type", self%type, unit_internal)
    call self%write_key_val("u0", self%u0, unit_internal)
    call self%write_key_val("uprim", self%uprim, unit_internal)
    call self%write_key_val("uprim2", self%uprim2, unit_internal)
    call self%write_key_val("vcprim", self%vcprim, unit_internal)
    call self%write_key_val("vcrit", self%vcrit, unit_internal)
    call self%write_key_val("vnewk", self%vnewk, unit_internal)
    call self%write_key_val("z", self%z, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_species_element_config

  !> Resets the config object to the initial empty state
  subroutine reset_species_element_config(self)
    class(species_element_config_type), intent(in out) :: self
    type(species_element_config_type) :: empty
    select type (self)
    type is (species_element_config_type)
       self = empty
    end select
  end subroutine reset_species_element_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_species_element_config(self)
    use mp, only: broadcast
    implicit none
    class(species_element_config_type), intent(in out) :: self
    call broadcast(self%bess_fac)
    call broadcast(self%dens)
    call broadcast(self%dens0)
    call broadcast(self%f0type)
    call broadcast(self%fprim)
    call broadcast(self%mass)
    call broadcast(self%nu)
    call broadcast(self%nu_h)
    call broadcast(self%nustar)
    call broadcast(self%temp)
    call broadcast(self%tpar0)
    call broadcast(self%tperp0)
    call broadcast(self%tprim)
    call broadcast(self%type)
    call broadcast(self%u0)
    call broadcast(self%uprim)
    call broadcast(self%uprim2)
    call broadcast(self%vcprim)
    call broadcast(self%vcrit)
    call broadcast(self%vnewk)
    call broadcast(self%z)

    call broadcast(self%exist)
  end subroutine broadcast_species_element_config 

  !> Gets the default name for this namelist
  function get_default_name_species_element_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_species_element_config
    get_default_name_species_element_config = "species_parameters"
  end function get_default_name_species_element_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_species_element_config()
    implicit none
    logical :: get_default_requires_index_species_element_config
    get_default_requires_index_species_element_config = .true.
  end function get_default_requires_index_species_element_config
end module species
