!> Provides random forcing for Alfven wave problem
!> Originally by Hammett, Dorland and Quataert, Dec. 5, 2001
!>
!> init_antenna should be called once per run
!> antenna_amplitudes provides the necessary info per time step
!>
!> added terms needed to calculate heating by antenna (apar_new)
!>
!> two kinds of namelists should be added to the input file:
!>
!> driver:
!>
!>    amplitude : RMS amplitude of | apar_antenna |
!>    w_antenna : frequency of driving, normalized to kpar v_Alfven
!>                Note: the imaginary part is the decorrelation rate
!>    nk_stir   : Number of k modes that should be driven
!>    write_antenna: .true. for direct antenna output to runname.antenna
!>                   default is .false.
!>
!> stir: (stir_1, stir_2, stir_3, etc., one for each k we want to drive)
!>
!>    kx, ky, kz : each of size nk_stir (not larger!)
!>               : list of k components that are to be driven; integers
!>    travel     : logical; choose false to make a standing wave
!>               : default value is .true.
module antenna
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use antenna_data, only: a_ant, b_ant
  implicit none

  private
  public :: init_antenna, dump_ant_amp, finish_antenna, reset_init
  public :: wnml_antenna, check_antenna, no_driver
  public :: antenna_w, antenna_apar, antenna_amplitudes, a_ext_data
  public :: get_jext !GGH

  public :: driver_config_type, stir_config_type, set_antenna_config
  public :: get_driver_config, get_stir_config
  
  complex :: w_antenna
  complex, dimension(:), allocatable :: w_stir
  complex, dimension(:,:,:), allocatable :: apar_new, apar_old
  integer, dimension(:), allocatable :: kx_stir, ky_stir, kz_stir
  real :: amplitude, t0, w_dot
  integer :: nk_stir, out_unit
  logical :: restarting = .false.
  logical :: no_driver = .false.
  logical :: write_antenna = .false.
  logical, dimension (:), allocatable :: trav
  real :: beta_s
  complex :: wtmp
  logical :: initialized = .false.

  !> Used to represent the input configuration of driver. See the
  !> documentation of the [[antenna]] module for more details.
  type, extends(abstract_config_type) :: driver_config_type
     ! namelist : driver
     ! indexed : false
     !> Amplitude of Langevin antenna.
     real :: amplitude = 0.0
     !> Overrides all other options and turns off antenna if `true`.
     logical :: ant_off = .false.
     !> Number of independent Fourier modes driven by antenna.
     integer :: nk_stir = 1
     !> If `true` then try to get initial antenna amplitudes from the
     !> restart file. If not found in restart file then will be set to
     !> 0.
     logical :: restarting = .false.
     !> Sets the time at which the antenna frequency changes from
     !> constant [[driver:w_antenna]] to linearly varying, with linear piece
     !> proprtional to [[driver:w_dot]] multiplied by `time - t0`.
     real :: t0 = -1.0
     !> Frequency of Langevin antenna. Sets the constant part of the
     !> complex antenna frequency.
     complex :: w_antenna = (1.0, 0.0)
     !> Sets the coefficient in front of the time varying antenna
     !> frequency (real) component activated when `time > t0`.
     real :: w_dot = 0.0
     !> Write antenna amplitudes to ASCII file for debugging.
     logical :: write_antenna = .false.
   contains
     procedure, public :: read => read_driver_config
     procedure, public :: write => write_driver_config
     procedure, public :: reset => reset_driver_config
     procedure, public :: broadcast => broadcast_driver_config
     procedure, public, nopass :: get_default_name => get_default_name_driver_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_driver_config
  end type driver_config_type

  type(driver_config_type) :: driver_config

  !> Used to represent the input configuration of stir. See the
  !> documentation of the [[antenna]] module for more details.
  type, extends(abstract_config_type) :: stir_config_type
     ! namelist : stir
     ! indexed : true
     !> Initial amplitude of right-moving component. It is not
     !> necessary to set [[stir:a]] and [[stir:b]] unless you are
     !> doing restarts, which are rather clunky at the moment with the
     !> antenna included.
     real :: a = -1.0
     !> Initial amplitude of left-moving component. It is not
     !> necessary to set [[stir:a]] and [[stir:b]] unless you are
     !> doing restarts, which are rather clunky at the moment with the
     !> antenna included.
     real :: b = -1.0
     !> Mode number for stirring.
     integer :: kx = 1
     !> Mode number for stirring.
     integer :: ky = 1
     !> Mode number for stirring.
     integer :: kz = 1
     !> Launches traveling wave if `true` (default) or standing wave
     !> if `false`.
     logical :: travel = .true.
   contains
     procedure, public :: read => read_stir_config
     procedure, public :: write => write_stir_config
     procedure, public :: reset => reset_stir_config
     procedure, public :: broadcast => broadcast_stir_config
     procedure, public, nopass :: get_default_name => get_default_name_stir_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_stir_config
  end type stir_config_type
  
  type(stir_config_type), dimension(:), allocatable :: stir_config
  
contains
  !> FIXME : Add documentation
  subroutine check_antenna(report_unit)
    use file_utils, only: run_name
    implicit none
    integer, intent(in) :: report_unit
    integer :: i
    if (no_driver) return
    if (amplitude == 0.0) then 
       write (report_unit, *) 
       write (report_unit, fmt="('No Langevin antenna included.')")
       write (report_unit, *) 
    else
       write (report_unit, *) 
       write (report_unit, fmt="('A Langevin antenna is included, with characteristics:')")
       write (report_unit, *) 
       write (report_unit, fmt="('Frequency:  (',e11.4,', ',e11.4,')')") w_antenna
       write (report_unit, fmt="('Number of independent k values: ',i3)") nk_stir
       if (write_antenna) then
          write (report_unit, *) 
          write (report_unit, fmt="('Antenna data will be written to ',a)") trim(run_name)//'.antenna'
          write (report_unit, *) 
       end if
       write (report_unit, fmt="('k values:')")
       do i=1,nk_stir
          if (trav(i)) then
             write (report_unit, fmt="('Travelling wave:')")
             write (report_unit, fmt="('   kx = ',i2,'    ky = ',i2,'    kz = ',i2)") &
                  & kx_stir(i), ky_stir(i), kz_stir(i)
          else
             write (report_unit, fmt="('Standing wave:')")
             write (report_unit, fmt="('   kx = ',i2,'    ky = ',i2,'    kz = ',i2)") &
                  & kx_stir(i), ky_stir(i), kz_stir(i)
          end if
       end do
    end if
  end subroutine check_antenna

  !> FIXME : Add documentation
  subroutine wnml_antenna(unit)
    implicit none
    integer, intent(in) :: unit
    integer :: i
    character (100) :: line
    if (no_driver) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "driver"
    write (unit, fmt="(' ant_off = ',L1)") no_driver
    write (unit, fmt="(' write_antenna = ',L1)") write_antenna
    write (unit, fmt="(' amplitude = ',e17.10)") amplitude
    write (unit, fmt="(' w_antenna = (',e17.10,', ',e17.10,')')") w_antenna
    write (unit, fmt="(' w_dot = ',e17.10)") w_dot
    write (unit, fmt="(' t0 = ',e17.10)") t0
    write (unit, fmt="(' nk_stir = ',i3)") nk_stir
    write (unit, fmt="(' /')")

    do i=1,nk_stir
       write (unit, *)
       write (line, *) i
       write(unit, fmt="(' &',a)") trim("stir_"//trim(adjustl(line)))
       write(unit, fmt="(' kx = ',i2,' ky = ',i2,' kz = ',i2)") &
            kx_stir(i), ky_stir(i), kz_stir(i)
       write(unit, fmt="(' travel = ',L1)") trav(i)
       write(unit, fmt="(' a = (',e20.13,',',e20.13,')')") a_ant(i)
       write(unit, fmt="(' b = (',e20.13,',',e20.13,') /')") b_ant(i)
    end do
  end subroutine wnml_antenna

  !> FIXME : Add documentation
  subroutine init_antenna(driver_config_in, stir_config_in)
    use species, only: spec, nspec, init_species
    use run_parameters, only: beta
    implicit none
    type(driver_config_type), intent(in), optional :: driver_config_in
    type(stir_config_type), intent(in), dimension(:), allocatable, optional :: stir_config_in
    integer :: i
    
    if (initialized) return
    initialized = .true.
    
    if (.not. allocated(w_stir)) then
       call read_parameters(driver_config_in, stir_config_in)
       allocate (w_stir(nk_stir))
    end if
    
    beta_s = 0.
    if (beta > epsilon(0.0)) then
       do i=1,nspec
          !GGH What does this mean?  Is this a normalization?
          !GGH This converts from input frequency (normalized to kpar vA)
          !    to the code frequency
          !GGH This seems to give beta_s = 2*n0*T0/v_A^2
          beta_s = beta_s + beta*spec(i)%dens*spec(i)%mass
       end do
    else
       beta_s = 2.
    end if
  end subroutine init_antenna

  !> FIXME : Add documentation
  subroutine read_parameters(driver_config_in, stir_config_in)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use file_utils, only: open_output_file, input_unit
    use mp, only: proc0
    use antenna_data, only: init_antenna_data
    use gs2_save, only: init_ant_amp
    implicit none
    type(driver_config_type), intent(in), optional :: driver_config_in
    type(stir_config_type), intent(in), dimension(:), allocatable, optional :: stir_config_in
    complex :: a, b
    integer :: kx, ky, kz, i, ierr
    logical :: exist, travel, ant_off

    if (present(driver_config_in)) driver_config = driver_config_in

    call driver_config%init(name = 'driver', requires_index = .false.)

    ! Copy out internal values into module level parameters
    amplitude = driver_config%amplitude
    ant_off = driver_config%ant_off
    nk_stir = driver_config%nk_stir
    restarting = driver_config%restarting
    t0 = driver_config%t0
    w_antenna = driver_config%w_antenna
    w_dot = driver_config%w_dot
    write_antenna = driver_config%write_antenna

    exist = driver_config%exist

    ! If the namelist wasn't present, or if [[driver_config:ant_off]]
    ! was explicitly set to true, then we won't be using the antenna,
    ! so for consistency, we turn it off in this case
    no_driver = driver_config%ant_off .or. .not. driver_config%exist
    driver_config%ant_off = no_driver

    if (no_driver) return

    call init_antenna_data (nk_stir)
    allocate (kx_stir(nk_stir))
    allocate (ky_stir(nk_stir))
    allocate (kz_stir(nk_stir))
    allocate (trav(nk_stir))

    ! BD
    ! get initial antenna amplitudes from restart file
    ! if none are found, a_ant = b_ant = 0 will be returned
    ! and ierr will be non-zero.
    !
    ! TO DO: need to know if there IS a restart file to check...
    !

    ierr = 1
    if (restarting) call init_ant_amp (a_ant, b_ant, nk_stir, ierr)

    if (present(stir_config_in)) stir_config = stir_config_in

    if (.not.allocated(stir_config)) allocate(stir_config(nk_stir))
    if (size(stir_config) .ne. nk_stir) then
      if (proc0) print*,"Warning: inconsistent config size in antenna"
    endif

    do i=1,nk_stir
      call stir_config(i)%init(name = 'stir', requires_index = .true., index = i)

      ! Copy out internal values into module level parameters
      a = stir_config(i)%a
      b = stir_config(i)%b
      kx = stir_config(i)%kx
      ky = stir_config(i)%ky
      kz = stir_config(i)%kz
      travel = stir_config(i)%travel

      kx_stir(i) = kx
      ky_stir(i) = ky
      kz_stir(i) = kz
      trav(i) = travel
      ! If a, b are not specified in the input file:
      if (a == -1.0 .and. b == -1.0) then
        ! And if a, b are not specified in the restart file
        ! (else use values from restart file by default)
        if (ierr /= 0) then
          a_ant(i) = amplitude*cmplx(1.,1.)/2.
          b_ant(i) = amplitude*cmplx(1.,1.)/2.
        end if
        ! Else if a, b ARE specified in the input file (ignore restart file):
      else
        a_ant(i) = a
        b_ant(i) = b
      end if
    end do

    if (proc0) then
       call open_output_file (out_unit, ".antenna")
    end if
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine finish_antenna
    use file_utils, only: close_output_file
    use antenna_data, only: finish_antenna_data
    use mp, only: proc0
    implicit none

    call finish_antenna_data
    if (allocated(w_stir)) deallocate (w_stir)
    if (allocated(kx_stir)) deallocate (kx_stir, ky_stir, kz_stir, trav)
    if (allocated(apar_new)) deallocate (apar_new, apar_old)
    if (proc0) call close_output_file(out_unit)
    initialized = .false.
    call driver_config%reset()
    if (allocated(stir_config)) deallocate(stir_config)
  end subroutine finish_antenna

  !> FIXME : Add documentation
  subroutine antenna_amplitudes (apar)
    use mp, only: broadcast, proc0
    use gs2_time, only: code_dt, code_time
    use kt_grids, only: naky, ntheta0, reality
    use theta_grid, only: theta, ntgrid, gradpar
    use ran, only: ranf
    use constants, only: zi
    
    complex, dimension (-ntgrid:,:,:), intent(out) :: apar
    complex :: force_a, force_b
    real :: dt, sigma, time
    integer :: i, it, ik

    apar=0.

    if (no_driver) return

    if (.not. allocated(apar_new)) then
       allocate (apar_new(-ntgrid:ntgrid,ntheta0,naky)) ; apar_new = 0.
       allocate (apar_old(-ntgrid:ntgrid,ntheta0,naky)) 
    end if

    apar_old = apar_new

    dt = code_dt
    time = code_time

    if (time > t0) then       
       wtmp = w_antenna + w_dot*(time-t0)
    else
       wtmp = w_antenna
    end if

! GGH fixed a bug with the frequency sweeping here.  11.17.05
    do i = 1, nk_stir
       w_stir(i) = gradpar(0)*abs(kz_stir(i))*sqrt(1./beta_s) * wtmp
    end do

    apar = 0.

    do i=1, nk_stir
       
       it = 1 + mod(kx_stir(i)+ntheta0,ntheta0)
       ik = 1 + mod(ky_stir(i)+naky,naky)

       if (proc0) then
          !GGH Decorrelation
          sigma = sqrt(12.*abs(aimag(w_stir(i)))/dt)*amplitude
          
          force_a = cmplx(ranf() - 0.5, ranf() - 0.5) * sigma
          if (trav(i)) then
             force_b = cmplx(ranf() - 0.5, ranf() - 0.5) * sigma
          else
             force_b = force_a
          end if
          
          if (trav(i) .and. abs(aimag(w_stir(i))) < epsilon(0.0)) then
             a_ant(i) = a_ant(i)*exp(-zi*w_stir(i)*dt)+force_a*dt
             b_ant(i) = 0.
          else
             a_ant(i) = a_ant(i)*exp(-zi*w_stir(i)*dt)+force_a*dt
             b_ant(i) = b_ant(i)*exp(zi*conjg(w_stir(i))*dt)+force_b*dt
          end if
       end if

! potential existed for a bug here, if different processors 
! got different random numbers; antennas driving different parts of
! velocity space at the same k could be driven differently.  BD 6.7.05

       call broadcast (a_ant)
       call broadcast (b_ant)

!       if (time < t0) then
!          apar(:,it,ik) = apar(:,it,ik) &
!               + (a_ant(i)+b_ant(i))/sqrt(2.)*exp(zi*kz_stir(i)*theta) &
!               * (0.5-0.5*cos(time*pi/t0)
!       else
! NOTE: Is the apar(:,it,ik) unnecessary on the RHS here (=0)?
! ANSWER: No, we are inside the nk_stir loop.  In general case, apar /= 0 here.
          apar(:,it,ik) = apar(:,it,ik) &
               + (a_ant(i)+b_ant(i))/sqrt(2.)*exp(zi*kz_stir(i)*theta) 
!       end if

       if (write_antenna) then
         if (proc0) write(out_unit, fmt='(8(1x,e13.6),1x,i2)') &
              & time, real((a_ant(i)+b_ant(i))/sqrt(2.)), &
              &     aimag((a_ant(i)+b_ant(i))/sqrt(2.)),real(a_ant(i)), aimag(a_ant(i)), &
              &     real(b_ant(i)), aimag(b_ant(i)), real(wtmp), i
       end if
    end do

    if (reality) then
       apar(:,1,1) = 0.0
       
       do it = 1, ntheta0/2
          apar(:,it+(ntheta0+1)/2,1) = conjg(apar(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    !GGH NOTE: Here apar_new is equal to apar+ apar_ext
    apar_new = apar

  end subroutine antenna_amplitudes

  !> FIXME : Add documentation
  subroutine antenna_apar (kperp2, j_ext)
!      GGH ERR? Is this correct? It uses apar_new rather than just the applied
!      current, so does this include the plasma response? 
!      BD: No error.  Here, apar and apar_new are local variables for the antenna only. 7.1.06
!
!      this routine returns space-centered (kperp**2 * A_ext) at the current time

    use kt_grids, only: naky, ntheta0
    use run_parameters, only: beta, fapar
    use theta_grid, only: ntgrid

    complex, dimension (-ntgrid:,:,:), intent(out) :: j_ext
    real, dimension (-ntgrid:,:,:), intent(in) :: kperp2
    integer :: ik, it, ig

    if (fapar > epsilon(0.0)) then

       if (.not. allocated(apar_new)) return
       
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid-1
                j_ext(ig,it,ik) = 0.5* &
                     ( &
!                  kperp2(ig+1,it,ik)*apar_old(ig+1,it,ik) &
!                 +kperp2(ig,  it,ik)*apar_old(ig,  it,ik) &
                     +kperp2(ig+1,it,ik)*apar_new(ig+1,it,ik) &
                     +kperp2(ig,  it,ik)*apar_new(ig,  it,ik))
             end do
          end do
       end do
     
! GGH Is this some normalization?  Yes.
       if (beta > 0.) j_ext = j_ext * 0.5 / beta

    else ! if apar itself is not being advanced in time, there should be no j_ext
       j_ext = 0.
    end if

  end subroutine antenna_apar

  !> Calculate the external current in the antenna
  subroutine get_jext(j_ext)
    use kt_grids, only: ntheta0, naky,aky, kperp2
    use theta_grid, only: jacob, delthet, ntgrid
    implicit none
    !Passed
    !Intent(in) only needed as initialised to zero at top level
    real, dimension(:,:), intent(in out) ::  j_ext
    !Local
    complex, dimension(:,:,:), allocatable :: j_extz
    integer :: ig,ik, it                             !Indices
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt

    !Get z-centered j_ext at current time
    allocate (j_extz(-ntgrid:ntgrid, ntheta0, naky)) ; j_extz = 0.
    call antenna_apar (kperp2, j_extz)

    !Set weighting factor for z-averages
    !    if (proc0) then
    allocate (wgt(-ntgrid:ntgrid))
    !GGH NOTE: Here wgt is 1/(2*ntgrid)
    wgt = 0.
    do ig=-ntgrid,ntgrid-1
       wgt(ig) = delthet(ig)*jacob(ig)
    end do
    wgt = wgt/sum(wgt)
    !    endif

    !Take average over z
    do ig=-ntgrid, ntgrid-1
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             j_ext(it,ik)=j_ext(it,ik)+real(j_extz(ig,it,ik))*wgt(ig)*fac2
          end do
       end do
    enddo

    !    if (proc0) then
    deallocate (wgt)
    !    endif
    deallocate (j_extz)

  end subroutine get_jext

  !> FIXME : Add documentation
  subroutine a_ext_data (A_ext_old, A_ext_new)
!      this routine returns current and previous A_ext vectors
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:,:,:), intent(out) :: A_ext_old, A_ext_new

    if (.not. allocated(apar_new)) then
       A_ext_old = 0.
       A_ext_new = 0.
       return
    else
       A_ext_old = apar_old
       A_ext_new = apar_new
    end if

  end subroutine a_ext_data

  !> FIXME : Add documentation  
  pure complex function antenna_w ()
    implicit none
    antenna_w = wtmp
  end function antenna_w

  !> Write antenna amplitudes to [[out_unit]]
  !>
  !> @warning This writes if [[driver:write_antenna]] is `.false.`
  !> while [[antenna_amplitudes]] writes if it is `.true.`!
  subroutine dump_ant_amp
    use gs2_time, only: user_time
    use mp, only: proc0
    implicit none
    real :: time
    integer :: i

    if (no_driver) return
    if (.not. write_antenna) then
       time = user_time
       do i=1,nk_stir
         if (proc0) write(out_unit, fmt='(7(1x,e13.6),1x,i2)') &
              & time, real((a_ant(i)+b_ant(i))/sqrt(2.)), &
              &     aimag((a_ant(i)+b_ant(i))/sqrt(2.)),real(a_ant(i)), aimag(a_ant(i)), &
              &     real(b_ant(i)), aimag(b_ant(i)),i
  
       end do
    end if
  end subroutine dump_ant_amp

  !> FIXME : Add documentation
  subroutine reset_init
    ! used when interfacing GS2 with Trinity -- MAB
    implicit none
    initialized = .false.
    call driver_config%reset()
    if(allocated(stir_config)) deallocate(stir_config)
  end subroutine reset_init

  !> Set the module level config types
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_antenna_config(driver_config_in, stir_config_in)
    use mp, only: mp_abort
    type(driver_config_type), intent(in), optional :: driver_config_in
    type(stir_config_type), intent(in), dimension(:), allocatable, optional :: stir_config_in
  
    if (initialized) then
       call mp_abort("Trying to set antenna config when already initialized.", to_screen = .true.)
    end if
    if (present(driver_config_in)) then
       driver_config = driver_config_in
    end if
    if (present(stir_config_in)) then
       stir_config = stir_config_in
    end if    
  end subroutine set_antenna_config

#include "driver_auto_gen.inc"
#include "stir_auto_gen.inc"
end module antenna

