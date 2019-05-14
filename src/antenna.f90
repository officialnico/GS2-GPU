!> Provides random forcing for Alfven wave problem
!! Originally by Hammett, Dorland and Quataert, Dec. 5, 2001
!!
!! init_antenna should be called once per run
!! antenna_amplitudes provides the necessary info per time step
!!
!! added terms needed to calculate heating by antenna (apar_new)
!!
!! two kinds of namelists should be added to the input file:
!! 
!! driver:
!!       
!!    amplitude : RMS amplitude of | apar_antenna |
!!    w_antenna : frequency of driving, normalized to kpar v_Alfven
!!                Note: the imaginary part is the decorrelation rate
!!    nk_stir   : Number of k modes that should be driven
!!    write_antenna: .true. for direct antenna output to runname.antenna
!!                   default is .false.
!!
!! stir: (stir_1, stir_2, stir_3, etc., one for each k we want to drive)
!!
!!    kx, ky, kz : each of size nk_stir (not larger!)
!!               : list of k components that are to be driven; integers
!!    travel     : logical; choose false to make a standing wave
!!               : default value is .true.
module antenna
  use abstract_config, only: abstract_config_type
  use antenna_data, only: nk_stir, a_ant, b_ant
  implicit none

  private
  public :: init_antenna, dump_ant_amp, amplitude, finish_antenna, reset_init
  public :: wnml_antenna, check_antenna, no_driver
  public :: antenna_w, antenna_apar, antenna_amplitudes, a_ext_data

  public :: init_antenna_level_1, init_antenna_level_2
  public :: finish_antenna_level_1, finish_antenna_level_2

  public :: driver_config_type
  public :: stir_config_type
  
  complex :: w_antenna
  complex, dimension(:), allocatable :: w_stir
  complex, dimension(:,:,:), allocatable :: apar_new, apar_old
  integer, dimension(:), allocatable :: kx_stir, ky_stir, kz_stir
  real :: amplitude, t0, w_dot
  integer :: out_unit
  logical :: restarting = .false.
  logical :: no_driver = .false.
  logical :: write_antenna = .false.
  logical, dimension (:), allocatable :: trav
  logical :: ant_off = .false.
  real :: beta_s
  complex :: wtmp
  logical :: initialized = .false.

  !> Used to represent the input configuration of driver
  type, extends(abstract_config_type) :: driver_config_type
     ! namelist : driver
     ! indexed : false
     !> FIXME: Add documentation
     real :: amplitude = 0.0
     !> FIXME: Add documentation
     logical :: ant_off = .false.
     !> FIXME: Add documentation
     integer :: nk_stir = 1
     !> FIXME: Add documentation
     logical :: restarting = .false.
     !> FIXME: Add documentation
     real :: t0 = -1.0
     !> FIXME: Add documentation
     complex :: w_antenna = (1.0, 0.0)
     !> FIXME: Add documentation
     real :: w_dot = 0.0
     !> FIXME: Add documentation
     logical :: write_antenna = .false.
   contains
     procedure, public :: read => read_driver_config
     procedure, public :: write => write_driver_config
     procedure, public :: reset => reset_driver_config
     procedure, public :: broadcast => broadcast_driver_config
  end type driver_config_type

  interface driver_config_type
     procedure :: driver_config_constructor
  end interface driver_config_type

  type(driver_config_type) :: driver_config

  !> Used to represent the input configuration of stir
  type, extends(abstract_config_type) :: stir_config_type
     ! namelist : stir
     ! indexed : true     
     !> FIXME: Add documentation
     real :: a = -1.0
     !> FIXME: Add documentation
     real :: b = -1.0
     !> FIXME: Add documentation
     integer :: kx = 1
     !> FIXME: Add documentation
     integer :: ky = 1
     !> FIXME: Add documentation
     integer :: kz = 1
     !> FIXME: Add documentation
     logical :: travel = .true.
   contains
     procedure, public :: read => read_stir_config
     procedure, public :: write => write_stir_config
     procedure, public :: reset => reset_stir_config
     procedure, public :: broadcast => broadcast_stir_config
  end type stir_config_type

  interface stir_config_type
     procedure :: stir_config_constructor
  end interface stir_config_type
  
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
    write (unit, fmt="(' ant_off = ',L1)") ant_off
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
  subroutine init_antenna_level_1
     call read_parameters
  end subroutine init_antenna_level_1

  !> FIXME : Add documentation
  subroutine finish_antenna_level_1
  end subroutine finish_antenna_level_1

  !> FIXME : Add documentation
  subroutine init_antenna_level_2
    call init_antenna
  end subroutine init_antenna_level_2

  !> FIXME : Add documentation
  subroutine finish_antenna_level_2
    call finish_antenna
  end subroutine finish_antenna_level_2

  !> FIXME : Add documentation
  subroutine init_antenna
    use species, only: spec, nspec, init_species
    use run_parameters, only: beta
    use antenna_data, only: init_antenna_data
    implicit none
    integer :: i
    
    if (initialized) return
    initialized = .true.
    
    if (.not. allocated(w_stir)) then
       !call init_species
       
       !call read_parameters
       if (no_driver) then
          i = -1
          call init_antenna_data (i)
          return
       end if
       
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
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit, get_unused_unit
    use file_utils, only: open_output_file, input_unit
    use mp, only: proc0
    use antenna_data, only: init_antenna_data
    use gs2_save, only: init_ant_amp
    implicit none
    type(driver_config_type), intent(in), optional :: driver_config_in
    type(stir_config_type), intent(in), dimension(:), optional :: stir_config_in
    complex :: a, b
    integer :: kx, ky, kz, i, in_file, ierr
    logical :: exist, travel

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

    no_driver = (.not. exist) .or. ant_off

    if (exist) then
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
    endif
    
    if (no_driver) return

    if (proc0) then
       call get_unused_unit (out_unit)
       call open_output_file (out_unit, ".antenna")
    end if
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine finish_antenna
    use file_utils, only: close_output_file
    use antenna_data, only: finish_antenna_data
    
    implicit none

    call finish_antenna_data
    if (allocated(w_stir)) deallocate (w_stir)
    if (allocated(kx_stir)) deallocate (kx_stir, ky_stir, kz_stir, trav)
    if (allocated(apar_new)) deallocate (apar_new, apar_old)
    call close_output_file(out_unit)
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
  function antenna_w ()
    implicit none
    complex :: antenna_w
    antenna_w = wtmp
  end function antenna_w

  !> FIXME : Add documentation
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

  !---------------------------------------
  ! Following is for the driver config_type
  !---------------------------------------
  
  !> Reads in the driver namelist and populates the member variables
  subroutine read_driver_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(driver_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: ant_off, restarting, write_antenna
    integer :: nk_stir
    complex :: w_antenna
    real :: amplitude, t0, w_dot

    namelist /driver/ amplitude, ant_off, nk_stir, restarting, t0, w_antenna, w_dot, write_antenna

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    amplitude = self%amplitude
    ant_off = self%ant_off
    nk_stir = self%nk_stir
    restarting = self%restarting
    t0 = self%t0
    w_antenna = self%w_antenna
    w_dot = self%w_dot
    write_antenna = self%write_antenna

    ! Now read in the main namelist
    if (self%requiresIndex) then
       call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
       read(in_file, nml = driver, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%name, exist)
       if (exist) read(in_file, nml = driver)
    endif

    ! Now copy from local variables into type members
    self%amplitude = amplitude
    self%ant_off = ant_off
    self%nk_stir = nk_stir
    self%restarting = restarting
    self%t0 = t0
    self%w_antenna = w_antenna
    self%w_dot = w_dot
    self%write_antenna = write_antenna

    self%exist = exist
  end subroutine read_driver_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_driver_config(self, unit)
    implicit none
    class(driver_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("amplitude", self%amplitude, unit_internal)
    call self%write_key_val("ant_off", self%ant_off, unit_internal)
    call self%write_key_val("nk_stir", self%nk_stir, unit_internal)
    call self%write_key_val("restarting", self%restarting, unit_internal)
    call self%write_key_val("t0", self%t0, unit_internal)
    call self%write_key_val("w_antenna", self%w_antenna, unit_internal)
    call self%write_key_val("w_dot", self%w_dot, unit_internal)
    call self%write_key_val("write_antenna", self%write_antenna, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_driver_config

  !> Resets the config object to the initial empty state
  subroutine reset_driver_config(self)
    class(driver_config_type), intent(in out) :: self
    type(driver_config_type) :: empty
    select type (self)
    type is (driver_config_type)
       self = empty
    end select
  end subroutine reset_driver_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_driver_config(self)
    use mp, only: broadcast
    implicit none
    class(driver_config_type), intent(in out) :: self
    call broadcast(self%amplitude)
    call broadcast(self%ant_off)
    call broadcast(self%nk_stir)
    call broadcast(self%restarting)
    call broadcast(self%t0)
    call broadcast(self%w_antenna)
    call broadcast(self%w_dot)
    call broadcast(self%write_antenna)

    call broadcast(self%exist)
  end subroutine broadcast_driver_config

  !> Sets the name and requiresIndex paremeters as appropriate
  function driver_config_constructor()
    implicit none
    type(driver_config_type) :: driver_config_constructor
    driver_config_constructor%name = "driver"
    driver_config_constructor%requiresIndex = .false.
  end function driver_config_constructor
  
  !---------------------------------------
  ! Following is for the stir config_type
  !---------------------------------------
  
  !> Reads in the stir namelist and populates the member variables
  subroutine read_stir_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(stir_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: travel
    real :: a, b
    integer :: kx, ky, kz

    namelist /stir/ a, b, kx, ky, kz, travel

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    a = self%a
    b = self%b
    kx = self%kx
    ky = self%ky
    kz = self%kz
    travel = self%travel

    ! Now read in the main namelist
    if (self%requiresIndex) then
       call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
       read(in_file, nml = stir, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%name, exist)
       if (exist) read(in_file, nml = stir)
    endif

    ! Now copy from local variables into type members
    self%a = a
    self%b = b
    self%kx = kx
    self%ky = ky
    self%kz = kz
    self%travel = travel

    self%exist = exist
  end subroutine read_stir_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_stir_config(self, unit)
    implicit none
    class(stir_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("a", self%a, unit_internal)
    call self%write_key_val("b", self%b, unit_internal)
    call self%write_key_val("kx", self%kx, unit_internal)
    call self%write_key_val("ky", self%ky, unit_internal)
    call self%write_key_val("kz", self%kz, unit_internal)
    call self%write_key_val("travel", self%travel, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_stir_config

  !> Resets the config object to the initial empty state
  subroutine reset_stir_config(self)
    class(stir_config_type), intent(in out) :: self
    type(stir_config_type) :: empty
    select type (self)
    type is (stir_config_type)
       self = empty
    end select
  end subroutine reset_stir_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_stir_config(self)
    use mp, only: broadcast
    implicit none
    class(stir_config_type), intent(in out) :: self
    call broadcast(self%a)
    call broadcast(self%b)
    call broadcast(self%kx)
    call broadcast(self%ky)
    call broadcast(self%kz)
    call broadcast(self%travel)

    call broadcast(self%exist)
  end subroutine broadcast_stir_config

  !> Sets the name and requiresIndex paremeters as appropriate
  function stir_config_constructor()
    implicit none
    type(stir_config_type) :: stir_config_constructor
    stir_config_constructor%name = "stir"
    stir_config_constructor%requiresIndex = .true.
  end function stir_config_constructor
end module antenna

