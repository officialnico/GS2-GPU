!> A small module used to calculate ideal ballooning stability.
!! Based on the original ball program in geo/ball.f90 of GS2
!! but tweaked to integrate closer to full GS2 runs etc.
module ballstab
  use abstract_config, only: abstract_config_type
  implicit none

  private

  !Routines
  public :: init_ballstab, finish_ballstab
  public :: run_stability_check, is_unstable
  public :: write_stability_ascii_1d, write_stability_ascii_2d

  !Vars
  public :: stability

  public :: ballstab_config_type
  
  logical :: make_salpha
  integer :: n_shat, n_beta
  real :: shat_min, shat_max
  real :: beta_mul, beta_div
  real, dimension(:), allocatable :: shat_arr, beta_arr, dbdrho_arr
  integer, dimension(:,:), allocatable :: stability
  real :: diff

  logical :: initialised = .false. 
  character(len=14) :: namelist_name="ballstab_knobs"

  !> Used to represent the input configuration of ballstab
  type, extends(abstract_config_type) :: ballstab_config_type
     ! namelist : ballstab_knobs
     ! indexed : false     
     !> FIXME: Add documentation
     real :: beta_div = 1.0
     !> FIXME: Add documentation
     real :: beta_mul = 1.0
     !> FIXME: Add documentation
     real :: diff = 0.0
     !> FIXME: Add documentation
     logical :: make_salpha = .false.
     !> FIXME: Add documentation
     integer :: n_beta = 1
     !> FIXME: Add documentation
     integer :: n_shat = 1
     !> @note gets a smart default
     real :: shat_max = 0.0
     !> @note gets a smart default
     real :: shat_min = 0.0
   contains
     procedure, public :: read => read_ballstab_config
     procedure, public :: write => write_ballstab_config
     procedure, public :: reset => reset_ballstab_config
     procedure, public :: broadcast => broadcast_ballstab_config
  end type ballstab_config_type

  interface ballstab_config_type
     procedure :: ballstab_config_constructor
  end interface ballstab_config_type

  type(ballstab_config_type) :: ballstab_config
  
contains
  
  !> Initialise grids, other modules etc.
  subroutine init_ballstab
    use theta_grid, only: init_theta_grid
    implicit none

    !Other modules
    call init_theta_grid

    !This module
    call real_init_ballstab
  end subroutine init_ballstab

  !> Finalise this module
  subroutine finish_ballstab
    implicit none
    !>Deallocate
    call dealloc_arrays

    initialised = .false.
    call ballstab_config%reset()
  end subroutine finish_ballstab

  !> Initialise this module
  subroutine real_init_ballstab
    use mp, only: proc0
    implicit none

    if(initialised) return
    initialised = .true.

    !Get settings
    call read_parameters

    !Note only proc0 will be allowed to do any work
    if (.not.proc0) return
    
    !Verify parameters
    call verify_parameters

    !Create arrays
    call setup_arrays

  end subroutine real_init_ballstab

  !> Read parameters
  subroutine read_parameters(ballstab_config_in)
    use theta_grid, only: shat
    implicit none
    type(ballstab_config_type), intent(in), optional :: ballstab_config_in    
    logical :: exist

    if (present(ballstab_config_in)) ballstab_config = ballstab_config_in

    !Smart defaults
    if (.not.ballstab_config%is_initialised()) then    
       ballstab_config%shat_min = shat
       ballstab_config%shat_max = shat
    end if
    
    call ballstab_config%init(name = 'ballstab_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    beta_div = ballstab_config%beta_div
    beta_mul = ballstab_config%beta_mul
    diff = ballstab_config%diff
    make_salpha = ballstab_config%make_salpha
    n_beta = ballstab_config%n_beta
    n_shat = ballstab_config%n_shat
    shat_max = ballstab_config%shat_max
    shat_min = ballstab_config%shat_min

    exist = ballstab_config%exist
  end subroutine read_parameters

  !> Check parameters are consistent etc.
  subroutine verify_parameters
    use geometry, only: bishop, iflux
    use file_utils, only: error_unit
    implicit none
    integer :: ierr

    ierr = error_unit()

    !Disable scan in s_hat when incompatible bishop used
    if(n_shat.gt.1)then
       select case(bishop)
       case(2,3,4,5,8,9)
          !These are valid values
       case default
          if((bishop.ne.0).and.(iflux.ne.1)) then
             !This is also ok
          else
             !This is not ok
             n_shat=1
             write(ierr,'("ERROR : Cannot scan in shat for current combination of bishop and iflux")')
          endif
       end select
    endif

    !Disable scan in beta_prime when incompatible bishop used
    if(n_beta.gt.1)then
       select case(bishop)
       case(2,3,4,5,6,7,8,9)
          !These are valid values
       case default
          !This is not ok
          n_beta=1
          write(ierr,'("ERROR : Cannot scan in beta_prime for current combination of bishop and iflux")')
       end select
    endif

  end subroutine verify_parameters

  !> Allocate+populate arrays
  subroutine setup_arrays
    implicit none
    integer :: is, ib
    real :: ds, db, b_min, b_max, bp
    
    allocate(shat_arr(n_shat))
    allocate(beta_arr(n_beta))
    allocate(dbdrho_arr(n_beta))
    dbdrho_arr = 0. !This will be populated later

    allocate(stability(n_shat,n_beta))
    stability = 0

    ds = 0
    if(n_shat.gt.1) ds = (shat_max-shat_min)/(n_shat-1)
    do is=1,n_shat
       shat_arr(is)=shat_min+ds*(is-1)
    enddo

    db = 0
    call get_beta_prime(bp)
    b_max = bp*beta_mul
    b_min = bp/beta_div
    if(n_beta.gt.1) then
       db = (b_max-b_min)/(n_beta-1)
    end if
    do ib=1,n_beta
       beta_arr(ib)=b_min+db*(ib-1)
    enddo
  end subroutine setup_arrays

  !> Deallocate arrays
  subroutine dealloc_arrays
    implicit none
    if(allocated(shat_arr)) deallocate(shat_arr)
    if(allocated(beta_arr)) deallocate(beta_arr)
    if(allocated(dbdrho_arr)) deallocate(dbdrho_arr)
    if(allocated(stability)) deallocate(stability)
  end subroutine dealloc_arrays

  !> Write a namelist based on current values
  subroutine wnml_ballstab(unit)
    implicit none
    integer, intent(in) :: unit

    !Header
    write(unit,'()')
    write(unit,'("&",A)') trim(adjustl(namelist_name))

    !Values
    write(unit,'(" make_salpha = ",L1)') make_salpha
    write(unit,'(" n_shat      = ",I0)') n_shat
    if(n_shat.gt.1) then
       write(unit,'(" shat_min    = ",F12.5)') shat_min
       write(unit,'(" shat_max    = ",F12.5)') shat_max
    endif
    write(unit,'(" n_beta      = ",I0)') n_beta

    !Footer
    write(unit,'("/")')
  end subroutine wnml_ballstab

  !> Check value consistency etc.
  subroutine check_ballstab(report_unit)
    implicit none
    integer, intent(in) :: report_unit
  end subroutine check_ballstab

  !> Gets the current value of shat
  subroutine get_shat(shat_out)
    use geometry, only: shat
    implicit none
    real, intent(out) :: shat_out
    shat_out = shat
  end subroutine get_shat

  !> Sets the value of shat
  subroutine set_shat(shat_in)
    use geometry, only: s_hat_input, shat
    implicit none
    real, intent(in) :: shat_in
    !Note for bishop = 6,7 or not in 2-9 and iflux==1 cannot change shat
    !Should flag this somewhere (during init)
    s_hat_input = shat_in
    shat = shat_in
  end subroutine set_shat

  !> Gets the current value of beta_prime (or equivalent var)
  !! Really just gets the variable that we can use to control 
  !! the pressure gradient.
  subroutine get_beta_prime(bp_out)
    use geometry, only: bishop, p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
    implicit none
    real, intent(out) :: bp_out

    select case (bishop)
    case(2)
       bp_out = p_prime_input
    case(3)
       bp_out = invLp_input
    case(4)
       bp_out = beta_prime_input
    case(5)
       bp_out = alpha_input
    case(6)
       bp_out = beta_prime_input
    case(7)
       bp_out = dp_mult
    case(8)
       bp_out = beta_prime_input
    case(9)
       bp_out = beta_prime_input
    case default
       bp_out = 0.0
       print*,"ERROR: For bishop not in 2-9 cannot modify the pressure gradient used"
    end select
  end subroutine get_beta_prime

  !> Sets the current value of beta_prime (or equivalent var)
  !! Really just sets the variable that we can use to control 
  !! the pressure gradient.
  subroutine set_beta_prime(bp_in)
    use geometry, only: bishop, p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
    implicit none
    real, intent(in) :: bp_in

    select case (bishop)
    case(2)
       p_prime_input = bp_in
    case(3)
       invLp_input = bp_in
    case(4)
       beta_prime_input = bp_in
    case(5)
       alpha_input = bp_in
    case(6)
       beta_prime_input = bp_in
    case(7)
       dp_mult = bp_in
    case(8)
       beta_prime_input = bp_in
    case(9)
       beta_prime_input = bp_in
    case default
       print*,"ERROR: For bishop not in 2-9 cannot modify the pressure gradient used"
    end select
  end subroutine set_beta_prime

  !> Runs the stability check scan
  subroutine run_stability_check
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: stab
    integer :: is, ib

    !Exit if not initialised
    if(.not.initialised) return

    !Check stability over range
    do ib=1,n_beta
       do is=1,n_shat
          call check_stability(shat_arr(is),beta_arr(ib),ib,stab)
          stability(is,ib)=stab
       enddo
    enddo
  end subroutine run_stability_check

  !> Routine to write out stability data to 1D ascii file
  subroutine write_stability_ascii_1d
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: is, ib, unit

    call open_output_file(unit,".ballstab_1d")
    write(unit,'(4(A12," "))') "shear", "beta_prime", "dbetadrho", "stability"
    do ib=1,n_beta
       do is=1,n_shat
          write(unit,'(3(F12.4," "),I12)') shat_arr(is), beta_arr(ib), dbdrho_arr(ib), stability(is,ib)
       enddo
    enddo
    call close_output_file(unit)
  end subroutine write_stability_ascii_1d

  !> Routine to write out stability data to 2D ascii file + 1d axis data
  subroutine write_stability_ascii_2d
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: is, ib, unit

    !/2D -- data
    call open_output_file(unit,".ballstab_2d")
    do ib=1,n_beta
       do is=1,n_shat
          write(unit,'(I0," ")',advance="no") stability(is,ib)
       enddo
       write(unit,'()') !Move to the next line
    enddo

    !/2D -- axis
    call open_output_file(unit,".ballstab_shat")
    !Commented lines here can be used to give row output
    do is=1,n_shat
       !write(unit,'(F9.2," ")',advance="no") shat_arr(is)
       write(unit,'(F9.2)') shat_arr(is)
    enddo
    !write(unit,'()')
    call close_output_file(unit)

    call open_output_file(unit,".ballstab_bp")
    do ib=1,n_beta
       write(unit,'(F12.4)') beta_arr(ib)
    enddo
    call close_output_file(unit)

    call open_output_file(unit,".ballstab_dbdr")
    do ib=1,n_beta
       write(unit,'(F12.4)') dbdrho_arr(ib)
    enddo
    call close_output_file(unit)
  end subroutine write_stability_ascii_2d

  !> Test if given shat/beta is unstable
  subroutine check_stability(shat_in, beta_prime_in, ib, iunstable, restore)
    !Note iunstable>0 means unstable - Could use a logical but size may be interesting
    use geometry, only: dbetadrho, eikcoefs
    implicit none
    real, intent(in) :: shat_in, beta_prime_in
    integer, intent(in) :: ib
    integer, intent(out) :: iunstable
    logical, intent(in), optional :: restore
    real :: shat_bak, beta_prime_bak
    
    !Store initial state
    call get_shat(shat_bak)
    call get_beta_prime(beta_prime_bak)

    !Setup geometry
    !NOTE: As how we can influence beta_prime (and shat) varies with bishop
    !we actually end up setting different parameters in different cases. The
    !physical parameter dbetadrho should really be the relevant parameter.
    call set_beta_prime(beta_prime_in)
    call set_shat(shat_in)
    call eikcoefs !Recalculate geometry terms

    !Here we store the value of dbetadrho, this should be bishop value independent
    dbdrho_arr(ib)=dbetadrho

    !Test stability
    iunstable=is_unstable()

    !Restore geometry if requested
    if(present(restore))then
       if(restore) then
          call set_beta_prime(beta_prime_bak)
          call set_shat(shat_bak)
          call eikcoefs
       endif
    endif
  end subroutine check_stability

  !> Test stability of current system -- return integer
  integer function is_unstable()
    !Note iunstable>0 means unstable - Could use a logical but size may be interesting
    !Assumes a call has been made to eikcoefs so that geometry terms are setup correctly
!<DD> NEED TO IMPROVE VARIABLE NAMES
    use geometry, only: gds2, gradpar, bmag, cvdrift, dbetadrho, theta
    implicit none
    real, dimension(:), allocatable :: g, c, psi, c1, c2, ch, g1, g2, gh, delthet
    integer :: ig, ntgrid
    real :: cflmax, cflmin, one_m_diff, psi_prime
    logical, parameter :: debug = .false.

    !Initialise
    is_unstable=0
    one_m_diff=1.-diff

    !Note we don't use theta_grid:ntgrid as in the case that
    !gridgen resizes the theta grid used by GS2 we don't end up
    !resizing any of the geometry module arrays (that we use directly here)
    !Really we would like to use everything from theta_grid directly but we'd
    !need some way to force theta_grid to update it's internal arrays when we
    !change parameters
    ntgrid=ubound(theta,1)
 
    !Calculate geometry arrays
    allocate(delthet(-ntgrid:ntgrid))
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)

    allocate(g(-ntgrid:ntgrid),c(-ntgrid:ntgrid))
    g = gds2*gradpar/bmag
    c = -0.5*dbetadrho*cvdrift/(bmag*gradpar)

    allocate(g1(-ntgrid:ntgrid), c1(-ntgrid:ntgrid))
    allocate(g2(-ntgrid:ntgrid), c2(-ntgrid:ntgrid))
    allocate(gh(-ntgrid:ntgrid), ch(-ntgrid:ntgrid))
    allocate(psi(-ntgrid:ntgrid))

    do ig=-ntgrid+1, ntgrid
       ch(ig)=0.5*(c(ig)+c(ig-1))
       gh(ig)=0.5*(g(ig)+g(ig-1))
    enddo

    if(debug) then
       cflmax = delthet(-ntgrid)**2
       cflmax = cflmax*(c(-ntgrid+1)+c(-ntgrid))/(g(-ntgrid+1)+g(-ntgrid))
       cflmin = cflmax
       do ig=-ntgrid+1,ntgrid
          cflmax=max(cflmax,delthet(ig-1)**2*ch(ig)/gh(ig))
          cflmin=min(cflmin,delthet(ig-1)**2*ch(ig)/gh(ig))
       enddo
       write(6,'("|cflmax| = ",F12.5)') abs(cflmax)
       write(6,'("|cflmin| = ",F12.5)') abs(cflmin)
       write(6,'("These should be much less than 1")') !Could actually do a quantitative check
    endif

    !Calculate coefficients
    do ig=-ntgrid+1,ntgrid-1
       c1(ig) = -delthet(ig)*(one_m_diff*c(ig)+0.5*diff*ch(ig+1))&
            -delthet(ig-1)*(one_m_diff*c(ig)+0.5*diff*ch(ig))&
            -delthet(ig-1)*0.5*diff*ch(ig)
       c1(ig)=0.5*c1(ig)
    enddo

    !Calculate coefficients
    do ig=-ntgrid+1,ntgrid
       c2(ig) = -0.25*diff*ch(ig)*delthet(ig-1)
       g1(ig) = gh(ig)/delthet(ig-1)
       g2(ig) = 1.0/(0.25*diff*ch(ig)*delthet(ig-1)+gh(ig)/delthet(ig-1))
    enddo

    !Set boundary values
    psi(-ntgrid)=0.
    psi(-ntgrid+1)=delthet(-ntgrid)
    psi_prime=psi(-ntgrid+1)/g2(-ntgrid+1)!-g1(-ntgrid+1)*psi(-ntgrid) !Second term zero by defn

    !Integrate along in theta
    do ig=-ntgrid+1,ntgrid-1
       psi_prime=psi_prime+c1(ig)*psi(ig)+c2(ig)*psi(ig-1)
       psi(ig+1)=(g1(ig+1)*psi(ig)+psi_prime)*g2(ig+1)
    enddo

    !Find how many times psi crosses axis || Actually how many times it hits zero
    !If we really do just want to count crossings then need to replace le with lt
    do ig=-ntgrid+1,ntgrid-1
       if(psi(ig)*psi(ig+1) .le. 0 ) is_unstable=is_unstable+1
    enddo

    !Probably want to think about writing out psi if debug

    !Clean up
    deallocate(g,c,g1,c1,g2,c2,gh,ch)
  end function is_unstable

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the ballstab_knobs namelist and populates the member variables
  subroutine read_ballstab_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(ballstab_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: n_beta, n_shat
    logical :: make_salpha
    real :: beta_div, beta_mul, diff, shat_max, shat_min

    namelist /ballstab_knobs/ beta_div, beta_mul, diff, make_salpha, n_beta, n_shat, shat_max, shat_min

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    beta_div = self%beta_div
    beta_mul = self%beta_mul
    diff = self%diff
    make_salpha = self%make_salpha
    n_beta = self%n_beta
    n_shat = self%n_shat
    shat_max = self%shat_max
    shat_min = self%shat_min

    ! Now read in the main namelist
    if (self%requiresIndex) then
       call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
       read(in_file, nml = ballstab_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%name, exist)
       if (exist) read(in_file, nml = ballstab_knobs)
    endif

    ! Now copy from local variables into type members
    self%beta_div = beta_div
    self%beta_mul = beta_mul
    self%diff = diff
    self%make_salpha = make_salpha
    self%n_beta = n_beta
    self%n_shat = n_shat
    self%shat_max = shat_max
    self%shat_min = shat_min

    self%exist = exist
  end subroutine read_ballstab_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_ballstab_config(self, unit)
    implicit none
    class(ballstab_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("beta_div", self%beta_div, unit_internal)
    call self%write_key_val("beta_mul", self%beta_mul, unit_internal)
    call self%write_key_val("diff", self%diff, unit_internal)
    call self%write_key_val("make_salpha", self%make_salpha, unit_internal)
    call self%write_key_val("n_beta", self%n_beta, unit_internal)
    call self%write_key_val("n_shat", self%n_shat, unit_internal)
    call self%write_key_val("shat_max", self%shat_max, unit_internal)
    call self%write_key_val("shat_min", self%shat_min, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_ballstab_config

  !> Resets the config object to the initial empty state
  subroutine reset_ballstab_config(self)
    class(ballstab_config_type), intent(in out) :: self
    type(ballstab_config_type) :: empty
    select type (self)
    type is (ballstab_config_type)
       self = empty
    end select
  end subroutine reset_ballstab_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_ballstab_config(self)
    use mp, only: broadcast
    implicit none
    class(ballstab_config_type), intent(in out) :: self
    call broadcast(self%beta_div)
    call broadcast(self%beta_mul)
    call broadcast(self%diff)
    call broadcast(self%make_salpha)
    call broadcast(self%n_beta)
    call broadcast(self%n_shat)
    call broadcast(self%shat_max)
    call broadcast(self%shat_min)

    call broadcast(self%exist)
  end subroutine broadcast_ballstab_config
  
  !> Sets the name and requiresIndex paremeters as appropriate
  function ballstab_config_constructor()
    implicit none
    type(ballstab_config_type) :: ballstab_config_constructor
    ballstab_config_constructor%name = "ballstab_knobs"
    ballstab_config_constructor%requiresIndex = .false.
  end function ballstab_config_constructor
end module ballstab
