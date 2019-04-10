!> FIXME : Add documentation
module gs2_time

  implicit none

  private

  public :: user_dt, code_dt, update_time, code_dt_old
  public :: time_history_available
  public :: code_dt_prev1, code_dt_prev2
  public :: user_time, code_time, check_time_step_too_large
  public :: save_dt_min, save_dt, save_dt_cfl, write_dt
  public :: init_tstart, init_delt
  public :: code_dt_cfl, code_dt_min, user2code

  real, parameter :: dt_not_set = -1.0
  real :: user_dt, code_dt 
  real :: code_dt_prev1, code_dt_prev2
  real :: user_dt_cfl = dt_not_set
  real :: code_dt_cfl = dt_not_set

  real :: code_dt_min, user_dt_min

  ! added May 18, 2009 to take care of problems
  ! in exb_shear calculation after change in time step size
  real :: code_dt_old = 0.

! GGH points out that this initialization is not necessary (we think)
  real :: user_time = 0.
  real :: code_time = 0.

contains

  !> FIXME : Add documentation  
  subroutine init_tstart (tstart)
    implicit none
    real, intent (in) :: tstart

    user_time = tstart
    code_time = tstart

  end subroutine init_tstart

  !> Set the full timestep history.
  !! Inputs set the current, previous and last but one timesteps
  subroutine init_delt (delt, delt1, delt2)
    implicit none
    real, intent (in) :: delt, delt1, delt2
!
! delt* are user input, from the run_parameters module.
! In a perfect world, we could have a gs2_time namelist. 
! 
    user_dt = delt
    code_dt = delt
    code_dt_prev1 = delt1
    code_dt_prev2 = delt2
  end subroutine init_delt

  !> Check if the current time step is too big
  subroutine check_time_step_too_large(reset)
    implicit none
    logical, intent(out) :: reset
    reset=(code_dt>code_dt_cfl)
  end subroutine check_time_step_too_large

  !> FIXME : Add documentation  
  subroutine update_time
    implicit none
! for using a changing-timestep AB3 algorithm for the nonlinear term  
    code_dt_prev2 = code_dt_prev1
    code_dt_prev1 = code_dt


! MAB+CMR, 21/5/09: set code_dt_old to code_dt BEFORE any changes in timestep
    code_dt_old = code_dt
    code_time = code_time + code_dt
    user_time = user_time + user_dt
  

  end subroutine update_time

  !> FIXME : Add documentation
  subroutine save_dt_cfl (delt_cfl)
    implicit none
    real, intent (in) :: delt_cfl

    code_dt_cfl = delt_cfl
    user_dt_cfl = delt_cfl

  end subroutine save_dt_cfl

  !> FIXME : Add documentation  
  subroutine save_dt_min (dt_min)
    implicit none
    real, intent (in) :: dt_min

    user_dt_min = dt_min
    code_dt_min = dt_min

  end subroutine save_dt_min

  !> FIXME : Add documentation
  subroutine save_dt(delt)
    implicit none
    real, intent (in) :: delt

    code_dt = delt
    user_dt = delt
    
  end subroutine save_dt

  !> FIXME : Add documentation  
  subroutine write_dt
    implicit none
    if (user_dt_cfl > 0. .and. user_dt_cfl < 1.e7) &
         write(*,*) user_dt_cfl,' : ',user_dt
  end subroutine write_dt

  !> FIXME : Add documentation  
  subroutine user2code (usertime, codetime)
    implicit none
    real, intent (in) :: usertime
    real, intent (out) :: codetime
    codetime = usertime
  end subroutine user2code

  !> Reports how many previous time steps are available
  !! by checking if the previous code_dt values have been
  !! set to anything.
  function time_history_available()
    implicit none
    integer :: time_history_available
    if (code_dt_prev2 .ne. dt_not_set) then
       return 3
    else if (code_dt_prev1 .ne. dt_not_set) then
       return 2
    else
       return 1
    endif
  end function time_history_available
end module gs2_time
