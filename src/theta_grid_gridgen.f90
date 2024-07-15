!> FIXME : Add documentation
module theta_grid_gridgen
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: theta_grid_gridgen_init, finish_theta_grid_gridgen
  public :: gridgen_get_grids
  public :: wnml_theta_grid_gridgen

  public :: theta_grid_gridgen_config_type
  public :: set_theta_grid_gridgen_config
  public :: get_theta_grid_gridgen_config

  ! knobs
  integer :: npadd
  real :: alknob, epsknob, bpknob, extrknob, regrid_tension, tension
  real :: thetamax, deltaw, widthw
  logical :: skip_gridgen
  logical :: exist, initialized=.false.

  !> Used to represent the input configuration of theta_grid_gridgen
  type, extends(abstract_config_type) :: theta_grid_gridgen_config_type
     ! namelist : theta_grid_gridgen_knobs
     ! indexed : false
     !> Relative weighting of pitch-angle metric to \(\theta\) metric
     real :: alknob = 0.0
     !> Consider when the right grid point is equal to the target bmag.
     !>
     !> FIXME: What does this mean?
     real :: bpknob = 1.e-8
     !> Parameter for weighted resolution in theta. Each theta grid point
     !> contributes to the resolution metric according to the function
     !> $$ 1 + deltaw *[ 1 / (1 + [theta-thetamax]**2/widthw**2) +
     !>                  1 / (1 + [theta+thetamax]**2/widthw**2)] $$
     !> which has peaks at theta = +/- thetamax and deltaw controls the
     !> relative weighting of the theta dependent contribution.
     real :: deltaw = 0.0
     !> Maximum difference between neighbouring points for determining
     !> if a point is an extremum.
     real :: epsknob = 1e-5
     !> Used to set a "bonus" contribtion to resolution at B extrema with an
     !> even number of theta grid points. Those with an odd number of points
     !> and the assumed extrema at -pi have a metric of 1e20. Here extrknob
     !> can be used to bias the algorithm towards keeping extrema with an
     !> even number of points.
     real :: extrknob = 0.0
     !> Number of points between original grid points to oversample \(\theta\) by
     integer :: npadd = 2
     !> Tension to use in interpolating splines for regrid of geometrical quantities
     !> Defaults to [[theta_grid_gridgen_config_type:tension]] if not set.
     real :: regrid_tension = -1.0
     !> If true then skip gridgen call and instead just use the existing grid.
     !>
     !> @note This is primarily for debugging and testing. When active we are
     !> effectively forced to assumed that the input bmag is symmetric and
     !> monotonic. If this is not true then we should not trust the results.
     logical :: skip_gridgen = .false.
     !> Tension for \(B\) spline
     real :: tension = 1.0
     !> Parameter for weighted resolution in theta. Each theta grid point
     !> contributes to the resolution metric according to the function
     !> $$ 1 + deltaw *[ 1 / (1 + [theta-thetamax]**2/widthw**2) +
     !>                  1 / (1 + [theta+thetamax]**2/widthw**2)] $$
     !> which has peaks at theta = +/- thetamax.
     real :: thetamax = 0.0
     !> Parameter for weighted resolution in theta. Each theta grid point
     !> contributes to the resolution metric according to the function
     !> $$ 1 + deltaw *[ 1 / (1 + [theta-thetamax]**2/widthw**2) +
     !>                  1 / (1 + [theta+thetamax]**2/widthw**2)] $$
     !> which has peaks at theta = +/- thetamax and widthw controls the
     !> scale length of the peaks.
     real :: widthw = 1.0
   contains
     procedure, public :: read => read_theta_grid_gridgen_config
     procedure, public :: write => write_theta_grid_gridgen_config
     procedure, public :: reset => reset_theta_grid_gridgen_config
     procedure, public :: broadcast => broadcast_theta_grid_gridgen_config
     procedure, public, nopass :: get_default_name => get_default_name_theta_grid_gridgen_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_theta_grid_gridgen_config
  end type theta_grid_gridgen_config_type

  type(theta_grid_gridgen_config_type) :: theta_grid_gridgen_config

contains

  !> FIXME : Add documentation
  subroutine wnml_theta_grid_gridgen(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_gridgen_knobs"
    write (unit, fmt="(' npadd =    ',i4)") npadd
    write (unit, fmt="(' alknob =   ',e17.10)") alknob
    write (unit, fmt="(' epsknob =  ',e17.10)") epsknob
    write (unit, fmt="(' bpknob =   ',e17.10)") bpknob
    write (unit, fmt="(' extrknob = ',e17.10)") extrknob
    write (unit, fmt="(' regrid_tension =  ',e17.10)") regrid_tension
    write (unit, fmt="(' skip_gridgen =  ',L1)") skip_gridgen
    write (unit, fmt="(' tension =  ',e17.10)") tension
    write (unit, fmt="(' thetamax = ',e17.10)") thetamax
    write (unit, fmt="(' deltaw =   ',e17.10)") deltaw
    write (unit, fmt="(' widthw =   ',e17.10)") widthw
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_gridgen

  !> FIXME : Add documentation
  subroutine theta_grid_gridgen_init(theta_grid_gridgen_config_in)
    implicit none
    type(theta_grid_gridgen_config_type), intent(in), optional :: theta_grid_gridgen_config_in
    if (initialized) return
    initialized = .true.
    call read_parameters(theta_grid_gridgen_config_in)
  end subroutine theta_grid_gridgen_init

  !> FIXME : Add documentation
  subroutine finish_theta_grid_gridgen
    implicit none
    initialized = .false.
    call theta_grid_gridgen_config%reset()
  end subroutine finish_theta_grid_gridgen

  !> FIXME : Add documentation
  subroutine read_parameters(theta_grid_gridgen_config_in)
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    type(theta_grid_gridgen_config_type), intent(in), optional :: theta_grid_gridgen_config_in

    if (present(theta_grid_gridgen_config_in)) theta_grid_gridgen_config = theta_grid_gridgen_config_in

    call theta_grid_gridgen_config%init(name = 'theta_grid_gridgen_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    npadd  = theta_grid_gridgen_config%npadd
    alknob  = theta_grid_gridgen_config%alknob
    epsknob  = theta_grid_gridgen_config%epsknob
    bpknob  = theta_grid_gridgen_config%bpknob
    extrknob  = theta_grid_gridgen_config%extrknob
    regrid_tension  = theta_grid_gridgen_config%regrid_tension
    skip_gridgen = theta_grid_gridgen_config%skip_gridgen
    tension  = theta_grid_gridgen_config%tension
    thetamax  = theta_grid_gridgen_config%thetamax
    deltaw  = theta_grid_gridgen_config%deltaw
    widthw  = theta_grid_gridgen_config%widthw

    exist = theta_grid_gridgen_config%exist

    ! Handle special treatment
    if (regrid_tension < 0) regrid_tension = tension
  end subroutine read_parameters

  !> FIXME : Add documentation
  subroutine gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
       theta, bset, bmag, gradpar, gbdrift, gbdrift0, cvdrift, &
       cvdrift0, cdrift, cdrift0, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    use gridgen4mod, only: gridgen4_2
    use constants, only: pi
    use mp, only: mp_abort
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: theta
    real, dimension (nbset), intent (in out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    integer :: ntheta_old, ntgrid_old, nbset_old
    real, dimension (-ntgrid:ntgrid) :: thetasave
    real, dimension (ntheta+1) :: thetaold, thetanew
    real, dimension (ntheta+1) :: bmagold, bmagnew
    integer :: i
    logical, parameter :: debug=.false.
    if (debug) write(6,*) 'gridgen_get_grids'

    ! If we're skipping gridgen then we'd better
    ! make sure that we set the outputs. As we're not
    ! changing the theta grid the only thing we need to
    ! set is `bset`, which is just the unique bmag values
    ! on our grid. Note here we're assuming bmag is symmetric,
    ! ordered etc. By skipping gridgen we don't guarantee that
    ! this will be the case.
    if (skip_gridgen) then
       bset = bmag(-ntheta/2:0)
       return
    end if

    ntheta_old = ntheta
    ntgrid_old = ntgrid
    nbset_old = nbset

    thetasave = theta
    thetaold = theta(-ntheta/2:ntheta/2)
    bmagold = bmag(-ntheta/2:ntheta/2)

if (debug) write(6,*) 'gridgen_get_grids: call gridgen4_2'
    call gridgen4_2 (1,ntheta_old+1,thetaold,bmagold, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nbset,thetanew,bmagnew,bset)

    if (ntheta_old /= ntheta) then
       write(*,*) 'Error in theta_grid_gridgen?'
       write(*,*) 'ntheta_old = ',ntheta_old
       write(*,*) 'ntheta_new = ',ntheta
       write(*,*) 'Stopping this run would be wise so will now abort.'
       write(*,*) 'Try again with ntheta = ',ntheta_old + 2
       if(ntheta_old<ntheta)then
          write(*,*) 'ntheta_old<ntheta but code assumes ntheta<ntheta_old.'
       endif
       call mp_abort("Bad behaviour spotted in theta_grid_gridgen. Consider changing ntheta or set skip_gridgen = .true.", to_screen = .true.)
    end if

    ! interpolate to new grid
    ntgrid = ntheta/2 + (nperiod-1)*ntheta

    theta(-ntheta/2:ntheta/2-1) = thetanew(1:ntheta)
    theta(ntheta/2) = thetanew(1) + real(2)*pi
    bmag(-ntheta/2:ntheta/2-1) = bmagnew(1:ntheta)
    bmag(ntheta/2) = bmagnew(1)
    do i = 1, nperiod-1
       theta(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) &
            = thetanew(1:ntheta) + real(2*i)*pi
       theta(ntheta/2+i*ntheta) = thetanew(1) + real(2*(i+1))*pi
       theta(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) &
            = thetanew(1:ntheta) - real(2*i)*pi
       bmag(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) = bmagnew(1:ntheta)
       bmag( ntheta/2+i*ntheta) = bmagnew(1)
       bmag(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) = bmagnew(1:ntheta)
    end do

if (debug) write(6,*) 'gridgen_get_grids: call regrid'
!<DD>NOTE: Regrid assumes nnew<nold but doesn't check it. Do we need to?
!          This only packs the new (splined) data into the old array, it
!          doesn't actually resize the array. This resizing currently takes
!          place during finish_init call (look for eik_save). Would be safer
!          if we could make regrid actually reallocate/resize the array.
!<DD>NOTE: We only resize the arrays internal to theta_grid. This is what should be used
!          by the rest of the code, but anywhere where the geometry arrays are used directly
!          could lead to an error if these arrays are resized as we don't resize the geometry arrays.
    call regrid (ntgrid_old, thetasave, gradpar, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds2, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds21, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds22, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds23, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds24, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds24_noq, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, grho, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Rplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Zplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, aplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Rprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Zprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, aprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Bpol, ntgrid, theta)

if (debug) write(6,*) 'gridgen_get_grids: end'
  end subroutine gridgen_get_grids

  !> FIXME : Add documentation
  !!
  !! @note This routine assumes nnew<nold
  subroutine regrid (nold, x, y, nnew, xnew)
    use splines, only: new_spline, splint, delete_spline, spline
    implicit none
    integer, intent (in) :: nold
    real, dimension (-nold:nold), intent (in) :: x
    real, dimension (-nold:nold), intent (in out) :: y
    integer, intent (in) :: nnew
    real, dimension (-nnew:nnew), intent (in) :: xnew
    type (spline) :: spl
    integer :: i

    spl = new_spline (x(-nold:nold), y(-nold:nold), tension = regrid_tension)

    do i = -nnew, nnew
       y(i) = splint(xnew(i), spl)
    end do

    call delete_spline (spl)
  end subroutine regrid

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_theta_grid_gridgen_config(theta_grid_gridgen_config_in)
    use mp, only: mp_abort
    type(theta_grid_gridgen_config_type), intent(in), optional :: theta_grid_gridgen_config_in
    if (initialized) then
       call mp_abort("Trying to set theta_grid_gridgen_config when already initialized.", to_screen = .true.)
    end if
    if (present(theta_grid_gridgen_config_in)) then
       theta_grid_gridgen_config = theta_grid_gridgen_config_in
    end if
  end subroutine set_theta_grid_gridgen_config

#include "theta_grid_gridgen_auto_gen.inc"
end module theta_grid_gridgen
