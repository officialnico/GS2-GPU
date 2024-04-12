!> Given B{theta_i} attempt to produce an optimal set of theta values
!> with the same size as the input {theta_i}. The output {theta} set
!> should include all extrema (global and local) and as many other
!> grid points as required to meet the requested number.
!>
!> The general approach is to:
!>    1. Generate a periodic interpolating spline of B(theta).
!>    2. Find all extrema (essential points).
!>    3. Evaluate the spline on a finer theta mesh, controlled by npadd.
!>    4. Construct a set of B values consisted of all extrema and values on the fine mesh.
!>    5. For each unique B value find all values of theta for which this B appears (setset).
!>    6. Eliminate B values, removing the corresponding values of theta, until we have
!>       the required number of theta grid points. We choose which B values to remove based
!>       on (user controllable) theta and lambda grid resolution metrics. We strongly bias
!>       towards keeping extrema.
!>
!> @note Instead of splining and evaluating on a finer mesh we could instead pass a higher
!> resolution B(theta) calculated in geometry when ntheta_geometry > ntheta.
!>
!> @note Instead of packing extra theta points in step 3, it might make more sense to
!> finely sample B between the global min and max values?
!>
!> MRH 15/11/2017
!> Previous versions of gridgen4_2  produced different theta grids on different machines
!> In this version of gridgen4mod I have corrected this bug by making changes in the following subroutines (in gridgen4_2):
!> heapsort- modified the inequalities in the sorting algorithm for reasons discussed below.
!> get_thetares-corrected this routine so that thetares is 0 when thetares is computed using
!> (theta (i) - theta (j)) and theta(i) equals theta (j) to machine error, and the output should be identically 0.
!> get_lambdares-corrected this routine in the same manner as get_thetares, the output is 0 when it should be identically 0.
!> gg4remove- modification of the inequalities used to choose which points to delete from the superset of theta points.
!>
!> The function of this module seems to be to optimise the choice of theta grid points.
!> It does this by taking the input magnetic field and theta grid (usually generated in the geometry module)
!> and spline fitting the magnetic field on a fine-mesh theta grid.
!>
!> By using the metrics defined in get_lambdares and get_thetares
!> the subroutine gg4remove removes points from the fine-mesh theta grid iteratively until the required number of points remain.
!> The theta grid pair with the lowest value for the metric is deleted at each iteration.
!> The reason for the bug was that these metrics were calculated in a sloppy manner,
!> and the routine for finding the lowest value in the metric did not allow for the
!> case that the metrics were identical up to machine precision:
!> -In the 1st iteration of gg4remove all of the metrics should be identically 0 but were not
!> -The 1st theta grid pair were deleted based on a decision made with machine precision sized values of the metric-this led to random behaviour
!> -When all duplicate values of the theta grid pairs
!> ----- (theta grid superset contains duplicate values of the theta grid pairs for reasons unbeknownst to me)
!> had all been deleted, because the fine-mesh theta grid had almost equally spaced points
!> the values for the metrics for the theta grid pairs near theta =0 were almost exactly equal.
!> -Depending on the machine, the metric would be slightly lower for one pair of points than another
!> and therefore the choice of which pair of points to delete 1st from the non-duplicate set was random.
!> This behaviour led to different theta grids between different machines.
!>
!> I fixed the above bug by biasing the gg4remove procedure to remove the 1st theta grid pair in array element order
!> (of the sorted theta grid superset array) in the case that the metrics were equal up to machine error.
!> The variable precision_bound controls this behaviour, precision_bound equal 0.0 returns this module to its original behaviour.
!> Finite but small precision_bound is enough to remove machine dependence.
!> I biased heapsort in a similar manner to enable debugging - the sorting algorithm
!> was also machine dependent because of the duplicate theta grid points in the superset.
!>
!> run_name.gridgen.200 contains the debug output from this module
!> To see verbose output ensure that debug_output =.true.and verbose_debug_output =.true.
!>
!> @todo : This file/module needs tidying/pruning. More consistent usage of tolerance comparisons would help
!> Testing may be tricky but would be very helpful.
module gridgen4mod
  implicit none
  private
  public :: gridgen4_1, gridgen4_2
contains

  !> FIXME : Add documentation
  !> @note Currently unused
  subroutine gridgen4read (filename,ntheta,nperiod,ntgrid,nlambda,theta,alambda, &
       gbdrift,gradpar,cvdrift,gds2,bmag,gds21,gds22,cvdrift0,gbdrift0)
    use mp, only: mp_abort
    use constants, only: run_name_size
    implicit none
    character(*), intent (in) :: filename
    integer, intent (in out) :: ntheta, nperiod, ntgrid, nlambda
    real, dimension (-ntgrid:ntgrid), intent(out) :: theta
    real, dimension (nlambda), intent(out) :: alambda
    real, dimension (-ntgrid:ntgrid), intent(out) :: gbdrift, gradpar
    real, dimension (-ntgrid:ntgrid), intent(out) :: cvdrift, gds2, bmag
    real, dimension (-ntgrid:ntgrid), intent(out) :: gds21, gds22
    real, dimension (-ntgrid:ntgrid), intent(out) :: cvdrift0, gbdrift0
    integer :: nthetain, nperiodin, ntgridin, nlambdain, i, unit
    character(run_name_size) :: line

    open (newunit=unit, file=filename, status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) nlambdain
    if (nlambdain > nlambda) then
       print *, "grid.out:nlambda > nlambda: ", nlambdain, nlambda
       call mp_abort("grid.out:nlambda > nlambda: ")
    end if
    nlambda = nlambdain
    read (unit=unit, fmt="(a)") line
    do i = 1, nlambda
       read (unit=unit, fmt=*) alambda(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) ntgridin, nperiodin, nthetain
    if (ntgridin > ntgrid) then
       print *, "grid.out:ntgrid > ntgrid: ", ntgridin, ntgrid
       call mp_abort("grid.out:ntgrid > ntgrid: ")
    end if
    ntgrid = ntgridin
    nperiod = nperiodin
    ntheta = nthetain
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gbdrift(i), gradpar(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gds21(i), gds22(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
    end do
    
    close (unit=unit)
  end subroutine gridgen4read

  !> FIXME : Add documentation  
  subroutine gridgen4 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw, &
       ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    implicit none
    integer, intent (in) :: nbmag, n
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw
    integer, intent (in out) :: ntheta, nlambda
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nlambda), intent (out) :: alambdagrid
    
    call gridgen4_1 (n,nbmag,thetain,bmagin, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,1.0, &
         ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
  end subroutine gridgen4

  !> FIXME : Add documentation  
  subroutine gridgen4_1 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
       ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: nbmag
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw, tension
    integer, intent (in out) :: ntheta, nlambda
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nlambda), intent (out) :: alambdagrid
    
    call gridgen4_2 (n,nbmag,thetain,bmagin, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    alambdagrid(1:nlambda) = 1.0/alambdagrid(1:nlambda)
  end subroutine gridgen4_1

  !> FIXME : Add documentation  
  subroutine gridgen4_2 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
       ntheta,nbset,thetagrid,bmaggrid,bset)
    use constants, only: pi, twopi
    use mp, only: proc0
    use splines, only: periodic_spline
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: nbmag
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw, tension
    integer, intent (in out) :: ntheta, nbset
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nbset), intent (out) :: bset
    real :: npi
    integer, parameter :: maxiter = 100
    !> MRH verbose_debug_output prints out the resolution_metric to
    !> debug_unit for every iteration of gg4remove to debug the point
    !> using procedure
    logical :: debug_output = .true., verbose_debug_output = .false.
    !> MRH a number used to take into account machine error in
    !> floatingpoint arithmetic
    real, parameter :: precision_bound = 1.0e3*epsilon (0.0)
    type(periodic_spline) :: bmagspl

    integer :: nstart
    real, allocatable, dimension (:) :: thetastart, bmagstart
    logical, allocatable, dimension (:) :: essentialstart
    real, allocatable, dimension (:) :: bextr, thetaextr, bprime

    integer :: nset, nsetset
    real, allocatable, dimension (:) :: thetasetset, bmagset
    !> MRH for debugging
    real, allocatable, dimension (:) :: thetasetset_sorted, bmagset_sorted
    integer, allocatable, dimension (:) :: ibmagsetset, icollsetset
    logical, allocatable, dimension (:) :: essentialset
    
    integer, allocatable, dimension (:) :: ithetasort, ibmagsort
    real, allocatable, dimension (:) :: thetares, alambdares
    
    integer :: nthetaout, nlambdaout, itheta, debug_unit
    logical, parameter :: debug=.false.!.true.
    npi = n * twopi

    ! Disable debug output when iproc /= 0 and make sure debug_output
    ! active if verbose_debug_output
    verbose_debug_output = verbose_debug_output .and. proc0
    debug_output = (debug_output .or. verbose_debug_output) .and. proc0

if (debug) write(6,*) "gridgen4_2: call gg4init"
    call gg4init
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nbmag, thetain, bmagin, "gg4init - input grid - thetain bmagin ")

! 1 Set up starting grid.
if (debug) write(6,*) "gridgen4_2: call gg4start"
    call gg4start (n)
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nstart, thetastart, thetastart, "gg4start - starting grid - thetastart thetastart")
! MRH thetastart is the array containing the theta grid on a finer grid, after a spline fit of the 
! MRH magnetic field has been done in gg4init and some obscure manipulations in gg4start
! MRH quirkily the 1st element is - pi and the 2nd element is 0, the rest of the elements increase monotonically from - pi to pi

! MRH bmagstart Contains a spline fit of the magnetic field, although a different one to gg4init
! MRH and is on the same grid as thetastart with the -pi and 0 elements as the 1st 2 elements. 

! MRH essentialstart is an array containing logicals indicating which elements of the thetastart grid are essential (it's the 1st 2 only)
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nstart, bmagstart, bmagstart, "gg4start - starting field - bmagstart bmagstart")
    call gg4debugl (nstart, essentialstart, "gg4start - essentialstart")
!MRH
! 2 Collect all bounce points associated with starting grid.
if (debug) write(6,*) "gridgen4_2: call gg4collect"
    call gg4collect
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nset, bmagset, bmagset, " gg4collect - bmagset bmagset")
if (debug) write(6,*) "gridgen4_2: call gg4debugi"
    call gg4debugi (nsetset,thetasetset,ibmagsetset,icollsetset," gg4collect - thetasetset ibmagsetset icollsetset")

! MRH bmagset seems not to differ from bmagstart

! MRH thetasetset contains paired elements of the theta grid with the same magnetic field strength except 
! MRH for the 1st 2 elements which are still -pi and 0
! MRH note that this grid contains some redundancy, every theta point appears 
! MRH in a (- theta, theta) pair and these pairs appear twice.

! MRH ibmagsetset is an array that is paired with thetasetset and contains 
! MRH the index of bmagset corresponding to the theta value stored in thetasetset
! MRH so thetasetset (i), bmagset(ibmagset(i)) are the paired values of (theta, B(theta)) 

! MRH icollsetset seems to be a debugging label which is not used outside of gg4collect

! 3 Sort collected grids.
if (debug) write(6,*) "gridgen4_2: call gg4sort"
    call gg4sort
! MRH ithetasort is an array which contains an index such that for a given index s 
! MRH of the monotonically increasing theta grid, the grid value for 
! MRH the monotonically increasing theta growth can be attained by 
! MRH thetasetset_sorted(s) = thetasetset(ithetasort(s))

! MRH ibmagsort is the corresponding array for the magnetic field sort
! MRH thetasetset_sorted and bmagset_sorted are the monotonically increasing sorted arrays for debugging
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debugi (nsetset, thetasetset, ithetasort, ithetasort, "gg4sort - thetasetset ithetasort ithetasort")
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debugi (nset, bmagset, ibmagsort, ibmagsort, "gg4sort - bmagset ibmagsort ibmagsort")
if (debug) then
    allocate (thetasetset_sorted(nsetset))
    allocate (bmagset_sorted(nset))
    do itheta=1,nsetset
        thetasetset_sorted(itheta) = thetasetset(ithetasort(itheta))
    enddo
    do itheta=1,nset
        bmagset_sorted(itheta) = bmagset(ibmagsort(itheta))
    enddo

    call gg4debug (nset, bmagset_sorted, bmagset, "gg4sort - bmagset_sorted bmagset")

    call gg4debug (nsetset,thetasetset_sorted,thetasetset,"gg4sort - thetasetset_sorted thetasetset")
    deallocate(thetasetset_sorted)
    deallocate(bmagset_sorted)
! MRH essentialset is an array containing logicals indicating which elements of the thetastart grid are essential (it's the 1st 2 only)    
    call gg4debugl (nstart, essentialset, "gg4sort - essentialset")
endif
!MRH
! 4 Calculate spacing and resolution metrics.
if (debug) write(6,*) "gridgen4_2: call gg4metrics"
    call gg4metrics
!MRH
! MRH are the metrics used to determine which points of the larger fine grid (thetasetset-with redundant points)
! MRH to keep as one reduces to the required grid size
! MRH the larger the value of the metric the more valuable the point is
! MRH A pair of points is deleted when it is the least valuable member of the set
! MRH the process of reducing to less fine grid is iterative and pairs of points are deleted one at a time
! MRH The metrics are recalculated after a pair of points is deleted within gg4remove
! MRH Points are deleted by setting their corresponding magnetic field strength bmagset = 0.0 
! MRH so that it is ignored when reconstructing the grid and field in gg4results
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nset, thetares, alambdares, "gg4metrics - thetares alambdares")
!MRH
! 5 Remove sets until the number of points is small enough.
if (debug) write(6,*) "gridgen4_2: call gg4metrics"
    call gg4remove

! 6 Build output grids.
if (debug) write(6,*) "gridgen4_2: call gg4results"
    call gg4results (n)
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nthetaout, thetagrid, bmaggrid, "gg4results - output grid - thetagrid thetagrid")
if (debug) write(6,*) "gridgen4_2: call gg4debugi"
    call gg4debug (nlambdaout, bset, bset, "gg4results - output 1/lambda - bset bset")

    ntheta = 2*(nthetaout/2)
    nbset = nlambdaout

! 7 Clean up.
if (debug) write(6,*) "gridgen4_2: call gg4finish"
    call gg4finish

  contains
    
    !> FIXME : Add documentation  
    subroutine gg4init
      use splines, only: new_periodic_spline
      use file_utils, only: open_output_file
      implicit none
      integer :: i
      bmagspl = new_periodic_spline(thetain, bmagin, npi, &
           drop_last_point = .true., tension = tension)

      nstart = 0
      nset = 0
      nsetset = 0

      if (debug_output) then
         call open_output_file(debug_unit,".gridgen.200")
         write (unit=debug_unit, fmt=*) "nbmag=", nbmag
         write (unit=debug_unit, fmt=*) "thetain,bmagin="
         do i=1,nbmag
            write (unit=debug_unit, fmt=*) thetain(i), bmagin(i)
         end do
         write (unit=debug_unit, fmt=*) "alknob=", alknob
         write (unit=debug_unit, fmt=*) "epsknob=", epsknob
         write (unit=debug_unit, fmt=*) "bpknob=", bpknob
         write (unit=debug_unit, fmt=*) "extrknob=", extrknob
         write (unit=debug_unit, fmt=*) "thetamax=", thetamax
         write (unit=debug_unit, fmt=*) "deltaw,widthw=", deltaw, widthw
         write (unit=debug_unit, fmt=*) "tension=", tension
         write (unit=debug_unit, fmt=*) "ntheta,nlambda=", ntheta, nbset
      end if
      if (debug) write(6,*) 'gg4init: end'
    end subroutine gg4init

    !> FIXME : Add documentation    
    subroutine gg4finish
      use mp, only: mp_abort
      use file_utils, only: close_output_file
      use splines, only: delete_periodic_spline
      implicit none
      
      if (debug_output) call close_output_file(debug_unit)
      call delete_periodic_spline(bmagspl)
      deallocate (thetastart, bmagstart, essentialstart)
      deallocate (bextr,thetaextr,bprime)

      deallocate (thetasetset, bmagset, essentialset)
      deallocate (ibmagsetset, icollsetset)
      deallocate (ithetasort, ibmagsort)
      deallocate (thetares, alambdares)

      if (mod(nthetaout,2) /= 1) then
         print *, "gridgen4_1:gg4results:nthetaout=", nthetaout
         print *, "nthetaout is not odd, so there is a problem with gridgen"
         print *, "This can sometimes happen if this file is not compiled with promotion to double precision."
         call mp_abort('gridgen has finished with an odd ntheta. This suggests something has gone wrong so abort. Consider changing the gridgen input flags or ntheta.', .true.)
      end if
    end subroutine gg4finish

    !> FIXME : Add documentation    
    subroutine gg4start (n)
      implicit none
      integer, intent (in) :: n
      integer :: i, j, nextr
      real :: thetal, thetar, theta0, bprime0, dtheta
      logical :: skip_point
      allocate(bextr(nbmag-1), thetaextr(nbmag-1), bprime(nbmag))
      bextr = 0.0 ; thetaextr = 0.0

      ! Evaluate dB/dtheta on the input theta grid from the spline of B
      do i = 1, nbmag
         bprime(i) = bmagp(thetain(i))
      end do

      ! Find all extrema. - Note we skip the first point and later
      ! assume that this is an extrema (at theta = -pi*n). We _could_
      ! just make the loop go over i = 1, nbmag - 1 instead.
      ! Here we start with one extrema corresponding to the assumed
      ! extrema at the start of the grid
      nextr = 1
      skip_point = .false.
      do i = 2, nbmag-2
         if (skip_point) then
            skip_point = .false.
            cycle
         end if

         if (bprime(i) == 0.0) then
            ! Exactly zero derivative on grid
            bextr(i) = bmagin(i)
            thetaextr(i) = thetain(i)
            nextr = nextr + 1
         else if (bprime(i) * bprime(i+1) < 0.0 ) then
            ! Change in derivative - extrema between grid points
            thetal = thetain(i) ; thetar = thetain(i+1)
            !Bisect to find extrema. Note we don't check actually found
            do j = 1, maxiter
               theta0 = 0.5 * (thetal + thetar)
               bprime0 = bmagp(theta0)
               if (bprime0 * bprime(i) > 0.0) then
                  !If no B' sign change set left theta bound to mid-point
                  thetal = theta0
               else if (bprime0 * bprime(i) < 0.0) then
                  !If B' sign change set right theta bound to mid-point
                  thetar = theta0
               else
                  exit
               end if
            end do
            !If extrema close enough to bounding grid points just
            !say derivative is zero on the grid point.
            if (theta0 - thetain(i) <= epsknob) then
               theta0 = thetain(i)
               bprime(i) = 0.0
            else if (thetain(i+1)-theta0 <= epsknob) then
               ! @Note: If this triggers then on the next iteration the bprime(i)==0 test
               ! should pass and we will end up adding this extremum twice, so mark as
               ! needing to skip next iteration.
               skip_point = .true.
               theta0 = thetain(i+1)
               bprime(i+1) = 0.0
            end if
            thetaextr(i) = theta0
            bextr(i) = bmagint(theta0)
            nextr = nextr + 1
         end if
      end do

      ! Collect -pi, all local extrema, original grid points, points between
      ! original grid points
      nstart = nextr + (nbmag - 1) * (1 + npadd)

      allocate (thetastart(nstart), essentialstart(nstart), bmagstart(nstart))
      essentialstart(:nextr) = .true.
      essentialstart(nextr+1:) = .false.
      thetastart(1) = -pi*n !Should this be thetain(1)?
      thetastart(2:nextr) = pack(thetaextr,bextr /= 0.0)
      
      nstart = nextr
      do i = 1, nbmag-1
         dtheta = (thetain(i+1)-thetain(i))
         do j = 0, npadd
            ! For j=0 add thetagrid(i) point, for j > 0 try to add extra theta grid
            ! points inbetween thetagrid(i) and thetagrid(i+1). Only actually keep
            ! the points if they are more than epsknob far away from all previously
            ! added points. As such epsknob sets the minimum theta grid spacing
            ! allowed here.
            thetastart(nstart+1) &
                 = thetain(i) + dtheta * real(j)/real(npadd+1)
            if (all(abs(thetastart(nstart+1)-thetastart(1:nstart)) > epsknob)) then
               nstart = nstart + 1
            end if
         end do
      end do

      ! Get B on the new theta grid point through spline interpolation
      do i = 1, nstart
         bmagstart(i) = bmagint(thetastart(i))
      end do
    end subroutine gg4start

    subroutine gg4collect
      use mp, only: mp_abort
      implicit none
      integer :: i, iset, j
      real :: thetai, bmagi

      ! Allocate
      allocate (thetasetset(0))
      allocate (ibmagsetset(0))
      allocate (icollsetset(0))

      nsetset = 0
      starting_points: do iset = 1, nstart
         thetai = thetastart(iset)
         bmagi = bmagstart(iset)
         
! 1 For extrema, check previous sets to attach to.
         if (essentialstart(iset)) then
            do i = 1, iset-1
               if (abs(bmagstart(i)-bmagi) <= epsknob) then
! 1.1.1 If the extremum does belong in this set, eliminate points
!       near the extremum from this set.
                  where (ibmagsetset(1:nsetset) == i &
                       .and. abs(thetasetset(1:nsetset)-thetai) <= epsknob)
                     ibmagsetset(1:nsetset) = 0
                  end where
! 1.1.2 Attach the extremum to this set.
                  call add_setset (thetai, i, 112)
                  bmagstart(iset) = 0.0
                  cycle starting_points
               end if
            end do
         end if

! 2 Start a new set.
         call add_setset (thetai, iset, 2)

! 3 Check each original grid interval for matching bmag.
         grid_interval: do i = 1, nbmag-1
! 3.0.1 Stoopid problems near -pi.
! 'original' coding:
!            if (iset == 1 .and. i == 1) cycle grid_interval
! modification 5.16.02:
            if (iset == 1 .and. (i == 1 .or. i == nbmag-1)) cycle grid_interval
            if (bmagin(i) > bmagi .and. thetain(i) /= thetai) then
! 3.1 Consider when the left grid point is greater than the target bmag.
!    Then, there are three cases in which there are matching points.
!    (1) The right grid point is equal to the target bmag, and the slope
!        at the right point is positive.
!    (2) The right gridpoint is less than the target bmag.
!    (3) The right grid point is greater than the target bmag,
!        and the interval is concave down, and the minimum in
!        this interval, which is guaranteed to exist, is less than
!        the target bmag.
               if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                    .and. bprime(i+1) > 0.0) then
! 3.1.1 Consider when the right grid point is equal to the target bmag.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,311)
                  cycle grid_interval
               end if
               if (bmagin(i+1) < bmagi) then
! 3.1.2 Consider when the right grid point is less than the target bmag.
! 3.1.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.1.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,3122)
                  cycle grid_interval
               end if
! 3.1.3 Check if the grid interval is concave down.
! 3.1.3.1 If not, skip to the next interval.
               if (bprime(i) >= 0.0 .or. bprime(i+1) <= 0.0) cycle grid_interval
! 3.1.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.1.3.2.2 Otherwise, the other target point is right of the starting
!           point is the slope is negative, and left if positive.
                  if (bprime(i) < 0.0) then
                     call add_setset_root (thetai,thetain(i+1),bmagi,iset,313221)
                  else
                     call add_setset_root (thetai,thetain(i),bmagi,iset,313222)
                  end if
                  cycle grid_interval
               end if
! 3.1.3.3 If the minimum within this interval is less than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the minimum.
! 3.1.3.3.1 If this interval is on the edge, there will not be any minimum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.1.3.3.2 Find the minimum in this interval.
               if (bextr(i) == 0.0) then
                  print *, "gridgen4.f90:gg4collect:3.1.3.3.2:"," missing extremum"
                  print *, "iset,i:",iset,i
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  call mp_abort("gridgen4.f90:gg4collect:3.1.3.3.2: missing extremum")
               end if
! 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
!             the next interval.
               if (bextr(i) > bmagi) cycle grid_interval
! 3.1.3.3.2.2 Collect the point left of the minimum.
               call add_setset_root (thetaextr(i),thetain(i),bmagi,iset,313322)
! 3.1.3.3.2.3 Collect the point right of the minimum.
               call add_setset_root (thetaextr(i),thetain(i+1),bmagi,iset,313323)
               cycle grid_interval
            else if (bmagin(i) < bmagi .and. thetain(i) /= thetai) then
! 3.2 Consider then the left grid point is less than the target bmag.
!     Then, there are three cases in which there are matching points.
!     (1) The right grid point is equal to the target bmag, and the
!         slope at the right point is negative.
!     (2) The right grid point is greater than the target bmag.
!     (3) The right grid point is less than the target bmag,
!         and the interval is concave up, and the maximum in
!         this interval, which is guaranteed to exist, is greater
!         than the target bmag.
! 3.2.1 Consider when the right grid point is equal to the target bmag.
               if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                    .and. bprime(i+1) < -bpknob) then
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,321)
                  cycle grid_interval
               end if
               if (bmagin(i+1) > bmagi) then
! 3.2.2 Consider when the right grid point is greater than the target bmag.
! 3.2.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.2.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,3222)
                  cycle grid_interval
               end if
! 3.2.3 Check if the grid interval is concave up.
! 3.2.3.1 If not, skip to the next interval.
               if (bprime(i) <= 0.0 .or. bprime(i+1) >= 0.0) cycle grid_interval
! 3.2.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.2.3.2.2 Otherwise, the other target point is right of the starting
!           point if the slope is positive, and left if negative.
                  if (bprime(i) > 0.0) then
                     call add_setset_root (thetain(i+1),thetai,bmagi,iset,323221)
                  else
                     call add_setset_root (thetain(i),thetai,bmagi,iset,323222)
                  end if
                  cycle grid_interval
               end if
! 3.2.3.3 If the maximum within this interval is greater than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the maximum.
! 3.2.3.3.1 If this interval is on the edge, there will not be any maximum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.2.3.3.2 Find the maximum in this interval.
               if (bextr(i) == 0.0) then
                  print *, "gridgen4.f90:gg4collect:3.2.3.3.2:"," missing extremum"
                  print *, "iset,i:",iset,i
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  call mp_abort("gridgen4.f90:gg4collect:3.2.3.3.2: missing extremum")
               end if
! 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
!             the next interval.
               if (bextr(i) <= bmagi) cycle grid_interval
! 3.2.3.3.2.2 Collect the point left of the maximum.
               call add_setset_root (thetain(i),thetaextr(i),bmagi,iset,323322)
! 3.2.3.3.2.3 Collect the point right of the maximum.
               call add_setset_root (thetain(i+1),thetaextr(i),bmagi,iset,323323)
               cycle grid_interval
            else if (bmagin(i) == bmagi .or. thetain(i) == thetai) then
! 3.3 Consider when then left grid point is equal to the target bmag.
! 3.3.1 Add the point if it is not the starting grid point.
               if (thetai /= thetain(i)) then
                  call add_setset (thetain(i), iset, 331)
               end if
! 3.3.2 Check if there is another matching target bmag in the interval.
               if (bprime(i) > 0.0 .and. bmagin(i+1) < bmagi &
                    .and. i /= 1 .and. i /= nbmag-1) then
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,3321)
                  cycle grid_interval
               else if (bprime(i) < 0.0 .and. bmagin(i+1) > bmagi &
                    .and. i /= 1 .and. i /= nbmag-1) then
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,3322)
                  cycle grid_interval
               end if
            end if

         end do grid_interval
      end do starting_points

! 4 Compatibility with the old gg4collect.
      allocate (bmagset(nstart))
      allocate (essentialset(nstart))

      nset = nstart
      bmagset = bmagstart(1:nstart)
      essentialset = essentialstart(1:nstart)

    end subroutine gg4collect

    !> FIXME : Add documentation    
    subroutine add_setset (theta, ibmag, icoll)
      implicit none
      real, intent (in) :: theta
      integer, intent (in) :: ibmag, icoll
      thetasetset = [thetasetset, theta]
      ibmagsetset = [ibmagsetset, ibmag]
      icollsetset = [icollsetset, icoll]
      nsetset = size(thetasetset)
      if (any(ibmagsetset(1:nsetset-1) == ibmag .and. &
           abs(thetasetset(1:nsetset-1)-theta) <= epsknob)) &
           then
         ibmagsetset(nsetset) = 0
      end if
    end subroutine add_setset

    !> FIXME : Add documentation    
    subroutine add_setset_root (thetal, thetag, bmagi, ibmag, icoll)
      implicit none
      real, intent (in) :: thetal, thetag, bmagi
      integer, intent (in) :: ibmag, icoll
      real :: thl, thg, theta0, bmag0
      integer :: i
      !Bisect between thl and thg to find theta such that B(theta0) = bmagi
      thl = thetal ; thg = thetag
      do i = 1, maxiter
         theta0 = 0.5 * (thl + thg)
         bmag0 = bmagint(theta0)
         if (bmag0 > bmagi) then
            thg = theta0
         else if (bmag0 < bmagi) then
            thl = theta0
         else
            exit
         end if
      end do
      call add_setset (theta0, ibmag, icoll)
    end subroutine add_setset_root

    !> FIXME : Add documentation    
    subroutine gg4sort
      implicit none
      allocate (ithetasort(nsetset), ibmagsort(nset))
      call heapsort (nsetset, thetasetset(1:nsetset), ithetasort)
      call heapsort (nset, bmagset(1:nset), ibmagsort)
    end subroutine gg4sort

    !MRH
    !bset(n) = bmagset(ibmagsort(i))
    !The above line comes from GG4results
    ! For that line to make sense index here must be:
    ! A map which takes the new position of an element in the sorted list
    ! Returns the index of the old position of the same element in the unsorted list
    ! MRH

    !> FIXME : Add documentation    
    subroutine heapsort (n, v, index)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: v
      integer, dimension (n), intent (out) :: index
      integer :: i, j, l, ir, ira
      
      index = (/ (i, i=1,n) /)
      !if (debug) write(6,*) "heapsort index", index !MRH
      l = n/2+1
      ir = n
      
      do
         if (l > 1) then
            l = l - 1
            ira = index(l)
         else
            ira = index(ir)
            index(ir) = index(1)
            ir = ir - 1
            if (ir == 1) then
               index(1) = ira
               return
            end if
         end if
         i = l
         j = l + l
         do while (j <= ir)
            if (j < ir) then
               if (v(index(j)) < v(index(j+1))-precision_bound) j = j + 1 ! precision_bound =1.0e3*epsilon(0.0)
            ! MRH the epsilon number in the inequality has been inserted so that the logical operator gives a 
            ! MRH well-defined result across machines when v(index(j))==v(index(j+1))
               !if (v(index(j+1)) - v(index(j)) >= 0.0) j = j + 1 !MRH
            end if
            if (v(ira) < v(index(j))-precision_bound) then ! precision_bound = 1.0e3*epsilon(0.0)
            ! MRH the epsilon number in the inequality has been inserted so that the logical operator gives a 
            ! MRH well-defined result across machines when v(ira)==v(index(j))
            !if (v(index(j))-v(ira) >=0.0) then !MRH
               index(i) = index(j)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
         end do
         index(i) = ira
      end do
    end subroutine heapsort

    !> FIXME : Add documentation    
    subroutine gg4metrics
      implicit none
      integer :: i
      
      allocate (thetares(nset), alambdares(nset))
      
      do i = 1, nset
         if (bmagset(i) /= 0.0) then
            call get_thetares (i, thetares(i))
            call get_lambdares (i, alambdares(i))
         else
            thetares(i) = 0.0
            alambdares(i) = 0.0
         end if
      end do
      
    end subroutine gg4metrics

    !> FIXME : Add documentation    
    subroutine get_thetares (iset, thetares)
      implicit none
      integer, intent (in) :: iset
      real, intent (out) :: thetares
      integer :: i, ileft, iright, npts
      real :: res, dthetal, dthetar, dtheta_tmp !MRH

! 1 Return large value for essential sets with odd number of points or
!   for set containing -pi.
      if (essentialset(iset)) then
! 1.1 Set containing -pi.
         if (iset == 1) then
            thetares = 2e20
            return
         end if
! 1.2 Sets with odd number of points.
! 1.2.1 Count points.
!       npts = count(ibmagsetset == iset)
         npts = count(ibmagsetset(1:nsetset) == iset)
         if (mod(npts,2) == 1) then
            thetares = 1e20
            return
         end if
         res = extrknob/real(npts*npts)
      else
         res = 0.0
      end if

! 2 Look through all points for points in the set.
      npts = 0
      points_in_set: do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) /= iset) cycle
!MRH The above line forces i to be such that we are looking at 
!MRH an element in the unsorted lists with the iset label we are considering (labels pairs of bmag values)

! 2.1 Look for the point to the left.
! MRH the purpose of the below loops seems to be
! MRH to find the index "ileft" of the closest location to the left with a nonzero magnetic field "bmagset" (In the sorted thetasetset grid)
! MRH to find the index "iright" of the closest location to the right with a nonzero magnetic field "bmagset" (In the sorted thetasetset grid)
! MRH to see which theta values correspond to ibmagsetset values, (and hence which bmagset values)
! MRH look at the "# gg4collect - thetasetset ibmagsetset icollsetset" debug output
         ileft = i
         do
            ileft = ileft - 1
            if (ileft < 1) cycle points_in_set
            if (ibmagsetset(ithetasort(ileft)) == 0) cycle
            if (bmagset(ibmagsetset(ithetasort(ileft))) /= 0.0) exit
         end do
! 2.2 Look for the point to the right.
         iright = i
         do
            iright = iright + 1
            if (iright > nsetset) cycle points_in_set
            if (ibmagsetset(ithetasort(iright)) == 0) cycle
            if (bmagset(ibmagsetset(ithetasort(iright))) /= 0.0) exit
         end do
! 2.3 Add contribution from this interval.
! MRH the following seems to compute some kind of weighted resolution in theta
! MRH all points with the same Theta value (or possibly magnetic field strength value) contribute to this resolution

      !MRH The point of the following modifications to avoid numbers
      !MRH of the size of machine precision appearing in thetares
      !MRH which is used to choose which theta points to delete
      !MRH machine precision size numbers differ from machine to machine 
      !MRH and so lead to different data grids on different machines 
         dtheta_tmp = abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(ileft)))
         if (dtheta_tmp < precision_bound) then !1.0e3*epsilon(0.0)
            dthetal =0.0
         else 
            dthetal=dtheta_tmp
         endif
         
         dtheta_tmp = abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(iright)))
         if (dtheta_tmp < precision_bound) then !1.0e3*epsilon(0.0)
            dthetar=0.0
         else
            dthetar=dtheta_tmp
         endif
        !MRH
         !dthetal=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(ileft)))
         !dthetar=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(iright)))
         res = res + tfunc(thetasetset(ithetasort(i)))*dthetal*dthetar &
              /(dthetal+dthetar+1e-20)
         npts = npts + 1
      end do points_in_set
      if (npts > 0) then
         thetares = res/real(npts)
      else
         thetares = res
      end if
    end subroutine get_thetares
    
    !> FIXME : Add documentation
    subroutine get_lambdares (iset, alambdares)
      implicit none
      integer, intent (in) :: iset
      real, intent (out) :: alambdares
      integer :: i, iplus, iminus, npts
      real :: res, al, alplus, alminus, dalplus, dalminus, dal_tmp !MRH
      
      al = 1.0/bmagset(iset)
      do i = 1, nset
! 1 Look for target lambda
         if (ibmagsort(i) == iset) then
! 2 Look for bordering lambdas.
            alplus = 0.0
            do iplus = i-1, 1, -1
               if (ibmagsort(iplus) == 0) cycle
               if (bmagset(ibmagsort(iplus)) /= 0.0) then
                  alplus = 1.0/bmagset(ibmagsort(iplus))
                  exit
               end if
            end do
            alminus = 0.0
            do iminus = i+1, nset
               if (ibmagsort(iminus) == 0) cycle
               if (bmagset(ibmagsort(iminus)) /= 0.0) then
                  alminus = 1.0/bmagset(ibmagsort(iminus))
                  exit
               end if
            end do
            exit
         end if
      end do

! 3 Add up contributions to the result.
      res = 0.0
      npts = 0
      do i = 1, nbmag
      !MRH The point of the following modifications to avoid numbers
      !MRH of the size of machine precision appearing in lambdares
      !MRH which is used to choose which theta points to delete
      !MRH machine precision size numbers differ from machine to machine 
      !MRH and so lead to different data grids on different machines 
         dal_tmp = abs(sqrt(max(1.0-alplus*bmagin(i),0.0))-sqrt(max(1.0-al*bmagin(i),0.0))) 
         if (dal_tmp< precision_bound) then !1.0e5*epsilon(0.0)
            dalplus = 0.0 
         else
            dalplus=dal_tmp
         endif
         dal_tmp = abs(sqrt(max(1.0-alminus*bmagin(i),0.0))-sqrt(max(1.0-al*bmagin(i),0.0)))
         if (dal_tmp < precision_bound) then !1.0e5*epsilon(0.0)
            dalminus = 0.0
         else
            dalminus = dal_tmp 
         endif
         !MRH
         !dalplus = abs(sqrt(max(1.0-alplus*bmagin(i),0.0)) &
         !     -sqrt(max(1.0-al*bmagin(i),0.0)))
         !dalminus = abs(sqrt(max(1.0-alminus*bmagin(i),0.0)) &
         !     -sqrt(max(1.0-al*bmagin(i),0.0)))
         if (dalplus+dalminus /= 0.0) then
            npts = npts + 1
            res = res + dalplus*dalminus/(dalplus+dalminus+1e-20)
         end if
      end do
      if (npts /= 0) then
         alambdares = res/real(npts)
      else
         alambdares = res
      end if
    end subroutine get_lambdares

    !> FIXME : Add documentation    
    subroutine gg4remove
      use mp, only: mp_abort
      implicit none
      integer :: idel, i
      integer :: ntheta_left, nlambda_left
      integer, dimension (nset) :: work
      real, dimension (nset) :: resolution_metric ! MRH to address machine dependence
      
      nlambda_left = nset
      ntheta_left = nsetset
      do i = 1, nset
         if (.not. any(ibmagsetset(1:nsetset) == i)) then
            nlambda_left = nlambda_left - 1
         end if
      end do
! 1 Find the set with the minimum resolution metric. !MRH Actually we are looking for a set of maximal resolution_metric values
      do while (ntheta_left > ntheta .or. nlambda_left > nbset)
         
         ! ! ! resolution_metric = thetares+alknob*alambdares 
        ! ! ! ! MRH since we want to delete the iset value for which the resolution_metric is minimum
        ! ! ! ! MRH we do not change the result by subtracting the minimum value of the resolution_metric array from
        ! ! ! ! MRH from itself for all array values. We need to do this because resolution_metric is filled with either zeros
        ! ! ! ! MRH or pairs of values that identical to machine precision
        ! ! ! ! MRH the Fortran intrinsic minloc should pick the 1st minimum value in array element order when there are multiple minimum values
        ! ! ! ! MRH however depending on the compiler it seems to have trouble doing this for nonzero values of the minimum.
        ! ! ! ! MRH This is likely because for finite R,L L==R (when displayed) L-R ~ epsilon(0.0) the machine precision 
        ! ! ! ! MRH machine precision size numbers differ from machine to machine 
        ! ! ! ! MRH and so this leads to different theta grids on different machines 
        ! ! ! ! MRH ==> therefore we insure here that the minimum value is always exactly 0.
        
         ! ! ! ! MRH the following line subtracts the minimum value 
         ! ! ! ! MRH for which the mask bmagset(1:nset) /= 0.0 is true (i.e. we have not already deleted this value of iset)
         ! ! ! resolution_metric = resolution_metric - minval(resolution_metric,bmagset(1:nset) /= 0.0)
         ! ! ! ! MRH we then ensure that for every value which is minimum, the entry is exactly 0
         ! ! ! ! MRH this is so that minloc selects the 1st minimum value in array element order
         ! ! ! do i =1,nset
            ! ! ! if(resolution_metric(i)< precision_bound ) then !1.0e5*epsilon(0.0)
                ! ! ! resolution_metric(i)=0.0
            ! ! ! endif
         ! ! ! enddo
         ! MRH this line then gives us the next element to delete
         
         ! ! !work(1:1) = minloc(resolution_metric(1:nset),bmagset(1:nset) /= 0.0)
         
         ! MRH I have improved the above fix by making it more transparent-
        ! MRH I have written my own version of the intrinsic function minloc
        ! MRH minimum_index has all the documented features of minloc for a 1D array
        ! MRH with the added feature of the final machine_error variable
        ! MRH this variable gives a definite truth value to (array(i) < array (j) - machine_error) 
        ! MRH in the case that array(i) is identical to array(j).
        ! MRH This ensures the minimum_index returns the 1st minimum index in array order
        ! MRH this turns out to be necessary because in the 1st iteration in this loopafter all of the duplicate 
        ! MRH theta grid pairs have been deleted the resolution metric is calculated for the 
        ! MRH remaining theta grid which is made of almost equally spaced points, and as a
        ! MRH consequence the resolution_metric contains values which are almost
        ! MRH identical-differing only by machine precision- 
        ! MRH this leads to a random choice in which theta grid points to keep
        ! MRH in order to make this choice consistently across machines
        ! MRH we always remove the theta pair corresponding to the 1st value (in array
        ! MRH element order) of resolution_metric when the values of resolution_metric
        ! MRH are identical up to machine precision
         resolution_metric = thetares+alknob*alambdares
         work(1:1) = minimum_index (nset, resolution_metric, bmagset(1:nset) /= 0.0, precision_bound)
         !work(1:1) = minloc(thetares(1:nset)+alknob*alambdares(1:nset), &
         !     bmagset(1:nset) /= 0.0)
         idel = work(1)
         if (debug) then
            write(6,*) "gridgen4_2: call gg4debug"
            ! MRH if one follows this debug output one sees that:
! MRH whilst there are duplicate theta values remaining to be considered in the thetasetset grid
! MRH the values of the resolution_metric are initially all 0 (whilst every grid point has a duplicate)
! MRH and as the duplicates are deleted (every element of the list that is deleted is a duplicate until we reach midway)
! MRH more values of resolution_metric steadily become non 0.
! MRH finally when all the duplicates are deleted the values of the resolution_metric almost identical and differ
! MRH only at machine precision
! MRH this is because after all of the duplicate pairs have been deleted but before any more pairs have been deleted
! MRH the resolution metric is calculated based on the difference between almost equally spaced theta grid points
! MRH (see "# gg4collect - thetasetset ibmagsetset icollsetset" in .gridgen.200 to see the pairs of theta grid points
! MRH and their magnetic field strength label ibmagsetset - idel is always an element of ibmagsetset
! MRH - in this way one can keep track of which pairs are being deleted)
! MRH the bug that I have fixed was that depending on the machine, the resolution_metric at this stage in the computation
! MRH had entries which differed at the machine precision level, and as a consequence the next theta grid pair to be
! MRH deleted was effectively chosen at random, which biased the rest of the procedure so that the theta grid differed
! MRH from machine to machine.
! MRH By making the comparison in minimum_index with a bias, we avoid a machine dependent theta grid
            if (debug_output) write (unit=debug_unit, fmt=*) "gg4remove idel ",idel
            call gg4debug4 (nset, thetares, alambdares,bmagset,resolution_metric, "gg4metrics thetares alambdares bmagset(mask F if 0.0) resolution_metric")
         end if
! 2 Delete the set just found.
! MRH the subroutine delete_set deletes the iset value that the above lines just found
! MRH it does this by setting bmagset (idel) =0.0 so that when gg4results computes the final field and theta grid
! MRH this value of iset is ignored. Finally the metrics are computed taking into account that the point has been deleted
         if (idel == 0) then
            call mp_abort("gridgen4.f:gg4remove:2: This cannot happen.", .true.)
         end if
         call delete_set (idel, work, ntheta_left)
         nlambda_left = nlambda_left - 1
      end do
    end subroutine gg4remove

    !> FIXME : Add documentation
    pure integer function minimum_index (n, array, mask, machine_error)
      ! MRH the function minimum_index has been written to have the same features as minloc
      ! MRH n is the length of the array and mask
      ! MRH array is the array from which we want the index of the minimum value in the array
      ! MRH mask is an array of logicals of the same dimension as array.
      ! MRH we consider values in array for which mask is true. If all values of mask are false minimum_index is 0
      ! MRH machine_error is a control parameter in case we are comparing numbers which are identical in size
      ! MRH this function returns the 1st minimum in array element order if there are multiple minimum of the same size

      ! MRH minimum_index has all the documented features of minloc for a 1D array
      ! MRH with the added feature of the final machine_error variable
      ! MRH this variable gives a definite truth value to (array(i) < array (j) - machine_error)
      ! MRH in the case that array(i) is identical to array(j).
      ! MRH This ensures the minimum_index returns the 1st minimum index in array order
      implicit none
        real, dimension (:), intent (in) :: array
        logical, dimension (:), intent (in) :: mask
        integer, intent (in) :: n
        real, intent (in) :: machine_error
        integer :: i
        real :: minimum_value

        minimum_index = 0
        if (all (.not. mask)) return ! MRH if all of mask is false there are no values to consider in the array and so we return with 0 as the minimum index
        ! MRH get a starting value for the minimum-above we confirmed above there is at least one value in the array
        do i=1,n
            if (mask(i)) then
                minimum_value = array(i)
                minimum_index = i
                exit
            end if
        end do
        ! MRH now go through all the values in the array to find the minimum

        do i=1,n
            if (mask (i) .and. (array (i) < minimum_value- machine_error )) then !
                minimum_value = array (i)
                minimum_index =i
            end if
        end do
    end function minimum_index

    !> FIXME : Add documentation    
    subroutine delete_set (idel, work, ntheta_left)
      implicit none
      integer, intent (in) :: idel
      integer, intent (in out), dimension (nset) :: work
      integer, intent (out) :: ntheta_left
      integer :: i, j
      
      work = 0
! 1 Mark neighboring lambda sets to be recalculated.
      do i = 1, nset
         if (ibmagsort(i) == idel) then
            do j = i-1, 1, -1
               if (bmagset(ibmagsort(j)) /= 0.0) then
                  work(ibmagsort(j)) = ibmagsort(j)
                  exit
               end if
            end do
            do j = i+1, nset
               if (bmagset(ibmagsort(j)) /= 0.0) then
                  work(ibmagsort(j)) = ibmagsort(j)
                  exit
               end if
            end do
            exit
         end if
      end do

! 2 Mark lambda sets with neighboring theta points to have their resolution
!   metrics recalculated, counting the remaining theta points.
      ntheta_left = 0
      do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) == 0) cycle
         if (bmagset(ibmagsetset(ithetasort(i))) == 0.0) cycle
         if (ibmagsetset(ithetasort(i)) /= idel) then
            ntheta_left = ntheta_left + 1
            cycle
         end if
! 2.1 Found point to be deleted.
! 2.1.1 Mark set of neighboring point to the left to be recalculated.
         do j = i-1, 1, -1
            if (ibmagsetset(ithetasort(j)) == idel) cycle
            if (ibmagsetset(ithetasort(j)) == 0) cycle
            if (bmagset(ibmagsetset(ithetasort(j))) == 0.0) cycle
            work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
            exit
         end do
! 2.1.2 Mark set of neighboring point to the right to be recalculated.
         do j = i+1, nsetset
            if (ibmagsetset(ithetasort(j)) == idel) cycle
            if (ibmagsetset(ithetasort(j)) == 0) cycle
            if (bmagset(ibmagsetset(ithetasort(j))) == 0.0) cycle
            work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
            exit
         end do
      end do

! 3 Delete this set.
      bmagset(idel) = 0.0

! 4 Recalculate resolution metric for affected sets.
      do i = 1, nset
         if (work(i) /= 0) then
            call get_thetares (work(i), thetares(work(i)))
            call get_lambdares (work(i), alambdares(work(i)))
         end if
      end do
    end subroutine delete_set
    
    !> FIXME : Add documentation
    subroutine gg4results (iperiod)
      implicit none
      integer, intent (in) :: iperiod
      integer :: i, n
      
      n = 0
      do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) == 0) cycle
         if (bmagset(ibmagsetset(ithetasort(i))) == 0.0) cycle
         n = n + 1
         thetagrid(n) = thetasetset(ithetasort(i))
         bmaggrid(n) = bmagset(ibmagsetset(ithetasort(i)))
      end do
      
! Check point at +pi*iperiod
      if (bmaggrid(n) /= bmaggrid(1)) then
         n = n + 1
         thetagrid(n) = pi*iperiod
         bmaggrid(n) = bmaggrid(1)
      end if
      nthetaout = n

      n = 0
      do i = nset, 1, -1
         if (bmagset(ibmagsort(i)) == 0.0) cycle
         n = n + 1
         bset(n) = bmagset(ibmagsort(i))
      end do
      nlambdaout = n
    end subroutine gg4results

    !> FIXME : Add documentation    
    subroutine gg4debug (n, x, y, label)
      implicit none
      integer, intent (in) :: n
      real, dimension (:), intent (in) :: x, y
      character(*), intent (in) :: label
      integer :: i
      if (.not. debug_output) return
      write (unit = debug_unit, fmt = "('#',a)") label
      do i = 1, n
         write (unit = debug_unit, fmt = "(i5,1x,2(g25.18,1x))") i, x(i), y(i)
      end do
    end subroutine gg4debug

    !> FIXME : Add documentation    
    subroutine gg4debugi (n, x, i1, i2, label)
      implicit none
      integer, intent (in) :: n
      real, dimension (:), intent (in) :: x
      integer, dimension (:), intent (in) :: i1, i2
      character(*), intent (in) :: label
      integer :: i
      if (.not. debug_output) return
      write (unit = debug_unit, fmt = "('#',a)") label
      do i = 1, n
         write (unit = debug_unit, fmt = "(i5,1x,g25.18,1x,i4,1x,i10)") i, x(i), i1(i), i2(i)
      end do
    end subroutine gg4debugi

    !> FIXME : Add documentation    
    subroutine gg4debugl (n, x, label)
      implicit none
      integer, intent (in) :: n
      logical, dimension (:), intent (in) :: x
      character(*), intent (in) :: label
      integer :: i
      if (.not. debug_output) return
      write (unit = debug_unit, fmt = "('#',a)") label
      do i = 1, n
         write (unit = debug_unit, fmt = *) i, x(i)
      end do
    end subroutine gg4debugl

    !> FIXME : Add documentation    
    subroutine gg4debug4 (n, x, y, p, q, label)
      implicit none
      integer, intent (in) :: n
      real, dimension (:), intent (in) :: x, y, p, q
      character(*), intent (in) :: label
      integer :: i
      if (.not. verbose_debug_output) return
      write (unit = debug_unit, fmt = "('#',a)") label
      do i = 1, n
         write (unit = debug_unit, fmt = "(i5,1x,4(g22.15,1x))") i, x(i), y(i), p(i), q(i)
      end do
    end subroutine gg4debug4

    !> FIXME : Add documentation    
    real function bmagint (theta)
      use splines, only: periodic_splint
      implicit none
      real, intent (in) :: theta
      bmagint = periodic_splint(theta, bmagspl)
    end function bmagint

    !> FIXME : Add documentation    
    real function bmagp (theta)
      use splines, only: periodic_dsplint
      implicit none
      real, intent (in) :: theta
      bmagp = periodic_dsplint(theta, bmagspl)
    end function bmagp

    !> FIXME : Add documentation    
    elemental real function tfunc (theta)
      implicit none
      real, intent(in) :: theta
      ! The brackets multiplied by deltaw looks like a smooth function
      ! with peaks located at theta = +/- thetamax and peaks widths
      ! controlled by widthw.  As the return of tfunc is used to
      ! determine this theta grid point's contribution to the
      ! resolution metric it seems that the intent of this form is to
      ! allow the user to bias the resolution metric to be larger at
      ! theta = +/- thetamax and hence to preferentially keep points
      ! in this area.
      tfunc = 1.0 + deltaw * ( &
           1.0 / (1.0 + (theta - thetamax)**2 / widthw**2) + &
           1.0 / (1.0 + (theta + thetamax)**2 / widthw**2))
    end function tfunc
    
  end subroutine gridgen4_2
end module gridgen4mod
