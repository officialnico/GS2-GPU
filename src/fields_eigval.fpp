#ifdef FIELDS_EIG
#include "unused_dummy.inc"
module fields_eigval
  use petscvec
  use petscmat
  use slepceps

  implicit none

!Allows us to use SLEPC_VERSION_MAJOR etc. to automatically
!disable unsupported features
#include <slepcversion.h>
#include <petscversion.h>

#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>
#elif PETSC_VERSION_LT(3,8,0)
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#endif

#if SLEPC_VERSION_LT(3,6,0)
#include <finclude/slepcepsdef.h>
#elif SLEPC_VERSION_LT(3,8,0)
#include <slepc/finclude/slepcepsdef.h>
#else
#include <slepc/finclude/slepceps.h>
#endif

  private

  public :: eigval_functional
  public :: init_eigval, finish_eigval, BasicSolve

  logical :: initialized = .false.

  !Number of GS2 advance steps per SLEPc call to advance_eigen
  !May be useful when looking at marginal modes etc.
  integer :: nadv = 500

contains

  !> Returns true if GS2 was compiled with WITH_EIG defined
  function eigval_functional()
    logical :: eigval_functional
    eigval_functional=.true.
  end function eigval_functional

  !> Initialise eigenvalue problem
  subroutine init_eigval!(eigval_config_in)
    use mp, only: mp_comm
    implicit none
    !type(eigval_config_type), intent(in), optional :: eigval_config_in
    PetscErrorCode :: ierr

    if (initialized) return
    initialized = .true.

    !Get the input parameters
    !call read_parameters(eigval_config_in)

    !Set PETSC_COMM_WORLD to be mp_comm so we can use whatever sub-comm
    !needed for list mode to work
    PETSC_COMM_WORLD=mp_comm

    !Initialise slepc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  end subroutine init_eigval

  !> Clean up eigenvalue problem
  subroutine finish_eigval
    implicit none
    PetscErrorCode :: ierr
    initialized = .false.
    !Finish up slepc
    call SlepcFinalize(ierr)
  end subroutine finish_eigval

  !> Main routine that does the setup, starts the solver and
  !! coordinates reporting of the results
  subroutine BasicSolve
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: has_phi, has_apar, has_bpar
    use mp, only: proc0
    implicit none
    EPS :: my_solver
    Mat :: my_operator
    PetscInt :: d1_size, d2_size, d3_size, d4_size
    PetscInt :: n_loc, n_glob
    PetscErrorCode :: ierr
    integer :: nfield

    nfield = count([has_phi, has_apar, has_bpar])

    !Define dimensions for problem
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=ntheta0    !Size of theta0 grid -- has to be 1
    d3_size=naky       !Size of ky grid -- has to be
    d4_size=nfield

    !How big is the operator
    n_loc  = d1_size*d2_size*d3_size*d4_size
    n_glob = n_loc
    if (.not. proc0) n_loc = 0

    !Initialise the matrix operator

    !Make a shell matrix operator (AdvMat)
    call MatCreateShell(PETSC_COMM_WORLD, n_loc, n_loc,&
         n_glob, n_glob, PETSC_NULL_INTEGER, my_operator, ierr)

    !Set the shell MATOP_MULT operation, i.e. which routine returns the result
    !of a AdvMat.x
    call MatShellSetOperation(my_operator, MATOP_MULT, advance_eigen, ierr)

    !Create the solver

    !Now create the eps solver
    call EPSCreate(PETSC_COMM_WORLD, my_solver, ierr)

    !Set the matrix operator we want to find eigenvalues of
#if PETSC_VERSION_LT(3,8,0)
    call EPSSetOperators(my_solver, my_operator, PETSC_NULL_OBJECT, ierr)
#else
    call EPSSetOperators(my_solver, my_operator, PETSC_NULL_MAT, ierr)
#endif

    !Set the type of eigenproblem -->Always non-hermittian for us
    call EPSSetProblemType(my_solver,EPS_NHEP,ierr)

    !Now setup from options
    call EPSSetFromOptions(my_solver,ierr)

    !Now we're ready to solve
    call EPSSolve(my_solver,ierr)

    !Now do the reporting and diagnostic output
    call ReportResults(my_solver, my_operator)

    !Now destroy objects
    call EPSDestroy(my_solver,ierr)
    call MatDestroy(my_operator,ierr)
  end subroutine BasicSolve

  !> Call back routine to represent linear advance operator
  subroutine advance_eigen(MatOperator, VecIn, Res, ierr)
    use fields, only: set_init_fields, advance
    use dist_fn_arrays, only: gnew, g
    use fields_arrays, only: phinew
    use array_utils, only: zero_array, copy
    use gs2_layouts, only: g_lo
    implicit none
    Mat, intent(in) :: MatOperator
    Vec, intent(in out) :: VecIn, Res
    PetscErrorCode, intent(inout) :: ierr
    integer, parameter :: istep=2
    integer :: is, iglo
    UNUSED_DUMMY(MatOperator); UNUSED_DUMMY(ierr)

    !First unpack input vector into gnew
    call VecToFields(VecIn)
!    call zero_array(gnew)
!    call zero_array(g)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       gnew(:, 1, iglo) = phinew(:, 1, 1)
       gnew(:, 2, iglo) = phinew(:, 1, 1)
    end do
    call copy(gnew, g)

    !Now do a number of time steps
    do is=1,nadv
       !Note by using a fixed istep we
       !force the explicit terms to be excluded
       !except for the very first call.
       call advance(istep)
    enddo

    !Now pack fields into output
    call FieldsToVec(Res)
  end subroutine advance_eigen

  !> Unpacks a SLEPc vector into fields
  subroutine VecToFields(VecIn)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: has_phi, has_apar, has_bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use mp, only: proc0, broadcast
    implicit none
    Vec, intent(inout) :: VecIn
    PetscScalar, pointer :: VecInArr(:)
    PetscInt :: d1_size, d2_size, d3_size, d4_size
    PetscErrorCode :: ierr
    integer :: nfield, ig, ik, it, local_index

    nfield = count([has_phi, has_apar, has_bpar])

    !Define dimensions
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=ntheta0    !Size of theta0 grid -- has to be 1
    d3_size=naky       !Size of ky grid -- has to be 1
    d4_size=nfield

    !Get a pointer to the data
    call VecGetArrayReadF90(VecIn,VecInArr,ierr)

    ! Force iproc /= 0 to skip any work
    if (.not. proc0) d3_size = 0
    !Extract
    local_index = 0
    if (has_phi) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                phinew(ig, it, ik) = VecInArr(local_index)
             end do
          end do
       end do
       call broadcast(phinew)
       phi = phinew
    end if

    if (has_apar) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                aparnew(ig, it, ik) = VecInArr(local_index)
             end do
          end do
       end do
       call broadcast(aparnew)
       apar = aparnew
    end if

    if (has_bpar) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                bparnew(ig, it, ik) = VecInArr(local_index)
             end do
          end do
       end do
       call broadcast(bparnew)
       bpar = bparnew
    end if

    if (proc0) then
       if (local_index /= d1_size*d2_size*d3_size*d4_size) print*,'Inconsistent sizes V2F'
    end if

    !Get rid of pointer
    call VecRestoreArrayReadF90(VecIn,VecInArr,ierr)
  end subroutine VecToFields

  !> Packs fields into a SLEPc vector
  subroutine FieldsToVec(VecIn)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: has_phi, has_apar, has_bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0
    implicit none
    Vec, intent(in out) :: VecIn
    PetscScalar, pointer :: VecInArr(:)
    PetscInt :: d1_size, d2_size, d3_size, d4_size
    PetscErrorCode :: ierr
    integer :: nfield, ig, ik, it, local_index

    ! Only proc0 has any data
    if (.not. proc0) return

    nfield = count([has_phi, has_apar, has_bpar])

    !Define dimensions
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=ntheta0    !Size of theta0 grid -- has to be 1
    d3_size=naky       !Size of ky grid -- has to be 1
    d4_size=nfield

    !Get a pointer to the data
    call VecGetArrayF90(VecIn,VecInArr,ierr)

    !Extract
    local_index = 0
    if (has_phi) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                VecInArr(local_index) = phinew(ig, it, ik)
             end do
          end do
       end do
    end if

    if (has_apar) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                VecInArr(local_index) = aparnew(ig, it, ik)
             end do
          end do
       end do
    end if

    if (has_bpar) then
       do ik = 1, d3_size
          do it = 1, d2_size
             do ig = -ntgrid, ntgrid
                local_index = local_index + 1
                VecInArr(local_index) = bparnew(ig, it, ik)
             end do
          end do
       end do
    end if

    if (proc0) then
       if (local_index /= d1_size*d2_size*d3_size*d4_size) print*,'Inconsistent sizes F2V'
    end if

    !Get rid of pointer
    call VecRestoreArrayF90(VecIn,VecInArr,ierr)
  end subroutine FieldsToVec

  !> Routine to write results to screen and netcdf
  subroutine ReportResults(my_solver,my_operator)
    use mp, only: proc0
    use constants, only: zi
    use gs2_time, only: code_dt
    use run_parameters, only: has_phi, has_apar, has_bpar
    use fields, only: set_init_fields
    use fields_arrays, only: phi, phinew, apar, aparnew, bpar, bparnew
    use gs2_save, only: EigNetcdfID, init_eigenfunc_file, finish_eigenfunc_file
    use gs2_save, only: add_eigenpair_to_file, finish_gs2_save
    use file_utils, only: run_name
    use gs2_diagnostics, only: loop_diagnostics, nwrite, navg
    use run_parameters, only: use_old_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: run_diagnostics, gnostics
#endif
    implicit none
    EPS, intent(in) :: my_solver
    Mat, intent(in) :: my_operator
    PetscErrorCode :: ierr
    PetscInt :: iteration_count, n_converged
    Vec :: eig_vec_r, eig_vec_i
    PetscScalar :: eig_val_r, eig_val_i
    complex :: EigVal, previous_field_factor
    logical :: diagnostics_exit
    PetscInt :: ieig
    type(EigNetcdfID) :: io_ids
    integer :: original_navg

    !Find out how many iterations were performed
    call EPSGetIterationNumber(my_solver,iteration_count,ierr)

    !Find out how many converged eigenpairs where found
    call EPSGetConverged(my_solver,n_converged,ierr)
    if(proc0) write(6,'("Found ",I0," eigenvalues in ",I0," iterations.")') n_converged,iteration_count

    if(n_converged.gt.0)then
       !Initialise the eigenvector arrays
#if PETSC_VERSION_LT(3,6,0)
       call MatGetVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_r,ierr)
       call MatGetVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_i,ierr)
#elif PETSC_VERSION_LT(3,8,0)
       call MatCreateVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_r,ierr)
       call MatCreateVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_i,ierr)
#else
       call MatCreateVecs(my_operator,PETSC_NULL_VEC,eig_vec_r,ierr)
       call MatCreateVecs(my_operator,PETSC_NULL_VEC,eig_vec_i,ierr)
#endif
       if(proc0)call init_eigenfunc_file(trim(run_name)//"_eig.out.nc", io_ids)

       !Now loop over converged values
       do ieig=0,n_converged-1
          !Get eigenpair
          call EPSGetEigenpair(my_solver,ieig,eig_val_r,&
               eig_val_i,eig_vec_r,eig_vec_i,ierr)

          EigVal=cmplx(PetscRealPart(eig_val_r),PetscImaginaryPart(eig_val_r))
          print*,'Eigval --',ieig,EigVal,'|',Eig_SlepcToGs2(EigVal)
          !Convert to GS2 eigenvalue
          Eigval=Eig_SlepcToGs2(EigVal)

          !Get field eigenmodes
          call VecToFields(eig_vec_r)

          !Add to file
          if(proc0)call add_eigenpair_to_file(EigVal, has_phi, has_apar, has_bpar, io_ids)

          !/Calculate what the fields at the previous time step would
          !/look like given current state and Eigval. This is intended
          !/to allow diagnostics to recover Eigval from the usual
          !/calculation of omega
          previous_field_factor = exp(zi * Eigval * code_dt)
          if (has_phi) phi = phinew * previous_field_factor
          if (has_apar) apar = aparnew * previous_field_factor
          if (has_bpar) bpar = bparnew * previous_field_factor

          if (use_old_diagnostics) then
             original_navg = navg
             navg = 1
             call loop_diagnostics ((ieig+1) * nwrite, diagnostics_exit)
             navg = original_navg
#ifdef NEW_DIAG
          else
             original_navg = gnostics%navg
             gnostics%navg = 1
             call run_diagnostics ((ieig+1) * gnostics%nwrite, diagnostics_exit)
             gnostics%navg = original_navg
#endif
          end if
       end do

       !Close netcdf file
       if(proc0)call finish_eigenfunc_file(io_ids)

       call VecDestroy(eig_vec_r,ierr)
       call VecDestroy(eig_vec_i,ierr)
    else
       if(proc0) write(6,'("No converged eigenvalues found.")')
    endif
  end subroutine ReportResults

  !> Function to convert SLEPc eigenvalue to a GS2 one
  elemental complex function Eig_SlepcToGs2(eig)
    use gs2_time, only: code_dt
    use constants, only: zi
    implicit none
    complex, intent(in) :: eig
    Eig_SlepcToGs2=log(eig)*zi/(code_dt*nadv)
  end function Eig_SlepcToGs2

end module fields_eigval

#else

!> A stub module representing the eigenvalue solver
module fields_eigval
  implicit none
  private
  public :: eigval_functional
contains
  !> Returns true if GS2 was compiled with WITH_EIG defined -- here forced to .false.
  !! as compiled without SLEPc support
  function eigval_functional()
    logical :: eigval_functional
    eigval_functional=.false.
  end function eigval_functional
end module fields_eigval

#endif
