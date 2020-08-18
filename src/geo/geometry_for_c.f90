!> THE FOLLOWING FREE ROUTINES SHOULD BE PLACED IN A MODULE IF POSSIBLE

!> This subroutine implements a general interface to geometry
!! which can be used outside GS2; in particular, it can be
!! called from C/C++ programs. 
!!
!! Equilibrium types:
!!     1: Local (Miller)
!!     2: EFIT
!!     3: Chease
!!     4: Transp
!!
!! Some equilibria, those that are created from a grid of psi(R,Z),
!! which means EFIT and Chease, change ntheta. 
!! This should then be used as the number of 
!! gridpoints along the fieldline. 
!!
!! Massive note to self, check where iflux and local_eq are set
!! in GS2
subroutine geometry_set_inputs(equilibrium_type,&
                               !eqfile_in, &
                               irho_in, &
                               rhoc_in, &
                               bishop_in, &
                               nperiod_in, &
                               ntheta)
  use geometry, only: rhoc, nperiod, eqfile, irho, iflux, local_eq, gen_eq
  use geometry, only: ppl_eq, transp_eq, chs_eq, efit_eq, equal_arc, bishop
  use geometry, only: dp_mult, delrho, rmin, rmax, isym, in_nt, writelots, itor
  use geometry, only: verb
  use mp, only: mp_abort
  implicit none
  integer, intent(in) :: equilibrium_type, nperiod_in, irho_in, bishop_in
  real, intent(in) :: rhoc_in
  integer, intent(inout) :: ntheta
  !character(len=800), intent(in) :: eqfile_in

  write (*,*) 'Setting inputs'
  write (*,*)

  rhoc = rhoc_in
  nperiod = nperiod_in
  !eqfile = eqfile_in
  irho = irho_in

  iflux = 1 ! Numerical unless changed for miller

  local_eq = .false.
  gen_eq = .false.
  ppl_eq = .false.
  transp_eq = .false.
  chs_eq = .false.
  efit_eq = .false.

  equal_arc = .true. ! JRB - are we sure equal_arc works properly?
  bishop = bishop_in
  dp_mult = 1.0
  delrho = 1e-3
  rmin = 1e-3
  rmax = 1.0
  isym = 0
  in_nt = .false.
  writelots = .true.
  !! Advanced use only
  itor = 1

  verb = 4

  
  !ntheta_out = -1

  select case (equilibrium_type)
  case (1)  ! Miller
    local_eq = .true.
    iflux = 0
    bishop = 4
  case (2)  ! EFIT
    efit_eq = .true.
  case (3)  ! CHEASE
    chs_eq = .true.
  case default
    write(*,*) "Whatever geometry you are using hasn't "
    write(*,*) "been tested in a long time!"
    call mp_abort("Whatever geometry you are using hasn't been tested in a long time!",.true.)
  end select

end subroutine geometry_set_inputs  


!> This subroutine gets default values for the advanced switches,
!! returning a derived type (advanced_parameters_type) containing
!! all the parameters. It does not modify the geometry variables. 
subroutine geometry_get_default_advanced_parameters(advanced_parameters_out)
  use geometry, only: advanced_parameters_type
  implicit none
  type(advanced_parameters_type), intent(out) :: advanced_parameters_out
 
  advanced_parameters_out%equal_arc = .true.  ! JRB - are we sure equal_arc works properly?
  advanced_parameters_out%dp_mult = 1.0
  advanced_parameters_out%delrho = 1e-3
  advanced_parameters_out%rmin = 1e-3
  advanced_parameters_out%rmax = 1.0
  advanced_parameters_out%isym = 0
  advanced_parameters_out%in_nt = .false.
  advanced_parameters_out%writelots = .false.
  advanced_parameters_out%itor =  1
end subroutine geometry_get_default_advanced_parameters

!> Returns a derived type containing values of the advanced geometry
!! parameters. 
subroutine geometry_get_advanced_parameters(advanced_parameters_out)
  use geometry, only: equal_arc, dp_mult, delrho, rmin, rmax, isym, in_nt, writelots, itor
  use geometry, only: advanced_parameters_type, advanced_parameters
  implicit none
  type(advanced_parameters_type), intent(out) :: advanced_parameters_out

  advanced_parameters%equal_arc = equal_arc
  advanced_parameters%dp_mult = dp_mult
  advanced_parameters%delrho = delrho
  advanced_parameters%rmin = rmin
  advanced_parameters%rmax = rmax
  advanced_parameters%isym = isym
  advanced_parameters%in_nt = in_nt
  advanced_parameters%writelots = writelots
  advanced_parameters%itor = itor

  advanced_parameters_out = advanced_parameters

end subroutine geometry_get_advanced_parameters

!> Accepts a derived type containing values of the advanced geometry
!! parameters and sets the variables in the geometry module accordingly. 
subroutine geometry_set_advanced_parameters(advanced_parameters_in)
  use geometry, only: equal_arc, dp_mult, delrho, rmin, rmax, isym, in_nt, writelots, itor
  use geometry, only: advanced_parameters_type, advanced_parameters
  implicit none
  type(advanced_parameters_type), intent(in) :: advanced_parameters_in

  write(*,*) 'Setting advanced parameters'
  advanced_parameters = advanced_parameters_in

  equal_arc = advanced_parameters%equal_arc
  dp_mult = advanced_parameters%dp_mult
  delrho = advanced_parameters%delrho
  rmin = advanced_parameters%rmin
  rmax = advanced_parameters%rmax
  isym = advanced_parameters%isym
  in_nt = advanced_parameters%in_nt
  writelots = advanced_parameters%writelots
  itor = advanced_parameters%itor

  write(*,*) 'delrho was set to', delrho

end subroutine geometry_set_advanced_parameters

!> FIXME : Add documentation
subroutine geometry_get_miller_parameters(miller_parameters_out)
  use geometry, only: geoType, rmaj, r_geo, aSurf, shift, shiftVert, mMode, nMode
  use geometry, only: deltam, deltan, deltampri, deltanpri, thetam, thetan, qinp, shat
  use geometry, only: miller_parameters_type, miller_parameters
  implicit none
  type(miller_parameters_type), intent(out) :: miller_parameters_out
 
  miller_parameters%geoType = geoType
  miller_parameters%rmaj = rmaj
  miller_parameters%r_geo = r_geo
  miller_parameters%aSurf = aSurf
  miller_parameters%shift = shift
  miller_parameters%shiftVert = shiftVert
  miller_parameters%mMode = mMode
  miller_parameters%nMode = nMode
  miller_parameters%deltam = deltam
  miller_parameters%deltan = deltan
  miller_parameters%deltampri = deltampri
  miller_parameters%deltanpri = deltanpri
  miller_parameters%thetam = thetam
  miller_parameters%thetan = thetan
  miller_parameters%qinp = qinp
  miller_parameters%shat = shat
  
  miller_parameters_out = miller_parameters

end subroutine geometry_get_miller_parameters

!> FIXME : Add documentation
subroutine geometry_set_miller_parameters(miller_parameters_in)
  use geometry, only: geoType, rmaj, r_geo, aSurf, shift, shiftVert, mMode, nMode
  use geometry, only: deltam, deltan, deltampri, deltanpri, thetam, thetan, qinp, shat
  use geometry, only: miller_parameters_type, miller_parameters
  implicit none
  type(miller_parameters_type), intent(in) :: miller_parameters_in

  miller_parameters = miller_parameters_in

  geoType = miller_parameters_in%geoType
  rmaj = miller_parameters_in%rmaj
  R_geo = miller_parameters_in%R_geo
  aSurf = miller_parameters_in%aSurf
  shift = miller_parameters_in%shift
  shiftVert = miller_parameters_in%shiftVert
  mMode = miller_parameters_in%mMode
  nMode = miller_parameters_in%nMode
  deltam = miller_parameters_in%deltam
  deltan = miller_parameters_in%deltan
  deltampri = miller_parameters_in%deltampri
  deltanpri = miller_parameters_in%deltanpri
  thetam = miller_parameters_in%thetam
  thetan = miller_parameters_in%thetan
  qinp = miller_parameters_in%qinp
  shat = miller_parameters_in%shat

  write(*,*) 's_hat was set to', shat
  write(*,*) 'deltam was set to', deltam
  write(*,*) 'deltan was set to', deltan
  write(*,*) 'R_geo was set to', R_geo
  write(*,*) 'rmaj was set to', rmaj
  write(*,*) 'qinp was set to', qinp
  !write(*,*) 'ntheta was set to', ntheta

end subroutine geometry_set_miller_parameters

!> Vary the magnetic shear and pressure gradient 
!! while maintaing a solution to the Grad-Shrafranov eqn,
!! as described variously by Greene & Chance, Bishop and Miller.
subroutine geometry_vary_s_alpha(s_hat_input_in, beta_prime_input_in)
  use geometry, only: s_hat_input, beta_prime_input
  implicit none
  real, intent(in) :: s_hat_input_in, beta_prime_input_in
  s_hat_input = s_hat_input_in
  beta_prime_input = beta_prime_input_in
end subroutine geometry_vary_s_alpha

!> FIXME : Add documentation
subroutine geometry_calculate_coefficients(grid_size)
  use theta_grid_params, only: ntheta
  use geometry, only: eikcoefs, nperiod, equal_arc
  implicit none
  integer, intent(out) :: grid_size
  ! We calculate the equal_arc grids manually,
  ! and override equal_arc = .false.
  equal_arc = .false.
  call eikcoefs(ntheta)
  write (*,*) 'ntheta out is ', ntheta
  grid_size = (2*nperiod - 1)*ntheta + 1 ! = 2*ntgrid + 1
  write (*,*) 'grid_size out is ', grid_size
end subroutine geometry_calculate_coefficients

!> Get the geometric coefficients calculated by the geometry module.
subroutine geometry_get_coefficients(grid_size, coefficients_out)
  use geometry, only: grho, bmag, gradpar, cvdrift, cvdrift0, coefficients_type
  use geometry, only: gbdrift, gbdrift0, cdrift, cdrift0, gbdrift_th
  use geometry, only: cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq
  use geometry, only: jacob, Rplot, Zplot, aplot, Rprime, Zprime, aprime, Uk1, Uk2, Bpol
  use geometry, only: theta, theta_eqarc, theta_prime, theta_prime_eqarc
  use geometry, only: gradpar_prime
  implicit none
  integer, intent(in) :: grid_size
  type(coefficients_type), dimension(grid_size) :: coefficients_out
  integer ::ntgrid, i
  real, dimension(:,:),allocatable :: interp_matrix
   
  allocate(interp_matrix(grid_size,grid_size))

  !ntgrid = (2*nperiod - 1) * ntheta/2
  ntgrid = (grid_size - 1)/2

  call create_interp_matrix

  write (*,*) 'HERE'
  write (*,*) 'Grid size should be ', 2*ntgrid + 1
  do i = -ntgrid,ntgrid
     write (*,*) 'i', i
     coefficients_out(i+ntgrid+1)%theta        = theta(i)   
     coefficients_out(i+ntgrid+1)%theta_eqarc = theta_eqarc(i)
     coefficients_out(i+ntgrid+1)%theta_prime        = theta_prime(i)   
     coefficients_out(i+ntgrid+1)%theta_prime_eqarc = theta_prime_eqarc(i)
     coefficients_out(i+ntgrid+1)%grho        = grho(i)   
     coefficients_out(i+ntgrid+1)%grho_eqarc = interp(grho,i)
     coefficients_out(i+ntgrid+1)%bmag        = bmag(i)       
     coefficients_out(i+ntgrid+1)%bmag_eqarc = interp(bmag,i)
     coefficients_out(i+ntgrid+1)%gradpar     = gradpar(i)    
     coefficients_out(i+ntgrid+1)%gradpar_eqarc = interp(gradpar,i)
     coefficients_out(i+ntgrid+1)%gradpar_prime     = gradpar_prime(i)    
     coefficients_out(i+ntgrid+1)%gradpar_prime_eqarc = interp(gradpar_prime,i)
     coefficients_out(i+ntgrid+1)%cvdrift     = cvdrift(i)    
     coefficients_out(i+ntgrid+1)%cvdrift_eqarc = interp(cvdrift,i)
     coefficients_out(i+ntgrid+1)%cvdrift0    = cvdrift0(i)   
     coefficients_out(i+ntgrid+1)%cvdrift0_eqarc = interp(cvdrift0,i)
     coefficients_out(i+ntgrid+1)%gbdrift     = gbdrift(i)    
     coefficients_out(i+ntgrid+1)%gbdrift_eqarc = interp(gbdrift,i)
     coefficients_out(i+ntgrid+1)%gbdrift0    = gbdrift0(i)   
     coefficients_out(i+ntgrid+1)%gbdrift0_eqarc = interp(gbdrift0,i)
     coefficients_out(i+ntgrid+1)%cdrift     = cdrift(i)    
     coefficients_out(i+ntgrid+1)%cdrift_eqarc = interp(cdrift,i)
     coefficients_out(i+ntgrid+1)%cdrift0     = cdrift0(i)    
     coefficients_out(i+ntgrid+1)%cdrift0_eqarc = interp(cdrift0,i)
     coefficients_out(i+ntgrid+1)%gbdrift_th  = gbdrift_th(i) 
     coefficients_out(i+ntgrid+1)%gbdrift_th_eqarc = interp(gbdrift_th,i)
     coefficients_out(i+ntgrid+1)%cvdrift_th  = cvdrift_th(i) 
     coefficients_out(i+ntgrid+1)%cvdrift_th_eqarc = interp(cvdrift_th,i)
     coefficients_out(i+ntgrid+1)%gds2        = gds2(i)       
     coefficients_out(i+ntgrid+1)%gds2_eqarc = interp(gds2,i)
     coefficients_out(i+ntgrid+1)%gds21       = gds21(i)      
     coefficients_out(i+ntgrid+1)%gds21_eqarc = interp(gds21,i)
     coefficients_out(i+ntgrid+1)%gds22       = gds22(i)      
     coefficients_out(i+ntgrid+1)%gds22_eqarc = interp(gds22,i)
     coefficients_out(i+ntgrid+1)%gds23       = gds23(i)      
     coefficients_out(i+ntgrid+1)%gds23_eqarc = interp(gds23,i)
     coefficients_out(i+ntgrid+1)%gds24       = gds24(i)      
     coefficients_out(i+ntgrid+1)%gds24_eqarc = interp(gds24,i)
     coefficients_out(i+ntgrid+1)%gds24_noq   = gds24_noq(i)  
     coefficients_out(i+ntgrid+1)%gds24_noq_eqarc = interp(gds24_noq,i)
     coefficients_out(i+ntgrid+1)%jacob       = jacob(i)      
     coefficients_out(i+ntgrid+1)%jacob_eqarc = interp(jacob,i)
     coefficients_out(i+ntgrid+1)%Rplot       = Rplot(i)      
     coefficients_out(i+ntgrid+1)%Rplot_eqarc = interp(Rplot,i)
     coefficients_out(i+ntgrid+1)%Zplot       = Zplot(i)      
     coefficients_out(i+ntgrid+1)%Zplot_eqarc = interp(Zplot,i)
     coefficients_out(i+ntgrid+1)%aplot       = aplot(i)      
     coefficients_out(i+ntgrid+1)%aplot_eqarc = interp(aplot,i)
     coefficients_out(i+ntgrid+1)%Rprime      = Rprime(i)     
     coefficients_out(i+ntgrid+1)%Rprime_eqarc = interp(Rprime,i)
     coefficients_out(i+ntgrid+1)%Zprime      = Zprime(i)     
     coefficients_out(i+ntgrid+1)%Zprime_eqarc = interp(Zprime,i)
     coefficients_out(i+ntgrid+1)%aprime      = aprime(i)     
     coefficients_out(i+ntgrid+1)%aprime_eqarc = interp(aprime,i)
     coefficients_out(i+ntgrid+1)%Uk1         = Uk1(i)        
     coefficients_out(i+ntgrid+1)%Uk1_eqarc = interp(Uk1,i)
     coefficients_out(i+ntgrid+1)%Uk2         = Uk2(i)        
     coefficients_out(i+ntgrid+1)%Uk2_eqarc = interp(Uk2,i)
     coefficients_out(i+ntgrid+1)%Bpol        = Bpol(i)       
     coefficients_out(i+ntgrid+1)%Bpol_eqarc = interp(Bpol,i)
  end do

  deallocate(interp_matrix)

  write (*,*) 'Returning....'
  contains
    subroutine create_interp_matrix
      use geometry, only: theta, theta_eqarc
      use splines, only: inter_cspl
      integer :: j
      real, dimension(grid_size) :: delta_array
      do j = 1,grid_size
        delta_array = 0.
        delta_array(j) = 1.
        call inter_cspl(grid_size, theta, delta_array, &
                        grid_size, theta_eqarc, interp_matrix(:,j))
      end do

    end subroutine create_interp_matrix

    function interp(array, ig)
      real, dimension(:), intent(in) :: array
      integer, intent(in) :: ig
      real :: interp
      integer :: j
      j = ig + ntgrid + 1
      interp = dot_product(array(:),interp_matrix(j,:))
      !interp = array(j)
    end function interp

end subroutine geometry_get_coefficients

!> FIXME : Add documentation
subroutine geometry_get_constant_coefficients(constant_coefficients_out)
  use geometry, only: qsf, rmaj, shat, kxfac, aminor, drhodpsin
  use geometry, only: bi_out
  use geometry, only: constant_coefficients_type, constant_coefficients
  implicit none
  type(constant_coefficients_type), intent(out) :: constant_coefficients_out

  constant_coefficients%qsf = qsf
  constant_coefficients%rmaj = rmaj
  constant_coefficients%shat = shat
  constant_coefficients%kxfac = kxfac
  constant_coefficients%aminor = aminor
  constant_coefficients%drhodpsin = drhodpsin
  constant_coefficients%bi = bi_out
  constant_coefficients_out = constant_coefficients
  
end subroutine geometry_get_constant_coefficients
