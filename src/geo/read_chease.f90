!> CHEASE Output Reader
!!
!! A module to read in datafiles from CHEASE.
!!
!! Written by Edmund Highcock
!! edmundhighcock@sourceforge.net
!!
!! This is free software released under the MIT Licence
!!
!! Available quantities are:
!! 
!! Zero D:
!! r0exp_chease,b0exp_chease
!! One D:
!! rgeom_chease,ageom_chease,q_chease,dqdpsi_chease,
!! d2qdpsi2_chease,p_chease,dpdpsi_chease,f_chease,
!! fdfdpsi_chease,v_chease,rho_t_chease,shear_chease,
!! dsheardpsi_chease,kappa_chease,delta_lower_chease,
!! delta_upper_chease,dvdpsi_chease,dpsidrhotor_chease,
!! gdpsi_av_chease,radius_av_chease,r_av_chease,te_chease,
!! dtedpsi_chease,ne_chease,dnedpsi_chease,ti_chease,
!! dtidpsi_chease,ni_chease,dnidpsi_chease,zeff_chease,
!! signeo_chease,jbsbav_chease
!! Two D:
!! g11_chease,g12_chease,g22_chease,g33_chease,b_chease,
!! dbdpsi_chease,dbdchi_chease,dpsidr_chease,dpsidz_chease,
!! dchidr_chease,dchidz_chease,jacobian_chease,r_chease,
!! z_chease
!!
!! This module is generated automatically 
!! using 
!!    $ ruby generate_read_chease.rb
!!
!! DO NOT EDIT!
!! YOUR CHANGES WILL BE LOST

! An example chease_namelist will be on the wiki

module read_chease
  implicit none

  private

  public :: read_infile
  public :: finish

  public :: npsi_chease, nchi_chease 
  public :: psi_chease, chi_chease
  public :: r0exp_chease
  public :: b0exp_chease
  public :: rgeom_chease
  public :: ageom_chease
  public :: q_chease
  public :: dqdpsi_chease
  public :: d2qdpsi2_chease
  public :: p_chease
  public :: dpdpsi_chease
  public :: f_chease
  public :: fdfdpsi_chease
  public :: v_chease
  public :: rho_t_chease
  public :: shear_chease
  public :: dsheardpsi_chease
  public :: kappa_chease
  public :: delta_lower_chease
  public :: delta_upper_chease
  public :: dvdpsi_chease
  public :: dpsidrhotor_chease
  public :: gdpsi_av_chease
  public :: radius_av_chease
  public :: r_av_chease
  public :: te_chease
  public :: dtedpsi_chease
  public :: ne_chease
  public :: dnedpsi_chease
  public :: ti_chease
  public :: dtidpsi_chease
  public :: ni_chease
  public :: dnidpsi_chease
  public :: zeff_chease
  public :: signeo_chease
  public :: jbsbav_chease
  public :: g11_chease
  public :: g12_chease
  public :: g22_chease
  public :: g33_chease
  public :: b_chease
  public :: dbdpsi_chease
  public :: dbdchi_chease
  public :: dpsidr_chease
  public :: dpsidz_chease
  public :: dchidr_chease
  public :: dchidz_chease
  public :: jacobian_chease
  public :: r_chease
  public :: z_chease

  integer :: npsi_chease, nchi_chease 
  real, dimension (:), allocatable :: psi_chease,chi_chease  
  integer :: infile=1212

  real :: r0exp_chease
  real :: b0exp_chease
  real, dimension(:), allocatable :: rgeom_chease
  real, dimension(:), allocatable :: ageom_chease
  real, dimension(:), allocatable :: q_chease
  real, dimension(:), allocatable :: dqdpsi_chease
  real, dimension(:), allocatable :: d2qdpsi2_chease
  real, dimension(:), allocatable :: p_chease
  real, dimension(:), allocatable :: dpdpsi_chease
  real, dimension(:), allocatable :: f_chease
  real, dimension(:), allocatable :: fdfdpsi_chease
  real, dimension(:), allocatable :: v_chease
  real, dimension(:), allocatable :: rho_t_chease
  real, dimension(:), allocatable :: shear_chease
  real, dimension(:), allocatable :: dsheardpsi_chease
  real, dimension(:), allocatable :: kappa_chease
  real, dimension(:), allocatable :: delta_lower_chease
  real, dimension(:), allocatable :: delta_upper_chease
  real, dimension(:), allocatable :: dvdpsi_chease
  real, dimension(:), allocatable :: dpsidrhotor_chease
  real, dimension(:), allocatable :: gdpsi_av_chease
  real, dimension(:), allocatable :: radius_av_chease
  real, dimension(:), allocatable :: r_av_chease
  real, dimension(:), allocatable :: te_chease
  real, dimension(:), allocatable :: dtedpsi_chease
  real, dimension(:), allocatable :: ne_chease
  real, dimension(:), allocatable :: dnedpsi_chease
  real, dimension(:), allocatable :: ti_chease
  real, dimension(:), allocatable :: dtidpsi_chease
  real, dimension(:), allocatable :: ni_chease
  real, dimension(:), allocatable :: dnidpsi_chease
  real, dimension(:), allocatable :: zeff_chease
  real, dimension(:), allocatable :: signeo_chease
  real, dimension(:), allocatable :: jbsbav_chease
  real, dimension(:,:), allocatable :: g11_chease
  real, dimension(:,:), allocatable :: g12_chease
  real, dimension(:,:), allocatable :: g22_chease
  real, dimension(:,:), allocatable :: g33_chease
  real, dimension(:,:), allocatable :: b_chease
  real, dimension(:,:), allocatable :: dbdpsi_chease
  real, dimension(:,:), allocatable :: dbdchi_chease
  real, dimension(:,:), allocatable :: dpsidr_chease
  real, dimension(:,:), allocatable :: dpsidz_chease
  real, dimension(:,:), allocatable :: dchidr_chease
  real, dimension(:,:), allocatable :: dchidz_chease
  real, dimension(:,:), allocatable :: jacobian_chease
  real, dimension(:,:), allocatable :: r_chease
  real, dimension(:,:), allocatable :: z_chease

  logical, parameter :: debug = .true.

contains
  subroutine read_infile(filename)
    implicit none
    character(len=*), intent(in) :: filename
    open(infile,file= filename)

    read(infile, *)
    read(infile, *) npsi_chease
    if (debug) write(*, *) npsi_chease, "<---npsi_chease"

    read(infile, *)
    read(infile, *) nchi_chease
    if (debug) write(*, *) nchi_chease, "<---nchi_chease"

    read(infile, *)
    read(infile, *) r0exp_chease

    read(infile, *)
    read(infile, *) b0exp_chease

    allocate(psi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) psi_chease

    allocate(chi_chease(nchi_chease))
    read(infile, *)
    read(infile, *) chi_chease

    allocate(rgeom_chease(npsi_chease))
    read(infile, *)
    read(infile, *) rgeom_chease

    allocate(ageom_chease(npsi_chease))
    read(infile, *)
    read(infile, *) ageom_chease

    allocate(q_chease(npsi_chease))
    read(infile, *)
    read(infile, *) q_chease

    allocate(dqdpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dqdpsi_chease

    allocate(d2qdpsi2_chease(npsi_chease))
    read(infile, *)
    read(infile, *) d2qdpsi2_chease

    allocate(p_chease(npsi_chease))
    read(infile, *)
    read(infile, *) p_chease

    allocate(dpdpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dpdpsi_chease

    allocate(f_chease(npsi_chease))
    read(infile, *)
    read(infile, *) f_chease

    allocate(fdfdpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) fdfdpsi_chease

    allocate(v_chease(npsi_chease))
    read(infile, *)
    read(infile, *) v_chease

    allocate(rho_t_chease(npsi_chease))
    read(infile, *)
    read(infile, *) rho_t_chease

    allocate(shear_chease(npsi_chease))
    read(infile, *)
    read(infile, *) shear_chease

    allocate(dsheardpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dsheardpsi_chease

    allocate(kappa_chease(npsi_chease))
    read(infile, *)
    read(infile, *) kappa_chease

    allocate(delta_lower_chease(npsi_chease))
    read(infile, *)
    read(infile, *) delta_lower_chease

    allocate(delta_upper_chease(npsi_chease))
    read(infile, *)
    read(infile, *) delta_upper_chease

    allocate(dvdpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dvdpsi_chease

    allocate(dpsidrhotor_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dpsidrhotor_chease

    allocate(gdpsi_av_chease(npsi_chease))
    read(infile, *)
    read(infile, *) gdpsi_av_chease

    allocate(radius_av_chease(npsi_chease))
    read(infile, *)
    read(infile, *) radius_av_chease

    allocate(r_av_chease(npsi_chease))
    read(infile, *)
    read(infile, *) r_av_chease

    allocate(te_chease(npsi_chease))
    read(infile, *)
    read(infile, *) te_chease

    allocate(dtedpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dtedpsi_chease

    allocate(ne_chease(npsi_chease))
    read(infile, *)
    read(infile, *) ne_chease

    allocate(dnedpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dnedpsi_chease

    allocate(ti_chease(npsi_chease))
    read(infile, *)
    read(infile, *) ti_chease

    allocate(dtidpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dtidpsi_chease

    allocate(ni_chease(npsi_chease))
    read(infile, *)
    read(infile, *) ni_chease

    allocate(dnidpsi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) dnidpsi_chease

    allocate(zeff_chease(npsi_chease))
    read(infile, *)
    read(infile, *) zeff_chease

    allocate(signeo_chease(npsi_chease))
    read(infile, *)
    read(infile, *) signeo_chease

    allocate(jbsbav_chease(npsi_chease))
    read(infile, *)
    read(infile, *) jbsbav_chease

    allocate(g11_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) g11_chease

    allocate(g12_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) g12_chease

    allocate(g22_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) g22_chease

    allocate(g33_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) g33_chease

    allocate(b_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) b_chease

    allocate(dbdpsi_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dbdpsi_chease

    allocate(dbdchi_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dbdchi_chease

    allocate(dpsidr_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dpsidr_chease

    allocate(dpsidz_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dpsidz_chease

    allocate(dchidr_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dchidr_chease

    allocate(dchidz_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) dchidz_chease

    allocate(jacobian_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) jacobian_chease

    allocate(r_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) r_chease

    allocate(z_chease(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) z_chease
  end subroutine read_infile

  subroutine finish
    if (allocated(psi_chease)) deallocate(psi_chease)
    if (allocated(chi_chease)) deallocate(chi_chease)
    if (allocated(rgeom_chease)) deallocate(rgeom_chease)
    if (allocated(ageom_chease)) deallocate(ageom_chease)
    if (allocated(q_chease)) deallocate(q_chease)
    if (allocated(dqdpsi_chease)) deallocate(dqdpsi_chease)
    if (allocated(d2qdpsi2_chease)) deallocate(d2qdpsi2_chease)
    if (allocated(p_chease)) deallocate(p_chease)
    if (allocated(dpdpsi_chease)) deallocate(dpdpsi_chease)
    if (allocated(f_chease)) deallocate(f_chease)
    if (allocated(fdfdpsi_chease)) deallocate(fdfdpsi_chease)
    if (allocated(v_chease)) deallocate(v_chease)
    if (allocated(rho_t_chease)) deallocate(rho_t_chease)
    if (allocated(shear_chease)) deallocate(shear_chease)
    if (allocated(dsheardpsi_chease)) deallocate(dsheardpsi_chease)
    if (allocated(kappa_chease)) deallocate(kappa_chease)
    if (allocated(delta_lower_chease)) deallocate(delta_lower_chease)
    if (allocated(delta_upper_chease)) deallocate(delta_upper_chease)
    if (allocated(dvdpsi_chease)) deallocate(dvdpsi_chease)
    if (allocated(dpsidrhotor_chease)) deallocate(dpsidrhotor_chease)
    if (allocated(gdpsi_av_chease)) deallocate(gdpsi_av_chease)
    if (allocated(radius_av_chease)) deallocate(radius_av_chease)
    if (allocated(r_av_chease)) deallocate(r_av_chease)
    if (allocated(te_chease)) deallocate(te_chease)
    if (allocated(dtedpsi_chease)) deallocate(dtedpsi_chease)
    if (allocated(ne_chease)) deallocate(ne_chease)
    if (allocated(dnedpsi_chease)) deallocate(dnedpsi_chease)
    if (allocated(ti_chease)) deallocate(ti_chease)
    if (allocated(dtidpsi_chease)) deallocate(dtidpsi_chease)
    if (allocated(ni_chease)) deallocate(ni_chease)
    if (allocated(dnidpsi_chease)) deallocate(dnidpsi_chease)
    if (allocated(zeff_chease)) deallocate(zeff_chease)
    if (allocated(signeo_chease)) deallocate(signeo_chease)
    if (allocated(jbsbav_chease)) deallocate(jbsbav_chease)
    if (allocated(g11_chease)) deallocate(g11_chease)
    if (allocated(g12_chease)) deallocate(g12_chease)
    if (allocated(g22_chease)) deallocate(g22_chease)
    if (allocated(g33_chease)) deallocate(g33_chease)
    if (allocated(b_chease)) deallocate(b_chease)
    if (allocated(dbdpsi_chease)) deallocate(dbdpsi_chease)
    if (allocated(dbdchi_chease)) deallocate(dbdchi_chease)
    if (allocated(dpsidr_chease)) deallocate(dpsidr_chease)
    if (allocated(dpsidz_chease)) deallocate(dpsidz_chease)
    if (allocated(dchidr_chease)) deallocate(dchidr_chease)
    if (allocated(dchidz_chease)) deallocate(dchidz_chease)
    if (allocated(jacobian_chease)) deallocate(jacobian_chease)
    if (allocated(r_chease)) deallocate(r_chease)
    if (allocated(z_chease)) deallocate(z_chease)
  end subroutine finish
end module read_chease

