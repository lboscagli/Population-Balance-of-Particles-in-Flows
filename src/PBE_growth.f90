!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ni,index,growth_source)

!**********************************************************************************************
!
! Growth for finite volume method
! Size can be diameter, volume etc.
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(out)              :: growth_source

double precision :: G_terml,G_termr,phi
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Number density in left cell
double precision :: nll               !< Number density in left-left cell
double precision :: nc                !< Number density in this cell
double precision :: nr                !< Number density in right cell
double precision :: nrr               !< Number density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio
double precision :: rl,rr             !< r+ at left and right surface

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************

!Only growth to the right present at nucleation interval

if (growth_function==1) then

  !Size-independent growth
  g_termr = g_coeff1
  g_terml = g_coeff1

else if (growth_function==2) then

  !Power-law growth
  g_termr = g_coeff1*(v(index)**g_coeff2)
  g_terml = g_coeff1*(v(index-1)**g_coeff2)

else if (growth_function==4) then
  ! Ice growth model (Karcher et al 1999)
  ! Luca Boscagli 25/04/205: adapted for CPMOD 

  ! compute g_coeff1 and g_coeff2 from kinetic growth model
  call pbe_growth_ice(index,g_coeff1,g_coeff2)

!  For debug  
!  if (index<5) then
!    write(*,*) 'g_coeff1: ',g_coeff1
!    write(*,*) 'g_coeff2: ',g_coeff2
!  endif
  
  g_termr = g_coeff1*(v(index)**g_coeff2)
  g_terml = g_coeff1*(v(index-1)**g_coeff2)  

end if

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if (g_coeff1>0.D0) then

  ! growth rate is along positive direction
  if (index==1) then

    gnl = 0.0D0
    gnr = g_termr * 0.5 * (ni(1)+ni(2))

  else if (index==m) then

    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1) + 0.5 * phi * (ni(m-1) - ni(m-2)))
    gnr = g_termr * (ni(m) + 0.5*(ni(m) - ni(m-1)))

  else

    ! Fluxes at cell right surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))

    ! Fluxes at cell left surface
    if (index==2 ) then
      gnl = g_terml * 0.5 * (ni(1)+ni(2))
    else
      nl = ni(index-2)
      nc = ni(index-1)
      nr = ni(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
    end if
  end if

else

  ! growth rate is along negative direction
  if (index==1) then

    gnl = g_terml * (ni(1) + 0.5 * (ni(1) - ni(2)))
    rr = (ni(1) - ni(2) + eps) / (ni(2) - ni(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ni(2) + 0.5 * phi * (ni(2) - ni(3)))

  else if (index==m) then

    gnr = 0
    gnl = g_terml * 0.5 * (ni(m)+ni(m-1))

  else

    ! Fluxes at cell right surface
    if (index==m-1) then
      gnr = g_termr * 0.5 * (ni(m)+ni(m-1))
    else
      nl = ni(index)
      nc = ni(index+1)
      nr = ni(index+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      gnr = g_termr * (nc + 0.5 * phi * (nc - nr))
    end if

    ! Fluxes at cell left surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc + 0.5 * phi * (nc - nr))

  end if

end if

if (i_gm==1) then
  ! For mass-conservative growth scheme, apply it after the first interval
  if (index>1) then
    growth_source = (v(index-1)*gnl - v(index)*gnr) / (0.5*(v(index)**2-v(index-1)**2))
  else
    growth_source = (gnl - gnr) / dv(index)
  end if
else
  growth_source = (gnl - gnr) / dv(index)
end if

end subroutine growth_tvd

!**********************************************************************************************

!**********************************************************************************************

subroutine pbe_growth_ice(index,g_coeff1,g_coeff2)

!**********************************************************************************************
!
! Computation of g_coeff1 from kinetic growth model
!
! Based on deposition model by Karcher et al. 1996:
! "The initial composition of jet condensation trails"
! and also based on Ferreira et al. 2024: 
! "Developing a numerical framework for high-fidelity simulation of contrails: sensitivity analysis for conventional contrails"
!    
! Daniel Fredrich 13/01/2025
! Luca Boscagli 25/04/205: adapted for CPMOD 
!
!**********************************************************************************************

  !use arrays, only : p,ajc
  !use chemistry, only : nsp,fsc,temp,sumn,names,wm
  !use euler_part_interface
  use pbe_mod, only :v0, v_m !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
  use pbe_mod, only :amb_temp, amb_p, RH, part_den_l 
  implicit none

  !class(pbe_growth) :: this
  ! integer, intent(in)  :: ijk
  integer, intent(in)  :: index
  double precision, intent(out)              :: g_coeff1, g_coeff2
  !integer :: isp
  double precision :: p_water,p_water_sat_ice,amb_rho,p_water_sat_liq !,RH,amb_temp,amb_p
  double precision :: M_water !,M_air, X_water
  double precision :: r_part,r_nuc,den_ice
  double precision :: dif_water,lambda_water
  double precision :: coll_factor,Kn,alpha
  double precision :: fornow,drdt
  double precision :: part_den !, part_den_l, part_den_r
  real,parameter :: pi = 3.141592653589793E+00
  double precision :: gascon=8314.3
  !real, allocatable, dimension(:)                 :: part_den(:) ! Density of particles in an arbitrary bin of the PSD (interpolated value)

  
  !----------------------------------------------------------------------------------------------

  !do isp = 1,nsp          
  !  if (trim(names(isp)) .eq. 'H2O') then
  !    X_water = fsc(isp,ijk) / sumn(ijk) ! water molecular fraction
  !    M_water = wm(isp) ! water molecular weigth
  !  endif
  !enddo

  ! X_water = 0.0815 ! water molecular fraction
  !RH = 1.2 ! Relative humidity
  M_water = 18.016 ! water molecular weigth
  ! M_air = 28.96 ! air molecular weigth
  !amb_temp = 208.15 ! ambient temperature in kelvin
  !amb_p = 16235.70 ! ambient pressure in Pascal
  !amb_rho = amb_p/(gascon/M_air)/amb_temp ! ambient air density
  !part_den_r = 1550.0				! Density of particles on the right side of the PSD (kg/m^3)
  !part_den_l = 1550.0				! Density of particles on the left side of the PSD (kg/m^3)
  !v0 = v_nuc ! 3.35103e-23 
  
  p_water_sat_liq = exp(54.842763 - 6763.22 / amb_temp - 4.21 * log(amb_temp) &
  + 0.000367 * amb_temp + tanh(0.0415 * (amb_temp - 218.8)) &
  * (53.878 - 1331.22 / amb_temp - 9.44523 * log(amb_temp) &
  + 0.014025 * amb_temp))
  
  !write(*,*) 'p_water_sat_liq: ',p_water_sat_liq

  ! water vapor partial pressure
  !p_water = p * X_water
  p_water = p_water_sat_liq * RH
  
  !write(*,*) 'p_water: ',p_water

  ! saturated (relative to ice) water vapour partial pressure 
  p_water_sat_ice = exp(9.550426 - 5723.265 / amb_temp + 3.53068 * log(amb_temp) &
                    - 0.00728332 * amb_temp) 
  
  !write(*,*) 'p_water_sat_ice: ',p_water_sat_ice

  den_ice = 917.0 ! constant ice particle density (neglecting soot core) 

  r_part = (3.0 / (4.0 * pi) * v_m(index))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
  r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        ! radius of the nuclei (samllest particle volume) - this is constant

  !write(*,*) 'r_nuc: ',r_nuc
  part_den = (part_den_l * r_nuc**3.0 + den_ice * &
                          (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density

  !write(*,*) 'part_den: ',r_nuc

  dif_water = 2.11D-5 * 101325.0 / amb_p * (amb_temp / 273.15)**1.94  ! diffusion coefficient of water vapor molecules in air 
  lambda_water = 6.15D-8 * 101325.0 / amb_p * amb_temp / 273.15       ! water vapor mean free path

  Kn = lambda_water / r_part  ! knudsen number 
  alpha = 0.1 ! deposition coefficient. Higher values of alpha speed up the deposition rate

  coll_factor = 1.0 / (1.0 / (1.0 + Kn) + 4.0 / 3.0 * Kn / alpha) ! collision factor (G), accounts for transition from gas kinetic energy (G->1 for Kn->0) to continuum regime (G->0 for Kn->1)

  fornow = dif_water * coll_factor * M_water * (p_water - p_water_sat_ice) &
          / (gascon * amb_temp)   ! nominator of dr/dt 
                                  ! [ M_water * p_water_sat_ice / (gascon * amb_temp)]: saturated (relative to ice) water vapor density (considering compressibility factor approx = 1)
                                  ! This is calling gascon = 8314.3 J/K/kg (in module_chemistry.f90)
                                  ! gascon/M_water is the specific gas constant for water vapor, i.e., Rv = 461.52 J/K/kg
                                  

  drdt = fornow / (part_den * r_part) ! change in particle radius over time

  g_coeff2 = 2.0 / 3.0 
  g_coeff1 = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt

  

!! Luca: the section below is not needed for CPMOD standalone (no coupling with BOFFIN)
!  dmdt = 4.0 * pi * r_part * fornow

!  m_source(ijk) = m_source(ijk) + dmdt * this%ni_part_pbe(index) &
!                * this%dv_part(index) * ajc(ijk)                       ! [kg/s] !ajc is the jacobian of the coordinate transformation (BOFFIN manual p.16 and p.41)
                                                                                !Presumably if the computational grid is already orthogonal, then this is an identity matrix?


end subroutine pbe_growth_ice

!-------------------------------------------------------------------------------