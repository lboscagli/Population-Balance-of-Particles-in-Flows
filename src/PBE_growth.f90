!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ni,index,growth_source,G_term_m)

!**********************************************************************************************
!
! Growth for finite volume method
! Size can be diameter, volume etc.
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
!
!**********************************************************************************************

use pbe_mod
use ice_microphys_mod
use thermo

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(out)              :: growth_source, G_term_m

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

else if (growth_function>=4) then
  ! Ice growth model
  ! compute g_coeff1 and g_coeff2 from kinetic growth model

  if (growth_function==4) then
    !Depositional growth model activation (Karcher et al. 1996)
    if (p_water .ge. p_sat_liq) then !SAC criterion
      call pbe_depositional_growth_ice(index,ni,g_coeff1,g_coeff2)  
    else
      g_coeff1 =0.0
      g_coeff2 =0.0
    endif
  elseif (growth_function==5) then
    !Droplet activation and growth based on Ponsonby et al. 2025
    if (activation_logical) then 
      if ((p_water .ge. p_sat_liq) .and. (p_water .ge. p_sat_ice)) then
        call pbe_droplet_growth(index, ni, g_coeff1,g_coeff2)
      !! TODO 
      !! Here is where I need to check if freezing-relaxation starts based on freezing temperature
      !! Note that freezing temperature depends on liquid volume available for freezing (i.e., depends on (v(index)-v_0))  
      elseif ((p_water .le. p_sat_liq) .and. (p_water .ge. p_sat_ice)) then 
        call pbe_depositional_growth_ice(index, ni, g_coeff1,g_coeff2) 
      else
        g_coeff1 = 0.0
        g_coeff2 = 0.0  
      endif
    else
      call kohler_crit(current_temp, (3.0 / (4.0 * 3.141592653589793E+00) * v0)**(1.0/3.0), kappa, .true., r_vc, S_vc)
      S_vc = S_vc + 1.0
      if (Smw_time_series(step_update) .ge. S_vc) then
        activation_logical = .true.
      else
        g_coeff1 = 0.0
        g_coeff2 = 0.0
      endif
    endif
  endif
  
  g_termr = g_coeff1*(v(index)**g_coeff2)
  g_terml = g_coeff1*(v(index-1)**g_coeff2)  

  g_term_m = 0.5*(g_termr+g_terml)

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




        
    
    