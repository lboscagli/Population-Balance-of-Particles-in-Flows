!**********************************************************************************************
!
! PBE solver subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ(ni,dt)

!**********************************************************************************************
!
! Temporal integration
!
! Stelios Rigopoulos (02/06/2019)
!
! Modified 14/06/06
! Modified 19/12/2017
! Modified 13/05/2018
! Modified 02/06/2019
! Modified 25/06/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(inout) :: ni
double precision, intent(inout)                  :: dt

integer i

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------
if (growth_function>=4) then

  !Eliminate the negative values of number density which are not physical - Luca 27/04/2025
  do i = 1,m
    if (ni(i) < 0.0) then
      ni(i) = 0.0
    endif
  enddo 

endif

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot(ni,dt,niprime)
  ni = ni + niprime * dt 

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot(ni,dt,niprime)
  nitemp = ni + 0.5D0 * niprime * dt
  call pbe_ydot(nitemp,dt,niprime)
  ni = ni + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot(ni,dt,niprime)
  k1 = niprime * dt
  nitemp = ni + 0.5D0 * k1
  call pbe_ydot(nitemp,dt,niprime)
  k2 = niprime * dt
  nitemp = ni + 0.5D0 * k2
  call pbe_ydot(nitemp,dt,niprime)
  k3 = niprime * dt
  nitemp = ni + k3
  call pbe_ydot(nitemp,dt,niprime)
  k4 = niprime * dt
  ni = ni + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

if (growth_function>=4) then 
  
  !Eliminate the negative values of number density which are not physical - Luca 27/04/2025
  do i = 1,m
    if (ni(i) < 0.0) then
      ni(i) = 0.0
    endif
  enddo 

endif
 

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot(ni,dt,niprime)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! 14/01/2002
! Modified 04/05/2017
! Modified 23/06/2020
! Modified 25/06/2020
! Modified 05/07/2020
!
!**********************************************************************************************

use pbe_mod
use thermo

implicit none

double precision, dimension(m), intent(inout)  :: ni
double precision, intent(inout)  :: dt
double precision, dimension(m), intent(out) :: niprime

double precision dn(m)
double precision, dimension(m)  :: ni_tmp
double precision :: courant_max, Gdv_max

double precision growth_source,growth_mass_source,params(1),G_m,tau_g_num,tau_g_den,p_water_sat_liq

integer index, jj

real,parameter :: pi = 3.141592653589793E+00

!----------------------------------------------------------------------------------------------

niprime = 0.
params(1) = 0.

!Initialize kinetics timescale
tau_g_num = 0. ! numerator of growth timescale 
tau_g_den = 0. ! denominator of growth timescale

!Initialize Supersaturation consumption
Loss_Sw = 0.
g_coeff1_l_prev = 0.
gnl_prev = 0.
Loss_Sw_bins = 0.
m_source_bins = 0.
m_source_pbe = 0.
Loss_Sw_bins = 0.
if (step_update .eq. 1) then
  courant_max = 0.0
  Gdv_max = 0.0
  ratio_max = 0.0
endif

!Initialize ni_tmp
!ni_tmp = ni

!Nucleation
if (max_nuc>0) then
! Note: expressions for i_gm==0 and i_gm==1 are equivalent but written separately for clarity
  if (i_gm==0) then
    do index = 1,max_nuc
      niprime(index) = nuc(index)/dv(index)
    end do
  else if (i_gm==1) then
    do index = 1,max_nuc
      niprime(index) = nuc(index)*v_m(index)/(0.5*(v(index)**2-v(index-1)**2))
    end do
  end if
end if

!Growth
if (growth_function>0) then

  if (growth_function .eq. 7) then
   !Check k-kholer condition for each index if it is a nuclei
    do index = 1,m
      if (nuclei_logical(index) .and. (.not. activation_logical_bins(index))) then
        call kohler_crit(current_temp, (3.0 / (4.0 * 3.141592653589793E+00) * v0_bins(index) )**(1.0/3.0), kappa_bins(index), .false., r_vc, S_vc)
        S_vc = S_vc + 1.0
        if ((Smw_time_series(step_update) .ge. S_vc)) then  
          !activation_logical = .true.
          activation_logical_bins(index) = .true.
          S_vc_bins(index) = S_vc
          write(*,*) 'S_vc [-]',S_vc
          write(*,*) 'r_vc [nm]',r_vc*1E9
          
          ! if (((3.0 / (4.0 * 3.141592653589793E+00) * v(index) )**(1.0/3.0)) < r_vc) then
          !Nucleation
          ! write(*,*) 'WARNING: nucleation'
          ! niprime(index) = niprime(index) + (((4.0 * 3.141592653589793E+00) / 3.0 * (r_vc**3)) - v_m(index))/dt/dv(index)
          ! write(*,*) 'growth due to nucleation',niprime(index)
          ! endif

          if (((3.0 / (4.0 * 3.141592653589793E+00) * v(index) )**(1.0/3.0)) < r_vc) then
            write(*,*) 'WARNING: nuclei required shit from bin index',index
            do jj = index,m
              if (((3.0 / (4.0 * 3.141592653589793E+00) * v(jj) )**(1.0/3.0)) >= r_vc) then
                ni(jj) = ni(jj) + ni(index)*(dv(index)/dv(jj))!*rho_w(index)/part_rho_bins(index))
                ni(index) = 0.0
                exit
              endif
            end do
            write(*,*) 'WARNING: nuclei shifted to bin index ',jj
          endif

        endif
      endif  
    end do
  endif 

  call p_sat_liq_murphy_koop(p_water_sat_liq)
  do index = 1,m
    call growth_tvd(ni,index,growth_source,G_m)
    ! update max growth rate
    courant_max = max(courant_max,G_m*dt/dv(index))
    if (ni(index)>0) then
      Gdv_max = max(Gdv_max,G_m/dv(index))!growth_source/ni(index))!G_m/dv(index))
      ! write(*,*) 'Gdv_max',Gdv_max
      ! write(*,*) 'courant_max',G_m/dv(index)
    endif

    niprime(index) = niprime(index) + growth_source 
    
    if (ni_type(index) .eq. 1.0) then
      m_source_pbe(index) = rho_w * niprime(index) * v_m(index) * dv(index)
    elseif (ni_type(index) .eq. 2.0) then
      m_source_pbe(index) =  rho_ice * niprime(index) * v_m(index) * dv(index)    
    else
      m_source_pbe(index) = 0.0
      ! if (inps_distribution_logical) then
      !   m_source_pbe(index) = part_rho_bins(index) * niprime(index) * v_m(index) * dv(index)   
      ! else 
      !   m_source_pbe(index) = part_den_l * niprime(index) * v_m(index) * dv(index) 
      ! endif
    endif
    !Loss_Sw_bins(index) = (amb_p/p_water_sat_liq/epsilon_fluid/amb_rho) * m_source_pbe(index)
    if (m_source_pbe(index) .ge. 0) then
      Loss_Sw_bins(index) = (amb_p/p_water_sat_liq) * (Rv/Rd) * (1.0/amb_rho) * m_source_pbe(index)
    else
      m_source_pbe(index) = 0.0
      Loss_Sw_bins(index) = 0.0
    endif
    tau_g_num = tau_g_num + 3 * ((4/3*pi)**(2/3)) * (v_m(index)*niprime(index)*dv(index))
    tau_g_den = tau_g_den + (4*pi) * (G_m * (v_m(index))**(-1/3)*niprime(index)*dv(index))
  end do
  Loss_Sw = sum(Loss_Sw_bins)

  if (Gdv_max*dt_input > 0.25) then !.and. Smw_time_series(step_update)>1.5
    !if 
    write(*,*) 'CLF old',Gdv_max*dt_input
    dt = min(0.25/Gdv_max,dt_input) !(dt*0.25)/(courant_max)
    write(*,*) 'CLF new',Gdv_max*dt
    ! write(*,*) 'dt',dt
    ! if ((Smw_time_series(step_update)) > 1.0) then
    !   ratio_max = max(ratio_max, abs(((Smw_time_series(step_update) - Smw_time_series(step_update-1))/dt_input) / Gdv_max ))
    ! endif 
    ! write(*,*) 'Gdv_max',Gdv_max
    ! write(*,*) 'dSmw/dt',(Smw_time_series(step_update) - Smw_time_series(step_update-1))/dt_input
  else
    ! if (Gdv_max > 0.0) then
    !   dt = min(0.95/Gdv_max,dt_input)
    ! else
      dt = dt_input
    ! endif
  endif

  
  if (i_gm==1) then
    ! For mass-conservative growth scheme, apply growth source term after the first interval
    do index=2,m
      niprime(index) = niprime(index) + (ni(index)*g_coeff1 &
      /((g_coeff2+1.)*0.5*(v(index)**2-v(index-1)**2))) & 
      * (v(index)**(g_coeff2+1.)-v(index-1)**(g_coeff2+1.))
    end do
  end if
end if

if (growth_function>=4) then
  if (tau_g_den>0) then
    tau_g = tau_g_num / tau_g_den
  else
    tau_g = 0.
  endif
endif

!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni,niprime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime,ni)
end if

end subroutine pbe_ydot

subroutine p_sat_liq_murphy_koop(p_sat_liq)

!**********************************************************************************************
!
! Determines the saturation vapor pressure of water as a function of temperature using the Murphy-Koop equation 
!
! Murphy, D. M., & Koop, T. (2005)
!
!**********************************************************************************************
  use pbe_mod, only :current_temp

  implicit none

  double precision, intent(out)                  :: p_sat_liq
  
  !----------------------------------------------------------------------------------------------
  

  p_sat_liq = (EXP(54.842763 - 6763.22/current_temp - 4.210 * LOG(current_temp) + 0.000367 * current_temp + TANH(0.0415 * (current_temp - 218.8)) * \
  (53.878 - 1331.22/current_temp - 9.44523 * LOG(current_temp) + 0.014025 * current_temp)))
  
  
end subroutine p_sat_liq_murphy_koop

!**********************************************************************************************