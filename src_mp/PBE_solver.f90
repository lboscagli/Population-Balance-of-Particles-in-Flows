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
double precision, intent(in)                  :: dt

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

  !do i = 1,m

  !  ni(i) = ni(i) * current_rho
  
  !enddo

endif

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot(ni,niprime)
  ni = ni + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot(ni,niprime)
  nitemp = ni + 0.5D0 * niprime * dt
  call pbe_ydot(nitemp,niprime)
  ni = ni + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot(ni,niprime)
  k1 = niprime * dt
  nitemp = ni + 0.5D0 * k1
  call pbe_ydot(nitemp,niprime)
  k2 = niprime * dt
  nitemp = ni + 0.5D0 * k2
  call pbe_ydot(nitemp,niprime)
  k3 = niprime * dt
  nitemp = ni + k3
  call pbe_ydot(nitemp,niprime)
  k4 = niprime * dt
  ni = ni + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

if (growth_function>=4) then
  !do i = 1,m

  !  ni(i) = ni(i) / current_rho
  
  !enddo  
  
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

subroutine pbe_ydot(ni,niprime)

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

double precision, dimension(m), intent(in)  :: ni
double precision, dimension(m), intent(out) :: niprime

double precision dn(m)

double precision growth_source,growth_mass_source,params(1),G_m,tau_g_num,tau_g_den

integer index

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
Loss_Sw_bins = 0.
m_source_bins = 0.

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
        endif
      endif  
    end do
  endif 

  do index = 1,m
    call growth_tvd(ni,index,growth_source,G_m)
    niprime(index) = niprime(index) + growth_source
    tau_g_num = tau_g_num + 3 * ((4/3*pi)**(2/3)) * (v_m(index)*niprime(index)*dv(index))
    tau_g_den = tau_g_den + (4*pi) * (G_m * (v_m(index))**(-1/3)*niprime(index)*dv(index))
  end do
  
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

!**********************************************************************************************