!**********************************************************************************************
!
! PBE general subroutines
! Stelios Rigopoulos, Anxiong Liu, Binxuan Sun, Daniel O'Sullivan
! Modified 19/12/2017
! Modified 15/01/2018
! Modified 09/05/2018
!
!**********************************************************************************************



!**********************************************************************************************

module pbe_mod

!**********************************************************************************************
!
! Declaration of common variables related to grid and kernels
!
! by Stelios Rigopoulos
! 17/12/2001
! 19/04/2002 new version
!
! Modified 07/05/2017
! Modified 31/05/2017
! Modified 14/07/2020
!
!**********************************************************************************************

implicit none

save

double precision, allocatable, dimension(:) :: v
double precision, allocatable, dimension(:) :: dv
double precision, allocatable, dimension(:) :: v_m
double precision, allocatable, dimension(:) :: nuc

double precision v0,grid_lb,grid_rb
double precision agg_kernel_const
double precision break_const
double precision g_coeff1,g_coeff2
double precision nuc1
double precision N0

double precision :: amb_temp, amb_p, amb_rho, G_mixing_line, part_den_l, alpha_ice, p_sat_liq, p_sat_ice
double precision :: jet_cl_model, diameter_jet, u_0j, T_0j, current_temp, current_rho, p_water, current_XH2O
double precision :: tau_g
double precision :: kappa !hygroscopicity
real :: Loss_Sw !Saturation consumption
real :: Smw !Saturation ratio along the mixingline with no particles
real(kind=8) :: S_vc !Critical saturation ratio for aerosol droplet with specific dry radius and hygroscopicity
real(kind=8) :: r_vc !Critical radius for aerosol particle to activate - this is dummy as we actually specify a wet radius for the nuclei 
real, allocatable :: Smw_time_series(:)

logical :: activation_logical !Water droplet activation flag (initialize as .false. and then check saturation)
logical :: consumption_logical !Flag for user activation of consumption of supersaturation based on K15 model

integer m,grid_type
integer i_gm,solver_pbe
integer initdis,n_th1,n_th2
integer agg_kernel
integer growth_function
integer max_nuc
integer order_of_gq
integer step_update

end module pbe_mod

!**********************************************************************************************



!**********************************************************************************************

module agg_cfv_mod

implicit none

save

integer :: ID_Ajk, ID_xj, ID_xk, ID_dgj, ID_dgk
integer :: ID_j, ID_k, ID_i
integer :: NDou_AggSec, NInt_AggSec
integer :: Nbeta_Agg
integer :: N_AggSec
integer :: index_Sec

double precision, allocatable, dimension(:) :: Dou_AggSec
integer, allocatable, dimension(:)          :: Int_AggSec
double precision, allocatable, dimension(:) :: beta_AggSec
double precision, allocatable, dimension(:) :: beta_AggSecOrig

double precision :: Pbe_AggrInfo(3)

parameter (ID_Ajk = 1, ID_xj = 2, ID_xk = 3, ID_dgj=4, ID_dgk=5)
parameter (ID_j = 1, ID_k = 2, ID_i = 3)
parameter (NDou_AggSec = 5, NInt_AggSec = 3)

end module agg_cfv_mod

!**********************************************************************************************



!**********************************************************************************************

module frag_cfv_mod

implicit none

save

double precision :: break_a, break_b, break_c, h, delta, lambda ! break dist parameters
integer :: break_section_count, break_global_cntr
double precision, allocatable :: break_kernel_store(:,:)

end module frag_cfv_mod

!**********************************************************************************************



!**********************************************************************************************

module gauss_mod

! Data for Gaussian integration

implicit none

save

double precision, allocatable, dimension(:,:) :: gw, gp

! Gauss data 2-point
double precision, parameter, dimension(4,2) :: gp_2 &
& = reshape((/-1.D0/(3.D0**0.5D0),-1.D0/(3.D0**0.5D0), &
&1.D0/(3.D0**0.5D0),1.D0/(3.D0**0.5D0),-1.D0/(3.D0**0.5D0),1.D0/3.D0**0.5D0,-1.D0/(3.D0**0.5D0),1.D0/(3.D0**0.5D0)/),shape(gp_2))
double precision, parameter :: gw_2 = 1.D0

! Gauss data 3 point
double precision, parameter, dimension(9,2) :: gp_3 &
& = reshape((/ -(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,&
&0.D0,0.D0,0.D0,(6.D0/10.D0)**0.5D0,(6.D0/10.D0)**0.5D0,(6.D0/10.D0)**0.5D0,&
&-(6.D0/10.D0)**0.5D0,0.D0,(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,&
& 0.D0,(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,0.D0,(6.D0/10.D0)**0.5D0/),shape(gp_3))
double precision, parameter, dimension(9,2) :: gw_3 &
& = reshape((/ 5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 , 8.D0/9.D0 , 8.D0/9.D0,&
&5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 , 5.D0/9.D0,&
& 5.D0/9.D0 , 8.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 ,5.D0/9.D0/),shape(gw_3))

end module gauss_mod

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_read(n_pbe)

!**********************************************************************************************
!
! Initialises PBE data
!
! Stelios Rigopoulos 03/12/2014
! Modified 06/05/2017
! Modified 25/06/2020
! Modified 14/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, intent(out) :: n_pbe

integer i

!----------------------------------------------------------------------------------------------

! Read PBE parameters

open(30,file='pbe/pbe.in')
do i=1,2
  read(30,*)
end do
read(30,*) agg_kernel
read(30,*) agg_kernel_const
read(30,*) max_nuc
read(30,*) nuc1
read(30,*) growth_function
read(30,*) g_coeff1
read(30,*) g_coeff2
read(30,*) i_gm
read(30,*) break_const
read(30,*) order_of_gq
read(30,*) grid_type
read(30,*) m
read(30,*) grid_lb
read(30,*) grid_rb
read(30,*) v0
read(30,*) initdis
read(30,*) n_th1
read(30,*) n_th2
read(30,*) N0
read(30,*) solver_pbe

n_pbe = m

if ((agg_kernel>0).and.(growth_function==1)) then
  write(*,*) 'mass-conservative discretisation should not be used with size-independent &
  growth, stopping'
  stop
end if

end subroutine pbe_read

!**********************************************************************************************


!**********************************************************************************************

subroutine pbe_ice_read()

!**********************************************************************************************
!
! Initialises ICE data
!
! Luca Boscagli 29/04/2025
!
!**********************************************************************************************
  use pbe_mod

  implicit none
  
  double precision :: gascon=8314.3
  double precision :: M_air=28.96

  integer i
  
  !----------------------------------------------------------------------------------------------
  
  ! Read ICE input data
  open(30,file='psr/ice.in')
  do i=1,2
    read(30,*)
  end do
  read(30,*) amb_p
  read(30,*) amb_temp
  read(30,*) G_mixing_line
  read(30,*) part_den_l
  read(30,*) alpha_ice
  read(30,*) jet_cl_model
  read(30,*) diameter_jet
  read(30,*) u_0j
  read(30,*) T_0j
  read(30,*) kappa
  read(30,*) consumption_logical
  close(30)

  amb_rho = amb_p / amb_temp / (gascon / M_air) 
  
end subroutine pbe_ice_read

subroutine pbe_ice_update(time, jet_temp, jet_rho)

  !**********************************************************************************************
  !
  ! Update ICE data thermodynamic state based on Witze1974 + Karcher2015 jet centerline 
  ! velocity and cooling rate model
  !
  ! Luca Boscagli 06/05/2025
  !
  !**********************************************************************************************
    use pbe_mod
  
    implicit none

    double precision, intent(in)                  :: time
    double precision, intent(out)                  :: jet_temp, jet_rho
    double precision :: epsilon_t, beta, x_m, r_0j, tau_m, p_water_amb

    double precision :: gascon=8314.3
    double precision :: M_air = 28.96
    
    !----------------------------------------------------------------------------------------------
    
    ! Read ICE input data
    
    !Initial Jet radius
    r_0j = 0.5*diameter_jet 
    
    !Model constants
    epsilon_t = 0.0285 ! dimensioless turbulent diffusivity (Tollmien, 1926)
    beta = 0.9 ! Dilution parameter (Karcher, 1999)
    
    !Maximum distance over which central jet region is unaffected by entrainement
    x_m = r_0j * (2.0/epsilon_t)**0.5
    
    !Mixing timescale
    tau_m = x_m / u_0j
    
    !Jet cooling rate
    if (time.ge.tau_m) then
      jet_temp = amb_temp + (T_0j-amb_temp) * (tau_m/time)**beta
      !write(*,*) 'amb_temp'
    else
      jet_temp = T_0j
    endif


    ! Compute saturation ratio based on user defined slope of mixing line (G_mixing_line) 
    p_sat_ice = EXP(9.550426 - 5723.265/jet_temp + 3.53068 * LOG(jet_temp) -  0.00728332 * jet_temp )
    p_sat_liq = (EXP(54.842763 - 6763.22/jet_temp - 4.210 * LOG(jet_temp) + 0.000367 * jet_temp + TANH(0.0415 * (jet_temp - 218.8)) * \
    (53.878 - 1331.22/jet_temp - 9.44523 * LOG(jet_temp) + 0.014025 * jet_temp)))
    
    ! Water vapor partial pressure at ambient condition
    p_water_amb = 1.0 * (EXP(9.550426 - 5723.265/amb_temp + 3.53068 * LOG(amb_temp) -  0.00728332 * amb_temp)) ! this is assuming that at ambient condition we are saturated relative to ice, hence 1.0
    
    !Saturation ratio along the mixing line
    Smw = ((p_water_amb + G_mixing_line * (jet_temp - amb_temp)) / p_sat_liq)

    if (step_update .eq. 0) then
      !Water vapor partial pressure along the mixing line
      p_water = Smw*p_sat_liq
    else
      p_water = Smw_time_series(step_update)*p_sat_liq
    endif
  
    ! Jet moist air density
    jet_rho = amb_p/(gascon/M_air*(jet_temp*(1.0_8+0.61_8*((p_water/p_sat_liq)*0.622_8*((611.2_8*exp(17.67_8*(jet_temp-273.15)/((jet_temp-273.15)+243.5_8)))/amb_p)))))
    !write(*,*) 'Moist air density',jet_rho
    ! Jet dry air density
    !jet_rho = amb_p / jet_temp / (gascon / M_air)
    !write(*,*) 'Dry air density',jet_rho  

  end subroutine pbe_ice_update

!**********************************************************************************************

  subroutine pbe_ice_update_LES(time, jet_temp, jet_rho, jet_XH2O)

    !**********************************************************************************************
    !
    ! Update ICE data thermodynamic state and water vapor molar fraction based on LES data 
    !
    ! Luca Boscagli 29/05/2025
    !
    !**********************************************************************************************
    use pbe_mod
  
    implicit none

    double precision, intent(in)                  :: time
    double precision, intent(out)                  :: jet_temp, jet_rho, jet_XH2O
    double precision :: a_T, b_T, a_XH2O, b_XH2O, c_XH2O, d_XH2O

    double precision :: gascon=8314.3
    double precision :: M_air = 28.96
    
    !----------------------------------------------------------------------------------------------
      
    ! Read ICE input data
    
    !Parameters from fitting of LES data - this are case dependent and not general
    a_T = 0.00128922
    b_T = 0.77394931
    a_XH2O = 1.15187008e-03
    b_XH2O = 7.22690952e-01
    c_XH2O = 2.37237876e-02
    d_XH2O = -1.14574449e-04
    
    !temperature
    jet_temp = amb_temp + (T_0j - amb_temp)*(a_T/(time+a_T))**b_T 
    
    !Water vapor molar fraction
    jet_XH2O = d_XH2O + (c_XH2O - d_XH2O)*(a_XH2O/(time+a_XH2O))**b_XH2O
    
    !Water vapor partial pressure
    p_water = amb_p*jet_XH2O

    !Saturation pressure relative to ice and liquid
    p_sat_ice = EXP(9.550426 - 5723.265/jet_temp + 3.53068 * LOG(jet_temp) -  0.00728332 * jet_temp )
    p_sat_liq = (EXP(54.842763 - 6763.22/jet_temp - 4.210 * LOG(jet_temp) + 0.000367 * jet_temp + TANH(0.0415 * (jet_temp - 218.8)) * \
    (53.878 - 1331.22/jet_temp - 9.44523 * LOG(jet_temp) + 0.014025 * jet_temp)))

    !Saturation ratio and Water vapor partial pressure along the mixing line
    Smw = p_water / p_sat_liq

    if (consumption_logical) then
      if (step_update .eq. 0) then
        p_water = Smw*p_sat_liq
      else
        p_water = Smw_time_series(step_update)*p_sat_liq
        Smw = p_water/p_sat_liq
      endif     
    endif

    ! Jet moist air density
    jet_rho = amb_p/(gascon/M_air*(jet_temp*(1.0_8+0.61_8*((p_water/p_sat_liq)*0.622_8*((611.2_8*exp(17.67_8*(jet_temp-273.15)/((jet_temp-273.15)+243.5_8)))/amb_p)))))

  
  end subroutine pbe_ice_update_LES
  
  !**********************************************************************************************
  


!**********************************************************************************************

subroutine pbe_init(ni)

!**********************************************************************************************
!
! Initialises PBE data
!
! Stelios Rigopoulos 03/12/2014
! Modified 06/05/2017
! Modified 25/06/2020
! Modified 14/07/2020
!
!**********************************************************************************************

use pbe_mod
use frag_cfv_mod
use gauss_mod

implicit none

double precision, dimension(m), intent(inout) :: ni


integer i

!----------------------------------------------------------------------------------------------

! Initialise nucleation
if (max_nuc>0) then
  nuc = 0.D0
  nuc(1:max_nuc) = nuc1
end if

! Initialise aggregation
if (agg_kernel>0) then
  call PBE_agg_fvLocate(v,1)
  call PBE_agg_fvLocate(v,2)
  call PBE_agg_beta(1)
  call PBE_agg_beta(2)
end if

if (break_const>0.) then
  call gauss_init(order_of_gq)
  call pbe_breakage_calc(1)
  allocate(break_kernel_store(break_global_cntr,3))
  call pbe_breakage_calc(2)
end if

! Initialise distribution
if (initdis==1) then
  ! Zero
  ni = 0.D0
else if (initdis==2) then
  ! Top hat
  ni = 0.D0
  ni(n_th1:n_th2) = N0
end if

end subroutine pbe_init

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_grid()

!**********************************************************************************************
!
! Subroutine grid
! Generates/reads the grid
!
! By Stelios Rigopoulos
! 31/10/2001 - initial version
! Modified 01/06/2017
! Modified 12/05/2018
! Modified 14/07/2020
!
! m                     number of points
! ni                    population density
! v                     independent variable - volume
! dv                    lengths of intervals
! v_m                   interval mid-points
! v0                    size of nuclei
! grid_lb               left boundary
! grid_rb               right boundary
!
!**********************************************************************************************e

use pbe_mod

implicit none

double precision alpha,v1,v2
double precision :: alpha_l,alpha_mid,alpha_r
double precision :: grid_mid_l,grid_mid_r, v3, v4

integer i
integer :: m_l
integer :: m_r

!----------------------------------------------------------------------------------------------

! Allocate arrays
allocate(v(0:m),dv(m),v_m(m),nuc(m))

if (grid_type==1) then

  !Option 1: geometric grid
  v(0) = grid_lb
  v(1) = v0 + (v0 - grid_lb)
  v1 = v(1) - grid_lb
  v2 = grid_rb - grid_lb
  call inc_ratio(v1,v2,m,alpha)
  do i=2,m
    v(i) = v(0)+(v(1)-v(0))*(1-alpha**real(i))/(1-alpha)
  end do
  write(*,*) 'left boundary: ',v(0)
  write(*,*) 'first node: ',v(1)
  write(*,*) 'right boundary: ',v(m)
  write(*,*) 'ratio: ',alpha

else if (grid_type==2) then

  !Option 2: uniform grid
  alpha = (grid_rb - grid_lb)/m
  v(0) = grid_lb
  do i=1,m
    v(i) = v(i-1) + alpha
  end do
  write(*,*) 'left boundary: ',v(0)
  write(*,*) 'first node: ',v(1)
  write(*,*) 'right boundary: ',v(m)
  write(*,*) 'increment: ',alpha

else if (grid_type==3) then
  !Option 3: composite grid

  ! Size the grid
  m_l = int(m*0.3)
  m_r = int(m*0.5)

  ! Define intermediate points
  grid_mid_l = (v0 + (v0 - grid_lb)) * 50.0
  grid_mid_r = (grid_mid_l + grid_rb) * 0.5

  write(*,*) 'grid_lb: ',grid_lb
  write(*,*) 'grid_mid_l: ',grid_mid_l
  write(*,*) 'grid_mid_r: ',grid_mid_r
  write(*,*) 'grid_rb: ',grid_rb

  ! First part: geometric progression
  v(0) = grid_lb
  v(1) = v0 + (v0 - grid_lb)
  v1 = v(1) - grid_lb
  v2 = grid_mid_l - grid_lb
  call inc_ratio(v1,v2,m_l,alpha_l)
  do i=2,m_l
    v(i) = v(0)+(v(1)-v(0))*(1-alpha_l**real(i))/(1-alpha_l)
    !write(*,*) 'i: ',i 
    !write(*,*) 'v(i): ',v(i)
  end do
  
  !write(*,*) 'left boundary: ',v(0)
  !write(*,*) 'first node: ',v(1)
  !write(*,*) 'right boundary mid-1: ',v(m_l)

  ! Second (middle) part: uniform
  alpha_mid = (grid_mid_r - grid_mid_l)/m_r
  do i=m_l+1,m_l+m_r
    v(i) = v(i-1) + alpha_mid
    !write(*,*) 'i: ',i 
    !write(*,*) 'v(i): ',v(i)    
  end do

  write(*,*) 'right boundary mid-2: ',v(m_l+m_r)   

  ! Third (final) part: geometric progression
  v3 = grid_mid_r - v(m_l+m_r-1)
  v4 = grid_rb - v(m_l+m_r-1)
  call inc_ratio(v3,v4,m-(m_l+m_r),alpha_r)
  do i=m_l+m_r+1,m
    v(i) = v(m_l+m_r-1)+(grid_mid_r-v(m_l+m_r-1))*(1-alpha_r**real(i-m_l-m_r))/(1-alpha_r)
    !write(*,*) 'i: ',i 
    !write(*,*) 'v(i): ',v(i)    
  end do
 
  !write(*,*) 'alpha_l: ',alpha_l
  !write(*,*) 'alpha_r: ',alpha_r


end if

!Interval length
do i=1,m
  dv(i) = v(i)-v(i-1)
end do

!Determine mid-points
do i=1,m
  v_m(i) = v(i-1)+0.5D0*dv(i)
end do

end subroutine pbe_grid

!**********************************************************************************************



!**********************************************************************************************

subroutine inc_ratio(lb,rb,k,q)

!**********************************************************************************************
!
! Calculation of common ratio in exponential grid
!
! By Binxuan Sun
! 28/02/2018
! lb: left boundary, which is also the first term a1 in geometric progression
! rb: right boundary
! k:  number of total elements
! q:  common ratio
!
! Sum a(1,2,3,....,m) = a1 * ( 1 - q^m )/ ( 1 - q )
! Sum = the length of the interval [0,rb]
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)    :: lb
double precision, intent(in)    :: rb
integer, intent(in)             :: k
double precision, intent(inout) :: q

double precision s,r1,r2,rm,rt
double precision t1,t2,a1,t
double precision abs_error

integer i,IMAX

!**********************************************************************************************

IMAX = 1000 !maximum iteration times
r1 = 1.0001
r2 = 20
abs_error= 1.0D-6

! a1 is the first number in the series
a1 = lb
!s = ( rb - lb )/ a1
s = rb / a1
t1 = r1**k - s * r1 + s -1
t2 = r2**k - s * r2 + s -1

if ((t1.GE.0.0).OR.(t2.LE.0.0)) then
  write(*,*) "Error in the right grid boundary and the number of nodes, re-select them"
  stop
end if
rt = 0

if ((t1.LT.0.0).AND.(t2.GT.0.0)) then
  do i=1,IMAX
    rm = ( r1+ r2 )*0.5
    if (abs(rm-rt).LE.abs_error) then
      go to 666
    end if
    rt = rm
    q = rm
    t = q**k - s * q + s - 1
    if (t.GT.0) then
      r2 = rm
    else
      r1 = rm
    end if
  end do
end if

666 q = rm

write(*,*) "alpha = ", q

end subroutine inc_ratio

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_moments(ni,moment,meansize)

!**********************************************************************************************
!
! Calculation of zeroth and first moments
!
! By Stelios Rigopoulos
! Modified 06/05/2017
! Modified 04/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)    :: ni
double precision, dimension(0:1), intent(out) :: moment
double precision, intent(out)                 :: meansize

double precision M1_lp,lp

integer i

!----------------------------------------------------------------------------------------------

moment(0) = 0.0
moment(1) = 0.0

do i=1,m
  moment(0) = moment(0) + ni(i)*dv(i)
  moment(1) = moment(1) + 0.5D0*(v(i-1)+v(i))*ni(i)*dv(i)
end do

M1_lp = 0.D0
do i=m-5,m
  M1_lp = M1_lp + 0.5D0*(v(i-1)+v(i))*ni(i)*dv(i)
end do

lp = M1_lp/moment(1)
if (lp.gt.0.001) then
  write(*,*) 'warning, more than 0.1% of mass in the last five nodes'
end if

if (moment(0).gt.1.D-10) then
  meansize = moment(1)/moment(0)
else
  meansize = 0.D0
end if

end subroutine pbe_moments

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output(ni,i_writesp,i_step)

!**********************************************************************************************
!
! Calculation of zeroth and first moments
!
! By Stelios Rigopoulos
! Modified 06/05/2017
! Modified 04/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in) :: i_writesp
integer, intent(in) :: i_step

double precision :: nitemp(m),eta(m),psi(m)
double precision, dimension(0:1) :: moment

double precision meansize

integer i

character(len=10) :: i_step_str


!----------------------------------------------------------------------------------------------

call pbe_moments(ni,moment,meansize)

do i=1,m
  if (abs(ni(i))<1.D-16) then
    nitemp(i) = 0.D0
  else
    nitemp(i) = ni(i)
  end if
end do
if (growth_function==4) then
  if (jet_cl_model>0) then
    ! Convert integer to string
    write(i_step_str, '(I0)') i_step
    open(99,file='pbe/psd_T'//trim(i_step_str)//'.out')    
  else
    open(99,file='pbe/psd.out')
  endif
else
  open(99,file='pbe/psd.out')
endif
do i=1,m
  write(99,1001) v_m(i),(6.D0/3.14159*v_m(i))**(1.D0/3.D0),nitemp(i), &
  & nitemp(i)*dv(i)/moment(0),v_m(i)*nitemp(i),v_m(i)*nitemp(i)*dv(i)/moment(1),dv(i)
end do
close(99)
if (i_writesp==1) then
  open(99,file='pbe/psd_sp.out')
  do i=1,m
    eta(i) = moment(0)*v_m(i)/moment(1)
    psi(i) = nitemp(i)*moment(1)/moment(0)**2
    write(99,1002) eta(i),psi(i)
  end do
end if

1001 format(7E20.10)
1002 format(2E20.10)

end subroutine pbe_output

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_deallocate()

!**********************************************************************************************
!
! Deallocate PBE variables
!
! By Stelios Rigopoulos
! Modified 06/05/2017
! Modified 04/07/2020
!
!**********************************************************************************************

use pbe_mod
 
deallocate(v,dv,v_m,nuc)

end subroutine pbe_deallocate

!**********************************************************************************************
