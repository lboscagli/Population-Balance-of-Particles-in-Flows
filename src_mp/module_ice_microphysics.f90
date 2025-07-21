module ice_microphys_mod
  implicit none
  private
  public :: append_scalar, saturation_ratio, p_sat_ice_murphy_koop, p_sat_liq_murphy_koop
  public :: pbe_condensational_droplet_growth_Bier, pbe_depositional_growth_ice_Bier, pbe_depositional_growth_ice, pbe_freezing_temperature
contains

!---------------------------------------------------------------------------------------------------
! Auxiliary functions for extended K15 model implementation
!---------------------------------------------------------------------------------------------------

  subroutine p_sat_ice_murphy_koop(p_sat_ice)

  !**********************************************************************************************
  !
  ! Determines the saturation vapor pressure of ice as a function of temperature using the Murphy-Koop equation 
  !
  ! Murphy, D. M., & Koop, T. (2005)
  !
  !**********************************************************************************************
    use pbe_mod, only :current_temp
  
    implicit none

    double precision, intent(out)                  :: p_sat_ice
    
    !----------------------------------------------------------------------------------------------
    
  
    p_sat_ice = EXP(9.550426 - 5723.265/current_temp + 3.53068 * LOG(current_temp) -  0.00728332 * current_temp )


  end subroutine p_sat_ice_murphy_koop


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

  subroutine saturation_ratio(Smw_time_series, k, dt, Smw, Loss_Sw)
  !**********************************************************************************************
  !
  ! This function is to compute the consumption of supersaturation based either on water droplet 
  ! or ice crystal growth. The saturation ratio (S) is appended to an array such that at every step
  ! after the first one (k>1) the derivative of S wrt time can be computed and S can be updated at
  ! the next step
  !
  ! Based on Karcher et al. 2015 and Ponsonby et al. 2025
  !    
  ! Luca Boscagli 07/07/2025
  !
  !**********************************************************************************************
    use pbe_mod, only :current_temp, amb_p, G_mixing_line, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model
    use thermo

    implicit none
    ! Input/output
    real, allocatable, intent(inout) :: Smw_time_series(:)
    !real, intent(in) :: Smw_time_series(:)
    !real, allocatable, intent(inout) :: Smw_time_derivative(:)
    integer, intent(in) :: k
    real, intent(in) :: Smw, Loss_Sw
    double precision, intent(in) :: dt
    real :: dSmwdt, Smw_new 
    !double precision :: p_sat_liq
    
    if (k .le. 1) then
      ! Append new value for Smw time series
      call append_scalar(Smw_time_series, Smw) 

    elseif (k > 1) then
      !Compute saturation ratio derivative wrt time and subtract supersaturation consumption (note that if there is no growth then Loss_Sw will be zero)
      dSmwdt = (Smw-Smw_time_series(k-1))/dt - Loss_Sw

      ! Adjust next predicted scalar if applicable
      Smw_new = Smw + dSmwdt*dt

      ! Append new value for Smw time series
      call append_scalar(Smw_time_series, Smw_new) 
      
    endif
  
  end subroutine saturation_ratio

! Append a new scalar to a 1D allocatable array
  subroutine append_scalar(arr, val)
    implicit none
    real, allocatable, intent(inout) :: arr(:)
    real, intent(in) :: val
    real, allocatable :: tmp(:)
  
    allocate(tmp(size(arr)+1))
    if (size(arr) > 0) tmp(1:size(arr)) = arr
    tmp(size(tmp)) = val
    call move_alloc(tmp, arr)
  end subroutine append_scalar

  subroutine pbe_condensational_droplet_growth_Bier(index, ni, g_coeff1_l,g_coeff1_r,g_coeff2)

  !**********************************************************************************************
  !
  ! Computation of g_coeff1 from kinetic growth model. The model applies between activation-relaxation
  ! and freezing-relaxation point where the 
  !
  ! Based on Bier et al. 2021 - Eq. 2
  !    
  ! Luca Boscagli 14/07/2025
  !
  !**********************************************************************************************
  
    use pbe_mod, only :v0, v_m, m, dv, v !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
    use pbe_mod, only :current_temp, amb_temp, amb_p, G_mixing_line, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model, kappa, Loss_Sw, current_rho
    use pbe_mod, only :Smw_time_series, step_update, r_vc, S_vc
    use thermo

    implicit none
  
    integer, intent(in)  :: index
    double precision, dimension(m), intent(in) :: ni
    double precision, intent(out)              :: g_coeff1_l,g_coeff1_r, g_coeff2

    double precision :: p_water_sat_ice,p_water_sat_liq, S_v, S_e !,RH,amb_temp,amb_p,amb_rho
    double precision :: r_part,r_nuc
    double precision :: drdt, dmdt
    double precision :: F_M, F_H, lambda_v, Kn_v, beta_M, lambda_a, Kn_a, beta_H
    double precision :: part_den !, part_den_l, part_den_r
    double precision :: Fd, Fk
    real,parameter :: pi = 3.141592653589793E+00
    real(kind=8) :: accom=1.0, r_crit, s_crit
    double precision :: gascon=8314.3
    double precision :: M_water = 18.016 ! molar weight of water [kg/mol]
    double precision :: M_air = 28.96 ! molar weight of dry air [kg/mol]
    
    !Right side

    !Compute droplet particle radius and nucleri (dry) radius from the volume
    r_part = (3.0 / (4.0 * pi) * v(index))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        ! radius of the nuclei (samllest particle volume) - this is constant
 
    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + rho_w * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density

    !Compute equilibrium saturation ratio over the particle
    if (r_part .eq. r_nuc) then
      S_v = S_vc
    else  
      S_v = Seq(r_part, r_nuc, current_temp, kappa) ! this computes the supersaturation, so we need to add 1.0 to get the saturatiom
      S_v = S_v + 1.0   
    endif

    ! Saturation presure over liquid
    call p_sat_liq_murphy_koop(p_water_sat_liq)

    
    !Mass (M) diffusion term
    F_M = (gascon/M_water)*current_temp/dv_cont(current_temp, amb_p)
    !Heat (H) diffusion term
    F_H = (S_v*p_water_sat_liq) * (latent_heat_cond_evap(current_temp))**2 / (gascon/M_water) / kair_conductivity(current_temp) / current_temp**2
    !correction for mass diffusion term to account for non continuum effects (large knudsen)
    lambda_v = 2 * dv_cont(current_temp, amb_p) * SQRT(1.0 / (2 * (gascon/M_water) * current_temp))
    Kn_v = lambda_v / r_part
    beta_M = (Kn_v + 1.0) / ((4.0/3.0/1.0 + 0.377) * Kn_v  + (4.0/3.0/1.0) * Kn_v**2)
    !correction for heat diffusion term to account for non continuum effects (large knudsen)
    lambda_a = 0.8 * kair_conductivity(current_temp) * current_temp / amb_p * SQRT(1.0 / (2 * (gascon/M_air) * current_temp)) 
    Kn_a = lambda_a / r_part    
    beta_H = (Kn_a + 1.0) / ((4.0/3.0/1.0 + 0.377) * Kn_a  + (4.0/3.0/1.0) * Kn_a**2)
    
    !Compute droplet mass growth 
    dmdt = (4*pi*r_part) * (Smw_time_series(step_update)*p_water_sat_liq - S_v*p_water_sat_liq) / (F_m/beta_M + F_H*beta_H)

    !Compute droplet radial growth 
    drdt = dmdt / (4 * pi * part_den * r_part**2)

    !write(*,*) 'Condensational growth'
    !write(*,*) 'drdt',drdt

    g_coeff2 = 2.0 / 3.0 
    if (drdt .ge. 0) then
       g_coeff1_r = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
       Loss_Sw = Loss_Sw + (amb_p/p_water_sat_liq/epsilon_fluid) * pi * part_den / current_rho * (ni(index)*dv(index) * r_part**2 * drdt)
    else
       g_coeff1_r = 0.0
    endif

    !Left side

    !Compute droplet particle radius from the volume
    r_part = (3.0 / (4.0 * pi) * v(index-1))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
  
    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + rho_w * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density

    !Compute equilibrium saturation ratio over the particle
    if (r_part .eq. r_nuc) then
      S_v = S_vc
    else  
      S_v = Seq(r_part, r_nuc, current_temp, kappa) ! this computes the supersaturation, so we need to add 1.0 to get the saturatiom
      S_v = S_v + 1.0   
    endif

    
    !Mass (M) diffusion term
    F_M = (gascon/M_water)*current_temp/dv_cont(current_temp, amb_p)
    !Heat (H) diffusion term
    F_H = (S_v*p_water_sat_liq) * (latent_heat_cond_evap(current_temp))**2 / (gascon/M_water) / kair_conductivity(current_temp) / current_temp**2
    !correction for mass diffusion term to account for non continuum effects (large knudsen)
    lambda_v = 2 * dv_cont(current_temp, amb_p) * SQRT(1.0 / (2 * (gascon/M_water) * current_temp))
    Kn_v = lambda_v / r_part
    beta_M = (Kn_v + 1.0) / ((4.0/3.0/1.0 + 0.377) * Kn_v  + (4.0/3.0/1.0) * Kn_v**2)
    !correction for heat diffusion term to account for non continuum effects (large knudsen)
    lambda_a = 0.8 * kair_conductivity(current_temp) * current_temp / amb_p * SQRT(1.0 / (2 * (gascon/M_air) * current_temp)) 
    Kn_a = lambda_a / r_part    
    beta_H = (Kn_a + 1.0) / ((4.0/3.0/1.0 + 0.377) * Kn_a  + (4.0/3.0/1.0) * Kn_a**2)
    
    !Compute droplet mass growth 
    dmdt = (4*pi*r_part) * (Smw_time_series(step_update)*p_water_sat_liq - S_v*p_water_sat_liq) / (F_m/beta_M + F_H*beta_H)

    !Compute droplet radial growth 
    drdt = dmdt / (4 * pi * part_den * r_part**2)

    !write(*,*) 'Condensational growth'
    !write(*,*) 'drdt',drdt

    if (drdt .ge. 0) then
       g_coeff1_l = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
    else
       g_coeff1_l = 0.0
    endif
    
  
  end subroutine pbe_condensational_droplet_growth_Bier   

  subroutine pbe_depositional_growth_ice_Bier(index, ni, g_coeff1_l,g_coeff1_r,g_coeff2)

  !**********************************************************************************************
  !
  ! Computation of g_coeff1 from kinetic growth model. The model applies after
  ! freezing-relaxation point where the 
  !
  ! Based on Bier et al. 2021 - Eq. 7
  !    
  ! Luca Boscagli 14/07/2025
  !
  !**********************************************************************************************
  
    use pbe_mod, only :v0, v_m, m, dv, v !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
    use pbe_mod, only :current_temp, amb_temp, amb_p, G_mixing_line, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model, kappa, Loss_Sw, current_rho
    use pbe_mod, only :Smw_time_series, step_update, r_vc, S_vc
    use thermo

    implicit none
  
    integer, intent(in)  :: index
    double precision, dimension(m), intent(in) :: ni
    double precision, intent(out)              :: g_coeff1_l,g_coeff1_r, g_coeff2

    double precision :: p_water_sat_ice,p_water_sat_liq, S_v !,RH,amb_temp,amb_p,amb_rho
    double precision :: r_part,r_nuc,den_ice
    double precision :: drdt, dmdt, denom_dmdt
    double precision :: lambda_v, lambda_a, beta_v, beta_k, rho_air_dry
    double precision :: part_den !, part_den_l, part_den_r
    double precision :: Fd, Fk
    real,parameter :: pi = 3.141592653589793E+00
    real(kind=8) :: accom=1.0, r_crit, s_crit
    double precision :: gascon=8314.3
    double precision :: M_water = 18.016 ! molar weight of water [kg/mol]
    double precision :: M_air = 28.96 ! molar weight of dry air [kg/mol]
    !double precision :: Cp = 1004.0 !! specific heat of dry air at const. pressure [J/kg/K]
    
    !Right side

    !Compute droplet particle radius and nucleri (dry) radius from the volume
    r_part = (3.0 / (4.0 * pi) * v(index))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        ! radius of the nuclei (samllest particle volume) - this is constant

    ! Constant ice particle density 
    den_ice = 917.0  
    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + den_ice * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density

    !Compute equilibrium saturation ratio over the particle
    if (r_part .eq. r_nuc) then
      S_v = S_vc
    else  
      S_v = Seq(r_part, r_nuc, current_temp, kappa) ! this computes the supersaturation, so we need to add 1.0 to get the saturatiom
      S_v = S_v + 1.0   
    endif

    ! Saturation presure over ice and liquid
    call p_sat_ice_murphy_koop(p_water_sat_ice)
    call p_sat_liq_murphy_koop(p_water_sat_liq)

    !correction for mass diffusion term to account for non continuum effects (large knudsen)
    lambda_v = 2 * dv_cont(current_temp, amb_p) * SQRT(1.0 / (2 * (gascon/M_water) * current_temp))
    beta_v = r_part / (r_part + lambda_v) + (4 * 1.0 * dv_cont(current_temp, amb_p))/(0.5 * r_part * mean_thermal_speed(current_temp))  
    !correction for heat diffusion term to account for non continuum effects (large knudsen)
    lambda_a = 0.8 * kair_conductivity(current_temp) * current_temp / amb_p * SQRT(1.0 / (2 * (gascon/M_air) * current_temp)) 
    rho_air_dry = amb_p / current_temp / (gascon/M_air)
    beta_k = r_part / (r_part + lambda_a) + (4 * kair_conductivity(current_temp))/(0.7 * Cp * mean_thermal_speed(current_temp) * rho_air_dry)
      

    !Compute droplet mass growth - spherical assumption (C=1.0)
    denom_dmdt = ((dv_cont(current_temp, amb_p)*(beta_v**(-1))*latent_heat_dep_sub(current_temp)*S_v*p_water_sat_ice) / (kair_conductivity(current_temp) * (beta_k**(-1)) * (beta_v**(-1)) * current_temp))  * (latent_heat_dep_sub(current_temp) / ((gascon/M_water) * current_temp) - 1.0) + (gascon/M_water) * current_temp
    dmdt = (4*pi*r_part * 1.0) * dv_cont(current_temp, amb_p) * (1.0/beta_v) * (Smw_time_series(step_update)*p_water_sat_liq - S_v*p_water_sat_ice) / denom_dmdt

    !Compute droplet radial growth 
    drdt = dmdt / (4 * pi * part_den * r_part**2)
    !write(*,*) 'Depositional growth'
    !write(*,*) 'drdt',drdt

    g_coeff2 = 2.0 / 3.0 
    if (drdt .ge. 0) then
       g_coeff1_r = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
       Loss_Sw = Loss_Sw + (amb_p/p_water_sat_liq/epsilon_fluid) * pi * part_den / current_rho * (ni(index)*dv(index) * r_part**2 * drdt)
    else
       g_coeff1_r = 0.0
    endif

    !Left side

    !Compute droplet particle radius  from the volume
    r_part = (3.0 / (4.0 * pi) * v(index-1))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m

    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + den_ice * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density

    !Compute equilibrium saturation ratio over the particle
    if (r_part .eq. r_nuc) then
      S_v = S_vc
    else  
      S_v = Seq(r_part, r_nuc, current_temp, kappa) ! this computes the supersaturation, so we need to add 1.0 to get the saturatiom
      S_v = S_v + 1.0   
    endif


    !correction for mass diffusion term to account for non continuum effects (large knudsen)
    lambda_v = 2 * dv_cont(current_temp, amb_p) * SQRT(1.0 / (2 * (gascon/M_water) * current_temp))
    beta_v = r_part / (r_part + lambda_v) + (4 * 1.0 * dv_cont(current_temp, amb_p))/(0.5 * r_part * mean_thermal_speed(current_temp))  
    !correction for heat diffusion term to account for non continuum effects (large knudsen)
    lambda_a = 0.8 * kair_conductivity(current_temp) * current_temp / amb_p * SQRT(1.0 / (2 * (gascon/M_air) * current_temp)) 
    rho_air_dry = amb_p / current_temp / (gascon/M_air)
    beta_k = r_part / (r_part + lambda_a) + (4 * kair_conductivity(current_temp))/(0.7 * Cp * mean_thermal_speed(current_temp) * rho_air_dry)
      

    !Compute droplet mass growth - spherical assumption (C=1.0)
    denom_dmdt = ((dv_cont(current_temp, amb_p)*(beta_v**(-1))*latent_heat_dep_sub(current_temp)*S_v*p_water_sat_ice) / (kair_conductivity(current_temp) * (beta_k**(-1)) * (beta_v**(-1)) * current_temp))  * (latent_heat_dep_sub(current_temp) / ((gascon/M_water) * current_temp) - 1.0) + (gascon/M_water) * current_temp
    dmdt = (4*pi*r_part * 1.0) * dv_cont(current_temp, amb_p) * (1.0/beta_v) * (Smw_time_series(step_update)*p_water_sat_liq - S_v*p_water_sat_ice) / denom_dmdt

    !Compute droplet radial growth 
    drdt = dmdt / (4 * pi * part_den * r_part**2)
    !write(*,*) 'Depositional growth'
    !write(*,*) 'drdt',drdt

    if (drdt .ge. 0) then
       g_coeff1_l = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
    else
       g_coeff1_l = 0.0
    endif    
  
  end subroutine pbe_depositional_growth_ice_Bier  

  subroutine pbe_depositional_growth_ice(index, ni, g_coeff1_l,g_coeff1_r,g_coeff2)

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
    use pbe_mod, only :v0, v_m, m, dv, v !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
    use pbe_mod, only :current_temp, amb_p, G_mixing_line, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model, Loss_Sw, current_rho  
    use thermo

    implicit none
  
  !class(pbe_growth) :: this
    integer, intent(in)  :: index
    double precision, dimension(m), intent(in) :: ni
    double precision, intent(out)              :: g_coeff1_l,g_coeff1_r, g_coeff2
    !integer :: isp
    double precision :: p_water_sat_ice,p_water_sat_liq !,RH,amb_temp,amb_p,amb_rho
    double precision :: M_water !,M_air, X_water
    double precision :: r_part,r_nuc,den_ice
    double precision :: dif_water,lambda_water
    double precision :: coll_factor,Kn!,alpha_ice
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
  
    ! X_water =  water molecular fraction
    ! RH = Relative humidity
    M_water = 18.016 ! water molecular weigth
    ! M_air = 28.96 ! air molecular weigth
    ! amb_temp = ambient temperature in kelvin
    ! amb_p = ambient pressure in Pascal
    ! amb_rho = amb_p/(gascon/M_air)/amb_temp ! ambient air density
    ! part_den_l = 1550.0				! Density of particles on the left side of the PSD (kg/m^3)
    ! v0 = v_nuc ! Nuclei volume  
    
    ! saturated (relative to liquid) water vapour partial pressure 
    call p_sat_liq_murphy_koop(p_water_sat_liq)

    ! saturated (relative to ice) water vapour partial pressure 
    call p_sat_ice_murphy_koop(p_water_sat_ice)
     
    ! Constant ice particle density 
    den_ice = 917.0  
  
    ! radius of the nuclei (samllest particle volume) - this is constant  
    r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        

    ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    r_part = (3.0 / (4.0 * pi) * v(index))**(1.0/3.0)    
    
    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + den_ice * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density
  
    ! diffusion coefficient of water vapor molecules in air
    dif_water = 2.11D-5 * 101325.0 / amb_p * (current_temp / 273.15)**1.94   
    ! water vapor mean free path
    lambda_water = 6.15D-8 * 101325.0 / amb_p * current_temp / 273.15       
  
    ! knudsen number
    Kn = lambda_water / r_part   
    ! deposition coefficient. Higher values of alpha speed up the deposition rate - user input in CPMOD
    ! alpha_ice = 0.1 
  
    ! collision factor (G), accounts for transition from gas kinetic energy (G->1 for Kn->0) to continuum regime (G->0 for Kn->1)
    coll_factor = 1.0 / (1.0 / (1.0 + Kn) + 4.0 / 3.0 * Kn / alpha_ice) 
  
    ! nominator of dr/dt 
    ! [ M_water * p_water_sat_ice / (gascon * current_temp)]: saturated (relative to ice) water vapor density 
    ! (considering compressibility factor approx = 1)
    ! This is calling gascon = 8314.3 J/K/kg (in module_chemistry.f90)
    ! gascon/M_water is the specific gas constant for water vapor, i.e., Rv = 461.52 J/K/kg    

    fornow = dif_water * coll_factor * M_water * (p_water - p_water_sat_ice) &
            / (gascon * current_temp)   

                                    
    ! Change in particle radius over time
    drdt = fornow / (part_den * r_part) 
  
    ! Compute coefficients needed for growth and supersaturation consumption
    g_coeff2 = 2.0 / 3.0 
    if (drdt .ge. 0) then
      g_coeff1_r = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
      Loss_Sw = Loss_Sw + (amb_p/p_water_sat_liq/epsilon_fluid) * pi * part_den / current_rho * (ni(index)*dv(index) * r_part**2 * drdt)
    else
       g_coeff1_r = 0.0
    endif    

    !Left side
    ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    r_part = (3.0 / (4.0 * pi) * v(index-1))**(1.0/3.0)    
    
    ! particle density
    part_den = (part_den_l * r_nuc**3.0 + den_ice * &
                            (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density
  
    ! knudsen number
    Kn = lambda_water / r_part   
    ! deposition coefficient. Higher values of alpha speed up the deposition rate - user input in CPMOD
    ! alpha_ice = 0.1 
  
    ! collision factor (G), accounts for transition from gas kinetic energy (G->1 for Kn->0) to continuum regime (G->0 for Kn->1)
    coll_factor = 1.0 / (1.0 / (1.0 + Kn) + 4.0 / 3.0 * Kn / alpha_ice) 
  
    ! nominator of dr/dt 
    ! [ M_water * p_water_sat_ice / (gascon * current_temp)]: saturated (relative to ice) water vapor density 
    ! (considering compressibility factor approx = 1)
    ! This is calling gascon = 8314.3 J/K/kg (in module_chemistry.f90)
    ! gascon/M_water is the specific gas constant for water vapor, i.e., Rv = 461.52 J/K/kg    

    fornow = dif_water * coll_factor * M_water * (p_water - p_water_sat_ice) &
            / (gascon * current_temp)   

                                    
    ! Change in particle radius over time
    drdt = fornow / (part_den * r_part)     

    if (drdt .ge. 0) then
      g_coeff1_l = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
    else
      g_coeff1_l = 0.0
    endif    
  !! Luca: the section below is not needed for CPMOD standalone (no coupling with BOFFIN)
  !  dmdt = 4.0 * pi * r_part * fornow
  
  !  m_source(ijk) = m_source(ijk) + dmdt * this%ni_part_pbe(index) &
  !                * this%dv_part(index) * ajc(ijk)                       ! [kg/s] !ajc is the jacobian of the coordinate transformation (BOFFIN manual p.16 and p.41)
                                                                                  !Presumably if the computational grid is already orthogonal, then this is an identity matrix?
  
    
  end subroutine pbe_depositional_growth_ice

  subroutine pbe_freezing_temperature(index, T_frz)

  !**********************************************************************************************
  !
  ! Computation of freezing temperature. This forms the criterion to witch from condensational
  ! droplet growth to depositional ice growth.
  !
  ! Based on Karcher et al. 2015, Bier et al. 2021 and Ponsonby et al. 2025
  !    
  ! Luca Boscagli 14/07/2025
  !
  !**********************************************************************************************

    use pbe_mod, only :v0, v_m, m, dv, v !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
    use pbe_mod, only :current_temp, plume_cooling_rate
    use thermo

    implicit none
  
    integer, intent(in)  :: index
    double precision, intent(out)  :: T_frz

    double precision :: LWV, J_freez_rate_coeff
    real,parameter :: pi = 3.141592653589793E+00
    double precision :: gascon=8314.3
    double precision :: a_1=-3.5714 !(1/K) - Karcher et al
    double precision :: a_2=858.719
    

    !Compute liquid water volume (LWV)
    LWV = v_m(index) - v0

    !COmpute freezing rate coefficient
    J_freez_rate_coeff = 1E6 * exp(a_1*current_temp + a_2)

    !COmpute freezing temperature
    if (LWV .eq. 0.) then
      T_frz = 0.0
    else  
      T_frz = (1.0/a_1) * (LOG((1.0E-6 * a_1 * plume_cooling_rate)/LWV) - a_2) ! Expression based on eq. 36 in (Karcher et al. 2015) and eq. 6 in (Bier et al. 2021)
    endif
    
  end subroutine pbe_freezing_temperature

end module ice_microphys_mod

