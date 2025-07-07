module ice_microphys_mod
  implicit none
  private
  public :: append_scalar, saturation_ratio

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

  !NOTE THE FOLLOWING ARE JUST PLACE-HOLDERS THAT MAY BE USEFUL TO UPDATE THE SUPERSATURATION DUE TO CONSUMPTION of WATER VAPOR
  subroutine saturation_ratio(Smw_time_series, Smw_time_derivative, k, dt, Loss_Sw, SAC_flag, AR_flag, FR_flag)
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
    use pbe_mod, only :current_temp, amb_p, RH, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model
    use thermo

    implicit none
    ! Input/output
    real, allocatable, intent(inout) :: Smw_time_series(:)
    real, allocatable, intent(inout) :: Smw_time_derivative(:)
    integer, intent(in) :: k
    real, intent(in) :: dt, Loss_Sw
    logical, intent(in) :: SAC_flag ! Check Schmidt and Appleman criterio
    logical, intent(in) :: AR_flag ! Check Activation-Relaxation point
    logical, intent(in) :: FR_flag ! Check Relaxation-Freezing point


    real :: Smw
    double precision :: p_sat_liq 

    if (jet_cl_model==1) then
      Smw = RH
    elseif (jet_cl_model==2) then
      call p_sat_liq_murphy_koop(p_sat_liq)
      Smw = amb_p * current_XH2O / p_sat_liq
    endif 

    ! Compute saturation ratio at this step
    if (k .le. 1) then

      ! Append new value for Smw time series
      call append_scalar(Smw_time_series, Smw)  
  
    ! Compute and append derivative if k >= 2
    else
      
      !Check if SAC satisfied
      if (SAC_flag) then
        
        if ((AR_flag) .or. (FR_flag)) then
          Smw_time_derivative(k) = Smw_time_derivative(k-1) - Loss_Sw

          ! Adjust next predicted scalar if applicable
          !if (k < 1000) then  ! or whatever max nt you have
          Smw_time_series(k+1) = Smw_time_series(k) + Smw_time_derivative(k-1)*dt
          !end if

          ! Append new value for Smw time series
          call append_scalar(Smw_time_series, Smw) 
          
          ! APpend new values for dSmw/dt time series
          call append_scalar(Smw_time_derivative, (Smw_time_series(k)-Smw_time_series(k-1))/dt)

        else
          ! If activation-relaxation is not satisfied, then simply append the data for saturation ratio from the model or from the LES dataset
          if (jet_cl_model==1) then
            Smw = RH
          elseif (jet_cl_model==2) then
            call p_sat_liq_murphy_koop(p_sat_liq)
            Smw = amb_p * current_XH2O / p_sat_liq
          endif  
    
          ! Append new value for Smw time series
          call append_scalar(Smw_time_series, Smw) 
          
          ! APpend new values for dSmw/dt time series
          call append_scalar(Smw_time_derivative, (Smw_time_series(k)-Smw_time_series(k-1))/dt)
        
        endif

      endif

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

subroutine pbe_droplet_growth(index, S_e, ni, g_coeff1,g_coeff2)

  !**********************************************************************************************
  !
  ! Computation of g_coeff1 from kinetic growth model. The model applies between activation-relaxation
  ! and freezing-relaxation point where the 
  !
  ! Based on Karcher et al. 2015 and Ponsonby et al. 2025
  !    
  ! Luca Boscagli 04/07/2025
  !
  !**********************************************************************************************
  
    use pbe_mod, only :v0, v_m, m !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
    use pbe_mod, only :current_temp, amb_p, RH, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model, kappa, Loss_Sw
    use thermo

    implicit none
  
    integer, intent(in)  :: index, S_e
    double precision, dimension(m), intent(in) :: ni
    double precision, intent(out)              :: g_coeff1, g_coeff2

    double precision :: p_water_sat_ice,p_water_sat_liq, S_v !,RH,amb_temp,amb_p,amb_rho
    double precision :: r_part,r_nuc
    double precision :: drdt
    double precision :: part_den !, part_den_l, part_den_r
    double precision :: Fd, Fk
    real,parameter :: pi = 3.141592653589793E+00
    real(kind=8) :: accom=1.0
    double precision :: gascon=8314.3
    

    !Compute droplet particle radius and nucleri (dry) radius from the volume
    r_part = (3.0 / (4.0 * pi) * v_m(index))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        ! radius of the nuclei (samllest particle volume) - this is constant

    !Compute equilibrium saturation ratio over the particle
    S_v = Seq(r_part, r_nuc, current_temp, kappa) ! this computes the supersaturation, so we need to add 1.0 to get the saturatiom
    S_v = S_v + 1.0

    ! Compute saturation ratio at this step - this is an input
    !call saturation_ratio(Smw_time_series, Smw_time_derivative, k, dt)
    !S_e = Smw_time_series(size(Smw_time_series))    

    ! Saturation presure over liquid
    call p_sat_liq_murphy_koop(p_water_sat_liq)

    !Compute growth coefficient G=1/(Fk+Fd) based on Rogers and Yau (1996)
    Fk = Lw*rho_w/(4 * ka_corr(current_temp, part_den, r_part) * current_temp) * (Lw*Mw/(gascon*current_temp) - 1.0) !thermodynamic term associated to heat conduction
    Fd = rho_w * gascon * current_temp / (4 * Mw * dv_corr(current_temp, r_part, amb_p, accom) * p_water_sat_liq) !term associated with vapor diffusion

    drdt = (1.0 / (Fk + Fd)) * (1 / r_part) * (S_e - S_v) 
  
    g_coeff2 = 2.0 / 3.0 
    if (drdt .ge. 0) then
       g_coeff1 = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
       Loss_Sw = Loss_Sw + amb_p * pi * rho_w / rho_air(current_temp, amb_p, amb_p/p_water_sat_liq) * (ni(index) * r_part**2 * drdt)
      else
       g_coeff1 = 0.0
    endif

    

  
  end subroutine pbe_droplet_growth  

  subroutine pbe_ice_growth(index, g_coeff1,g_coeff2)

    !**********************************************************************************************
    !
    ! Computation of g_coeff1 from kinetic growth model. The model applies fter freezing-relaxation point
    !
    ! Based on deposition model by Karcher et al. 2015 and Ponsonby et al. 2025
    !    
    ! Luca Boscagli 04/07/2025
    !
    !**********************************************************************************************
    
      use pbe_mod, only :v0, v_m !v0 = nuclei volume (named v_nuc in BOFFIN+PBE)
      use pbe_mod, only :current_temp, amb_p, RH, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model 
      
      implicit none
    
      integer, intent(in)  :: index
      double precision, intent(out)              :: g_coeff1, g_coeff2
  
      double precision :: p_water_sat_ice,p_water_sat_liq !,RH,amb_temp,amb_p,amb_rho
      double precision :: M_water !,M_air, X_water
      double precision :: r_part,r_nuc,den_ice
      double precision :: dif_water,lambda_water
      double precision :: coll_factor,Kn!,alpha_ice
      double precision :: fornow,drdt
      double precision :: part_den !, part_den_l, part_den_r
      real,parameter :: pi = 3.141592653589793E+00
      double precision :: gascon=8314.3
     
    !  M_water = 18.016 ! water molecular weigth
    !  M_air = 28.96 ! air molecular weigth
  
    
    !  r_part = (3.0 / (4.0 * pi) * v_m(index))**(1.0/3.0)  ! radius of spherically assumed particles computed from the volume at the middle of each bin v_m
    !  r_nuc = (3.0 / (4.0 * pi) * v0)**(1.0/3.0)        ! radius of the nuclei (samllest particle volume) - this is constant
    
    !  part_den = (part_den_l * r_nuc**3.0 + den_ice * &
    !                          (r_part**3.0 - r_nuc**3.0)) / r_part**3.0  ! particle density
    
    
    !  dif_water = 2.11D-5 * 101325.0 / amb_p * (current_temp / 273.15)**1.94  ! diffusion coefficient of water vapor molecules in air 
    ! lambda_water = 6.15D-8 * 101325.0 / amb_p * current_temp / 273.15       ! water vapor mean free path
    
    !  Kn = lambda_water / r_part  ! knudsen number 
    
    !  coll_factor = 1.0 / (1.0 / (1.0 + Kn) + 4.0 / 3.0 * Kn / alpha_ice) ! collision factor (G), accounts for transition from gas kinetic energy (G->1 for Kn->0) to continuum regime (G->0 for Kn->1)
    
    !  fornow = dif_water * coll_factor * M_water * (p_water - p_water_sat_ice) &
    !          / (gascon * current_temp)   ! nominator of dr/dt 
                                      ! [ M_water * p_water_sat_ice / (gascon * current_temp)]: saturated (relative to ice) water vapor density (considering compressibility factor approx = 1)
                                      ! This is calling gascon = 8314.3 J/K/kg (in module_chemistry.f90)
                                      ! gascon/M_water is the specific gas constant for water vapor, i.e., Rv = 461.52 J/K/kg
                                      
    
    !  drdt = fornow / (part_den * r_part) ! change in particle radius over time
    
    !  g_coeff2 = 2.0 / 3.0 
    !  if (p_water.ge.p_water_sat_liq) then 
    !    g_coeff1 = 4.0 * pi * (3.0 / (4.0 * pi))**g_coeff2 * drdt ! Equivalent to  3 * (4/3 pi)**(1/3) * dr/dt
    !  else
    !    g_coeff1 = 0.0
    !  endif  
      
    
    end subroutine pbe_ice_growth  

end module ice_microphys_mod

