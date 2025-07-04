module ice_microphys_mod
  implicit none
  private
  public :: update_time_series, append_scalar, some_scalar_function

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
  subroutine saturation_ratio(Smw_time_series, Smw_time_derivative, k, dt)

    use pbe_mod, only :current_temp, amb_p, RH, part_den_l, alpha_ice, p_water, current_XH2O, jet_cl_model
    
    implicit none
    ! Input/output
    real, allocatable, intent(inout) :: Smw_time_series(:)
    real, allocatable, intent(inout) :: Smw_time_derivative(:)
    integer, intent(in) :: k
    real, intent(in) :: dt
  
    real :: Smw, p_sat_liq 
    logical :: update_condition
  
    ! Compute saturation ratio at this step
    if (jet_cl_model==1) then
      Smw = RH
    elseif (jet_cl_model==2) then
      call p_sat_liq_murphy_koop(p_sat_liq)
      Smw = amb_p * current_XH2O / p_sat_liq
    endif  
  
    ! Append new value
    call append_scalar(Smw_time_series, Smw)
  
    ! Compute and append derivative if k >= 2
    if (k > 1) then
        call append_scalar(Smw_time_derivative, (Smw_time_series(k)-Smw_time_series(k-1))/dt)
  
        ! Check condition
        update_condition = (Smw_time_derivative(k-1) > 0.5) !0.5 is dummy value, we should check is consumption of saturation based on K15 is greater than 0 and subtract if so
  
        if (update_condition) then
          ! Adjust the most recent derivative
          Smw_time_derivative(k-1) = Smw_time_derivative(k-1) * 0.9 !0.9 is dummy value, we should subtract the consumption of saturation based on K15
  
          ! Adjust next predicted scalar if applicable
          !if (k < 1000) then  ! or whatever max nt you have
          Smw_time_series(k+1) = Smw_time_series(k) + Smw_time_derivative(k-1)*dt
          !end if
        end if
    end if
  
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
  
  ! Dummy placeholder for your model's scalar value
  subroutine saturation_ratio(step)
    implicit none
    
    integer, intent(in) :: step
    
    ! water vapor partial pressure
    if (jet_cl_model==1) then
      p_water = p_water_sat_liq * RH
    elseif (jet_cl_model==2) then
      p_water = amb_p * current_XH2O
    endif  
  end subroutine some_scalar_function


subroutine pbe_droplet_growth(index, g_coeff1,g_coeff2)

  !**********************************************************************************************
  !
  ! Computation of g_coeff1 from kinetic growth model. The model applies between activation-relaxation
  ! and freezing-relaxation point where the 
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
    
  
  end subroutine pbe_droplet_growth  

end module ice_microphys_mod

