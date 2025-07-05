module constants
  !! Physical and model constants.
  implicit none
  public
  real(kind=8), parameter :: g = 9.81_8     !! gravitational acceleration [m/s^2]
  real(kind=8), parameter :: Cp = 1004.0_8  !! specific heat of dry air at const. pressure [J/kg/K]
  real(kind=8), parameter :: L = 2.25e6_8   !! latent heat of vaporization [J/kg]
  real(kind=8), parameter :: rho_w = 1000.0_8 !! density of water [kg/m^3]
  real(kind=8), parameter :: R = 8.314_8    !! universal gas constant [J/mol/K]
  real(kind=8), parameter :: Mw = 18.016_8/1.0e3_8 !! molar weight of water [kg/mol]
  real(kind=8), parameter :: Ma = 28.96_8/1.0e3_8 !! molar weight of dry air [kg/mol]
  real(kind=8), parameter :: Rd = R/Ma      !! gas constant for dry air [J/(kg/K)]
  real(kind=8), parameter :: Rv = R/Mw      !! gas constant for water vapor [J/(kg/K)]
  real(kind=8), parameter :: Dv = 3.0e-5_8  !! water vapor diffusivity in air [m^2/s]
  real(kind=8), parameter :: ac = 1.0_8     !! condensation coefficient [-]
  real(kind=8), parameter :: Ka = 2.0e-2_8  !! thermal conductivity of air [J/m/s/K]
  real(kind=8), parameter :: at = 0.96_8    !! thermal accommodation coefficient [-]
  real(kind=8), parameter :: epsilon = 0.622_8 !! Mw/Ma [-]

  integer, parameter :: N_STATE_VARS = 7
end module constants


module optimizer
  !! Provides a simple golden-section minimizer on [a,b].
  implicit none
  private
  public :: fminbound
contains
  function fminbound(f, a, b, tol) result(x_min)
    !! Find the minimizer of the function f(x) on the interval [a,b].
    !! Uses the golden-section search algorithm.
    !!
    !! Arguments:
    !!   f   : external real function of one variable
    !!   a,b : real, the lower and upper bounds
    !!   tol : real, the convergence tolerance
    !!
    !! Returns:
    !!   x_min : real, approximate minimizer of f on [a,b].
    implicit none
    interface
      function f(x) result(fval)
        real(kind=8), intent(in) :: x
        real(kind=8) :: fval
      end function f
    end interface
    real(kind=8), intent(in) :: a, b, tol
    real(kind=8) :: x_min
    real(kind=8) :: gr, c, d, fc, fd, left, right

    gr = (sqrt(5.0_8)-1.0_8)/2.0_8
    left = a; right = b
    c = right - (right-left)*gr
    d = left + (right-left)*gr
    fc = f(c); fd = f(d)

    do while (abs(right-left) > tol)
      if (fc < fd) then
        right = d
        d = c
        fd = fc
        c = right - (right-left)*gr
        fc = f(c)
      else
        left = c
        c = d
        fc = fd
        d = left + (right-left)*gr
        fd = f(d)
      end if
    end do
    x_min = (left+right)/2.0_8
  end function fminbound
end module optimizer


module thermo
  !! Thermodynamics & microphysics functions.
  use constants
  use optimizer
  implicit none
  public :: dv_cont, dv_corr, rho_air, es, ka_cont, ka_corr, sigma_w
  public :: Seq, Seq_approx, kohler_crit, critical_curve, r_eff

  ! Temporary variables for minimization:
  real(kind=8), save :: tmp_r_dry, tmp_T, tmp_kappa

contains

  function dv_cont(T, P) result(Dv_c)
    !! Compute water vapor diffusivity (continuum regime) at temperature T and pressure P.
    !! T [K], P [Pa]; returns Dv_c [m^2/s].
    real(kind=8), intent(in) :: T, P
    real(kind=8) :: Dv_c, P_atm
    P_atm = P * 1.01325e-5_8
    Dv_c = 1.0e-4_8*(0.211_8/P_atm)*( (T/273.0_8)**1.94_8 )
  end function dv_cont

  function dv_corr(T, r, P, accom) result(Dv_nc)
    !! Compute water vapor diffusivity corrected for non-continuum effects.
    !! T [K], r [m], P [Pa], accom: condensation coefficient.
    !! Returns Dv_nc [m^2/s].
    real(kind=8), intent(in) :: T, r, P, accom
    real(kind=8) :: Dv_nc, Dv_c, denom
    Dv_c = dv_cont(T, P)
    denom = 1.0_8+(Dv_c/(accom*r))*sqrt((2.0_8*3.141592653589793_8*Mw)/(R*T))
    Dv_nc = Dv_c/denom
  end function dv_corr

  function es(T_c) result(e_sat)
    !! Compute saturation vapor pressure over water at temperature T_c [C].
    !! Returns e_sat [Pa].
    real(kind=8), intent(in) :: T_c
    real(kind=8) :: e_sat
    e_sat = 611.2_8*exp(17.67_8*T_c/(T_c+243.5_8))
  end function es

  function rho_air(T, P, RH) result(rho_a)
    !! Compute moist air density at temperature T [K], pressure P [Pa], relative humidity RH [-].
    !! Uses ideal gas law with virtual temperature.
    real(kind=8), intent(in) :: T, P, RH
    real(kind=8) :: rho_a, qsat, Tv
    qsat = RH*0.622_8*(es(T-273.15_8)/P)
    Tv = T*(1.0_8+0.61_8*qsat)
    rho_a = P/(Rd*Tv)
  end function rho_air

  function ka_cont(T) result(ka_c)
    !! Compute thermal conductivity of air (continuum), T [K].
    !! Returns ka_c [J/m/s/K].
    real(kind=8), intent(in) :: T
    real(kind=8) :: ka_c
    ka_c = 1.0e-3_8*(4.39_8+0.071_8*T)
  end function ka_cont

  function ka_corr(T, rho, r) result(ka_nc)
    !! Compute thermal conductivity of air corrected for non-continuum effects.
    !! T [K], rho [kg/m^3], r [m].
    !! Returns ka_nc [J/m/s/K].
    real(kind=8), intent(in) :: T, rho, r
    real(kind=8) :: ka_nc, ka_c, denom
    ka_c = ka_cont(T)
    denom = 1.0_8+(ka_c/(at*r*rho*Cp))*sqrt((2.0_8*3.141592653589793_8*Ma)/(R*T))
    ka_nc = ka_c/denom
  end function ka_corr

  function sigma_w(T) result(sigma)
    !! Compute surface tension of water at temperature T [K].
    !! Returns sigma [J/m^2].
    real(kind=8), intent(in) :: T
    real(kind=8) :: sigma
    sigma = 0.0761_8-1.55e-4_8*(T-273.15_8)
  end function sigma_w

  function Seq(r, r_dry, T, kappa) result(S_eq)
    !! Compute equilibrium supersaturation over a particle per κ-Köhler theory.
    !! r,r_dry [m]; T [K]; kappa [-].
    !! Returns S_eq [-].
    real(kind=8), intent(in) :: r, r_dry, T, kappa
    real(kind=8) :: S_eq, A, B
    A = (2.0_8*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (r**3-r_dry**3)/(r**3-(r_dry**3)*(1.0_8-kappa))
    S_eq = exp(A)*B - 1.0_8
  end function Seq

  function Seq_approx(r, r_dry, T, kappa) result(S_eq_approx)
    !! Compute approximate equilibrium supersaturation using Taylor-expansion.
    !! r,r_dry [m]; T [K]; kappa [-].
    !! Returns S_eq_approx [-].
    real(kind=8), intent(in) :: r, r_dry, T, kappa
    real(kind=8) :: S_eq_approx, A
    A = (2.0_8*Mw*sigma_w(T))/(R*T*rho_w*r)
    S_eq_approx = A - kappa*(r_dry**3)/(r**3)
  end function Seq_approx

  function neg_Seq_fixed(r) result(fval)
    real(kind=8), intent(in) :: r
    real(kind=8) :: fval
    fval = -1.0_8*Seq(r, tmp_r_dry, tmp_T, tmp_kappa)
  end function neg_Seq_fixed

  subroutine kohler_crit(T, r_dry, kappa, approx, r_crit, s_crit)
    !! Compute Köhler critical radius and supersaturation for dry particle.
    !! T [K], r_dry [m], kappa [-].
    !! approx: logical - use approximate eqn.
    !! r_crit, s_crit are returned.    
    real(kind=8), intent(in) :: T, r_dry, kappa
    logical, intent(in) :: approx
    real(kind=8), intent(out) :: r_crit, s_crit
    real(kind=8) :: A

    if (approx) then
      A = (2.0_8*Mw*sigma_w(T))/(R*T*rho_w)
      s_crit = sqrt((4.0_8*A**3)/(27.0_8*kappa*(r_dry**3)))
      r_crit = sqrt((3.0_8*kappa*(r_dry**3))/A)
    else
      ! Set the module variables for the function
      tmp_r_dry = r_dry
      tmp_T = T
      tmp_kappa = kappa
      ! Call minimizer on neg_Seq_fixed
      r_crit = fminbound(neg_Seq_fixed, r_dry, r_dry*1.0e4_8, 1.0e-10_8)
      s_crit = Seq(r_crit, r_dry, T, kappa)
    end if
  end subroutine kohler_crit

  subroutine critical_curve(T, r_a, r_b, kappa, approx, rs, rcrits, scrits, n)
    !! Generate Köhler critical radii and supersaturations across dry sizes.
    !! T [K], r_a,r_b [m] dry radius bounds; kappa [-]; approx logical.
    !! rs(n), rcrits(n), scrits(n) returned.
    real(kind=8), intent(in) :: T, r_a, r_b, kappa
    logical, intent(in) :: approx
    integer, intent(in) :: n
    real(kind=8), intent(out) :: rs(n), rcrits(n), scrits(n)
    integer :: i
    real(kind=8) :: r_crit, s_crit
    do i = 1,n
       rs(i) = 10.0_8**( log10(r_a) + (i-1)*(log10(r_b)-log10(r_a))/(n-1) )
       call kohler_crit(T, rs(i), kappa, approx, r_crit, s_crit)
       rcrits(i) = r_crit
       scrits(i) = s_crit
    end do
  end subroutine critical_curve

  function r_eff(rho, wc, Ni) result(r_eff_val)
    !! Compute effective radius given liquid water mixing ratio.
    !! rho [kg/m^3], wc [-], Ni [m^-3].
    !! Returns r_eff_val [m].
    real(kind=8), intent(in) :: rho, wc, Ni
    real(kind=8) :: r_eff_val
    r_eff_val = (3.0_8*rho*wc/(4.0_8*3.141592653589793_8*rho_w*Ni))**(1.0_8/3.0_8)
  end function r_eff

end module thermo
