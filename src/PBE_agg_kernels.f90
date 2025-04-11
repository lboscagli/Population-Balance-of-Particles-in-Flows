!**********************************************************************************************
!
! PBE kernels
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine agg_kernel_compute(v1,v2,beta)

! This subroutine contains the volume dependent part of aggregation kernels
! Other dependencies (e.g. on temperature) should be implemented in the agg_kernel_const
! factor and updated if needed
!
! The kernel options are as follows:
! 1. Constant
! 2. Sum
! 3. Product
! 4. Free molecular
! 5. Brownian motion
! 6. Shear
! 7. Gravitational (differential sedimentation)

!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)  :: v1
double precision, intent(in)  :: v2
double precision, intent(out) :: beta

!**********************************************************************************************

if (agg_kernel==1) then
  beta = 1.D0
else if (agg_kernel==2) then
  beta = v1+v2
else if (agg_kernel==3) then
  beta = v1*v2
else if (agg_kernel==4) then
  if ((v1==0).or.(v2==0)) then
    beta = 0.D0
  else
    beta = (1.D0/v1+1.D0/v2)**0.5D0*(v1**(1.D0/3.D0)+v2**(1.D0/3.D0)) &
          *(v1**(1.D0/3.D0)+v2**(1.D0/3.D0))
  end if
else if (agg_kernel==5) then
  if ((v1==0).or.(v2==0)) then
    beta = 0.D0
  else
    beta = (1.D0/v1**(1.D0/3.D0)+1.D0/v2**(1.D0/3.D0))*(v1**(1.D0/3.D0)+v2**(1.D0/3.D0))
  end if
else if (agg_kernel==6) then
  beta = (v1**(1.D0/3.D0)+v2**(1.D0/3.D0))*(v1**(1.D0/3.D0)+v2**(1.D0/3.D0)) &
        *(v1**(1.D0/3.D0)+v2**(1.D0/3.D0))
else if (agg_kernel==7) then
  beta = (v1**(1.D0/3.D0)+v2**(1.D0/3.D0))*(v1**(1.D0/3.D0)+v2**(1.D0/3.D0)) &
        *dabs(v1**(2.D0/3.D0)-v2**(2.D0/3.D0))
end if

end subroutine agg_kernel_compute

!**********************************************************************************************