!**********************************************************************************************
!
! Discrete PBE
! Stelios Rigopoulos 1/10/2023
!
!**********************************************************************************************



!**********************************************************************************************

module pbe_discrete_mod

!**********************************************************************************************
!
! Declaration of common variables related to grid and kernels
!
! by Stelios Rigopoulos
!
!**********************************************************************************************

implicit none

save

double precision, allocatable, dimension(:)   :: v
double precision, allocatable, dimension(:,:) :: agg_discrete_factor

double precision agg_kernel_const,v0

integer m,agg_kernel,integer,initdis,n_th1,n_th2,solver_pbe_discrete,N0

end module pbe_discrete_mod

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_discrete_integ(ni,dt)

!**********************************************************************************************
!
! Temporal integration
!
! By Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

external pbe_discrete_ydot

double precision, dimension(m), intent(inout) :: ni
double precision, intent(in)                  :: dt

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------

if (solver_pbe_discrete == 1) then

  !Euler explicit
  call pbe_discrete_ydot(ni,niprime)
  ni = ni + niprime * dt

else if (solver_pbe_discrete == 2) then

  !Runge-Kutta 2nd order
  call pbe_discrete_ydot(ni,niprime)
  nitemp = ni + 0.5D0 * niprime * dt
  call pbe_discrete_ydot(nitemp,niprime)
  ni = ni + niprime * dt

else if (solver_pbe_discrete == 3) then

  !Runge-Kutta 4th order
  call pbe_discrete_ydot(ni,niprime)
  k1 = niprime * dt
  nitemp = ni + 0.5D0 * k1
  call pbe_discrete_ydot(nitemp,niprime)
  k2 = niprime * dt
  nitemp = ni + 0.5D0 * k2
  call pbe_discrete_ydot(nitemp,niprime)
  k3 = niprime * dt
  nitemp = ni + k3
  call pbe_discrete_ydot(nitemp,niprime)
  k4 = niprime * dt
  ni = ni + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

end subroutine pbe_discrete_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_discrete_ydot(y,yprime)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

double precision, dimension(m), intent(in)  :: y
double precision, dimension(m), intent(out) :: yprime

double precision agg_source,agg_sink

double precision agg_total_source,agg_total_sink

integer index

!----------------------------------------------------------------------------------------------

agg_source = 0.D0
agg_sink = 0.D0

agg_total_source = 0.D0
agg_total_sink = 0.D0

do index=1,m

  !Aggregation
  if (agg_kernel>0) then
    if (index>1) then
      call agg_discrete_source(y,index,agg_source)
    end if
    call agg_discrete_sink(y,index,agg_sink)
    agg_total_source = agg_total_source + agg_source
    agg_total_sink = agg_total_sink + agg_sink
  end if

  yprime(index) = agg_source + agg_sink

end do

end subroutine pbe_discrete_ydot

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_discrete_moments(ni,moment,meansize)

!**********************************************************************************************
!
! Calculation of zeroth and first moments
!
! Stelios Rigopoulos
! Modified 1/10/2023
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

double precision, dimension(m), intent(in)    :: ni
double precision, dimension(0:1), intent(out) :: moment
double precision, intent(out)                 :: meansize

double precision m1_lp,lp

integer i

!----------------------------------------------------------------------------------------------

moment(0) = 0.D0
moment(1) = 0.D0

do i=1,m
  moment(0) = moment(0) + ni(i)
  moment(1) = moment(1) + v(i)*ni(i)
end do

m1_lp = 0.D0
do i=m-5,m
  m1_lp = m1_lp + v(i)*ni(i)
end do

lp = m1_lp/moment(1)
if (lp.gt.0.001) then
  write(*,*) 'warning, more than 0.1% of mass in the last five nodes'
end if

if (moment(0).gt.1.D-10) then
  meansize = moment(1)/moment(0)
else
  meansize = 0.D0
end if

end subroutine pbe_discrete_moments

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_discrete_agg_init()

! Tabulation of aggregation operations to speed up aggregation tabulations
! Tabulation depends on mesh, which must be already generated
! Determines complimentaries
!
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

integer i,j

!----------------------------------------------------------------------------------------------

write(*,*) 'Tabulating aggregation terms'

allocate(agg_discrete_factor(m,m))
agg_discrete_factor = 0.D0
do i=1,m
  do j=1,m
    call agg_discrete_kernel_compute(v(i),v(j),agg_discrete_factor(i,j))
  end do
end do

!----------------------------------------------------------------------------------------------

write(*,*) 'preliminary aggregation operations complete'

end subroutine pbe_discrete_agg_init

!**********************************************************************************************



!**********************************************************************************************

subroutine agg_discrete_source(ni,index,agg_source)

!**********************************************************************************************
!
! Calculates the aggregation birth term
!
! By Stelios Rigopoulos
! 28/06/2020
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(out)              :: agg_source

integer i

!----------------------------------------------------------------------------------------------

agg_source = 0.D0

if (index - 2 * floor(real(index) / 2) == 0) then
  !For even numbers, add aggregation of index/2 halved (aggregation between equal volumes)
  do i=1,floor(real(index)/2)-1
    agg_source = agg_source + agg_discrete_factor(i,index-i)*ni(i)*ni(index-i)
  end do
  i = floor(real(index)/2)
  agg_source = agg_source + 0.5D0 * agg_discrete_factor(i,i)*ni(i)*ni(i)
else
  !For odd numbers, do not halve aggregation of index/2
  do i=1,floor(real(index)/2)
    agg_source = agg_source + agg_discrete_factor(i,index-i)*ni(i)*ni(index-i)
  end do
end if

agg_source = agg_kernel_const * agg_source

end subroutine agg_discrete_source

!**********************************************************************************************



!**********************************************************************************************

subroutine agg_discrete_sink(ni,index,agg_sink)

!**********************************************************************************************
!
! Calculates the aggregation death term
!
! By Stelios Rigopoulos
! 28/06/2020
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(out)              :: agg_sink

integer i

!----------------------------------------------------------------------------------------------

agg_sink = 0.D0
do i=1,m
  agg_sink = agg_sink + agg_discrete_factor(index,i)*ni(i)
end do
agg_sink = -agg_kernel_const*agg_sink*ni(index)

end subroutine agg_discrete_sink

!**********************************************************************************************



!**********************************************************************************************

subroutine agg_discrete_kernel_compute(v1,v2,beta)

!**********************************************************************************************

use pbe_discrete_mod

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

end subroutine agg_discrete_kernel_compute

!**********************************************************************************************