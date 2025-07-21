!**********************************************************************************************
!
! PBE test subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_test_agg_const()

!**********************************************************************************************
!
! Testing of the PBE solver
! Test case: aggregation with constant kernel
! Reference for analytical solution: form of Seinfeld, original reference is
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, allocatable :: ni(:),nia(:)

double precision moment(0:1),meansize
double precision int_time,dt,vmean_init
double precision m0,m1,tau

integer i,i_step,n_steps,i_writesp

!**********************************************************************************************

! Parameters
agg_kernel = 1
agg_kernel_const = 1.
max_nuc = 0
growth_function = 0
break_const = 0.
grid_type = 1
m = 80
grid_lb = 0.
grid_rb = 1.E3
v0 = 0.1
initdis = 0 ! To be set in this subroutine
N0 = 1.
solver_pbe = 1
int_time = 10.
dt = 1.E-2
i_writesp = 0
vmean_init = 1.

! Initialise PBE
allocate(ni(m))
call pbe_grid()
call pbe_init(ni)
allocate(nia(m))

! Exponential initial distribution
do i=1,m
  ni(i) = (N0/vmean_init)*exp(-v_m(i)/vmean_init)
  if (abs(ni(i))<1.D-16) then
    ni(i) = 0.D0
  end if
end do

! Write initial moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of initial distribution', moment(0)
write(*,*) 'First moment of initial distribution', moment(1)
write(*,*) 'Mean size of initial distribution:', meansize

! Initialise PSR integration
n_steps = int_time/dt

! Integration
do i_step = 1,n_steps
  call pbe_integ(ni,dt)
end do

! Write final moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of final distribution', moment(0)
write(*,*) 'First moment of final distribution', moment(1)
write(*,*) 'Mean size of final distribution:',meansize

! Write final PSD
call pbe_output(ni,i_writesp)

! Calculate and write analytical solution
open(99,file='pbe/psd_analytical.out')
tau = agg_kernel_const*N0*int_time
do i=1,m
  nia(i) = agg_kernel_const*(4.*N0/(tau+2.)**2)*exp((-v_m(i)/vmean_init)*(2./(tau+2.)))
  if (abs(nia(i))<1.D-16) then
    nia(i) = 0.D0
  end if
  write(99,1001) v_m(i),nia(i),v_m(i)*nia(i)
end do
close(99)

1001 format(3E20.10)

end subroutine pbe_test_agg_const

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_test_agg_Brownian()

!**********************************************************************************************
!
! Testing of the PBE solver
! Test case: Brownian aggregation
!
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, allocatable :: ni(:)

double precision moment(0:1),meansize
double precision int_time,dt,vmean_init
double precision m0,m1

integer i_step,n_steps,i_writesp

!**********************************************************************************************

! Parameters
agg_kernel = 5
agg_kernel_const = 1.
max_nuc = 0
growth_function = 0
break_const = 0.
grid_type = 1
m = 80
grid_lb = 0.05
grid_rb = 1.E4
v0 = 1.
initdis = 2
n_th1 = 1
n_th2 = 1
N0 = 1.
solver_pbe = 1
int_time = 100.
dt = 1.E-2
i_writesp = 1

! Initialise PBE
allocate(ni(m))
call pbe_grid()
call pbe_init(ni)

! Initialise PSR integration
n_steps = int_time/dt

! Write initial moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of initial distribution', moment(0)
write(*,*) 'First moment of initial distribution', moment(1)
write(*,*) 'Mean size of initial distribution:', meansize

! Integration
do i_step = 1,n_steps
  call pbe_integ(ni,dt)
end do

! Write final moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of final distribution', moment(0)
write(*,*) 'First moment of final distribution', moment(1)
write(*,*) 'Mean size of final distribution:',meansize

! Write final PSD
call pbe_output(ni,i_writesp)

end subroutine pbe_test_agg_Brownian

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_test_agg_growth()

!**********************************************************************************************
!
! Testing of the PBE solver
! Test case: aggregation with constant kernel and size-dependent growth
! Reference for analytical solution:
! Ramabhadran, T. E., Peterson, T. W., and Seinfeld, J. H. 1976. Dynamics of aerosol 
! coagulation and condensation. AIChE Journal, 22(5), 840–851.
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, allocatable :: ni(:),nia(:),CN(:)

double precision moment(0:1),meansize
double precision int_time,dt,vmean_init
double precision m0,m1

integer i,i_step,n_steps,i_writesp

!**********************************************************************************************

! Parameters
agg_kernel = 1
agg_kernel_const = 1.
max_nuc = 0
growth_function = 2
break_const = 0.
g_coeff1 = 0.5
g_coeff2 = 1.
i_gm = 1
grid_type = 1
m = 80
grid_lb = 0.
grid_rb = 1.E5
v0 = 0.0001
initdis = 0 ! To be set in this subroutine
N0 = 1.
solver_pbe = 1
int_time = 1.
dt = 1.E-2
i_writesp = 0
vmean_init = 1.

! Initialise PBE
allocate(ni(m))
call pbe_grid()
call pbe_init(ni)
allocate(nia(m),CN(m))

! Exponential initial distribution
do i=1,m
  ni(i) = (N0/vmean_init)*exp(-v_m(i)/vmean_init)
  if (abs(ni(i))<1.D-16) then
    ni(i) = 0.D0
  end if
end do

! Write initial moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of initial distribution', moment(0)
write(*,*) 'First moment of initial distribution', moment(1)
write(*,*) 'Mean size of initial distribution:',meansize

! Initialise PSR integration
n_steps = int_time/dt

do i=1,m
  CN(i) = G_coeff1*v(i)**G_coeff2*dt/dv(i)
end do
write (*,*) 'Maximum Courant number: ',maxval(CN)

! Integration
do i_step = 1,n_steps
  call pbe_integ(ni,dt)
end do

! Write final moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of final distribution', moment(0)
write(*,*) 'First moment of final distribution', moment(1)
write(*,*) 'Mean size of final distribution:',meansize

! Write final PSD
call pbe_output(ni,i_writesp)

! Calculate and write analytical solution
open(99,file='pbe/psd_analytical.out')
m0 = 2.D0 * N0 / (2.D0 + agg_kernel_const * N0 * int_time)
m1 = N0 * vmean_init * exp(g_coeff1 * int_time)
do i=1,m
  nia(i) = (m0**2/m1)*exp(-m0/m1*v_m(i))
  if (abs(nia(i))<1.D-16) then
    nia(i) = 0.D0
  end if
  write(99,1001) v_m(i),nia(i),v_m(i)*nia(i)
end do
close(99)

1001 format(3E20.10)

end subroutine pbe_test_agg_growth

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_test_nuc_growth()

!**********************************************************************************************
!
! Testing of the PBE solver
! Test case: nucleation and size-independent growth
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, allocatable :: ni(:),nia(:),CN(:)

double precision moment(0:1),meansize
double precision int_time,dt
double precision m0,m1

integer i,i_step,n_steps,i_writesp

!**********************************************************************************************

! Parameters
agg_kernel = 0
max_nuc = 1
nuc1 = 1.
growth_function = 1
break_const = 0.
g_coeff1 = 1.
grid_type = 2
m = 2000
grid_lb = 0.
grid_rb = 2.
initdis = 1
i_gm = 0
n_th1 = 1
n_th2 = m/10
N0 = 1.
solver_pbe = 1
int_time = 1.
dt = 1.E-5
i_writesp = 0

! Initialise PBE
allocate(ni(m))
call pbe_grid()
call pbe_init(ni)
allocate(nia(m),CN(m))

! Initialise PSR integration
n_steps = int_time/dt

do i=1,m
  CN(i) = G_coeff1*dt/dv(i)
end do
write (*,*) 'Maximum Courant number: ',maxval(CN)

! Integration
do i_step = 1,n_steps
  call pbe_integ(ni,dt)
end do

! Write final PSD
call pbe_output(ni,i_writesp)

1001 format(3E20.10)

end subroutine pbe_test_nuc_growth

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_test_frag_un_bin()

!**********************************************************************************************
!
! Testing of the PBE solver
! Test case: uniform binary fragmentation
! Stelios Rigopoulos
! 9/10/2024
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, allocatable :: ni(:),nia(:)

double precision moment(0:1),meansize
double precision int_time,dt,vmean_init
double precision m0,m1

integer i,i_step,n_steps,i_writesp

!**********************************************************************************************

! Parameters
agg_kernel = 0
max_nuc = 0
growth_function = 0
break_const = 1.
order_of_gq = 3
grid_type = 1
m = 80
grid_lb = 0.
grid_rb = 1.E1
v0 = 1.E-6 !0.1
initdis = 0 ! To be set in this subroutine
N0 = 1.
solver_pbe = 1
int_time = 1.
dt = 1.E-2
i_writesp = 0
vmean_init = 1.

! Initialise PBE
allocate(ni(m))
call pbe_grid()
call pbe_init(ni)
allocate(nia(m))

! Exponential initial distribution
do i=1,m
  ni(i) = (N0/vmean_init)*exp(-v_m(i)/vmean_init)
  if (abs(ni(i))<1.D-16) then
    ni(i) = 0.D0
  end if
end do

! Write initial moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of initial distribution', moment(0)
write(*,*) 'First moment of initial distribution', moment(1)
write(*,*) 'Mean size of initial distribution:', meansize

! Initialise PSR integration
n_steps = int_time/dt

! Integration
do i_step = 1,n_steps
  call pbe_integ(ni,dt)
end do

! Write final moments
call pbe_moments(ni,moment,meansize)
write(*,*)
write(*,*) 'Zeroth moment of final distribution', moment(0)
write(*,*) 'First moment of final distribution', moment(1)
write(*,*) 'Mean size of final distribution:',meansize

! Write final PSD
call pbe_output(ni,i_writesp)

! Calculate and write analytical solution
open(99,file='pbe/psd_analytical.out')
do i=1,m
  nia(i) = (1+int_time)**2*exp(-v_m(i)*(1+int_time))
  if (abs(nia(i))<1.D-16) then
    nia(i) = 0.D0
  end if
  write(99,1001) v_m(i),nia(i),v_m(i)*nia(i)
end do
close(99)

1001 format(3E20.10)

end subroutine pbe_test_frag_un_bin

!**********************************************************************************************



!**********************************************************************************************

subroutine discrete_pbe_test_agg_Brownian()

!**********************************************************************************************
!
! Testing of the discrete PBE solver
! Test case: self-preserving distribution for Brownian kernel
!
! Stelios Rigopoulos
! 1/10/2023
!
!**********************************************************************************************

use pbe_discrete_mod

implicit none

double precision, allocatable :: ni(:),eta(:),psi(:)

double precision moment(0:1)
double precision int_time,dt,current_time,meansize

integer i,i_step,n_steps

!**********************************************************************************************

!Parameters
agg_kernel = 5
agg_kernel_const = 1.
m = 4000
N0 = 1000
v0 = 1.
solver_pbe_discrete = 1
int_time = 0.1
dt = 1.E-4

!Initialise PBE
allocate(v(m))
do i=1,m
  v(i) = i*v0
end do
call pbe_discrete_agg_init()
allocate(ni(m),eta(m),psi(m))

!Top hat initial distribution
ni = 0.D0
ni(1) = N0

! Initialise PSR integration
n_steps = int_time/dt

!Calculate initial moments
call pbe_discrete_moments(ni,moment,meansize)
write(*,'(F10.5,E20.5,E20.5)') current_time,moment(0),moment(1)

!Integration
dt = int_time/n_steps
current_time= 0.D0
do i_step = 1,n_steps
  call pbe_discrete_integ(ni,dt)
  call pbe_discrete_moments(ni,moment,meansize)
  current_time = current_time + dt
  write(*,'(F10.5,E20.5,E20.5)') current_time,moment(0),moment(1)
end do

!Write-up
do i=1,m
  if (abs(ni(i))<1.D-16) then
    ni(i) = 0.D0
  end if
end do
open(99,file='pbe/psd_discrete.out')
do i=1,m
  write(99,1134) v(i),ni(i)
end do
close(99)
open(99,file='pbe/psd_sp_discrete.out')
do i=1,m
  eta(i) = moment(0)*v(i)/moment(1)
  psi(i) = ni(i)*moment(1)/moment(0)**2
  write(99,1134) eta(i),psi(i)
end do

!Finish
deallocate(ni,eta,psi)

1134 format(2E20.10)

end subroutine discrete_pbe_test_agg_Brownian

!**********************************************************************************************