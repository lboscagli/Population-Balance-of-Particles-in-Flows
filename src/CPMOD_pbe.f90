!**********************************************************************************************

program CPMOD_pbe

!**********************************************************************************************
!
! CPMOD - Chemistry and Particulates MODelling
!
! Stelios Rigopoulos
! 13/12/2015
!
! Modified 23/10/2016
!
!**********************************************************************************************

implicit none

integer :: i,problem

double precision :: time1,time2,time

!**********************************************************************************************

call cpu_time(time1)

! Choose problem
open(30,file='problem.in')
do i=1,2
  read(30,*)
end do
read(30,*) problem
close(30)

if (problem==1) then
  call psr_pbe()
else if (problem==2) then
  call pbe_test_agg_const()
else if (problem==3) then
  call pbe_test_agg_Brownian()
else if (problem==4) then
  call pbe_test_agg_growth()
else if (problem==5) then
  call pbe_test_nuc_growth()
else if (problem==6) then
  call pbe_test_frag_un_bin()
else if (problem==7) then
  call discrete_pbe_test_agg_Brownian()
else
  write(*,*) 'incorrect problem specification'
  stop
end if

call cpu_time(time2)
time = time2-time1
write(*,*)
write(*,*) 'CPU time = ',time

end program CPMOD_pbe

!**********************************************************************************************