!**********************************************************************************************
!
! Implementation of the conservative finite volume for fragmentation
! Initial code by Daniel O'Sullivan
! Modified by Stelios Rigopoulos
! 10/03/2025
!
!**********************************************************************************************



!**********************************************************************************************

subroutine breakage_fcn(r,s,break_kernel)

use pbe_mod

implicit none

integer :: i,j
integer :: cell_i,cell_j

double precision, intent(in) :: r,s ! dummy particle volumes
double precision, intent(inout) :: break_kernel
double precision :: brate

break_kernel = 2./s

brate = break_const*s
break_kernel = brate*break_kernel

end subroutine

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_breakage_calc(mode_select)

use pbe_mod
use frag_cfv_mod

implicit none

integer, intent(in) :: mode_select
integer :: cntr, i, j
double precision :: brk, blank

cntr = 0
do i = 1,m-1
  do j = i+1, m
    if (mode_select.eq.1) then
      cntr = cntr +1
    else if(mode_select.eq.2) then
      cntr = cntr +1
      call gauss_quad_beta(v(i-1),v(i),v(j-1),v(j),brk,blank)
      call breakage_store2(i,j,brk,cntr)
    end if
  end do
end do

if (mode_select.eq.1) then
  break_global_cntr = cntr
end if

end subroutine pbe_breakage_calc

!**********************************************************************************************



!**********************************************************************************************

subroutine breakage_store2(i,j,brkern, cntr)

use frag_cfv_mod

implicit none

integer, intent(in) :: i,j
integer, intent(in) :: cntr
double precision, intent(in) :: brkern

break_kernel_store(cntr,1) = i
break_kernel_store(cntr,2) = j
break_kernel_store(cntr,3) = brkern

end subroutine breakage_store2

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_breakage_cfv(dn,n)
! Danny O'Sullivan Feb 2020
! executes breakage fractional step

use pbe_mod
use frag_cfv_mod

implicit none

double precision, intent(inout) :: dn(m),n(m)
double precision :: dmb(m),dmd(m),factor, brkern, w_m, u_m,blank
double precision :: d_mass_check
double precision :: bmass_j, bmass_k, d_mass
double precision :: dnb(m), dnd(m)
integer :: i,k,cntr,j,end_break

cntr = 1
dmb = 0.D0
dmd = 0.D0

dnd = 0.D0
dnb = 0.D0

end_break = size(break_kernel_store,1)

factor = 1.D0

do cntr = 1,break_global_cntr

  i = break_kernel_store(cntr,1)
  j = break_kernel_store(cntr,2)
  brkern = break_kernel_store(cntr,3)

  if (isnan(brkern)) then
    write(*,*) 'nan found case 1'
    read(*,*)
  end if

  ! in reworking i is birth cell and j is death cell

  bmass_j = factor*n(j)*brkern

  if (brkern.lt.0.) then
    write(*,*) 'brkern -ve'
    write(*,*) brkern,cntr
    read(*,*)
  end if

  if (isnan(bmass_j)) then
    write(*,*) 'nan found, case 2'
    write(*,*) factor, n(j),'j: ',j
    read(*,*)
  end if

  dmb(j) = dmb(j) - bmass_j
  dmd(i) = dmd(i) + bmass_j

end do

do i = 1, m
  dnb(i) = (dmb(i)) / (dv(i)*v_m(i))
  dnd(i) = (dmd(i))/(dv(i)*v_m(i))
  if (isnan(dnb(i)))then
    write(*,*) 'nan found case 3'
    read(*,*)
  end if
  if (isnan(dnd(i)))then
    write(*,*) 'nan found case 4'
    read(*,*)
  end if
end do

dn = dn + dnb + dnd

end subroutine pbe_breakage_cfv

!**********************************************************************************************



!**********************************************************************************************

subroutine gauss_init(order)

use gauss_mod

implicit none
integer, intent(in) :: order

select case(order)
  case(2)
    allocate(gw(4,2))
    gw = gw_2
    allocate(gp(4,2))
    gp = gp_2
  case(3)
    allocate(gw(9,2))
    gw = gw_3
    allocate(gp(9,2))
    gp = gp_3
  case default
    write(*,*) 'error in gauss quad'
end select

end subroutine gauss_init

!**********************************************************************************************



!**********************************************************************************************

subroutine gauss_quad_beta(w_lb,w_ub,y_lb,y_ub,beta_u,beta_w)
! implements 2 point gaus quad in 2 directions therefore requiring 4 gps.
! a summation over the 4 gps gps obtained through module array gp(4,2)

use pbe_mod
use gauss_mod

implicit none

double precision, intent(inout) :: beta_u,beta_w ! to each k and j death terms
double precision, intent(in) :: w_lb, w_ub,y_lb,y_ub
double precision :: b, x_u, x_w
double precision :: v_lb, v_ub ! helps to rename y_ub/y_lb 
double precision :: delta_w, w_mid, delta_v, v_mid
integer :: i !increment through gp pairs

delta_w = w_ub - w_lb
w_mid = 0.5D0*(w_ub + w_lb)

beta_u = 0.D0
beta_w = 0.D0

delta_v = y_ub - y_lb
v_mid = 0.5D0*(y_ub + y_lb)

do i=1,order_of_gq*order_of_gq
  ! daughter cell
   x_w = 0.5*delta_w*gp(i,2) + w_mid
  ! parent cell
  x_u = 0.5*delta_v*gp(i,1) + v_mid
  call breakage_fcn(x_w,x_u,b)
  beta_u = beta_u + 0.25D0*gw(i,1)*gw(i,2)*delta_w*delta_v*(0.5D0*delta_w*gp(i,2)+ w_mid)*b
end do

end subroutine gauss_quad_beta

!**********************************************************************************************