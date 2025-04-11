!**********************************************************************************************
!
! Implementation of the conservative finite volume for aggregation
! Initial code by Anxiong Liu
! Modified by Stelios Rigopoulos
! 16/09/2024
!
!**********************************************************************************************



!-------------------------------------------------------------------------------
!
!> \brief This subroutine is to locate all the subsections during the coagulation 
!> process for the finite volume scheme. The double array Dou_Array will save 
!> the surface area 'Ajk' for eahc subsection combination and weighted point value 
!> x1 and x2. The integral array Int_Array will save the index number j, k, i for
!> each subsection combination.

! INPUT
!> \param[in]     X            Finite volume grid cell [m^3]
!> \param[in]     CountStore   Count for number of sections or Store information 

! OUTPUT
!> Dou_Array [Ajk, x1, x2]
!> Int_Array [j, k, i]

!-------------------------------------------------------------------------------
subroutine PBE_agg_fvLocate(x, CountStore)

  ! Modules
  use pbe_mod
  use agg_cfv_mod

  ! Header
  implicit none

  ! Declarations
  double precision, intent(in) :: x(0: m)
  integer, intent(in) :: CountStore

  ! Doubles

  double precision :: xiL, xiR            !< I cell left / right boundary
  double precision :: xiHalfL, xiHalfR    !< Half of I cell, left boundary, right boundary
  double precision :: xjL, xjR, xjM       !< J cell(section), left boundary, right boundary, middle point
  double precision :: xkR_old, xkL_old, xkR, xkL  ! Corresponding K section positions
  double precision :: xkL_TR              
  double precision :: xkSupp, xkSupp_old  !< Supplementary K starting position and ending position for coagulation death
  double precision :: xjSecL, xjSecR      !< Left and right boundary for the subsetion in cell J 
  double precision :: xkSecR, xkSecL      !< Right and left boundary for the subsection in block K         
  double precision :: x1, x2              !< Weighted points for subsections in cell J and block K
  double precision :: hj, hk, hjRect      !< Length of J section, Length of K subinterval, Displacement for K section

  double precision :: eps, epsConst       !< Infinitely small value for condition clauses
  double precision :: oneThird, twoThird  !< 1.0/3.0, 2.0/3.0
  double precision :: deltaX(m)           !< Cell length
  double precision :: dnb(m), dnd(m)      !< Number density increment / decrease
  double precision :: dmb(m), dmd(m)      !< First order moment increment / decrease
  double precision :: xMid(m)             !< Cell center value
  double precision :: Nvol(m)             !< number density times cell length
  double precision :: deltahi             !< I cell length
  double precision :: Ajk                 !< Surface are for any subsection combination during coagulation process
  
  integer :: kL, kR, kL_old, kR_old, kL_TR !< Cell index for the positions of K section
  integer :: kSupp, kSupp_old              !< Cell index of the supplermentary parts for coagulation death

  ! Integers
  integer :: i              !< Id of the current cell
  integer :: j              !< Id of daughter cell 1 (= right boundary node id)
  integer :: k              !< Id of daughter cell 2 (= left boundary node id)

  ! Flags
  logical :: ifiMid         !< flag for judging if a daughter section larger than half of cell I
  logical :: ifjNext        !< flag for whether integrating the next daughter section J

  ! Parameter doubles
  parameter (epsConst = 1.0D-15)
  parameter (oneThird = 1.0d0 / 3.0d0)
  parameter (twoThird = 2.0d0 / 3.0d0)

  ! Initialization
  if(CountStore == 2) then 
    allocate(Dou_AggSec(N_AggSec * Ndou_AggSec))
    allocate(Int_AggSec(N_AggSec * Nint_AggSec))
  end if

  index_Sec = 0

  deltaX = 0.0d0
  xMid = 0.0d0

  do i = 1,m      ! I cell coagulation birth

    xiL = x(i - 1)
    xiR = x(i)
    xMid(i) = 0.5d0 * (xiL + xiR)
    deltaX(i) = xiR - xiL
    eps = epsConst * xiR

    xiHalfL = 0.5d0 * xiL
    xiHalfR = 0.5d0 * xiR

    ! Initialize for J
    j = 1
    ifjNext = .true.
    ifiMid = .false.

    xjL = x(j - 1)
    xjR = x(j)

    if(xjL - xiHalfR > -1.0d0 * eps) then
      GO TO 20
    elseif(xjL < xiHalfL) then
      if(xjR > xiHalfL) then
        xjR = xiHalfL
        ifjNext = .false.
      end if
    else
      ifiMid = .true.
      if(xjR > xiHalfR) then
        xjR = xiHalfR
        ifjNext = .false.
       end if
    end if

    if(.Not. ifiMid) then
      xkL = xiL - xjR       !y1
      xkL_old = xiL - xjL   !y2
    else
      xkL = xjR        !y1
      xkL_old = xjL    !y2
    end if

    xkR = xiR - xjR    !y3
    xkR_old = xiR - xjL    !y4

    kR = max(1, i)
    do     ! Find kR_old
      if(xkR_old .le. x(kR - 1))  then     ! Find y4
        kR = kR - 1
      else
        exit
      end if
    end do

    kL = kR
    do     ! Find kL_old
      if(xkL_old .lt. x(kL - 1)) then      ! Find y2
         kL = kL - 1
      else
        exit
      end if
    end do

    do    ! Loop for all J daughter sections

      hj = xjR - xjL
      xjM = 0.5d0 * (xjR + xjL)
      kR_old = kR

      do     ! Find New kR
        if(xkR .le. x(kR - 1)) then     ! Find y3
          kR = kR - 1
        else
          exit
        end if
      end do

      if(.NOT.ifimid) then
        kL_old = kL
        do     ! Find New kL
          if(xkL .lt. x(kL - 1)) then     ! Find y1
            kL = kL - 1
          else
            exit
          end if
        end do
      else
        kL_old = j
        if(ifjNext) then
          kL = j + 1
        else
          kL = j
        end if
      end if

      !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Begin Finding all sections
      if(hj .le. deltaX(i)) then

        ! Right Wing Part
        if(kR_old == kR)  then   ! Whole Triangle

          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjR - twoThird * hj
            x2 = xkR_old - twoThird * hj
            Ajk = hj * hj * 0.5d0
            call Store_DouInt(j, kR_old, i, Ajk, x1, x2) 
            index_Sec = index_Sec + 1
          end if

        else  ! Top Triangle & Trapezoid - Triangle

          xjSecL = xjL
          xkSecR = xkR_old

          do k = kR_old, kR, -1

            if(k > kR) then
              hk = xkSecR - x(k - 1)
            else
              hk = xkSecR - xkR
            end if
          
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjSecL + oneThird * hk
            x2 = xkSecR - twoThird * hk
            Ajk = hk * hk * 0.5d0
            call Store_DouInt(j, k, i, Ajk, x1, x2) 
            index_Sec = index_Sec + 1
          end if

          if(k < kR_old) then
            !Inside Trapezoid - Rectangle

            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then
              Ajk = (xjSecL - xjL) * hk
              x1 = 0.5d0 * (xjL + xjSecL)
              x2 = xkSecR  - 0.5d0 * hk
              call Store_DouInt(j, k, i, Ajk, x1, x2) 
              index_Sec = index_Sec + 1
            end if
 
          end if

          if(k > kR) then
            xjSecL = xjSecL + hk
            xkSecR = xkSecR - hk
          end if

        end do
      end if

      !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Left Triangle Part

      if(.NOT.ifiMid) then     ! NOT IFMID

        if(kL == kL_old) then  ! Whole Triangle
          
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjR - oneThird * hj
            x2 = xkL + twoThird * hj
            Ajk = 0.5d0 * hj * hj
            call Store_DouInt(j, kL_old, i, Ajk, x1, x2) 
            index_Sec = index_Sec + 1
          end if
        
        else
          xjSecR = xjR
          xkSecL = xkL

          do k = kL, kL_old, 1
            if(k < kL_old) then
              hk = x(k) - xkSecL
            else
              hk = xkL_old - xkSecL
            end if

            ! Top Triangle & Trapezoid - Triangle
            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then
              x1 = xjSecR - oneThird * hk
              x2 = xkSecL + twoThird * hk
              Ajk = 0.5d0 * hk * hk
              call Store_DouInt(j, k, i, Ajk, x1, x2) 
              index_Sec = index_Sec + 1
            end if   

            if(k > kL) then
              !Inside Trapezoid - Rectangle
              if(CountStore == 1) then 
                index_Sec = index_Sec + 1
              elseif(CountStore == 2) then
                Ajk = (xjR - xjSecR) * hk
                x1 = 0.5d0 * (xjR + xjSecR)
                x2 = xkSecL  + 0.5d0 * hk
                call Store_DouInt(j, k, i, Ajk, x1, x2) 
                index_Sec = index_Sec + 1
              end if
            
            end if

            if(k < kL_old) then
              xjSecR = xjSecR - hk
              xkSecL = xkSecL + hk
            end if

          end do
        end if

      else   ! IFMid: Middle Left Triangle

        ! Left Triangle
        if(CountStore == 1) then 
          index_Sec = index_Sec + 1
        elseif(CountStore == 2) then
          x1 = xjL + oneThird * hj
          x2 = xkL_old + twoThird * hj
          Ajk = 0.5 * hj * hj
          call Store_DouInt(j, j, i, Ajk, x1, x2) 
          index_Sec = index_Sec + 1
        end if 
      end if  ! Whether IFMID

      !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Middle Rectangle Part
      if(.NOT. ifiMid) then
        kL_TR = kL_old
        xkL_TR = xkL_old
      else
        kL_TR = kL
        xkL_TR = xkL
      end if

      if(kR == kL_TR) then
        hk = xkR - xkL_TR
        if(hk > eps) then
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x2 = 0.5d0 * (xkR + xkL_TR)
            Ajk = hj * hk
            call Store_DouInt(j, kR, i, Ajk, xjM, x2) 
            index_Sec = index_Sec + 1
          end if  
        end if
      else
        do k = kR, kL_TR, -1

          if(k == kR) then
            hk = xkR - x(k - 1)
            x2 = 0.5d0 * (xkR + x(k - 1))
          elseif(k == kL_TR) then
            hk = x(k) - xkL_TR
            x2 = 0.5d0 * (x(k) + xkL_TR)
          else
            hk = deltaX(k)
            x2 = xMid(k)
          end if
         
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            Ajk = hj * hk        
            call Store_DouInt(j, k, i, Ajk, xjM, x2) 
            index_Sec = index_Sec + 1
          end if

        end do
      end if

  !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else !(hj > deltaX(i))
      !Right Wing Part
      if(kR_old == kL_old)  then   ! Whole Triangle
        
        if(CountStore == 1) then 
          index_Sec = index_Sec + 1
        elseif(CountStore == 2) then
          x1 = xjL + oneThird * deltaX(i)
          x2 = xkR_old - twoThird * deltaX(i)
          Ajk = 0.5d0 * deltaX(i) * deltaX(i)
          call Store_DouInt(j, kR_old, i, Ajk, x1, x2) 
          index_Sec = index_Sec + 1
        end if

      else  ! Top Triangle & Trapezoid - Triangle

        xjSecL = xjL
        xkSecR = xkR_old

        do k = kR_old, kL_old, -1

          if(k > kL_old) then
            hk = xkSecR - x(k - 1)
          else
            hk = xkSecR - xkL_old
          end if
        
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjSecL + oneThird * hk
            x2 = xkSecR - twoThird * hk
            Ajk = 0.5d0 * hk * hk
            call Store_DouInt(j, k, i, Ajk, x1, x2)  
            index_Sec = index_Sec + 1
          end if

          if(k < kR_old) then
            !Inside Trapezoid - Rectangle
        
            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then
              x1 = 0.5d0 * (xjSecL +  xjL)
              x2 = xkSecR  - 0.5d0 * hk
              Ajk = (xjSecL - xjL) * hk
              call Store_DouInt(j, k, i, Ajk, x1, x2)  
              index_Sec = index_Sec + 1
            end if

          end if

          if(k > kL_old) then
            xjSecL = xjSecL + hk
            xkSecR = xkSecR - hk
          end if
        end do
      end if

      ! Left Triangle Part

        if(kL == kR) then  ! Whole Triangle
         
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjR - oneThird * deltaX(i)
            x2 = xkL + twoThird * deltaX(i)
            Ajk = 0.5d0 * deltaX(i) * deltaX(i)
            call Store_DouInt(j, kL, i, Ajk, x1, x2)  
            index_Sec = index_Sec + 1
          end if

        else
          xjSecR = xjR
          xkSecL = xkL

          do k = kL, kR, 1
            if(k < kR) then
              hk = x(k) - xkSecL
            else
              hk = xkR - xkSecL
            end if

            ! Top Triangle & Trapezoid - Triangle
            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then 
              x1 = xjSecR - oneThird * hk
              x2 = xkSecL + twoThird * hk
              Ajk = 0.5d0 * hk * hk
              call Store_DouInt(j, k, i, Ajk, x1, x2)  
              index_Sec = index_Sec + 1
            end if

            if(k > kL) then
              !Inside Trapezoid - Rectangle
              if(CountStore == 1) then 
                index_Sec = index_Sec + 1
              elseif(CountStore == 2) then 
                x1 = 0.5d0 * (xjR + xjSecR)
                x2 = xkSecL  + 0.5d0 * hk
                Ajk = (xjR - xjSecR) * hk
                call Store_DouInt(j, k, i, Ajk, x1, x2)  
                index_Sec = index_Sec + 1
              end if

            end if

            if(k < kR) then
              xjSecR = xjSecR - hk
              xkSecL = xkSecL + hk
            end if

          end do
        end if

        ! Rectangle
        if(kR == kL_old) then
          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then 
            x2 = 0.5d0 * (xkR + xkL_old)
            x1 = xMid(i) - x2
            Ajk = deltaX(i) * (xkL_old - xkR)
            call Store_DouInt(j, kR, i, Ajk, x1, x2)  
            index_Sec = index_Sec + 1
          end if

        else

          do k = kL_old, kR, -1

            if(k == kL_old) then
              hk = xkL_old - x(k - 1)
              x2 = 0.5d0 * (xkL_old + x(k - 1))
            elseif(k == kR) then
              hk = x(k) - xkR
              x2 = 0.5d0 * (x(k) + xkR)
            else
              hk = deltaX(k)
              x2 = xMid(k)
            end if
           
            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then 
              x1 = xMid(i) - x2
              Ajk = deltaX(i) * hk
              call Store_DouInt(j, k, i, Ajk, x1, x2)  
              index_Sec = index_Sec + 1
            end if

          end do
        end if
      end if

      if(ifjNext) then
        j = j + 1
      end if

      ! New J Cell
      xjL = xjR
      xjR = x(j)

      ifjNext = .True.

      if(xjL - xiHalfR > -1.0D0 * eps) then
        exit
      elseif(xjL < xiHalfL) then
        if(xjR > xiHalfL) then
          xjR = xiHalfL
          ifjNext = .false.
        end if
      else
        ifiMid = .true.
        if(xjR > xiHalfR) then
          xjR = xiHalfR
          ifjNext = .false.
        end if
      end if

      xkL_old = xkL    !y2
      if(.Not. ifiMid) then
        xkL = xiL - xjR    !y1
      else
        xkL = xjR    !y1
      end if

      xkR_old = xkR      !y4
      xkR = xiR - xjR    !y3

    end do

20  CONTINUE

  end do

  ! Supplement for Coagulation Death
    kSupp = m
    xkSupp = x(m) - x(0)
    xiHalfR = 0.5d0 * x(m)
    do
      if(xkSupp.lt.x(kSupp - 1)) then
        kSupp = kSupp - 1
      else
        exit
      end if
    end do

    do j = 1, m
      kSupp_old = kSupp
      xkSupp_old = xkSupp

      if(x(j) .le. xiHalfR) then
        xjR = x(j)
        hj = deltaX(j)
        ifjNext = .True.
        x1 = xMid(j)
      else
        xjR = xiHalfR
        hj = xjR - x(j - 1)
        ifjNext = .false.
        x1 = 0.5d0 * (xjR + x(j - 1))
      end if

      ! Rectangle
      do k = m, kSupp_old + 1, -1
           
        if(CountStore == 1) then 
          index_Sec = index_Sec + 1
        elseif(CountStore == 2) then
          x2 = xMid(k)
          Ajk = hj * deltaX(k)
          call Store_DouInt(j, k, m + 1, Ajk, x1, x2)  
          index_Sec = index_Sec + 1
        end if
      end do
      
      if(CountStore == 1) then 
        index_Sec = index_Sec + 1
      elseif(CountStore == 2) then
        x2 = 0.5d0 * (xkSupp_old + x(kSupp_old))
        hk = x(kSupp_old) - xkSupp_old
        Ajk = hj * hk
        call Store_DouInt(j, kSupp_old, m + 1, Ajk, x1, x2)   
        index_Sec = index_Sec + 1
      end if
     
      ! Triangle
      xkSupp = x(m) - xjR
      do
        if(xkSupp.lt.x(kSupp - 1)) then
          kSupp = kSupp - 1
        else
          exit
        end if
      end do

      if(kSupp == kSupp_old) then
        
        if(CountStore == 1) then 
          index_Sec = index_Sec + 1
        elseif(CountStore == 2) then
          x1 = xjR - onethird * hj
          x2 = xkSupp_old - onethird * hj
          Ajk = 0.5d0 * deltaX(j) * hj
          call Store_DouInt(j, kSupp_old, m + 1, Ajk, x1, x2) 
          index_Sec = index_Sec + 1
        end if  

      else
        xjSecR = xjR
        xkSecL = xkSupp

        do k = kSupp, kSupp_old, 1
          if(k < kSupp_old) then
            hk = x(k) - xkSecL
          else
            hk = xkSupp_old - xkSecL
          end if
 
          ! Top Triangle & Trapezoid - Triangle      

          if(CountStore == 1) then 
            index_Sec = index_Sec + 1
          elseif(CountStore == 2) then
            x1 = xjSecR - oneThird * hk
            x2 = xkSecL + twoThird * hk
            Ajk = 0.5d0 * hk * hk
            call Store_DouInt(j, k, m + 1, Ajk, x1, x2)   
            index_Sec = index_Sec + 1
          end if

          if(k > kSupp) then
            !Inside Trapezoid - Rectangle
            if(CountStore == 1) then 
              index_Sec = index_Sec + 1
            elseif(CountStore == 2) then
              Ajk = (xjR - xjSecR) * hk
              x1 = 0.5d0 * (xjR + xjSecR)
              x2 = xkSecL  + 0.5d0 * hk
              call Store_DouInt(j, k, m + 1, Ajk, x1, x2)   
              index_Sec = index_Sec + 1
            end if
          end if

          xjSecR = xjSecR - hk
          xkSecL = xkSecL + hk

        end do
      end if

      if(.NOT.ifjNext) exit

    end do

  ! Save the total number of sbusection combinations
  if(CountStore == 1) then
    N_AggSec = index_Sec
    ! allocate kernel arrays
    allocate(beta_AggSecOrig(N_AggSec)) 
    allocate(beta_AggSec(N_AggSec)) 
  end if

end subroutine PBE_agg_fvLocate


!-------------------------------------------------------------------------------
!> \brief This subroutine is to calculate the coagulation kernel for each subsection 
!> combinations. 
!-------------------------------------------------------------------------------
subroutine PBE_agg_beta(iflag)
  
  ! Modules
  use pbe_mod
  use agg_cfv_mod

  ! Header
  implicit none
 
  integer, intent(in)                  :: iflag

  ! Doubles

  integer :: preID_dou, preID_int      !< temporary index of subsection combination
  double precision :: x1, x2, b        !< weighted J postion, K position and their colgulation kernel 
  integer :: i                         !< index number
  
  if(iflag == 1) then
    do i = 1, N_AggSec
      preID_dou = (i - 1) * NDou_AggSec
      ! Retrieve all weighted J position, K position
      x1 = Dou_AggSec(preID_dou + ID_xj) 
      x2 = Dou_AggSec(preID_dou + ID_xk)
      ! Calculate coagualtion kernel for each subsection combinations  
      call agg_kernel_compute(x1, x2, b)
      beta_AggSecOrig(i) = b  
    end do
  elseif(iflag == 2) then
    beta_AggSec = beta_AggSecOrig * agg_kernel_const
  end if 

end subroutine PBE_agg_beta


!-------------------------------------------------------------------------------
!> \brief This subroutine is to calculate the coagulation death and birth for each 
!> subsection combinations
!-------------------------------------------------------------------------------

subroutine pbe_agg_cfv(deltaX,xMid,N,dn)

  use pbe_mod
  use agg_cfv_mod

  implicit none

  double precision, intent(in) :: deltaX(m)
  double precision, intent(in) :: xMid(m)
  double precision, intent(in) :: n(m)
  double precision, intent(inout) :: dn(m)

  ! Integers & Doubles  
  integer :: index                    !< Index of a subsection combination
  integer :: preID_dou, preID_int     !< temporary index of subsection combination
  integer :: i, j, k                  !< Id of the current cell, the daughter cell J, the daughter block K
  double precision :: deltaN, deltaMj, deltaMk    !< Number density change due to coagulation, first order moment change in cell J and block K
  double precision :: dnb(m), dnd(m)  !< Number density increment / decrease
  double precision :: dmb(m), dmd(m)  !< First order moment increment / decrease
  double precision :: x1, x2, Ajk     !< Weight position for daughter cell J and daughter block K

  ! Initialization 
  dmb = 0.0
  dmd = 0.0
  dnb = 0.0
  dnd = 0.0 

  do index = 1, N_AggSec
    preID_dou = (index - 1) * NDou_AggSec
    preID_int = (index - 1) * NInt_AggSec
 
    ! Retrieve all indexes of corresponding cell for each subsection combination
    i = Int_AggSec(preID_int + ID_i)
    j = Int_AggSec(preID_int + ID_j)
    k = Int_AggSec(preID_int + ID_k)
    ! Retrieve all the double precision values for each subsection combination
    Ajk = Dou_AggSec(preID_dou + ID_Ajk)
    x1 = Dou_AggSec(preID_dou + ID_xj)
    x2 = Dou_AggSec(preID_dou + ID_xk)
    ! Calculate number density change and moment change in corresponding cell
    deltaN = Ajk * beta_AggSec(index) * n(j) * n(k)     
    deltaMj = deltaN * x1
    deltaMk = deltaN * x2
    dmd(j) = dmd(j) - deltaMj
    dmd(k) = dmd(k) - deltaMk

    if(i .le. m) then
      dmb(i) = dmb(i) + deltaMj + deltaMk
    end if
  end do

  do i = 1, m

!    dnb(i) = dnb(i) / deltahi
!    dnd(i) = dnd(i) / deltahi
    dnb(i) = dmb(i) / deltaX(i) / xMid(i)
    dnd(i) = dmd(i) / deltaX(i) / xMid(i)

  end do
  
! Update number density change in each cell
  dn = dn + dnb + dnd

end subroutine pbe_agg_cfv

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> \brief This subroutine is to calculate the coagulation death and birth for each 
!> subsection combinations

! INPUT
!> \param[in]     J, K, I  Indexes of the daughter cells and their parent cell 
!>                         for each subsection combination  
!> \param[in]     AJK      The surface area for each subsection combination
!> \param[in]     X1, X2   Weighted positions of the daughter particles for each subsection combination
!-------------------------------------------------------------------------------

subroutine Store_DouInt(j, k, i, Ajk, x1, x2)

  use pbe_mod
  use agg_cfv_mod

  implicit none
  
  integer, intent(in) :: j
  integer, intent(in) :: k
  integer, intent(in) :: i
  double precision, intent(in) :: Ajk
  double precision, intent(in) :: x1
  double precision, intent(in) :: x2  

  ! Integers & Doubles  
  integer :: preID_dou, preID_int

  preID_dou = index_Sec * Ndou_AggSec
  preID_int = index_Sec * Nint_AggSec
 
  ! Integers 
  INT_AggSec(preID_int + ID_j) = j
  INT_AggSec(preID_int + ID_i) = i
  INT_AggSec(preID_INT + ID_k) = k
  ! Doubles
  Dou_AggSec(preID_dou + ID_Ajk) = Ajk
  Dou_AggSec(preID_dou + ID_xj) = x1
  Dou_AggSec(preID_dou + ID_xk) = x2

end subroutine Store_DouInt