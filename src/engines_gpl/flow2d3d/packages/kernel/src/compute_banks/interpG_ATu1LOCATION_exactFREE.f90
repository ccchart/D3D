subroutine interpG_ATu1LOCATION_exactFREE(u     , v      , kcs, lunscr, mmax, nmax, nmaxus, kmax    , &
                                        & nst   , nlb    , nub, mlb   , mub , nmlb, nmub  , exitloop, &
                                        & iterFR, MAXaguu, gdp)
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  Ndryact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: My_intersec.f90 
!  $HeadURL:
!!--description-----------------------------------------------------------------
!
!   Function:   Interpolate the values at the Image points
!
!   NOTE: Vandermonde matrices are notoriously ill-conditioned. The inverses of large floating-point
! Vandermonde matrices are subject to severe round-off effects.fROM http://www.mathworks.com/help/symbolic/mupad_ref/linalg-invvandermonde.html
!
!   Author: Alberto Canestrelli 
!                        
!!--declarations----------------------------------------------------------------
!
    use globaldata
    use mathconsts, only: pi
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    real(fp)                  , pointer :: tolFREEexact
    integer                   , pointer :: subtypeTESTghost
    integer                   , pointer :: freeNONhomo
    integer                   , pointer :: totGHOSTu1
    integer                   , pointer :: typeEXTRAPstencil
    integer, dimension(:,:)   , pointer :: FROMmnTOghostU1
    integer, dimension(:)     , pointer :: nGPu1
    integer, dimension(:)     , pointer :: mGPu1
    integer, dimension(:)     , pointer :: mIPu1
    integer, dimension(:)     , pointer :: nIPu1
    integer, dimension(:)     , pointer :: mBIu1
    integer, dimension(:)     , pointer :: nBIu1
    integer, dimension(:,:)   , pointer :: GHOSTu1
    real(fp), dimension(:,:)  , pointer :: PSIx
    real(fp), dimension(:,:)  , pointer :: PSIy
    real(fp), dimension(:,:)  , pointer :: ETAx
    real(fp), dimension(:,:)  , pointer :: ETAy
    real(fp), dimension(:,:,:), pointer :: dutdn_U
    real(fp), dimension(:,:)  , pointer :: aguu
    real(fp), dimension(:,:)  , pointer :: xcorU1
    real(fp), dimension(:,:)  , pointer :: xcorV1
    real(fp), dimension(:,:)  , pointer :: ycorV1
    real(fp), dimension(:,:)  , pointer :: xG_V1
    real(fp), dimension(:,:)  , pointer :: xG_U1
    real(fp), dimension(:,:)  , pointer :: yG_U1
    real(fp), dimension(:)    , pointer :: xIPu1
    real(fp), dimension(:)    , pointer :: yIPu1
    real(fp), dimension(:)    , pointer :: xBIu1
    real(fp), dimension(:)    , pointer :: yBIu1
    real(fp), dimension(:,:)  , pointer :: u1IPu
    real(fp), dimension(:,:)  , pointer :: v1IPu
    real(fp), dimension(:)    , pointer :: nxG_U1
    real(fp), dimension(:)    , pointer :: nyG_U1
    real(fp), dimension(:,:)  , pointer :: shiftBIu_x
    real(fp), dimension(:,:)  , pointer :: shiftBIu_y
    logical                   , pointer :: periodSURFACE
    logical                   , pointer :: solvedU_ATu
    logical                   , pointer :: solvedV_ATu
!
! global variables
!
    real(fp), dimension(nlb:nub,mlb:mub,kmax)                           , intent(inout) :: u
    real(fp), dimension(nlb:nub,mlb:mub,kmax)                           , intent(inout) :: v
    real(fp)                                                            , intent(in)    :: MAXaguu !max value of aguu for which I prescribe ghost cells
    integer, dimension(nlb:nub,mlb:mub)                                 , intent(in)    :: kcs
    integer                                                             , intent(in)    :: mmax   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: nmax   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: nmaxus   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: kmax   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: lunscr   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: nst   !  Description and declaration in iidim.f90
    integer                                                             , intent(in)    :: nlb
    integer                                                             , intent(in)    :: nub
    integer                                                             , intent(in)    :: mlb
    integer                                                             , intent(in)    :: mub
    integer                                                             , intent(in)    :: nmlb
    integer                                                             , intent(in)    :: nmub
    integer                                                             , intent(inOUT) :: exitloop
    integer                                                             , intent(in)    :: iterFR
!
!
! local variables
!
  integer                    :: i
  integer                    :: j
  integer                    :: md
  integer                    :: nu
  integer                    :: m
  integer                    :: n
  integer                    :: mGP
  integer                    :: nGP
  integer                    :: k
  integer                    :: kk
  integer                    :: kkk
  integer                    :: kfl
  integer                    :: Kdry
  integer                    :: k1
  integer                    :: k2
  integer                    :: k3
  integer                    :: k4
  integer                    :: mBI
  integer                    :: nBI
  integer                    :: mBI1
  integer                    :: nBI1
  integer                    :: mBI2
  integer                    :: nBI2
  integer                    :: iG1
  integer                    :: iG2
  integer                    :: nL(4)
  integer                    :: mL(4)
  integer                    :: INFO
  integer                    :: contFLUID
  integer                    :: contGHOST
  integer                    :: contDRYnoGH   
  integer                    :: contNNghostWD
  integer                    :: nDRY   
  integer                    :: mDRY  
  integer                    :: signINT
  integer                    :: cont
!
  integer                    :: kOK(4)
  integer                    :: kGHOS(4)
  integer                    :: kDRYnoGH(4)
  integer                    :: kNNG(4)
  integer                    :: iwork(8)
  integer                    :: IPIV(8) 
!
  real(fp)                   :: nxDRY1
  real(fp)                   :: nyDRY1
  real(fp)                   :: nxDRY2
  real(fp)                   :: nyDRY2
  real(fp)                   :: signREAL
  real(fp)                   :: xBI_DRY1
  real(fp)                   :: yBI_DRY1
  real(fp)                   :: xBI_DRY2
  real(fp)                   :: yBI_DRY2
  real(fp)                   :: halfPSI
  real(fp)                   :: halfETA
  real(fp)                   :: RsignBC
  real(fp)                   :: angleADJ
  real(fp)                   :: dx
  real(fp)                   :: dy
  real(fp)                   :: x1
  real(fp)                   :: y1
  real(fp)                   :: rcond
  real(fp)                   :: colsum
  real(fp)                   :: anorm
!
  logical                    :: REDO
  logical                    :: DRYUadj
  logical                    :: setTOzero
!
  real(fp)                   :: butta
  real(fp)                   :: A(4,4)
  real(fp)                   :: A3(3,3)
  real(fp)                   :: A6(6,6)
  real(fp)                   :: A8(8,8)
  real(fp)                   :: B(4)
  real(fp)                   :: B3(3)
  real(fp)                   :: B6(6)
  real(fp)                   :: B8(8)
  real(fp)                   :: work(16)
  real(fp)                   :: work8(64)
  real(fp)                   :: work6(36)
  real(fp)                   :: iniVALUE(kmax)         
  real(fp)                   :: dutdy(kmax)  
  real(fp)                   :: dutdx(kmax)   
  real(fp)                   :: ut(4,kmax)   
  real(fp)                   :: xx(3)
  real(fp)                   :: yy(3)
  real(fp)                   :: DN       
!
  character, dimension(1)    :: norm
!
! executable statements -------------------------------------------------------
!  
    tolFREEexact      => gdp%gdimbound%tolFREEexact
    subtypeTESTghost  => gdp%gdimbound%subtypeTESTghost
    freeNONhomo       => gdp%gdimbound%freeNONhomo
    totGHOSTu1        => gdp%gdimbound%totGHOSTu1
    typeEXTRAPstencil => gdp%gdimbound%typeEXTRAPstencil
    FROMmnTOghostU1   => gdp%gdimbound%FROMmnTOghostU1
    nGPu1             => gdp%gdimbound%nGPu1
    mGPu1             => gdp%gdimbound%mGPu1
    mIPu1             => gdp%gdimbound%mIPu1
    nIPu1             => gdp%gdimbound%nIPu1
    mBIu1             => gdp%gdimbound%mBIu1
    nBIu1             => gdp%gdimbound%nBIu1
    GHOSTu1           => gdp%gdimbound%GHOSTu1
    PSIx              => gdp%gdimbound%PSIx
    PSIy              => gdp%gdimbound%PSIy
    ETAx              => gdp%gdimbound%ETAx
    ETAy              => gdp%gdimbound%ETAy
    dutdn_U           => gdp%gdimbound%dutdn_U
    aguu              => gdp%gdimbound%aguu
    xcorU1            => gdp%gdimbound%xcorU1
    xcorV1            => gdp%gdimbound%xcorV1
    ycorV1            => gdp%gdimbound%ycorV1
    xG_V1             => gdp%gdimbound%xG_V1
    xG_U1             => gdp%gdimbound%xG_U1
    yG_U1             => gdp%gdimbound%yG_U1
    xIPu1             => gdp%gdimbound%xIPu1
    yIPu1             => gdp%gdimbound%yIPu1
    xBIu1             => gdp%gdimbound%xBIu1
    yBIu1             => gdp%gdimbound%yBIu1
    u1IPu             => gdp%gdimbound%u1IPu
    v1IPu             => gdp%gdimbound%v1IPu
    nxG_U1            => gdp%gdimbound%nxG_U1
    nyG_U1            => gdp%gdimbound%nyG_U1
    shiftBIu_x        => gdp%gdimbound%shiftBIu_x
    shiftBIu_y        => gdp%gdimbound%shiftBIu_y
    periodSURFACE     => gdp%gdimbound%periodSURFACE
    solvedU_ATu       => gdp%gdimbound%solvedU_ATu
    solvedV_ATu       => gdp%gdimbound%solvedV_ATu
!
!  extrapolate at the boundary, in order to have some sort of tranmissive behaviour
!
    if (typeEXTRAPstencil.ge.1.and..not.periodSURFACE) then ! if periodic tangential velocity is already extrapolated
       do m = 1, mmax
          md = m-1
          do n = 1, nmaxus - 1
             nu = n + 1
   !         extrapolate vertically 
             if (kcs(n, m) == 2 .and. kcs(nu, m) == 1) then
                if (ghostu1(n , m )/=0) u(n, m, 1:kmax)  = u(nu, m, 1:kmax)
                if (ghostu1(n , md)/=0) u(n, md, 1:kmax) = u(nu, md, 1:kmax)
             elseif (kcs(n, m) == 1 .and. kcs(nu, m) == 2) then
                if (ghostu1(nu, m )/=0) u(nu, m, 1:kmax)  = u(n, m, 1:kmax)
                if (ghostu1(nu, md)/=0) u(nu, md, 1:kmax) = u(n, md, 1:kmax)
             endif
          enddo
       enddo
    endif
!
    do i = 1,totGHOSTu1       
!        
       m = mIPu1(i)
       n = nIPu1(i)
       mGP = mGPu1(i)
       nGP = nGPu1(i)
       if  (comparereal(aguu(nGP,mGP),MAXaguu).gt.0) cycle
       mBI = mBIu1(i)
       nBI = nBIu1(i)
       setTOzero = .false.
       iniVALUE(1:kmax) = u(nGP,mGP,1:kmax)
!
!      if m equal to 1 do an extrapolation from the stencil corresponding to n=2, otherwise 2 nodes of the stencil are outside the domain
!
       if (m.eq.1) then 
          m = m + 1
       elseif (m.eq.mmax) then
          m = m - 1
       endif
!
       contFLUID = 0
       contGHOST = 0
       contDRYnoGH = 0
       contNNghostWD = 0
!
      ! lower left  
       nL(1)=n-1
       mL(1)=m-1 
       ! lower right
       nL(2)=n-1
       mL(2)=m
       ! upper right
       nL(3)=n
       mL(3)=m  
       ! upper left  
       nL(4)=n 
       mL(4)=m-1

       Do K=1,4
         IF (GHOSTu1(nL(K)+1,mL(K)) == 0) then !fluid cell ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
            contFLUID = contFLUID + 1
            kOK(contFLUID) = K
         ELSEIF (GHOSTu1(nL(K)+1,mL(K)) == 1) then ! dry ghost cell ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
            contGHOST = contGHOST + 1
            kGHOS(contGHOST) = K
         ELSEIF (GHOSTu1(nL(K)+1,mL(K)) == 2) then ! dry NOT ghost cell 
            contDRYnoGH = contDRYnoGH + 1
            kDRYnoGH(contDRYnoGH) = K
         ELSEIF (GHOSTu1(nL(K)+1,mL(K)) >= 3) then ! bank/wet channel interface or orthogonal to it
            contNNghostWD = contNNghostWD + 1
            kNNG(contNNghostWD) = K
         endif
      ENDDO   
!
!     compute homogeneus part of free slip condition
!
      if (freeNONhomo==0) then
         dutdn_U(n,m,1:kmax)=0
         DN = 0
      elseif (freeNONhomo==1) then
         DN = sqrt( (xIPu1(i)-xG_U1(nGP,mGP))**2 + (yIPu1(i)-yG_U1(nGP,mGP))**2 )  !should be precomputed
         if (contFLUID+contGHOST==4) then
            !NOTE: VALID ONLY FOR GRID WITH M AND N ORIENTED IN X AND Y DIRECTION. OTHERWISE dutdx,dutdy has to be rotated and dx and dy has to be computed as distance (or from guu/gvv)
            do k=1,4
               ut(k,1:kmax) = - u(nL(k)+1,mL(k),1:kmax)*nyG_U1(i) + v(nL(k)+1,mL(k),1:kmax)*nxG_U1(i)
            enddo
            dx = xG_U1(nL(2)+1,mL(2))-xG_U1(nL(1)+1,mL(1)) !same as xG_U1(nL(3),mL(3))-xG_U1(nL(4),mL(4)) for orthogonal grid
            do k=1,kmax
               dutdx(k) = 0.5_fp * ( ( ut(2,k)-ut(1,k) )/dx ) + &  !lower right minus lower left 
                          0.5_fp * ( ( ut(3,k)-ut(4,k) )/dx )      !upper right minus upper left 
            enddo
            dy = yG_U1(nL(4)+1,mL(4))-yG_U1(nL(1)+1,mL(1)) !same as yG_U1(nL(3),mL(3))-yG_U1(nL(2),mL(2)) for orthogonal grid
            do k=1,kmax
               dutdy(k) = 0.5_fp * ( ( ut(4,k)-ut(1,k) )/dy ) + &  !lower right minus lower left 
                          0.5_fp * ( ( ut(3,k)-ut(2,k) )/dy )      !upper right minus upper left 
            enddo
            dutdn_U(n,m,1:kmax) = dutdx(1:kmax)*nxG_U1(i) + dutdy(1:kmax)*nyG_U1(i)
           ! write(*,*) 'check if normal goes toward dry otherwise change sign'
            write(1000000,'(4i8,f25.15)') iterFR,i,mGP,nGP,dutdn_U(n,m,1)
         else if(contFLUID+contGHOST==3) then
            !slopes given by plane between 3 points
            cont = 0
            do k=1,4
               if (GHOSTu1(nL(K)+1,mL(K))<=1) then !ghost or fluid
                  cont = cont + 1
                  ut(cont,1:kmax) = - u(nL(k)+1,mL(k),1:kmax)*nyG_U1(i) + v(nL(k)+1,mL(k),1:kmax)*nxG_U1(i)
                  xx(cont) = xG_U1(nL(k),mL(k))
                  yy(cont) = yG_U1(nL(k),mL(k))
               endif
            enddo
            !next two can be optimized
            dutdy(1:kmax) = ( (ut(2,1:kmax)-ut(1,1:kmax))*(xx(1)-xx(3))+(ut(3,1:kmax)-ut(1,1:kmax))*(xx(2)-xx(1)) )/ &
                            ( (yy(2)       -yy(1)       )*(xx(1)-xx(3))+(yy(3)       -yy(1)       )*(xx(2)-xx(1)) )  !slope along n
            dutdx(1:kmax) = ( (ut(2,1:kmax)-ut(1,1:kmax))-dutdy(1:kmax)*(yy(2)-yy(1)) )/(xx(2)-xx(1))!slope along m
            dutdn_U(n,m,1:kmax) =  dutdx(1:kmax)*nxG_U1(i) + dutdy(1:kmax)*nyG_U1(i)  !note:  nxG_U1,nyG_U1 points toward water so also dutdn_U points toward wet area (see my notes)
            write(1000000,'(4i8,f25.15)') iterFR,i,mGP,nGP,dutdn_U(n,m,1)
         else !it should never happen
            write(*,*) 'Error: cannot compute derivative at image point'
            call d3stop(1, gdp)
         endif
      elseif (freeNONhomo==5) then !exact solution circular channel !2D
         DN = sqrt( (xIPu1(i)-xG_U1(nGP,mGP))**2 + (yIPu1(i)-yG_U1(nGP,mGP))**2 )  !should be precomputed
         if (subtypeTESTghost==3) then ! anular channel 
            if (xBIu1(i)**2+yBIu1(i)**2>3600._fp ) then !R^2>60^2, where 60 is axis of the channel 
               dutdn_U(n,m,1:kmax) =  +0.008479956984545_fp !2D
            else !if<60
               dutdn_U(n,m,1:kmax) =  +0.012771766940083_fp !2D
            endif
         elseif (subtypeTESTghost==2) then ! straight channel 43/80, slope 0.001 orthogonal
            if ( yBIu1(i)>0.525_fp*xBIu1(i)+250._fp ) then !R^2>60^2, where 60 is axis of the channel 
               dutdn_U(n,m,1:kmax) =  -0.001_fp !2D
            else !if<60
               dutdn_U(n,m,1:kmax) =  -0.001_fp !2D
            endif   
         endif
      endif
       !
       ! note: Bilinear interpolation for a rectangle has a simple form (Nasr-Azadani, E. Meiburg 2009)
       ! to espress the inverse so  I think it should not be solved the linear system. However, 
       ! if the rectangle is not with edges along x and y teh form is not that simple (but actually it is if I use
       ! geometric interpretation)so I prefer to always
       ! solve the linear system. Moreover the form is not simple when I use the boundary condition itself for
       ! the interpolation, even for dirichlet BC (its not a rectangle anymore)
       !
      ! write(3434343,'(5i6)')  nst,n,m,contFLUID,contGHOST
       IF (contFLUID==4) then      
!
!         simplest case: all the four corner nodes of the interpolation cell are on the fluid
!         
          A(1,1) = xcorV1(nL(1),mL(1))*ycorV1(nL(1),mL(1)) 
          A(2,1) = xcorV1(nL(2),mL(2))*ycorV1(nL(2),mL(2))
          A(3,1) = xcorV1(nL(3),mL(3))*ycorV1(nL(3),mL(3))
          A(4,1) = xcorV1(nL(4),mL(4))*ycorV1(nL(4),mL(4))
          A(1,2) = xcorV1(nL(1),mL(1))  
          A(2,2) = xcorV1(nL(2),mL(2)) 
          A(3,2) = xcorV1(nL(3),mL(3)) 
          A(4,2) = xcorV1(nL(4),mL(4)) 
          A(1,3) = ycorV1(nL(1),mL(1)) 
          A(2,3) = ycorV1(nL(2),mL(2))
          A(3,3) = ycorV1(nL(3),mL(3))
          A(4,3) = ycorV1(nL(4),mL(4))
          A(1,4) = 1._fp
          A(2,4) = 1._fp
          A(3,4) = 1._fp
          A(4,4) = 1._fp
!
          CALL DGETRF( contFLUID, contFLUID, A, contFLUID, IPIV, INFO ) !compute LU
          CALL ERROR_DGETRF(INFO,contFLUID,mBI,nBI)   
          DO K=1,kmax
             !solve linear system for u
             if (.NOT.solvedU_ATu) THEN
                B(1) = u(nL(1)+1,mL(1),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(2) = u(nL(2)+1,mL(2),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(3) = u(nL(3)+1,mL(3),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(4) = u(nL(4)+1,mL(4),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                CALL DGETRS('N', contFLUID, 1, A, contFLUID, IPIV, B, contFLUID, INFO ) 
                CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)  
                u1IPu(i,k) = B(1)*xIPu1(i)*yIPu1(i)+B(2)*xIPu1(i)+B(3)*yIPu1(i)+B(4)              
             endif
             !solve linear system for v
             if (.NOT.solvedV_ATu) THEN
                B(1) = v(nL(1)+1,mL(1),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(2) = v(nL(2)+1,mL(2),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(3) = v(nL(3)+1,mL(3),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B(4) = v(nL(4)+1,mL(4),k) 
                CALL DGETRS('N', contFLUID, 1, A, contFLUID, IPIV, B, contFLUID, INFO ) 
                CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)  
                v1IPu(i,k) = B(1)*xIPu1(i)*yIPu1(i)+B(2)*xIPu1(i)+B(3)*yIPu1(i)+B(4)
             endif
             u(nGP,mGP,k) = - u1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,k)*nxG_U1(i)*nyG_U1(i) + DN*dutdn_U(n,m,k)*nyG_U1(i)
             if (freeNONhomo>0) v(nGP,mGP,k) =   v1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,k)*nxG_U1(i)*nyG_U1(i) - DN*dutdn_U(n,m,k)*nxG_U1(i)
             if ( u(nGP,mGP,k).lt.-1000.or.u(nGP,mGP,k).gt.1000 ) then
                write(*,*) ' wrong value of u in n,m,nst', n,m,nst
                call d3stop(1, gdp)
             endif  
          enddo
 !
       ELSEIF (contFLUID==3) then      
!         three of the  four corner nodes of the interpolation cell are on the fluid
!         
          k1 = kOK(1)  
          k2 = kOK(2) 
          k4 = kOK(3)  !! I put the BC on the third row so the diagonal term is not going to be zero (for zero-velocity prescribed) 
         
          if (contGHOST.eq.1) then
            !note its not necessarely the same ghost cell i!
             k3 = kGHOS(1)
             iG1  = FROMmnTOghostU1(nL(k3)+1,mL(k3))   ! note here m is increased by one because the location of xcorU1(n,m) is the location of xG_v1(n+1,m)
             mBI1 = mBIu1(iG1) !not needed for no slip
             nBI1 = nBIu1(iG1) !not needed for no slip
             nxDRY1 = nxG_U1(iG1)
             nyDRY1 = nyG_U1(iG1) 
             xBI_DRY1 = xBIu1(iG1) + shiftBIu_x(nL(k3)+1,mL(k3)) 
             yBI_DRY1 = yBIu1(iG1) + shiftBIu_y(nL(k3)+1,mL(k3)) 
!
             A8(1,1) = xcorV1(nL(k1),mL(k1))*ycorV1(nL(k1),mL(k1)) 
             A8(2,1) = 0._fp
             A8(3,1) = xcorV1(nL(k2),mL(k2))*ycorV1(nL(k2),mL(k2))
             A8(4,1) = 0._fp
             A8(5,1) = xBI_DRY1*yBI_DRY1*nxDRY1  
             A8(6,1) = - yBI_DRY1*nxDRY1*nyDRY1 - xBI_DRY1*nyDRY1**2
             A8(7,1) = xcorV1(nL(k4),mL(k4))*ycorV1(nL(k4),mL(k4))
             A8(8,1) = 0._fp
            !
             A8(1,2) = 0._fp
             A8(2,2) = xcorV1(nL(k1),mL(k1))*ycorV1(nL(k1),mL(k1)) 
             A8(3,2) = 0._fp
             A8(4,2) = xcorV1(nL(k2),mL(k2))*ycorV1(nL(k2),mL(k2))
             A8(5,2) = xBI_DRY1*yBI_DRY1*nyDRY1  
             A8(6,2) = yBI_DRY1*nxDRY1**2 + xBI_DRY1*nxDRY1*nyDRY1
             A8(7,2) = 0._fp
             A8(8,2) = xcorV1(nL(k4),mL(k4))*ycorV1(nL(k4),mL(k4))
             !
             A8(1,3) = xcorV1(nL(k1),mL(k1)) 
             A8(2,3) = 0._fp
             A8(3,3) = xcorV1(nL(k2),mL(k2)) 
             A8(4,3) = 0._fp
             A8(5,3) = xBI_DRY1*nxDRY1 
             A8(6,3) = - nxDRY1*nyDRY1
             A8(7,3) = xcorV1(nL(k4),mL(k4)) 
             A8(8,3) = 0._fp
             !
             A8(1,4) = 0._fp
             A8(2,4) = xcorV1(nL(k1),mL(k1)) 
             A8(3,4) = 0._fp
             A8(4,4) = xcorV1(nL(k2),mL(k2)) 
             A8(5,4) = xBI_DRY1*nyDRY1 
             A8(6,4) = nxDRY1**2
             A8(7,4) = 0._fp
             A8(8,4) = xcorV1(nL(k4),mL(k4)) 
             !
             A8(1,5) = ycorV1(nL(k1),mL(k1)) 
             A8(2,5) = 0._fp
             A8(3,5) = ycorV1(nL(k2),mL(k2))
             A8(4,5) = 0._fp
             A8(5,5) = yBI_DRY1*nxDRY1
             A8(6,5) = - nyDRY1**2
             A8(7,5) = ycorV1(nL(k4),mL(k4))
             A8(8,5) = 0._fp
             !
             A8(1,6) = 0._fp
             A8(2,6) = ycorV1(nL(k1),mL(k1)) 
             A8(3,6) = 0._fp
             A8(4,6) = ycorV1(nL(k2),mL(k2))
             A8(5,6) = yBI_DRY1*nyDRY1
             A8(6,6) = nxDRY1*nyDRY1
             A8(7,6) = 0._fp
             A8(8,6) = ycorV1(nL(k4),mL(k4))
             !
             A8(1,7) = 1._fp
             A8(2,7) = 0._fp
             A8(3,7) = 1._fp
             A8(4,7) = 0._fp
             A8(5,7) = nxDRY1
             A8(6,7) = 0._fp
             A8(7,7) = 1._fp
             A8(8,7) = 0._fp
             !
             A8(1,8) = 0._fp
             A8(2,8) = 1._fp
             A8(3,8) = 0._fp
             A8(4,8) = 1._fp
             A8(5,8) = nyDRY1
             A8(6,8) = 0._fp
             A8(7,8) = 0._fp
             A8(8,8) = 1._fp
!
             CALL DGETRF( 8, 8, A8, 8, IPIV, INFO ) !compute LU
             CALL ERROR_DGETRF(INFO,contFLUID,mBI,nBI)      
!            INSTEAD OF DOING A CYCLE CALL ONCE DGETRS CREATING B8(1:8,1:KMAX) AND calling CALL CALL DGETRS('N', 8, 8 ,.... instead of DGETRS('N', 8, 1
             do k=1,kmax
                B8(1) = u(nL(k1)+1,mL(k1),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B8(2) = v(nL(k1)+1,mL(k1),k)
                B8(3) = u(nL(k2)+1,mL(k2),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B8(4) = v(nL(k2)+1,mL(k2),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B8(5) = 0._fp
                B8(6) = dutdn_U(n,m,k)
                B8(7) = u(nL(k4)+1,mL(k4),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                B8(8) = v(nL(k4)+1,mL(k4),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
   !  
                CALL DGETRS('N', 8, 1, A8, 8, IPIV, B8, 8, INFO ) !solve linear system
                CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)     
                u1IPu(i,k) = B8(1)*xIPu1(i)*yIPu1(i)+B8(3)*xIPu1(i)+B8(5)*yIPu1(i)+B8(7)
                v1IPu(i,k) = B8(2)*xIPu1(i)*yIPu1(i)+B8(4)*xIPu1(i)+B8(6)*yIPu1(i)+B8(8)
                u(nGP,mGP,k) = - u1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,k)*nxG_U1(i)*nyG_U1(i) + DN*dutdn_U(n,m,k)*nyG_U1(i)
                if (freeNONhomo>0) v(nGP,mGP,k) =   v1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,k)*nxG_U1(i)*nyG_U1(i) - DN*dutdn_U(n,m,k)*nxG_U1(i)
   !
                if (u(nGP,mGP,k) .lt.-1000.or.u(nGP,mGP,k) .gt.1000) then
                  write(*,*) ' wrong value of u in n,m,nst', n,m,nst
                  call d3stop(1, gdp)
                endif        
             enddo      
!                  
          else !if contghost=0, i.e. contDRYnoGH=1  or contNNghostWD=1
!
            !i use only the three fluids for simplicity
            k1 = kOK(1)   
            k2 = kOK(2)   
            k3 = kOK(3)   
!
            A3(1,1) = xcorV1(nL(k1),mL(k1))
            A3(2,1) = xcorV1(nL(k2),mL(k2))
            A3(3,1) = xcorV1(nL(k3),mL(k3))
            A3(1,2) = ycorV1(nL(k1),mL(k1))
            A3(2,2) = ycorV1(nL(k2),mL(k2))
            A3(3,2) = ycorV1(nL(k3),mL(k3))
            A3(1,3) = 1._fp
            A3(2,3) = 1._fp
            A3(3,3) = 1._fp         
!
            CALL DGETRF( 3, 3, A3, 3, IPIV, INFO ) !compute LU
            CALL ERROR_DGETRF(INFO,contFLUID,mBI,nBI) 
!
            do k=1,kmax
               !solve linear system for u
               if (.NOT.solvedU_ATu) THEN      
                  B3(1) = u(nL(k1)+1,mL(k1),k)
                  B3(2) = u(nL(k2)+1,mL(k2),k)
                  B3(3) = u(nL(k3)+1,mL(k3),k)                      
                  CALL DGETRS('N', 3, 1, A3, 3, IPIV, B3, 3, INFO ) !solve linear system
                  CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)  
                  u1IPu(i,k) = B3(1)*xIPu1(i)+B3(2)*yIPu1(i)+B3(3) 
               endif
               !solve linear system for v
               if (.NOT.solvedV_ATu) THEN 
                  B3(1) = v(nL(k1)+1,mL(k1),k)
                  B3(2) = v(nL(k2)+1,mL(k2),k)
                  B3(3) = v(nL(k3)+1,mL(k3),k)                      
                  CALL DGETRS('N', 3, 1, A3, 3, IPIV, B3, 3, INFO ) !solve linear system
                  CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)  
                  v1IPu(i,k) = B3(1)*xIPu1(i)+B3(2)*yIPu1(i)+B3(3)                  
               endif
               u(nGP,mGP,k) = - u1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,k)*nxG_U1(i)*nyG_U1(i) + DN*dutdn_U(n,m,k)*nyG_U1(i)
               if (freeNONhomo>0) v(nGP,mGP,k) =   v1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,k)*nxG_U1(i)*nyG_U1(i) - DN*dutdn_U(n,m,k)*nxG_U1(i)
   !
               if (u(nGP,mGP,k) .lt.-1000.or.u(nGP,mGP,k) .gt.1000) then
                  write(*,*) ' wrong value of u in n,m,nst', n,m,nst
                  call d3stop(1, gdp)
               endif     
            enddo
         endif
                       
       ELSEIF (contFLUID==2) then      
!         two of the  four corner nodes of the interpolation stencil are on the same 
!         side of the cut-cell boundary  (or in general: i have only 2 fluid points)
!         
         k1 = kOK(1)  
         k4 = kOK(2)  !! I put the BCs on the second and third row so the diagonal term is not going to be zero (for zero-velocity prescribed) 
!       
         if (contGHOST.eq.2) then
 
          !first ghost (note iG1 and iG2 can in general be different from i)
            k2 = kGHOS(1)
            iG1  = FROMmnTOghostU1(nL(k2)+1,mL(k2))   ! note here m is increased by one because the location of xcorU1(n,m) is the location of xG_v1(n+1,m)
            mBI1 = mBIu1(iG1) !not needed for no slip
            nBI1 = nBIu1(iG1) !not needed for no slip
            nxDRY1 = nxG_U1(iG1)
            nyDRY1 = nyG_U1(iG1) 
            xBI_DRY1 = xBIu1(iG1) + shiftBIu_x(nL(k2)+1,mL(k2)) 
            yBI_DRY1 = yBIu1(iG1) + shiftBIu_y(nL(k2)+1,mL(k2)) 
          !second ghost
            k3 = kGHOS(2)
            iG2 = FROMmnTOghostU1(nL(k3)+1,mL(k3))   ! note here m is increased by one because the location of xcorU1(n,m) is the location of xG_v1(n+1,m)
            mBI2 = mBIu1(iG2) !not needed for no slip
            nBI2 = nBIu1(iG2) !not needed for no slip
            nxDRY2 = nxG_U1(iG2)
            nyDRY2 = nyG_U1(iG2) 
            xBI_DRY2 = xBIu1(iG2) + shiftBIu_x(nL(k3)+1,mL(k3)) 
            yBI_DRY2 = yBIu1(iG2) + shiftBIu_y(nL(k3)+1,mL(k3)) 
         elseif ((contGHOST.eq.1).and.(contNNghostWD.EQ.1)) then
          !first ghost (note iG1 and iG2 can in general be different from i)
            k2 = kGHOS(1)
            iG1  = FROMmnTOghostU1(nL(k2)+1,mL(k2))   ! note here m is increased by one because the location of xcorU1(n,m) is the location of xG_v1(n+1,m)
            mBI1 = mBIu1(iG1) !not needed for no slip
            nBI1 = nBIu1(iG1) !not needed for no slip
            nxDRY1 = nxG_U1(iG1)
            nyDRY1 = nyG_U1(iG1) 
            xBI_DRY1 = xBIu1(iG1) + shiftBIu_x(nL(k2)+1,mL(k2)) 
            yBI_DRY1 = yBIu1(iG1) + shiftBIu_y(nL(k2)+1,mL(k2)) 
          !contNNghostWD
            k3 = kNNG(1)
            nBI2 = nL(k3)+1      !n of the non ghost BC
            mBI2 = mL(k3)        !m of the non ghost BC            
            if (ghostU1(nL(k3)+1,mL(k3)).eq.3)  then !on the wet/dry interface
               nxDRY2 =  PSIx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  PSIy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction         
               xBI_DRY2 = xcorV1(nL(k3),mL(k3)) + shiftBIu_x(nL(k3)+1,mL(k3)) 
               yBI_DRY2 = ycorV1(nL(k3),mL(k3)) + shiftBIu_y(nL(k3)+1,mL(k3)) 
            else !if (ghostU1(nL(k3)+1,mL(k3).eq.4)  then !orthogonal to the wet/dry interface
               !check at which side is the image point (it could be a double-valued point, orthogonal to an exact wet-dry interface above and below tiself
               nxDRY2 =  ETAx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  ETAy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction               
               x1 = xcorV1(nL(k3),mL(k3))
               y1 = ycorV1(nL(k3),mL(k3))
               dx = xIPu1(i) - x1  
               dy = yIPu1(i) - y1    
               angleADJ= atan2(dx*PSIy(nBI2,mBI2)-dy*PSIx(nBI2,mBI2),dx*PSIx(nBI2,mBI2)+dy*PSIy(nBI2,mBI2)) 
               DRYUadj = (angleADJ.gt.-pi/2.0_fp.and.angleADJ.lt.pi/2.0_fp)
               if (DRYUadj) then  
                  xBI_DRY2 = xcorV1(nL(k3),mL(k3)) + nxDRY2 + shiftBIu_x(nL(k3)+1,mL(k3)) 
                  yBI_DRY2 = ycorV1(nL(k3),mL(k3)) + nyDRY2 + shiftBIu_y(nL(k3)+1,mL(k3)) 
               else                  
                  xBI_DRY2 = xcorV1(nL(k3),mL(k3)) - nxDRY2 + shiftBIu_x(nL(k3)+1,mL(k3)) 
                  yBI_DRY2 = ycorV1(nL(k3),mL(k3)) - nyDRY2 + shiftBIu_y(nL(k3)+1,mL(k3)) 
               endif              
            endif
         elseif (contNNghostWD.EQ.2) then ! unlucky case (a) on my scanned notes
           !first contNNghostWD  
            k2 = kNNG(1)
            nBI1 = nL(k2)+1      !n of the non ghost BC
            mBI1 = mL(k2)        !m of the non ghost BC
            if (ghostU1(nL(k2)+1,mL(k2)).eq.3)  then !on the wet/dry interface
               nxDRY1 =  PSIx(nBI1,mBI1) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY1 =  PSIy(nBI1,mBI1) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               xBI_DRY1 = xcorV1(nL(k2),mL(k2)) + shiftBIu_x(nL(k2)+1,mL(k2)) 
               yBI_DRY1 = ycorV1(nL(k2),mL(k2)) + shiftBIu_y(nL(k2)+1,mL(k2)) 
            else
               !check at which side is the image point (it could be a double-valued point, orthogonal to an exact wet-dry interface above and below tiself
               nxDRY1 =  ETAx(nBI1,mBI1) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY1 =  ETAy(nBI1,mBI1) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction               
               x1 = xcorV1(nL(k2),mL(k2))
               y1 = ycorV1(nL(k2),mL(k2))
               dx = xIPu1(i) - x1  
               dy = yIPu1(i) - y1    
               angleADJ= atan2(dx*ETAy(nBI1,mBI1)-dy*ETAx(nBI1,mBI1),dx*ETAx(nBI1,mBI1)+dy*ETAy(nBI1,mBI1)) 
               DRYUadj = (angleADJ.gt.-pi/2.0_fp.and.angleADJ.lt.pi/2.0_fp)
               if (DRYUadj) then  
                  xBI_DRY1 = xcorV1(nL(k2),mL(k2)) + nxDRY1 + shiftBIu_x(nL(k2)+1,mL(k2))
                  yBI_DRY1 = ycorV1(nL(k2),mL(k2)) + nyDRY1 + shiftBIu_y(nL(k2)+1,mL(k2))
               else                  
                  xBI_DRY1 = xcorV1(nL(k2),mL(k2)) - nxDRY1 + shiftBIu_x(nL(k2)+1,mL(k2))
                  yBI_DRY1 = ycorV1(nL(k2),mL(k2)) - nyDRY1 + shiftBIu_y(nL(k2)+1,mL(k2))
               endif    
            endif
           !second contNNghostWD
            k3 = kNNG(2)
            nBI2 = nL(k3)+1      !n of the non ghost BC
            mBI2 = mL(k3)        !m of the non ghost BC
            if (ghostU1(nL(k3)+1,mL(k3)).eq.3)  then !on the wet/dry interface
               nxDRY2 =  PSIx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  PSIy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction         
               xBI_DRY2 = xcorV1(nL(k3),mL(k3)) + shiftBIu_x(nL(k3)+1,mL(k3))
               yBI_DRY2 = ycorV1(nL(k3),mL(k3)) + shiftBIu_y(nL(k3)+1,mL(k3))
            else !if (ghostU1(nL(k3)+1,mL(k3).eq.4)  then !orthogonal to the wet/dry interface
               !check at which side is the image point of the i^th GP (it could be a double-valued point, orthogonal to an exact wet-dry interface above and below itself
               nxDRY2 =  ETAx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  ETAy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction               
               x1 = xcorV1(nL(k3),mL(k3))
               y1 = ycorV1(nL(k3),mL(k3))
               dx = xIPu1(i) - x1  
               dy = yIPu1(i) - y1    
               angleADJ= atan2(dx*ETAy(nBI2,mBI2)-dy*ETAx(nBI2,mBI2),dx*ETAx(nBI2,mBI2)+dy*ETAy(nBI2,mBI2)) 
               DRYUadj = (angleADJ.gt.-pi/2.0_fp.and.angleADJ.lt.pi/2.0_fp)
               if (DRYUadj) then  
                  xBI_DRY2 = xcorV1(nL(k3),mL(k3)) + nxDRY2 + shiftBIu_x(nL(k3)+1,mL(k3))
                  yBI_DRY2 = ycorV1(nL(k3),mL(k3)) + nyDRY2 + shiftBIu_y(nL(k3)+1,mL(k3))
               else                  
                  xBI_DRY2 = xcorV1(nL(k3),mL(k3)) - nxDRY2 + shiftBIu_x(nL(k3)+1,mL(k3))
                  yBI_DRY2 = ycorV1(nL(k3),mL(k3)) - nyDRY2 + shiftBIu_y(nL(k3)+1,mL(k3))
               endif              
            endif    
         else
            !if (mod(nst,80).eq.0) write(*,*) 'Vel set to 0: 2 cells in the stencil are dry non-veg/dry bank non-ghost cells for u1 at (m,n)', mGP,ngp  
            !SHOULD not occur often, only in the case depicted in singlurCASE_2FLUIDSoneGHOST.bmp
            setTOzero = .true.
            !call d3stop(1, gdp)
         endif
!
         if (.not.setTOzero) then
!
             A8(1,1) = xcorV1(nL(k1),mL(k1))*ycorV1(nL(k1),mL(k1)) 
             A8(2,1) = 0._fp
             A8(3,1) = xBI_DRY1*yBI_DRY1*nxDRY1  
             A8(4,1) = - yBI_DRY1*nxDRY1*nyDRY1 - xBI_DRY1*nyDRY1**2
             A8(5,1) = xBI_DRY2*yBI_DRY2*nxDRY2  
             A8(6,1) = - yBI_DRY2*nxDRY2*nyDRY2 - xBI_DRY2*nyDRY2**2
             A8(7,1) = xcorV1(nL(k4),mL(k4))*ycorV1(nL(k4),mL(k4))
             A8(8,1) = 0._fp
            !
             A8(1,2) = 0._fp
             A8(2,2) = xcorV1(nL(k1),mL(k1))*ycorV1(nL(k1),mL(k1)) 
             A8(3,2) = xBI_DRY1*yBI_DRY1*nyDRY1  
             A8(4,2) = yBI_DRY1*nxDRY1**2 + xBI_DRY1*nxDRY1*nyDRY1
             A8(5,2) = xBI_DRY2*yBI_DRY2*nyDRY2 
             A8(6,2) = yBI_DRY2*nxDRY2**2 + xBI_DRY2*nxDRY2*nyDRY2
             A8(7,2) = 0._fp
             A8(8,2) = xcorV1(nL(k4),mL(k4))*ycorV1(nL(k4),mL(k4))
             !
             A8(1,3) = xcorV1(nL(k1),mL(k1)) 
             A8(2,3) = 0._fp
             A8(3,3) = xBI_DRY1*nxDRY1 
             A8(4,3) = - nxDRY1*nyDRY1
             A8(5,3) = xBI_DRY2*nxDRY2 
             A8(6,3) = - nxDRY2*nyDRY2
             A8(7,3) = xcorV1(nL(k4),mL(k4)) 
             A8(8,3) = 0._fp
             !
             A8(1,4) = 0._fp
             A8(2,4) = xcorV1(nL(k1),mL(k1)) 
             A8(3,4) = xBI_DRY1*nyDRY1 
             A8(4,4) = nxDRY1**2
             A8(5,4) = xBI_DRY2*nyDRY2
             A8(6,4) = nxDRY2**2
             A8(7,4) = 0._fp
             A8(8,4) = xcorV1(nL(k4),mL(k4)) 
             !
             A8(1,5) = ycorV1(nL(k1),mL(k1)) 
             A8(2,5) = 0._fp
             A8(3,5) = yBI_DRY1*nxDRY1
             A8(4,5) = - nyDRY1**2
             A8(5,5) = yBI_DRY2*nxDRY2
             A8(6,5) = - nyDRY2**2
             A8(7,5) = ycorV1(nL(k4),mL(k4))
             A8(8,5) = 0._fp
             !
             A8(1,6) = 0._fp
             A8(2,6) = ycorV1(nL(k1),mL(k1)) 
             A8(3,6) = yBI_DRY1*nyDRY1
             A8(4,6) = nxDRY1*nyDRY1
             A8(5,6) = yBI_DRY2*nyDRY2
             A8(6,6) = nxDRY2*nyDRY2
             A8(7,6) = 0._fp
             A8(8,6) = ycorV1(nL(k4),mL(k4))
             !
             A8(1,7) = 1._fp
             A8(2,7) = 0._fp
             A8(3,7) = nxDRY1
             A8(4,7) = 0._fp
             A8(5,7) = nxDRY2
             A8(6,7) = 0._fp
             A8(7,7) = 1._fp
             A8(8,7) = 0._fp
             !
             A8(1,8) = 0._fp
             A8(2,8) = 1._fp
             A8(3,8) = nyDRY1
             A8(4,8) = 0._fp
             A8(5,8) = nyDRY2
             A8(6,8) = 0._fp
             A8(7,8) = 0._fp
             A8(8,8) = 1._fp
             !
            CALL DGETRF( 8, 8, A8, 8, IPIV, INFO ) !compute LU
            !compute 1-norm needed for condition number
            anorm = 0.d0
            do j=1,8
               colsum = 0.d0
               do k=1,8
                  colsum = colsum + abs(A8(k,j))
               enddo
               anorm = max(anorm, colsum)
            enddo
            norm = '1'  ! use 1-norm to compute condition number with dgecon
            CALL DGECON(norm,8,A8,8,anorm,rcond,work8,iwork,info) !http://www.mathworks.com/support/solutions/en/data/1-3KL67Y/?product=SL&solution=1-3KL67Y
            if (rcond.gt.0.000000000000001_fp) then ! 
               CALL ERROR_DGETRF(INFO,contFLUID,mBI,nBI)                  
               do k=1,kmax 
                  B8(1) = u(nL(k1)+1,mL(k1),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                  B8(2) = v(nL(k1)+1,mL(k1),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                  B8(3) = 0._fp                 
                  B8(4) = dutdn_U(n,m,k)     
                  B8(5) = 0._fp
                  B8(6) = dutdn_U(n,m,k)
                  B8(7) = u(nL(k4)+1,mL(k4),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                  B8(8) = v(nL(k4)+1,mL(k4),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)                  
   !                   
                  CALL DGETRS('N', 8, 1, A8, 8, IPIV, B8, 8, INFO ) !solve linear system
                  CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)     
                  u1IPu(i,k) = B8(1)*xIPu1(i)*yIPu1(i)+B8(3)*xIPu1(i)+B8(5)*yIPu1(i)+B8(7)
                  v1IPu(i,k) = B8(2)*xIPu1(i)*yIPu1(i)+B8(4)*xIPu1(i)+B8(6)*yIPu1(i)+B8(8)
                  u(nGP,mGP,k) = - u1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,k)*nxG_U1(i)*nyG_U1(i) + DN*dutdn_U(n,m,k)*nyG_U1(i)
                  if (freeNONhomo>0) v(nGP,mGP,k) =   v1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,k)*nxG_U1(i)*nyG_U1(i) - DN*dutdn_U(n,m,k)*nxG_U1(i)
   !                  
                  if (u(nGP,mGP,k).lt.-1000.or.u(nGP,mGP,k).gt.1000) then
                     write(*,*) ' wrong value of u in n,m,nst', n,m,nst
                     call d3stop(1, gdp)
                  endif   
               enddo
!
            else !The matrix for u1 is singular at machine precision just compute the average of the fluid cell points
               write(*,*) 'Matrix for u1 is singular, rcond = ',rcond
               u(nGP,mGP,1:kmax) = (u(nL(k1)+1,mL(k1),1:kmax) + u(nL(k4)+1,mL(k4),1:kmax))*0.5_fp  !(B8(1)+B8(7))*0.5_fp 
               if (freeNONhomo>0) v(nGP,mGP,1:kmax) = (v(nL(k1)+1,mL(k1),1:kmax) + v(nL(k4)+1,mL(k4),1:kmax) ) *0.5_fp              
            endif
 
         ELSE ! IF(setTOzero) 
           !average of fluid points
            u1IPu(i,1:kmax) = (u(nL(k1)+1,mL(k1),1:kmax) + u(nL(k4)+1,mL(k4),1:kmax) ) *0.5_fp
            v1IPu(i,1:kmax) = (v(nL(k1)+1,mL(k1),1:kmax) + v(nL(k4)+1,mL(k4),1:kmax) ) *0.5_fp
            u(nGP,mGP,1:kmax) = - u1IPu(i,1:kmax)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,1:kmax)*nxG_U1(i)*nyG_U1(i) !+ DN*dutdn_U(n,m,1:kmax) !if singular likely there is no plane 
            if (freeNONhomo>0) v(nGP,mGP,1:kmax) =   v1IPu(i,1:kmax)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,1:kmax)*nxG_U1(i)*nyG_U1(i) !- DN*dutdn_U(n,m,1:kmax) !if singular likely there is no plane 
         ENDIF
!
! 
       ELSEIF (contFLUID==1) then      
!         one of the  four corner nodes of the interpolation cell are on the fluid 
!         
          !first fluid cell
          k3 = kOK(1)  

          if (contGHOST.ge.2) then !even if its 3 I just use the first 2 (rare case, see CASEwith3GHOSTSand1FLUIDpoint.bmp mine)
!                       
             !first ghost (note iG1 and iG2 can in general be different from i)
             k1 = kGHOS(1)
             iG1  = FROMmnTOghostU1(nL(k1)+1,mL(k1)) 
             mBI1 = mBIu1(iG1)
             nBI1 = nBIu1(iG1)
             xBI_DRY1 = xBIu1(iG1) + shiftBIu_x(nL(k1)+1,mL(k1))
             yBI_DRY1 = yBIu1(iG1) + shiftBIu_y(nL(k1)+1,mL(k1)) 
             nxDRY1 = nxG_U1(iG1)
             nyDRY1 = nyG_U1(iG1)
             !second ghost
             k2 = kGHOS(2)  
             iG2 = FROMmnTOghostU1(nL(k2)+1,mL(k2)) 
             mBI2 = mBIu1(iG2)
             nBI2 = nBIu1(iG2) 
             xBI_DRY2 = xBIu1(iG2) + shiftBIu_x(nL(k2)+1,mL(k2)) 
             yBI_DRY2 = yBIu1(iG2) + shiftBIu_y(nL(k2)+1,mL(k2)) 
             nxDRY2 = nxG_U1(iG2)
             nyDRY2 = nyG_U1(iG2)
!
          elseif(contGHOST.eq.1.and.contNNghostWD.ge.1) then  
!   
            !first ghost (note iG1 and iG2 can in general be different from i)
             k1 = kGHOS(1) 
             iG1  = FROMmnTOghostU1(nL(k1)+1,mL(k1)) 
             mBI1 = mBIu1(iG1)
             nBI1 = nBIu1(iG1)
             xBI_DRY1 = xBIu1(iG1) + shiftBIu_x(nL(k1)+1,mL(k1)) 
             yBI_DRY1 = yBIu1(iG1) + shiftBIu_y(nL(k1)+1,mL(k1)) 
             nxDRY1 = nxG_U1(iG1)
             nyDRY1 = nyG_U1(iG1)
             !look for the second dry point at the wet/dry interface
             k2 = kNNG(1)
             nBI2 = nL(k2)+1      !n of the non ghost BC
             mBI2 = mL(k2)        !m of the non ghost BC
             if (ghostU1(nL(k2)+1,mL(k2)).eq.3)  then
               nxDRY2 =  PSIx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  PSIy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction         
               xBI_DRY2 = xcorV1(nL(k2),mL(k2)) + shiftBIu_x(nL(k2)+1,mL(k2))
               yBI_DRY2 = ycorV1(nL(k2),mL(k2)) + shiftBIu_y(nL(k2)+1,mL(k2))
             else !if (ghostU1(nL(k2)+1,mL(k2).eq.4)  then
               !check at which side is the image point (it could be a double-valued point, orthogonal to an exact wet-dry interface above and below tiself
               nxDRY2 =  ETAx(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction 
               nyDRY2 =  ETAy(nBI2,mBI2) ! I figured out that clearly the orientation (verso) does not matter to impose the derivative equal to zero, only the direction               
               x1 = xcorV1(nL(k2),mL(k2))
               y1 = ycorV1(nL(k2),mL(k2))
               dx = xIPu1(i) - x1  
               dy = yIPu1(i) - y1    
               angleADJ= atan2(dx*ETAy(nBI2,mBI2)-dy*ETAx(nBI2,mBI2),dx*ETAx(nBI2,mBI2)+dy*ETAy(nBI2,mBI2)) 
               DRYUadj = (angleADJ.gt.-pi/2.0_fp.and.angleADJ.lt.pi/2.0_fp)
               if (DRYUadj) then  
                  xBI_DRY2 = xcorV1(nL(k2),mL(k2)) + nxDRY2 + shiftBIu_x(nL(k2)+1,mL(k2))
                  yBI_DRY2 = ycorV1(nL(k2),mL(k2)) + nyDRY2 + shiftBIu_y(nL(k2)+1,mL(k2))
               else                  
                  xBI_DRY2 = xcorV1(nL(k2),mL(k2)) - nxDRY2 + shiftBIu_x(nL(k2)+1,mL(k2))
                  yBI_DRY2 = ycorV1(nL(k2),mL(k2)) - nyDRY2 + shiftBIu_y(nL(k2)+1,mL(k2))
               endif              
             endif
                  
          else 
             if (mod(nst,40).eq.0) WRITE(*,*) 'Error for interpolation for u1: only one fluid cell and less then two ghosts/boundary points. BI in cell',mBI,nBI
             setTOzero = .true.
             !    pause
             !stop
          endif    
!
          if (.not.setTOzero) then
             A6(1,1) = xBI_DRY1*nxDRY1 
             A6(2,1) = - nxDRY1*nyDRY1
             A6(3,1) = xBI_DRY2*nxDRY2 
             A6(4,1) = - nxDRY2*nyDRY2
             A6(5,1) = xcorV1(nL(k3),mL(k3)) 
             A6(6,1) = 0._fp
            !
             A6(1,2) = xBI_DRY1*nyDRY1 
             A6(2,2) = nxDRY1**2
             A6(3,2) = xBI_DRY2*nyDRY2
             A6(4,2) = nxDRY2**2
             A6(5,2) = 0._fp
             A6(6,2) = xcorV1(nL(k3),mL(k3)) 
            !
             A6(1,3) = yBI_DRY1*nxDRY1
             A6(2,3) = - nyDRY1**2
             A6(3,3) = yBI_DRY2*nxDRY2
             A6(4,3) = - nyDRY2**2
             A6(5,3) = ycorV1(nL(k3),mL(k3))
             A6(6,3) = 0._fp
            !
             A6(1,4) = yBI_DRY1*nyDRY1
             A6(2,4) = nxDRY1*nyDRY1
             A6(3,4) = yBI_DRY2*nyDRY2
             A6(4,4) = nxDRY2*nyDRY2
             A6(5,4) = 0._fp
             A6(6,4) = ycorV1(nL(k3),mL(k3))
            !
             A6(1,5) = nxDRY1
             A6(2,5) = 0._fp
             A6(3,5) = nxDRY2
             A6(4,5) = 0._fp
             A6(5,5) = 1._fp
             A6(6,5) = 0._fp
            !
             A6(1,6) = nyDRY1
             A6(2,6) = 0._fp
             A6(3,6) = nyDRY2
             A6(4,6) = 0._fp
             A6(5,6) = 0._fp
             A6(6,6) = 1._fp       
   !
             CALL DGETRF( 6, 6, A6, 6, IPIV, INFO ) !compute LU
             if (info==0) then
   !
                CALL ERROR_DGETRF(INFO,contFLUID,mBI,nBI)   
   !
                do k=1,kmax
                   B6(1) = 0._fp                 
                   B6(2) = dutdn_U(n,m,k)                
                   B6(3) = 0._fp
                   B6(4) = dutdn_U(n,m,k)
                   B6(5) = u(nL(k3)+1,mL(k3),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)
                   B6(6) = v(nL(k3)+1,mL(k3),k)   ! note here n is increased by one because the location of xcorV1(n,m) is the location of xG_u1(n+1,m)                  
   ! 
                   CALL DGETRS('N', 6, 1, A6, 6, IPIV, B6, 6, INFO ) !solve linear system
                   CALL ERROR_DGETRS(INFO,contFLUID,mBI,nBI)   
                   u1IPu(i,k) = B6(1)*xIPu1(i)+B6(3)*yIPu1(i)+B6(5) 
                   v1IPu(i,k) = B6(2)*xIPu1(i)+B6(4)*yIPu1(i)+B6(6) 
                   u(nGP,mGP,k) = - u1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*v1IPu(i,k)*nxG_U1(i)*nyG_U1(i) + DN*dutdn_U(n,m,k)*nyG_U1(i)
                   if (freeNONhomo>0) v(nGP,mGP,k) =   v1IPu(i,k)*(nxG_U1(i)**2-nyG_U1(i)**2)-2._fp*u1IPu(i,k)*nxG_U1(i)*nyG_U1(i) - DN*dutdn_U(n,m,k)*nxG_U1(i)
                   if (u1IPu(i,k).lt.-1000.or.u1IPu(i,k).gt.1000) then
                     write(*,*) ' wrong value of u in n,m,nst', n,m,nst
                     call d3stop(1, gdp)
                   endif  
                enddo
   !
             else  
                write(*,*) 'The two normals for u1 are coincident for ghost point (mGP,nGP)=,',mGP,nGP,'at time step ',nst
                 !call d3stop(1,gdp)  !uncomment this to set the value constant in the cell
                 u(nGP,mGP,1:kmax) = u(nL(k3)+1,mL(k3),1:kmax) 
                 if (freeNONhomo>0) v(nGP,mGP,1:kmax) = v(nL(k3)+1,mL(k3),1:kmax) 
             endif
          else !if settozero=true
             u(nGP,mGP,1:kmax) = u(nL(k3)+1,mL(k3),1:kmax) 
             if (freeNONhomo>0) v(nGP,mGP,1:kmax) = v(nL(k3)+1,mL(k3),1:kmax) 
          endif
       ELSE
          WRITE(515151,*) 'Error for interpolation: NO fluid cells found for u1. BI in cell',mBI,nBI,i
          WRITE(*,*) 'Error for interpolation: NO fluid cells found for u1.  BI in cell',mBI,nBI,i
          call d3stop(1, gdp)
       ENDIF
       !
       ! Check convergence!
       !
       loop_k_1: do k=1,kmax
         ! write(1232112,'(4i8,10f25.15)') iterFR,n,m,K,iniVALUE(k),u(nGP,mGP,k)
          if (comparereal(iniVALUE(k),0._fp)==0) then
             if (  abs(iniVALUE(k) - u(nGP,mGP,k)) > tolFREEexact) then
                exitloop = 0
                if(nst<3) write(2232112,'(5i8,10f25.15)') nst,iterFR,nGP,mGP,K,iniVALUE(k),u(nGP,mGP,k),iniVALUE(k)-u(nGP,mGP,k)
                exit loop_k_1
             endif
          !elseif (  abs((iniVALUE(k) - u(nGP,mGP,k)) / iniVALUE(k)) > tolFREEexact) then
          !   exitloop = 0
          !   if(nst<3) write(2232112,'(5i8,10f25.15)') nst,iterFR,nGP,mGP,K,iniVALUE(k),u(nGP,mGP,k),iniVALUE(k)-u(nGP,mGP,k),(iniVALUE(k) - u(nGP,mGP,k)) / iniVALUE(k)
          !   exit loop_k_1
          endif
       enddo loop_k_1
    enddo
    !
end subroutine interpG_ATu1LOCATION_exactFREE