!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id$
! $HeadURL$

!> set center related linkxy weights 
 subroutine setlinktocenterweights1d()                 

 use m_flow
 use m_flowgeom
 use m_longculverts
 implicit none
 !
 double precision       :: wud, wuL1, wuL2, cs, sn 
 integer                :: L 
 integer                :: k1, k2 
 integer                :: ilongc, L1dlink
 !
 double precision, external :: lin2nodx, lin2nody

 wc(ndx2D+1:ndx1Db) = 0d0 
 wcxy(1:2,ndx2D+1:ndx1Db) = 0d0
 
 do L = 1, lnx1d

    if (kcu(L) == 3) cycle ! no contribution from 1D2D internal links

    k1   = ln(1,L) ; k2 = ln(2,L) !left and right node
    wud  = wu(L)*dx(L) !flow surface area

    wuL1  = acl(L)*wud ! 2d center factor
    wcL  (1,L ) = wuL1
    wc     (k1) = wc(k1) + wuL1

    wuL2  = (1d0-acl(L))*wud
    wcL  (2, L) = wuL2
    wc     (k2) = wc(k2) + wuL2

    cs = lin2nodx(L,1,csu(L),snu(L))
    sn = lin2nody(L,1,csu(L),snu(L))
    wcx1(L)     = cs*wuL1
    wcy1(L)     = sn*wuL1

    cs = lin2nodx(L,2,csu(L),snu(L))
    sn = lin2nody(L,2,csu(L),snu(L))
    wcx2(L)     = cs*wuL2
    wcy2(L)     = sn*wuL2

    wcxy (1,k1) = wcxy (1,k1) + abs(wcx1(L))
    wcxy (2,k1) = wcxy (2,k1) + abs(wcy1(L))

    wcxy (1,k2) = wcxy (1,k2) + abs(wcx2(L))
    wcxy (2,k2) = wcxy (2,k2) + abs(wcy2(L))
 enddo

 do L = lnxi + 1, lnx1db

    if (kcu(L) == 3) cycle ! no contribution from 1D2D internal links

    k1   = ln(1,L) ; k2 = ln(2,L) !left and right node
    wud  = wu(L)*dx(L) !flow surface area

    wuL1  = acl(L)*wud ! 2d center factor
    wcL  (1,L ) = wuL1
    wc     (k1) = wc(k1) + wuL1

    wuL2  = (1d0-acl(L))*wud
    wcL  (2, L) = wuL2
    wc     (k2) = wc(k2) + wuL2

    cs = lin2nodx(L,1,csu(L),snu(L))
    sn = lin2nody(L,1,csu(L),snu(L))
    wcx1(L)     = cs*wuL1
    wcy1(L)     = sn*wuL1

    cs = lin2nodx(L,2,csu(L),snu(L))
    sn = lin2nody(L,2,csu(L),snu(L))
    wcx2(L)     = cs*wuL2
    wcy2(L)     = sn*wuL2

    wcxy (1,k1) = wcxy (1,k1) + abs(wcx1(L))
    wcxy (2,k1) = wcxy (2,k1) + abs(wcy1(L))

    wcxy (1,k2) = wcxy (1,k2) + abs(wcx2(L))
    wcxy (2,k2) = wcxy (2,k2) + abs(wcy2(L))
 enddo 
 
 ! This code recomputes weighting factors, but subsequently removes the contribution by the long culverts
 if(newculverts) then
    do ilongc = 1, nlongculvertsg
       L = abs(longculverts(ilongc)%flowlinks(1))
       L1Dlink = abs(longculverts(ilongc)%flowlinks(2))
       if (L > 0 .and. L1Dlink > 0) then
          k1   = ln(1,L) ; k2 = ln(2,L) !left and right node
          wud  = wu(L)*dx(L) !flow surface area
          wuL1  = acl(L)*wud ! 2d center factor
          wcL  (1,L ) = wuL1
          wuL2  = (1d0-acl(L))*wud
          wcL  (2, L) = wuL2
 
          !replace last addition of wcx1 etc.
          wcxy (1,k1) = wcxy (1,k1) - abs(wcx1(L))
          wcxy (2,k1) = wcxy (2,k1) - abs(wcy1(L))
 
          wcxy (1,k2) = wcxy (1,k2) - abs(wcx2(L))
          wcxy (2,k2) = wcxy (2,k2) - abs(wcy2(L))
 
          cs = lin2nodx(L1Dlink,1,csu(L1Dlink),snu(L1Dlink)) !L van buur 1D linkje
          sn = lin2nody(L1Dlink,1,csu(L1Dlink),snu(L1Dlink)) ! idem
          wcx1(L)     = cs*wuL1
          wcy1(L)     = sn*wuL1
 
          cs = lin2nodx(L1Dlink,2,csu(L1Dlink),snu(L1Dlink)) !L van buur 1D linkje
          sn = lin2nody(L1Dlink,2,csu(L1Dlink),snu(L1Dlink)) ! idem
          wcx2(L)     = cs*wuL2
          wcy2(L)     = sn*wuL2
 
          wcxy (1,k1) = wcxy (1,k1) + abs(wcx1(L))
          wcxy (2,k1) = wcxy (2,k1) + abs(wcy1(L))
 
          wcxy (1,k2) = wcxy (1,k2) + abs(wcx2(L))
          wcxy (2,k2) = wcxy (2,k2) + abs(wcy2(L))
       end if
 
       L = abs(longculverts(ilongc)%flowlinks(longculverts(ilongc)%numlinks))
       L1Dlink = abs(longculverts(ilongc)%flowlinks(longculverts(ilongc)%numlinks-1))
       if (L > 0 .and. L1Dlink > 0) then
          k1   = ln(1,L) ; k2 = ln(2,L) !left and right node
          wud  = wu(L)*dx(L) !flow surface area
          wuL1  = acl(L)*wud ! 2d center factor
          wcL  (1,L ) = wuL1
          wuL2  = (1d0-acl(L))*wud
          wcL  (2, L) = wuL2
 
          !replace last addition of wcx1 etc.
          wcxy (1,k1) = wcxy (1,k1) - abs(wcx1(L))
          wcxy (2,k1) = wcxy (2,k1) - abs(wcy1(L))
 
          wcxy (1,k2) = wcxy (1,k2) - abs(wcx2(L))
          wcxy (2,k2) = wcxy (2,k2) - abs(wcy2(L))
 
          cs = lin2nodx(L1Dlink,1,csu(L1Dlink),snu(L1Dlink)) !L van buur 1D linkje
          sn = lin2nody(L1Dlink,1,csu(L1Dlink),snu(L1Dlink)) ! idem
          wcx1(L)     = cs*wuL1
          wcy1(L)     = sn*wuL1
 
          cs = lin2nodx(L1Dlink,2,csu(L1Dlink),snu(L1Dlink)) !L van buur 1D linkje
          sn = lin2nody(L1Dlink,2,csu(L1Dlink),snu(L1Dlink)) ! idem
          wcx2(L)     = cs*wuL2
          wcy2(L)     = sn*wuL2
 
          wcxy (1,k1) = wcxy (1,k1) + abs(wcx1(L))
          wcxy (2,k1) = wcxy (2,k1) + abs(wcy1(L))
 
          wcxy (1,k2) = wcxy (1,k2) + abs(wcx2(L))
          wcxy (2,k2) = wcxy (2,k2) + abs(wcy2(L))
       end if
    enddo
 endif

 do L = 1, lnx1d
    k1 = ln(1,L) ; k2 = ln(2,L)
    wcx1(L) = wcx1(L)/a1(k1)
    wcy1(L) = wcy1(L)/a1(k1)
    wcx2(L) = wcx2(L)/a1(k2)
    wcy2(L) = wcy2(L)/a1(k2)
    if (wc(k1) > 0d0) wcL(1,L) = wcL(1,L) / wc(k1)
    if (wc(k2) > 0d0) wcL(2,L) = wcL(2,L) / wc(k2)

 enddo

 !boundaries
 do L = lnxi + 1, lnx1db
    k1 = ln(1,L) ; k2 = ln(2,L)
    wcx1(L) = wcx1(L)/a1(k1)
    wcy1(L) = wcy1(L)/a1(k1)
    wcx2(L) = wcx2(L)/a1(k2)
    wcy2(L) = wcy2(L)/a1(k2)
    if (wc(k1) > 0d0) wcL(1,L) = wcL(1,L) / wc(k1)
    if (wc(k2) > 0d0) wcL(2,L) = wcL(2,L) / wc(k2)
 enddo

 end subroutine setlinktocenterweights1d
