!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2023.                                
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

! 
! 

 subroutine setpressurehull()
 use m_ship
 use m_flowgeom
 use m_flow
 implicit none
 integer :: L, LL, k1, k2
 double precision :: r1,r2
 do LL = 1,Lnx
    k1 = ln(1,LL) ; k2 = ln(2,LL)
    if (zsp(k1) .ne. 0d0 .or. zsp(k2) .ne. 0d0) then
       do L = Lbot(LL), Ltop(LL)
          adve(L) = adve(L) + ag*( zsp(k2) - zsp(k1) )*dxi(LL)          ! impose ship hull
       enddo
    endif
 enddo
 end subroutine setpressurehull
