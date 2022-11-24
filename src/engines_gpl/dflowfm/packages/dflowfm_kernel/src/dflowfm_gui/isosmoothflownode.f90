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

 subroutine isosmoothflownode(k) ! smooth isolines in flow cells
 use m_flowgeom
 use m_flow
 use m_netw
 implicit none
 integer :: k

 integer          :: nn4, n
 double precision :: zz(10)

 nn4 = size(nd(k)%nod)
 do n = 1, nn4
    zz(n) = rnod( nd(k)%nod(n) )
 enddo
 nn4 = min(nn4, size(nd(k)%x) )
 call isofil(nd(k)%x, nd(k)%y, zz, nn4, 0)
 !call isocel(nd(k)%x, nd(k)%y, zz, nn4, 0)
 end subroutine isosmoothflownode
