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

! $Id$
! $HeadURL$

 subroutine lineinterp(xx, yy, ktx, x,y,n)
 implicit none
 integer          :: ktx, n, k, ip
 double precision :: xx(ktx), yy(ktx), x(n), y(n)
 double precision :: a, b


 ip = 1
 do k = 1, ktx
    do while ( xx(k) > x(ip+1) .and. ip < n-1 )
       ip = ip + 1
    enddo
    if ( xx(k) <= x(ip) ) then
       yy(k)  =  y(ip)
    else if ( xx(k) > x(ip) .and. xx(k) <= x(ip+1) ) then
       a     = ( xx(k) - x(ip) ) / max( 1d-4 , x(ip+1) - x(ip) ) ; b = 1d0 - a
       yy(k) = b*y(ip) + a*y(ip+1)
    else
       yy(k) = y(ip+1)
    endif
 enddo

 end subroutine lineinterp
