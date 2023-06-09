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

      SUBROUTINE WRIARCsam(MARC,DP,MMAX,NMAX,MC,NC,X0,Y0,DX,DY,dmiss)
      implicit none
      double precision :: dmiss
      double precision :: dp
      double precision :: dx, dy
      integer :: i
      integer :: j
      integer :: marc
      integer :: mc
      integer :: mmax
      integer :: nc
      integer :: nmax
      double precision :: x0
      double precision :: y0
!      DIMENSION DP  (MMAX,NMAX)
      DIMENSION DP  (NMAX,MMAX)   ! SPvdP: j-index is fastest running in sample arrays
 
      CALL WRITEARCINFOHEADER(MARC,MC,NC,X0,Y0,DX,DY,dmiss)
   
!      DO 10 J = NC,1,-1      ! SPvdP: j-index is already reverse-order in sample arrays
      DO 10 J = 1,NC
!         WRITE(MARC,'(5000F10.2)') ( DP(I,J),I = 1,MC)
      
         WRITE(MARC,'(99999999F10.2)') ( DP(J,I),I = 1,MC)   ! SPvdP: j-index is fastest running in sample arrays
   10 CONTINUE

      RETURN
      END
