!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
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

      subroutine DISPFORMscale0(value,fmt)
      implicit none
      integer :: n1
      integer :: n2
      integer :: ndec
      double precision :: scalesize
      double precision :: value
      double precision :: xsc
      double precision :: ysc
      character fmt*(*)
      COMMON /SCALEPOS/ XSC,YSC,SCALESIZE,NDEC

      fmt='(f9.3)'

      if (value .eq. 0.0) then
         fmt='(f3.1)'
         return
      endif

      n1 = int(log10(abs(value)))

      if (n1 .le. 6 .and. n1 .gt. 0) then
         n2 = min(9,n1 + 3)
         write (fmt(5:5),'(i1)') 9 - n2
      else if (n1 .ge. -5 .and. n1 .lt. 0) then
         write (fmt(5:5),'(i1)') 6
      else if ( n1 .eq. 0) then
         write (fmt(5:5),'(i1)') 6
      else
         fmt ='(e9.3)'
      endif
      IF (NDEC .GT. 0) write (fmt(5:5),'(i1)') NDEC
      return
      end
