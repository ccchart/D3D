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

      SUBROUTINE MAKEF1D(XR,YR,MNMAX) ! naar rekenvlak SUBROUTINE MAKEF1D
      USE M_SFERIC
      USE M_MISSING
      implicit none
      integer :: mnmax
      DOUBLE PRECISION :: XR(MNMAX), YR(MNMAX), FI2
      integer :: i

      DO I = 1,MNMAX
         IF (XR(I) .NE. DXYMIS) THEN
            FI2     = DG2RD*YR(I)
            YR(I) = ( 1D0 + SIN(FI2) ) / COS(FI2)
            YR(I) = LOG(YR(I))
            XR(I) = DG2RD*XR(I)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MAKEF1D
