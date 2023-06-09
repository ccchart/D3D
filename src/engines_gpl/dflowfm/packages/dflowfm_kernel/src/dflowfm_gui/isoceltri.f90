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

      SUBROUTINE ISOCELTRI(X,Y,P,NCOLR)
      implicit none
      integer          :: i, ih, ja, jaauto, ncolr, ncols, nh, nie, nis, nplus, nv
      double precision :: dv, p, p1, p2, val, vmax, vmin, vn, x, x1, x2, xh, xhit, y, y1, y2, yh, yhit


!     TEKENT ALLE NV ISOLIJNEN IN EEN CEL TEKAL-METHODE
      DIMENSION P(3),X(3),Y(3),XH(3),YH(3)
      COMMON /DEPMAX/ VMAX,VMIN,DV,VAL(256),NCOLS(256),NV,NIS,NIE,JAAUTO

      DO 10 I = 1,NV
         NPLUS = 1
         VN    = VAL(I)
         NH    = 0
         DO 20 IH = 1,3
            IF (IH .EQ. 3) NPLUS = -2
            P1 = P(IH)
            P2 = P(IH + NPLUS)
            X1 = X(IH)
            X2 = X(IH + NPLUS)
            Y1 = Y(IH)
            Y2 = Y(IH + NPLUS)
            CALL HITLIN(P1,P2,X1,Y1,X2,Y2,VN,XHIT,YHIT,JA)
            IF (JA .EQ. 1) THEN
               NH     = NH + 1
               XH(NH) = XHIT
               YH(NH) = YHIT
            ENDIF
   20    CONTINUE
!        IF (NH .GT. 1) CALL DISPF2(XH,YH,NH,3,NCOLS(I+1))
         IF (NH .GT. 1) CALL DISPF2(XH,YH,NH,3,0)
   10 CONTINUE

      IF (NCOLR .NE. 0) CALL DISPF2(X,Y,3,3,NCOLR)

      RETURN
      END
