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

      SUBROUTINE SWITCH(X, Y, mmax, nmax, JN, NUMPJ)
!      USE DIMENS
      implicit none
      integer :: mmax, nmax, jn, numpj
      double precision :: X(MMAX,NMAX),   Y(MMAX,NMAX)

       integer :: j
       double precision :: xh, yh

      DO 10 J = 1,NUMPJ/2
         XH              = X(JN,J)
         X(JN,J)         = X(JN,NUMPJ-J+1)
         X(JN,NUMPJ-J+1) = XH
         YH              = Y(JN,J)
         Y(JN,J)         = Y(JN,NUMPJ-J+1)
         Y(JN,NUMPJ-J+1) = YH
    10 CONTINUE
      RETURN
      END subroutine switch
