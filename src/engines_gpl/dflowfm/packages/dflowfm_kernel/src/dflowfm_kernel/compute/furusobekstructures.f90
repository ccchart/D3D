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

subroutine furusobekstructures()
use m_flow
use m_flowgeom
use m_strucs
implicit none
integer :: ng, n, L, Ls, LL, Lb, Lt
double precision :: zup, bup, a, fac

logical :: firstiter=.true. , jarea= .false.

firstiter = .true.
jarea     = .false.

do ng = 1, ncgensg      ! loop over generalstruc signals, sethu
   do n  = L1cgensg(ng), L2cgensg(ng)
      L  = kcgen(3,n)
      if (kcgen(1,n) == ln(2,L)) then
         Ls = -L ! Flow link has opposite orientation to structure's orientation.
      else
         Ls = L
      end if

      if (hu(L) > 0d0) then ! hu is above lowest sill
         call flgsfm( n, ng, Ls, firstiter , jarea )
      endif
      if (kmx > 0) then
         call getLbotLtop(L,Lb,Lt)
         do LL = Lb, Lt
            fu(LL) = fu(L) ; ru(LL) = ru(L)
            au(LL) = au(L)*( hu(LL)-hu(LL-1) ) / ( hu(Lt)-hu(Lb-1) )
         enddo
      endif
   enddo
enddo

end subroutine furusobekstructures
