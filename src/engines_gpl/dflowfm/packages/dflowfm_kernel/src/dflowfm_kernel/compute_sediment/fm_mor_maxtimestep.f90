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

   subroutine fm_mor_maxtimestep()
   use m_flowtimes
   use m_flow
   use m_flowgeom
   use m_sediment
   use m_partitioninfo
   use m_fm_erosed, only: sxtot, sytot, cdryb, morfac, lsedtot, ndx_mor, nd_mor, ln_mor

   implicit none

   integer           :: k, k1, k2, kk, L, ised, ac1, ac2
   double precision  :: dum, sx, sy, sL, dt, dtmaxmor, dhmax

   dtmaxmor = huge(0d0)

   do k = 1, ndx_mor
      dum = 0.d0
      if (kcsmor(k)==0) then
         cycle
      endif

      do ised = 1, lsedtot
         !
         do kk = 1, nd_mor(k)%lnx
            L = iabs(nd_mor(k)%ln(kk))
            k1 = ln_mor(1,L)
            k2 = ln_mor(2,L)
            ac1 = acl(L)
            ac2 = 1d0-ac1

            sx = (ac1*sxtot(k1,ised) + ac2*sxtot(k2,ised))/cdryb(ised)*max(morfac,1d0)
            sy = (ac1*sytot(k1,ised) + ac2*sytot(k2,ised))/cdryb(ised)*max(morfac,1d0)
            sL = csu(L)*sx + snu(L)*sy

            if (k2 .eq. k) sL = -sL

            if (sL .ge. 0.) then        ! outgoing transport fluxes only
               dum = dum + sL*wu(L)
            end if
         end do
      end do
      if (dum > tiny(0d0)) then
         dt = dzbdtmax*ba(k) / max(dum,eps10)   ! safety
         if ( dt.lt.dtmaxmor ) then
            dtmaxmor = dt
         end if
      end if
   end do

   if ( jampi.eq.1 ) then
      call reduce_double_min(dtmaxmor)
   end if

   if (dtmaxmor > dts) dtmaxmor = dts
   dtmaxmor = dts/ceiling(dts/dtmaxmor)
   !
   dts = dtmaxmor
   dti = 1d0/dts

   end subroutine fm_mor_maxtimestep
