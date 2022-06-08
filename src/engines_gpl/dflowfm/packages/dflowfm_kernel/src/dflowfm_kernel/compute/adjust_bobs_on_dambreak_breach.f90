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

 !> Calculate the links affected by the dam break and sets bobs accordingly
subroutine adjust_bobs_on_dambreak_breach(width, maxwidth, crl, startingLink, L1, L2, strucid)

   use m_flowgeom
   use m_flowexternalforcings
   use MessageHandling

   implicit none

   !input
   double precision, intent(in) :: width, crl, maxwidth
   integer, intent(in)          :: startingLink, L1, L2
   character(len=*), intent(in) :: strucid
   !local variables
   integer                      :: k, Lf
   double precision             :: leftBreachWidth, rightBreachWidth
   double precision             :: leftside, rightside
   double precision             :: remainder, leftfrac

   !nothing is open
   if (width<=0) return

   !something is open
   Lf = iabs(kdambreak(3,startingLink))
   if (Lf > 0) then
      bob(1,Lf) = max(bob0(1, Lf), crl)
      bob(2,Lf) = max(bob0(2, Lf), crl)
      activeDambreakLinks(startingLink) = 1
   endif
   if ((width - dambreakLinksEffectiveLength(startingLink))<= 0) then
      ! breach width still less than width of starting link
      dambreakLinksBreachLength(startingLink) = width
      return
   else
      !breach width larger than width of starting link
      dambreakLinksBreachLength(startingLink) = dambreakLinksEffectiveLength(startingLink)
      leftside  = sum(dambreakLinksEffectiveLength(L1:startingLink-1))
      rightside = sum(dambreakLinksEffectiveLength(startingLink+1:L2))
      if (dambreakWidening == DBW_SYMM) then
         ! original code
         leftBreachWidth = (width - dambreakLinksEffectiveLength(startingLink))/2.0d0
         rightBreachWidth = leftBreachWidth
      elseif (dambreakWidening == DBW_PROP) then
         ! proportional
         remainder = width - dambreakLinksEffectiveLength(startingLink)
         leftfrac  = leftside / (leftside +  rightside)
         leftBreachWidth = leftfrac * remainder
         rightBreachWidth = (1.0d0 - leftfrac) * remainder
      elseif (dambreakWidening == DBW_SYMM_ASYMM) then
         ! first symmetric, then asymmetric
         remainder = width - dambreakLinksEffectiveLength(startingLink)
         if (remainder < 2 * min(leftside,rightside)) then
            leftBreachWidth = remainder / 2.0d0
            rightBreachWidth = (1.0d0 - leftfrac) * remainder
         elseif (leftside <= rightside) then
            leftBreachWidth = leftside
            rightBreachWidth = remainder - leftside
         else
            rightBreachWidth = rightside
            leftBreachWidth = remainder - rightside
         endif
      endif
      do k = startingLink - 1, L1, -1
         Lf = iabs(kdambreak(3,k))
         if (Lf > 0) then
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
         endif
         activeDambreakLinks(k) = 1
         if (leftBreachWidth>=dambreakLinksEffectiveLength(k)) then
            dambreakLinksBreachLength(k) = dambreakLinksEffectiveLength(k)
            leftBreachWidth = leftBreachWidth - dambreakLinksEffectiveLength(k)
         else
            dambreakLinksBreachLength(k) = leftBreachWidth
            leftBreachWidth = 0d0
            exit
         endif
      enddo
      !right from the breach point
      do k = startingLink + 1, L2
         Lf = iabs(kdambreak(3,k))
         if (Lf > 0) then
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
         endif
         activeDambreakLinks(k) = 1
         if (rightBreachWidth>=dambreakLinksEffectiveLength(k)) then
            dambreakLinksBreachLength(k) = dambreakLinksEffectiveLength(k)
            rightBreachWidth = rightBreachWidth - dambreakLinksEffectiveLength(k)
         else
            dambreakLinksBreachLength(k) = rightBreachWidth
            rightBreachWidth = 0d0
            exit
         endif
      enddo
      if (leftBreachWidth > 1.0d-6 * maxwidth .or. rightBreachWidth > 1.0d-6 * maxwidth) then
         write (msgbuf, '(3a)' ) 'The breach width of dam ''', trim(strucid), ''' is wider than the actual dam width.'
         call SetMessage(LEVEL_WARN, msgbuf)
      end if
   endif

end subroutine adjust_bobs_on_dambreak_breach
