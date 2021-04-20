 !> Calculate the links affected by the dam break and sets bobs accordingly
subroutine adjust_bobs_on_dambreak_breach(width, crl, startingLink, L1, L2, strucid)

   use m_flowgeom
   use m_flowexternalforcings
   use MessageHandling

   implicit none

   !input
   double precision, intent(in) :: width, crl
   integer, intent(in)          :: startingLink, L1, L2
   character(len=*), intent(in) :: strucid
   !local variables
   integer                      :: k, Lf
   double precision             :: leftBreachWidth, rightBreachWidth

   !nothing is open
   if (width<=0) return;

   !something is open
   Lf = iabs(kdambreak(3,startingLink))
   bob(1,Lf) = max(bob0(1, Lf), crl)
   bob(2,Lf) = max(bob0(2, Lf), crl)
   activeDambreakLinks(startingLink) = 1
   if ((width - dambreakLinksEffectiveLength(startingLink))<= 0) then
      dambreakLinksActualLength(startingLink) = width
      return
   else
      !left from the breach point: the breach width is larger
      dambreakLinksActualLength(startingLink) = dambreakLinksEffectiveLength(startingLink)
      leftBreachWidth = (width - dambreakLinksEffectiveLength(startingLink))/2.0d0
      rightBreachWidth = leftBreachWidth
        do k = startingLink - 1, L1, -1
         Lf = iabs(kdambreak(3,k))
         if (leftBreachWidth>=dambreakLinksEffectiveLength(k)) then
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
            activeDambreakLinks(k) = 1
            dambreakLinksActualLength(k) = dambreakLinksEffectiveLength(k)
            leftBreachWidth = leftBreachWidth - dambreakLinksEffectiveLength(k)
         else
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
            activeDambreakLinks(k) = 1
            dambreakLinksActualLength(k) = leftBreachWidth
            leftBreachWidth = 0d0
            exit
         endif
      enddo
      !right from the breach point
      do k = startingLink + 1, L2
         Lf = iabs(kdambreak(3,k))
         if (rightBreachWidth>=dambreakLinksEffectiveLength(k)) then
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
            activeDambreakLinks(k) = 1
            dambreakLinksActualLength(k) = dambreakLinksEffectiveLength(k)
            rightBreachWidth = rightBreachWidth - dambreakLinksEffectiveLength(k)
         else
            bob(1,Lf) = max(bob0(1, Lf), crl)
            bob(2,Lf) = max(bob0(2, Lf), crl)
            activeDambreakLinks(k) = 1
            dambreakLinksActualLength(k) = rightBreachWidth
            rightBreachWidth = 0d0
            exit
         endif
      enddo
      if (leftBreachWidth /= 0d0 .or. rightBreachWidth /= 0d0) then
         write (msgbuf, '(3a)' ) 'The breach width of dam ''', trim(strucid), ''' is wider than the actual dam width.'
         call SetMessage(LEVEL_WARN, msgbuf)
      end if
   endif


end subroutine adjust_bobs_on_dambreak_breach
