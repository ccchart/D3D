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

subroutine update_dambreak_breach(startTime, deltaTime)

   use m_flowgeom
   use m_flow
   use m_missing
   use m_structures
   use unstruc_channel_flow
   use m_Dambreak
   use m_partitioninfo
   use m_meteo
   use m_flowexternalforcings
   use m_flowtimes
   use parallel_dambreaks, only: mpi_dambreak_exchange_max_avg

   implicit none

   !in-out
   double precision, intent(in)          :: startTime
   double precision, intent(in)          :: deltaTime

   !locals
   double precision                      :: tempValue, smax, smin, hmx, hmn
   integer                               :: indAverageUpStream(ndambreak)
   integer                               :: indAverageDownStream(ndambreak)
   integer                               :: nAverageUpStream, nAverageDownStream
   integer                               :: n, ierr, istru, indexLevelsAndWidths
   integer                               :: i
   integer                               :: k
   integer                               :: l
   integer                               :: Lf
   double precision                      :: e   !< energy height (m)
   double precision                      :: h   !< ratio of water depth over dambreak weir height (-)
   double precision                      :: s   !< water level (m)
   double precision                      :: u   !< flow velocity magnitude (m/s)
   double precision                      :: zb  !< bed level at dambreak (m)
   double precision                      :: zc  !< dambreak weir crest level (m)

   if (ndambreak > 0) then

      !
      ! Initialize
      !
      dambreakAveraging              = 0.0d0
      dambreakMaximum                = -huge(1.0d0)
      waterLevelsDambreakUpStream    = 0.0d0
      waterLevelsDambreakDownStream  = 0.0d0
      normalVelocityDambreak         = 0.0d0
      breachWidthDerivativeDambreak  = 0.0d0
      waterLevelJumpDambreak         = 0.0d0

      !
      ! Water level upstream of breach
      !
      if (nDambreakLocationsUpstream > 0) then
         dambreakAveraging(1,dambreakLocationsUpstreamMapping(1:nDambreakLocationsUpstream), IDB_S1U) = s1(dambreakLocationsUpstream(1:nDambreakLocationsUpstream))
         dambreakAveraging(2,dambreakLocationsUpstreamMapping, IDB_S1U) = 1d0
      endif
      do i = 1, nDambreakAveragingUpstream
         n = dambreakAveragingUpstreamMapping(i)
         do l = L1dambreaksg(n), L2dambreaksg(n)
            Lf = abs(kdambreak(3,l))
            if (Lf == 0) cycle
            if (hu(Lf) > dmiss .and. activeDambreakLinks(l) > 0 .and. (.not. db_ghost(l))) then
               dambreakAveraging(1,n, IDB_S1U) = dambreakAveraging(1,n, IDB_S1U) + wu(Lf) * s1(kdambreak(1,l))
               dambreakAveraging(2,n, IDB_S1U) = dambreakAveraging(2,n, IDB_S1U) + wu(Lf)
            endif
         enddo
         if (LStartBreach(n) > 0) then
            if ( kdambreak(1,LStartBreach(n)) > 0) then
               dambreakAveraging(3,n, IDB_S1U) = s1(kdambreak(1,LStartBreach(n)))
            endif
         endif
      enddo
      !
      ! Water level downstream of breach
      !
      if (nDambreakLocationsDownstream > 0) then
         dambreakAveraging(1,dambreakLocationsDownstreamMapping(1:nDambreakLocationsDownstream), IDB_S1D) = s1(dambreakLocationsDownstream(1:nDambreakLocationsDownstream))
         dambreakAveraging(2,dambreakLocationsDownstreamMapping, IDB_S1D) = 1d0
      endif
      do i = 1, nDambreakAveragingDownstream
         n = dambreakAveragingDownstreamMapping(i)
         do l = L1dambreaksg(n), L2dambreaksg(n)
            Lf = abs(kdambreak(3,l))
            if (Lf == 0) cycle
            if (hu(Lf) > dmiss .and. activeDambreakLinks(l) > 0 .and. (.not. db_ghost(l))) then
               dambreakAveraging(1,n, IDB_S1D) = dambreakAveraging(1,n, IDB_S1D) + wu(Lf) * s1(kdambreak(2,l))
               dambreakAveraging(2,n, IDB_S1D) = dambreakAveraging(2,n, IDB_S1D) + wu(Lf)
            endif
         enddo
         if (LStartBreach(n) > 0) then
            if (kdambreak(2,LStartBreach(n)) > 0) then
               dambreakAveraging(3,n, IDB_S1D) = s1(kdambreak(2,LStartBreach(n)))
            endif
         endif
      enddo
      !
      ! Velocity through breach
      !
      do n = 1, ndambreaksg
         do l = L1dambreaksg(n), L2dambreaksg(n)
            Lf = abs(kdambreak(3,l))
            if (Lf == 0) cycle
            if (hu(Lf) > dmiss .and. activeDambreakLinks(l) > 0 .and. (.not. db_ghost(l))) then
               dambreakAveraging(1,n, IDB_U1) = dambreakAveraging(1,n, IDB_U1) + au(Lf) * u1(Lf)
               dambreakAveraging(2,n, IDB_U1) = dambreakAveraging(2,n, IDB_U1) + au(Lf)
            endif
         enddo
      enddo
      !
      ! Maximum water level before or after the breach
      !
      do i = 1, nDambreakLocationsUpstream
         n = dambreakLocationsUpstreamMapping(i)
         ! use crest level zc and bed level zb at LStartBreach
         zc = dambreakCrestLevel(n)
         zb = dambreakUpstreamBedLevel(n)
         if (zc > zb) then
            s = s1(dambreakLocationsUpstream(i))
            u = ucmag(dambreakLocationsUpstream(i))
            e = s + 1.25*u**2/(2*ag)
            h = (s - zb)/(zc - zb)
            dambreakMaximum(n, ST_FC_WATERLEVEL) = s
            dambreakMaximum(n, ST_FC_VELMAG    ) = u
            dambreakMaximum(n, ST_FC_ENERGYHGHT) = e
            dambreakMaximum(n, ST_FC_RELDEPTH  ) = h
            dambreakMaximum(n, ST_FC_BEDLEVEL  ) = zb ! use max bed level along dam to determine failure crest level; subpar if dambreak spans multiple flow links
         endif
      enddo
      do i = 1, nDambreakLocationsDownstream
         n = dambreakLocationsDownstreamMapping(i)
         ! use crest level zc and bed level zb at LStartBreach
         zc = dambreakCrestLevel(n)
         zb = dambreakDownstreamBedLevel(n)
         if (zc > zb) then
            s = s1(dambreakLocationsDownstream(i))
            u = ucmag(dambreakLocationsDownstream(i))
            e = s + 1.25*u**2/(2*ag)
            h = (s - zb)/(zc - zb)
            dambreakMaximum(n, ST_FC_WATERLEVEL) = max(dambreakMaximum(n, ST_FC_WATERLEVEL), s)
            dambreakMaximum(n, ST_FC_VELMAG    ) = max(dambreakMaximum(n, ST_FC_VELMAG    ), u)
            dambreakMaximum(n, ST_FC_ENERGYHGHT) = max(dambreakMaximum(n, ST_FC_ENERGYHGHT), e)
            dambreakMaximum(n, ST_FC_RELDEPTH  ) = max(dambreakMaximum(n, ST_FC_RELDEPTH)  , h)
            dambreakMaximum(n, ST_FC_BEDLEVEL  ) = max(dambreakMaximum(n, ST_FC_BEDLEVEL  ), zb) ! use max bed level along dam to determine failure crest level
         endif
      enddo
      do i = 1, nDambreakAveragingUpstream
         n = dambreakAveragingUpstreamMapping(i)
         do k = L1dambreaksg(n), L2dambreaksg(n)
            if (.not. db_ghost(k)) then
               L = kdambreak(3,k)
               Lf = iabs(L)
               ! use local crest level zc and bed level zb
               if (L > 0) then
                  zc = bob(1,Lf)
                  zb = bob0(1,Lf)
               else
                  zc = bob(2,Lf)
                  zb = bob0(2,Lf)
               endif
               if (zc > zb) then
                  s = s1(kdambreak(1,k))
                  u = ucmag(kdambreak(1,k))
                  e = s + 1.25*u**2/(2*ag)
                  h = (s - zb)/(zc - zb)
                  dambreakMaximum(n, ST_FC_WATERLEVEL) = max(dambreakMaximum(n, ST_FC_WATERLEVEL), s)
                  dambreakMaximum(n, ST_FC_VELMAG    ) = max(dambreakMaximum(n, ST_FC_VELMAG    ), u)
                  dambreakMaximum(n, ST_FC_ENERGYHGHT) = max(dambreakMaximum(n, ST_FC_ENERGYHGHT), e)
                  dambreakMaximum(n, ST_FC_RELDEPTH  ) = max(dambreakMaximum(n, ST_FC_RELDEPTH  ), h)
                  dambreakMaximum(n, ST_FC_BEDLEVEL  ) = max(dambreakMaximum(n, ST_FC_BEDLEVEL  ), zb) ! use max bed level along dam to determine failure crest level
               endif
            endif
         enddo
      enddo
      do i = 1, nDambreakAveragingDownstream
         n = dambreakAveragingDownstreamMapping(i)
         do k = L1dambreaksg(n), L2dambreaksg(n)
            if (.not. db_ghost(k)) then
               L = kdambreak(3,k)
               Lf = iabs(L)
               ! use local crest level zc and bed level zb
               if (L > 0) then
                  zc = bob(2,Lf)
                  zb = bob0(2,Lf)
               else
                  zc = bob(1,Lf)
                  zb = bob0(1,Lf)
               endif
               if (zc > zb) then
                  s = s1(kdambreak(2,k))
                  u = ucmag(kdambreak(2,k))
                  e = s + 1.25*u**2/(2*ag)
                  h = (s - zb)/(zc - zb)
                  dambreakMaximum(n, ST_FC_WATERLEVEL) = max(dambreakMaximum(n, ST_FC_WATERLEVEL), s)
                  dambreakMaximum(n, ST_FC_VELMAG    ) = max(dambreakMaximum(n, ST_FC_VELMAG    ), u)
                  dambreakMaximum(n, ST_FC_ENERGYHGHT) = max(dambreakMaximum(n, ST_FC_ENERGYHGHT), e)
                  dambreakMaximum(n, ST_FC_RELDEPTH  ) = max(dambreakMaximum(n, ST_FC_RELDEPTH  ), h)
                  dambreakMaximum(n, ST_FC_BEDLEVEL  ) = max(dambreakMaximum(n, ST_FC_BEDLEVEL  ), zb) ! use max bed level along dam to determine failure crest level
               endif
            endif
         enddo
      enddo

      call mpi_dambreak_exchange_max_avg()

      !
      ! Water level upstream of breach
      !
      do n = 1, ndambreaksg
         if (dambreakAveraging(2,n, IDB_S1U) > 0.0d0) then
            waterLevelsDambreakUpStream(n)  = dambreakAveraging(1,n, IDB_S1U)/dambreakAveraging(2,n, IDB_S1U)
         else if (startTime >= network%sts%struct(dambreaks(n))%dambreak%t0) then
            waterLevelsDambreakUpStream(n) = dambreakAveraging(3,n, IDB_S1U)
         else
            continue
         endif
      enddo
      !
      ! Water level downstream of breach
      !
      do n = 1, ndambreaksg
         if (dambreakAveraging(2,n, IDB_S1D) > 0.0d0) then
            waterLevelsDambreakDownStream(n)  = dambreakAveraging(1,n, IDB_S1D)/dambreakAveraging(2,n, IDB_S1D)
         else if (startTime >= network%sts%struct(dambreaks(n))%dambreak%t0) then
            waterLevelsDambreakDownStream(n) = dambreakAveraging(3,n, IDB_S1D)
         else
            continue
         endif
      enddo
      !
      ! Velocity through breach
      !
      do n = 1, ndambreaksg
         if (dambreakAveraging(2,n, IDB_U1) > 0.0d0) then
            normalVelocityDambreak(n)  = dambreakAveraging(1,n, IDB_U1)/dambreakAveraging(2,n, IDB_U1)
         endif
      enddo

      ! Compute dambreak widths
      do n = 1, ndambreaksg
         istru = dambreaks(n)
         if (istru.ne.0) then
            if (network%sts%struct(istru)%dambreak%algorithm == ST_DB_VDKNAAP_00 .or. network%sts%struct(istru)%dambreak%algorithm == ST_DB_VERHEY_VDKNAAP_02) then
               ! Compute the breach width
               call prepareComputeDambreak(network%sts%struct(istru)%dambreak, waterLevelsDambreakUpStream(n), waterLevelsDambreakDownStream(n), normalVelocityDambreak(n), startTime, deltaTime, maximumDambreakWidths(n))

            else if (network%sts%struct(istru)%dambreak%algorithm == ST_DB_FRAGCURVE) then ! fragility curve
               if (dambreakMaximum(n, network%sts%struct(istru)%dambreak%failQuantity) > network%sts%struct(istru)%dambreak%failValue) then
                  ! if breaching condition is satisfied, break width completely, but lower structure by failFraction * (initial crest level - bed level)
                  ! dambreakCrestLevel = dambreakInitialCrestLevel - network%sts%struct(istru)%dambreak%failFraction * max(0d0, dambreakInitialCrestLevel - dambreakBedLevel)
                  network%sts%struct(istru)%dambreak%crl   = network%sts%struct(istru)%dambreak%crestLevelIni - network%sts%struct(istru)%dambreak%failFraction * max(0d0, network%sts%struct(istru)%dambreak%crestLevelIni - dambreakMaximum(n, ST_FC_BEDLEVEL))
                  network%sts%struct(istru)%dambreak%width = maximumDambreakWidths(n)
               endif

            else if(network%sts%struct(istru)%dambreak%algorithm == ST_DB_PRESCRIBED .and. startTime > network%sts%struct(istru)%dambreak%t0) then
               ! Time in the tim file is relative to the start time
               success = ec_gettimespacevalue_by_itemID(ecInstancePtr, item_dambreakLevelsAndWidthsFromTable, irefdate, tzone, tunit, startTime-network%sts%struct(istru)%dambreak%t0)
               ! NOTE: AvD: the code above works correctly, but is dangerous:
               ! the addtimespace for dambreak has added each dambreak separately with a targetoffset.
               ! The gettimespace above, however, gets the values for *all* dambreaks, but with the relative time
               ! of the *current* dambreak #n.
               ! This means that if t0 values for all dambreaks are different, then the dambreakLevelsAndWidthsFromTable(1:n-1) have become obsolete now.
               ! It works, because in the previous loop iterations the values that were then still correct
               ! have already been set into the %crl and %width values.
               if (success)  then
                  indexLevelsAndWidths = (n - 1) * 2 + 1
                  network%sts%struct(istru)%dambreak%crl   = dambreakLevelsAndWidthsFromTable(indexLevelsAndWidths)
                  network%sts%struct(istru)%dambreak%width = dambreakLevelsAndWidthsFromTable(indexLevelsAndWidths + 1 )
               else
                   return
               endif
            endif
            ! Store breach width derivatives
            tempValue = network%sts%struct(istru)%dambreak%breachWidthDerivative
            if (tempValue>0) then
               breachWidthDerivativeDambreak(n) = tempValue
            else
               breachWidthDerivativeDambreak(n) = &
                  (network%sts%struct(istru)%dambreak%width - breachWidthDambreak(n)) / deltaTime
            endif

            ! Store the current dambreak width
            breachWidthDambreak(n) = network%sts%struct(istru)%dambreak%width
            ! Store the current dambreak crest level
            breachDepthDambreak(n) = network%sts%struct(istru)%dambreak%crl

            ! Store water level jump
            tempValue = network%sts%struct(istru)%dambreak%waterLevelJumpDambreak
            if (tempValue>0) then
               ! Algo 1 or 2: from prepareComputeDambreak
               waterLevelJumpDambreak(n) = tempValue
            else
               ! Algo 3 (timeseries), compute here:
               smax = max(waterLevelsDambreakUpStream(n), waterLevelsDambreakDownStream(n))
               smin = min(waterLevelsDambreakUpStream(n), waterLevelsDambreakDownStream(n))
               hmx  = max(0d0,smax - network%sts%struct(istru)%dambreak%crl)
               hmn  = max(0d0,smin - network%sts%struct(istru)%dambreak%crl)
               waterLevelJumpDambreak(n) = hmx - hmn
            endif

         endif
      enddo
   endif
end subroutine update_dambreak_breach
