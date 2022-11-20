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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module parallel_dambreaks

implicit none
private

public mpi_dambreak1
public mpi_dambreak2
public mpi_dambreak3
public mpi_dambreak_exchange_max_avg

contains

!> subroutine to make dambreak information consistent across parallel partitions (call 1)
subroutine mpi_dambreak1(lftopol)
   use m_flowexternalforcings, only: ndambreaksg, L1dambreaksg, L2dambreaksg, ndambreak, kedb, db_ghost, db_ds, db_ndomains, ndambreak_glob
   use m_partitioninfo, only: jampi, ndomains, my_rank, DFM_COMM_DFMWORLD, idomain, link_ghostdata
   use m_timer, only: jatimer, starttimer, stoptimer, IMPIREDUCE
#ifdef HAVE_MPI
   use mpi, only: MPI_INTEGER, MPI_SUM, MPI_DOUBLE_PRECISION
#endif
   use sorting_algorithms, only: sort
   use m_flowgeom, only: ln
   
   implicit none
   
   integer        , allocatable, dimension(:)    :: lftopol            !< mapping of flow link to dambreak polyline segment
   
   double precision                              :: ds                 !< offset of a link along the dambreak polyline
   double precision, allocatable, dimension(:)   :: dbl_ds             !< offset of each link along the dambreak polyline (local only)
   double precision, allocatable, dimension(:)   :: dbl_ds_sorted      !< sorted offsets of each link along the dambreak polyline (one dambreak object only)
   double precision, allocatable, dimension(:)   :: dbl_ds_glob        !< offset of each link along the dambreak polyline (global collection)
   integer         , allocatable, dimension(:)   :: dbl_domain         !< domain number if link belongs to this domain, or zero if it is associated to another domain (local only)
   integer         , allocatable, dimension(:)   :: dbl_domain_glob    !< link domain number (global collection)
   integer         , allocatable, dimension(:)   :: domain_sorted      !< domain numbers for the sorted links (one dambreak object only)
   integer         , allocatable, dimension(:,:) :: ndbLinks           !< number of links per dambreak object for this partition (entered in (my_rank,:))
   integer         , allocatable, dimension(:,:) :: ndbLinks_glob      !< number of links per dambreak object per partition
   integer                                       :: idmn_link          !< domain number associated with selected link
   integer                                       :: ierr               !< communication error flag
   integer         , allocatable, dimension(:)   :: index_sorted       !< order of the links when sorted by offset (one dambreak object only)
   integer                                       :: jaghost            !< flow link is ghost link (1) or not (0)
   integer                                       :: k1                 !< index of upstream cell
   integer                                       :: k2                 !< index of downstream cell
   integer                                       :: ll                 !< flow link number
   integer                                       :: l                  !< local dambreak link index
   integer                                       :: lg                 !< global dambreak link index
   integer                                       :: my_p               !< parallel partition index for this domain (1-based)
   integer                                       :: n                  !< dambreak index
   integer                                       :: ndambreak_updated  !< updated value of ndambreak (the number of dambreak links)
   integer         , allocatable, dimension(:)   :: nlinks_glob        !< number of links per dambreak summed over all partitions
   integer         , allocatable, dimension(:)   :: nlinks_local       !< number of links per dambreak important for this partition (0 if dambreak is in one partition, but not this one; otherwise nlinks_unique)
   integer         , allocatable, dimension(:)   :: nlinks_unique      !< number of unique links per dambreak summed over all partitions
   integer                                       :: offset             !< offset in tmp_kedb array
   integer                                       :: p                  !< parallel partition index
   double precision, allocatable, dimension(:)   :: tmp_db_ds          !< temporary array for offset of links along dambreak polyline
   integer         , allocatable, dimension(:)   :: tmp_db_ghost       !< temporary array for ghost cell flags (0 if link not in this partition)
   integer         , allocatable, dimension(:)   :: tmp_kedb           !< temporary array for local links numbers (0 if link not in this partition)
   integer         , allocatable, dimension(:)   :: tmp_L1dambreaksg   !< temporary array for start index of local links numbers per dambreak object
   integer         , allocatable, dimension(:)   :: tmp_L2dambreaksg   !< temporary array for end index of local links numbers per dambreak object
   integer         , allocatable, dimension(:)   :: tmp_lftopol        !< temporary array for mapping of flow link to dambreak polyline segment
   integer                                       :: ulinks             !< counter for number of unique links associated with a particular dambreak object
   logical                                       :: mpi_needed         !< flag indicating whether MPI calls are actually needed
   
   if (.not.jampi) return
   
   allocate(ndbLinks(ndomains,ndambreaksg), ndbLinks_glob(ndomains,ndambreaksg), db_ndomains(ndambreaksg), nlinks_glob(ndambreaksg), nlinks_unique(ndambreaksg), nlinks_local(ndambreaksg))
   ndbLinks        = 0
   ndbLinks_glob   = 0
   db_ndomains     = 0
   my_p = my_rank + 1
   mpi_needed = .false.
   
   ! exchange the number of dambreak links (including ghost links) per dambreak object per partition
   do n = 1, ndambreaksg
      ndbLinks(my_p,n) = L2dambreaksg(n) - L1dambreaksg(n) + 1
   enddo
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(ndblinks,ndblinks_glob,ndomains*ndambreaksg,MPI_INTEGER,MPI_SUM,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   
   ! determine the global number of dambreak links per dambreak object (including double counting of ghost links)
   ! determine whether dambreak mpi calls are needed (at least one dambreak object partly in multiple partitions - including ghost area)
   ! determine total number of dambreak links associated with dambreak objects that are partly in multiple partitions
   ndambreak_glob = 0
   nlinks_glob = sum(ndblinks_glob,1)
   do n = 1, ndambreaksg
      do p = 1, ndomains
         if (ndblinks_glob(p,n) > 0) db_ndomains(n) = db_ndomains(n) + 1
      enddo
      if (db_ndomains(n) > 1) then
         mpi_needed = .true.
         ndambreak_glob = ndambreak_glob + nlinks_glob(n)
      endif
   enddo
   
   ! if all dambreak objects are inside one partition, we don't need to exchange any data
   ! TODO: note that this may not hold for the locations that influence the dambreak
   if (.not.mpi_needed) return
   
   ! for all dambreak objects that are in multiple partitions, exchange the distance along the defining polyline (dbl_ds)
   ! and for non ghost links, exchange the associated partition (dbl_domain)
   allocate(dbl_ds(ndambreak_glob), dbl_domain(ndambreak_glob))
   allocate(dbl_ds_glob(ndambreak_glob), dbl_domain_glob(ndambreak_glob))
   dbl_ds = 0d0
   dbl_domain = 0
   dbl_ds_glob = 0d0
   dbl_domain_glob = 0
   lg = 0
   do n = 1, ndambreaksg
      if (db_ndomains(n) <= 1) cycle
      do p = 1, my_p-1
         lg = lg + ndbLinks_glob(p,n)
      enddo
      do l = L1dambreaksg(n), L2dambreaksg(n)
         lg = lg+1
         dbl_ds(lg) = db_ds(l)
         ll = abs(kedb(l))
         k1 = ln(1,ll)
         k2 = ln(2,ll)
         call link_ghostdata(my_rank,idomain(k1),idomain(k2),jaghost,idmn_link)
         if (.not.jaghost) then
            dbl_domain(lg) = my_p
         endif
      enddo
      do p = my_p+1, ndomains
         lg = lg + ndbLinks_glob(p,n)
      enddo
   enddo
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(dbl_ds,dbl_ds_glob,ndambreak_glob,MPI_DOUBLE_PRECISION,MPI_SUM,DFM_COMM_DFMWORLD, ierr)
   call MPI_allreduce(dbl_domain,dbl_domain_glob,ndambreak_glob,MPI_INTEGER,MPI_SUM,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   
   ! sort the links per dambreak object based on their distance along the dambreak polyline
   lg = 0
   do n = 1, ndambreaksg
      if (db_ndomains(n) <= 1) then
          nlinks_unique(n) = nlinks_glob(n)
          nlinks_local(n) = ndbLinks(my_p,n)
          cycle
      endif
      !
      allocate(dbl_ds_sorted(nlinks_glob(n)), index_sorted(nlinks_glob(n)), domain_sorted(nlinks_glob(n)))
      call sort(nlinks_glob(n), dbl_ds_glob(lg+1:lg+nlinks_glob(n)), dbl_ds_sorted, index_sorted)
      domain_sorted = dbl_domain_glob(lg+index_sorted)
      dbl_ds_glob(lg+1:lg+nlinks_glob(n)) = dbl_ds_sorted
      dbl_domain_glob(lg+1:lg+nlinks_glob(n)) = domain_sorted
      !
      ulinks = 0
      do ll =1,nlinks_glob(n)
         if (domain_sorted(ll) > 0) ulinks = ulinks+1
      enddo
      nlinks_unique(n) = ulinks
      nlinks_local(n) = ulinks
      deallocate(dbl_ds_sorted, index_sorted, domain_sorted)
      !
      lg = lg + nlinks_glob(n)
   enddo
   ndambreak_updated = sum(nlinks_local)
   deallocate(nlinks_unique, ndbLinks_glob, dbl_domain, dbl_ds)
   
   ! build new versions of administration based on the gained knowledge:
   ! the administration includes the dambreaks that are only in this partition, and
   ! all dambreaks that are in multiple partitions (but not necessarily this partition)
   allocate(tmp_kedb(ndambreak_updated), tmp_db_ghost(ndambreak_updated), tmp_db_ds(ndambreak_updated), tmp_lftopol(ndambreak_updated), tmp_L1dambreaksg(ndambreaksg), tmp_L2dambreaksg(ndambreaksg))
   tmp_kedb = 0
   tmp_db_ghost = 0
   tmp_db_ds = 0d0
   tmp_lftopol = 0
   tmp_L1dambreaksg = 0
   tmp_L2dambreaksg = 0
   lg = 0
   offset = 0
   do n = 1, ndambreaksg
      tmp_L1dambreaksg(n) = offset + 1
      tmp_L2dambreaksg(n) = offset + nlinks_local(n)
      if (db_ndomains(n) <= 1) then
         ! dambreak in only one partition
         if (nlinks_local(n) > 0) then
            ! this dambreak onject is located only in this partition! Copy the links from kedb
            tmp_kedb(tmp_L1dambreaksg(n):tmp_L2dambreaksg(n)) = kedb(L1dambreaksg(n):L2dambreaksg(n))
            tmp_db_ghost(tmp_L1dambreaksg(n):tmp_L2dambreaksg(n)) = 0
            tmp_db_ds(tmp_L1dambreaksg(n):tmp_L2dambreaksg(n)) = db_ds(L1dambreaksg(n):L2dambreaksg(n))
            tmp_lftopol(tmp_L1dambreaksg(n):tmp_L2dambreaksg(n)) = lftopol(L1dambreaksg(n):L2dambreaksg(n))
            offset = offset + nlinks_local(n)
         endif
      else
         ! dambreak is located in multiple partitions (but it may not have any link inside this partition)
         do ll = 1,nlinks_glob(n)
            lg = lg+1
            if (dbl_domain_glob(lg) > 0) then
               ! this link is associated with a partition, so it wasn't a ghost link
               offset = offset + 1
               ds = dbl_ds_glob(lg)
               if (dbl_domain_glob(lg) == my_p) then
                  ! link is internal to this domain
                  tmp_db_ghost(offset) = 0
               else
                  tmp_db_ghost(offset) = 1
               endif
               tmp_db_ds(offset) = ds
               ! check if this partition has an entry at the same offset ...
               do l = L1dambreaksg(n), L2dambreaksg(n)
                  ! the next line allows for comparing doubles since ds has been copied from db_ds
                  if (ds == db_ds(l)) then
                     ! if yes ... copy the associated link number
                     tmp_kedb(offset) = kedb(l)
                     tmp_lftopol(offset) = lftopol(l)
                     exit
                  endif
               enddo
            endif
         enddo
      endif
   enddo
   deallocate(dbl_ds_glob, dbl_domain_glob, nlinks_local, nlinks_glob)
   
   ! replace the old administration
   deallocate(kedb, db_ghost, db_ds, lftopol)
   allocate(kedb(ndambreak_updated), db_ghost(ndambreak_updated), db_ds(ndambreak_updated), lftopol(ndambreak_updated))
   ndambreak = ndambreak_updated
   kedb = tmp_kedb
   db_ghost = tmp_db_ghost
   db_ds = tmp_db_ds
   L1dambreaksg = tmp_L1dambreaksg
   L2dambreaksg = tmp_L2dambreaksg
   lftopol = tmp_lftopol
   deallocate(tmp_kedb, tmp_db_ghost, tmp_db_ds, tmp_L1dambreaksg, tmp_L2dambreaksg, tmp_lftopol)
end subroutine mpi_dambreak1


!> subroutine to make dambreak information consistent across parallel partitions (call 2)
subroutine mpi_dambreak2()
   use m_flowexternalforcings, only: ndambreaksg, dsStartBreach, LStartBreach, L1dambreaksg, L2dambreaksg, db_ndomains, ndambreak_glob, dambreakLinksEffectiveLength, maximumDambreakWidths
   use m_partitioninfo, only: jampi, DFM_COMM_DFMWORLD
   use m_timer, only: jatimer, starttimer, stoptimer, IMPIREDUCE
   use mpi, only: MPI_MIN, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_INTEGER
   
   implicit none
   
   double precision, allocatable, dimension(:)   :: dblEffLen          !< dambreak link effective (=maximum) length filled with local values
   double precision, allocatable, dimension(:)   :: dblEffLen_glob     !< dambreak link effective (=maximum) length filled with all values (global)
   double precision, allocatable, dimension(:)   :: dsStartBreach_glob !< distance between specified breach point and nearest link
   integer         , allocatable, dimension(:)   :: LStartBreach_glob  !< index of dambreak link at which breach starts
   integer                                       :: ierr               !< communication error flag
   integer                                       :: l
   integer                                       :: lg
   integer                                       :: Lf
   integer                                       :: n
   
   if (.not.jampi) return
   
   allocate(dblEffLen(ndambreak_glob), dblEffLen_glob(ndambreak_glob))
   dblEffLen = 0d0
   dblEffLen_glob = 0d0
   lg = 0
   do n = 1,ndambreaksg
      if (db_ndomains(n) <= 1) cycle
      do l = L1dambreaksg(n), L2dambreaksg(n)
         lg = lg+1
         dblEffLen(lg) = dambreakLinksEffectiveLength(l)
      enddo
   enddo
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(dblEffLen,dblEffLen_glob,ndambreak_glob,MPI_DOUBLE_PRECISION,MPI_MAX,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   lg = 0
   do n = 1,ndambreaksg
      if (db_ndomains(n) <= 1) cycle
      maximumDambreakWidths(n) = 0d0
      do l = L1dambreaksg(n), L2dambreaksg(n)
         lg = lg+1
         dambreakLinksEffectiveLength(l) = dblEffLen_glob(lg)
         maximumDambreakWidths(n) = maximumDambreakWidths(n) + dambreakLinksEffectiveLength(l)
      enddo
   enddo
   deallocate(dblEffLen, dblEffLen_glob)
   
   allocate(dsStartBreach_glob(ndambreaksg), LStartBreach_glob(ndambreaksg))
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(dsStartBreach,dsStartBreach_glob,ndambreaksg,MPI_DOUBLE_PRECISION,MPI_MIN,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   
   do n = 1,ndambreaksg
      if (dsStartBreach_glob(n) < dsStartBreach(n)) then
         LStartBreach(n) = 0
      endif
   enddo
   
   LStartBreach = LStartBreach - L1dambreaksg + 1
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(LStartBreach,LStartBreach_glob,ndambreaksg,MPI_INTEGER,MPI_MAX,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   LStartBreach = LStartBreach_glob
   do n = 1, ndambreaksg
      LStartBreach(n) = LStartBreach(n) + L1dambreaksg(n) - 1
      if (LStartBreach(n) > L2dambreaksg(n)) then
         LStartBreach(n) = 0
      endif
   enddo
   deallocate(dsStartBreach_glob, LStartBreach_glob)
end subroutine mpi_dambreak2


!> subroutine to make dambreak information consistent across parallel partitions (call 3)
subroutine mpi_dambreak3()
   use m_flowexternalforcings, only: ndambreaksg, LStartBreach, kdambreak, dambreakCrestLevel, dambreakUpstreamBedLevel, dambreakDownstreamBedLevel
   use m_flowgeom, only: bob, bob0
   use m_partitioninfo, only: jampi, DFM_COMM_DFMWORLD
   use m_timer, only: jatimer, starttimer, stoptimer, IMPIREDUCE
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_MAX
   
   implicit none
   
   double precision, allocatable, dimension(:,:) :: dbLevels           !< characteristic levels at dambreak breach link filled with local values
   double precision, allocatable, dimension(:,:) :: dbLevels_glob      !< characteristic levels at dambreak breach link filled with all values (global)
   integer                                       :: ierr               !< communication error flag
   integer                                       :: l
   integer                                       :: Lf
   integer                                       :: n
   
   allocate(dbLevels(3,ndambreaksg), dbLevels_glob(3,ndambreaksg))
   dbLevels = 0d0
   dbLevels_glob = 0d0
   do n = 1, ndambreaksg
      if (LStartBreach(n) > 0) then
         L = kdambreak(3,LStartBreach(n))
         if (L == 0) cycle
         !
         Lf = iabs(L)
         if (L>0) then
            dbLevels(1,n) = bob0(1,Lf)
            dbLevels(2,n) = bob0(2,Lf)
         else
            dbLevels(1,n) = bob0(2,Lf)
            dbLevels(2,n) = bob0(1,Lf)
         endif
         dbLevels(3,n) = max(bob(1,Lf),bob(2,Lf))
      endif
   enddo
   if (jampi) then
      if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
      call MPI_allreduce(dbLevels,dbLevels_glob,3*ndambreaksg,MPI_DOUBLE_PRECISION,MPI_MAX,DFM_COMM_DFMWORLD, ierr)
#endif
      if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   else
      dbLevels_glob = dbLevels
   endif
   allocate(dambreakUpstreamBedLevel(ndambreaksg), dambreakDownstreamBedLevel(ndambreaksg),dambreakCrestLevel(ndambreaksg))
   do n = 1, ndambreaksg
      dambreakUpstreamBedLevel(n)   = dbLevels_glob(1,n)
      dambreakDownstreamBedLevel(n) = dbLevels_glob(2,n)
      dambreakCrestLevel(n)         = dbLevels_glob(3,n)
   enddo
   deallocate(dbLevels, dbLevels_glob)
end subroutine mpi_dambreak3


!> subroutine to make dambreak information consistent across parallel partitions (call to exchange dambreakMaximum and dambreakAveraging during time step)
subroutine mpi_dambreak_exchange_max_avg()
   use m_flowexternalforcings, only: ndambreaksg, dambreakAveraging, dambreakMaximum, IDB_NQUANT, IDBMAX_NQUANT
   use m_partitioninfo, only: jampi, DFM_COMM_DFMWORLD
   use m_timer, only: jatimer, starttimer, stoptimer, IMPIREDUCE
   use mpi, only: MPI_SUM, MPI_MAX, MPI_DOUBLE_PRECISION
   
   implicit none
   
   double precision, allocatable, dimension(:,:,:):: dambreakAveraging_glob !< array for global information on dambreak quantities (averaging)
   double precision, allocatable, dimension(:,:)  :: dambreakMaximum_glob   !< array for global information on dambreak quantities (maximum)
   integer                                        :: ierr                   !< communication error flag
   integer :: n
   
   if (.not.jampi) return
   
   allocate(dambreakAveraging_glob(3,ndambreaksg,IDB_NQUANT), dambreakMaximum_glob(ndambreaksg,IDBMAX_NQUANT))
   if ( jatimer == 1 ) call starttimer(IMPIREDUCE)
#ifdef HAVE_MPI
   call MPI_allreduce(dambreakAveraging,dambreakAveraging_glob,3*IDB_NQUANT*ndambreaksg,MPI_DOUBLE_PRECISION,MPI_SUM,DFM_COMM_DFMWORLD, ierr)
   call MPI_allreduce(dambreakMaximum,dambreakMaximum_glob,IDB_NQUANT*ndambreaksg,MPI_DOUBLE_PRECISION,MPI_MAX,DFM_COMM_DFMWORLD, ierr)
#endif
   if ( jatimer == 1 ) call stoptimer(IMPIREDUCE)
   dambreakAveraging = dambreakAveraging_glob
   dambreakMaximum   = dambreakMaximum_glob
   deallocate(dambreakAveraging_glob, dambreakMaximum_glob)
end subroutine mpi_dambreak_exchange_max_avg

end module parallel_dambreaks