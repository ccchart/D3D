 !> adjust bobs and iadvec for dams and structs
 subroutine adjust_bobs_for_dams_and_structs()
    use m_flowgeom
    use m_flow
    use m_netw
    use m_fixedweirs
    use unstruc_channel_flow
    use m_1d_structures
    use m_compound

    implicit none

    double precision :: zcdamn, minzcdamn, blmx
    type(t_structure), pointer :: pstru
    type(t_compound),  pointer :: pcompound

    integer :: L0
    integer          :: ng, k1, k2, L, n, istru, icompound, i

    do ng = 1,ncdamsg                                   ! loop over cdam signals, sethu
       zcdamn = zcdam(ng)
       do n   = L1cdamsg(ng), L2cdamsg(ng)
           k1       = kcdam(1,n)
           k2       = kcdam(2,n)
           L        = kcdam(3,n)
           blmx     = max(bl(k1), bl(k2))
           bob(1,L) = max(zcdamn,blmx)
           bob(2,L) = max(zcdamn,blmx)
        enddo
    enddo

    do ng = 1,ncgensg                                   ! loop over general structures signals, sethu
       zcdamn = zcgen(3*(ng-1)+1) ! TODO: actually, the crest/sill_width should be included here: not all flow links may be open
       do n   = L1cgensg(ng), L2cgensg(ng)
           k1       = kcgen(1,n)
           k2       = kcgen(2,n)
           L        = kcgen(3,n)
           blmx     = max(bl(k1), bl(k2))
           bob(1,L) = max(zcdamn,blmx)
           bob(2,L) = max(zcdamn,blmx)
           call switchiadvnearlink(L)
       enddo
    enddo

    do istru = 1, network%sts%count
        pstru => network%sts%struct(istru)
        zcdamn = get_crest_level(pstru)
        if (zcdamn == huge(1d0)) then
           ! Do not shut off structures that have no relevant crest (e.g. pumps)
           cycle
        end if

        do L0 = 1, pstru%numlinks
           L  = iabs(pstru%linknumbers(L0))
           k1 = ln(1,L)
           k2 = ln(2,L)
           bob(1,L) = max(zcdamn, bob0(1, L))
           bob(2,L) = max(zcdamn, bob0(2, L))
           iadv(L) = 22
           call switchiadvnearlink(L)
        enddo

    enddo

    ! correct BOBS for compound structures
    do icompound = 1, network%cmps%Count
       pcompound => network%cmps%compound(icompound)
       minzcdamn = huge(1d0)
       do i = 1, pcompound%numstructs
          istru = pcompound%structure_indices(i)
          pstru => network%sts%struct(istru)
          zcdamn = get_crest_level(pstru)
          if (zcdamn == huge(1d0)) then
             ! Obviously this is a pump. So do not adust the bob
             minzcdamn = huge(1d0)
             exit
          endif

          minzcdamn = min(minzcdamn, zcdamn)
       enddo
       if (minzcdamn < huge(1d0)) then
          do L0 = 1, pcompound%numlinks
             L  = iabs(pcompound%linknumbers(L0))
             k1 = ln(1,L)
             k2 = ln(2,L)
             bob(1,L) = max(minzcdamn, bob0(1, L))
             bob(2,L) = max(minzcdamn, bob0(2, L))
          enddo
       endif
    enddo


   !Adjust bobs for dambreak
   if (ndambreak > 0) then ! needed, because ndambreaksg may be > 0, but ndambreak==0, and then arrays are not available.
   do n = 1, ndambreaksg
      istru = dambreaks(n)
      if (istru.ne.0) then
         ! Update the bottom levels
         call adjust_bobs_on_dambreak_breach(network%sts%struct(istru)%dambreak%width, network%sts%struct(istru)%dambreak%crl,  LStartBreach(n), L1dambreaksg(n), L2dambreaksg(n), network%sts%struct(istru)%id)
      endif
   enddo
   end if

   return
   end subroutine adjust_bobs_for_dams_and_structs
