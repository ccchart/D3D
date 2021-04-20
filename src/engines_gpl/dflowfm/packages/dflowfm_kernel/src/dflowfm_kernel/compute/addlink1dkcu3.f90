 subroutine addlink1Dkcu3(L,japerim)                        ! and add area's and volumes of 1D link kcu3
 use m_flowgeom
 use m_flow
 use m_missing
 use unstruc_channel_flow
 use precision_basics

 implicit none

 integer           :: japerim, L, ja

 integer           :: k1, k2, K, calcConv
 double precision  :: ar1, wid1, cf1, ar2, wid2, cf2, dx1, dx2, widu, perim
 double precision  :: hpr

 if (japerim == 0) then
    calcConv = 0
    k1  = ln(1,L) ; k2 = ln(2,L)

    hpr = max(0d0,s1(k1)-bl(k1))    ! this statement is called most nr of times through waterlevel iteration
    if (hpr > 0) then                             !
       call getprof_1D(L, hpr, ar1, wid1, japerim, calcConv, perim)
       dx1  = dx(L)*acl(L)
       a1(k1) =   a1(k1) + dx1*wid1
       vol1(k1)   = vol1(k1)   + dx1*ar1
    endif

    hpr = max(0d0,s1(k2)-bl(k2))
    if (hpr > 0) then                             !
       call getprof_1D(L, hpr, ar2, wid2, japerim, calcConv, perim)
       dx2  = dx(L)*(1d0-acl(L))
       a1(k2) =   a1(k2) + dx2*wid2
       vol1(k2)   = vol1(k2)   + dx2*ar2
    endif

 else
    if (hu(L) > 0) then
       calcConv = 1
       call getprof_1D(L, hu(L), au(L), widu, japerim, calcConv, perim)
    endif

    calcConv = 0
    if(network%loaded) then
       ! Only in case of a 1d-network, vol1 and vol1_f can be different
       if (kcs(k1) == 1) then ! TODO: consider *also* adding storage area to the 2D side k1, if kcu(L)==5, maybe not for kcu(L)==7
          hpr = s1(k1)-bob0(1,L)
          if (hpr >= 0d0) then
             if (comparereal(hu(L), hpr)== 0) then
                vol1_f(k1) = vol1_f(k1) + dx1*au(L)
             else
                call getprof_1D(L, hpr, ar1, wid1, 1, calcConv, perim)
                vol1_f(k1) = vol1_f(k1) + dx1*ar1
             endif
          endif
       endif
       if (kcs(k2) == 1) then ! TODO: consider *also* adding storage area to the 2D side k2, if kcu(L)==5, maybe not for kcu(L)==7
          hpr = s1(k2)-bob0(2,L)
          if (hpr >= 0d0) then
             ! flow volume
             if (comparereal(hu(L), hpr)== 0) then
                vol1_f(k2) = vol1_f(k2) + dx2*au(L)
             else
                call getprof_1D(L, hpr, ar2, wid2, 1, calcConv, perim)
                vol1_f(k2) = vol1_f(k2) + dx2*ar2
             endif
          endif
       endif
    else

       vol1_f(k1) = vol1(k1)
       vol1_f(k2) = vol1(k2)
    endif

 endif
 end subroutine addlink1Dkcu3
