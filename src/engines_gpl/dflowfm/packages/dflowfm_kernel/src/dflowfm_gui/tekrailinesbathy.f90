 subroutine tekrailinesBATHY(ncol,jaall,ITYP)
 use m_flowgeom
 USE M_FLOW
 use m_flowtimes
 use m_sferic
 use unstruc_display
 use m_netw, only : xk,yk,zk
 use m_sediment

 implicit none
 integer          :: nx, ncol, jaall, ITYP
 integer          :: r, L, k1,k2, kk,k,n
 double precision :: zz1, zz2, xx1, xx2, yy1, yy2
 integer          :: ja, jg

 call setcol(ncol)
 do L = 1,lnx
    if (mod(L,200) == 0) then
       call halt2(ja)
       if (ja == 1) exit
    endif

    if (ityp == 1) then                                ! bottom layer
       if (ibedlevtyp == 1 .or. ibedlevtyp == 6) then  ! tegelen
          k1  = ln  (1,L)
          k2  = ln  (2,L)
          xx1 = xz(k1)
          xx2 = xz(k2)
          zz1 = bl(k1)
          zz2 = bl(k2)
          if (yfac > 0) then
             yy1 = yz(k1)
             yy2 = yz(k2)
          endif
       else
          k1  = lncn(1,L)
          k2  = lncn(2,L)
          xx1 = xk(k1)
          xx2 = xk(k2)
          zz1 = zk(k1)
          zz2 = zk(k2)
          if (yfac > 0) then
             yy1 = yk(k1)
             yy2 = yk(k2)
          endif
       endif
    else                                             ! non erodable layer
       jg = ityp - 1
       if (jaceneqtr == 1 ) then                     ! combined data not really drawn precise if no tegel
          k1  = ln (1,L)
          k2  = ln  (2,L)
          xx1 = xz(k1)
          xx2 = xz(k2)
          zz1 = bl(k1) - sum(grainlay(1:jg,k1) )
          zz2 = bl(k2) - sum(grainlay(1:jg,k2) )
          if (yfac > 0) then
             yy1 = yz(k1)
             yy2 = yz(k2)
          endif
       else
          k1  = lncn(1,L)
          k2  = lncn(2,L)
          xx1 = xk(k1)
          xx2 = xk(k2)
          zz1 = zk(k1) - sum(grainlay(1:jg,k1) )
          zz2 = zk(k2) - sum(grainlay(1:jg,k2) )
          if (yfac > 0) then
             yy1 = yk(k1)
             yy2 = yk(k2)
          endif
       endif
    endif

    if (yfac > 0) then
       zz1 = zz1 + (yy1 - ymn)*yfac
       zz2 = zz2 + (yy1 - ymn)*yfac
    endif

    if (jsferic == 1) then ! jglobe
       if (abs( xz(ln(1,L)) - xz(ln(2,L) ) ) > 10d0) cycle
    endif

    if (abs(zz1) < 1d-6) zz1 = 0d0  ! heh heh, eindelijk
    if (abs(zz2) < 1d-6) zz2 = 0d0

    call movabs(xx1, zz1 )
    call  lnabs(xx2, zz2 )
 enddo

 do L = 1,lnx1D
    k1  = ln(1,L)
    k2  = ln(2,L)
    xx1 = xk(k1)
    xx2 = xk(k2)
    zz1 = bl(k1)
    zz2 = bl(k2)
    if (yfac > 0) then
        yy1 = yk(k1)
        yy2 = yk(k2)
        zz1 = zz1 + (yy1 - ymn)*yfac
        zz2 = zz2 + (yy1 - ymn)*yfac
    endif
    call movabs(xx1, zz1 )
    call  lnabs(xx2, zz2 )
 enddo



 if (jaceneqtr == 2 .and. ityp ==4) then
     do kk  = 1,mxban
        n   = nban(1,kk)
        k   = nban(2,kk)
        xx1 = xz(k)
        yy1 = yz(k)
        zz1 = zmn + grainlay(1,kk) ! bl(k) - grainlay(1,kk)
        xx2 = xk(n)
        yy2 = yk(n)
        zz2 = zmn + grainlay(1,kk) ! zk(n) - grainlay(1,kk)

        call movabs(xx1, zz1 )
        call  lnabs(xx2, zz2 )
     enddo
 endif

 end subroutine tekrailinesBATHY
