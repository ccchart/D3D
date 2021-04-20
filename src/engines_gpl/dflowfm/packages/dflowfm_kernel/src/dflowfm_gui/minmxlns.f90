 subroutine minmxlns()

 use m_flowgeom
 use m_flow
 use m_missing

 implicit none
 double precision :: zlin
 double precision :: zn
 double precision :: rmin, rmax
 integer          :: i, l, n, k1, k2
 double precision :: VMAX,VMIN,DV,VAL(256)
 integer :: NCOLS(256),NIS,NIE,nv,JAAUTO
 common /depmax2/ vmax,vmin,dv,val,ncols,nv,nis,nie,jaauto
logical inview

 if (jaauto > 0) then
    rmin =  1d30; lnmin = 0
    rmax = -1d30; lnmax = 0
    do L = 1,lnx
       k1 = ln(1,L)
       k2 = ln(2,L)
       if (inview( xz(k1), yz(k1) ) .or. inview( xz(k2), yz(k2) ) ) then
           zn   = zlin(L)
           if ( zn.eq.DMISS ) cycle
           if (zn < rmin) then
              rmin = zn ; lnmin = L
           endif
           if (zn > rmax) then
              rmax = zn ; lnmax = L
           endif
       endif
    enddo
    vmax = rmax
    vmin = rmin
 endif

 dv   = vmax - vmin
 do i = 1,nv
    val(i) = vmin + (i-1)*dv/(nv-1)
 enddo

 return
 end subroutine minmxlns
