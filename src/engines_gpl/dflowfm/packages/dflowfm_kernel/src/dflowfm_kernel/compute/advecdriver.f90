 subroutine advecdriver()
 use m_flowtimes
 use m_flow
 use m_flowgeom
 implicit none

 double precision :: dta, das, ds
 integer          :: L, k1, k2, k

 if (itstep == 3) then

   if (.not. allocated(adve0) ) then
      allocate(adve0(lnkx))
   endif

   dta = 0.7D0*dts/cflmx
   das = dta/dts
   do k = 1,2

      adve = 0d0

      call advec()
      if (k == 1) then
         adve0 = adve
      endif
      do L = 1,lnx
         k1 = ln(1,L) ; k2 = ln(2,L)
         ds = ag*dxi(L)*(s0(k2) - s0(k1))
         u1(L) = ( u1(L)*(1d0 - das) + u0(L)*das - dta*(adve(L) + ds) ) / (1d0 + dta*advi(L))
      enddo
      call setucxucyucxuucyunew()

   enddo
   ! adve = teta0*adve + (1d0-teta0)*adve0
   ! u1 = u0

 else

   call advec()                                       ! advection term, must be called after set-umod and cell velocity updates

 endif

 call setextforcechkadvec()                           ! set external forcings and check explicit part adve
 end subroutine advecdriver
