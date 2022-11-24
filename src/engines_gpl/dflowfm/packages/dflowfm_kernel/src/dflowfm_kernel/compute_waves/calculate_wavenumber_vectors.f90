subroutine calculate_wavenumber_vectors()
   ! Get wave number vector in flowlink orientation, and cartesian in cellcentre
   ! Calculate all these once per timestep, instead of ad hoc
   use m_waves
   use m_sferic
   use m_flowgeom

   integer          :: i, L
   double precision :: ac1, ac2, rk

   do k=1,ndx
      if (hs(k)>epshu) then
         call getwavenr(hs(k), twav(k) ,rk)
         kw(k)  = rk
         kwx(k) = rk*dcos(phiwav(k))
         kwy(k) = rk*dsin(phiwav(k))
      else
         kwx = 0d0
         kwy = 0d0
      endif
   enddo
   !
   kwn = 0d0
   kwt = 0d0
   do i=1,wetlinkcount
      L      = onlyWetLinks(i)
      k1     = ln(1,L); k2=ln(2,L)
      ac1    = acl(L);    ac2 = 1d0-ac1
      kwn(L) = ac1*( csu(L)*kwx(k1) + snu(L)*kwy(k2)) + ac2*( csu(L)*kwx(k1) + snu(L)*kwy(k2))
      kwt(L) = ac1*(-snu(L)*kwx(k1) + csu(L)*kwy(k2)) + ac2*(-snu(L)*kwx(k1) + csu(L)*kwy(k2))
      kwL(L) = hypot(kwn(L), kwt(L))
   enddo

end subroutine calculate_wavenumber_vectors