!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
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

subroutine tauwavefetch(tim)               ! fetchlength and fetchdepth based significant wave height and period
 use m_sediment                             ! based on Hurdle, Stive formulae
 use m_netw                                 ! tauwave based on Swart
 use m_flowgeom                             ! taus = taubmx = taucur + tauwave, as in Delwaq
 use m_flow
 use m_waves, only: fetch, nwf, fetdp, uorb, twav, hwav
 use m_flowtimes
 use m_partitioninfo
 use m_sferic

 implicit none

 double precision :: tim

 integer            :: ierr, k, kb, ki, k1, k2, kkk, L, n, ndoner
 integer, external  :: initialise_fetch_proc_data
 logical, external  :: should_fetch_computation_be_stopped
 logical, parameter :: call_from_tauwavefetch=.true.
 double precision   :: U10, fetchL, fetchd, hsig, tsig, sqrt2, dum
 
 double precision, dimension(:), allocatable :: wxc, wyc

 integer :: ndraw
 COMMON /DRAWTHIS/ ndraw(50)

 if ( .not. allocated (fetch) .or. size (fetch,2) .ne. ndx) then
      nwf = 13

      if (  allocated (fetch) )  deallocate (fetch)
      allocate ( fetch(nwf, ndx) , stat = ierr)
      call aerr('fetch(nwf, ndx)', ierr ,  ndx*nwf)
      if (  allocated (fetdp) )  deallocate (fetdp)
      allocate ( fetdp(nwf, ndx) , stat = ierr)
      call aerr('fetdp(nwf, ndx)', ierr ,  ndx*nwf)

      ndx2dr = 0
      if (jampi == 1 .and. use_fetch_proc == 0 ) then
         allocate ( fett(2, ndx) , stat = ierr)
         call aerr('fett(2, ndx)', ierr ,  ndx*2)

         do k = 1,ndxi
            if (idomain(k) .ne. my_rank) cycle
            if (kcs(k) == 2) ndx2dr = ndx2dr + 1
         end do
         call reduce_int_sum(ndx2dr,ndoner)
         ndx2dr = ndoner
      else
         do k = 1,ndxi
             if (kcs(k) == 2) ndx2dr = ndx2dr + 1
         end do
      endif
             
      if ( use_fetch_proc > 0 ) then
         ierr = initialise_fetch_proc_data()
      endif
     
 endif

 if (tim >= time_fetch) then

    do !  inifinite loop for the fetch proc
        
        if ( use_fetch_proc > 0 ) then
            if ( should_fetch_computation_be_stopped(call_from_tauwavefetch) ) then
                return
            endif
        endif
                
        if ( use_fetch_proc == 1  ) then
           call send_s1_to_fetch_proc()
        endif
     
        time_fetch = max(tim, time_fetch + tifetch )
        if (tifetch == 0d0) time_fetch = 1d30

        if (use_fetch_proc == 0 .or. my_rank == fetch_proc_rank ) then
           call calculate_fetch_values()
        endif
      
        if ( use_fetch_proc == 1  ) then
           call get_fetch_values_from_fetch_proc()
        endif
        
        ! writing for testing purpose
        !call write_dp_data_over_cells("fetch", 1, 1, nwf, ndx, fetch)
        !call mpi_barrier(DFM_COMM_ALLWORLD,ierr)
        !stop
        
        if (use_fetch_proc == 0 .or. my_rank /= fetch_proc_rank) then
           exit
        endif
      enddo
 endif

 sqrt2 = sqrt(2d0)

 do k = 1,ndx2d
    Hwav(k)   = 0d0
    Twav(k)   = 0d0
    Uorb(k)   = 0d0
    rlabda(k) = 0d0

    if ( hs(k) > 0.01d0 ) then

       call getfetch(k,U10,FetchL,FetchD)
       if (FetchL > 0) then

          if (jawave == 1) then

             call hurdlestive (U10, fetchL, fetchD, Hsig, Tsig)

          else if (jawave == 2) then

             call ian_young_pt(U10, fetchL, fetchD, Hsig, Tsig)

          endif

          Hwav(k) = Hsig / sqrt2          ! Hwav === hrms
          Twav(k) = Tsig
          call tauwavehk(Hwav(k), Twav(k), hs(k), Uorb(k), rlabda(k), dum)      ! basically now just a dispersion function with 2DH stokes drift magnitude
       endif
    endif

    if (NDRAW(28) == 35) then
       plotlin(k) = fetchL
    else if (NDRAW(28) == 36) then
       plotlin(k) = fetchD
    else if (NDRAW(28) == 37) then
       plotlin(k) = Hsig
    else if (NDRAW(28) == 38) then
       plotlin(k) = Tsig
    else if (NDRAW(28) == 39) then
       !  plotlin(k) = Taucur
    else if (NDRAW(28) == 40) then
       plotlin(k) = uorb(k)
    endif

 enddo
 
 ! get phiwav
 call realloc(wxc,ndx, keepExisting=.true.)
 call realloc(wyc,ndx, keepExisting=.true.)
 wxc = 0d0; wyc = 0d0
 do L = 1, lnx
    k1 = ln(1,L); k2=ln(2,L)
    wxc(k1) = wxc(k1) + wcL(1,L)*wx(L)
    wxc(k2) = wxc(k2) + wcL(2,L)*wx(L)
    wyc(k1) = wyc(k1) + wcL(1,L)*wy(L)
    wyc(k2) = wyc(k2) + wcL(2,L)*wy(L)  
 enddo   
 phiwav = atan2(wyc,wxc)*180d0/pi
 
 ! Copy values to boundary nodes
 do n = 1, nbndz
    kb = kbndz(1,n)
    ki = kbndz(2,n)
    hwav(kb) = hwav(ki)
    twav(kb) = twav(ki)
    Uorb(kb) = uorb(ki)
    rlabda(kb) = rlabda(ki)
    phiwav(kb) = phiwav(ki)
 enddo
 
  do n = 1, nbndu
    kb = kbndu(1,n)
    ki = kbndu(2,n)
    hwav(kb) = hwav(ki)
    twav(kb) = twav(ki)
    Uorb(kb) = uorb(ki)
    rlabda(kb) = rlabda(ki)
    phiwav(kb) = phiwav(ki)    
 enddo  
 
end subroutine tauwavefetch
 
!> calculates fetch length and depth  
subroutine calculate_fetch_values()
 
 use m_sediment                             
 use m_netw                                 
 use m_flowgeom                             
 use m_flow
 use m_waves,         only: nwf, fetch, fetdp
 use m_partitioninfo
 use unstruc_display, only: jagui
 use geometry_module, only: getdx, getdy, dbdistance, cross, normalout, normalin
 use m_missing,       only: dmiss
 use m_sferic

 implicit none
 
 integer          :: k, kk, kkmin, l, k1, k2, kup, n, ndone, ndoner, ndoneprevcycle, nup, nupf, jaopen, jacros
 double precision :: dir, uwin, vwin, prin, dist, distmin, celsiz, wdep, xn, yn, crp, xkk1, ykk1, xkk2, ykk2
 double precision :: sl, sm, xcr, ycr, fetc, fetd, sumw, cs, sn, dsk2, www
 
 fetch = dmiss ; fetdp = dmiss
 do n  = 1, nwf
    if (jagui > 0) then
        call cls1()
        call setcol(221)
        ! numdots = 0
    endif
	
    dir   = twopi *real(n - 1) / real(nwf - 1)
    uwin  = cos(dir) ; vwin = sin(dir)
    ndone = 0

    do k = 1,ndxi
        if ( kcs(k) /= 2 ) cycle
        kkmin = 0 ; distmin = 1d10; celsiz = 0d0
        if ( jampi == 1  .and. use_fetch_proc == 0 ) then
            if ( idomain(k) /= my_rank) cycle
        endif
        do kk  = 1,netcell(k)%n
            L  = netcell(k)%lin(kk)
            k1 = netcell(k)%nod(kk)
            if ( kk == netcell(k)%n ) then
                k2 = netcell(k)%nod(1)
            else
                k2 = netcell(k)%nod(kk+1)
            endif
            celsiz = max(celsiz, dbdistance(xk(k1), yk(k1), xk(k2), yk(k2), jsferic, jasfer3D, dmiss) )
        enddo
        if (jsferic == 1) celsiz=celsiz*rd2dg/ra ! Herman questioned this line 17.11.2022

        jaopen = 0
        do kk  = 1,nd(k)%lnx
            L  = iabs( nd(k)%ln(kk) )
            if ( ln(1,L) > ndxi ) then
                jaopen = 1
            endif
        enddo

        do kk  = 1,netcell(k)%n
            L  = netcell(k)%lin(kk)
            k1 = netcell(k)%nod(kk)
            if ( kk == netcell(k)%n ) then
                k2 = netcell(k)%nod(1)
            else
                k2 = netcell(k)%nod(kk+1)
            endif

            wdep = s1(k) - min(zk(k1),zk(k2))
            if ( lnn(L) == 1 .or.  wdep < 0.5d0 .or. kn(3,L) == 0 .or. jaopen == 1 ) then    ! link shallow or closed => start fetch here
                call normalout(xk(k1), yk(k1), xk(k2), yk(k2), xn, yn, jsferic, jasfer3D, dmiss, dxymis)
                prin = uwin*xn + vwin*yn
                if ( prin < 0d0 ) then                   ! if upwind
                    crp  = xn ; xn  = -yn ; yn = crp
                    crp  = 0d0
                    xkk1 = xk(k1) - 2*celsiz*xn
                    ykk1 = yk(k1) - 2*celsiz*yn
                    xkk2 = xk(k2) + 2*celsiz*xn
                    ykk2 = yk(k2) + 2*celsiz*yn
                    call cross(xkk1,ykk1,xkk2,ykk2,xzw(k),yzw(k),xzw(k)-1d4*uwin,yzw(k)-1d4*vwin, &
                                jacros,sl,sm,xcr,ycr,crp,jsferic, dmiss)
                    if ( jacros == 1 ) then
                        dist = dbdistance(xz(k), yz(k), xcr, ycr, jsferic, jasfer3D, dmiss)
                        if ( dist < distmin ) then
                           distmin = dist ; kkmin = kk        ! closest crossed upwind edge
                        endif
                    endif
                endif
            endif
        enddo
		
        if ( kkmin > 0 ) then
            if ( jaopen == 1 ) then
                fetch(n,k) = 1d5
            else
                fetch(n,k) = min(distmin, celsiz)
            endif
            fetdp(n,k) = max( s1(k) - bl(k), .1d0)
            if ( jagui > 0 ) then
                   !CALL rCIRc(Xz(k),Yz(k) ) !, fetch(n,k))
                   !call adddot(Xz(k),Yz(k),1d0)
            endif
            ndone = ndone + 1
        endif

    enddo

    if ( jampi == 1 .and. use_fetch_proc == 0 ) then
        call reducefett(n)
        call reduce_int_sum(ndone,ndoner)
        ndone = ndoner
    endif

    if ( jagui > 0 ) call setcol(31)
	
    do while ( ndone < ndx2dr )

        ndoneprevcycle = ndone
        ndone          = 0

        do k = 1, ndxi
            if ( kcs(k) /= 2 ) cycle
            if ( jampi == 1 .and. use_fetch_proc == 0 ) then
                if ( idomain(k) /= my_rank ) cycle
            endif

            if ( fetch(n,k) == dmiss ) then
                kup = 0 ; fetc = 0; fetd = 0; sumw = 0; nup = 0; nupf = 0
                do kk  = 1,nd(k)%lnx
                    L  = iabs( nd(k)%ln(kk) )
                    k2 = ln(1,L) ; if (k2 == k) k2 = ln(2,L)
                    if ( kcs(k2) == 2 ) then  ! internal
                        !prin = uwin*getdx(xz(k2),yz(k2),xz(k),yz(k), jsferic) + vwin*getdy( xz(k2),yz(k2),xz(k),yz(k), jsferic)
                        !dsk2 = dbdistance(xz(k2),yz(k2),xz(k),yz(k), jsferic, jasfer3D, dmiss)
                        !cs   = min(max(prin/dsk2,-1d0),1d0)

                        cs   = uwin*csu(L) + vwin*snu(L)
                        if ( L /= nd(k)%ln(kk) ) cs = -1d0*cs
                        dsk2 = dx(L)
                        prin = dsk2*cs

                        if ( cs > 0 ) then ! internal upwind points
                            nup = nup + 1
                            if ( fetch(n,k2) .ne. dmiss ) then ! do not look at open boundaries
                               nupf = nupf + 1
                               sn   = sqrt( 1d0 - cs*cs)
                               ! www  = (1d0-sn)/dsk2               ! first attempt
                               www  = (cs   + 0.05d0*sn)*wu(L)/dsk2 ! some diffusion
                               fetc = fetc  + www*(fetch(n,k2) + prin)
                               fetd = fetd  + www*(fetch(n,k2) + prin)*max(.1d0, 0.8d0*fetdp(n,k2) + 0.2d0*(s1(k)-bl(k)) )
                               sumw = sumw  + www
                            endif
                        endif
                    endif
                enddo
                if ( nup == nupf .and. sumw > 0d0 ) then
                    fetch(n,k) = fetc/sumw
                    fetdp(n,k) = fetd/ ( sumw*fetch(n,k) )
                    ndone      = ndone + 1
                    if ( jagui > 0 ) then
                        !CALL rCIRc(Xz(k),Yz(k) )
                        !call adddot(Xz(k),Yz(k),2d0)
                        call KCIR(Xz(k),Yz(k),1d0)
                    endif
                endif
            else
                ndone = ndone + 1
            endif
        enddo ! k

        if ( jampi == 1 .and. use_fetch_proc == 0 ) then
            call reducefett(n)
            call reduce_int_sum(ndone,ndoner)
            ndone = ndoner
        endif

        if ( ndone == ndoneprevcycle ) then
            call qnerror('connectivity issue in fetch', ' ', ' ')
            return
        endif

    enddo

enddo 

end subroutine calculate_fetch_values


