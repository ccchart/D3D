subroutine sethu(jazws0)                            ! Set upwind waterdepth hu
 use m_flowgeom                                      ! Todo: higher order + limiter, see transport
 use m_flow
 use m_flowtimes
 use m_sediment
 use m_fixedweirs
 use m_sobekdfm
 use m_sferic
 use m_missing
 use m_netw, only: xk, yk
 use unstruc_model, only:md_restartfile
 use geometry_module, only: dbdistance

 implicit none

 ! locals
 integer           :: L, k1, k2, ku, kd, isg, LL, k, Ld, iup, nq, kk, ifrctyp, jazws0
 integer           :: n, kb, kb0, kt, itpbn, ng, jawet, Lb, nfw
 double precision  :: zb, hh, dtgh
 double precision  :: sup, bup, sk1, sk2, hs1, hs2, epsh, qdak
 double precision  :: hsav, hul, utp, hup
 double precision  :: huk1, huk2, hsku, sigm
 double precision  :: onet = 1d0/3d0
 double precision  :: twot = 2d0/3d0
 double precision  :: tgi  = 1d0/(2d0*9.81d0)
 double precision  :: ds, ds1, ds2, ds0, xx, hdl, bupmin, he, zcdamn, hcrest, uin, blmx, fdx
 double precision  :: h0, dzb, cz, sixth = 1d0/6d0, frcn, z00, sqcf, uuL, vhei, eup, ucxku, ucyku, vben, uLL, agwdxi

 double precision  :: Qweir_super, Qweir_sub, wu_orig
 double precision  :: weirrelax = 0.75d0
 double precision  :: Qrat, re, Edown, ucxkd, ucykd
 double precision  :: dp, d2, aa, hucrest, hunoweir, qweirsimple, ufac, efac

 double precision  :: avolk, hkruin, wsbov, wsben, d1, ewben, eweir, qvolk, qunit, hov, vov, vbov, hvolk, dte0, dtefri, qov, tol
 double precision  :: sl1, sl2, sl3, sku, hskub, hub, sigmd
 character (len=4) :: toest

 integer           :: k3, k4, itel, kuu, ku2, kku, ip , Lnu, kbd, ktd, kbd0, LLbup

 double precision, external :: dslim,  nod2linx, nod2liny

 ! SPvdP: s0 at the old time level already satisfies the boundary conditions at the old time level (see s1nod)
 !  Nevertheless, s0 at the boundary (at the old time-level) will now be filled with boundary conditions at the new time level
 if(jazws0==0 .or. len_trim(md_restartfile)==0) then
    ! if(jazws0==1 .and. len_trim(md_restartfile)>0) then s0 and s1 are read from restart file, in this case, no need to call the following subroutine
    call sets01zbnd(0, 0)                               ! set s0 on z-boundaries
 endif

 if (uniformhu > 0d0) then
    hu = uniformhu ; return
 endif

 !
 ! SPvdP: water-levels at the velocity boundaries do already satisfy the Neumann condition (see s1nod)
 !
 !do n  = 1, nbndu                                    ! velocity boundaries
 !   kb     = kbndu(1,n)
 !   k2     = kbndu(2,n)
 !   s0(kb) = s0(k2)
 !enddo

 ! adjust bobs for controllable dams
 call adjust_bobs_for_dams_and_structs()

 avolk = twot*sqrt(twot*ag)
 nfw   = 0

 do L = 1,lnx

!   for cut-cells
    if (wu(L).eq.0d0 ) then
       hu(L) = 0d0
       cycle
    end if

    k1 = ln(1,L) ; k2 = ln(2,L)

    if (jazws0 == 0) then
       uuL = u1(L)
    else
       uuL = u0(L)
    endif

    if ( uuL > 0) then           ! was q1 (halfway the timestep), jazws0 assigns if start of loop or loop itself is considered
       iup = 1  ; ku = k1 ; kd = k2 ; isg =  1
    else if (uuL < 0) then
       iup = 2  ; ku = k2 ; kd = k1 ; isg = -1
    else if ( s0(k1) > s0(k2) ) then
       iup = 1  ; ku = k1 ; kd = k2 ; isg =  1
    else
       iup = 2  ; ku = k2 ; kd = k1 ; isg = -1
    endif

    sup = s0(ku)

    if (limtyphu > 0 ) then
        if (limtyphu == 21) then      ! central
           sup = 0.5d0*( s0(k1) + s0(k2) )
        else if (limtyphu == 22) then ! perot alfa
           sup = acl(L)*s0(k1) + (1d0-acl(L))*s0(k2)
        else if (limtyphu == 23) then ! regular linear interpolation
           sup = acl(L)*s0(k2) + (1d0-acl(L))*s0(k1)
        else                          ! usual limiters except 6
           if (uuL > 0) then
              ip = 0
           else
              ip = 3
           endif
           ds2 = s0(kd) - s0(ku)

           kku = klnup(1+ip,L)
           kuu = abs(kku) ; sku = dmiss
           if (kku < 0) then
              sku = s0(kuu)
           else if (kuu > 0) then
              ku2 = iabs(klnup(2+ip,L))
              if ( ku2 > 0) then
                 sl1 = slnup(1+ip,L) ; sl2  = slnup(2+ip,L)
                 sku = s0(kuu)*sl1 + s0(ku2)*sl2
              endif
           endif
           if (sku .ne. dmiss) then
              sl3 = slnup(3+ip,L)
              ds1 = (s0(ku)  - sku)*sl3
              sup = sup + dslim(ds1, ds2, limtyphu)
           endif
        endif
    endif

    ! TODO: while documenting 1D2D code, we discovered the following undesirable bup:
    ! it should by default be the max(bob1/2), if conveyance2D < 1. Not yet changed.
    bup = bob(iup,L)

    if (L <= lnx1D) then      ! 1D
       if (kcu(L) == 4 .and. jaconveyance2D >= 1) then
           bup = min( bob(1,L), bob(2,L) )
       else if (kcu(L) == 5 .or. kcu(L) == 7) then
           bup = max( bob(1,L), bob(2,L) )
       else
           bup = max( bob(1,L), bob(2,L) )
       endif
       if (jagrounlay > 0) then
          bup = bup + grounlay(L)
       endif

    else if (kmx == 0 .and. jaconveyance2D >= 1) then
       bup = min( bob(1,L), bob(2,L) )
    endif

    huL = sup-bup

    if (huL  > epshu) then

       if (ncdamsg > 0 .or. ifixedweirscheme > 0) then                           ! sethu

          if (iadv(L) == 21 .or. iadv(L) >= 23 .and. iadv(L) <= 25) then         ! weir velocity point


             if (iadv(L) >= 23 .and. iadv(L) <= 25) then                         ! undisturbed velocity as if no weir present, WAQUA like
                hunoweir  = sup - blu(L) ! bob(1,L)                              ! 23 = Rajaratnam, 24 = Tabellenboek, 25 = Villemonte
                hunoweir  = max(hunoweir, huL)
                ucxku     = ucx(ku) ; ucyku = ucy(ku)
             else
                call getucxucynoweirs(ku, ucxku, ucyku, ifixedweirscheme )
             endif

             if (iadv(L) >= 23 .and. iadv(L) <= 25) then                         ! 23 = Rajaratnam, 24 = Tabellenboek, 25 = Villemonte
                 ! uin = ucxku*csfxw(nfw) + ucyku*snfxw(nfw)
                 uin = abs(u1(L))
             else
                 uin = ucxku*csu(L) + ucyku*snu(L)                               ! semi subgrid
             endif
             vhei = 0.5d0*uin*uin / ag
             eup  = hul + vhei

             if (iadv(L) == 21 .or. iadv(L) == 23) then

                hcrest= s0(kd) - bup

                if ( hcrest < hul ) then

                    huL = hcrest

                    if (hul < twot*eup) then ! supercritical

                       hul = twot*eup

                       hup = hcrest - hul

                       if (hup < 0) then
                          adve(L) = adve(L) - isg*hup*ag*dxi(L)
                       endif

                    endif


                endif ! hcrest< hul

             endif

             if (iadv(L) == -21) then !+-21

                if (jasfer3D == 1) then
                   uin = nod2linx(L,iup,ucxku,ucyku)*csu(L) + nod2liny(L,iup,ucxku,ucyku)*snu(L)
                endif
                fdx     = 0.5d0*dxi(L)*isg
                advi(L) = advi(L) + fdx*u0(L)
                adve(L) = adve(L) - fdx*uin*uin

             else if (iadv(L) == 23) then       ! simple Rajaratnam

                 ufac    = hunoweir / huL  ! compensates for undisturbed field velocity
                 efac    = 1d0 - (1d0/ufac**2)
                 advi(L) = advi(L) + 0.5d0*dxi(L)*abs(u1(L))*ufac*ufac*efac
                 huL     = hunoweir

             else if (iadv(L) == 24 .or. iadv(L) == 25) then  !  Tabellenboek or Villemonte from WAQUA

                 nfw    =  nfxwL(L)

                 wsbov  =  sup
                 wsben  =  s0(kd)
                 hkruin = -bup

                 ! d1     =  bup - blu(L)   !! old implementation

                 ! determine sill height downstream of weir
                 !
                 if (uin .ge. 0.0 ) then
                     d1 = shrxw(nfw)
                 else
                     d1 = shlxw(nfw)
                 endif

                 ! vbov   =  abs(u1(L))
                 ! vbov   =  sqrt(ucxku*ucxku + ucyku*ucyku)
                 vbov   =  abs(uin)
                 vhei   =  0.5d0*vbov*vbov / ag
                 eweir  =  max (0.000001d0, wsbov + hkruin) + vhei
                 qvolk  =  avolk*eweir**1.5d0
                 qunit  =  vbov*hunoweir

                 ! Compute energy height downstream (EWBEN)
                 vben   = qunit / max (0.000001d0,wsben - bl(kd))
                 vhei   =  0.5d0*vben*vben / ag
                 ewben  =  max (0.000001d0, wsben + hkruin) + vhei
                 ! limit downstream energy height EWBEN by upstream enegy height EWEIR
                 ewben = min(ewben, eweir)

                 ! Qunit  = abs(q1(L)) / wu(L)

                 hov    =  wsbov + hkruin
                 vov    =  qunit/hov
                 if (vov < 0.5d0 ) then
                    itel  = 0
                    hvolk = twot*eweir
                    tol   = 0.001d0 *max(0.0001d0, qunit)
                    qov   = 0d0
                    do while (itel < 100 .and. (abs(qunit - qov)) > tol )
                       itel = itel + 1
                       vov  = qunit / hov
                       hov  = max(hvolk, eweir - (vov**2)/(2d0*ag) )
                       qov  = vov*hov
                    enddo
                 endif
                 dte0   = weirdte(nfw)
                 dtefri = 0.0d0
                 call enloss(ag, d1, eweir, hkruin, hov,                   &
                             qunit, qvolk, toest, vov,                     &
                             ewben, wsbov, wsben, weirdte(nfw),              &
                             dtefri,iadv(L), crestlxw(nfw),                &
                             taludlxw(nfw), taludrxw(nfw), vegxw(nfw) )
                 weirdte(nfw) = (1d0 - waquaweirthetaw)*weirdte(nfw) + waquaweirthetaw*dte0
                 !
                 ! attention total waterdepth instead of water above crest
                 if ( toest == 'volk' ) then
                     vbov = qvolk/max(hunoweir, 1d-6 )
                 endif
                 if (vbov > 1d-4) then
                     agwdxi  = ag*weirdte(nfw)*dxi(L)
                     if (kmx == 0) then
                        advi(L) = advi(L) + agwdxi/vbov        ! 1/s
                     else
                        do LL = Lbot(L), Ltop(L)
                           uLL      = max(1d-4, abs(u1(LL)))
                           advi(LL) = advi(LL) + agwdxi/uLL
                        enddo
                    endif
                 endif
                 huL = hunoweir

             else


             endif

          endif ! kadepunt

       endif

       hu(L) = huL

    else
       hu(L) = 0d0
    endif

    if (kmx > 0) then
       Lb       = Lbot(L)
       if(hu(L) > 0d0) then
          kt      = ktop(ku)
          kb      = min ( ln0( iup,Lb ) , kt )  ! dickv, was ln

          kb0     = kb - 1                   ! kbot(ku) - 1
          Ltop(L) = Lb + kt - kb

          if (Ltop(L) > Lb + kmxL(L) - 1) then
             call qnerror('Ltop too large',' ',' ')
          endif

          hsku  = zws(kt) - zws(kb0)
          au(L) = 0d0
          hu(Lb-1) = 0d0

          if (layertype == 2 .and. keepzlayeringatbed == 3) then  ! split in a central sigma oriented part

             if ( Lb == Ltop(L) ) then                ! one layer
                LL     = Lb
                hu(LL) = hu(L)
                au(LL) = wu(L)*(hu(LL) - hu(LL-1))    ! this is only for now here, later move to addlink etc
                au(L)  = au(L) + au(LL)               ! add to integrated 2Dh layer

             else                                     ! two or more

                ktd  = ktop(kd)
                kbd  = min ( ln0(3-iup,Lb ) , ktd )
                kbd0 = kbd - 1


                hub  = 0d0
                do LL  = Lb+1, Ltop(L)                ! search upwind cell for first layer above local bob
                   hub = zws(kb+LL-Lb) - bup
                   if (hub > 0) then
                      LLbup = LL
                      exit
                   endif
                enddo

                do LL  = Lb, LLbup                    ! central in lower part
                   sigm   = ( zws(kb+LL-Lb)  - zws(kb0)  ) / ( zws(kb+LLbup-Lb)  - zws(kb0)  )
                   sigmd  = ( zws(kbd+LL-Lb) - zws(kbd0) ) / ( zws(kbd+LLbup-Lb) - zws(kbd0) )
                   sigm   = 0.5d0*(sigm + sigmd)
                   hu(LL) = sigm*hub
                   au(LL) = wu(L)*(hu(LL) - hu(LL-1)) ! this is only for now here, later move to addlink etc
                   au(L)  = au(L) + au(LL)            ! add to integrated 2Dh layer
                enddo

                hub = hu(L) - hub
                do LL  = LLbup+1, Ltop(L)             ! upwind in upper part
                   sigm   = ( zws(kb+LL-Lb) - zws(kb+LLbup-Lb) ) / ( zws(kt) - zws(kb+LLbup-Lb) )
                   hu(LL) = hu(LLbup) + sigm*hub
                   au(LL) = wu(L)*(hu(LL) - hu(LL-1)) ! this is only for now here, later move to addlink etc
                   au(L)  = au(L) + au(LL)            ! add to integrated 2Dh layer
                enddo

            endif

          else                                     ! default: upwind sigma oriented distribution of hu(L)

             do LL = Lb, Ltop(L)
                sigm   = (zws(kb+LL-Lb)-zws(kb0)) / hsku
                hu(LL) = sigm*hu(L)
                au(LL) = wu(L)*(hu(LL)-hu(LL-1))   ! this is only for now here, later move to addlink etc
                au(L)  = au(L) + au(LL)            ! add to integrated 2Dh layer
             enddo

          endif

       else
          Ltop(L) = 1 ! lb - 1 ! 1 ! flag dry
       endif

    endif

 enddo

 do L = 1,lnx
    k1 = ln(1,L) ; k2 = ln(2,L)
    hsav  = max(epshs, acl(L)*hs(k1) + (1d0-acl(L))*hs(k2) )
    huvli(L) = 1d0 / hsav
 enddo

 if (lincontin == 1) then
    do L = 1,lnx
       hu(L) = -0.5d0*( bob(1,L) + bob(2,L) )
    enddo
 endif


 if (nbnd1d2d > 0) then       ! 1d2d boundary check for closed boundaries
    call sethu_1d2d()
 endif

 if (javeg > 0) then
    call setveg()
 endif

end subroutine sethu
