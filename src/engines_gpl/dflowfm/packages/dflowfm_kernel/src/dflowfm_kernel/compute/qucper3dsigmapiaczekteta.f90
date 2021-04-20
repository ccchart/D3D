 subroutine QucPer3Dsigmapiaczekteta(LL,Lb,Lt,cs,sn,ae,ai)  ! Piaczekteta in 3D

 use m_flow                                                 ! advect the cell center velocities (dimension: m4/s2)
 use m_flowgeom
 use m_flowtimes, only : dts                                !
 use m_sferic
 implicit none

 integer,          intent(in) :: LL,Lb,Lt                   ! working for basis link LL
 double precision, intent(in) :: cs, sn
 double precision, intent(out):: ae(Lt-Lb+1)                ! explicit part
 double precision, intent(out):: ai(Lt-Lb+1)                ! implicit part


 ! locals
 integer          :: La, LLL, LLLL, Lb2, Lt2, Lk            ! for links LL,
 integer          :: k12, n12, k1, k2                       ! relevant node, 1 or 2, L/R
 double precision :: ucin, cfl, tet, volu, ac, acq          ! velocity surplus

 double precision ::  ucinx, uciny
 integer          :: nn12

 double precision, external:: lin2nodx, lin2nody, nod2linx, nod2liny

 ae = 0d0 ; ai = 0d0

 do n12 = 1,2
    if (n12 ==1) then
        ac = acL(LL)
    else
        ac = 1d0-acL(LL)
    endif
    k12  = ln(n12,LL)
    do La   = 1, nd(k12)%lnx                                ! loop over all attached links
       LLL  = nd(k12)%ln(La)
       nn12 = 1; if ( LLL.gt.0 ) nn12 = 2
       LLLL = iabs(LLL)

       Lb2  = Lbot(LLLL) ; Lt2 = Ltop(LLLL)
       do Lk = LB2, LT2

          if ( qa(Lk) .ne. 0) then                          ! include own link
             k1   = ln(1,Lb+Lk-Lb2)  ;  k2 = ln(2,Lb+Lk-Lb2)
             volu = acL(LL)*vol1(k1) + (1d0-acl(LL))*vol1(k2)
             if (volu > 0d0) then
                cfl  = abs(qa(Lk))*dts/volu
                if (nd(k12)%lnx ==3) cfl=1.4d0*cfl
                if (cfl > 0d0) then
                   tet  = max(0d0, 1d0 - 1d0/cfl  )
                   if (jasfer3D == 0) then
                      ucin = ucxu(Lk)*cs + ucyu(Lk)*sn  - (1d0-tet)*u1(Lb + Lk - Lb2)
                   else
                      ucinx = lin2nodx(LLLL,nn12,ucxu(Lk),ucyu(Lk))
                      uciny = lin2nody(LLLL,nn12,ucxu(Lk),ucyu(Lk))
                      ucin  = nod2linx(LL,n12,ucinx,uciny)*cs + nod2liny(LL,n12,ucinx,uciny)*sn  - (1d0-tet)*u1(Lb + Lk - Lb2)
                   endif

                   acq  = ac*qa(Lk)/volu
                   if (LLL > 0) then                        ! incoming link
                      ae(Lk-Lb2+1) = ae(Lk-Lb2+1) - acq*ucin
                      ai(Lk-Lb2+1) = ai(Lk-Lb2+1) + acq*tet
                   else
                      ae(Lk-Lb2+1) = ae(Lk-Lb2+1) + acq*ucin
                      ai(Lk-Lb2+1) = ai(Lk-Lb2+1) - acq*tet
                   endif
                endif
             endif

          endif

       enddo

    enddo

 enddo

 end subroutine QucPer3Dsigmapiaczekteta
