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

 subroutine addbarocL(LL,Lb,Lt)
 use m_flowgeom
 use m_flow
 use m_flowtimes

 implicit none
 integer, intent(in) :: LL,Lb,Lt

 integer             :: L, k1, k2, k1t, k2t, k, kt, kz ! ,%ktz
 double precision    :: gradpu(kmxx), rhovol(kmxx), gr3, barocL, ft, dum
 double precision    :: rv1, rv2, gr1, gr2, rvk, grk, saw0, saw1, tmw0, tmw1, fzu, fzd, dzz, rv0, rhow0, rhow1, pdb, p0d

 do L = Lb,Lt
    k1 = ln(1,L) ; k1t = k1
    k2 = ln(2,L) ; k2t = k2
    if (L == Lt) then
       k1t = ktop(ln(1,LL)) ; k2t = ktop(ln(2,LL))
    endif

    rhovol(L-Lb+1) = 0.5d0*( zws(k1t) - zws(k1-1) + zws(k2t) - zws(k2-1) )*0.5d0*(rho(k1) + rho(k2)) ! write in rvdn

    ! rhovol(L-Lb+1) = 0.5d0*( (zws(k1t) - zws(k1-1))*rho(k1) + (zws(k2t) - zws(k2-1))*rho(k2) )
    if (jarhoxu > 0) then
       rhou(L) = rhovol(L-Lb+1) /  0.5d0*( zws(k1t) - zws(k1-1) + zws(k2t) - zws(k2-1) )
    endif
    rhovol(L-Lb+1) = rhovol(L-Lb+1)*dx(LL)

    rv1 = rvdn(k1)
    rv2 = rvdn(k2)
    gr1 = grn (k1)
    gr2 = grn (k2)

   if (jabaroczlaybed > 0) then ! keyword will vanish later
    if (L == Lb) then
       if (kmxn(ln(1,LL)) > kmxn(ln(2,LL)) .or. kmxn(ln(1,LL)) < kmxn(ln(2,LL))) then ! extrapolate at 'bed' layer of deepest side

          if ( kmxn(ln(1,LL)) > kmxn(ln(2,LL)) ) then ! k1=deepest
             k  = k1 ; kt  = ktop(ln(1,LL))
             kz = k2 !; ktz = ktop(ln(2,LL))
          else
             k  = k2 ; kt  = ktop(ln(2,LL))
             kz = k1 !; ktz = ktop(ln(1,LL))
          endif

          fzu   = (zws(kz+1) - zws(kz)) / (zws(kz+1) - zws(kz-1)) ; fzd = 1d0 - fzu
          dzz   =  zws(kz) - zws(kz-1)

          rhow1 = fzu*rho(k+1) + fzd*rho(k)
          rhow0 = 2d0*rho(k) - rhow1
          rhow1 = rhow1 - rhomean
          rhow0 = rhow0 - rhomean
          rvk   =   rvdn(k+1) + 0.5d0*dzz*    ( rhow1 + rhow0 )
          grk   = ( rvdn(k+1) + 0.5d0*dzz*( 2d0*rhow1 + rhow0 )/3d0 )*dzz

          if ( kmxn(ln(1,LL)) > kmxn(ln(2,LL)) ) then ! k1=deepest
             rv1 = rvk ; gr1 = grk
          else
             rv2 = rvk ; gr2 = grk
          endif
       endif
    endif
   endif

    gr3 = 0.5d0*( rv1 + rv2 )*( zws(k1-1) - zws(k2-1) )

    gradpu(L-Lb+1) = gr1 - gr2 + gr3
    if (L > Lb ) then
       gradpu(L-Lb) = gradpu(L-Lb) - gr3   ! ceiling of ff# downstairs neighbours
    endif
 enddo

 ! this last piece is identical to addbaroc2, that will be removed at some moment
 if (jabaroctimeint == 3) then                                           ! original AB implementation

    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl    = ag*gradpu(L-Lb+1) /rhovol(L-Lb+1)                  !
           adve(L)   = adve(L) - 1.5d0*barocl + 0.5d0*dpbdx0(L)
           dpbdx0(L) = barocL
        endif
    enddo

 else if (abs(jabaroctimeint) == 4) then                                 ! AB + better drying flooding

    ft = 0.5d0*dts/dtprev
    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1) / rhovol(L-Lb+1)

           if (dpbdx0(L) .ne. 0d0) then
               adve(L)   = adve(L) - (1d0+ft)*barocl + ft*dpbdx0(L)
           else
               adve(L)   = adve(L) - barocl
           endif
           dpbdx0(L) = barocL
       endif
    enddo

    do L = Lt+1,Lb+kmxL(LL)-1
       dpbdx0(L) = 0d0
    enddo

 else

    do L = Lb, Lt
        if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)                     !  Explicit
           adve(L) = adve(L) - barocl
       endif
    enddo

 endif

 end subroutine addbarocL

 subroutine addbarocLrho_w(LL,Lb,Lt)
 use m_flowgeom
 use m_flow
 use m_flowtimes
 use m_transport, only : NUMCONST, ISALT, ITEMP, ISED1, ISEDN, ITRA1, ITRAN, ITRAN0, constituents
 use m_physcoef , only : rhomean


 implicit none
 integer, intent(in) :: LL,Lb,Lt

 integer             :: L, k1, k2, k1t, k2t, k, kt, i, kz, ktz
 double precision    :: gradpu(kmxx), rhovol(kmxx), gr3, barocL, ft
 double precision    :: rv1, rv2, gr1, gr2, rvk, grk, saw0, saw1, tmw0, tmw1, fzu, fzd, dzz, rv0, rhow0, pdb, p0d
 double precision  , external :: densfm


 do L = Lb,Lt
    k1 = ln(1,L) ; k1t = k1
    k2 = ln(2,L) ; k2t = k2
    if (L == Lt) then
       k1t = ktop(ln(1,LL)) ; k2t = ktop(ln(2,LL))
    endif

    rhovol(L-Lb+1) = 0.5d0*( rvdn(k1t) - rvdn(k1-1) + (zws(k1t) - zws(k1-1))*rhomean  +    &
                             rvdn(k2t) - rvdn(k2-1) + (zws(k2t) - zws(k2-1))*rhomean  )
    if (jarhoxu > 0) then
       rhou(L) = rhovol(L-Lb+1) /  0.5d0*( zws(k1t) - zws(k1-1) + zws(k2t) - zws(k2-1) )
    endif
    rhovol(L-Lb+1) = rhovol(L-Lb+1)*dx(LL)

    rv1 = rvdn(k1)
    rv2 = rvdn(k2)
    gr1 = grn (k1)
    gr2 = grn (k2)

   if (jabaroczlaybed > 0) then ! keyword will vanish later
    if (L == Lb) then
       if (kmxn(ln(1,LL)) > kmxn(ln(2,LL)) .or. kmxn(ln(1,LL)) < kmxn(ln(2,LL))) then ! extrapolate at 'bed' layer of deepest side
          if ( kmxn(ln(1,LL)) > kmxn(ln(2,LL)) ) then ! k1=deepest
             k  = k1 ; kt  = ktop(ln(1,LL))
             kz = k2 ; ktz = ktop(ln(2,LL))
          else
             k  = k2 ; kt  = ktop(ln(2,LL))
             kz = k1 ; ktz = ktop(ln(1,LL))
          endif

          fzu  = (zws(kz+1) - zws(kz)) / (zws(kz+1) - zws(kz-1)) ; fzd = 1d0 - fzu
          dzz  =  zws(kz) - zws(kz-1)

          saw1 = fzu*constituents(isalt,k+1) + fzd*constituents(isalt,k)
          tmw1 = fzu*constituents(itemp,k+1) + fzd*constituents(itemp,k)
          saw0 = 2d0*constituents(isalt,k) - saw1
          tmw0 = 2d0*constituents(itemp,k) - tmw1

          if (idensform < 10) then
             rhow0 = densfm(saw0,tmw0,0d0) - rhomean
          else
             pdb  = ( zws(ktz) - zws(kz-1) )*rhomean
             rvk  = rvdn(k+1) + 0.5d0*dzz*( rhosww(k) + rhosww(k-1) )
             do i = 1,maxitpresdens
                p0d   = ag*(rvk + pdb)                                            ! total pressure
                rhow0 = densfm(saw0,tmw0,p0d) - rhomean
                rvk   =  rvdn(k+1) + 0.5d0*dzz*    ( rhosww(k) + rhow0 )
             enddo
          endif
          rvk   =   rvdn(k+1) + 0.5d0*dzz*    ( rhosww(k) + rhow0 )
          grk   = ( rvdn(k+1) + 0.5d0*dzz*( 2d0*rhosww(k) + rhow0 )/3d0 )*dzz

          if ( kmxn(ln(1,LL)) > kmxn(ln(2,LL)) ) then ! k1=deepest
             rv1 = rvk ; gr1 = grk
          else
             rv2 = rvk ; gr2 = grk
          endif
       endif
    endif
   endif

    gr3            = 0.5d0*( rv1 + rv2 )*(zws(k1-1) - zws(k2-1))
    gradpu(L-Lb+1) = gr1 - gr2 + gr3
    if (L > Lb ) then
       gradpu(L-Lb) = gradpu(L-Lb)     - gr3            ! ceiling of ff# downstairs neighbours
    endif
 enddo

 ! this last piece is identical to addbaroc2, that will be removed at some moment
 if (jabaroctimeint == 3) then                                           ! original AB implementation

    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl    = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)                  !
           adve(L)   = adve(L) - 1.5d0*barocl + 0.5d0*dpbdx0(L)
           dpbdx0(L) = barocL
        endif
    enddo

 else if (abs(jabaroctimeint) == 4) then                                 ! AB + better drying flooding

    ft = 0.5d0*dts/dtprev
    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)
           if (dpbdx0(L) .ne. 0d0) then
               adve(L)   = adve(L) - (1d0+ft)*barocl + ft*dpbdx0(L)
           else
               adve(L)   = adve(L) - barocl
           endif
           dpbdx0(L) = barocL
       endif
    enddo

    do L = Lt+1,Lb+kmxL(LL)-1
       dpbdx0(L) = 0d0
    enddo

 else

    do L = Lb, Lt
        if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)                     !  Explicit
           adve(L) = adve(L) - barocl
       endif
    enddo

 endif

 end subroutine addbarocLrho_w

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

 subroutine addbarocLorg(LL,Lb,Lt)
 use m_flowgeom
 use m_flow
 use m_flowtimes

 implicit none
 integer, intent(in) :: LL,Lb,Lt

 integer             :: L, k1, k2, k1t, k2t
 double precision    :: gradpu(kmxx), rhovol(kmxx), gr3, barocL, ft

 do L = Lb,Lt
    k1 = ln(1,L) ; k1t = k1
    k2 = ln(2,L) ; k2t = k2
    if (L == Lt) then
       k1t = ktop(ln(1,LL)) ; k2t = ktop(ln(2,LL))
    endif
    rhovol(L-Lb+1) = 0.5d0*( zws(k1t) - zws(k1-1) + zws(k2t) - zws(k2-1) )*dx(LL)*0.5d0*(rho(k1) + rho(k2)) ! write in rvdn
    gr3            = 0.5d0*( rvdn(k1) + rvdn(k2) )*(zws(k1-1) - zws(k2-1))
    gradpu(L-Lb+1) = grn(k1) - grn(k2) + gr3
    if (L > Lb ) then
       gradpu(L-Lb) = gradpu(L-Lb)     - gr3            ! ceiling of ff# downstairs neighbours
    endif
 enddo

 ! this last piece is identical to addbaroc2, that will be removed at some moment
 if (jabaroctimeint == 3) then                                           ! original AB implementation

    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl    = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)                  !
           adve(L)   = adve(L) - 1.5d0*barocl + 0.5d0*dpbdx0(L)
           dpbdx0(L) = barocL
        endif
    enddo

 else if (abs(jabaroctimeint) == 4) then                                 ! AB + better drying flooding

    ft = 0.5d0*dts/dtprev
    do L = Lb, Lt
       if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)
           if (dpbdx0(L) .ne. 0d0) then
               adve(L)   = adve(L) - (1d0+ft)*barocl + ft*dpbdx0(L)
           else
               adve(L)   = adve(L) - barocl
           endif
           dpbdx0(L) = barocL
       endif
    enddo

    do L = Lt+1,Lb+kmxL(LL)-1
       dpbdx0(L) = 0d0
    enddo

 else

    do L = Lb, Lt
        if (rhovol(L-Lb+1) > 0d0) then
           barocl  = ag*gradpu(L-Lb+1)/rhovol(L-Lb+1)                     !  Explicit
           adve(L) = adve(L) - barocl
       endif
    enddo

 endif

 end subroutine addbarocLorg
