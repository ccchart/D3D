!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
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

subroutine flgsfm( n, ng, L, firstiter, jarea)
use m_flowgeom
!!--description-----------------------------------------------------------------
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    ! use cpluv
    ! use m_strucs
    ! use ident

    use m_strucs
    use m_flow

    implicit none
!
! Global variables
!
    integer, intent(in)  :: n          !< general structure point n
    integer, intent(in)  :: ng         !< is a member of general structure sigal ng
    integer, intent(in)  :: L          !< Flow link number, signed! If L < 0 then flow link is in opposite direction than structure left-right orientation.
    logical, intent(in)  :: firstiter
    logical              :: jarea


!
!
! Local variables
!
    integer                        :: il, ir, k1, k2, kL, kR, m, Lf
    integer                        :: L0
    logical                        :: velheight
    double precision               :: cgd
    double precision               :: cgf
    double precision               :: crest
    double precision               :: cwd
    double precision               :: cwf
    double precision               :: dg
    double precision               :: ds
    double precision               :: ds1
    double precision               :: ds2
    double precision               :: hdsb
    double precision               :: husb
    double precision               :: lambda
    double precision               :: mugf
    double precision               :: relax
    double precision               :: rholeft
    double precision               :: rhoright
    double precision               :: strdamf
    double precision               :: teken, tekenstr
    double precision               :: ud
    double precision               :: uu
    double precision               :: w2
    double precision               :: wsd
    double precision               :: wstr
    double precision               :: zb2
    double precision               :: zs, gateloweredgelevel, gatedoorheight
    double precision               :: DsL
    double precision               :: gatefraction, fu_sav, ru_sav, au_sav

!
!! executable statements -------------------------------------------------------
!
    !
    !=======================================================================
    !                      Deltares
    !                One-Two Dimensional Modelling System
    !                           S O B E K
    !
    ! Subsystem:          Flow Module
    !
    ! Programmer:         J.Kuipers
    !
    ! Module:             FLGS (FLow General Structure)
    !
    ! Module description: In subroutine FLGS the QH-relationship for a
    !                     general structure will be transformed to a
    !                     linearized equation
    !
    !
    ! Parameters:
    ! NR NAME              IO DESCRIPTION
    !  2 il                I  Grid point on left side of structure (lower
    !                         index).
    !  3 ir                I  Grid point on right side of structure (upper
    !                         index).
    !  8 jarea             I  If True then claculate only area
    !  1 m                 I  Grid index of structure
    !  4 istru             I  Number of structure.
    !  7 firstiter         I  True in case of first iteration step.

    ! Subprogram calls:
    ! NAME    DESCRIPTION
    ! flgtar  FLow get General sTructure ARguments
    ! flupdg  FLow UP/Downstream near General structure
    ! errmsg  generate ERRer MeSsaGe to log file
    ! flqhgs  FLow QH relation for General Structure
    !=======================================================================
    !     Include Pluvius data space
    !     Include identifiers of objects
    !
    !     Declaration of parameters:
    !
    !
    !     Declaration of local variables:
    !
    !
    !
    !TEM  WRITE (11,*) 'Call structure',istru,'(',m,il,ir,istru,')'

    Lf = abs(L)
    ! NOTE: Since a single general structure may be crossed by multiple flow links,
    ! pay attention to proper directions: structure parameters are typically determined
    ! by the structure's left-right direction, whereas upwinding and furu-computations
    ! are typically in the flow link's 1-2 direction.
    k1 = ln(1,Lf) ; k2 = ln(2,Lf)      ! 1 -> 2 flow link direction
    kL = kcgen(1,n) ; kR = kcgen(2,n)  ! L -> R structure direction

    m  = L
    il = k1
    ir = k2
    L0 = n - L1cgensg(ng)+1

    zs                 = min  ( bob(1,Lf), bob(2,Lf) )         ! == zcgen(3*ng - 2) crest/silllevel
    gateloweredgelevel = generalstruc(ng)%gateheightonlink(L0) ! == zcgen(3*ng - 1) under gate door and infinity in open part.
    gatefraction = generalstruc(ng)%gateclosedfractiononlink(L0)

    ! TODO: RTC: AvD/Herman: hier ook wu'tjes en zb1-tjes etc gaan zetten, voordat we de flupd/flgtar-subroutines gaan callen?
    ! Velheight is always true for river structures
    ! velheight = istrtyp(7, istru)==1
    velheight = .true.

    relax = 1.0D0
    au(Lf) = 0d0 ; fu(Lf) = 0d0 ; ru(Lf) = 0d0
    !
    ! ng instead of istru

    dg = gateloweredgelevel - zs

    call flupdofm(m, il, ir, ng, velheight, rholeft, rhoright, crest, husb, hdsb,     &
                  uu, ud, teken, relax)

    gatedoorheight = 0d0

    tekenstr = teken*sign(1, L) ! if flow link abs(L) is in opposite orientation to the structure's orientation, then negate the just computed upwind (flow) teken.

    if (husb > zs) then
       call flgtarfm(ng, L0, wu(Lf), bl(kL), bl(kR), tekenstr, zs, wstr, w2, wsd, zb2, dg, ds1, ds2, cgf, cgd,   &
                     cwf, cwd, mugf, lambda, strdamf, gatedoorheight)

       DsL   = s1(k2) - s1(k1)
       u1(Lf) = rusav(1,n) - fusav(1,n)*DsL ; u0(Lf) = u1(Lf) ; q1(Lf) = ausav(1,n)*u1(Lf)
       call flqhgsfm(Lf, teken, husb, hdsb, uu, zs, gatefraction*wstr, w2, wsd, zb2, ds1, ds2, dg,  &
                     cgf, cgd, cwf, cwd, mugf, lambda, strdamf, jarea, ds)
       fusav(1,n) = fu(Lf) ; rusav(1,n) = ru(Lf) ; ausav(1,n) = au(Lf)
    else
      fusav(1,n) = 0d0
      rusav(1,n) = 0d0
      ausav(1,n) = 0d0
    endif

    if (gatedoorheight > 0d0) then  ! now add water overflowing top of gate
       zs = gateloweredgelevel + gatedoorheight
       if (husb > zs) then          ! husb = upwind waterlevel instead of height
          dg    = 1d9               ! sky is the limit, this gate fully open
          u1(Lf) = rusav(2,n) - fusav(2,n)*dsL ; u0(Lf) = u1(Lf) ; q1(Lf) = ausav(2,n)*u1(Lf)
          call flgtarfm(ng, L0, wu(Lf), bl(kL), bl(kR), tekenstr, zs, wstr, w2, wsd, zb2, dg, ds1, ds2, cgf, cgd,   &
                        cwf, cwd, mugf, lambda, strdamf, gatedoorheight)
          call flqhgsfm(Lf, teken, husb, hdsb, uu, zs, wstr, w2, wsd, zb2, ds1, ds2, dg,  &
                        cgf, cgd, cwf, cwd, mugf, lambda, strdamf, jarea, ds)
          fusav(2,n) = fu(Lf) ; rusav(2,n) = ru(Lf) ; ausav(2,n) = au(Lf)

       else
          fusav(2,n) = 0d0
          rusav(2,n) = 0d0
          ausav(2,n) = 0d0
       endif
    else
       fusav(2,n) = 0d0
       rusav(2,n) = 0d0
       ausav(2,n) = 0d0
    endif

    zs                 = min  ( bob(1,Lf), bob(2,Lf) )         ! == zcgen(3*ng - 2) crest/silllevel

    if ( husb >= zs .and. gatefraction < 1d0) then
       fu_sav = fu(Lf)
       ru_sav = ru(Lf)
       au_sav = au(Lf)

       zs =  min  ( bob(1,Lf), bob(2,Lf) )
       dg = huge(1d0)
       u1(Lf) = rusav(3,n) - fusav(3,n)*dsL ; u0(Lf) = u1(Lf) ; q1(Lf) = ausav(3,n)*u1(Lf)
       call flgtarfm(ng, L0, wu(Lf), bl(kL), bl(kR), tekenstr, zs, wstr, w2, wsd, zb2, dg, ds1, ds2, cgf, cgd,   &
                     cwf, cwd, mugf, lambda, strdamf, gatedoorheight)
       call flqhgsfm(Lf, teken, husb, hdsb, uu, zs, (1d0-gatefraction)*wstr, w2, wsd, zb2, ds1, ds2, dg,  &
                     cgf, cgd, cwf, cwd, mugf, lambda, strdamf, jarea, ds)
       fusav(3,n) = fu(Lf) ; rusav(3,n) = ru(Lf) ; ausav(3,n) = au(Lf)
    else
       fusav(3,n) = 0d0
       rusav(3,n) = 0d0
       ausav(3,n) = 0d0
    end if

    au(Lf) =  ausav(1,n) + ausav(2,n) + ausav(3,n)
    if (au(Lf) > 0d0) then
       fu(Lf) = (fusav(1, n)*ausav(1, n) + fusav(2, n)*ausav(2, n) + fusav(3, n)*ausav(3, n))/au(Lf)
       ru(Lf) = (rusav(1, n)*ausav(1, n) + rusav(2, n)*ausav(2, n) + rusav(3, n)*ausav(3, n))/au(Lf)
    else
       fu(Lf) = 0d0
       ru(Lf) = 0d0
    endif

    if (au(Lf) == 0d0) then
        hu(Lf) =  0d0
    endif

    ! TEMP = laatste statement
    ! strhis(15, istru) = ds + crest     ! waterlevel on crest
end subroutine flgsfm
