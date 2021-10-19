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

  ! =================================================================================================
  ! =================================================================================================
  subroutine setpillars()
    use m_flowgeom            , only: ndx, lnx, ba, wu, nd
    use m_flowexternalforcings, only: pillar, Cpil
    use m_vegetation          , only: rnveg, diaveg, stemheight
    use gridoperations
    use m_flowparameters      , only: japillar
    use m_crspath
    implicit none
    integer                                     :: i, j, k, L, Lf, La, m, n
    double precision                            :: pi
    integer         , dimension(:), allocatable :: npil
    double precision, dimension(:), allocatable :: cdeq
    double precision, dimension(:), allocatable :: Aeff
    integer         , dimension(:), allocatable :: linktype
    integer                                     :: nPath
    type (tcrspath) , dimension(:), allocatable :: Path
    integer         , dimension(:), allocatable :: idum

    if (allocated(Cpil)) deallocate( Cpil )
    if (japillar == 1) then
      allocate( Cpil(ndx) )
    else if (japillar == 3) then
      allocate( Cpil(lnx) )
    endif

    if (allocated(idum)) deallocate(idum)
    allocate(idum(1))
    idum = 0

    pi = 4.0d0 * atan( 1d0 )

    if( japillar == 2 ) then
      if( allocated( cdeq ) ) deallocate( cdeq, npil )
      allocate( cdeq(ndx), npil(ndx) )
      cdeq = 0d0
      npil = 0
      do m = 1,size(pillar)
        do i = 1,pillar(m)%np
          if( pillar(m)%dia(i) == -999d0 .or. pillar(m)%cd(i) == -999d0 ) cycle
          call incells( pillar(m)%xcor(i), pillar(m)%ycor(i), j )
          if( j == 0 ) cycle
          rnveg(j) = rnveg(j) + pillar(m)%dia(i)**2 * pi * 0.25d0 / ba(j)
          cdeq(j)  = cdeq(j)  + pillar(m)%cd(i) * pillar(m)%dia(i)
          npil(j) = npil(j) + 1
        enddo
      enddo
      do j = 1,ndx
        if( npil(j) == 0 ) cycle
        diaveg(j)  = diaveg(j) + cdeq(j) / npil(j)
        stemheight(j) = 1d30
      enddo
      deallocate( cdeq )
      deallocate( npil )

    elseif( japillar == 1 ) then   ! Delft3D implimentation, but modified version on flow cells
      if (allocated(Aeff) ) deallocate( Aeff, cdeq )
      allocate( Aeff(ndx), cdeq(ndx) )
      do j = 1,ndx
        Aeff(j) = ba(j)
      enddo
      cdeq = 0d0
      do m = 1,size(pillar)
        do i = 1,pillar(m)%np
          if (pillar(m)%dia(i) == -999d0 .or. pillar(m)%cd(i) == -999d0) cycle
          call incells( pillar(m)%xcor(i), pillar(m)%ycor(i), j )
          if (j == 0) cycle
          cdeq(j) = cdeq(j) + pillar(m)%cd(i) * pillar(m)%dia(i)
          Aeff(j) = Aeff(j) - pillar(m)%dia(i)**2 * pi * 0.25d0
        enddo
      enddo
      Cpil = 0d0
      do j = 1,ndx
        if( cdeq(j) == 0 ) cycle
        if( Aeff(j) <= 0d0 ) then
          Cpil(j) = 1d30
          cycle
        endif
        Cpil(j) = cdeq(j) * 0.25d0 / Aeff(j) * sqrt( ba(j) * pi )
      enddo
      deallocate( Aeff )
      deallocate( cdeq )

    else if (japillar == 3) then       ! Based on D3D approach on flow links
      if (allocated(Aeff) ) deallocate( Aeff, cdeq, linktype )
      allocate( Aeff(lnx), cdeq(lnx), linktype(lnx) )
      linktype = 0
      Aeff = wu
      cdeq = 0d0
      do m = 1,size(pillar)
        call pol_to_flowlinks(pillar(m)%xcor, pillar(m)%ycor, pillar(m)%xcor*0d0, pillar(m)%np, nPath, Path)
        do n = 1,nPath
          call crspath_on_flowgeom(Path(n),0,0,1,idum,0,1)
          do L = 1,Path(n)%lnx
            Lf = Path(n)%ln(L)
            La = iabs(Lf)
            linktype(La) = 1
          enddo
        enddo
        do i = 1,pillar(m)%np
          if (pillar(m)%dia(i) == -999d0 .or. pillar(m)%cd(i) == -999d0) cycle
          call incells( pillar(m)%xcor(i), pillar(m)%ycor(i), k )
          if( k == 0 ) cycle
          do L = 1,nd(k)%lnx
            Lf = nd(k)%ln(L)
            La = iabs(Lf)
            if (linktype(La) /= 1) cycle
            cdeq(La) = cdeq(La) + pillar(m)%cd(i) * pillar(m)%dia(i)
            Aeff(La) = Aeff(La) - pillar(m)%dia(i)
          enddo
        enddo
      enddo
      Cpil = 0d0
      do L = 1,lnx
        if( cdeq(L) == 0 ) cycle
        if( Aeff(L) <= 0d0 ) then
          Cpil(L) = 1d30
          cycle
        endif
        Cpil(L)  = cdeq(L) * 0.5d0 * wu(L) / Aeff(L)**2
      enddo
      deallocate( Aeff )
      deallocate( cdeq )
      deallocate( linktype )
    endif

  end subroutine setpillars
  
  ! =================================================================================================
  ! =================================================================================================
  subroutine init_sealock ()
     use m_flowexternalforcings, only: nsealocksg, sealock, numsrc, ksealock, L1sealocksg, L2sealocksg
     use m_flowgeom            , only: bl
     use m_flow                , only: hs, s1
     use gridoperations        , only: incells
     use network_data, only: xzw, yzw
     implicit none
     integer                        :: m, k1, k2, ksea, klake, L, ierr
     double precision               :: area
     double precision, dimension(2) :: xpin, ypin, zpin, dzlin

     interface
        subroutine addsorsin(filename, area, ierr, xpin, ypin, zpin, dzlin, sorsin_name)
            character (len=*), intent(in)  :: filename
            double precision,  intent(in)  :: area
            integer,           intent(out) :: ierr
            double precision, dimension(:), optional, intent(in) :: xpin, ypin, zpin, dzlin
            character (len=*), optional, intent(in)  :: sorsin_name !< Optional custom name for source-sink, if filename is not used, but xpin is used instead.
        end subroutine addsorsin
     end interface
     
     ! use m_flowexternalforcings, only: *sealockblabla
     ! do i=1,numslf
     !
     !
     !   ... call add_sorsin( ..., (/ sealock(i)%xsea_transport, xlake_transport /), (/ .. y... /) ..)
     
     ! end do
    do m = 1,nsealocksg
       call incells (sealock(m)%xsea_probe , sealock(m)%ysea_probe , k1)
       call incells (sealock(m)%xlake_probe, sealock(m)%ylake_probe, k2)
       sealock(m)%ksea_probe  = k1
       sealock(m)%klake_probe = k2
       
       L = L1sealocksg(m)
       ksea = ksealock(1,L) ! for now assume that a sealock is always on 1 gridcell, so use L from L1sealocksg.
       klake = ksealock(2,L)
       
       area = 0d0
       ierr = 0
       xpin(1) = xzw(ksea)
       ypin(1) = yzw(ksea)
       zpin(1)  = bl(ksea)
       dzlin(1) = bl(ksea) + hs(ksea) * 0.3d0   !TODO: hardcoded 30% needs to be as input
       
       xpin(2) = xzw(klake)
       ypin(2) = yzw(klake)
       zpin(2)  = bl(klake)
       dzlin(2) = bl(klake) + hs(klake) * 0.3d0
       
       call addsorsin('', area, ierr, xpin, ypin, zpin, dzlin, trim(sealock(m)%id)//'_s2l')  !TODO: buttom source sink from sea to lake
       sealock(m)%sorsin_index(1) = numsrc ! Store source-sink index sea->lake for this sealock

       xpin(1) = xzw(klake)
       ypin(1) = yzw(klake)
       zpin(1) = bl(klake) + hs(klake) * 0.7d0
       dzlin(1) = s1(klake)
       
       xpin(2) = xzw(ksea)
       ypin(2) = yzw(ksea)
       zpin(2) = bl(ksea) + hs(ksea) * 0.7d0
       dzlin(2) = s1(ksea)
       
       call addsorsin('', area, ierr, xpin, ypin, zpin, dzlin, trim(sealock(m)%id)//'_l2s')
       sealock(m)%sorsin_index(2) = numsrc ! Store source-sink index lake->sea for this sealock

    enddo
  
  end subroutine init_sealock
