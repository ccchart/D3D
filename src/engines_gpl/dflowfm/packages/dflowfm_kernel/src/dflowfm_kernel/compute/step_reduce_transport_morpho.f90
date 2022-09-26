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

 subroutine step_reduce_transport_morpho(key)        ! do a flow timestep dts guus, reduce once, then elimin conjugate grad substi
 use m_flow                                          ! when entering this subroutine, s1=s0, u1=u0, etc
 use m_flowgeom
 use m_sediment, only: stm_included, stmpar, mtd
 use Timers
 use m_flowtimes
 use m_sferic
 use m_wind
 use m_reduce
 use m_ship
 use m_partitioninfo
 use m_timer
 use m_xbeach_data
 use MessageHandling
 use m_sobekdfm
 use unstruc_display
 use m_waves
 use m_subsidence

 implicit none

 integer :: ndraw
 COMMON /DRAWTHIS/  ndraw(50)

 integer            :: key, jposhchk_sav, LL, L, k1,k2, itype
 integer            :: ja, k, ierror, n, kt, num, js1, noddifmaxlevm, nsiz
 character (len=40) :: tex
 logical            :: firstnniteration
 double precision   :: wave_tnow, wave_tstop, t0, t1, dif, difmaxlevm
 double precision   :: hw,tw, uorbi,rkw,ustt,hh,cs,sn,thresh

 character(len=128) :: msg

!-----------------------------------------------------------------------------------------------
 numnodneg = 0
 if (wrwaqon.and.allocated(qsrcwaq)) then
    qsrcwaq0 = qsrcwaq ! store current cumulative qsrc for waq at the beginning of this time step
 end if

 111 continue

!-----------------------------------------------------------------------------------------------
! TODO: AvD: consider moving everything below to flow_finalize single_timestep?
 call setkbotktop(0)                                 ! bottom and top layer indices and new sigma distribution

 call u1q1()                                         ! the vertical flux qw depends on new sigma => after setkbotktop
 call compute_q_total_1d2d()

 !if ( jacheckmonitor.eq.1 ) then
 !   call comp_checkmonitor()
 !end if

 if ( itstep.eq.4 ) then   ! explicit time-step
    call update_s_explicit()
 end if
 hs = s1-bl
 hs = max(hs,0d0)


  if ((jawave==3 .or. jawave==6) .and. .not. flowWithoutWaves) then
    ! Normal situation: use wave info in FLOW
    ! 3D not implementend
    if( kmx == 0 ) then
       call wave_comp_stokes_velocities()
       call wave_uorbrlabda()                                          ! hwav gets depth-limited here
       call tauwave()
    end if
    call setwavfu()
    call setwavmubnd()
 end if

  if ((jawave==3 .or. jawave==6) .and. flowWithoutWaves) then
    ! Exceptional situation: use wave info not in FLOW, only in WAQ
    ! Only compute uorb
    ! Works both for 2D and 3D
    call wave_uorbrlabda()                       ! hwav gets depth-limited here
  end if

  if (jawave.eq.4 .and. jajre.eq.1 .and. .not. flowWithoutWaves) then
    if (swave.eq.1 ) then
       call xbeach_waves(ierror)
    endif
    call tauwave()
    if ( jaGUI.eq.1 ) then                                          ! this part is for online visualisation
       if (ntek > 0) then
          if (mod(int(dnt_user),ntek) .eq. 0) then
             call wave_makeplotvars()                                ! Potentially only at ntek interval
          end if
       endif
    endif
    if (jamombal==1) then
       call xbeach_mombalance()
    endif
  end if

  if (jawave==5 .and. .not. flowWithoutWaves) then
    if (kmx==0) then
       do k=1,ndx
          hwav(k) = min(hwavuni, gammax*(s1(k)-bl(k)))
       enddo
       do L=1,lnx
          k1=ln(1,L); k2=ln(2,L)
          hh = hu(L); hw=0.5d0*(hwav(k1)+hwav(k2));tw=.5d0*(twav(k1)+twav(k2))
          cs = 0.5*(cos(phiwav(k1)*dg2rd)+cos(phiwav(k2)*dg2rd))
          sn = 0.5*(sin(phiwav(k1)*dg2rd)+sin(phiwav(k2)*dg2rd))
          call tauwavehk(hw, tw, hh, uorbi, rkw, ustt)
          ustokes(L) = ustt*(csu(L)*cs + snu(L)*sn)
          vstokes(L) = ustt*(-snu(L)*cs + csu(L)*sn)
       enddo
       do k=1,ndx
          call tauwavehk(hwav(k), twav(k), hs(k), uorbi, rkw, ustt)
          rlabda(k) = rkw; uorb(k) = uorbi
       enddo
       call tauwave()
    endif
 endif

 if (jased > 0 .and. stm_included) then
    if ( jatimer.eq.1 ) call starttimer(IEROSED)
    if (jawave==0) then
       call settaubxu_nowave()         ! set taubxu for no wave conditions BEFORE erosed
    endif
    !
!    call setucxucyucxuucyu()
    call setucxucy_mor (u1)
    call fm_fallve()                   ! update fall velocities
    call fm_erosed()                   ! source/sink, bedload/total load
    if ( jatimer.eq.1 ) call stoptimer(IEROSED)
 end if

 ! secondary flow
 if ( jasecflow > 0 .and. kmx == 0 ) then
    call get_curvature()
    if( jaequili > 0 ) then
       call equili_spiralintensity()
    endif
 end if

 !SPvdP: timestep is now based on u0, q0
 !       transport is with u1,q1 with timestep based on u0,q0
 if ( jatimer.eq.1 ) call starttimer(ITRANSPORT)
 call transport()
 if ( jatimer.eq.1 ) call stoptimer (ITRANSPORT)

 !update particles
 call update_part()

 if (jased > 0 .and. stm_included) then
    call fm_bott3d() ! bottom update
 endif

 if (jasubsupl>0) then
    call apply_subsupl()
 endif

 if ((jased > 0 .and. stm_included).or.(jasubsupl>0)) then
    call setbobs()   ! adjust administration - This option only works for ibedlevtyp = 1, otherwise original bed level [bl] is overwritten to original value
    if (jasubsupl>0) then
       call subsupl_update_s1()
    end if
    call volsur()                     ! update volumes 2d
    if (kmx>0) then
       call setkbotktop(0)            ! and 3D for cell volumes
    endif
 end if

 ! Moved to flow_finalize_single_timestep: call flow_f0isf1()                                  ! mass balance and vol0 = vol1

 if (layertype > 1 .and. kmx.gt.0 ) then

     ! ln = ln0 ! was ok.

 endif

 end subroutine step_reduce_transport_morpho
