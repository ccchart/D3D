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

!> Updates the variables of flow1d implicit solver that
!change every time step
subroutine fm1dimp_update_network(iresult)

!use m_flowparameters
use m_f1dimp
use unstruc_channel_flow, only: network
use m_CrossSections, only: createTablesForTabulatedProfile
use m_fm_erosed, only: ndx_mor, lnx_mor
use m_oned_functions, only: gridpoint2cross
use m_flow, only: hu
use m_flowgeom, only: lnx 

!use unstruc_messages

implicit none


!
!pointer
!

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: kd
integer :: k2
integer :: idx_crs, idx_sre

!pointers
integer                                  , pointer :: ngrid   

integer, dimension(:)                    , pointer :: grd_ghost_link_closest
integer, dimension(:)                    , pointer :: grd_sre_cs 

real, dimension(:,:)                     , pointer :: wft
real, dimension(:,:)                     , pointer :: aft
real, dimension(:,:)                     , pointer :: wtt
real, dimension(:,:)                     , pointer :: att
real, dimension(:,:)                     , pointer :: of

double precision, dimension(:)           , pointer :: bedlevel

double precision, dimension(:,:)         , pointer :: hlev

!!
!! POINT
!!

ngrid                  => f1dimppar%ngrid

grd_ghost_link_closest => f1dimppar%grd_ghost_link_closest
grd_sre_cs             => f1dimppar%grd_sre_cs

wft                    => f1dimppar%wft
aft                    => f1dimppar%aft
wtt                    => f1dimppar%wtt
att                    => f1dimppar%att
of                     => f1dimppar%of

bedlevel               => f1dimppar%bedlevel

hlev                   => f1dimppar%hlev




!!
!! CALC
!!

iresult=0

!do kd=1,ndx_mor
do kd=1,ngrid
    
    !idx_fm=f1dimppar%grd_sre_fm(k) !index of the global grid point in fm for the global gridpoint <k> in SRE
    !idx_crs=network%ADM%gpnt2cross(idx_fm)%C1 !< index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node 		
    
    !It is not possible to use <gridpoint2cross> via the fm index to get the cross-section index
    !because for the case of a junction with just two branches there are two cross-sections and
    !two SRE gridpoints but just one FM index. In this situation the SRE gridpoint and associated
    !cross-section are not updated. As a consequence, an independent mapping variable (<grd_sre_cs>)
    !is needed
    
    !skip the boundary-ghost flownodes, for which there is no cross-section
    !idx_sre=f1dimppar%grd_fm_sre(kd)
    !if (gridpoint2cross(kd)%NUM_CROSS_SECTIONS.eq.0) then
    !    cycle
    !endif
    !idx_crs=gridpoint2cross(kd)%cross(1)
    
    idx_sre=kd
    idx_crs=grd_sre_cs(idx_sre)
    
    
    !update cross-section flow variables after bed level changes
    !FM1DIMP2DO: remove debug
    !if (idx_crs.eq.1) then
    !write(42,*) network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(5), network%CRS%CROSS(idx_crs)%TABDEF%height(1)
    !endif
    
    call createTablesForTabulatedProfile(network%CRS%CROSS(idx_crs)%TABDEF)
    
    !FM1DIMP2DO: remove debug
    !if (idx_crs.eq.1) then
    !write(42,*) network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(5), network%CRS%CROSS(idx_crs)%TABDEF%height(1)
    !endif
    
    !copy values to <fm1dimppar>
    !FM1DIMP2DO: move to a subroutine?
    
    !nlev
    f1dimppar%nlev(idx_sre)=network%CRS%CROSS(idx_crs)%TABDEF%LEVELSCOUNT
    do k2=1,f1dimppar%nlev(idx_sre)
        wft (idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWWIDTH(k2) 
        aft (idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(k2)  
        wtt (idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALWIDTH(k2) 
        att (idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALAREA(k2)    
        of  (idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%WETPERIMETER(k2)     
        hlev(idx_sre,k2)=network%CRS%CROSS(idx_crs)%TABDEF%HEIGHT(k2)
    end do !k2

    !bed level
    bedlevel(idx_sre)=network%crs%cross(idx_crs)%bedLevel
    
        !FM1DIMP2DO: remove debug
    !write(42,*) idx_sre, idx_crs
    !write(42,*) f1dimppar%wft (idx_sre,:)
    !write(42,*) f1dimppar%aft (idx_sre,:)
    !write(42,*) f1dimppar%wtt (idx_sre,:)
    !write(42,*) f1dimppar%att (idx_sre,:)
    !write(42,*) f1dimppar%of (idx_sre,:)
    !write(42,*) f1dimppar%hlev (idx_sre,:)
    !write(42,*) f1dimppar%bedlevel (idx_sre)
    
enddo !kd

!in <fm_upwbed> there is a threshold check on `hu(Lf)>epshu`. Here we
!copy the value of the closest link to the ghost link for passing
!this threshold check. 
call reallocate_fill(hu      ,grd_ghost_link_closest,lnx,lnx_mor)

end subroutine fm1dimp_update_network

