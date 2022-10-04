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
use m_CrossSections, only: CalcConveyance
!use unstruc_messages

implicit none


!
!pointer
!

!logical                                  , pointer :: lconv                   
!logical                                  , pointer :: steady    
!                                         
!integer                                  , pointer :: flitmx                 
!integer                                  , pointer :: iterbc                 
!integer                                  , pointer :: ngrid   
!integer                                  , pointer :: ngridm   
!integer                                  , pointer :: nbran   
!integer                                  , pointer :: maxlev
!integer                                  , pointer :: nnode
!integer                                  , pointer :: nhstat
!integer                                  , pointer :: nqstat
!integer                                  , pointer :: maxtab
!integer                                  , pointer :: ntabm
!integer                                  , pointer :: nbrnod
!
!integer, dimension(:)                    , pointer :: nlev
!integer, dimension(:)                    , pointer :: numnod
!
!integer, dimension(:,:)                  , pointer :: branch
!integer, dimension(:,:)                  , pointer :: bfrict
!integer, dimension(:,:)                  , pointer :: hbdpar
!integer, dimension(:,:)                  , pointer :: qbdpar
!integer, dimension(:,:)                  , pointer :: ntab
!integer, dimension(:,:)                  , pointer :: node
!integer, dimension(:,:)                  , pointer :: nodnod
!
!real                                     , pointer :: g
!real                                     , pointer :: psi                    
!real                                     , pointer :: theta                  
!real                                     , pointer :: epsh                   
!real                                     , pointer :: epsq                   
!real                                     , pointer :: rhow                   
!real                                     , pointer :: omega                  
!real                                     , pointer :: epsqrl                 
!real                                     , pointer :: lambda                 
!real                                     , pointer :: relstr                 
!real                                     , pointer :: dhstru                 
!real                                     , pointer :: cflpse                               
!real                                     , pointer :: overlp                 
!real                                     , pointer :: omcfl                  
!real                                     , pointer :: dhtyp                  
!real                                     , pointer :: exrstp     
!
!real, dimension(:)                       , pointer :: table
!real, dimension(:)                       , pointer :: x
!
!real, dimension(:,:)                     , pointer :: bfricp
!real, dimension(:,:)                     , pointer :: wft
!real, dimension(:,:)                     , pointer :: aft
!real, dimension(:,:)                     , pointer :: wtt
!real, dimension(:,:)                     , pointer :: att
!real, dimension(:,:)                     , pointer :: of
!real, dimension(:,:)                     , pointer :: waoft
!
!double precision                         , pointer :: time
!double precision                         , pointer :: dtf
!double precision                         , pointer :: resid
!
!double precision, dimension(:,:)         , pointer :: hpack
!double precision, dimension(:,:)         , pointer :: qpack
!double precision, dimension(:,:)         , pointer :: hlev

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: k
integer :: k2
integer :: idx_crs

real :: swaoft

!!
!! CALC
!!

iresult=0

do k=1,f1dimppar%ngrid
    
    idx_crs=network%ADM%gpnt2cross(k)%C1 !< index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node 		
    
    !update cross-section flow variables after bed level changes
    call CalcConveyance(network%crs%cross(idx_crs))

    !copy values to <fm1dimppar>
    !FM1DIMP2DO: move to a subroutine?
    
    !nlev
    f1dimppar%nlev(k)=network%CRS%CROSS(idx_crs)%TABDEF%LEVELSCOUNT
    do k2=1,f1dimppar%nlev(k)
        f1dimppar%wft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWWIDTH(k2) 
        f1dimppar%aft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(k2)  
        f1dimppar%wtt(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALWIDTH(k2) 
        !FM1DIMP2DO: deal with (at least error) case of rectangular cross-section with only two elevation points. 
        !In this case, the area for the low point is 0. This causes a huge interpolation error in SRE. 
        f1dimppar%att(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALAREA(k2)    
        f1dimppar%of(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%WETPERIMETER(k2)     
        f1dimppar%hlev(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%HEIGHT(k2)
    end do !k2

end do !k

end subroutine fm1dimp_update_network

