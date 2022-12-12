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
!use unstruc_messages

implicit none


!
!pointer
!

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: k
integer :: k2
integer :: idx_crs!, idx_fm

!!
!! CALC
!!

iresult=0

do k=1,f1dimppar%ngrid
    
    !idx_fm=f1dimppar%grd_sre_fm(k) !index of the global grid point in fm for the global gridpoint <k> in SRE
    !idx_crs=network%ADM%gpnt2cross(idx_fm)%C1 !< index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node 		
    
    idx_crs=f1dimppar%idx_cs(k)
    
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
    f1dimppar%nlev(k)=network%CRS%CROSS(idx_crs)%TABDEF%LEVELSCOUNT
    do k2=1,f1dimppar%nlev(k)
        f1dimppar%wft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWWIDTH(k2) 
        f1dimppar%aft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(k2)  
        f1dimppar%wtt(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALWIDTH(k2) 
        f1dimppar%att(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALAREA(k2)    
        f1dimppar%of(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%WETPERIMETER(k2)     
        f1dimppar%hlev(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%HEIGHT(k2)
    end do !k2

    !bed level
    f1dimppar%bedlevel(k)=network%crs%cross(idx_crs)%bedLevel
    
end do !k

end subroutine fm1dimp_update_network

