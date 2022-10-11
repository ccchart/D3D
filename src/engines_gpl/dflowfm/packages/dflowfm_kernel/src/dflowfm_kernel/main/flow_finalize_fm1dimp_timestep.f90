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

!> Updates the flow variables of FM when the Sobek-RE kernel is used. 

!This function cannot be in the <fm1dimp> module because it uses FM variables and
!while the kernel of FM depends on the <fm1dimp> module, the opposite is not true.

subroutine flow_finalize_fm1dimp_timestep()

!
!MODULES
!

use m_flow, only: s1, u1
use m_flowgeom, only: lnx1d, ndxi, ln, lnx1Db, lnxi
!use unstruc_channel_flow, only: network
!use m_CrossSections, only: CalcConveyance
!use m_flowgeom
!use m_fm_erosed, only: link1, link1sign
use m_f1dimp 

implicit none

!
!DECLARATION
!

!pointer

real, dimension(:)                       , pointer :: x

real, dimension(:,:)                     , pointer :: waoft

double precision, dimension(:,:)         , pointer :: hpack
double precision, dimension(:,:)         , pointer :: qpack

!locals

integer :: L, N, n1, n2, nint, nout

!
!SET POINTERS
!

!dependent on gridpoints 
x      => f1dimppar%x 
waoft  => f1dimppar%waoft 
hpack  => f1dimppar%hpack
qpack  => f1dimppar%qpack

 ! 1:ndx2D, ndx2D+1:ndxi, ndxi+1:ndx1Db, ndx1Db+1:ndx
 ! ^ 2D int ^ 1D int      ^ 1D bnd       ^ 2D bnd ^ total
do N=1,ndxi !internal cell centres
    s1(N)=hpack(N,3)
enddo
   
do L=1,lnx1d !internal links
    n1 = ln(1,L) 
    n2 = ln(2,L)
    u1(L)=0.5*qpack(n1,3)/waoft(n1,3)+0.5*qpack(n2,3)/waoft(n2,3)
enddo
            
do L=lnxi+1,lnx1Db !boundary links
    !FM1DIMP2DO: we could create a variable with this mapping to prevent computation of max, min and x(nint) every timestep
    n1 = ln(1,L) 
    n2 = ln(2,L)
    nint=min(n1,n2) !from the two cells that this link connects, the minimum is internal, and hence we have data
    nout=max(n1,n2) !from the two cells that this link connects, the maximum is extrernal, and it is the one in which we have to set the water level
    
    !we could have a better reconstruction of <u1(L)> with the slope of the previous value rather than just copying the value.
    
    !not sure if x=0 is enough to order the branch. Another option is to get the next internal cell connected to the identified
    !internal cell and compute <dx>. If positive it is in the direction of the flow and viceversa.
    
    !FM1DIMP2DO: this is prone to error and tricky. I am assuming that the upstream end, where the discharge is specified, is a chainage 0
    !The difficulty is to set a direction of the branch that defines what it means that the velocity of SRE is positive. 
    if (x(nint).eq.0) then !upstream
        u1(L)=qpack(nint,3)/waoft(nint,3) 
    else ! downstream
        u1(L)=-qpack(nint,3)/waoft(nint,3) !we could have a better reconstruction with the slope of the previous value
    endif
    
    s1(nout)=s1(nint)
enddo

end subroutine flow_finalize_fm1dimp_timestep
