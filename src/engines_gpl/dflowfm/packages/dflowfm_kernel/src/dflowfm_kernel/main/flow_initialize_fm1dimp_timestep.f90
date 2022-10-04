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

subroutine flow_initialize_fm1dimp_timestep()

!
!MODULES
!

use m_flowgeom, only: ndxi
use unstruc_channel_flow, only: network
use m_CrossSections, only: CalcConveyance
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


end subroutine flow_initialize_fm1dimp_timestep
