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

module m_f1dimp
   use m_f1dimp_data
   !
   ! flow 1d implicit
   !
   logical                           :: f1dimp_initialized=.false.   
   type(f1dimppar_type), target      :: f1dimppar         !< flow 1d implicit parameters
   !
   
    contains
!
!BEGIN reallocate_fill
!
    
subroutine reallocate_fill(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

double precision, dimension(:), allocatable, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call realloc(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill
    
!
!END reallocate_fill_int
!

!
!BEGIN reallocate_fill_int
!

!FM1DIMP2DO: There must be a better way to do it. This is just the same but for integers
    
subroutine reallocate_fill_int(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

integer, dimension(:), allocatable, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call realloc(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill_int
    
!
!END reallocate_fill_int
!

!
!BEGIN reallocate_fill_pointer
!

!FM1DIMP2DO: There must be a better way to do it. This is just the same but for integers
    
subroutine reallocate_fill_pointer(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

real(fp), dimension(:), pointer, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call reallocp(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill_pointer
    
!
!END reallocate_fill_pointer
!

!
!BEGIN reallocate_fill_int_pointer
!
    
subroutine reallocate_fill_int_pointer(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

integer, dimension(:), pointer, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call reallocp(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill_int_pointer
    
!
!END reallocate_fill_int_pointer
!

end module m_f1dimp
