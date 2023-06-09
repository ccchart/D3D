!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2023.                                
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

! 
! 

!> Highlights the 'string' field of a user-input form field.
!! Input fields are highlighted automatically, but the string label isn't.
!! This assumes that string field number is always input field number minus 1.
!! Only use this subroutine as the FMUSER argument to IFormEditUser(.., .., FMUSER).
subroutine highlight_form_line(ifield, iexitk)
implicit none
integer, intent(in) :: ifield !< Form field number that lost focus (infoform(3) contains 'next' field).
integer, intent(in) :: iexitk !< 'Exit' key that was used to leave this form field.

integer :: ifieldnext
integer, external :: InfoForm

ifieldnext = InfoForm(3)

! Reset the 'current' field back to defaults (no highlights)
if (ifield > 1) then
    call iformattributen(ifield-1, 0, -1, -1)
    call iformshowfield(ifield-1)
end if
if (ifieldnext > 1) then
    call iformattribute(ifieldnext-1, 'UB', ' ', ' ')
    call iformshowfield(ifieldnext-1)
end if

end subroutine highlight_form_line
