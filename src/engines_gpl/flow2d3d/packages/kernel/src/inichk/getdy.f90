function getdy(sferic,x1,y1,x2,y2,gdp)
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!!--description-----------------------------------------------------------------
!
!    Function: Calculates "north-south" distance between two points on earth
! Method used: Circular distance when sferic is true,
!              Euclidic distance when sferic is false
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    real(hp) , pointer :: ddegrad
    real(hp) , pointer :: dearthrad
!
! Global variables
!
    real(fp) :: x1
    real(fp) :: y1
    real(fp) :: x2
    real(fp) :: y2
    real(fp) :: getdy
    logical  :: sferic
!
! Local variables
!
    real(hp) :: yy1
    real(hp) :: yy2
!
!! executable statements -------------------------------------------------------
!
    ddegrad    => gdp%gdconstd%ddegrad
    dearthrad  => gdp%gdconstd%dearthrad
    !
    if (sferic) then
       yy1   = real(y1,hp)*ddegrad
       yy2   = real(y2,hp)*ddegrad
       getdy = real(dearthrad*(yy2-yy1),fp)
    else
       getdy = y2-y1
    endif
end function getdy
