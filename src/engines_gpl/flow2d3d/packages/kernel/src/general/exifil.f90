function exifil(file      ,lundia    ,errornr    ,gdp       )
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
!    Function: A logical function which checks the existence of a
!              specified file (path may be included). Set to TRUE
!              when the file is found, FALSE otherwise.
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
!
! Global variables
!
    integer         :: lundia   !  Description and declaration in inout.igs
    logical         :: exifil   !!  Logical function which tests the existence of a specified file in a
                                !!  the assigned directory (usually the work directory)
    character(*)    :: file     !!  Help var.
    character(4)    :: errornr  !!  error number (stored as characters) in case error occurred
!
!
! Local variables
!
    integer    :: ipos   ! Help var. 
    logical    :: ex     ! Help flag = TRUE when file is found 
!
!
!! executable statements -------------------------------------------------------
!
    !
    !
    !
    call noextspaces(file      ,ipos      )
    !
    inquire (file = file(:ipos), exist = ex)
    if (.not.ex) then
       call prterr(lundia    ,errornr    ,file(:ipos)          ,gdp       )
       !
       exifil = .false.
    else
       exifil = .true.
    endif
end function exifil
