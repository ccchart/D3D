subroutine wait_until_finished (filename,gdp)
!----- GPL ---------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2011-2016.
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
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
!    Function: Compute the jet trajectory in "world" coordinates from cortim results
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
!
! Global variables
!
    character*256, intent(in)    :: filename
!
! Local variables
!
    integer , external      :: newlun
    integer                 :: lun
    integer                 :: ios
    integer                 :: numlines
    logical                 :: ex_file
    logical                 :: opend
!
!! executable statements -------------------------------------------------------
!
    ex_file = .false.
    !
    ! Examine if the file exists
    ! This will cost CPU time, but there is nothing else to do (Cosumo/Cormix run on another machine)
    !
    write(*,'(3a)') "Waiting for file '", trim(filename), "' to appear ..."
    do while (.not. ex_file)
        inquire (file=filename,exist=ex_file)
    enddo
    write(*,'(a)') "Scanning    file ..."
    !
    ! File found: open file and read until you find eof
    !
    lun = newlun(gdp)
    inquire(lun, iostat=ios, opened=opend)
 10 if (ios==0 .and. opend) close(lun, iostat=ios)
    open (lun,file=filename)
    rewind (lun)
    ios      = 0
    numlines = 0
    do while (ios == 0)
       read (lun,'(a1)',iostat=ios,err=10,end=20)
       numlines = numlines + 1
    enddo
    write(*,'(a)') "ERROR: This line should not be reached."
    goto 10
20  close(lun)
    write(*,'(a,i0)') "numlines: ", numlines
end subroutine wait_until_finished
