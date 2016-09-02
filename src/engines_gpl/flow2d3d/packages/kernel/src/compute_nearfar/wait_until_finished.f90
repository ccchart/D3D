subroutine wait_until_finished (no_dis, waitfiles, idis, filename, waitlog, gdp)
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
    use flow2d3d_timers
    !
    implicit none
    !
    type(globdat),target :: gdp
!
! Global variables
!
    integer                        , intent(in)  :: no_dis
    character(*), dimension(no_dis)              :: waitfiles
    integer                        , intent(out) :: idis
    character(*)                   , intent(out) :: filename
    logical                        , intent(in)  :: waitlog
!
! Local variables
!
    integer , external      :: newlun
    integer                 :: lun
    integer                 :: ios
    integer                 :: i
    integer                 :: numlines
    integer                 :: sleeptime
    logical                 :: ex_file
    logical                 :: opend
!
!! executable statements -------------------------------------------------------
!
    ! Check for how many files we are waiting to appear
    idis = 0
    do i=1, no_dis
       if (waitfiles(i) /= ' ') then
          idis = i
          if (waitlog) then
             write(*,'(3a)') "Waiting for file '", trim(waitfiles(i)), "' to appear ..."
          endif
       endif
    enddo
    !
    ! Return when all files did appear (and waitfiles is empty). idis must be 0.
    if (idis == 0) return
    !
    ex_file = .false.
    !
    ! Examine if one of the files exists
    ! This will cost CPU time, but there is nothing else to do (Cosumo/Cormix run on another machine)
    !
    call timer_start(timer_wait, gdp)
    idis     = 0
    filename = ' '
    do while (.not. ex_file)
       do i=1, no_dis
          if (waitfiles(i) /= ' ') then
             inquire (file=waitfiles(i), exist=ex_file)
             if (ex_file) then
                filename     = waitfiles(i)
                idis         = i
                waitfiles(i) =  ' '
                exit
             endif
          endif
       enddo
    enddo
    call timer_stop(timer_wait, gdp)
    write(*,'(3a)') "Scanning    file '", trim(filename), "' (after 1 second) ..."
    !
    ! sleep time in milliseconds
    ! Ad Hoc solution, needed for Cosumo to finish writing the file
    !
    sleeptime = 1000
    call CUTIL_SLEEP(sleeptime);
    !
    ! File found: open file and read until you find eof
    !
    lun = newlun(gdp)
    inquire(lun, iostat=ios, opened=opend)
 10 if (ios==0 .and. opend) close(lun, iostat=ios)
    open (lun,file=filename,err=10)
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
end subroutine wait_until_finished
