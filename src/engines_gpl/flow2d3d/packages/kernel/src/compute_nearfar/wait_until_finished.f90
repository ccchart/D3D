subroutine wait_until_finished (no_dis, waitfiles, idis, filename, waitlog, error, gdp)
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
    logical                        , intent(out) :: error
!
! Local variables
!
    integer , external      :: newlun
    integer                 :: lun
    integer                 :: ios
    integer                 :: i
    integer                 :: sleeptime
    logical                 :: ex_file    ! true: file exists
    logical                 :: fileok     ! true: file contains the end tag </NF2FF>
    logical                 :: opend
    character(300)          :: line
!
!! executable statements -------------------------------------------------------
!
    error = .false.
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
    fileok = .false.
    do while (.not.fileok)
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
                   !
                   ! Do not remove the found file from the waitfiles here:
                   ! It will be removed in near_field, when the full reading of the file finished successfully
                   !
                   exit
                endif
             endif
          enddo
       enddo
       call timer_stop(timer_wait, gdp)
       !
       ! File found: open file and search for the end tag </NF2FF>
       !
       lun = newlun(gdp)
       inquire(lun, iostat=ios, opened=opend)
       if (ios /= 0) then
          ! try again
          cycle
       endif
       if (opend) close(lun, iostat=ios)
       open (lun,file=filename,iostat=ios)
       if (ios /= 0) then
          ! try again
          cycle
       endif
       ios    = 0
       do while (ios == 0)
          read (lun,'(a)',iostat=ios) line
          if (index(line,'</NF2FF>') >= 1) then
             fileok = .true.
          endif
       enddo
       close(lun)
    enddo
    write(*,'(3a)') "Scanned    file '", trim(filename), "'"
end subroutine wait_until_finished
