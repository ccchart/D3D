subroutine wait_until_finished (filename,gdp)

!!--copyright-------------------------------------------------------------------
! Copyright (c) 2009, WL | Delft Hydraulics. All rights reserved.
!!--disclaimer------------------------------------------------------------------
! This code is part of the Delft3D software system. WL|Delft Hydraulics has
! developed c.q. manufactured this code to its best ability and according to the
! state of the art. Nevertheless, there is no express or implied warranty as to
! this software whether tangible or intangible. In particular, there is no
! express or implied warranty as to the fitness for a particular purpose of this
! software, whether tangible or intangible. The intellectual property rights
! related to this software code remain with WL|Delft Hydraulics at all times.
! For details on the licensing agreement, we refer to the Delft3D software
! license and any modifications to this license, if applicable. These documents
! are available upon request.
!!--version information---------------------------------------------------------
! $Author$
! $Date$
! $Revision$
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
    character*256                                                   , intent(in)    :: filename
!
! Local variables
!
    integer                                                         , external      :: newlun
    integer                                                                         :: lun
    integer                                                                         :: ios  

    logical                                                                         :: ex_file

!
!! executable statements -------------------------------------------------------
!
    lun     = newlun(gdp)
    ex_file = .false.

!   examine if the file exists

    do while (.not. ex_file)
        inquire (file=filename,exist = ex_file)
    end do

!   open file and read until you find eof

    open (lun,file=filename)

 10 rewind (lun)

    ios = 0
    do while (ios == 0)
       read (lun,'(a1)',iostat=ios,err=10,end=20)
    end do
    goto 10

 20 close(lun)
end subroutine wait_until_finished
