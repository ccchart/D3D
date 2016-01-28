subroutine wait_a_sec (no_seconds)

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
use precision
!
! Global variables
!
    integer                                                     , intent(in)    :: no_seconds
!
! Local variables
!
    real(sp)                                                                        :: timnow
    real(sp)                                                                        :: timstart

!
!! executable statements -------------------------------------------------------
!
    call  cpu_time(timstart)
    call  cpu_time(timnow)

    do while (timnow - timstart < no_seconds)
       call  cpu_time(timnow)
    enddo
    
end subroutine wait_a_sec
