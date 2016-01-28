subroutine det_num_dis(no_dis, gdp)
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
!    Function: Reads input needed for the coupling of Corjet/Cortime/Cormix
!              with Delft3d-Flow
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
!
    use precision
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
!
! Global variables
!
    integer                :: no_dis
!
! Local variables
!
    integer                :: luntmp
    integer, external      :: newlun
!
!! executable statements -------------------------------------------------------
!
    !
    luntmp = newlun(gdp)
    open (luntmp,file='corinp.dat')

    !
    ! Read base_path were to put the cortiem_linkinp file
    !

    call skipstarlines (luntmp)
    read (luntmp,*) gdp%gdnfl%base_path

    !
    ! Read number of discharges
    !

    call skipstarlines (luntmp)
    read (luntmp,*) no_dis

    !
    ! Close the general cormix input file
    !
    close (luntmp)
    !
end subroutine det_num_dis
