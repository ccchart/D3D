subroutine det_num_dis(no_dis, gdp)
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
    read (luntmp,'(a)') gdp%gdnfl%base_path

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
