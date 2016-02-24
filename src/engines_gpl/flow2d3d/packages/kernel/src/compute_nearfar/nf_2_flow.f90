subroutine nf_2_flow( filename,x_jet,y_jet,z_jet,s_jet,h_jet,b_jet, no_val)
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
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    integer, intent(inout)                   :: no_val
    real(fp), dimension(no_val), intent(out) :: x_jet, y_jet, z_jet, s_jet, h_jet, b_jet
    character(len=*), intent(in)             :: filename
    
    integer                 :: i, lun, ierr
    
    open( newunit=lun, file = filename )
    
    do i = 1,no_val
        read( lun, *, iostat = ierr ) x_jet(i), y_jet(i), z_jet(i), s_jet(i), h_jet(i), b_jet(i)
        if ( ierr /= 0 ) then
            no_val = i - 1
            exit
        endif
    enddo
    
    close( lun )
end subroutine nf_2_flow