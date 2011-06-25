subroutine trisim (numdom, nummap, context_id, fsm_flags, fsm_tracefile, runid)
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
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
!    Function: Main routine for the 2d / 3d program
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use mod_trisim
    use precision
    use dfparall
    !
    ! To raise floating-point invalid, divide-by-zero, and overflow exceptions:
    ! Activate the following line
    ! See also statements below
    !
    ! use ifcore
    !
    ! global data declaration; compare with include 'globdat.igd'
    !
    use globaldata
    implicit none
    type(globDat),pointer :: gdp
    !
    include 'fsm.i'
!
! Parameters
!
    integer       , intent(in)  :: context_id
    integer       , intent(in)  :: fsm_flags
    integer       , intent(in)  :: numdom        ! Number of subdomains (0 means single domain)
                                                 ! as detected by hydra
    integer       , intent(in)  :: nummap        ! Number of mappers (one for each DD boundaries connected with this subdomain)
                                                 ! as detected by hydra
    character(*)  , intent(in)  :: fsm_tracefile
    character(256)              :: runid
!
! Local variables
!
    integer :: ierr
    integer :: retval
    !
    ! To raise floating-point invalid, divide-by-zero, and overflow exceptions:
    ! Activate the following line
    ! See also statements below
    !
    ! INTEGER*4 OLD_FPE_FLAGS, NEW_FPE_FLAGS
!
!! executable statements -------------------------------------------------------
!
    ! To raise floating-point invalid, divide-by-zero, and overflow exceptions:
    ! Activate the following two lines
    ! See also use statement above
    !
    ! NEW_FPE_FLAGS = FPE_M_TRAP_OVF + FPE_M_TRAP_DIV0 + FPE_M_TRAP_INV
    ! OLD_FPE_FLAGS = FOR_SET_FPE (NEW_FPE_FLAGS)
    !
    ! create and initialize GDP structure
    !
    allocate(gdp)
    !
    ! is this necessary?
    !
    nullify(gdp%runid)
    !
    retval = trisim_init(numdom, nummap, context_id, fsm_flags, fsm_tracefile, runid, gdp)
    if (retval /= 0) then
       return
    endif
    !
    retval = trisim_step(gdp)
    if (retval /= 0) then
       return
    endif
    !
    retval = trisim_finish(gdp)
    if (retval /= 0) then
       return
    endif
    !
    !! TODORE call gdp_dealloc(gdp)
    !! TODORE deallocate(gdp, stat=ierr)
    !
    ! Finalizes MPI
    !
    if (parll) then
       call dfexitmpi(0)
    endif
end subroutine trisim
