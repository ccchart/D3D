subroutine wrmorm2(lundia    ,error     ,mmax      ,nmaxus    ,lsedtot   , &
                 & nlyr      ,irequest  ,fds       ,grpnam    ,lyrfrac   , &
                 & thlyr     ,gdp  )
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
!    Function: Writes the time varying data for the morphological under layers
!              to the sediment group on the NEFIS FLOW MAP file
!
! Method used:
!
!!--declarations----------------------------------------------------------------
use precision
use bedcomposition_module
    use sp_buffer
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer              , pointer :: celidt
    type (nefiselement)  , pointer :: nefiselem
!
! Global variables
!
    integer                                                                :: fds
    integer                                                                :: irequest
    integer                                                                :: lsedtot
    integer                                                                :: lundia
    integer                                                                :: mmax
    integer                                                                :: nlyr
    integer                                                                :: nmaxus
    logical                                                                :: error
    character(16)                                                          :: grpnam
    real(fp)           , dimension(gdp%d%nmlb:gdp%d%nmub,1:nlyr,1:lsedtot) :: lyrfrac
    real(fp)           , dimension(gdp%d%nmlb:gdp%d%nmub,1:nlyr)           :: thlyr
    type(bedcomp_data)                                                     :: gdmorlyr
!
! Local variables
!
    integer                 :: ierror      ! Local errorflag for NEFIS files
    integer                 :: i
    integer                 :: k
    integer                 :: l
    integer                 :: m
    integer                 :: n
    integer                 :: nm
    integer, external       :: neferr
    integer, external       :: putelt
    integer, dimension(3,5) :: uindex
    character(256)          :: errmsg
!
!! executable statements -------------------------------------------------------
!
    nefiselem => gdp%nefisio%nefiselem(nefiswrsedm)
    celidt  => nefiselem%celidt
    !
    select case (irequest)
    case (1)
       !
       ! Define elements
       !
       call addelm(nefiswrsedm,'LYRFRAC',' ','[   -   ]','REAL',4, &
          & 'Volume fraction of sediment in layer'               , &
          & 4      ,nmaxus ,mmax     ,nlyr     ,lsedtot  ,0      , &
          & lundia ,gdp    )
       call addelm(nefiswrsedm,'THLYR',' ','[   M   ]','REAL',4, &
          & 'Thickness of sediment layer'                      , &
          & 3      ,nmaxus ,mmax     ,nlyr    ,0       ,0      , &
          & lundia ,gdp    )
    case (2)
       !
       ! Write data to file
       !
       uindex (1,1) = celidt
       uindex (2,1) = celidt
       uindex (3,1) = 1 ! increment in time
       !
       ! element 'LYRFRAC'
       !
       call sbuff_checksize(lsedtot*nlyr*mmax*nmaxus)
       i = 0
       do l = 1, lsedtot
          do k = 1, nlyr
             do m = 1, mmax
                do n = 1, nmaxus
                   i        = i+1
                   call n_and_m_to_nm(n, m, nm, gdp)
                   sbuff(i) = real(lyrfrac(nm, k, l),sp)
                enddo
             enddo
          enddo
       enddo
       ierror = putelt(fds, grpnam, 'LYRFRAC', uindex, 1, sbuff)
       if (ierror/= 0) goto 9999
       !
       ! element 'THLYR'
       !
       i = 0
       do k = 1, nlyr
          do m = 1, mmax
             do n = 1, nmaxus
                i        = i+1
                call n_and_m_to_nm(n, m, nm, gdp)
                sbuff(i) = real(thlyr(nm, k),sp)
             enddo
          enddo
       enddo
       ierror = putelt(fds, grpnam, 'THLYR', uindex, 1, sbuff)
       if (ierror/= 0) goto 9999
       !
       ! write error message if error occurred and set error = .true.
       ! the files will be closed in clsnef (called in triend)
       !
 9999  continue
       if (ierror/= 0) then
          ierror = neferr(0, errmsg)
          call prterr(lundia, 'P004', errmsg)
          error = .true.
       endif
    endselect
end subroutine wrmorm2
