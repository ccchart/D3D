subroutine findnmk(xz     ,yz     ,dps    ,s1    ,kcs    ,nmmax  , &
                 & thick  ,kmax   ,x_jet  ,y_jet ,z_jet  ,nm_jet , &
                 & k_jet  ,kfsmn0 ,kfsmx0 ,dzs0  ,zmodel ,gdp    )
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
!    Function: Finds n,m and k coordinates of "jet" points
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
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
    integer                                             , intent(in)  :: kmax   ! Description and declaration in tricom.igs
    integer                                             , intent(in)  :: nmmax  ! Description and declaration in tricom.igs
    integer                                             , intent(out) :: nm_jet
    integer                                             , intent(out) :: k_jet
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: kcs    ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: kfsmn0 ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: kfsmx0 ! Description and declaration in esm_alloc_int.f90
    real(fp)                                            , intent(in)  :: x_jet
    real(fp)                                            , intent(in)  :: y_jet
    real(fp)                                            , intent(in)  :: z_jet
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: s1     ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: xz     ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: yz     ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(kmax)                        , intent(in)  :: thick  ! Description and declaration in esm_alloc_real.f90
    real(prec) , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: dps    ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax) , intent(in)  :: dzs0   ! Description and declaration in esm_alloc_real.f90
    logical                                             , intent(in)  :: zmodel
!
! Local variables
!
    integer       :: k
    integer       :: nm
    real(fp)      :: dist
    real(fp)      :: distmin
    real(fp)      :: r_onder
    real(fp)      :: r_boven
!
!! executable statements -------------------------------------------------------
!
    !
    ! Find the horizontal nm coordinate of the "Jet" trajectory point
    !
    nm_jet  = 0
    distmin = 1.0e+30_fp
    do nm = 1, nmmax
       if (kcs(nm) == 1) then
          dist = sqrt((xz(nm) - x_jet)**2 + (yz(nm) - y_jet)**2)
          if (dist < distmin) then
             nm_jet  = nm
             distmin = dist
          endif
       endif
    enddo
    !
    ! Find the vertical position
    !
    k_jet   = 0
    !
    ! Is this correct?? Seems as if z_jet is relative to the bed
    !
    r_boven = real(dps(nm_jet),fp) + s1(nm_jet)
    if (.not. zmodel) then
       if (z_jet >= r_boven) then
          k_jet = 1
       else
          do k = 1, kmax - 1
             r_onder = r_boven - thick(k)*(real(dps(nm_jet),fp) + s1(nm_jet))
             if (z_jet <= r_boven .and. z_jet >= r_onder) then
                k_jet = k
                exit
             endif
             r_boven = r_onder
          enddo
          if (k_jet == 0) k_jet = kmax
       endif
    else
       if (z_jet >= r_boven) then
          k_jet = kfsmx0(nm_jet)
       else
          do k = kfsmx0(nm_jet), kfsmn0(nm_jet) + 1,-1
             r_onder = r_boven - dzs0(nm_jet,k)
             if (z_jet < r_boven .and. z_jet >= r_onder) then
                k_jet = k
                exit
             endif
             r_boven = r_onder
          enddo
          if (k_jet == 0) k_jet = kfsmn0(nm_jet)
       endif
    endif
    !
end subroutine findnmk
