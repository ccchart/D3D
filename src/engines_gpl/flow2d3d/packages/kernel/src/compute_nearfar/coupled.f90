subroutine coupled (add     , r0, kmax, lstsci, lcon  , thick , m_intake, n_intake, &
                  & k_intake, s0, dps , dzs0  , kfsmn0, kfsmx0, zmodel  , gdp     )
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
!    Function: Coupled intakes/outlets for cortime to flow coupling
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
    integer                                                    , intent(in)  :: kmax
    integer                                                    , intent(in)  :: lstsci
    integer                                                    , intent(in)  :: lcon
    integer                                                    , intent(in)  :: m_intake
    integer                                                    , intent(in)  :: n_intake
    integer                                                    , intent(in)  :: k_intake
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmx0     ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmn0     ! Description and declaration in esm_alloc_int.f90
    real(fp)                                                   , intent(out) :: add
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: s0
    real(prec) , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: dps
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: dzs0       ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax,lstsci) , intent(in)  :: r0         ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(kmax)                               , intent(in)  :: thick      ! Description and declaration in esm_alloc_real.f90
    logical                                                    , intent(in)  :: zmodel
!
! Local variables
!
    integer                                                                  :: nm_intake
    integer                                                                  :: k
!
!! executable statements -------------------------------------------------------
!
    add = 0.0_fp

    if (m_intake > 0) then
       call n_and_m_to_nm(n_intake, m_intake, nm_intake, gdp)
       if (k_intake == 0) then
          if (.not. zmodel) then
             do k = 1, kmax
                add = add + thick(k)*r0(nm_intake,k,lcon)
             enddo
          else
             do k = kfsmn0(nm_intake), kfsmx0(nm_intake)
                add = add + dzs0(nm_intake,k)*r0(nm_intake,k,lcon)/(s0(nm_intake)+real(dps(nm_intake),fp))
             enddo
          endif
       else
          add = r0(nm_intake,k_intake,lcon)
       endif
    endif
end subroutine coupled
