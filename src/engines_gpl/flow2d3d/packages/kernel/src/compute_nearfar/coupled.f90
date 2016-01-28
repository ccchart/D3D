subroutine coupled (add,r0,kmax,lstsci,lcon,thick,m_intake,n_intake,k_intake,gdp)

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
    real(fp)                                                   , intent(out) :: add
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax,lstsci) , intent(in)  :: r0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(kmax)                               , intent(in)  :: thick    !  Description and declaration in esm_alloc_real.f90
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
          do k = 1, kmax
             add = add + thick(k)*r0(nm_intake,k,lcon)
          enddo
       else
          add = r0(nm_intake,k_intake,lcon)
       endif
    endif
end subroutine coupled
