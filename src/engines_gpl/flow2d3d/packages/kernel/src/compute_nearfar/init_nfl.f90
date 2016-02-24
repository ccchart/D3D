subroutine init_nfl  (kmax, lstsci, gdp   )
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
!    Function: Converts flow results to cormix input
!              nog wel doen:
!              1) interpoleren naar equidistante laagverdeling jet3d
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
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !

    integer                         , pointer :: idensform
    character(256)                  , pointer :: nflmod
    integer                         , pointer :: no_dis

    integer                                   :: kmax
    integer                                   :: lstsci
    integer                                   :: istat
!
!! executable statements -------------------------------------------------------
!

    idensform      => gdp%gdphysco%idensform
    nflmod         => gdp%gdnfl%nflmod
    no_dis         => gdp%gdnfl%no_dis

    !
    ! Initialisation
    !

    select case (nflmod)
       case ('corjet','cortime','generic')
          call det_num_dis(no_dis, gdp) !FIXME if generic actually is using xml input
       case ('jet3d')
          no_dis = 1
    end select

                  allocate (gdp%gdnfl%m_diff    (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%n_diff    (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%m_amb     (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%n_amb     (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%m_intake  (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%n_intake  (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%k_intake  (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%q_diff    (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%t0_diff   (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%s0_diff   (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%rho0_diff (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%d0        (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%h0        (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%sigma0    (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%theta0    (no_dis)                                    , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%basecase  (no_dis,2)                                  , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%disnf     (gdp%d%nmlb:gdp%d%nmub, kmax,no_dis)        , stat = istat)
    if(istat==0)  allocate (gdp%gdnfl%sournf    (gdp%d%nmlb:gdp%d%nmub, kmax,lstsci,no_dis) , stat = istat)

    !
    ! Read corjet/cortime input data from file corinp.dat
    !

!   select case (nflmod)
!      case ('corjet','cortime')
!         call corinp_gen(idensform,gdp)
!   end select

    !
gdp%gdnfl%m_diff   = 0.0_fp
gdp%gdnfl%n_diff   = 0.0_fp
gdp%gdnfl%m_amb    = 0.0_fp
gdp%gdnfl%n_amb    = 0.0_fp
gdp%gdnfl%m_intake = 0.0_fp
gdp%gdnfl%n_intake = 0.0_fp
gdp%gdnfl%k_intake = 0.0_fp
gdp%gdnfl%q_diff   = 0.0_fp
gdp%gdnfl%t0_diff  = 0.0_fp
gdp%gdnfl%s0_diff  = 0.0_fp
gdp%gdnfl%rho0_diff= 0.0_fp
gdp%gdnfl%d0       = 0.0_fp
gdp%gdnfl%h0       = 0.0_fp
gdp%gdnfl%sigma0   = 0.0_fp
gdp%gdnfl%theta0   = 0.0_fp
gdp%gdnfl%basecase = ' '
gdp%gdnfl%disnf    = 0.0_fp
gdp%gdnfl%sournf   = 0.0_fp

end subroutine init_nfl
