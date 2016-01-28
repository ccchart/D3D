subroutine init_nfl  (kmax, lstsci, gdp   )

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
       case ('corjet','cortime')
          call det_num_dis(no_dis, gdp)
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
end subroutine init_nfl
