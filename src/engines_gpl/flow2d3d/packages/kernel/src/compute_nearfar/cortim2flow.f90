subroutine cortim2flow(thick  ,kmax  ,dps   ,s0    ,r0       ,          &
                     & lstsci ,lsal  ,ltem  ,xz    ,yz       ,nmmax   , &
                     & kcs    ,filename     ,taua            ,idis    , &
                     & linkinf,gdp     )
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
!    Function: Converts cormix output to delft3d sources
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
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    integer ,dimension(:)          , pointer :: m_diff
    integer ,dimension(:)          , pointer :: n_diff
!
! Global variables
!
    integer                                                      , intent(in)  :: idis
    integer                                                      , intent(in)  :: kmax     !  Description and declaration in tricom.igs
    integer                                                      , intent(in)  :: lstsci   !  Description and declaration in tricom.igs
    integer                                                      , intent(in)  :: lsal     !  Description and declaration in tricom.igs
    integer                                                      , intent(in)  :: ltem     !  Description and declaration in tricom.igs
    integer                                                      , intent(in)  :: nmmax    !  Description and declaration in tricom.igs
    integer      , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kcs      !  Description and declaration in
    real(fp)                                                     , intent(in)  :: taua
    real(fp)     , dimension(8)                                  , intent(in)  :: linkinf
    real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: s0       !  Description and declaration in rjdim.f90 gs
    real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: xz       !  Description and declaration in rjdim.f90 gs
    real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: yz       !  Description and declaration in rjdim.f90
    real(fp)     , dimension(gdp%d%nmlb:gdp%d%nmub,kmax,lstsci)  , intent(in)  :: r0       !  Description and declaration in rjdim.f90
    real(fp)     , dimension(kmax)                               , intent(in)  :: thick    !  Description and declaration in rjdim.f90 gs
    real(prec)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: dps      !  Description and declaration in rjdim.f90
    character*256,dimension(3)                                   , intent(in)  :: filename
!
! Local variables
!
    integer                                     :: nm_diff
    integer                      , external     :: newlun
    integer                                     :: no_modules
    integer                                     :: luntmp
    integer                                     :: iocond
    integer                                     :: irow
    integer                                     :: itel
    integer                                     :: nrow
    integer                                     :: no_val
    integer       , dimension(:) , allocatable  :: modstart
    real(fp)                                    :: xx
    real(fp)                                    :: yy
    real(fp)                                    :: zz
    real(fp)                                    :: ss
    real(fp)                                    :: xstart
    real(fp)                                    :: xend
    real(fp)                                    :: ystart
    real(fp)                                    :: yend
    real(fp)      , dimension(:) , allocatable  :: x_cor
    real(fp)      , dimension(:) , allocatable  :: y_cor
    real(fp)      , dimension(:) , allocatable  :: z_cor
    real(fp)      , dimension(:) , allocatable  :: s_cor
    real(fp)      , dimension(:) , allocatable  :: bv_cor
    real(fp)      , dimension(:) , allocatable  :: bh_cor
    real(fp)      , dimension(:) , allocatable  :: v_cor
    real(fp)      , dimension(:) , allocatable  :: x_jet
    real(fp)      , dimension(:) , allocatable  :: y_jet
    real(fp)      , dimension(:) , allocatable  :: z_jet
    real(fp)      , dimension(:) , allocatable  :: s_jet
    real(fp)      , dimension(:) , allocatable  :: bv_jet
    real(fp)      , dimension(:) , allocatable  :: bh_jet
    real(fp)      , dimension(:) , allocatable  :: v_jet

    character*256 , dimension(:) , allocatable  :: modules
!
!! executable statements -------------------------------------------------------
!
    m_diff         => gdp%gdnfl%m_diff
    n_diff         => gdp%gdnfl%n_diff

    call n_and_m_to_nm(n_diff(idis), m_diff(idis), nm_diff, gdp)

    !
    ! Open cormix output file
    !

    luntmp = newlun (gdp)
    open (luntmp,file=filename(2),status='old')

    !
    ! Determine the number of modules used by Cormix and the number of jut/plume trajectory values
    ! Allocate temporary arrays
    !

    call cortim_no_modules (luntmp,no_modules,nrow)
    rewind luntmp

    allocate (x_cor   (nrow))
    allocate (y_cor   (nrow))
    allocate (z_cor   (nrow))
    allocate (s_cor   (nrow))
    allocate (bv_cor  (nrow))
    allocate (bh_cor  (nrow))
    allocate (v_cor   (nrow))
    allocate (x_jet   (nrow))
    allocate (y_jet   (nrow))
    allocate (z_jet   (nrow))
    allocate (s_jet   (nrow))
    allocate (bv_jet  (nrow))
    allocate (bh_jet  (nrow))
    allocate (v_jet   (nrow))
    allocate (modules (no_modules))
    allocate (modstart(no_modules))

    !
    ! Read cormix jet/plume trajectory and belonging characteristics
    ! Allocate temporary arrays
    !

    call cortim_xyzs (luntmp,no_modules, nrow, modules,modstart,x_cor,y_cor,z_cor,s_cor,bv_cor,bh_cor,v_cor)

    close (luntmp)

    !
    ! Cortim results to jet trajectory (in real world coordinates)
    !

    call ct2jettraj (no_modules, nrow, modules, modstart, x_cor  , y_cor  , z_cor , s_cor  , bv_cor , bh_cor , &
                                                          x_jet  , y_jet  , z_jet , s_jet  , bv_jet , bh_jet , &
                                                          xz     , yz     , dps   , s0     , nm_diff,&
                                                          no_val , taua           , idis   , gdp    )

    !
    ! Fill sources and sinks following the Desa Method of Prof. Lee
    !



    call desa(x_jet   ,y_jet    ,z_jet   ,s_jet   ,no_val  , &
            & kcs     ,xz       ,yz      ,dps     ,s0      , &
            & nmmax   ,thick    ,kmax    ,lstsci  ,lsal    , &
            & ltem    ,bv_jet  ,bh_jet   ,v_jet   ,idis    , &
            & xstart  ,xend    ,ystart   ,yend    ,r0      , &
            & linkinf ,gdp     )

    !
    ! Temporarily, write cormix trajectory to tekal file for postprocessing
    !

    call wri_tek (x_jet  , y_jet  , z_jet , no_val  , xstart  , xend  , ystart  , yend  , filename(3), &
                & bv_jet , bh_jet , s_jet , linkinf , gdp     )


    !
    ! Deallocate temporary arrays
    !

    deallocate (x_cor)
    deallocate (y_cor)
    deallocate (z_cor)
    deallocate (s_cor)
    deallocate (bv_cor)
    deallocate (bh_cor)
    deallocate (v_cor)
    deallocate (x_jet)
    deallocate (y_jet)
    deallocate (z_jet)
    deallocate (s_jet)
    deallocate (bv_jet)
    deallocate (bh_jet)
    deallocate (v_jet)
    deallocate (modules)
    deallocate (modstart)
    !
end subroutine cortim2flow
