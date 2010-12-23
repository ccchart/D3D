subroutine z_adi(stage     ,j         ,nmmaxj    ,nmmax     ,kmax      , &
               & mmax      ,nmax      ,zmodel    ,nst       ,nfltyp    , &
               & nocol     ,norow     ,nsrc      ,dismmt    ,irocol    , &
               & mnksrc    ,kfu       ,kfv       ,kfs       ,            &
               & kspu      ,kspv      ,kadu      ,kadv      ,kcs       , &
               & kcu       ,kcv       ,kfsmin    ,kfsmax    ,kfsmx0    , &
               & kfumin    ,kfumax    ,kfumx0    ,kfvmin    ,kfvmax    , &
               & kfvmx0    ,kfuz1     ,kfvz1     ,kfsz1     ,kcu45     , &
               & kcv45     ,kcscut    ,porosu    ,porosv    ,areau     , &
               & areav     ,volum1    ,s0        ,s1        ,w1        , &
               & u0        ,u1        ,v0        ,v1        ,hu        , &
               & hv        ,thick     ,umean     ,ubrlsu    ,ubrlsv    , &
               & vmean     ,dpu       ,dpv       ,dps       ,dzu0      , &
               & dzv0      ,dzs0      ,dzu1      ,dzv1      ,dzs1      , &
               & qxk       ,qyk       ,qzk       ,evap      ,circ2d    , &
               & circ3d    ,drhodx    ,drhody    ,disch     ,umdis     , &
               & vmdis     ,wsu       ,wsv       ,gud       ,gvd       , &
               & guu       ,guv       ,gvv       ,gvu       ,guz       , &
               & gvz       ,gsqs      ,gsqiu     ,gsqiv     ,taubpu    , &
               & taubpv    ,taubsu    ,taubsv    ,vicuv     ,vnu2d     , &
               & vicww     ,rxx       ,rxy       ,ryy       ,windu     , &
               & windv     ,patm      ,fcorio    ,tgfsep    ,wrka1     , &
               & wrka2     ,wrka3     ,wrka4     ,wrka5     ,wrka6     , &
               & wrka7     ,wrka8     ,wrka15    ,wrka16    ,wrkb1     , &
               & wrkb2     ,wrkb3     ,wrkb4     ,wrkb5     ,wrkb6     , &
               & wrkb7     ,wrkb8     ,wrkb9     ,wrkb10    ,zk        , &
               & p0        ,crbc      ,hu0       ,hv0       ,hu2       , &
               & dzu2      ,wrkb13    ,wrkb14    ,pship     ,diapl     , &
               & rnpl      ,sbkol     ,cfurou    ,cfvrou    ,gdp       )
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
!    Function: ADI performs one time step of the Alternating
!              Direction Implicit (ADI) method
! Method used: A.D.I. method is used.
!              Upwind-approach for wet cross section in shallow
!              areas or if the model area contains structures.
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use timers
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
    include 'flow_steps_f.inc'
!
! Global variables
!
    integer                                               :: j      !!  Begin pointer for arrays which have
                                                                    !!  been transformed into 1D arrays.
                                                                    !!  Due to the shift in the 2nd (M-)
                                                                    !!  index, J = -2*NMAX + 1
    integer                                               :: kmax   !  Description and declaration in iidim.f90
    integer                                               :: mmax   !  Description and declaration in iidim.f90
    integer                                               :: nfltyp !  Description and declaration in iidim.f90
    integer                                               :: nmax   !  Description and declaration in iidim.f90
    integer                                               :: nmmax  !  Description and declaration in dimens.igs
    integer                                               :: nmmaxj !  Description and declaration in dimens.igs
    integer                                               :: nocol  !  Description and declaration in iidim.f90
    integer                                               :: norow  !  Description and declaration in iidim.f90
    integer                                               :: nsrc   !  Description and declaration in iidim.f90
    integer                                               :: nst    !!  Time step number
    integer, dimension(5, norow + nocol)                  :: irocol !  Description and declaration in iidim.f90
    integer, dimension(7, nsrc)                           :: mnksrc !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kcs    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kcu    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kcv    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfs    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfsmax !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfsmin !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfsmx0 !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfu    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfumax !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfumin !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfumx0 !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfv    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfvmax !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfvmin !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)             :: kfvmx0 !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)     :: kspu   !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)     :: kspv   !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kadu   !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kadv   !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kcscut !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kcu45  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kcv45  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kfsz1  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kfuz1  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)       :: kfvz1  !  Description and declaration in iidim.f90
    logical                                               :: sbkol  !  Description and declaration in procs.igs
    logical                                               :: zmodel !  Description and declaration in procs.igs
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, 3)         :: cfurou !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, 3)         :: cfvrou !  Description and declaration in rjdim.f90
    real(fp), dimension(12, norow + nocol)                :: crbc   !  Description and declaration in rjdim.f90
    real(fp), dimension(4, norow + nocol)                 :: circ2d !  Description and declaration in rjdim.f90
    real(prec), dimension(gdp%d%nmlb:gdp%d%nmub)          :: dps    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: dpu    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: dpv    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: evap   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: fcorio !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gsqiu  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gsqiv  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gsqs   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gud    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: guu    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: guv    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: guz    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gvd    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gvu    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gvv    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: gvz    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: hu     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: hv     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: hu0    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: hu2    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: hv0    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: patm   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: pship  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: s0     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: s1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: taubpu !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: taubpv !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: taubsu !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: taubsv !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: tgfsep !!  Water elev. induced by tide gen.force
                                                                    !!  Internal work array WRKB17 used
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: umean  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: vmean  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: vnu2d  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: windu  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: windv  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka1  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka15 !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka16 !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka2  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka3  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka4  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka5  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka6  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka7  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wrka8  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wsu    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)            :: wsv    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)    :: qzk    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)    :: vicww  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)    :: w1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax+2)    :: vicuv  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: areau  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: areav  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: diapl  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzs0   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzs1   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzu0   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzu1   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzu2   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzv0   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: dzv1   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: p0     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: porosu !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: porosv !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: qxk    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: qyk    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: rnpl   !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: rxx    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: rxy    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: ryy    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: drhodx !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: drhody !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: u0     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: u1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: ubrlsu !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: ubrlsv !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: v0     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: v1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: volum1 !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb1  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb2  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb3  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb4  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb5  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb6  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb7  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb8  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb9  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb10 !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb13 !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)      :: wrkb14 !  Description and declaration in rjdim.f90
    real(fp), dimension(kmax)                             :: thick  !  Description and declaration in rjdim.f90
    real(fp), dimension(0:kmax)                           :: zk
    real(fp), dimension(kmax, 2, norow + nocol)           :: circ3d !  Description and declaration in rjdim.f90
    real(fp), dimension(nsrc)                             :: disch  !  Description and declaration in rjdim.f90
    real(fp), dimension(nsrc)                             :: umdis  !  Description and declaration in rjdim.f90
    real(fp), dimension(nsrc)                             :: vmdis  !  Description and declaration in rjdim.f90
    character(1), dimension(nsrc)                         :: dismmt !  Description and declaration in ckdim.f90
    character(8)                            , intent(in)  :: stage  !!  First or Second half time step
!
! Local variables
!
    integer :: icx
    integer :: icy
    integer :: idry
    integer :: nhystp
    integer :: nmaxddb
    logical :: flood   ! Flag for activating flooding part of checku subroutine 
!
!! executable statements -------------------------------------------------------
!
    nmaxddb = nmax + 2*gdp%d%ddbound
    !
    !
    ! =====================================
    ! COMPUTATION OF STAGE 1 FOR ADI METHOD
    ! =====================================
    !
    if (stage=='stage1') then
       !
       ! Computation of V1, i.e. evaluate momentum equation for one half timest
       !     calculate hv and set kfv = 0 for hv < htrsh (.5*dryflc)
       !
       !
       ! checku before uzd not needed
       !
       ! dummy z_uzd call (for synchronisation with domain decomposition
       ! implementation of the sigma-model)
       !
       call timer_start(timer_uzd, gdp)
       !
       gdp%dd%uzditer = 0
       icx = 1
       icy = nmaxddb
       call timer_start(timer_1stuzd, gdp)
       call z_uzd(j         ,nmmaxj    ,nmmax     ,kmax      ,icx       , &
                & icy       ,nsrc      ,kcs       ,kcv45     ,kcscut    , &
                & kcv       ,kfv       ,kfvz1     ,kfvmin    ,kfvmax    , &
                & kfu       ,kfuz1     ,kfs       ,kfsz1     ,kfsmin    , kfsmax    , &
                & v0        ,u1        ,w1        ,hv        ,dzv1      , &
                & gvv       ,guu       ,guv       ,gvu       ,gsqs      , &
                & gvd       ,gud       ,gvz       ,guz       ,gsqiv     , &
                & disch     ,umdis     ,kspv      ,mnksrc    ,dismmt    , &
                & wrkb1     ,wrkb2     ,wrkb3     ,wrkb4     ,wrkb5     , &
                & wrkb6     ,wrkb7     ,wrkb8     , &
                & vicuv     ,vnu2d     ,vicww     ,tgfsep    , &
                & drhody    ,wsv       ,taubpv    ,taubsv    ,ryy       , &
                & rxy       ,windv     ,patm      ,fcorio    ,p0        , &
                & ubrlsv    ,pship     ,diapl     ,rnpl      ,nst       , &
                & v1        ,s0        ,cfvrou    ,wrkb13    ,wrkb14    , &
                & nocol     ,irocol(1, norow + 1) ,gdp       )
       call timer_stop(timer_1stuzd, gdp)
       call timer_stop(timer_uzd, gdp)
       !
       !     computation proceeds in X direction
       !
       ! CHECK FOR FLOODING AND DRYING IN "U" POINTS
       !
       ! Update HU, DZU1. DZS1 is updated in Z_DRYCHK after Z_SUD based
       ! on the computed water level at new time step
       !
       flood = .true.
       idry  = 0
       icx   = nmaxddb
       icy   = 1
       call timer_start(timer_checku, gdp)
       call z_checku(j         ,nmmaxj    ,nmmax     ,icx       ,kmax      , &
                   & flood     , &
                   & kfu       ,kcs       ,kcu       ,kspu      , &
                   & kfuz1     ,kfumin    ,kfumax    ,kfumx0    ,hu        , &
                   & s0        ,dpu       ,dps       ,umean     ,u0        , &
                   & u1        ,dzu1      ,zk        ,gdp       )
       call timer_stop(timer_checku, gdp)
       !
       ! Computation of U1 and S1, i.e. evaluation of coupled momentum and
       ! continuity equation for one half time step. 
       !
       ! HU is updated by routine UPWHU in CHECKU
       !
11111  continue
       call timer_start(timer_sud, gdp)
       !
       gdp%dd%suditer = 0
       icx = nmaxddb
       icy = 1
       call timer_start(timer_1stsud, gdp)
       call z_sud(j         ,nmmaxj    ,nmmax     ,kmax      ,mmax      , &
                & nmax      ,nsrc      ,nst       ,icx       ,icy       , &
                & flood     ,norow     ,irocol(1, 1)         ,mnksrc    , &
                & kfu       ,kfv       ,kfs       ,kcs       ,kcu       , &
                & kfuz1     ,kfvz1     ,kfsz1     ,kspu      ,kcu45     , &
                & kcscut    ,kfumin    ,kfumax    ,kfsmin    ,kfsmax    , &
                & kfumx0    ,thick     ,circ2d(1, 1)         ,circ3d(1, 1, 1)      ,s0        , &
                & s1        ,u0        ,u1        ,v1        ,w1        , &
                & qxk       ,qyk       ,qzk       ,guu       ,gvv       , &
                & guv       ,gvu       ,gud       ,gvd       ,guz       , &
                & gvz       ,gsqiu     ,gsqs      ,disch     ,umdis     , &
                & dismmt    ,umean     ,evap      ,hu        ,dps       , &
                & dpu       ,dzu1      ,dzs1      ,dzu0      ,wrka1     , &
                & wrka2     ,wrka3     ,wrka4     ,wrka5     ,wrka6     , &
                & wrka7     ,wrka8     ,wrka15    ,wrkb1     ,wrkb2     , &
                & wrkb3     ,wrkb4     ,wrkb5     ,wrkb6     ,wrkb7     , &
                & wrkb8     ,wsu       ,taubpu    ,taubsu    ,vicuv     , &
                & vnu2d     ,vicww     ,rxx       ,rxy       ,windu     , &
                & patm      ,fcorio    ,tgfsep    ,drhodx    ,zk        , &
                & p0        ,crbc(1, 1),idry      ,porosu    ,ubrlsu    , &
                & pship     ,diapl     ,rnpl      ,cfurou    ,gdp)
       !
       call timer_stop(timer_1stsud, gdp)
       call timer_stop(timer_sud, gdp)
       !
       ! Check for drying in waterlevel points in the X-direction
       ! Update DZS1
       !
       icx = nmaxddb
       icy = 1
       call timer_start(timer_drychk, gdp)
       call z_drychk(idry      ,j         ,nmmaxj    ,nmmax     ,kmax      , &
                   & nfltyp    ,icx       ,icy       ,kfu       ,kfv       , &
                   & kfs       ,kcs       ,kfuz1     ,kfvz1     ,kfsz1     , &
                   & kfumin    ,kfumax    ,kfvmin    ,kfvmax    ,kfsmin    , &
                   & kfsmax    ,kfsmx0    ,s1        ,s0        ,dps       , &
                   & qxk       ,qyk       ,w1        ,dzu1      ,dzv1      , &
                   & dzs1      ,zk        ,gdp       )
       call timer_stop(timer_drychk, gdp)
       !
       ! Compute Volume and Areas to be used in routines that computes 
       ! the transport of matter (consistency with WAQ)
       !
       ! For areas use the DZU array from previous time step (DZU0)
       !
       call timer_start(timer_comvol, gdp)
       call comvol(nmmax     ,kmax      ,zmodel    ,kcs       ,kcu       , &
                 & thick     ,guu       ,gsqs      ,dps       ,s1        , &
                 & dzs1      ,dzu0      ,hu        ,porosu    ,volum1    , &
                 & areau     ,gdp       )
       call timer_stop(timer_comvol, gdp)
       !
       !
       ! DD code added:
       !
       !
       ! Synchronize on Dry Point
       !
       if (nfltyp==0) then
          nhystp = nxtdry(d3dflow_check_adi_dry, 0, gdp)
       else
          nhystp = nxtdry(d3dflow_check_adi_dry, idry, gdp)
       endif
       !
       if (nfltyp/=0) then
          !
          ! If waterlevel is below bottom then isolate waterlevel point by setting
          ! the surrounding velocities to zero and repeat the computation of SUD
          !
          if (nhystp==d3dflow_build_adi_zeta .or.                               &
            & (nhystp==noneighbors .and. idry==1)) goto 11111
       !
       ! DD code added end
       !
       endif
    !
    ! END OF COMPUTATION OF STAGE 1 FOR ADI METHOD
    !
    endif
    !
    ! =====================================
    ! COMPUTATION OF STAGE 2 FOR ADI METHOD
    ! =====================================
    !
    if (stage=='stage2') then
       !
       ! Computation of U1, i.e. evaluate momentum equation for one half timest
       !     calculate hu and set kfu = 0 for hu < htrsh (.5*dryflc)
       !
       call timer_start(timer_uzd, gdp)
       !
       ! dummy z_uzd call (for synchronisation with domain decomposition
       ! implementation of the sigma-model)
       !
       gdp%dd%uzditer = 0
       icx = nmaxddb
       icy = 1
       call timer_start(timer_2nduzd, gdp)
       call z_uzd(j         ,nmmaxj    ,nmmax     ,kmax      ,icx       , &
                & icy       ,nsrc      ,kcs       ,kcu45     ,kcscut    , &
                & kcu       ,kfu       ,kfuz1     ,kfumin    ,kfumax    , &
                & kfv       ,kfvz1     ,kfs       ,kfsz1     ,kfsmin    , kfsmax    , &
                & u0        ,v1        ,w1        ,hu        ,dzu1      , &
                & guu       ,gvv       ,gvu       ,guv       ,gsqs      , &
                & gud       ,gvd       ,guz       ,gvz       ,gsqiu     , &
                & disch     ,umdis     ,kspu      ,mnksrc    ,dismmt    , &
                & wrkb1     ,wrkb2     ,wrkb3     ,wrkb4     ,wrkb5     , &
                & wrkb6     ,wrkb7     ,wrkb8     , &
                & vicuv     ,vnu2d     ,vicww     ,tgfsep    , &
                & drhodx    ,wsu       ,taubpu    ,taubsu    ,rxx       , &
                & rxy       ,windu     ,patm      ,fcorio    ,p0        , &
                & ubrlsu    ,pship     ,diapl     ,rnpl      ,nst       , &
                & u1        ,s0        ,cfurou    ,wrkb13    ,wrkb14    , &
                & norow     ,irocol(1, 1)         ,gdp       )                
       call timer_stop(timer_2nduzd, gdp)
       call timer_stop(timer_uzd, gdp)
       !
       ! CHECK FOR FLOODING AND DRYING IN "V" POINTS
       !
       flood = .true.
       idry  = 0
       icx   = 1
       icy   = nmaxddb
       call timer_start(timer_checku, gdp)
       call z_checku(j         ,nmmaxj    ,nmmax     ,icx       ,kmax      , &
                   & flood     , &
                   & kfv       ,kcs       ,kcv       ,kspv      , &
                   & kfvz1     ,kfvmin    ,kfvmax    ,kfvmx0    ,hv        , &
                   & s0        ,dpv       ,dps       ,vmean     ,v0        , &
                   & v1        ,dzv1      ,zk        ,gdp       )
       call timer_stop(timer_checku, gdp)
       !
       ! Computation of V1 and S1, i.e. evaluation of coupled momentum and
       ! continuity equation for one half time step
       !
22222  continue
       call timer_start(timer_sud, gdp)
       !
       gdp%dd%suditer = 0
       icx = 1
       icy = nmaxddb
       call timer_start(timer_2ndsud, gdp)
       call z_sud(j         ,nmmaxj    ,nmmax     ,kmax      ,nmax      , &
                & mmax      ,nsrc      ,nst       ,icx       ,icy       , &
                & flood     ,nocol     ,irocol(1, norow + 1) ,mnksrc    , &
                & kfv       ,kfu       ,kfs       ,kcs       ,kcv       , &
                & kfvz1     ,kfuz1     ,kfsz1     ,kspv      ,kcv45     , &
                & kcscut    ,kfvmin    ,kfvmax    ,kfsmin    ,kfsmax    , &
                & kfvmx0    ,thick     ,circ2d(1, norow + 1) ,circ3d(1, 1, norow + 1)         ,s0        , &
                & s1        ,v0        ,v1        ,u1        ,w1        , &
                & qyk       ,qxk       ,qzk       ,gvv       ,guu       , &
                & gvu       ,guv       ,gvd       ,gud       ,gvz       , &
                & guz       ,gsqiv     ,gsqs      ,disch     ,vmdis     , &
                & dismmt    ,vmean     ,evap      ,hv        ,dps       , &
                & dpv       ,dzv1      ,dzs1      ,dzv0      ,wrka1     , &
                & wrka2     ,wrka3     ,wrka4     ,wrka5     ,wrka6     , &
                & wrka7     ,wrka8     ,wrka16    ,wrkb1     ,wrkb2     , &
                & wrkb3     ,wrkb4     ,wrkb5     ,wrkb6     ,wrkb7     , &
                & wrkb8     ,wsv       ,taubpv    ,taubsv    ,vicuv     , &
                & vnu2d     ,vicww     ,ryy       ,rxy       ,windv     , &
                & patm      ,fcorio    ,tgfsep    ,drhody    ,zk        , &
                & p0        ,crbc(1, norow + 1)   ,idry      ,porosv    ,ubrlsv    , &
                & pship     ,diapl     ,rnpl      ,cfvrou    ,gdp)
       call timer_stop(timer_2ndsud, gdp)
       call timer_stop(timer_sud, gdp)
       !
       ! Check for drying in waterlevel points in the X-direction
       !
       icx = nmaxddb
       icy = 1
       call timer_start(timer_drychk, gdp)
       call z_drychk(idry      ,j         ,nmmaxj    ,nmmax     ,kmax      , &
                   & nfltyp    ,icx       ,icy       ,kfu       ,kfv       , &
                   & kfs       ,kcs       ,kfuz1     ,kfvz1     ,kfsz1     , &
                   & kfumin    ,kfumax    ,kfvmin    ,kfvmax    ,kfsmin    , &
                   & kfsmax    ,kfsmx0    ,s1        ,s0        ,dps       , &
                   & qxk       ,qyk       ,w1        ,dzu1      ,dzv1      , &
                   & dzs1      ,zk        ,gdp       )
       call timer_stop(timer_drychk, gdp)
       !
       ! Compute Volume and Areas to be used in routines that computes 
       ! the transport of matter (consistency with WAQ)
       !
       call timer_start(timer_comvol, gdp)
       call comvol(nmmax     ,kmax      ,zmodel    ,kcs       ,kcv       , &
                 & thick     ,gvv       ,gsqs      ,dps       ,s1        , &
                 & dzs1      ,dzv0      ,hv        ,porosv    ,volum1    , &
                 & areav     ,gdp       )
       call timer_stop(timer_comvol, gdp)
       !
       !
       ! DD code added:
       !
       !
       ! Synchronize on Dry Point
       !
       if (nfltyp==0) then
          nhystp = nxtdry(d3dflow_check_adi_dry, 0, gdp)
       else
          nhystp = nxtdry(d3dflow_check_adi_dry, idry, gdp)
       endif
       !
       if (nfltyp/=0) then
          !
          ! If waterlevel is below bottom then isolate waterlevel point by setting
          ! the surrounding velocities to zero and repeat the computation of SUD
          !
          if (nhystp==d3dflow_build_adi_zeta .or.                               &
            & (nhystp==noneighbors .and. idry==1)) goto 22222
       !
       ! DD code added end
       !
       endif
    !
    ! END OF COMPUTATION OF STAGE 2 FOR ADI METHOD
    !
    endif
end subroutine z_adi
