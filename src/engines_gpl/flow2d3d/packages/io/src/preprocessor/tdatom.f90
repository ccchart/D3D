subroutine tdatom(runid, filmrs, nuerr, alone, gdp) 
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2012.                                
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
!    Function: Reads and writes the time dependent data from the 
!              MD-file and or attribute file to standard unfor- 
!              matted file for the simulation 
!              Optional extra function for astronomical tide 
!              Create extra intermediate file which smallest 
!              end time of all time dependent data 
! Method used: 
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
    ! parameters (hard coded, fixed dimensions)
    !
    include 'pardef.igd'         
    !
    type(globdat),target :: gdp 
    ! 
    ! The following list of pointer parameters is used to point inside the gdp structure 
    ! 
    integer                         , pointer :: nmax 
    integer                         , pointer :: mmax 
    integer                         , pointer :: nmaxus 
    integer                         , pointer :: kmax 
    integer                         , pointer :: lsts 
    integer                         , pointer :: lstsc 
    integer                         , pointer :: lstsci 
    integer                         , pointer :: lsal 
    integer                         , pointer :: ltem 
    integer                         , pointer :: lsec 
    integer                         , pointer :: ltur 
    integer                         , pointer :: ltur2d 
    integer                         , pointer :: nostat 
    integer                         , pointer :: nto 
    integer                         , pointer :: ntof 
    integer                         , pointer :: ntoq 
    integer                         , pointer :: ntot 
    integer                         , pointer :: ntruv 
    integer                         , pointer :: kc 
    integer                         , pointer :: nsrc 
    integer                         , pointer :: nsrcd 
    integer                         , pointer :: nsluv 
    integer                         , pointer :: upwsrc 
    integer                         , pointer :: itdate 
    real(fp)                        , pointer :: tstart 
    real(fp)                        , pointer :: tstop 
    real(fp)                        , pointer :: dt 
    real(fp)                        , pointer :: tunit 
    real(fp)                        , pointer :: tzone 
    integer                         , pointer :: lunmd 
    integer                         , pointer :: lundia 
    integer                         , pointer :: lunscr 
    integer                         , pointer :: itstrt 
    integer                         , pointer :: itfinish 
    integer                         , pointer :: julday 
    character(256)                  , pointer :: nflmod 
    logical                         , pointer :: wind 
    logical                         , pointer :: salin 
    logical                         , pointer :: temp 
    logical                         , pointer :: const 
    logical                         , pointer :: culvert 
    logical                         , pointer :: dredge 
    logical                         , pointer :: drogue 
    logical                         , pointer :: wave 
    logical                         , pointer :: waveol 
    logical                         , pointer :: threed 
    logical                         , pointer :: secflo 
    logical                         , pointer :: iweflg 
    logical                         , pointer :: struct 
    logical                         , pointer :: cdwstruct 
    logical                         , pointer :: sedim 
    logical                         , pointer :: htur2d 
    logical                         , pointer :: flmd2l 
    logical                         , pointer :: mudlay 
    logical                         , pointer :: couplemod 
    logical                         , pointer :: zmodel 
    logical                         , pointer :: nonhyd 
    logical                         , pointer :: roller 
    logical                         , pointer :: wavcmp 
    logical                         , pointer :: cnstwv 
    logical                         , pointer :: dpmveg 
    logical                         , pointer :: snelli 
    logical                         , pointer :: lrdamp 
    logical                         , pointer :: sbkol 
    logical                         , pointer :: bubble 
    logical                         , pointer :: nfl 
    integer                         , pointer :: rtcmod 
    logical                         , pointer :: reusetmp !!  TRUE when temporary files will be reused 
                                                          !!  if possible 
    real(fp)                        , pointer :: zbot 
    real(fp)                        , pointer :: ztop 
! 
! Local parameters 
! 
    integer, parameter :: mxtime = 1 
! 
! Global variables 
! 
    integer        , intent(out) :: nuerr  !! Exit code: 0 := ok, < 0 then error
    logical        , intent(in)  :: alone  !! TRUE when flow runs stand-alone, 
                                           !! FALSE when flow is part of morsys     
    character(12)  , intent(in)  :: filmrs !! File name for DELFT3D_MOR FLOW 
                                           !! input file (MD-flow.xxx) 
    character(*)   , intent(in)  :: runid  !! Run identification code for the current simulation
                                           !! (used to determine the names of the in/output files used by the system)
! 
! Local variables 
! 
    integer                                             :: ibarco      ! Barocline pressure correction for open boundary points 0 = no 1 = yes Default 
    integer                                             :: icreep      ! Identification for special approach horizontal gradients 0 = no anti creep 1 = yes anti creep  
    integer                                             :: iend        ! Help var.  
    integer                                             :: irov 
    integer                                             :: idensform 
    integer                                             :: it 
    integer                                             :: iter1       ! Loop parameter for iteration over continuity equation  
    integer                                             :: ivapop      ! Flag specifying whether vapour pres- sure data are to be computed or specified Only for KTEMP=4 IVAPOP can be 1  
    integer                                             :: j           ! Help var.  
    integer                                             :: k 
    integer                                             :: keva        ! Option flag for rain/evaporation module =1: Yes =0: No  
    integer                                             :: ktemp       ! Option flag for temperature module 
                                                                       !    =1: T-abs. model, radiation from eq. 
                                                                       !    =2: T-abs. model, radiation pre- scribed by time series 
                                                                       !    =3: Excess temperature model (1 and 2 not yet implemented)  
    integer                                             :: l           ! Help var.  
    integer                                             :: lconc       ! Number of constituents defined by user (excl. Salinity, Temperature, Secondary flow and quantities for the Turb. models)  
    integer                                             :: lkw         ! Length of char. str (usually the KEYWRD or RECNAM)  
    integer                                             :: lrid        ! Length of character string runid  
    integer                                             :: lunbch      ! Unit number for temporary harmonic or astronomical components data file  
    integer                                             :: luntdi      ! Unit number containing start time for all time dependent data if requested  
    integer                                             :: luntri      ! Unit number for trigger file for TRISIM for running programs simultaniously all time dependent data if requested  
    integer                                             :: merr        ! Help var. (nr. of 'ERROR'   string found in the diagnostic file)  
    integer                                             :: mwarn       ! Help var. (nr. of 'WARNING' string found in the diagnostic file)  
    integer                                             :: mxdnto      ! Array dimension for opening boundary definition arrays, which are dynamic declared (MXNTO in call to RDBNDD)  
    integer                                             :: n 
    integer                                             :: nn          ! Help var.  
    integer                                             :: nrrec       ! Pointer to the record number in the MD-file  
    integer                                             :: nrver       ! Version number (240 - ....)  
    integer                                             :: ntimtm      ! Actual number of times for which time varying data is allowed in the Md-file (dummy value)  
    integer                                             :: ntot0       ! Offset for time series bnd. in arrays NTOT0 = NTOF + NTOQ  
    integer                                             :: ntrec       ! Help. var to keep track of NRREC  
    integer, dimension(:,:), pointer                    :: mnbnd       ! Array containing the coordinates of the open boundary sections 
                                                                       !    MNBND(1,K)=X-Coor. of the begin pnt. 
                                                                       !    MNBND(2,K)=X-Coor. of the end   pnt. 
                                                                       !    MNBND(3,K)=Y-Coor. of the begin pnt. 
                                                                       !    MNBND(4,K)=Y-Coor. of the end   pnt. 
                                                                       !    K = 1,.....,NOPEN  
    integer, dimension(:,:), pointer                    :: mnksrc      ! MNK-coord. for discharges  
    integer, dimension(:,:,:), pointer                  :: kspu        ! NT+2,0:MXKMAX Mask array for total water depth upwind in special (u-)points 
                                                                       !    In KSPU(NM,0) the special point definition is given 
                                                                       !       = 1 Discharge location 
                                                                       !       = 2 Floating structure 
                                                                       !       = 3 Local weir 
                                                                       !       = 4 Gate 
                                                                       !       = 5 Rigid sheet 
                                                                       !       = 6 Porous plate 
                                                                       !       = 7 Bridge 
                                                                       !       = 8 Barrier 
                                                                       !       = 9 2D Weir 
                                                                       !    For type 1-3,5-8 the negative equivalence implice no upwind  
    integer, dimension(:,:,:), pointer                  :: kspv        ! NT+2,0:MXKMAX Mask array for total water depth upwind in special (v-)points 
                                                                       !    In KSPV(NM,0) the special point definition is given  
                                                                       !       = 1 Discharge location 
                                                                       !       = 2 Floating structure 
                                                                       !       = 3 Local weir 
                                                                       !       = 4 Gate 
                                                                       !       = 5 Rigid sheet 
                                                                       !       = 6 Porous plate 
                                                                       !       = 7 Bridge 
                                                                       !       = 8 Barrier 
                                                                       !       = 9 2D Weir 
                                                                       !    For type 1-3,5-8 the negative equivalence implice no upwind  
    integer, dimension(:), pointer                      :: nhsub       ! integer array to store sequence numbers of harmonic boundary condition in own subdomain 
    integer, external                                   :: newlun 
    logical                                             :: found       ! If FOUND = TRUE then recnam in the MD-file was found  
    logical                                             :: error       ! Flag=TRUE if an error is encountered  
    logical                                             :: solrad_read ! Flag=TRUE means Nett Solar Radiation is to be read from .tem file 
    logical                                             :: lexist      ! Logical to determine file existence  
    logical                                             :: newkw       ! Logical var. specifying whether a new recnam should be read from the MD-file or just new data in the continuation line  
    logical                                             :: noui        ! Flag for reading from User Interface  
    logical                                             :: sferic      ! Flag for spherical coordinates (TRUE or FALSE)  
    logical                                             :: tmpexist    ! Flag for call from TDATOM (.true.) for time varying data  
    logical                                             :: verify      ! Always FALSE, to be removed; was used for program=MD-VER
    logical                                             :: yestdd      ! Flag for call from TDATOM (.true.) for time varying data  
    real(fp)                                            :: anglat 
    real(fp)                                            :: betac       ! Coupling coefficient between intensity and impuls equations for secondary flows  
    real(fp)                                            :: chzmin 
    real(fp)                                            :: dco         ! Minimal depth for upwind approach in shallow areas [m]  
    real(fp)                                            :: dgcuni 
    real(fp)                                            :: dryflc      ! Drying & Flooding criterion [m]  
    real(fp)                                            :: dx 
    real(fp)                                            :: dy 
    real(fp)                                            :: fwfac 
    real(fp)                                            :: gammax 
    real(fp)                                            :: grdang 
    real(fp)                                            :: rhoa        ! Density of air (kg/m3)  
    real(fp)                                            :: rmincf 
    real(fp)                                            :: saleqs      ! Salinity value used in eq. of state which will be applied uniformly over the vertical.  
    real(fp)                                            :: temeqs      ! Temperature value used in the eq. of state which will applied uniformly over the vertical.  
    real(fp)                                            :: z0v 
    real(fp), dimension(:), pointer                     :: wstcof      ! Wind stress Coefficients (constant)
                                                                       !    Space varying: 1 - wstcof at wspeed1
                                                                       !                   2 - wspeed1
                                                                       !                   3 - wstcof at wspeed2
                                                                       !                   4 - wspeed2
                                                                       !                   5 - wstcof at wspeed3
                                                                       !                   6 - wspeed3
                                                                       !    Uniform: 1 - wstcof at wspeed1
                                                                       !             2 - 0.0 m/s
                                                                       !             3 - wstcof(2)
                                                                       !             4 - 50.0 m/s
                                                                       !             5 - wstcof(3)
                                                                       !             6 - 100.0 m/s

    real(fp), dimension(:,:,:), pointer                 :: hydrbc 
    real(fp), dimension(:,:,:,:), pointer               :: rval        ! ,MXLLLL Help array can have at most 4 array subscripts  
    real(fp), dimension(:), pointer                     :: omega 
    real(fp), dimension(:), pointer                     :: thick       ! Relative layer thickness  
    real(fp), dimension(:), pointer                     :: alpha 
    real(fp), dimension(:), pointer                     :: rtime 
    character(1)                                        :: ascon       ! 'Y' if open bnd. contains type 'A'  
    character(1)                                        :: ctunit      ! Time scale for time parameters, currently set to 'M'(inute - fixed).  
    character(1)                                        :: eol         ! ASCII code for End-Of-Line (^J)  
    character(1)                                        :: equili      ! Equilibrium or advection and diffusion default = no equilibrium ('N') which means lsec = 1  
    character(1)                                        :: evaint      ! Interpolation option for the rainfall / evaporation  
    character(1)                                        :: forfuv      ! Forester filter option for UV  
    character(1)                                        :: forfww      ! Forester filter option for W  
    character(1)                                        :: sphere      ! Flag Yes / No spherical coordinates  
    character(1)                                        :: temint      ! Interpolation option for the tempe- rature  
    character(1), dimension(:), pointer                 :: disint      ! Interpolation option for the disch.  
    character(1), dimension(:), pointer                 :: datbnd      ! Type of open boundary: -'A'(stronomical) -'H'(armonic/Tide) -'T'(ime series/time dependent)  
    character(1), dimension(:), pointer                 :: typbnd      ! Type of open boundary prescribed 
    character(10)                                       :: citdat      ! Reference date for the simulation times. Format: "DD MMM 'YY"  
    character(10)                                       :: date        ! Date to be filled in the header  
    character(30), dimension(:), pointer                :: runtxt      ! Textual description of model input  
    character(10), dimension(:,:),pointer               :: cval        ! Help array for character variables  
    character(11)                                       :: fmtfil      ! File format for the time varying data file  
    character(12)                                       :: filsim      ! Name for trigger file for running FLOW and RTC simultaniously  
    character(12), dimension(:,:),pointer               :: statns      ! References to tidal stations at boundary support points  
    character(13)                                       :: trasol      ! Transport scheme option  
    character(6)                                        :: momsol 
    character(16)                                       :: simdat      ! Dummy string 
    character(20)                                       :: rundat      ! Current date and time containing a combination of DATE and TIME  
    character(20), dimension(:),pointer                 :: namcon      ! Names of the constituents  
    character(20), dimension(:),pointer                 :: namsrc      ! Names of discharge points the open boundary section  
    character(20), dimension(:),pointer                 :: nambnd      ! Names of open boundaries  
    character(20), dimension(:),pointer                 :: tprofu 
    character(256)                                      :: filnam      ! File name for the time varying data file  
    character(256)                                      :: filmd       ! File name of MD-Flow file 
    character(300)                                      :: mdfrec      ! Standard rec. length in MD-file (300) 300 = 256 + a bit (field, =, ##, etc.)  
    character(36)                                       :: tgfcmp      ! Character containing names of tidal constant names to be included in the Tide generating forces. If blank then No forces will be included.  
    character(8)                                        :: dpsopt 
    character(8)                                        :: dpuopt 
    character(55)                                       :: txtput 
    character(6)                                        :: soort       ! Help var. determining the prog. name currently active  
    character(9)                                        :: keyw        ! Name of record to look for in the MD-file 
    character(5)                                        :: versio 
! 
!! executable statements ------------------------------------------------------- 
! 
    nmax        => gdp%d%nmax 
    mmax        => gdp%d%mmax 
    nmaxus      => gdp%d%nmaxus 
    kmax        => gdp%d%kmax 
    lsts        => gdp%d%lsts 
    lstsc       => gdp%d%lstsc 
    lstsci      => gdp%d%lstsci 
    lsal        => gdp%d%lsal 
    ltem        => gdp%d%ltem 
    lsec        => gdp%d%lsec 
    ltur        => gdp%d%ltur 
    ltur2d      => gdp%d%ltur2d 
    nostat      => gdp%d%nostat 
    nto         => gdp%d%nto 
    ntof        => gdp%d%ntof 
    ntoq        => gdp%d%ntoq 
    ntot        => gdp%d%ntot 
    ntruv       => gdp%d%ntruv 
    kc          => gdp%d%kc 
    nsrc        => gdp%d%nsrc 
    nsrcd       => gdp%d%nsrcd 
    nsluv       => gdp%d%nsluv 
    upwsrc      => gdp%d%upwsrc 
    itdate      => gdp%gdexttim%itdate 
    tstart      => gdp%gdexttim%tstart 
    tstop       => gdp%gdexttim%tstop 
    dt          => gdp%gdexttim%dt 
    tunit       => gdp%gdexttim%tunit 
    tzone       => gdp%gdexttim%tzone 
    lunmd       => gdp%gdinout%lunmd 
    lundia      => gdp%gdinout%lundia 
    lunscr      => gdp%gdinout%lunscr 
    itstrt      => gdp%gdinttim%itstrt 
    itfinish    => gdp%gdinttim%itfinish 
    julday      => gdp%gdinttim%julday 
    nflmod      => gdp%gdnfl%nflmod 
    wind        => gdp%gdprocs%wind 
    salin       => gdp%gdprocs%salin 
    temp        => gdp%gdprocs%temp 
    const       => gdp%gdprocs%const 
    culvert     => gdp%gdprocs%culvert 
    dredge      => gdp%gdprocs%dredge 
    drogue      => gdp%gdprocs%drogue 
    wave        => gdp%gdprocs%wave 
    waveol      => gdp%gdprocs%waveol 
    threed      => gdp%gdprocs%threed 
    secflo      => gdp%gdprocs%secflo 
    iweflg      => gdp%gdprocs%iweflg 
    struct      => gdp%gdprocs%struct 
    cdwstruct   => gdp%gdprocs%cdwstruct 
    sedim       => gdp%gdprocs%sedim 
    htur2d      => gdp%gdprocs%htur2d 
    flmd2l      => gdp%gdprocs%flmd2l 
    mudlay      => gdp%gdprocs%mudlay 
    couplemod   => gdp%gdprocs%couplemod 
    zmodel      => gdp%gdprocs%zmodel 
    nonhyd      => gdp%gdprocs%nonhyd 
    roller      => gdp%gdprocs%roller 
    wavcmp      => gdp%gdprocs%wavcmp 
    cnstwv      => gdp%gdprocs%cnstwv 
    dpmveg      => gdp%gdprocs%dpmveg 
    snelli      => gdp%gdprocs%snelli 
    lrdamp      => gdp%gdprocs%lrdamp 
    sbkol       => gdp%gdprocs%sbkol 
    bubble      => gdp%gdprocs%bubble 
    nfl         => gdp%gdprocs%nfl 
    rtcmod      => gdp%gdrtc%rtcmod 
    reusetmp    => gdp%gdtmpfil%reusetmp 
    zbot        => gdp%gdzmodel%zbot 
    ztop        => gdp%gdzmodel%ztop 
    ! 
    dpsopt      = ' ' 
    dpuopt      = ' ' 
    nflmod      = ' ' 
    ! 
    error       = .false. 
    noui        = .true. 
    yestdd      = .true. 
    found       = .false. 
    verify      = .false. 
    ! 
    soort       = 'tdatom' 
    cdwstruct   = .false. 
    cnstwv      = .false. 
    const       = .false. 
    couplemod   = .false. 
    culvert     = .false. 
    dpmveg      = .false. 
    dredge      = .false. 
    drogue      = .false. 
    htur2d      = .false. 
    iweflg      = .false. 
    lrdamp      = .false. 
    roller      = .false. 
    rtcmod      =  noRTC 
    salin       = .false. 
    secflo      = .false. 
    sedim       = .false. 
    sferic      = .false. 
    snelli      = .false. 
    struct      = .false. 
    temp        = .false. 
    threed      = .false. 
    wavcmp      = .false. 
    wave        = .false. 
    waveol      = .false. 
    wind        = .false. 
    nfl         = .false. 
    ! 
    ! Initializing tdatom part of FLOW simulation program 
    ! 
    call sysini(error     ,runid     ,filmrs    ,alone     ,soort     , & 
              & verify    ,versio    ,filmd     ,gdp       ) 
    ! 
    if (error) goto 9990
    !    
    ! Check if all necessary temporary files exist if they are reused 
    ! 
    call tmpcheck( runid,  reusetmp,  tmpexist,  gdp ) 
    ! 
    ! Run TDATOM 
    ! But not if temporary files are reused and they all exist 
    ! 
    if (.not.tmpexist) then 
        soort = 'tdatom' 
        ! 
        ! Define EOL 
        ! 
        eol = char(10) 
        ! 
        ! Define length of runid 
        ! 
        call noextspaces(runid     ,lrid      ) 
        ! 
        ! Read processes and dimensions 
        ! 
        call dimrd(lunmd     ,lundia    ,error     ,runid     ,nrver     , & 
                 & soort     ,wind      ,salin     ,temp      ,sedim     , & 
                 & const     ,secflo    ,drogue    ,wave      ,iweflg    , & 
                 & htur2d    ,mudlay    , & 
                 & flmd2l    ,zmodel    ,nonhyd    ,roller    ,wavcmp    , & 
                 & culvert   ,dredge    ,cdwstruct ,snelli    ,cnstwv    , & 
                 & dpmveg    ,waveol    ,lrdamp    ,sbkol     ,bubble    , & 
                 & nfl       ,nflmod    ,gdp       ) 
        if (error) goto 9990 
        ! 
        allocate(hydrbc(4, mxnto, mxkc))
        allocate(rval(4, mxtime, mxnto, lstsc)) 
        allocate(cval(mxnto, lstsc)) 
        allocate(namcon(lstsc)) 
        !
        allocate(mnbnd (7, mxnto))
        allocate(mnksrc(7, mxnsrc))
        allocate(kspu  (mxnpnt, -1:mxnpnt + 2, 0:mxkmax))
        allocate(kspv  (mxnpnt, -1:mxnpnt + 2, 0:mxkmax))
        allocate(nhsub (mxnto))
        allocate(wstcof(6))
        allocate(omega (mxkc))
        allocate(thick (mxkmax))
        allocate(alpha (mxnto ))
        allocate(rtime (mxtime))
        allocate(disint(mxnsrc))
        allocate(datbnd(mxnto ))
        allocate(typbnd(mxnto ))
        allocate(runtxt(10    ))
        allocate(statns(mxnto,2))
        allocate(namsrc(mxnsrc))
        allocate(nambnd(mxnto))
        allocate(tprofu(mxnto))
        ! 
        do j = 1, 4 
           do l = 1, lstsc 
              do it = 1, mxtime 
                 do nn = 1, mxnto 
                    cval(nn, l) = ' ' 
                    rval(j, it, nn, l) = 0.0 
                 enddo 
              enddo 
           enddo 
        enddo 
        ! 
        ! Test local dimensions, and define LCONC 
        ! 
        call chklod(lundia    ,error     ,nto       ,kmax      ,nsrc      , & 
                  & gdp       ) 
        if (error) goto 999 
        ! 
        lconc = lstsc - lsal - min(1, ltem) 
        ! 
        ! Define Constituent names depending on SALIN, TEMP and CONST 
        ! 
        rewind (lunmd) 
        read (lunmd, '(a300)') mdfrec 
        nrrec = 1 
        ! 
        call rdnamc(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                  & noui      ,salin     ,temp      ,lconc     ,lstsc     , & 
                  & namcon    ,gdp       ) 
        if (error) goto 999 
        ! 
        ! X, Y, Z and Orientation, initialize global parameters 
        ! Initialize global and local parameters in subroutine 
        ! 
        !call rdxyzo(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
        !          & noui      ,kmax      ,zbot      ,ztop      , & 
        !          & dx        ,dy        ,filnam    ,fmtfil    ,thick     , & 
        !          & anglat    ,grdang    ,sphere    ,sferic    ,zmodel    , & 
        !          & gdp       ) 
        !if (error) goto 999 
        ! 
        ! Read the simulation times 
        ! initialize global and local parameters in subroutine 
        ! 
        call rdirt(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                 & noui      ,citdat    ,tstart    ,tstop     ,tzone     , & 
                 & itdate    ,julday    ,itstrt    ,itfinish  ,dt        , & 
                 & ctunit    ,tunit     ,gdp       ) 
        if (error) goto 999 
        ! 
        ! Read KTEMP (Process heat module) 
        ! initialize global and local parameters in subroutine 
        ! 
        call rdproc(error    ,nrrec     ,mdfrec   ,noui        ,htur2d      , & 
                  & salin    ,temp      ,wind     ,ktemp       , & 
                  & keva     ,ivapop    ,irov     ,ctunit      , &
                  & z0v      ,sferic    ,tgfcmp   ,temeqs      ,saleqs      , & 
                  & wstcof   ,rhoa      ,secflo   ,betac       ,equili      , & 
                  & lsec     ,chzmin    ,rmincf   ,rtcmod      ,couplemod   , & 
                  & nonhyd   ,mmax      ,nmax     ,nmaxus      ,sedim       , & 
                  & idensform ,solrad_read , gdp) 
        if (error) goto 999 
        ! 
        ! Initialize global and local parameters in subroutine 
        ! 
        call rdnum(lunmd     ,lundia    ,nrrec     ,mdfrec    , & 
                 & iter1     ,dryflc    ,dco       ,ibarco    ,kmax      , & 
                 & lstsci    ,icreep    ,trasol    ,momsol    ,dgcuni    , & 
                 & forfuv    ,forfww    ,ktemp     ,temint    ,            & 
                 & keva      ,evaint    , & 
                 & dpsopt    ,dpuopt    ,zmodel    ,gammax    ,fwfac     , & 
                 & gdp       ) 
        ! 
        ! Read boundary definition. If specified in a separate file, 
        ! then check if TRIANA file is specified. If so, start tidals, the 
        ! subroutine that reads tidal components in TRIANA file, translates 
        ! cmp names to frequencies, performs ascon correction, and writes 
        ! the freq., ampl. and phas. at their proper positions. 
        ! Initialisation of global and local parameters in subroutine 
        ! 
        if (nto > 0) then 
           mxdnto = mxnto 
           call rdbndd(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                     & noui      ,nto       ,ntof      ,ntoq      ,mmax      , & 
                     & nmaxus    ,kmax      ,mxdnto    ,mxnto     ,filnam    , & 
                     & fmtfil    ,ascon     ,nambnd    ,typbnd    ,datbnd    , & 
                     & mnbnd     ,alpha     ,tprofu    ,statns    ,nhsub     , & 
                     & yestdd    , gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Special points, discharge sources 
        ! initialize global and local parameters in subroutine 
        ! 
        call rdspec(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                  & noui      ,yestdd    ,filnam    ,fmtfil    ,nsrcd     , & 
                  & mmax      ,nmax      ,nmaxus    ,mnksrc    ,namsrc    , & 
                  & disint    ,upwsrc    ,gdp       ) 
        if (error) goto 999 
        ! 
        ! Put header on the screen 
        ! 
        txtput = 'Part II   - Creating intermediate files...' 
        write (lunscr, '(a)') txtput 
        ! 
        ! Periodic boundary conditions 
        ! only if NTOF > 0 
        ! 
        if (ntof > 0) then 
           filnam = ' ' 
           fmtfil = 'FR' 
           ntimtm = 0 
           ! 
           do k = 1, mxkc 
              omega(k) = 0. 
              do n = 1, mxnto 
                 do j = 1, 4 
                    hydrbc(j, n, k) = 0.0 
                 enddo 
              enddo 
           enddo 
           ! 
           call rdbch(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & noui      ,filnam    ,fmtfil    ,ntof      ,mxnto     , & 
                    & kc        ,mxkc      ,omega     ,hydrbc    ,ascon     , & 
                    & gdp       ) 
           if (error) goto 999 
           ! 
           if (ascon /= 'Y') then 
              ! Write an unformatted intermediate file with 
              ! OMEGA, AMPLITUDES and PHASES 
              ! 
              lunbch = newlun(gdp) 
              filnam = 'TMP_' // runid(:lrid) // '.bch' 
              open (lunbch, file = filnam(:8 + lrid), form = 'unformatted',            & 
                   & status = 'unknown') 
              ! 
              ! Write file for DATBND values of NTOF are 'H' 
              ! 
              call wrtbch(lunbch    ,ntof      ,kc        ,omega     ,hydrbc    , & 
                        & mxnto     ,mxkc      ) 
              close (lunbch) 
           endif 
        endif 
        ! 
        ! QH boundary conditions 
        ! only if NTOQ > 0 
        ! 
        if (ntoq > 0) then 
           filnam = ' ' 
           ! 
           ! Read input file and write an unformatted intermediate file with 
           ! QH relations 
           ! 
           call rdbcq(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & runid     ,filnam    ,eol       ,nambnd    ,nto       , & 
                    & ntof      ,ntoq      ,bubble    ,gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Re-read MD-file  
        ! 
        rewind (lunmd) 
        read (lunmd, '(a300)') mdfrec 
        nrrec  = 1 
        ! 
        ! Boundary conditions Time series 
        ! only if NTOT > 0 
        ! 
        if (ntot > 0) then 
           filnam = ' ' 
           ntimtm = 0 
           ! 
           ntot0 = ntof + ntoq 
           call rdbct(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & noui      ,nrver     ,runid     ,filnam    ,eol       , & 
                    & nambnd    ,typbnd    ,tprofu    ,nto       ,ntot      , & 
                    & ntot0     ,kmax      ,rtime     ,itstrt    ,itfinish  , & 
                    & mxtime    ,ntimtm    ,rval      ,bubble    ,gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Boundary conditions Time series for the transported quantity 
        ! only if NTOT > 0 and LSTSC > 0 
        ! 
        if (nto>0 .and. lstsc>0) then 
           filnam = ' ' 
           ntimtm = 0 
           call rdbcc(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & noui      ,nrver     ,runid     ,filnam    ,eol       , & 
                    & nambnd    ,namcon    ,nto       ,lstsc     ,kmax      , & 
                    & rtime     ,itstrt    ,itfinish  ,mxtime    ,ntimtm    , & 
                    & salin     ,temp      ,const     ,lconc     ,rval      , & 
                    & rval      ,rval      ,rval      ,cval      ,bubble    , & 
                    & gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Discharge sources Time series 
        ! only if NSRC > 0 
        ! 
        if (nsrc > 0) then 
           filnam = ' ' 
           ntimtm = 0 
           call rddis(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & noui      ,nrver     ,runid     ,filnam    ,eol       , & 
                    & namsrc    ,disint    ,namcon    ,nsrcd     ,rtime     , & 
                    & itstrt    ,itfinish  ,mxnsrc    ,mxtime    ,ntimtm    , & 
                    & lstsc     ,salin     ,temp      ,const     ,lconc     , & 
                    & rval      ,rval      ,rval      ,rval      ,bubble    , & 
                    & gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Discharge sources Time series 
        ! only if NSLUV > 0 and RTC-coupling for Barrier heights to FLOW 
        ! 
        if (rtcmod == dataFromRTCToFLOW .and. nsluv>0) then 
           filnam = ' ' 
           ntimtm = 0 
           ! 
           call rdbcb(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & runid     ,filnam    ,itstrt    ,itfinish  ,gdp       ) 
           if (error) goto 999 
           ! 
           ! Create trigger file for starting RTC
           ! 
           luntri = newlun(gdp) 
           filsim = 'TMP_SYNC.RUN' 
           inquire (file = filsim, exist = lexist) 
           if (lexist) then 
              open (luntri, file = filsim) 
              close (luntri, status = 'delete') 
           endif 
           open (luntri, file = filsim, form = 'unformatted', status = 'unknown') 
           ! 
           ! Write 'RUNRTC' by telephone 
           ! 
           write (luntri) 786782 
           close (luntri) 
        endif 
        ! 
        ! Temperature Time series for HEAT modules 
        ! only if KTEMP > 0 
        ! 
        if (ktemp > 0) then 
           filnam = ' ' 
           fmtfil = 'FR' 
           ntimtm = 0 
           call rdheat(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec      , & 
                     & noui      ,runid     ,filnam    ,fmtfil    ,ktemp       , & 
                     & rval      ,dt        ,itstrt    ,itfinish  ,mxtime      , & 
                     & ntimtm    ,ivapop    ,rval      ,rval      ,rval        , & 
                     & rval      ,rval      ,rval      ,rval      ,solrad_read , & 
                     & gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! Rainfall / Evaporation Time series 
        ! only if KEVA > 0 
        ! 
        if (keva > 0) then 
           filnam = ' ' 
           fmtfil = 'FR' 
           ntimtm = 0 
           ! 
           call rdeva(lunmd     ,lundia    ,error     ,nrrec     ,mdfrec    , & 
                    & noui      ,runid     ,filnam    ,fmtfil    ,rval      , & 
                    & dt        ,itstrt    ,itfinish  ,mxtime    ,ntimtm    , & 
                    & rval      ,rval      ,rval      ,gdp       ) 
           if (error) goto 999 
        endif 
        ! 
        ! End of reading MD-file and attribute files or error found 
        ! 
        ! 
        ! Date and time 
        ! 
      999 continue 
        ! 
        deallocate(hydrbc) 
        deallocate(rval) 
        deallocate(cval) 
        deallocate(namcon) 

        deallocate(mnbnd )
        deallocate(mnksrc)
        deallocate(kspu  )
        deallocate(kspv  )
        deallocate(nhsub )
        deallocate(wstcof)
        deallocate(omega )
        deallocate(thick )
        deallocate(alpha )
        deallocate(rtime )
        deallocate(disint)
        deallocate(datbnd)
        deallocate(typbnd)
        deallocate(runtxt)
        deallocate(statns)
        deallocate(namsrc)
        deallocate(nambnd)
        deallocate(tprofu)
        ! 
     9990 continue 
        ! 
        call dattim(rundat    ) 
        date(1:4)  = rundat(1:4) 
        date(5:5)  = '-' 
        date(6:7)  = rundat(6:7) 
        date(8:8)  = '-' 
        date(9:10) = rundat(9:10) 
        write (lundia, *) 
        write (lundia, '(a,t77,a,t32,a)') '*** End   of tdatom for model:', '***',  & 
                                        & runid(:lrid) 
        write (lundia, '(a,t11,a,a,t77,a)') '*** ', date, rundat(11:19), '***' 
        ! 
        ! Calculate number of errors 
        ! 
        rewind (lundia) 
        ! 
        keyw  = '*** ERROR' 
        nrrec = 1 
        ntrec = 1 
        iend  = 20 
        lkw   = 9 
        merr  = 0 
        newkw = .true. 
        call search(lundia    ,error     ,newkw     ,nrrec     ,found     , & 
                  & ntrec     ,mdfrec    ,iend      ,keyw      ,lkw       , & 
                  & 'NOPR'    ) 
        if (found) then 
           merr = 1 
        endif 
        ! 
        ! Calculate number of warnings 
        ! 
        rewind (lundia) 
        ! 
        keyw  = '*** WARNI' 
        nrrec = 1 
        ntrec = 1 
        lkw   = 9 
        mwarn = 0 
        newkw = .true. 
        call search(lundia    ,error     ,newkw     ,nrrec     ,found     , & 
                  & ntrec     ,mdfrec    ,iend      ,keyw      ,lkw       , & 
                  & 'NOPR'    ) 
        if (found) then 
           mwarn = 1 
        endif 
        ! 
        ! Close files depending on number of errors and warnings 
        ! 
        if (merr+mwarn == 0) then 
           close (lundia, status = 'delete') 
        else 
           close (lundia) 
        endif 
    endif 
    ! 
    ! Close MD-file 
    ! 
    close (lunmd) 
    ! 
    if (error) nuerr = 1 
end subroutine tdatom 

