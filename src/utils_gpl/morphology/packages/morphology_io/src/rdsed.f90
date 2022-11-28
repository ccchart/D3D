module m_rdsed
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2022.                                
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
!-------------------------------------------------------------------------------
use m_depfil_stm

private

!
! functions and subroutines
!
public rdsed
public count_sed
public echosed

contains

subroutine rdsed(lundia    ,error     ,lsal      ,ltem      ,lsed      , &
               & lsedtot   ,lstsci    ,ltur      ,namcon    ,iopsus    , &
               & nmlb      ,nmub      ,kmax      ,filsed    ,sed_ptr   , &
               & slu_ptr   ,sedpar    ,trapar    ,griddim   )
!!--description-----------------------------------------------------------------
!
! Read sediment parameters from an input file
! File type:
!    ASCII-file if file version below 02.00
!    INI  -file if file version is 02.00 or higher
!
!!--declarations----------------------------------------------------------------
    use precision
    use properties
    use string_module
    use message_module
    use morphology_data_module
    use sediment_basics_module
    use system_utils, only:SHARED_LIB_PREFIX, SHARED_LIB_EXTENSION
    use grid_dimens_module, only: griddimtype
    !
    implicit none
    !
    ! The following list of pointer parameters is used to point inside sedpar and trapar
    !
    real(fp)                           , pointer :: csoil
    real(fp)                           , pointer :: mdcuni
    real(fp)                           , pointer :: kssilt
    real(fp)                           , pointer :: kssand
    real(fp)                           , pointer :: sc_cmf1
    real(fp)                           , pointer :: sc_cmf2
    real(fp)                           , pointer :: sc_flcf
    integer                            , pointer :: nmudfrac
    integer                            , pointer :: sc_mudfac
    integer                            , pointer :: rheo
    real(fp)         , dimension(:)    , pointer :: tpsnumber
    real(fp)         , dimension(:)    , pointer :: rhosol
    real(fp)         , dimension(:,:,:), pointer :: logseddia
    real(fp)         , dimension(:)    , pointer :: logsedsig
    real(fp)         , dimension(:)    , pointer :: sedd10
    real(fp)         , dimension(:)    , pointer :: sedd50
    real(fp)         , dimension(:)    , pointer :: sedd50fld
    real(fp)         , dimension(:)    , pointer :: seddm
    real(fp)         , dimension(:)    , pointer :: sedd90
    real(fp)         , dimension(:)    , pointer :: cdryb
    real(fp)         , dimension(:)    , pointer :: dstar
    real(fp)         , dimension(:)    , pointer :: taucr
    real(fp)         , dimension(:)    , pointer :: tetacr
    real(fp)         , dimension(:)    , pointer :: facdss
    real(fp)         , dimension(:)    , pointer :: sdbuni
    real(fp)         , dimension(:)    , pointer :: sedtrcfac
    real(fp)         , dimension(:)    , pointer :: thcmud
    real(fp)         , dimension(:)    , pointer :: tcguni
    real(fp)         , dimension(:)    , pointer :: mudcnt
    real(fp)         , dimension(:)    , pointer :: pmcrit
    real(fp)         , dimension(:,:)  , pointer :: rhocf
    real(fp)         , dimension(:,:)  , pointer :: cfvic
    real(fp)         , dimension(:,:)  , pointer :: cfmu
    real(fp)         , dimension(:,:)  , pointer :: xmu
    real(fp)         , dimension(:,:)  , pointer :: cfty
    real(fp)         , dimension(:,:)  , pointer :: cftau
    real(fp)         , dimension(:,:)  , pointer :: tyield
    real(fp)         , dimension(:,:)  , pointer :: taubh
    integer          , dimension(:)    , pointer :: nseddia
    integer          , dimension(:)    , pointer :: sedtyp
    character(10)    , dimension(:)    , pointer :: inisedunit
    character(20)    , dimension(:)    , pointer :: namsed
    character(256)   , dimension(:)    , pointer :: flsdbd
    character(256)   , dimension(:)    , pointer :: flstcg
    logical                            , pointer :: anymud
    logical                            , pointer :: bsskin
    character(256)                     , pointer :: flsdia
    character(256)                     , pointer :: flsmdc
    character(256)                     , pointer :: flspmc
    character(256)   , dimension(:)    , pointer :: dll_function_settle
    character(256)   , dimension(:)    , pointer :: dll_name_settle
    integer(pntrsize), dimension(:)    , pointer :: dll_handle_settle
    character(256)   , dimension(:)    , pointer :: dll_usrfil_settle
    integer          , dimension(:)    , pointer :: iform_settle
    real(fp)         , dimension(:,:)  , pointer :: par_settle
    integer          , dimension(:)    , pointer :: iform
    character(256)   , dimension(:)    , pointer :: flstrn
!
! Arguments
!
    integer                                  , intent(in)  :: lsal    !  Description and declaration in dimens.igs
    integer                                  , intent(in)  :: lsed    !  Description and declaration in esm_alloc_int.f90
    integer                                  , intent(in)  :: lsedtot !  Description and declaration in esm_alloc_int.f90
    integer                                  , intent(in)  :: lstsci  !  Description and declaration in esm_alloc_int.f90
    integer                                  , intent(in)  :: ltem    !  Description and declaration in dimens.igs
    integer                                  , intent(in)  :: ltur    !  Description and declaration in esm_alloc_int.f90
    integer                                  , intent(in)  :: lundia  !  Description and declaration in inout.igs
    logical                                  , intent(out) :: error   !!  Flag=TRUE if an error is encountered
    character(20) , dimension(lstsci + ltur) , intent(in)  :: namcon  !  Description and declaration in esm_alloc_char.f90
    integer                                  , intent(out) :: iopsus
    integer                                  , intent(in)  :: nmlb
    integer                                  , intent(in)  :: nmub
    integer                                  , intent(in)  :: kmax
    character(len=*)                         , intent(in)  :: filsed
    type(tree_data)                          , pointer     :: sed_ptr
    type(tree_data)                          , pointer     :: slu_ptr
    type(sedpar_type)                        , pointer     :: sedpar
    type(trapar_type)                        , pointer     :: trapar
    type(griddimtype)             , target   , intent(in)  :: griddim
!
! Local variables
!
    integer                     :: i
    integer                     :: iocond
    integer                     :: istat
    integer(pntrsize)           :: istat_ptr
    integer                     :: j
    integer                     :: l
    integer                     :: lbl                 ! bedload fraction number: lbl = l - lsed
    integer                     :: lenc                ! Help var. (length of character var.) 
    integer                     :: lfile
    integer                     :: luninp
    integer                     :: n                   ! Temporary storage for nseddia(l)
    integer                     :: nm
    integer                     :: version
    integer(pntrsize), external :: open_shared_library
    real(fp)                    :: rmissval
    real(fp)                    :: nclay
    real(fp)                    :: plasin
    real(fp)                    :: seddxx              ! Temporary storage for sediment diameter
    real(fp)                    :: sedsg               ! Temporary storage for geometric standard deviation
    real(fp)                    :: tpsmud
    logical                     :: ex
    logical                     :: success
    character(11)               :: fmttmp ! Format file ('formatted  ') 
    character(20)               :: sedname
    character(256)              :: filtrn
    character(256)              :: rec
    character(300)              :: message
    character(80)               :: parname
    character(20)               :: sc_type
    character(20)               :: sedtype             ! Local variable for sediment type
    character(78)               :: string
    character(10)               :: versionstring
    character(6)                :: seddxxstring
    character(256)              :: errmsg
    type(tree_data), pointer    :: sedblock_ptr
!
!! executable statements -------------------------------------------------------
!
    csoil                => sedpar%csoil
    mdcuni               => sedpar%mdcuni
    rheo                 => sedpar%rheologymodel
    kssilt               => sedpar%kssilt
    kssand               => sedpar%kssand
    sc_cmf1              => sedpar%sc_cmf1
    sc_cmf2              => sedpar%sc_cmf2
    sc_flcf              => sedpar%sc_flcf
    nmudfrac             => sedpar%nmudfrac
    sc_mudfac            => sedpar%sc_mudfac
    tpsnumber            => sedpar%tpsnumber
    rhosol               => sedpar%rhosol
    logseddia            => sedpar%logseddia
    logsedsig            => sedpar%logsedsig
    sedd10               => sedpar%sedd10
    sedd50               => sedpar%sedd50
    sedd50fld            => sedpar%sedd50fld
    seddm                => sedpar%seddm
    sedd90               => sedpar%sedd90
    cdryb                => sedpar%cdryb
    dstar                => sedpar%dstar
    taucr                => sedpar%taucr
    tetacr               => sedpar%tetacr
    sdbuni               => sedpar%sdbuni
    thcmud               => sedpar%thcmud
    mudcnt               => sedpar%mudcnt
    pmcrit               => sedpar%pmcrit
    nseddia              => sedpar%nseddia
    sedtyp               => sedpar%sedtyp
    inisedunit           => sedpar%inisedunit
    namsed               => sedpar%namsed
    flsdbd               => sedpar%flsdbd
    anymud               => sedpar%anymud
    bsskin               => sedpar%bsskin
    flsdia               => sedpar%flsdia
    flsmdc               => sedpar%flsmdc
    flspmc               => sedpar%flspmc
    dll_function_settle  => trapar%dll_function_settle
    dll_name_settle      => trapar%dll_name_settle
    dll_handle_settle    => trapar%dll_handle_settle
    dll_usrfil_settle    => trapar%dll_usrfil_settle
    iform_settle         => trapar%iform_settle
    par_settle           => trapar%par_settle
    iform                => trapar%iform
    flstrn               => trapar%flstrn
    !
    rmissval = -999.0_fp
    fmttmp   = 'formatted'
    !
    istat = 0
    !
    ! Read flag stressStrainRelation first:
    ! It influences the allocate statements
    !
    sedpar%stressStrainRelation = .false.
    call prop_get_logical(slu_ptr, 'Slurry', 'stressStrainRelation', sedpar%stressStrainRelation)
    !
    if (.not. associated(sedpar%sedd50)) then
       !
       ! allocation of namsed, rhosol and sedtyp have been allocated in count_sed routine
       !
       if (istat==0) allocate (sedpar%tpsnumber (                          lsedtot), stat = istat)
       !
       if (istat==0) allocate (sedpar%sedblock  (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%nseddia   (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%logseddia (2, 101,                   lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%logsedsig (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%sedd10    (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%sedd50    (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%seddm     (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%sedd90    (                          lsedtot), stat = istat)
       !
       if (istat==0) allocate (sedpar%cdryb     (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%dstar     (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%taucr     (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%tetacr    (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%sdbuni    (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%sedtrcfac (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%flsdbd    (                          lsedtot), stat = istat)
       if (istat==0) allocate (sedpar%inisedunit(                          lsedtot), stat = istat)
       !
       if (istat==0) allocate (sedpar%dss       (nmlb:nmub            ,max(1,lsed)), stat = istat)
       if (istat==0) allocate (sedpar%facdss    (                      max(1,lsed)), stat = istat)
       if (istat==0) allocate (sedpar%tcguni    (                      max(1,lsed)), stat = istat)
       !
       if (istat==0) allocate (sedpar%thcmud    (nmlb:nmub            ), stat = istat)
       if (istat==0) allocate (sedpar%flstcg    (                      max(1,lsed)), stat = istat)
       !
       if (istat==0) allocate (sedpar%mudcnt    (nmlb:nmub            ), stat = istat)
       if (istat==0) allocate (sedpar%pmcrit    (nmlb:nmub            ), stat = istat)
       if (istat==0) allocate (sedpar%sedd50fld (nmlb:nmub            ), stat = istat)
       !
       if (sedpar%stressStrainRelation) then
          if (istat==0) allocate (sedpar%phiclay   (nmlb:nmub            ,kmax), stat = istat)
          if (istat==0) allocate (sedpar%phisand   (nmlb:nmub            ,kmax), stat = istat)
       endif
       !
       ! Allocation of arrays rhocf and cfvic is after reading the value of parameter stressStrainRelation
       !
       if (istat/=0) then
          errmsg = 'RDSED: memory alloc error'
          call write_error(errmsg, unit=lundia)
          error = .true.
          return
       endif
       !
       ! update local pointers
       !
       tpsnumber     => sedpar%tpsnumber
       !
       nseddia       => sedpar%nseddia
       logseddia     => sedpar%logseddia
       logsedsig     => sedpar%logsedsig
       sedd10        => sedpar%sedd10
       sedd50        => sedpar%sedd50
       seddm         => sedpar%seddm
       sedd90        => sedpar%sedd90
       !
       cdryb         => sedpar%cdryb
       dstar         => sedpar%dstar
       taucr         => sedpar%taucr
       tetacr        => sedpar%tetacr
       sdbuni        => sedpar%sdbuni
       sedtrcfac     => sedpar%sedtrcfac
       flsdbd        => sedpar%flsdbd
       inisedunit    => sedpar%inisedunit
       !
       facdss        => sedpar%facdss
       thcmud        => sedpar%thcmud
       !
       mudcnt        => sedpar%mudcnt
       pmcrit        => sedpar%pmcrit
       sedd50fld     => sedpar%sedd50fld
       tcguni        => sedpar%tcguni
       flstcg        => sedpar%flstcg
       !
       !
       ! end check on assocation of sedpar%sedd50
       !
    endif 
    !
    ! Initialization of the just allocated arrays
    !
    do i = 1,lsedtot
       sedpar%sedblock(i)%node_name => null()
       if (sedtyp(i) == SEDTYP_COHESIVE) then
           tpsnumber(i) = 0.7_fp
       else
           tpsnumber(i) = 1.0_fp
       endif
    enddo
    flsdbd              = ' '
    flsmdc              = ' '
    flspmc              = ' '
    flsdia              = ' '
    dll_function_settle = ' '
    dll_usrfil_settle   = ' '
    flstcg              = ' '
    flstrn              = ' '
    !
    nseddia      = 0        ! nseddia counts relevant data
    logseddia    = rmissval
    logsedsig    = rmissval
    sedd10       = rmissval
    sedd50       = rmissval
    seddm        = rmissval
    sedd90       = rmissval
    !
    sedtrcfac    = rmissval
    !
    dstar        = rmissval
    taucr        = rmissval
    tetacr       = rmissval
    inisedunit   = 'kg/m2'
    !
    thcmud       = rmissval
    !
    mudcnt       =  0.0
    pmcrit       = -1.0
    sedd50fld    = rmissval
    !
    tcguni       = 1.5
    !
    ! Initialization of local parameters/arrays
    !
    version      = 0
    error        = .false.
    lenc         = 4
    !       
    facdss       = rmissval
    !
    seddxxstring = 'SedDXX'
    anymud       = .false.
    !
    ! Defaults fall velocity formulation of winterwerp 2004
    !
    sedpar%falflc   = .false.
    sedpar%flck3    =  0.8
    sedpar%flck4    =  0.0130
    sedpar%flcd0    = 10.0e-06
    sedpar%flcr     =  0.44
    sedpar%flcnf    = 2.15_fp
    sedpar%flwsmin  =  0.3e-3
    !
    ! Defaults Schelde erosion parameters formulation Winterwerp en
    !          van Kesteren 2004
    sedpar%eroschel = .false.
    sedpar%erosk1   =   1.0
    sedpar%erosk2   =   1.0
    sedpar%erosns   =   0.0
    sedpar%eroscv   =   3.0e-9
    sedpar%eroscu   = 400.0
    sedpar%erosd50  =  20.0e-6
    !
    ! Defaults consolidation
    !
    sedpar%cons_mud     = .false.
    sedpar%cons_kk      = 2.0e-14
    sedpar%cons_ksig    = 2.0e8
    sedpar%cons_eta     = 1.0e5
    sedpar%flcnf_cons   = 2.8
    !
    ! Defaults Bingham
    !
    sedpar%bin_cvisco    = 1.0e-3
    sedpar%bin_cnvisco   = 4.0
    sedpar%bin_cyield    = 5.0e8
    sedpar%bin_abingh    = 0.02
    sedpar%pow_bng_mix   = 2.35
    sedpar%pow_bng_silt  = 2.35
    sedpar%pow_rich_zaki = 3.1   
    !
    ! Defaults erosion interface
    !
    sedpar%ero_intfc     = .false.
    sedpar%nstress_intfc = 10.000
    sedpar%erosk1_int    =  1.0
    sedpar%erosk2_int    =  1.0
    !
    ! Slurry
    !
    sedpar%shearsettling = .false.
    sedpar%vmudToVicuv   = .false.
    sedpar%vicThresh     = 1.0e10_fp
    !
    ! Check version number of slu input file
    !
    versionstring = ' '
    call prop_get_string(slu_ptr, 'SlurryFileInformation', 'FileVersion', versionstring)
    if (versionstring /= '01.00') then
       errmsg = 'Unexpected version number of Slurry file. Expecting "01.00"'
       call write_error(errmsg, unit=lundia)
       error = .true.
       return
    endif
    !
    ! Check version number of sed input file
    !
    versionstring = ' '
    call prop_get_string(sed_ptr, 'SedimentFileInformation', 'FileVersion', versionstring)
    if (versionstring == '02.00' .or. versionstring == '03.00') then
       if (versionstring == '03.00') sedpar%version = 3.0_fp
       error  = .false.
       !
       csoil  = 1.0e4_fp
       call prop_get(sed_ptr, 'SedimentOverall', 'Cref', csoil)
       !
       tpsmud  = 0.7_fp
       call prop_get(sed_ptr, 'SedimentOverall', 'MudTPS', tpsmud)
       do i = 1,lsed
          if (sedtyp(i) == SEDTYP_COHESIVE) then
              tpsnumber(i) = tpsmud
          endif
       enddo
       !
       iopsus = 0
       call prop_get_integer(sed_ptr, 'SedimentOverall', 'IopSus', iopsus)
       !
       call prop_get_string(sed_ptr, 'SedimentOverall', 'MudCnt', flsmdc)
       !
       sedpar%flnrd(0) = ' '
       call prop_get_string(sed_ptr, 'SedimentOverall', 'NodeRelations', sedpar%flnrd(0))
       if (sedpar%flnrd(0) .ne. ' ') then
          call combinepaths(filsed, sedpar%flnrd(0))
       endif
       !
       !
       ! Intel 7.0 crashes on an inquire statement when file = ' '
       !
       if (flsmdc == ' ') then
          ex = .false.
       else
          call combinepaths(filsed, flsmdc)
          inquire (file = flsmdc, exist = ex)
       endif
       if (ex) then
          !
          ! Space varying data has been specified
          ! Use routine that also read the depth file to read the data
          !
          call depfil_stm(lundia    ,error     ,flsmdc    ,fmttmp    , &
                        & mudcnt    ,1         ,1         ,griddim   , errmsg)
          if (error) then
              call write_error(errmsg, unit=lundia)
              return
          endif
          do nm = 1, griddim%nmmax
             mudcnt(nm) = max(0.0_fp, min(mudcnt(nm), 1.0_fp))
          enddo
       else
          flsmdc = ' '
          mdcuni = 0.0_fp
          call prop_get(sed_ptr, 'SedimentOverall', 'MudCnt', mdcuni)
          !
          ! Uniform data has been specified
          !
          mudcnt = max(0.0_fp,min(mdcuni,1.0_fp))
       endif
       !
       if ( .not. associated(sed_ptr%child_nodes) ) then
          errmsg = 'Unable to read sediment information'
          call write_error(errmsg, unit=lundia)
          error = .true.
          return
       endif
       !
       ! Sand-mud interaction parameters
       !
       call prop_get_string(sed_ptr, 'SedimentOverall', 'PmCrit', flspmc)
       !
       ! Intel 7.0 crashes on an inquire statement when file = ' '
       !
       if (flspmc == ' ') then
          ex = .false.
       else
          call combinepaths(filsed, flspmc)
          inquire (file = flspmc, exist = ex)
       endif
       if (ex) then
          !
          ! Space varying data has been specified
          ! Use routine that also read the depth file to read the data
          !
          call depfil_stm(lundia    ,error     ,flspmc    ,fmttmp    , &
                        & pmcrit    ,1         ,1         ,griddim   )
          if (error) return
          do nm = 1, griddim%nmmax
             pmcrit(nm) = min(pmcrit(nm), 1.0_fp)
          enddo
       else
          flspmc = ' '
          call prop_get(sed_ptr, 'SedimentOverall', 'PmCrit', pmcrit(1))
          !
          ! Explicit loop because of stack overflow
          !
          do nm = 1, griddim%nmmax
             pmcrit(nm) = min(pmcrit(1), 1.0_fp)
          enddo
       endif
       !
       ! Get bed shear skin stress parameters
       !
       bsskin = .false.
       call prop_get_logical(sed_ptr, 'SedimentOverall', 'BsSkin', bsskin)
       if (bsskin) then
          call prop_get(sed_ptr, 'SedimentOverall', 'KsSilt', kssilt)
          call prop_get(sed_ptr, 'SedimentOverall', 'KsSand', kssand)
          !
          sc_type = 'thickness'
          call prop_get(sed_ptr, 'SedimentOverall', 'SC_mudfactor', sc_type)
          call str_lower(sc_type)
          select case (sc_type)
          case ('fraction')
             sc_mudfac = SC_MUDFRAC
          case ('thickness')
             sc_mudfac = SC_MUDTHC
          case default
             errmsg = 'Unknown option for sc_mudfactor. Expecting ''fraction'' or ''thickness''.'
             call write_error(errmsg, unit=lundia)
             error = .true.
             return
          end select
          !
          call prop_get(sed_ptr, 'SedimentOverall', 'SC_cmf1', sc_cmf1)
          call prop_get(sed_ptr, 'SedimentOverall', 'SC_cmf2', sc_cmf2)
          if (sc_mudfac == SC_MUDFRAC) then
             sc_cmf1 = max(0.0_fp , min(sc_cmf1, 1.0_fp))
             sc_cmf2 = max(sc_cmf1, min(sc_cmf2, 1.0_fp))
          else
             sc_cmf1 = max(0.0_fp , sc_cmf1)
             sc_cmf2 = max(sc_cmf1, sc_cmf2)
          endif
          !
          call prop_get(sed_ptr, 'SedimentOverall', 'CritFluffFactor', sc_flcf)
       endif
       !
       ! Get Fall velocity winterwerp parameters
       !
       call prop_get_logical(sed_ptr, 'SedimentOverall', 'FalFlc', sedpar%falflc)
       if (sedpar%falflc) then
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcNf_floc'   , sedpar%flcnf_floc )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcK2'        , sedpar%flck2 )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcK3'        , sedpar%flck3 )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcK4'        , sedpar%flck4 )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcD0'        , sedpar%flcd0 )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcR'         , sedpar%flcr  )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlWsmin'      , sedpar%flwsmin)
       endif
       !
       ! Schelde Erosion formulation parameters
       !
       call prop_get_logical(sed_ptr, 'SedimentOverall', 'EroSchel', sedpar%eroschel)
       if (sedpar%eroschel) then
          call prop_get(sed_ptr, 'SedimentOverall', 'EroSk1'    , sedpar%erosk1 )
          call prop_get(sed_ptr, 'SedimentOverall', 'EroSk2'    , sedpar%erosk2 )
          call prop_get(sed_ptr, 'SedimentOverall', 'EroSns'    , sedpar%erosns )
          call prop_get(sed_ptr, 'SedimentOverall', 'EroScv'    , sedpar%eroscv )
          call prop_get(sed_ptr, 'SedimentOverall', 'EroScu'    , sedpar%eroscu )
          call prop_get(sed_ptr, 'SedimentOverall', 'EroSd50'   , sedpar%erosd50)
          nclay     = (1.0 - sedpar%erosns) / 2.74
          plasin    = 254.0 * (nclay - 0.101)
       endif
       !
       ! Schelde consolidation parameters
       !
       call prop_get_logical(sed_ptr, 'SedimentOverall', 'cons_mud', sedpar%cons_mud)
       if (sedpar%cons_mud) then
          call prop_get(sed_ptr, 'SedimentOverall', 'cons_kk   '  , sedpar%cons_kk  )
          call prop_get(sed_ptr, 'SedimentOverall', 'cons_ksig '  , sedpar%cons_ksig)
          call prop_get(sed_ptr, 'SedimentOverall', 'cons_eta  '  , sedpar%cons_eta )
          call prop_get(sed_ptr, 'SedimentOverall', 'FlcNf_cons'  , sedpar%flcnf_cons )
       endif
       !
       ! Bingham model parameters
       !
       if (sedpar%stressStrainRelation) then
          call prop_get(sed_ptr, 'SedimentOverall', 'bin_cvisco'   , sedpar%bin_cvisco  )
          call prop_get(sed_ptr, 'SedimentOverall', 'bin_cnvisco'  , sedpar%bin_cnvisco )
          call prop_get(sed_ptr, 'SedimentOverall', 'bin_cyield'   , sedpar%bin_cyield  )
          call prop_get(sed_ptr, 'SedimentOverall', 'bin_abingh'   , sedpar%bin_abingh  )
          call prop_get(sed_ptr, 'SedimentOverall', 'pow_bng_mix'  , sedpar%pow_bng_mix )
          call prop_get(sed_ptr, 'SedimentOverall', 'pow_bng_silt' , sedpar%pow_bng_silt) 
          call prop_get(sed_ptr, 'SedimentOverall', 'pow_rich_zaki' , sedpar%pow_rich_zaki)          
          if (.not. sedpar%cons_mud) then
             call prop_get(sed_ptr, 'SedimentOverall', 'FlcNf_cons'   , sedpar%flcnf_cons )
          endif
       endif
       !
       ! Erosion of interface?
       !
       call prop_get_logical(sed_ptr, 'SedimentOverall', 'ero_intfc', sedpar%ero_intfc)
       if (sedpar%ero_intfc) then
          call prop_get(sed_ptr, 'SedimentOverall', 'nstress_intfc'   , sedpar%nstress_intfc  )
          call prop_get(sed_ptr, 'SedimentOverall', 'erosk1_int'      , sedpar%erosk1_int     )
          call prop_get(sed_ptr, 'SedimentOverall', 'erosk2_int'      , sedpar%erosk2_int     )
          if (.not. sedpar%eroschel) then
             call prop_get(sed_ptr, 'SedimentOverall', 'eroSns'       , sedpar%erosns         )
             call prop_get(sed_ptr, 'SedimentOverall', 'eroSd50'      , sedpar%erosd50        )
          endif
       endif
       !
       ! Some output parameters
       !
       allocate (sedpar%xmu   (nmlb:nmub,kmax), stat = istat)
       allocate (sedpar%tyield(nmlb:nmub,kmax), stat = istat)
       allocate (sedpar%taubh (nmlb:nmub,kmax), stat = istat)
       xmu           => sedpar%xmu
       tyield        => sedpar%tyield
       taubh         => sedpar%taubh
       xmu   (:,:) = -999.0
       tyield(:,:) = -999.0
       taubh (:,:) = -999.0
       !
       !
       if (sedpar%stressStrainRelation) then
          !
          ! Slurry
          !
          call prop_get_logical(slu_ptr, 'Slurry', 'shearSettling', sedpar%shearsettling)
          call prop_get        (slu_ptr, 'Slurry', 'vicThresh'    , sedpar%vicThresh)
          call prop_get_logical(slu_ptr, 'Slurry', 'vmudToVicuv'  , sedpar%vmudToVicuv)
          !
          ! Rheology (Jill Hanssen's work)
          !
          allocate (sedpar%rhocf (nmlb:nmub,kmax), stat = istat)
          allocate (sedpar%cfvic (nmlb:nmub,kmax), stat = istat)
          allocate (sedpar%cfmu  (nmlb:nmub,kmax), stat = istat)
          allocate (sedpar%cfty  (nmlb:nmub,kmax), stat = istat)
          allocate (sedpar%cftau (nmlb:nmub,kmax), stat = istat)
          rhocf         => sedpar%rhocf
          cfvic         => sedpar%cfvic
          cfmu          => sedpar%cfmu
          cfty          => sedpar%cfty
          cftau         => sedpar%cftau
          rhocf (:,:) = -999.0
          cfvic (:,:) = -999.0
          cfmu  (:,:) = -999.0
          cfty  (:,:) = -999.0
          cftau (:,:) = -999.0
          !
          ! Shearsettling: default: not set (rheo=-1)
          !
          sedpar%SluSettParam1 = 1.0_fp
          call prop_get(slu_ptr, 'Slurry', 'SluSettParam1', sedpar%SluSettParam1)
          sedpar%SluSettParam2 = 2.0_fp
          call prop_get(slu_ptr, 'Slurry', 'SluSettParam2', sedpar%SluSettParam2)
          sedpar%Shearsettle_w_opt = 1.0_fp  !% Arno's formula 1, Han's formula 2;
          call prop_get(slu_ptr, 'Slurry', 'Shearsettle_w_opt', sedpar%Shearsettle_w_opt)
          !
          ! Rheology: default: not set (rheo=-1)
          !
          string = 'dummy'
          call prop_get(slu_ptr, 'Slurry', 'rheology', string)
          call str_lower(string)
          if (string /= 'dummy') then
             if (string == 'winterwerp_kranenburg') then
                rheo = RHEOLOGY_WINTERWERP_KRANENBURG
             elseif (string == 'jacobs_vankesteren') then
                rheo = RHEOLOGY_JACOBS_VANKESTEREN
             elseif (string == 'thomas') then
                rheo = RHEOLOGY_THOMAS
             else
                errmsg = 'Value of keyword "rheology" must be equal to "Winterwerp_Kranenburg", "Jacobs_vanKesteren" or "Thomas"'
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             endif
          endif
          !
          ! Rheo= Winterwerp_Kranenburg: Read related parameters
          !
          if (rheo == RHEOLOGY_WINTERWERP_KRANENBURG) then
             sedpar%rheo_phisim = 0.6_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'phisim', sedpar%rheo_phisim)
             sedpar%rheo_ayield = 729884.0_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'ayield', sedpar%rheo_ayield)
             sedpar%rheo_frcdim = 2.6426_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'frcdim', sedpar%rheo_frcdim)
             sedpar%rheo_bety = 0.2752_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'bety', sedpar%rheo_bety)
             sedpar%rheo_avic = 931.86_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'avic', sedpar%rheo_avic)
             sedpar%rheo_powa = 3.65_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'powa', sedpar%rheo_powa)
             sedpar%rheo_betv = 0.2752_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'betv', sedpar%rheo_betv)
             sedpar%rheo_shrco = 5.0e6_fp
             call prop_get(slu_ptr, 'Winterwerp_Kranenburg', 'shrco', sedpar%rheo_shrco)
          endif
          !
          ! Rheo= Jacobs_vanKesteren: Read related parameters
          !
          if (rheo == RHEOLOGY_JACOBS_VANKESTEREN) then
             sedpar%rheo_phisim = 0.6_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'phisim', sedpar%rheo_phisim)
             sedpar%rheo_ayield = 4.167_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'ayield', sedpar%rheo_ayield)
             sedpar%rheo_powyie = -0.14_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'powyie', sedpar%rheo_powyie)
             sedpar%rheo_bety = 0.2752_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'bety', sedpar%rheo_bety)
             sedpar%rheo_avic = 0.043_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'avic', sedpar%rheo_avic)
             sedpar%rheo_powvic = -0.73_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'powvic', sedpar%rheo_powvic)
             sedpar%rheo_betv = 0.2752_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'betv', sedpar%rheo_betv)
             sedpar%rheo_shrco = 5000.0_fp
             call prop_get(slu_ptr, 'Jacobs_vanKesteren', 'shrco', sedpar%rheo_shrco)
          endif
          !
          ! Rheo= Thomas: Read related parameters
          !
          if (rheo == RHEOLOGY_THOMAS) then
             sedpar%rheo_phisim = 0.6_fp
             call prop_get(slu_ptr, 'Thomas', 'phisim', sedpar%rheo_phisim)
             sedpar%rheo_ayield = 7.45E5_fp
             call prop_get(slu_ptr, 'Thomas', 'ayield', sedpar%rheo_ayield)
             sedpar%rheo_powyie = 5.61_fp
             call prop_get(slu_ptr, 'Thomas', 'powyie', sedpar%rheo_powyie)
             sedpar%rheo_yieldk = 1.5_fp
             call prop_get(slu_ptr, 'Thomas', 'yieldk', sedpar%rheo_yieldk)
             sedpar%rheo_bvic = 17.7_fp
             call prop_get(slu_ptr, 'Thomas', 'bvic', sedpar%rheo_bvic)
             sedpar%rheo_visck = 1.25_fp
             call prop_get(slu_ptr, 'Thomas', 'visck', sedpar%rheo_visck)
             sedpar%rheo_shrco = 5000.0_fp
             call prop_get(slu_ptr, 'Thomas', 'shrco', sedpar%rheo_shrco)
          endif
       endif
       !
       do l = 1, lsedtot
          !
          ! namsed pre-filled in count-sed routine
          !
          sedname = namsed(l)
          write(lundia,'(a,i2,2a)') 'Sediment fraction ', l, ': ', sedname
          !
          ! find associated block (note that name order and block order may not correspond since
          ! suspended fractions need to be associated with the lower indices and e.g. Delft3D-FLOW
          ! requires suspended fraction order to match order of constituents).
          !
          do i = 1, size(sed_ptr%child_nodes)
             !
             ! Does sed_ptr contain a child with name 'Sediment' (converted to lower case)?
             !
             sedblock_ptr => sed_ptr%child_nodes(i)%node_ptr
             parname = tree_get_name( sedblock_ptr )
             call str_lower(parname)
             if ( trim(parname) /= 'sediment') cycle
             !
             parname = ' '
             call prop_get_string(sedblock_ptr, '*', 'Name', parname)
             if (.not. strcmpi(parname, sedname)) cycle
             !
             ! sediment fraction found
             !
             sedpar%flnrd(l) = ' '
             call prop_get_string(sedblock_ptr, '*', 'NodeRelations', sedpar%flnrd(l))
             if (sedpar%flnrd(l) .ne. ' ') then
                call combinepaths(filsed, sedpar%flnrd(l))
             endif
             !
             exit
          enddo
          sedpar%sedblock(l) = sedblock_ptr
          !
          rhosol(l) = rmissval
          call prop_get(sedblock_ptr, '*', 'RhoSol', rhosol(l))
          !
          ! Get the geometric standard deviation of the sediment fraction
          !
          sedsg = rmissval
          call prop_get(sedblock_ptr, '*', 'SedSg', sedsg)
          if (comparereal(sedsg,rmissval) /= 0) then
             logsedsig(l) = log(sedsg)
          endif
          !
          ! In case of one sediment fraction, it is possible to use a spatially
          ! varying grain size.
          !
          if (lsedtot == 1) then
             call prop_get_string(sedblock_ptr, '*', 'SedD50', flsdia)
             if (flsdia == ' ') then
                !
                ! Alternative for SedD50 is SedDia (backward compatibility)
                !
                call prop_get_string(sedblock_ptr, '*', 'SedDia', flsdia)
             endif
             !
             ! Intel 7.0 crashes on an inquire statement when file = ' '
             !
             if (flsdia == ' ') then
                ex = .false.
             else
                call combinepaths(filsed, flsdia)
                inquire (file = flsdia, exist = ex)
             endif
             if (ex) then
                !
                !  File with space varying data has been specified, read it now.
                !
                call depfil_stm(lundia    ,error     ,flsdia    ,fmttmp    , &
                              & sedd50fld ,1         ,1         ,griddim   , errmsg)
                if (error) then 
                    call write_error(errmsg, unit=lundia)
                    return
                endif      
             else
                flsdia = ' '
             endif
          else
             ex = .false.
          endif
          !
          ! If there are multiple sediment fractions, or if no spatially
          ! varying grain size file was specified, read the sediment size
          ! properties.
          !
          if (.not. ex) then
             do j = 0, 100
                seddxx = rmissval
                if (j == 0) then
                   call prop_get(sedblock_ptr, '*', 'SedMinDia', seddxx)
                elseif (j == 100) then
                   call prop_get(sedblock_ptr, '*', 'SedMaxDia', seddxx)
                else
                   write(seddxxstring(5:6),'(i2.2)') j
                   call prop_get(sedblock_ptr, '*', seddxxstring, seddxx)
                   if (j == 50 .and. comparereal(seddxx,rmissval) == 0) then
                      !
                      ! Alternative for SedD50 is SedDia (backward compatibility)
                      !
                      call prop_get(sedblock_ptr, '*', 'SedDia', seddxx)
                   endif
                endif
                if (comparereal(seddxx,rmissval) /= 0) then
                   if (seddxx <= 0.0_fp) then
                      !
                      ! error: sediment diameter less than or equal to zero!!
                      !
                      errmsg = 'Sediment diameters must be positive!'
                      call write_error(errmsg, unit=lundia)
                      error = .true.
                      return
                   endif
                   !
                   nseddia(l)       = nseddia(l) + 1
                   n                = nseddia(l)
                   logseddia(1,n,l) = real(j,fp)
                   logseddia(2,n,l) = log(seddxx)
                   if (n > 1) then
                      if (logseddia(2,n,l) <= logseddia(2,n-1,l)) then
                         !
                         ! error: sediment diameters not increasing!!
                         !
                         errmsg = 'Sediment diameters must be increasing!'
                         call write_error(errmsg, unit=lundia)
                         error = .true.
                         return
                      endif
                   endif
                endif
             enddo
          endif
          !
          if (l <= lsed) then
             rec = ' '
             call prop_get(sedblock_ptr, '*', 'SettleLib', rec)
             dll_name_settle(l) = rec
             if (rec /= ' ') then
                write(rec,'(3a)') SHARED_LIB_PREFIX, trim(dll_name_settle(l)), SHARED_LIB_EXTENSION
                dll_name_settle(l) = rec
                istat_ptr = 0
                istat_ptr = open_shared_library(dll_handle_settle(l), dll_name_settle(l))
                if (istat_ptr /= 0) then
                   write(errmsg,'(a,a)') 'Can not open shared library ', trim(dll_name_settle(l))
                   errmsg = FILE_NOT_FOUND // trim(errmsg)
                   call write_error(errmsg, unit=lundia)
                   error = .true.
                   return
                endif
                !
                call prop_get_string(sedblock_ptr, '*', 'SettleFunction', dll_function_settle(l))
                call prop_get_string(sedblock_ptr, '*', 'SettleInput'   , dll_usrfil_settle(l))
                iform_settle(l) = 15
             elseif (sedtyp(l) == SEDTYP_COHESIVE) then
                if (sedpar%falflc) then
                   iform_settle(l) = 3
                else
                   iform_settle(l) = 1
                endif
             elseif (sedtyp(l) == SEDTYP_NONCOHESIVE_SUSPENDED) then
                iform_settle(l) = 2
             endif
             !
             par_settle(:,l) = rmissval
             if (iform_settle(l) == 1) then
                call prop_get(sedblock_ptr, '*', 'SalMax', par_settle(1,l))
                call prop_get(sedblock_ptr, '*', 'WS0'   , par_settle(2,l))
                call prop_get(sedblock_ptr, '*', 'WSM'   , par_settle(3,l))
             elseif (iform_settle(l) == 2) then
                !
                ! These parameters will only be used for iform = -2, but unfortunately iform hasn't been determined yet.
                ! In the future we may have to read the parameters in a different order.
                !
                call prop_get(sedblock_ptr, '*', 'SalMax', par_settle(1,l))
                par_settle(2,l) = 1.0_fp
                call prop_get(sedblock_ptr, '*', 'GamFloc', par_settle(2,l))
             endif
             !
             ! Tracer calibration factor
             !
             call prop_get(sedblock_ptr, '*', 'TracerCalibrationFactor', sedtrcfac(l))
          endif
          !
          cdryb(l) = rmissval
          call prop_get(sedblock_ptr, '*', 'CDryB', cdryb(l))
          !
          ! First assume that 'IniSedThick'/'SdBUni' contains a filename
          ! If the file does not exist, assume that 'SdBUni' contains a uniform value (real)
          !
          call prop_get_string(sedblock_ptr, '*', 'IniSedThick', flsdbd(l))
          if (flsdbd(l) /= ' ') then
             inisedunit(l) = 'm'
          else
             inisedunit(l) = 'kg/m2'
             call prop_get_string(sedblock_ptr, '*', 'SdBUni', flsdbd(l))
          endif
          !
          ! Intel 7.0 crashes on an inquire statement when file = ' '
          !
          if (flsdbd(l) == ' ') then
             ex = .false.
          else
             call combinepaths(filsed, flsdbd(l))
             inquire (file = flsdbd(l), exist = ex)
          endif
          if (.not. ex) then
             sdbuni(l) = rmissval
             if (inisedunit(l) == 'm') then
                call prop_get(sedblock_ptr, '*', 'IniSedThick', sdbuni(l), success)
             else
                call prop_get(sedblock_ptr, '*', 'SdBUni', sdbuni(l), success)
             endif
             if (.not. success) then
                if (inisedunit(l) == 'm') then
                   errmsg = 'Error in IniSedThick: ' // trim(flsdbd(l)) // ' is not a file and not a value.'
                else
                   errmsg = 'Error in SdBUni.' // trim(flsdbd(l))
                endif
                errmsg = FILE_NOT_FOUND // trim(errmsg)
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             endif
             flsdbd(l) = ' '
          endif
          !
          if (l <= lsed) then
             if (sedpar%eroschel) then
                call write_error('eroschel code needs migration', unit=lundia)
                error = .true.
                return
                !erouni(l) = sedpar%erosk1*(sedpar%eroscv/sedpar%eroscu)* &
                !                  (cdryb(l)*cdryb(l)/rhosol(l)/(10.*sedpar%erosd50))
                !!
                !! Schelde formulations: Erosion = M*(taub-taucr)
                !!         for consistanct with general formulation (taub/taucr - 1)
                !!         M = Mschel*taucr
                !!
                !erouni(l) = erouni(l)*tceuni(l)
             endif
             call prop_get(sedblock_ptr, '*', 'FacDSS', facdss(l))
          endif
          !
          filtrn = ' '
          call prop_get(sedblock_ptr, '*', 'TraFrm', filtrn)
          !
          ! Intel 7.0 crashes on an inquire statement when file = ' '
          !
          if (filtrn /= ' ') then
             call combinepaths(filsed, filtrn)
             inquire (file = filtrn, exist = ex)
             if (ex) flstrn(l) = filtrn
             !
             call prop_get(sedblock_ptr, '*', 'TraFrm', iform(l))
          endif
       enddo
    else
       !
       ! sediment input version is 0 or 1:
       ! No keywords
       !
       call opensedfil(lundia    ,error     ,filsed    ,luninp    ,version  )
       if (error) return
       !write (versionstring, '(i4)') version
       !
       ! version 0 and 1 files don't contain sediment names.
       ! copy the names from the constituent array.
       !
       do l = 1, lsed
          namsed(l) = namcon(max(0, lsal, ltem) + l)
       enddo
       call rdsed01(lsed      ,luninp    ,lundia    ,csoil     ,iopsus    , &
                  & facdss    ,sedtyp    ,rhosol    ,sedd50    ,par_settle, &
                  & sdbuni    ,flsdbd    ,cdryb     ,sedpar%sedblock      , &
                  & version   ,error     )
       !
       ! nodal relations are not supported in older version sed files 
       !
       do l = 0, lsedtot
           sedpar%flnrd(l) = ' '
       enddo
       !
       close (luninp)
       if (error) return
       !
       do l = 1, lsed
          if (sedtyp(l) == SEDTYP_COHESIVE) then
             iform_settle(l) = 1
          elseif (sedtyp(l) == SEDTYP_NONCOHESIVE_SUSPENDED) then
             iform_settle(l) = 2
             par_settle(2,l) = 1.0_fp ! gamflc default
          endif
       enddo
       !
       ! Old format: SedDia specified, but stored in "wrong" array.
       !
       do l = 1, lsedtot
          if (sedd50(l) > 0.0_fp) then
             nseddia(l)       = 1
             logseddia(1,1,l) = 50.0_fp
             logseddia(2,1,l) = log(sedd50(l))
          endif
       enddo
    endif
    if (error) goto 9999
    do l = 1, lsedtot
       if (sedtyp(l) == SEDTYP_COHESIVE) then
          anymud   = .true.
          nmudfrac = nmudfrac + 1
       endif
    enddo
 9999 continue
end subroutine rdsed


subroutine rdsed01(lsed      ,luninp    ,lundia    ,csoil     ,iopsus    , &
                 & facdss    ,sedtyp    ,rhosol    ,seddia    ,parsettle , &
                 & sdbuni    ,flsdbd    ,cdryb     ,sedblock  ,version   , &
                 & error     )
!!--description-----------------------------------------------------------------
!
! Reads sediment input version 0 (or no version number found) or 1
!
!!--declarations----------------------------------------------------------------
    use precision
    use properties
    use sediment_basics_module
    use message_module, only: write_error
    !
    implicit none
!
! Arguments
!
    integer                        , intent(out):: iopsus
    integer                        , intent(in) :: lundia !  Description and declaration in inout.igs
    integer                        , intent(in) :: luninp
    integer                        , intent(in) :: version
    integer                                     :: lsed
    logical                        , intent(out):: error
    real(fp)                       , intent(out):: csoil
    real(fp)       , dimension(:)  , intent(out):: cdryb
    real(fp)       , dimension(:)  , intent(out):: facdss
    real(fp)       , dimension(:)  , intent(out):: rhosol
    real(fp)       , dimension(:,:), intent(out):: parsettle
    real(fp)       , dimension(:)  , intent(out):: sdbuni
    real(fp)       , dimension(:)  , intent(out):: seddia
    integer        , dimension(:)  , intent(out):: sedtyp !  sediment type: 0=total/1=noncoh/2=coh
    character(256) , dimension(:)  , intent(out):: flsdbd
    type(tree_data), dimension(:)  , intent(out):: sedblock
!
! Local variables
!
    integer                  :: iocond
    integer                  :: l
    integer                  :: lenc
    integer                  :: lsedlc
    character(4)             :: sedtype
    character(256)           :: errmsg
    character(256)           :: line
    type(tree_data), pointer :: sedblock_ptr
    type(tree_data), pointer :: anode
!
!! executable statements -------------------------------------------------------
!
    error  = .false.
    !
    if (version==0) then
       read (luninp, *, iostat = iocond) lsedlc, csoil
       iopsus = 0
    else
       read (luninp, *, iostat = iocond) lsedlc, csoil, iopsus
    endif
    !
    if (lsedlc /= lsed) then
       errmsg = 'Inconsistent number of sediment fractions'
       call write_error(errmsg, unit=lundia)
       error = .true.
       return
    endif
    !
    lenc = 4
    do l = 1, lsed
       call tree_create('Sediment',sedblock_ptr)
       !
       sedtype   = ' '
       read (luninp, '(a)', iostat = iocond) sedtype
       if (iocond == 0) then
          call small(sedtype, lenc)
          if (index(sedtype, 'sand') == 1) then
              sedtyp(l) = SEDTYP_NONCOHESIVE_SUSPENDED
          elseif (index(sedtype, 'mud') == 1) then
              sedtyp(l) = SEDTYP_COHESIVE
          else
             errmsg = 'Invalid suspended sediment type (must start with sand or mud)'
             call write_error(errmsg, unit=lundia)
             error = .true.
             return
          endif
       endif
       read (luninp, *, iostat = iocond) rhosol(l)
       if (iocond == 0) read (luninp, *, iostat = iocond) seddia(l)
       if (iocond == 0) read (luninp, *, iostat = iocond) parsettle(1,l) ! salmax
       if (iocond == 0) read (luninp, *, iostat = iocond) parsettle(2,l) ! ws0
       if (iocond == 0) read (luninp, *, iostat = iocond) parsettle(3,l) ! wsm
       if (iocond == 0) read (luninp, '(a)', iostat = iocond) line ! tcd
       call tree_create_node(sedblock_ptr,'TcrSed',anode)
       call tree_put_data(anode,transfer(trim(line),node_value),"STRING")
       if (iocond == 0) read (luninp, '(a)', iostat = iocond) line ! tce
       call tree_create_node(sedblock_ptr,'TcrEro',anode)
       call tree_put_data(anode,transfer(trim(line),node_value),"STRING")
       if (iocond == 0) read (luninp, '(a)', iostat = iocond) line ! ero
       call tree_create_node(sedblock_ptr,'EroPar',anode)
       call tree_put_data(anode,transfer(trim(line),node_value),"STRING")
       if (iocond == 0) read (luninp, *, iostat = iocond) cdryb(l)
       if (iocond == 0) read (luninp, *, iostat = iocond) sdbuni(l)
       if (iocond /= 0) then
          backspace (luninp)
          read (luninp, '(a)', iostat = iocond) flsdbd(l)
       endif
       if (version==0) then
          facdss(l) = 1.0_fp
       else
          if (iocond == 0) read (luninp, *, iostat = iocond) facdss(l)
       endif
       !
       sedblock(l) = sedblock_ptr
    enddo
    !
    deallocate(sedblock_ptr)
end subroutine rdsed01


subroutine opensedfil(lundia    ,error     ,filsed    ,luninp    ,version  )
!!--description-----------------------------------------------------------------
!
! Read sediment parameters from an input file
! File type:
!    ASCII-file if file version below 02.00
!    INI  -file if file version is 02.00 or higher
!
!!--declarations----------------------------------------------------------------
    use message_module
    use string_module
    implicit none
!
! Arguments
!
    integer                                  , intent(in)  :: lundia  ! unit of diagnostic file
    integer                                  , intent(out) :: luninp  ! unit of opened .sed file
    integer                                  , intent(out) :: version ! version of .sed file
    logical                                  , intent(out) :: error   ! .true. if error occured
    character(*)                                           :: filsed  ! name of .sed file
!
! Local variables
!
    integer                                                :: i
    integer                                                :: iocond
    character(256)                                         :: string
    character(256)                                         :: errmsg
!
!! executable statements -------------------------------------------------------
!
    error   = .false.
    version = 0
    !
    call remove_leading_spaces(filsed)
    open (newunit = luninp, file = filsed, form = 'formatted', status = 'old',            &
        & iostat = iocond)
    if (iocond /= 0) then
       errmsg = ERROR_FILE_OPEN // trim(filsed)
       call write_error(errmsg, unit=lundia)
       error = .true.
       return
    endif
    !
    ! Free formatted file, skip lines starting with a '*'
    !
    read (luninp, '(a)') string
    do while (string(:1) == '*')
       call str_lower(string)
       i = index(string, 'version')
       if (i /= 0) then
          read (string(i + 8:), '(i2)') version
       else
          write (lundia, '(a)') 'No version number found'
       endif
       read (luninp, '(a)') string
    enddo
    !
    rewind (luninp)
    call skipstarlines(luninp)
end subroutine opensedfil


subroutine echosed(lundia    ,error     ,lsed      ,lsedtot   , &
                 & iopsus    ,rhow      ,ag        ,sedpar    ,trapar    )
!!--description-----------------------------------------------------------------
!
! Report sediment parameter to diag file
!
!!--declarations----------------------------------------------------------------
    use precision
    use message_module
    use morphology_data_module
    use sediment_basics_module
    use m_rdtrafrm, only:echotrafrm
    !
    implicit none
!
! Arguments
!
    integer                                  , intent(in)  :: iopsus
    integer                                  , intent(in)  :: lsed    !  Description and declaration in iidim.f90
    integer                                  , intent(in)  :: lsedtot !  Description and declaration in iidim.f90
    integer                                                :: lundia  !  Description and declaration in inout.igs
    logical                                  , intent(out) :: error   !!  Flag=TRUE if an error is encountered
    type(sedpar_type)                        , pointer     :: sedpar
    type(trapar_type)                        , pointer     :: trapar
!
! Local variables
!
    real(fp)                          , pointer :: csoil
    real(fp)                          , pointer :: mdcuni
    real(fp)                          , pointer :: kssilt
    real(fp)                          , pointer :: kssand
    real(fp)                          , pointer :: sc_cmf1
    real(fp)                          , pointer :: sc_cmf2
    real(fp)                          , pointer :: sc_flcf
    integer                           , pointer :: sc_mudfac
    real(fp)        , dimension(:)    , pointer :: tpsnumber
    real(fp)        , dimension(:)    , pointer :: rhosol
    real(fp)        , dimension(:,:,:), pointer :: logseddia
    real(fp)        , dimension(:)    , pointer :: logsedsig
    real(fp)        , dimension(:)    , pointer :: sedd10
    real(fp)        , dimension(:)    , pointer :: sedd50
    real(fp)        , dimension(:)    , pointer :: sedd50fld
    real(fp)        , dimension(:)    , pointer :: seddm
    real(fp)        , dimension(:)    , pointer :: sedd90
    real(fp)        , dimension(:)    , pointer :: cdryb
    real(fp)        , dimension(:,:)  , pointer :: dss
    real(fp)        , dimension(:)    , pointer :: facdss
    real(fp)        , dimension(:)    , pointer :: sdbuni
    real(fp)        , dimension(:)    , pointer :: sedtrcfac
    real(fp)        , dimension(:)    , pointer :: tcguni
    real(fp)        , dimension(:)    , pointer :: pmcrit
    integer         , dimension(:)    , pointer :: nseddia
    integer         , dimension(:)    , pointer :: sedtyp
    character(10)   , dimension(:)    , pointer :: inisedunit
    character(20)   , dimension(:)    , pointer :: namsed
    character(256)  , dimension(:)    , pointer :: flsdbd
    character(256)  , dimension(:)    , pointer :: flstcg
    logical                           , pointer :: anymud
    logical                           , pointer :: bsskin
    character(256)                    , pointer :: flsdia
    character(256)                    , pointer :: flsmdc
    character(256)                    , pointer :: flspmc
    character(256)  , dimension(:)    , pointer :: dll_function_settle
    character(256)  , dimension(:)    , pointer :: dll_name_settle
    character(256)  , dimension(:)    , pointer :: dll_usrfil_settle
    integer         , dimension(:)    , pointer :: iform_settle
    real(fp)        , dimension(:,:)  , pointer :: par_settle
    integer         , dimension(:)    , pointer :: iform
    !
    integer                   :: i
    integer                   :: l
    integer                   :: n                   ! Temporary storage for nseddia(l)
    integer                   :: nm
    real(fp)                  :: logsedd50
    real(fp)                  :: rmissval
    real(fp)                  :: a_cu        ! coefficient
    real(fp)                  :: b_cu        ! exponent
    real(fp)                  :: taubin
    real(fp)                  :: watcon      ! water content
    real(fp)                  :: voidr       ! 1 + e/e with void ratio e
    real(fp)                  :: sedsg       ! Temporary storage for geometric standard deviation
    real(fp)                  :: xxinv               ! Help var. [1/xx or 1/(1-xx) in log unif distrib.]
    real(fp)                  :: nclay       ! Clay content
    real(fp)                  :: plasin      ! Plasticity index
    real(fp)                  :: xm
    logical        , external :: stringsequalinsens
    character(45)             :: txtput1
    character(10)             :: txtput2
    character(256)            :: errmsg
    real(fp), intent(in) :: rhow
    real(fp), intent(in) :: ag
!
!! executable statements -------------------------------------------------------
!
    csoil                => sedpar%csoil
    mdcuni               => sedpar%mdcuni
    kssilt               => sedpar%kssilt
    kssand               => sedpar%kssand
    sc_cmf1              => sedpar%sc_cmf1
    sc_cmf2              => sedpar%sc_cmf2
    sc_flcf              => sedpar%sc_flcf
    sc_mudfac            => sedpar%sc_mudfac
    tpsnumber            => sedpar%tpsnumber
    rhosol               => sedpar%rhosol
    logseddia            => sedpar%logseddia
    logsedsig            => sedpar%logsedsig
    sedd10               => sedpar%sedd10
    sedd50               => sedpar%sedd50
    sedd50fld            => sedpar%sedd50fld
    seddm                => sedpar%seddm
    sedd90               => sedpar%sedd90
    cdryb                => sedpar%cdryb
    dss                  => sedpar%dss
    facdss               => sedpar%facdss
    sdbuni               => sedpar%sdbuni
    sedtrcfac            => sedpar%sedtrcfac
    tcguni               => sedpar%tcguni
    pmcrit               => sedpar%pmcrit
    nseddia              => sedpar%nseddia
    sedtyp               => sedpar%sedtyp
    inisedunit           => sedpar%inisedunit
    namsed               => sedpar%namsed
    flsdbd               => sedpar%flsdbd
    flstcg               => sedpar%flstcg
    anymud               => sedpar%anymud
    bsskin               => sedpar%bsskin
    flsdia               => sedpar%flsdia
    flsmdc               => sedpar%flsmdc
    flspmc               => sedpar%flspmc
    dll_function_settle  => trapar%dll_function_settle
    dll_name_settle      => trapar%dll_name_settle
    dll_usrfil_settle    => trapar%dll_usrfil_settle
    iform_settle         => trapar%iform_settle
    par_settle           => trapar%par_settle
    iform                => trapar%iform
    !
    rmissval = -999.0_fp
    !
    ! echo input in diagnose-file
    !
    write (lundia, '(a)')   '*** Start  of sediment input'
    txtput1 = 'Ref concentration'
    write (lundia, '(2a,e12.4)') txtput1, ':', csoil
    if (csoil <= 0.0_fp) then
       errmsg = 'Reference concentration should be positive.'
       call write_error(errmsg, unit=lundia)
       error = .true.
       return
    endif
    txtput1 = 'Option Dss'
    write (lundia, '(2a,i12)') txtput1, ':', iopsus
    if (anymud) then
       if (flsmdc /= ' ' .or. comparereal(mdcuni,0.0_fp) /= 0) then
          errmsg = 'User defined mud content ignored: mud fraction simulated.'
          call write_warning(errmsg, unit=lundia)
       endif
       flsmdc = ' '
       mdcuni = 0.0_fp
       !
       ! Sand-mud interaction
       !
       txtput1 = 'Crit mud frac for sand mud interaction'
       if (flspmc /= ' ') then
          write (lundia, '(3a)') txtput1, ':  ', trim(flspmc)
       else
          write (lundia, '(2a,e12.4)') txtput1, ':', pmcrit(1)
       endif
    else
       if (flsmdc /= ' ') then
          txtput1 = 'File mud content'
          write (lundia, '(3a)') txtput1, ':  ', trim(flsmdc)
       else
          txtput1 = 'Uniform mud content'
          write (lundia, '(2a,e12.4)') txtput1, ':', mdcuni
       endif
       !
       ! Sand-mud interaction
       !
       if (flspmc /= ' ' .or. comparereal(pmcrit(1),0.0_fp) > 0) then
          errmsg = 'Sand mud interaction ignored: no mud fraction simulated.'
          call write_warning(errmsg, unit=lundia)
       endif
       flspmc = ' '
       pmcrit = -1.0_fp
    endif
    if (bsskin) then
       txtput1 = 'Skin friction Soulsby 2004'
       write (lundia, '(a)') txtput1
       !
       select case (sc_mudfac)
       case (SC_MUDFRAC)
          txtput2 = 'fraction'
       case (SC_MUDTHC)
          txtput2 = 'thickness'
       end select
       txtput1 = 'Lower crit mud '//txtput2
       write (lundia, '(2a,f12.6)') txtput1,':', sc_cmf1
       txtput1 = 'Upper crit mud '//txtput2
       write (lundia, '(2a,f12.6)') txtput1,':', sc_cmf2
       !
       txtput1 = 'Kssilt '
       write (lundia, '(2a,f12.6)') txtput1,':', kssilt
       txtput1 = 'Kssand '
       write (lundia, '(2a,f12.6)') txtput1,':', kssand
       !
       txtput1 = 'Critical fluff layer coverage factor'
       write (lundia, '(2a,f12.6)') txtput1,':', sc_flcf
    endif
    !
    if (sedpar%falflc) then
       txtput1 = 'Fall velocity Winterwerp 2004'
       write (lundia, '(a)') txtput1
       txtput1 = 'FlcNf_floc   '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flcnf_floc
       txtput1 = 'FlcK2   '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flck2
       txtput1 = 'FlcK3   '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flck3
       txtput1 = 'FlcK4   '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flck4
       txtput1 = 'FlcD0   '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flcd0
       txtput1 = 'FlcR    '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flcR
       txtput1 = 'FlWsmin '
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flwsmin
    endif
    !
    if (sedpar%eroschel) then
       txtput1 = 'New erosion formulation Schelde'
       write (lundia, '(a)') txtput1
       txtput1 = 'EroSk1  '
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%erosk1
       txtput1 = 'EroSk2  '
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%erosk2
       txtput1 = 'EroSns  '
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%erosns
       txtput1 = 'EroScu  '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%eroscu
       txtput1 = 'EroScv  '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%eroscv
       txtput1 = 'EroSd50 '
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%erosd50
    endif

    if (sedpar%cons_mud) then
       txtput1 = 'Consolidation of mud included'
       write (lundia, '(a)') txtput1
       txtput1 = 'Cons_kk   '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%cons_kk
       txtput1 = 'Cons_ksig '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%cons_ksig
       txtput1 = 'Cons_eta  '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%cons_eta
       txtput1 = 'FlcNf_cons'
       write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flcnf_cons
    endif

    if (sedpar%stressStrainRelation) then
       txtput1 = 'Non-Newtonian rheological model used'
       write (lundia, '(a)') txtput1
       txtput1 = 'bin_cvisco'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%bin_cvisco
       txtput1 = 'bin_cnvisco'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%bin_cnvisco
       txtput1 = 'bin_cyield'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%bin_cyield
       txtput1 = 'bin_abingh'
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%bin_abingh
       txtput1 = 'pow_bng_mix'
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%pow_bng_mix
       txtput1 = 'pow_bng_silt'
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%pow_bng_silt
       txtput1 = 'pow_rich_zaki'
       write (lundia, '(a,a,f14.6)') txtput1,':',sedpar%pow_rich_zaki       
       if (.not. sedpar%cons_mud) then
          txtput1 = 'FlcNf_cons'
          write (lundia, '(a,a,f12.6)') txtput1,':',sedpar%flcnf_cons
       endif
    endif
    !
    sedpar%power   = (3. - sedpar%flcnf_cons)/ 2.
    if (sedpar%ero_intfc) then
       !
       !  Inversion of normal stress criterion for erosion
       !
       sedpar%sedcint = rhosol(1)*(sedpar%nstress_intfc/sedpar%cons_ksig)**sedpar%power ! mass concentration g/l
       !
       !  Undrained shear strength C_U related to water content
       !
       a_cu    = 15.468
       b_cu    = -3.24
       watcon  = rhow/sedpar%sedcint - rhow / rhosol(1)
       nclay     = (1. - sedpar%erosns)/2.74
       if (watcon <= 0.0 ) then
         write(lundia,'(a)') 'Negative water content, check RHOSOL and CDRYB '
       else
         sedpar%eroscu  = a_cu * (watcon / nclay)**b_cu
       endif
       !
       ! Critical shear stress from erosion formulations for the
       ! Scheldt river
       !
       plasin    = 254.*(nclay - .101)
       !sedpar%tcrint(1) = 0.163*sedpar%erosk2_int*plasin**0.84
       !
       ! Consolidation coefficient CV
       !
       voidr          = 1./(1 - sedpar%sedcint/rhosol(1)) ! volume concentration needed
       sedpar%eroscv  = voidr * sedpar%cons_kk * sedpar%cons_ksig/(ag*rhow*sedpar%power)
       ! Parameter for internal erosion
       sedpar%xme_int = sedpar%erosk1_int*(sedpar%eroscv/sedpar%eroscu)* &
                    (sedpar%sedcint*sedpar%sedcint/rhosol(1)/(10.*sedpar%erosd50))
       !
       txtput1 = 'Erosion of interface'
       write (lundia, '(a)') txtput1
       txtput1 = 'nstress_intfc'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%nstress_intfc
       txtput1 = 'erosk1_int   '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%erosk1_int
       txtput1 = 'erosk2_int   '
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%erosk2_int
       if (.not.sedpar%eroschel) then
          txtput1 = 'erosns       '
          write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%erosns
          txtput1 = 'erosd50      '
          write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%erosd50
       endif

       txtput1 = 'sedcint'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%sedcint
       txtput1 = 'watcon'
       write (lundia, '(a,a,e14.6)') txtput1,':',watcon
       txtput1 = 'eroscu'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%eroscu
       txtput1 = 'eroscv'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%eroscv
       txtput1 = 'xme_int'
       write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%xme_int
    endif
    !
    if (sedpar%flnrd(0) /= ' ') then
       txtput1 = '1D nodal relations for bed/total load'
       write (lundia, '(3a)') txtput1, ':  ', trim(sedpar%flnrd(0))
    endif
    !
    do l = 1, lsedtot
       txtput1 = 'Sediment number'
       write (lundia, '(2a,i12)') txtput1, ':', l
       txtput1 = '  Name'
       write (lundia, '(3a)') txtput1, ': ', trim(namsed(l))
       txtput1 = '  Type'
       select case (sedtyp(l))
          case (SEDTYP_NONCOHESIVE_TOTALLOAD)
             write (lundia, '(2a,a12)') txtput1, ':', 'bedload'
          case (SEDTYP_NONCOHESIVE_SUSPENDED)
             write (lundia, '(2a,a12)') txtput1, ':', 'sand'
          case (SEDTYP_COHESIVE)
             write (lundia, '(2a,a12)') txtput1, ':', 'mud'
       end select
       if (sedtrcfac(l)>0.0_fp) then
           txtput1 = '  Tracer calibration factor '
           write (lundia, '(2a,e12.4)') txtput1, ':', sedtrcfac(l)
       endif
       if (l <= lsed) then
          txtput1 = '  Turbulent Prandtl-Schmidt number'
          write (lundia, '(2a,e12.4)') txtput1, ':', tpsnumber(l)
       endif
       txtput1 = '  Solid density (RHOSOL)'
       write (lundia, '(2a,e12.4)') txtput1, ':', rhosol(l)
       if (flsdia /= ' ') then
          !
          ! One sediment diameter between 0 and 100%
          !
          write (lundia, '(3a)') txtput1, ':  ', 'lognormal'
          if (comparereal(logsedsig(l),rmissval) == 0) then
             !
             ! no standard deviation specified: use default geometric
             ! standard deviation of 1.34 which is in the middle of the
             ! range of sigma_g values representing "well sorted" mixtures
             ! indicated in Blott & Pye, 2001. Earth Surface Processes
             ! and Landforms 26, p. 1237-1248.
             !
             logsedsig(l) = log(1.34)
          endif
          !
          nseddia(l) = -999
          txtput1 = '  geom. st. dev.'
          write (lundia, '(2a,e12.4)') txtput1, ':', exp(logsedsig(l))
          txtput1 = '  SedD50'
          write (lundia, '(3a)') txtput1, ':  ', trim(flsdia)
       elseif (sedtyp(l) /= SEDTYP_COHESIVE) then
          !
          ! Determine various sediment diameters in case of
          ! sand or bedload.
          !
          txtput1 = '  sed. distribution'
          if (nseddia(l) == 0) then
             !
             ! error: no sediment diameter specified!
             !
             errmsg = 'Missing sediment diameter data'
             call write_error(errmsg, unit=lundia)
             error = .true.
             return
          elseif (nseddia(l) == 1) then
             !
             ! Just one sediment diameter
             !
             if (nint(logseddia(1,1,l)) == 0) then
                !
                ! error: only minimum sediment diameter is insufficient
                !
                errmsg = 'Missing maximum diameter data'
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             elseif (nint(logseddia(1,1,l)) == 100) then
                !
                ! error: only maximum sediment diameter is insufficient
                !
                errmsg = 'Missing minimum diameter data'
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             endif
             !
             ! One sediment diameter between 0 and 100%
             !
             if (comparereal(sedpar%version,3.0_fp) >= 0) then
                !
                ! New behaviour: lognormal distribution
                !
                write (lundia, '(3a)') txtput1, ':  ', 'lognormal'
                if (comparereal(logsedsig(l),rmissval) == 0) then
                   !
                   ! no standard deviation specified: use default geometric
                   ! standard deviation of 1.34 which is in the middle of the
                   ! range of sigma_g values representing "well sorted" mixtures
                   ! indicated in Blott & Pye, 2001. Earth Surface Processes
                   ! and Landforms 26, p. 1237-1248.
                   !
                   logsedsig(l) = log(1.34)
                endif
                txtput1 = '  geom. st. dev.'
                write (lundia, '(2a,e12.4)') txtput1, ':', exp(logsedsig(l))
                !
                ! Approximate lognormal distribution using the following
                ! percentiles:
                ! 0.1 (set to 0) --- ilognormal --- 99.9 (set to 100)
                !
                logsedd50        = logseddia(2,1,l) - logsedsig(l)*lognormal(nint(logseddia(1,1,l)))
                logseddia(1,1,l) = 0.0_fp
                logseddia(2,1,l) = logsedd50 - 3.0_fp*logsedsig(l)
                n = 1
                do i = 1, size(ilognormal)
                   n = n + 1
                   logseddia(1,n,l) = real(ilognormal(i),fp)
                   logseddia(2,n,l) = logsedd50 + lognormal(ilognormal(i))*logsedsig(l)
                enddo
                n = n + 1
                logseddia(1,n,l) = 100.0_fp
                logseddia(2,n,l) = logsedd50 + 3.0_fp*logsedsig(l)
                nseddia(l)       = n
                !
                ! Compute characteristic sediment diameters
                !
                sedd50(l) = exp(logsedd50)
                sedd10(l) = exp(logsedd50 + lognormal(10)*logsedsig(l))
                sedd90(l) = exp(logsedd50 + lognormal(90)*logsedsig(l))
                seddm(l)  = exp(logsedd50 + 0.5_fp*logsedsig(l)*logsedsig(l))
             else
                !
                ! Old behaviour: D10 = 0.75 * D50, D90 = 1.5 * D50
                ! Piecewise loguniform approach: requires D50
                !
                if (nint(logseddia(1,1,l)) /= 50) then
                   !
                   ! error: old approach requires D50
                   !
                   errmsg = 'Missing median diameter data'
                   call write_error(errmsg, unit=lundia)
                   error = .true.
                   return
                endif
                !
                write (lundia, '(3a)') txtput1, ':  ', 'piecewise loguniform'
                !
                nseddia(l)       = 3
                logseddia(1,2,l) = 50.0_fp
                logseddia(2,2,l) = logseddia(2,1,l)
                !
                logseddia(1,1,l) = 0.0_fp
                logseddia(2,1,l) = logseddia(2,2,l) + 5.0_fp*log(0.75_fp)/4.0_fp
                !
                logseddia(1,3,l) = 100.0_fp
                logseddia(2,3,l) = logseddia(2,2,l) + 5.0_fp*log(1.5_fp)/4.0_fp
                !
                ! Compute characteristic sediment diameters
                !
                sedd50(l) = exp(logseddia(2,2,l))
                sedd10(l) = 0.75_fp * sedd50(l)
                sedd90(l) = 1.50_fp * sedd50(l)
                !
                seddm(l)  = 0.0_fp
                logsedsig(l) = 0.0_fp
                xm = 0.0_fp
                do n = 2, nseddia(l)
                   xxinv    = logseddia(1,n,l) - logseddia(1,n-1,l)
                   seddm(l) = seddm(l) &
                            & + xxinv * (exp(logseddia(2,n,l)) - exp(logseddia(2,n-1,l))) &
                            & / (logseddia(2,n,l) - logseddia(2,n-1,l))
                   logsedsig(l) = logsedsig(l) + xxinv * (&
                                & (logseddia(2,n,l)-logseddia(2,n-1,l))/sqrt(3.0_fp)/2.0_fp &
                                & + ((logseddia(2,n,l)+logseddia(2,n-1,l))/2.0_fp)**2 &
                                & )
                   xm = xm + xxinv * (logseddia(2,n,l)+logseddia(2,n-1,l))/2.0_fp
                enddo
                seddm(l) = seddm(l) / 100.0_fp
                xm = xm / 100.0_fp
                logsedsig(l) = logsedsig(l) / 100.0_fp - xm**2
             endif
          elseif (nseddia(l) == 2) then
             !
             ! Two sediment diameters specified
             !
             if (comparereal(logsedsig(l),rmissval) /= 0) then
                !
                ! standard deviation specified
                ! error: not allowed in combination with multiple
                ! sediment diameters
                !
                errmsg = 'Geom. st. dev. not allowed in combination with multiple sediment diameters'
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             endif
             if (nint(logseddia(1,1,l))==0 .or. nint(logseddia(1,2,l))==100) then
                !
                ! Minimum or maximum sediment diameter given:
                ! loguniform distribution
                !
                write (lundia, '(3a)') txtput1, ':  ', 'loguniform'
                !
                ! Compute characteristic sediment diameters
                !
                if (nint(logseddia(1,1,l)) > 0) then
                   !
                   ! Only maximum sediment diameter given and some
                   ! percentile: compute minimum sediment diameter.
                   !
                   xxinv            = 1.0_fp / (1.0_fp - (real(logseddia(1,1,l),fp)/100.0_fp))
                   logseddia(2,1,l) = logseddia(2,1,l)*xxinv + logseddia(2,2,l)*(1.0_fp-xxinv)
                   logseddia(1,1,l) = 0.0_fp
                elseif (nint(logseddia(1,2,l)) < 100) then
                   !
                   ! Only minimum sediment diameter given and some
                   ! percentile: compute maximum sediment diameter.
                   !
                   xxinv            = 100.0_fp / real(logseddia(1,2,l),fp)
                   logseddia(2,2,l) = logseddia(2,2,l)*xxinv + logseddia(2,1,l)*(1.0_fp-xxinv)
                   logseddia(1,2,l) = 100.0_fp
                endif
                !
                ! Both minimum and maximum sediment diameters given
                !
                sedd50(l) = exp(0.5_fp*logseddia(2,1,l) + 0.5_fp*logseddia(2,2,l))
                sedd10(l) = exp(0.9_fp*logseddia(2,1,l) + 0.1_fp*logseddia(2,2,l))
                sedd90(l) = exp(0.1_fp*logseddia(2,1,l) + 0.9_fp*logseddia(2,2,l))
                seddm(l)  = (exp(logseddia(2,2,l)) - exp(logseddia(2,1,l))) / &
                          & (logseddia(2,2,l) - logseddia(2,1,l))
                logsedsig(l) = (logseddia(2,2,l) - logseddia(2,1,l))/sqrt(3.0_fp)/2.0_fp
             else
                !
                ! Neither minimum nor maximum sediment diameter given:
                ! lognormal distribution
                !
                write (lundia, '(3a)') txtput1, ':  ', 'lognormal'
                !
                ! Compute geometric standard deviation
                !
                logsedsig(l) = (logseddia(2,2,l) - logseddia(2,1,l))/ &
                             & (lognormal(nint(logseddia(1,2,l))) -   &
                             &  lognormal(nint(logseddia(1,1,l))))
                txtput1 = ' geom. stand. dev.'
                write (lundia, '(2a,e12.4)') txtput1, ':', exp(logsedsig(l))
                !
                ! Approximate lognormal distribution using the following
                ! percentiles:
                ! 0.1 (set to 0) --- ilognormal --- 99.9 (set to 100)
                !
                logsedd50        = logseddia(2,1,l) - logsedsig(l)*lognormal(nint(logseddia(1,1,l)))
                logseddia(1,1,l) = 0.0_fp
                logseddia(2,1,l) = logsedd50 - 3.0_fp*logsedsig(l)
                n = 1
                do i = 1, size(ilognormal)
                   n = n + 1
                   logseddia(1,n,l) = real(ilognormal(i),fp)
                   logseddia(2,n,l) = logsedd50 + lognormal(ilognormal(i))*logsedsig(l)
                enddo
                n = n + 1
                logseddia(1,n,l) = 100.0_fp
                logseddia(2,n,l) = logsedd50 + 3.0_fp*logsedsig(l)
                nseddia(l)       = n
                !
                ! Compute characteristic sediment diameters
                !
                sedd50(l) = exp(logsedd50)
                sedd10(l) = exp(logsedd50 + lognormal(10)*logsedsig(l))
                sedd90(l) = exp(logsedd50 + lognormal(90)*logsedsig(l))
                seddm(l)  = exp(logsedd50 + 0.5_fp*logsedsig(l)*logsedsig(l))
             endif
          else
             !
             ! More than two sediment diameters specified
             !
             if (comparereal(logsedsig(l),rmissval) /= 0) then
                !
                ! standard deviation specified
                ! error: not allowed in combination with multiple
                ! sediment diameters
                !
                errmsg = 'Geom. std. dev. not allowed in combination with multiple sediment diameters'
                call write_error(errmsg, unit=lundia)
                error = .true.
                return
             endif
             !
             ! Always piecewise loguniform distribution
             !
             write (lundia, '(3a)') txtput1, ':  ', 'piecewise loguniform'
             !
             ! Compute characteristic sediment diameters
             !
             if (nint(logseddia(1,1,l)) > 0) then
                !
                ! sediment diameter table does not include sedmindia
                ! replace the first entry by extending the first loguniform range
                ! by keeping the density constant
                !
                xxinv            =  real(logseddia(1,2,l),fp) / &
                                 & (real(logseddia(1,2,l),fp) - real(logseddia(1,1,l),fp))
                logseddia(2,1,l) = logseddia(2,1,l)*xxinv + logseddia(2,2,l)*(1.0_fp-xxinv)
                logseddia(1,1,l) = 0.0_fp
             endif
             n = nseddia(l)
             if (nint(logseddia(1,n,l)) < 100) then
                !
                ! sediment diameter table does not include sedmaxdia
                ! replace the last entry by extending the last loguniform range
                ! by keeping the density constant
                !
                xxinv            = (100.0_fp - real(logseddia(1,n-1,l),fp)) / &
                                 & (real(logseddia(1,n,l),fp) - real(logseddia(1,n-1,l),fp))
                logseddia(2,n,l) = logseddia(2,n,l)*xxinv + logseddia(2,n-1,l)*(1.0_fp-xxinv)
                logseddia(1,n,l) = 100.0_fp
             endif
             !
             seddm(l)  = 0.0_fp
             do n = 2, nseddia(l)
                xxinv    = logseddia(1,n,l) - logseddia(1,n-1,l)
                seddm(l) = seddm(l) &
                         & + xxinv * (exp(logseddia(2,n,l)) - exp(logseddia(2,n-1,l))) &
                         & / (logseddia(2,n,l) - logseddia(2,n-1,l))
                xxinv    = 1.0_fp / xxinv
                if (logseddia(1,n-1,l) < 10.0_fp .and. logseddia(1,n,l) >= 10.0_fp) then
                   sedd10(l) = exp(logseddia(2,n-1,l)*((logseddia(1,n,l)-10.0_fp)*xxinv) + &
                                 & logseddia(2,n,l)*((10.0_fp-logseddia(1,n-1,l))*xxinv))
                endif
                if (logseddia(1,n-1,l) < 50.0_fp .and. logseddia(1,n,l) >= 50.0_fp) then
                   sedd50(l) = exp(logseddia(2,n-1,l)*((logseddia(1,n,l)-50.0_fp)*xxinv) + &
                                 & logseddia(2,n,l)*((50.0_fp-logseddia(1,n-1,l))*xxinv))
                endif
                if (logseddia(1,n-1,l) < 90.0_fp .and. logseddia(1,n,l) >= 90.0_fp) then
                   sedd90(l) = exp(logseddia(2,n-1,l)*((logseddia(1,n,l)-90.0_fp)*xxinv) + &
                                 & logseddia(2,n,l)*((90.0_fp-logseddia(1,n-1,l))*xxinv))
                endif
             enddo
             seddm(l) = seddm(l) / 100.0_fp
          endif
          !
          ! convert percentages to fractions
          !
          do n = 1, nseddia(l)
             logseddia(1,n,l) = logseddia(1,n,l) / 100.0_fp
          enddo
          !
          txtput1 = '  SedD10'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedd10(l)
          txtput1 = '  SedD50'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedd50(l)
          txtput1 = '  SedDM'
          write (lundia, '(2a,e12.4)') txtput1, ':', seddm(l)
          txtput1 = '  SedD90'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedd90(l)
       endif
       txtput1 = '  Dry bed (bulk) density (CDRYB)'
       write (lundia, '(2a,e12.4)') txtput1, ':', cdryb(l)
       if (flsdbd(l) /= ' ') then
          if (inisedunit(l) == 'kg/m2') then
             txtput1 = '  File IniCon'
          else
             txtput1 = '  File IniThick'
          endif
          write (lundia, '(3a)') txtput1, ':  ', trim(flsdbd(l))
       else
          if (inisedunit(l) == 'kg/m2') then
             txtput1 = '  Uniform IniCon'
             txtput2 = ' [kg/m2]'
          else
             txtput1 = '  Uniform IniThick'
             txtput2 = ' [m]'
          endif
          write (lundia, '(2a,e12.4,a)') txtput1, ':', sdbuni(l), trim(txtput2)
       endif
       if (l <= lsed) then
          txtput1 = '  FACDSS'
          write (lundia, '(2a,e12.4)') txtput1, ':', facdss(l)
          if (facdss(l) <= 0.0_fp) then
             errmsg = 'FACDSS <= 0.0'
             call write_error(errmsg, unit=lundia)
             error = .true.
             return
          endif
          !
          if (lsedtot==1 .and. flsdia/=' ') then ! l=1
             !
             ! Wrong: do nm=1,size(sedd50fld)
             ! nm goes from nmlb (negative!) to nmub
             !
             do nm = lbound(sedd50fld,1), ubound(sedd50fld,1)
                dss(nm, 1) = sedd50fld(nm)*facdss(1)
             enddo
          else
             !
             ! Explicit loop because of stack overflow
             !
             do nm = lbound(sedd50fld,1), ubound(sedd50fld,1)
                dss(nm, l) = sedd50(l)*facdss(l)
             enddo
          endif
       endif
       !
       call echotrafrm(lundia      ,trapar     ,l         )
       !
       if (iform_settle(l) == 1) then
          txtput1 = '  SALMAX'
          write (lundia, '(2a,e12.4)') txtput1, ':', par_settle(1,l)
          txtput1 = '  WS0'
          write (lundia, '(2a,e12.4)') txtput1, ':', par_settle(2,l)
          txtput1 = '  WSM'
          write (lundia, '(2a,e12.4)') txtput1, ':', par_settle(3,l)
       elseif (iform_settle(l) == 2) then
          if (iform(l) == -2 .or. iform(l) == -4) then
             iform_settle(l) = -2
             txtput1 = '  SALMAX'
             write (lundia, '(2a,e12.4)') txtput1, ':', par_settle(1,l)
             txtput1 = '  Flocculation factor GamFloc'
             write (lundia, '(2a,e12.4)') txtput1, ':', par_settle(2,l)
          endif
       elseif (iform_settle(l) == 15) then
          !
          ! User defined settling velocity function
          !
          txtput1 = '  Settle library'
          write (lundia, '(3a)') txtput1, ': ', trim(dll_name_settle(l))
          txtput1 = '  Function in Settle lib'
          write (lundia, '(3a)') txtput1, ': ', trim(dll_function_settle(l))
          if (dll_usrfil_settle(l) /= ' ') then
             txtput1 = '  Input for Settle function'
             write (lundia, '(3a)') txtput1, ': ', trim(dll_usrfil_settle(l))
          endif
       endif
       !
       if (sedpar%ero_intfc) then
          txtput1 = '  tcrint'
          write(lundia,'(A)') 'ERROR: SEDPAR%TCRINT not yet defined'
          !write (lundia, '(a,a,e14.6)') txtput1,':',sedpar%tcrint(l)
       endif
       !
       if (sedpar%flnrd(l) /= ' ') then
          txtput1 = '  1D nodal relations for bed/total load'
          write (lundia, '(3a)') txtput1, ':  ', trim(sedpar%flnrd(l))
       endif
    enddo
    if (sedpar%shearsettling) then
       write (lundia, '(a)') 'SLURRY: Fall velocity of sand is computed according to the shear settling formulation (dependent on mud concentration)'
    endif
    if (sedpar%shearsettling) then
       write (lundia, '(a)') 'SLURRY: Fall velocity of sand is computed according to the shear settling formulation (dependent on mud concentration)'
    endif
    if (sedpar%shearsettling) then
       write (lundia, '(a)') 'SLURRY: Fall velocity of sand is computed according to the shear settling formulation (dependent on mud concentration)'
    endif
    if (sedpar%stressStrainRelation) then
       !
       ! Check on maximum Bingham stress
       !
       taubin = sedpar%bin_cyield*2.0**(1.0/sedpar%power) / 1.e+3
       write (lundia, *)
       write (lundia, '(a,e12.4)') 'Indication of maximum Bingham stress : ', &
       &  taubin
       !
       if ( taubin > 0.2) then 
          write (lundia,'(a)') 'Is the Bingham stress possibly too large?'
          write (lundia,'(a)') 'Maximum value for the Bingham stress is 1 Pascal'
       endif
       write (lundia, *)
    endif
    !
    ! Rheology
    !
    if (sedpar%stressStrainRelation) then
       txtput1 = 'rheology'
       txtput2 = 'not specified'
       if (sedpar%rheologymodel == RHEOLOGY_WINTERWERP_KRANENBURG) then
          txtput2 = 'Winterwerp_Kranenburg'
       elseif (sedpar%rheologymodel == RHEOLOGY_JACOBS_VANKESTEREN) then
          txtput2 = 'Jacobs_vanKesteren'
       elseif (sedpar%rheologymodel == RHEOLOGY_THOMAS) then
          txtput2 = 'Thomas'
       endif
       write (lundia, '(3a)') txtput1, ': ', trim(txtput2)
       if (sedpar%rheologymodel == RHEOLOGY_WINTERWERP_KRANENBURG) then
          txtput1 = '  phisim'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_phisim
          txtput1 = '  ayield'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_ayield
          txtput1 = '  frcdim'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_frcdim
          txtput1 = '  bety'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_bety
          txtput1 = '  avic'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_avic
          txtput1 = '  powa'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_powa
          txtput1 = '  betv'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_betv
          txtput1 = '  shrco'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_shrco
       elseif (sedpar%rheologymodel == RHEOLOGY_JACOBS_VANKESTEREN) then
          txtput1 = '  phisim'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_phisim
          txtput1 = '  ayield'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_ayield
          txtput1 = '  powyie'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_powyie
          txtput1 = '  bety'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_bety
          txtput1 = '  avic'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_avic
          txtput1 = '  powvic'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_powvic
          txtput1 = '  betv'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_betv
          txtput1 = '  shrco'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_shrco
       elseif (sedpar%rheologymodel == RHEOLOGY_THOMAS) then
          txtput1 = '  phisim'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_phisim
          txtput1 = '  ayield'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_ayield
          txtput1 = '  powyie'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_powyie
          txtput1 = '  yieldk'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_yieldk
          txtput1 = '  bvic'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_bvic
          txtput1 = '  visck'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_visck
          txtput1 = '  shrco'
          write (lundia, '(2a,e12.4)') txtput1, ':', sedpar%rheo_shrco
       endif
    endif
    !
    write (lundia, '(a)') '*** End    of sediment input'
    write (lundia, *)
end subroutine echosed


subroutine count_sed(lundia    ,error     ,lsed      ,lsedtot   , &
                   & filsed    ,sedpar    ,sed_ptr   )
!!--description-----------------------------------------------------------------
!
! - Determines number of sediment fractions from sediment input file
!
!!--declarations----------------------------------------------------------------
    use properties
    use string_module
    use message_module
    use sediment_basics_module
    use morphology_data_module
    !
    implicit none
!
! Arguments
!
    integer                         , intent(out) :: lsed    ! Number of suspended sediment fractions
    integer                         , intent(out) :: lsedtot ! Total number of sediment fractions
    integer                         , intent(in)  :: lundia  ! Unit of diagnostic file
    logical                         , intent(out) :: error
    character(len=*)                , intent(in)  :: filsed
    type(sedpar_type)               , pointer     :: sedpar
    type(tree_data)                 , pointer     :: sed_ptr
!
! Local variables
!
    integer                                                :: i
    integer                                                :: ibl
    integer                                                :: istat
    integer                                                :: j
    integer                                                :: lsedbl        ! Number of bedload fractions
    integer                                                :: sedtypnr
    integer         , dimension(:) , allocatable           :: typsedim      ! Type of the sediments
    logical                                                :: found
    character(20)                                          :: versionstring
    character(20)   , dimension(:) , allocatable           :: namsedim      ! Names of the sediments as read from sed-file
    character(20)                                          :: sedtyptmp     ! Sediment type in sed-file
    character(80)                                          :: parname
    character(300)                                         :: message
    type(tree_data)                            , pointer   :: asedblock_ptr
!
!! executable statements -------------------------------------------------------
!
    lsedbl        = 0
    lsed          = 0
    lsedtot       = 0
    istat         = 0
    error         = .false.
    !
    ! Check 'Filsed' record for attribute file containing parameters
    ! for sediment transport computation
    !
    if (filsed == ' ') then
       !
       ! file does not exist
       !
       message = 'No sediment file name specified.'
       call write_error(message,unit=lundia)
       error = .true.
       return
    endif 
    !
    ! Read sed-file into tree data structure
    !
    call prop_file('ini', trim(filsed), sed_ptr, istat)
    if (istat /= 0) then
       select case (istat)
          case(1)
             call write_error(FILE_NOT_FOUND//trim(filsed), unit=lundia)
          case(3)
             call write_error(PREMATURE_EOF//trim(filsed), unit=lundia)
          case default
             call write_error(FILE_READ_ERROR//trim(filsed), unit=lundia)
       endselect
       error = .true.
       return
    endif
    !
    if ( .not.associated(sed_ptr%child_nodes) ) then
       call write_error(FILE_READ_ERROR//trim(filsed), unit=lundia)
       error = .true.
       return
    endif
    !
    ! Check version number of sed input file
    !
    versionstring = ' '
    call prop_get_string(sed_ptr, 'SedimentFileInformation', 'FileVersion', versionstring)
    if (trim(versionstring) == '02.00') then
       !
       ! allocate temporary arrays with length equal to the number of data blocks in the file
       ! the sediment blocks will be a subset of that
       !
       allocate(namsedim(size(sed_ptr%child_nodes)))
       allocate(typsedim(size(sed_ptr%child_nodes)))
       namsedim = ' '
       typsedim = -999
       !
       do j = 1, size(sed_ptr%child_nodes)
          !
          ! Does sed_ptr contain a child with name 'Sediment' (converted to lower case)?
          !
          asedblock_ptr => sed_ptr%child_nodes(j)%node_ptr
          parname = tree_get_name( asedblock_ptr )
          if (parname == 'sediment') then
             parname = ' '
             call prop_get_string(asedblock_ptr, '*', 'Name', parname)
             !
             ! Check if the same sediment name was used before
             !
             found = .false.
             do i = 1, j-1
                if (strcmpi(parname,namsedim(i))) then
                   found = .true.
                   exit
                endif
             enddo
             if (found) then
                message = 'Sediment #'//trim(parname)//'# is specified more than once in sediment file '//trim(filsed)
                call write_error(message, unit=lundia)
                error = .true.
                return
             endif
             !
             ! Determine sediment type
             !
             sedtyptmp = ' '
             call prop_get_string(asedblock_ptr, '*', 'SedTyp', sedtyptmp)
             call small(sedtyptmp, 999)
             !
             if (index(sedtyptmp, 'mud') == 1) then
                sedtypnr = SEDTYP_COHESIVE
                lsed = lsed+1
             elseif (index(sedtyptmp, 'sand') == 1) then
                sedtypnr = SEDTYP_NONCOHESIVE_SUSPENDED
                lsed = lsed+1
             elseif (index(sedtyptmp,'bedload') == 1) then
                sedtypnr = SEDTYP_NONCOHESIVE_TOTALLOAD
                lsedbl = lsedbl+1
             else
                message = 'Sediment type of '//trim(parname)//' invalid: '//trim(sedtyptmp)
                call write_error(message,unit=lundia)
                error = .true.
                return
             endif
             !
             namsedim(j) = parname
             typsedim(j) = sedtypnr
          endif
       enddo
    else
       !
       ! sediment input version is 0 or 1:
       ! No keywords and no bedload option
       !
       lsedbl = 0
       call count_sed01(lundia    ,error     ,lsed      ,filsed    , &
                      & namsedim  ,typsedim  )
    endif
    lsedtot = lsed + lsedbl
    !
    ! rhosol, namsed and sedtyp must always be allocated
    !
    if (lsedtot == 0) then
       message = 'No sediment fractions defined in ' // trim(filsed)
       call write_error(message, unit=lundia)
       error = .true.
       return
    endif
    
                    allocate (sedpar%rhosol(lsedtot), stat = istat)
    if (istat == 0) allocate (sedpar%namsed(lsedtot), stat = istat)
    if (istat == 0) allocate (sedpar%sedtyp(lsedtot), stat = istat)
    if (istat == 0) allocate (sedpar%flnrd(0:lsedtot), stat = istat)
    if (istat /= 0) then
       call write_error('Memory allocation error in COUNT_SED', unit=lundia)
    endif
    !
    i   = 1
    ibl = lsed+1
    do j = 1,size(typsedim)
       select case (typsedim(j)) 
          case (SEDTYP_NONCOHESIVE_TOTALLOAD)
             sedpar%namsed(ibl) = namsedim(j)
             sedpar%sedtyp(ibl) = typsedim(j)
             ibl = ibl+1
          case (SEDTYP_NONCOHESIVE_SUSPENDED, SEDTYP_COHESIVE)
             sedpar%namsed(i) = namsedim(j)
             sedpar%sedtyp(i) = typsedim(j)
             i = i+1
          case default
             ! not a sediment block, so continue
       end select
    enddo
    !
    deallocate(namsedim)
    deallocate(typsedim)
end subroutine count_sed


subroutine count_sed01(lundia    ,error     ,lsed      ,filsed    , &
                     & namsedim  ,typsedim  )
!!--description-----------------------------------------------------------------
!
! - Determines number of sediment fractions from sediment input file
!   supports both VERSION 0 and VERSION 1
!
!!--declarations----------------------------------------------------------------
    use string_module
    use message_module
    use sediment_basics_module
    !
    implicit none
!
! Arguments
!
    integer                         , intent(out) :: lsed     ! Number of suspended sediment fractions
    integer                         , intent(in)  :: lundia   ! Unit of diagnostic file
    character(*)                                  :: filsed
    logical                         , intent(out) :: error    ! error flag
    integer         , dimension(:) , allocatable  :: typsedim ! Type of the sediments
    character(20)   , dimension(:) , allocatable  :: namsedim ! Names of the sediments as read from sed-file
!
! Local variables
!
    integer                                                :: i
    integer                                                :: iocond
    integer                                                :: l
    integer                                                :: luninp   ! Unit of .sed file
    integer                                                :: npar
    integer                                                :: version  ! version number (0 or 1)
    character(20)                                          :: sedtyptmp     ! Sediment type in sed-file
    character(300)                                         :: line
    character(300)                                         :: errmsg
!
!! executable statements -------------------------------------------------------
!
    error  = .false.
    call opensedfil(lundia    ,error     ,filsed    ,luninp    ,version  )
    if (error) return
    !
    read (luninp, *, iostat = iocond) lsed ! skip other parameters on line
    !
    allocate(namsedim(lsed))
    allocate(typsedim(lsed))
    !
    do l = 1, lsed
       write(namsedim(l),'(a,1i0)') 'sed.fraction ',l
       read (luninp, '(a)', iostat = iocond) sedtyptmp
       if (iocond == 0) then
          call str_lower(sedtyptmp)
          if (index(sedtyptmp, 'sand') == 1) then
              typsedim(l) = SEDTYP_NONCOHESIVE_SUSPENDED
          elseif (index(sedtyptmp, 'mud') == 1) then
              typsedim(l) = SEDTYP_COHESIVE
          else
             errmsg = 'Invalid suspended sediment type (must start with sand or mud)'
             call write_error(errmsg, unit=lundia)
             error = .true.
             goto 999
          endif
       endif
       !
       ! number of parameter values is 10 for version 0 and 11 for version 1 files
       !
       npar = 10+version
       do i = 1,npar
          read (luninp, *, iostat = iocond) line
          if (iocond /= 0) then
             write(errmsg,'(A,I0)') 'Error reading parameters for sediment fraction',l
             call write_error(errmsg, unit=lundia)
             error = .true.
             goto 999
          endif
       enddo
    enddo
    !
999 continue
    close (luninp)
end subroutine count_sed01


end module m_rdsed
