subroutine inised(lundia    ,error     ,nmax      ,mmax      ,nmaxus    , &
                & nmmax     ,lsed      ,lsedtot   , &
                & facdss    ,dss       ,kcs       ,gdp       )
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
!    Function: - Initialisation total sediment at bed in each
!                horizontal point
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    use globaldata
    use bedcomposition_module
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    real(fp)                             , pointer :: eps
    real(fp)                             , pointer :: eps_plt
    integer,        dimension(:)         , pointer :: iform
    real(fp),       dimension(:,:)       , pointer :: par
    real(fp), dimension(:)               , pointer :: dm
    real(fp), dimension(:)               , pointer :: dg
    real(fp), dimension(:,:)             , pointer :: dxx
    real(fp), dimension(:,:)             , pointer :: frac
    real(fp), dimension(:)               , pointer :: mudfrac
    real(fp), dimension(:,:)             , pointer :: hidexp
    real(fp), dimension(:,:)             , pointer :: sbuuc
    real(fp), dimension(:,:)             , pointer :: sbvvc
    real(fp), dimension(:,:)             , pointer :: ssuuc
    real(fp), dimension(:,:)             , pointer :: ssvvc
    real(fp), dimension(:,:)             , pointer :: sucor
    real(fp), dimension(:,:)             , pointer :: svcor
    type (sv_erosed)                     , pointer :: sverosed
    real(fp)                             , pointer :: dsand
    real(fp)                             , pointer :: dgravel
    real(fp)                             , pointer :: sus
    real(fp)                             , pointer :: bed
    integer                              , pointer :: nxx
    real(fp)              , dimension(:) , pointer :: xx
    real(fp)                             , pointer :: rhow
    real(fp)                             , pointer :: ag
    real(fp)                             , pointer :: vicmol
    real(fp)                             , pointer :: mdcuni
    real(fp)      , dimension(:)         , pointer :: rhosol
    real(fp)      , dimension(:,:,:)     , pointer :: logseddia
    real(fp)      , dimension(:)         , pointer :: logsedsig
    real(fp)      , dimension(:)         , pointer :: sedd50
    real(fp)      , dimension(:)         , pointer :: sedd50fld
    real(fp)      , dimension(:)         , pointer :: cdryb
    real(fp)      , dimension(:)         , pointer :: dstar
    real(fp)      , dimension(:)         , pointer :: taucr
    real(fp)      , dimension(:)         , pointer :: tetacr
    real(fp)      , dimension(:)         , pointer :: ws0
    real(fp)      , dimension(:)         , pointer :: sdbuni
    real(fp)      , dimension(:,:)       , pointer :: tcrdep
    real(fp)      , dimension(:)         , pointer :: tcduni
    real(fp)      , dimension(:,:)       , pointer :: tcrero
    real(fp)      , dimension(:)         , pointer :: tceuni
    real(fp)      , dimension(:,:)       , pointer :: eropar
    real(fp)      , dimension(:)         , pointer :: erouni
    real(fp)      , dimension(:)         , pointer :: mudcnt
    integer       , dimension(:)         , pointer :: nseddia
    character(10) , dimension(:)         , pointer :: inisedunit
    character(4)  , dimension(:)         , pointer :: sedtyp
    character(256), dimension(:)         , pointer :: flsdbd
    character(256), dimension(:)         , pointer :: flstcd
    character(256), dimension(:)         , pointer :: flstce
    character(256), dimension(:)         , pointer :: flsero
    logical                              , pointer :: anymud
    character(256)                       , pointer :: flsdia
    character(256)                       , pointer :: flsmdc
    type (gd_sedpar)                     , pointer :: gdsedpar
    real(fp)                             , pointer :: factcr
!
! Global variables
!
    integer                                                         , intent(in)  :: lsed    !  Description and declaration in iidim.f90
    integer                                                         , intent(in)  :: lsedtot !  Description and declaration in iidim.f90
    integer                                                         , intent(in)  :: lundia  !  Description and declaration in inout.igs
    integer                                                         , intent(in)  :: mmax    !  Description and declaration in iidim.f90
    integer                                                         , intent(in)  :: nmax    !  Description and declaration in iidim.f90
    integer                                                         , intent(in)  :: nmaxus  !  Description and declaration in iidim.f90
    integer                                                         , intent(in)  :: nmmax   !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                    , intent(in)  :: kcs     !  Description and declaration in iidim.f90
    logical                                                                       :: error   !!  Flag=TRUE if an error is encountered
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, lsed)              , intent(out) :: dss     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(lsed)                                     , intent(in)  :: facdss  !  Description and declaration in rjdim.f90
!
! Local variables
!
    integer           :: icx
    integer           :: icy
    integer           :: istat
    integer           :: ll
    integer           :: nm
    integer           :: nmaxddb
    integer           :: nmlb
    integer           :: nmub
    real(fp)          :: drho
    real(fp)          :: s
    character(11)     :: fmttmp ! Format file ('formatted  ') 
!
!! executable statements -------------------------------------------------------
!
    eps                 => gdp%gdconst%eps
    eps_plt             => gdp%gdconst%eps_plt
    iform               => gdp%gdeqtran%iform
    par                 => gdp%gdeqtran%par
    dm                  => gdp%gderosed%dm
    dg                  => gdp%gderosed%dg
    dxx                 => gdp%gderosed%dxx
    frac                => gdp%gderosed%frac
    mudfrac             => gdp%gderosed%mudfrac
    hidexp              => gdp%gderosed%hidexp
    sbuuc               => gdp%gderosed%sbuuc
    sbvvc               => gdp%gderosed%sbvvc
    ssuuc               => gdp%gderosed%ssuuc
    ssvvc               => gdp%gderosed%ssvvc
    sucor               => gdp%gderosed%sucor
    svcor               => gdp%gderosed%svcor
    sverosed            => gdp%gderosed
    dsand               => gdp%gdmorpar%dsand
    dgravel             => gdp%gdmorpar%dgravel
    sus                 => gdp%gdmorpar%sus
    bed                 => gdp%gdmorpar%bed
    nxx                 => gdp%gdmorpar%nxx
    xx                  => gdp%gdmorpar%xx
    rhow                => gdp%gdphysco%rhow
    ag                  => gdp%gdphysco%ag
    vicmol              => gdp%gdphysco%vicmol
    mdcuni              => gdp%gdsedpar%mdcuni
    rhosol              => gdp%gdsedpar%rhosol
    logseddia           => gdp%gdsedpar%logseddia
    logsedsig           => gdp%gdsedpar%logsedsig
    sedd50              => gdp%gdsedpar%sedd50
    sedd50fld           => gdp%gdsedpar%sedd50fld
    cdryb               => gdp%gdsedpar%cdryb
    dstar               => gdp%gdsedpar%dstar
    taucr               => gdp%gdsedpar%taucr
    tetacr              => gdp%gdsedpar%tetacr
    ws0                 => gdp%gdsedpar%ws0
    sdbuni              => gdp%gdsedpar%sdbuni
    tcrdep              => gdp%gdsedpar%tcrdep
    tcduni              => gdp%gdsedpar%tcduni
    tcrero              => gdp%gdsedpar%tcrero
    tceuni              => gdp%gdsedpar%tceuni
    eropar              => gdp%gdsedpar%eropar
    erouni              => gdp%gdsedpar%erouni
    mudcnt              => gdp%gdsedpar%mudcnt
    nseddia             => gdp%gdsedpar%nseddia
    inisedunit          => gdp%gdsedpar%inisedunit
    sedtyp              => gdp%gdsedpar%sedtyp
    flsdbd              => gdp%gdsedpar%flsdbd
    flstcd              => gdp%gdsedpar%flstcd
    flstce              => gdp%gdsedpar%flstce
    flsero              => gdp%gdsedpar%flsero
    anymud              => gdp%gdsedpar%anymud
    flsdia              => gdp%gdsedpar%flsdia
    flsmdc              => gdp%gdsedpar%flsmdc
    gdsedpar            => gdp%gdsedpar
    factcr              => gdp%gdmorpar%factcr
    !
    nmlb    = gdp%d%nmlb
    nmub    = gdp%d%nmub
    nmaxddb = gdp%d%nub - gdp%d%nlb + 1
    !
    icx = 1
    icy = nmaxddb
    !
    fmttmp = 'formatted'
    istat  = 0
    if (.not. associated(sverosed%dm)) then
                     allocate (sverosed%dm     (gdp%d%nmlb:gdp%d%nmub)        , stat = istat)
       if (istat==0) allocate (sverosed%dg     (gdp%d%nmlb:gdp%d%nmub)        , stat = istat)
       if (istat==0) allocate (sverosed%dxx    (gdp%d%nmlb:gdp%d%nmub,nxx)    , stat = istat)
       if (istat==0) allocate (sverosed%frac   (gdp%d%nmlb:gdp%d%nmub,lsedtot), stat = istat)
       if (istat==0) allocate (sverosed%mudfrac(gdp%d%nmlb:gdp%d%nmub)        , stat = istat)
       if (istat==0) allocate (sverosed%hidexp (gdp%d%nmlb:gdp%d%nmub,lsedtot), stat = istat)
       if (istat==0) allocate (sverosed%sbuuc  (gdp%d%nmlb:gdp%d%nmub,lsedtot), stat = istat)
       if (istat==0) allocate (sverosed%sbvvc  (gdp%d%nmlb:gdp%d%nmub,lsedtot), stat = istat)
       if (istat==0) allocate (sverosed%ssuuc  (gdp%d%nmlb:gdp%d%nmub,lsed)   , stat = istat)
       if (istat==0) allocate (sverosed%ssvvc  (gdp%d%nmlb:gdp%d%nmub,lsed)   , stat = istat)
       if (istat==0) allocate (sverosed%sucor  (gdp%d%nmlb:gdp%d%nmub,lsed)   , stat = istat)
       if (istat==0) allocate (sverosed%svcor  (gdp%d%nmlb:gdp%d%nmub,lsed)   , stat = istat)
       if (istat/=0) then
          call prterr(lundia, 'U021', 'Inised: memory alloc error', gdp)
          call d3stop(1, gdp)
       endif
    dm                  => gdp%gderosed%dm
    dg                  => gdp%gderosed%dg
    dxx                 => gdp%gderosed%dxx
    frac                => gdp%gderosed%frac
    mudfrac             => gdp%gderosed%mudfrac
    hidexp              => gdp%gderosed%hidexp
    sbuuc               => gdp%gderosed%sbuuc
    sbvvc               => gdp%gderosed%sbvvc
    ssuuc               => gdp%gderosed%ssuuc
    ssvvc               => gdp%gderosed%ssvvc
    sucor               => gdp%gderosed%sucor
    svcor               => gdp%gderosed%svcor
    sverosed            => gdp%gderosed
       !
       dm      = 0.0
       dg      = 0.0
       dxx     = 0.0
       mudfrac = 0.0
       frac    = 0.0
       sucor   = 0.0
       svcor   = 0.0
    endif
    !
    ! Initialise cumulative sediment transport arrays
    !
    sbuuc = 0.0
    sbvvc = 0.0
    ssuuc = 0.0
    ssvvc = 0.0
    !
    ! Start filling array SEDD50FLD
    !
    if (lsedtot==1 .and. flsdia/=' ') then
       !
       !  Space varying data has been specified
       !  Use routine that also read the depth file to read the data
       !
       allocate (gdsedpar%sedd50fld(gdp%d%nmlb:gdp%d%nmub), stat = istat)
       if (istat /= 0) then
          call prterr(lundia, 'U021', 'Inised: memory alloc error', gdp)
          call d3stop(1, gdp)
       endif
    mdcuni              => gdp%gdsedpar%mdcuni
    rhosol              => gdp%gdsedpar%rhosol
    logseddia           => gdp%gdsedpar%logseddia
    logsedsig           => gdp%gdsedpar%logsedsig
    sedd50              => gdp%gdsedpar%sedd50
    sedd50fld           => gdp%gdsedpar%sedd50fld
    cdryb               => gdp%gdsedpar%cdryb
    dstar               => gdp%gdsedpar%dstar
    taucr               => gdp%gdsedpar%taucr
    tetacr              => gdp%gdsedpar%tetacr
    ws0                 => gdp%gdsedpar%ws0
    sdbuni              => gdp%gdsedpar%sdbuni
    tcrdep              => gdp%gdsedpar%tcrdep
    tcduni              => gdp%gdsedpar%tcduni
    tcrero              => gdp%gdsedpar%tcrero
    tceuni              => gdp%gdsedpar%tceuni
    eropar              => gdp%gdsedpar%eropar
    erouni              => gdp%gdsedpar%erouni
    mudcnt              => gdp%gdsedpar%mudcnt
    nseddia             => gdp%gdsedpar%nseddia
    inisedunit          => gdp%gdsedpar%inisedunit
    sedtyp              => gdp%gdsedpar%sedtyp
    flsdbd              => gdp%gdsedpar%flsdbd
    flstcd              => gdp%gdsedpar%flstcd
    flstce              => gdp%gdsedpar%flstce
    flsero              => gdp%gdsedpar%flsero
    anymud              => gdp%gdsedpar%anymud
    flsdia              => gdp%gdsedpar%flsdia
    flsmdc              => gdp%gdsedpar%flsmdc
    gdsedpar            => gdp%gdsedpar
       !
       call depfil(lundia    ,error     ,flsdia    ,fmttmp    ,mmax      , &
                 & nmax      ,nmaxus    ,sedd50fld ,gdp       )
       if (error) goto 9999
       !
       call mirror_bnd(icx       ,icy       ,nmmax     , &
                     & kcs       ,sedd50fld ,nmlb      ,nmub      )
       !
    endif
    !
    ! Fill sediment dependent arrays
    !
    do ll = 1, lsed
       !
       ! tcrdep; only for mud
       !
       if (sedtyp(ll) == 'mud') then
          if (flstcd(ll) == ' ') then
             !
             ! Uniform data has been specified
             !
             tcrdep(:, ll) = tcduni(ll)
          else
             !
             ! Space varying data has been specified
             ! Use routine that also reads the depth file to read the data
             !
             call depfil(lundia    ,error     ,flstcd(ll),fmttmp    ,mmax      , &
                       & nmax      ,nmaxus    ,tcrdep(nmlb, ll)     ,gdp       )
             if (error) goto 9999
          endif
          !
          ! Check whether taucr for deposition is zero somewhere
          ! If so, give a warning for the first cell where that occurs
          !
          do nm = 1, nmmax
             if (kcs(nm) == 1) then
                if (comparereal(tcrdep(nm,ll), eps_plt) /= 1) then
                   call prterr(lundia, 'G051', 'Critical shear stress for deposition of mud is 0.0 in at least one cell', gdp)
                   exit
                endif
             endif
          enddo
       endif
       !
       ! tcrero; only for mud
       !
       if (sedtyp(ll) == 'mud') then
          if (flstce(ll) == ' ') then
             !
             ! Uniform data has been specified
             !
             tcrero(:, ll) = tceuni(ll)
          else
             !
             ! Space varying data has been specified
             ! Use routine that also reads the depth file to read the data
             !
             call depfil(lundia    ,error     ,flstce(ll),fmttmp    ,mmax      , &
                       & nmax      ,nmaxus    ,tcrero(nmlb, ll)     ,gdp       )
             if (error) goto 9999
          endif
          do nm = 1, nmmax
             if (kcs(nm) == 1) then
                if (comparereal(tcrero(nm,ll), eps_plt) /= 1) then
                   call prterr(lundia, 'U021', 'Critical shear stress for erosion of mud must be > 0.0', gdp)
                   call d3stop(1, gdp)
                endif
             endif
          enddo
       endif
       !
       ! eropar
       !
       if (flsero(ll) == ' ') then
          !
          ! Uniform data has been specified
          !
          eropar(:, ll) = erouni(ll)
       else
          !
          ! Space varying data has been specified
          ! Use routine that also read the depth file to read the data
          !
          call depfil(lundia    ,error     ,flsero(ll),fmttmp    ,mmax      , &
                    & nmax      ,nmaxus    ,eropar(nmlb, ll)     ,gdp       )
          if (error) goto 9999
       endif
    enddo
    !
    ! Start filling array MUDCNT
    !
    if (flsmdc == ' ') then
       !
       ! Uniform data has been specified
       !
       mudcnt = mdcuni
    else
       !
       ! Space varying data has been specified
       ! Use routine that also read the depth file to read the data
       !
       call depfil(lundia    ,error     ,flsmdc    ,fmttmp    ,mmax      , &
                 & nmax      ,nmaxus    ,mudcnt    ,gdp       )
       if (error) goto 9999
    endif
    do nm = 1, nmmax
       mudcnt(nm) = max(0.0_fp, min(1.0_fp, mudcnt(nm)))
    enddo
    !
    ! Initialise suspended sediment diameter
    !
    if (lsedtot==1 .and. lsed==1 .and. flsdia/=' ') then
       do nm = 1, nmmax
          dss(nm, 1) = sedd50fld(nm)*facdss(1)
       enddo
    else
       do ll = 1, lsed
          dss(:, ll) = sedd50(ll)*facdss(ll)
       enddo
    endif
    !
    ! Calculation of dimensionless grain size and critical shear stress
    ! Only for uniform sedd50
    ! For space varying sedd50:
    ! - this is done every time step, for every nm in erosed
    ! - the do-loop below may not be executed due to the usage of the uninitialised sedd50(ll)
    !
    if (lsedtot/=1 .or. lsed/=1 .or. flsdia==' ') then
       do ll = 1, lsedtot
          if (sedtyp(ll)=='sand' .or. sedtyp(ll)=='bedl') then
              drho      = (rhosol(ll)-rhow) / rhow
              dstar(ll) = sedd50(ll) * (drho*ag/vicmol**2)**0.3333_fp
              if (dstar(ll) < 1.0_fp) then
                 if (iform(ll) == -2) then
                    tetacr(ll) = 0.115_fp / (dstar(ll)**0.5_fp)
                 else
                    tetacr(ll) = 0.24_fp / dstar(ll)
                 endif
              elseif (dstar(ll) <= 4.0_fp) then
                 if (iform(ll) == -2) then
                    tetacr(ll) = 0.115_fp / (dstar(ll)**0.5_fp)
                 else
                    tetacr(ll) = 0.24_fp / dstar(ll)
                 endif
              elseif (dstar(ll)>4.0_fp .and. dstar(ll)<=10.0_fp) then
                 tetacr(ll) = 0.14_fp  / (dstar(ll)**0.64_fp)
              elseif (dstar(ll)>10.0_fp .and. dstar(ll)<=20.0_fp) then
                 tetacr(ll) = 0.04_fp  / (dstar(ll)**0.1_fp)
              elseif (dstar(ll)>20.0_fp .and. dstar(ll)<=150.0_fp) then
                 tetacr(ll) = 0.013_fp * (dstar(ll)**0.29_fp)
              else
                 tetacr(ll) = 0.055_fp
              endif
              taucr(ll) = factcr * (rhosol(ll)-rhow) * ag * sedd50(ll) * tetacr(ll)
           else
              dstar(ll)  = 0.0_fp
              tetacr(ll) = 0.0_fp
              taucr(ll)  = 0.0_fp
           endif
       enddo
    endif
    !
    ! Initialise morphology layers
    !
    call inimorlyr(flsdbd    ,sdbuni    ,inisedunit,cdryb     , &
                 & lsedtot   ,mmax      ,nmax      ,nmaxus    ,nmmax     , &
                 & lundia    ,error     ,kcs       ,gdp       )
    !
    ! Inilialise fractions
    !
    call getfrac(gdp%gdmorlyr, cdryb ,frac     , &
               & sedtyp    ,anymud   ,mudcnt   ,mudfrac   )
    !
    ! Calculate arithmetic mean sediment diameter Dm
    !
    call compdmean(frac      ,sedd50    ,nmmax     ,lsedtot   , &
                 & sedtyp    ,dm        ,sedd50fld ,logsedsig , &
                 & gdp       )
    !
    ! Calculate geometric mean sediment diameter Dg
    !
    call compdgeomean(frac      ,sedd50    ,nmmax     ,lsedtot   , &
                    & sedtyp    ,dg        ,sedd50fld ,gdp       )
    !
    ! Calculate percentiles Dxx
    !
    call compdxx(frac      ,nseddia   ,logseddia ,logsedsig , &
               & nmmax     ,lsedtot   ,sedtyp    ,dxx       , &
               & xx        ,nxx       ,sedd50fld ,gdp       )
    !
    ! Determine hiding & exposure factors
    !
    if (lsedtot > 1) then
       call comphidexp(frac     ,dm        ,nmmax     ,lsedtot   , &
                     & sedd50   ,hidexp    ,gdp       )
    else
       hidexp = 1.0
    endif
    !
    ! Initialise settling velocity for sand
    !
    do ll = 1, lsed
       if (sedtyp(ll) == 'sand') then
          s = rhosol(ll)/rhow
          !
          if (sedd50(ll) < 1.5*dsand) then
             ws0(ll) = (s - 1.0)*ag*sedd50(ll)**2/(18.0*vicmol)
          elseif (sedd50(ll) < 0.5*dgravel) then
             ws0(ll) = 10.0*vicmol/sedd50(ll)                      &
                & *(sqrt(1.0 + (s - 1.0)*ag*sedd50(ll)**3          &
                &                      /(100.0*vicmol**2)) - 1.0)
          else
             ws0(ll) = 1.1*sqrt((s - 1.0)*ag*sedd50(ll))
          endif
       endif
    enddo
 9999 continue
end subroutine inised
