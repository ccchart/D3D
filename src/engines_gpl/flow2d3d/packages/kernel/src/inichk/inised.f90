subroutine inised(lundia    ,error     ,nmax      ,mmax      ,nmaxus    , &
                & nmmax     ,lsed      ,lsedtot   , &
                & facdss    ,dss       ,kcs       ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2013.                                
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
!    Function: - Initialisation total sediment at bed in each
!                horizontal point
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use mathconsts
    !
    use globaldata
    use bedcomposition_module
    use morphology_data_module, only: allocsedtra
    use sediment_basics_module
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer       , dimension(:)         , pointer :: iform
    real(fp)      , dimension(:)         , pointer :: dm
    real(fp)      , dimension(:)         , pointer :: dg
    real(fp)      , dimension(:)         , pointer :: dgsd
    real(fp)      , dimension(:,:)       , pointer :: dxx
    real(fp)      , dimension(:,:)       , pointer :: frac
    real(fp)      , dimension(:)         , pointer :: mudfrac
    real(fp)      , dimension(:)         , pointer :: sandfrac
    real(fp)      , dimension(:,:)       , pointer :: hidexp
    real(fp)      , dimension(:,:)       , pointer :: sbuuc
    real(fp)      , dimension(:,:)       , pointer :: sbvvc
    real(fp)      , dimension(:,:)       , pointer :: ssuuc
    real(fp)      , dimension(:,:)       , pointer :: ssvvc
    real(fp)      , dimension(:,:)       , pointer :: sucor
    real(fp)      , dimension(:,:)       , pointer :: svcor
    integer                              , pointer :: nxx
    real(fp)      , dimension(:)         , pointer :: xx
    integer                              , pointer :: ihidexp
    real(fp)                             , pointer :: asklhe
    real(fp)                             , pointer :: mwwjhe
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
    real(fp)      , dimension(:)         , pointer :: sedtrcfac
    real(fp)      , dimension(:)         , pointer :: mudcnt
    real(fp)      , dimension(:)         , pointer :: pmcrit
    real(fp)      , dimension(:)         , pointer :: tcguni
    integer       , dimension(:)         , pointer :: nseddia
    integer       , dimension(:)         , pointer :: sedtyp
    character(10) , dimension(:)         , pointer :: inisedunit
    character(256), dimension(:)         , pointer :: flsdbd
    character(256), dimension(:)         , pointer :: flstcg
    logical                              , pointer :: anymud
    character(256)                       , pointer :: flsdia
    character(256)                       , pointer :: flsmdc
    character(256)                       , pointer :: flspmc
    real(fp)                             , pointer :: factcr
!
! Global variables
!
    integer                                            , intent(in)  :: lsed    !  Description and declaration in esm_alloc_int.f90
    integer                                            , intent(in)  :: lsedtot !  Description and declaration in esm_alloc_int.f90
    integer                                            , intent(in)  :: lundia  !  Description and declaration in inout.igs
    integer                                            , intent(in)  :: mmax    !  Description and declaration in esm_alloc_int.f90
    integer                                            , intent(in)  :: nmax    !  Description and declaration in esm_alloc_int.f90
    integer                                            , intent(in)  :: nmaxus  !  Description and declaration in esm_alloc_int.f90
    integer                                            , intent(in)  :: nmmax   !  Description and declaration in esm_alloc_int.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)       , intent(in)  :: kcs     !  Description and declaration in esm_alloc_int.f90
    logical                                                          :: error   !  Flag=TRUE if an error is encountered
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, lsed) , intent(out) :: dss     !  Description and declaration in esm_alloc_real.f90
    real(fp)  , dimension(lsed)                        , intent(in)  :: facdss  !  Description and declaration in esm_alloc_real.f90
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
    iform               => gdp%gdtrapar%iform
    nxx                 => gdp%gdmorpar%nxx
    xx                  => gdp%gdmorpar%xx
    ihidexp             => gdp%gdmorpar%ihidexp
    asklhe              => gdp%gdmorpar%asklhe
    mwwjhe              => gdp%gdmorpar%mwwjhe
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
    sedtrcfac           => gdp%gdsedpar%sedtrcfac
    tcguni              => gdp%gdsedpar%tcguni
    mudcnt              => gdp%gdsedpar%mudcnt
    pmcrit              => gdp%gdsedpar%pmcrit
    nseddia             => gdp%gdsedpar%nseddia
    sedtyp              => gdp%gdsedpar%sedtyp
    inisedunit          => gdp%gdsedpar%inisedunit
    flsdbd              => gdp%gdsedpar%flsdbd
    flstcg              => gdp%gdsedpar%flstcg
    anymud              => gdp%gdsedpar%anymud
    flsdia              => gdp%gdsedpar%flsdia
    flsmdc              => gdp%gdsedpar%flsmdc
    flspmc              => gdp%gdsedpar%flspmc
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
    call allocsedtra(gdp%gderosed, gdp%d%kmax, lsed, lsedtot, &
                   & gdp%d%nmlb, gdp%d%nmub, gdp%d%nmlb, gdp%d%nmub, nxx)
    !
    dm                  => gdp%gderosed%dm
    dg                  => gdp%gderosed%dg
    dgsd                => gdp%gderosed%dgsd
    dxx                 => gdp%gderosed%dxx
    frac                => gdp%gderosed%frac
    mudfrac             => gdp%gderosed%mudfrac
    sandfrac            => gdp%gderosed%sandfrac
    hidexp              => gdp%gderosed%hidexp
    sbuuc               => gdp%gderosed%e_sbnc
    sbvvc               => gdp%gderosed%e_sbtc
    ssuuc               => gdp%gderosed%e_ssnc
    ssvvc               => gdp%gderosed%e_sstc
    sucor               => gdp%gderosed%e_scrn
    svcor               => gdp%gderosed%e_scrt
    !
    ! Start filling array SEDD50FLD
    !
    if (lsedtot==1 .and. flsdia/=' ') then
       !
       !  Space varying data has been specified
       !  Use routine that also read the depth file to read the data
       !
       allocate (gdp%gdsedpar%sedd50fld(gdp%d%nmlb:gdp%d%nmub), stat = istat)
       if (istat /= 0) then
          call prterr(lundia, 'U021', 'Inised: memory alloc error')
          call d3stop(1, gdp)
       endif
       sedd50fld           => gdp%gdsedpar%sedd50fld
       !
       call depfil(lundia    ,error     ,flsdia    ,fmttmp    , &
                 & sedd50fld ,1         ,1         ,gdp%griddim)
       if (error) goto 9999
       !
       call mirror_bnd(icx       ,icy       ,nmmax     , &
                     & kcs       ,sedd50fld ,nmlb      ,nmub      )
       !
    endif
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
       call depfil(lundia    ,error     ,flsmdc    ,fmttmp    , &
                 & mudcnt    ,1         ,1         ,gdp%griddim)
       if (error) goto 9999
    endif
    do nm = 1, nmmax
       mudcnt(nm) = max(0.0_fp, min(mudcnt(nm), 1.0_fp))
    enddo
    !
    ! Read array PMCRIT
    !
    if (flspmc /= ' ') then
       !
       ! Space varying data has been specified
       ! Use routine that also read the depth file to read the data
       !
       call depfil(lundia    ,error     ,flspmc    ,fmttmp    , &
                 & pmcrit    ,1         ,1         ,gdp%griddim)
       if (error) goto 9999
    endif
    do nm = 1, nmmax
       pmcrit(nm) = min(pmcrit(nm), 1.0_fp)
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
          if (sedtyp(ll) /= SEDTYP_COHESIVE) then
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
    call getfrac(gdp%gdmorlyr, frac     ,anymud    ,mudcnt    , &
               & mudfrac     ,gdp%d%nmlb,gdp%d%nmub)
    !
    ! Calculate arithmetic mean sediment diameter Dm
    ! Calculate geometric mean sediment diameter Dg
    ! Calculate percentiles Dxx
    !
    call compdiam(frac      ,sedd50    ,sedd50    ,sedtyp    ,lsedtot   , &
                & logsedsig ,nseddia   ,logseddia ,nmmax     ,nmlb      , &
                & nmub      ,xx        ,nxx       ,sedd50fld ,dm        , &
                & dg        ,dxx       ,dgsd      )
    !
    ! Determine hiding & exposure factors
    !
    if (lsedtot > 1) then
       call comphidexp(frac     ,dm        ,nmmax     ,lsedtot   , &
                     & sedd50   ,hidexp    ,ihidexp   ,asklhe    , &
                     & mwwjhe   ,nmlb      ,nmub      )
    else
       hidexp = 1.0
    endif
    !
    ! Initialise settling velocity for sand
    !
    do ll = 1, lsed
       if (sedtyp(ll) == SEDTYP_NONCOHESIVE_SUSPENDED) then
          s = rhosol(ll)/rhow
          !
          if (sedd50(ll) < 1.5_fp*dsand) then
             ws0(ll) = (s - 1.0_fp)*ag*sedd50(ll)**2/(18.0_fp*vicmol)
          elseif (sedd50(ll) < 0.5_fp*dgravel) then
             ws0(ll) = 10.0_fp*vicmol/sedd50(ll)                      &
                & *(sqrt(1.0_fp + (s - 1.0_fp)*ag*sedd50(ll)**3          &
                &                      /(100.0_fp*vicmol**2)) - 1.0_fp)
          else
             ws0(ll) = 1.1_fp*sqrt((s - 1.0_fp)*ag*sedd50(ll))
          endif
       endif
    enddo
 9999 continue
end subroutine inised
