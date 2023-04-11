!----- AGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2017-2023.
!
!  This file is part of Delft3D (D-Flow Flexible Mesh component).
!
!  Delft3D is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Affero General Public License as
!  published by the Free Software Foundation version 3.
!
!  Delft3D  is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Affero General Public License for more details.
!
!  You should have received a copy of the GNU Affero General Public License
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.
!
!  contact: delft3d.support@deltares.nl
!  Stichting Deltares
!  P.O. Box 177
!  2600 MH Delft, The Netherlands
!
!  All indications and logos of, and references to, "Delft3D",
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!
!-------------------------------------------------------------------------------

! 
! 

 !> Initialise flow model time dependent parameters
 !! @return Integer error status (0) if succesful.
 integer function flow_flowinit() result(iresult)
 use timers
 use m_netw
 use m_flowgeom
 use m_flow
 use m_flowtimes
 use m_sferic
 use unstruc_model
 use unstruc_files
 use m_reduce, only : nodtot, lintot
 use m_samples
 use m_missing
 use m_fixedweirs
 use m_partitioninfo
 use m_sediment
 use m_transport
 use dfm_error
 use m_sobekdfm
 use m_monitoring_crosssections, only : crs, ReallocCrosssectionSums
 use string_module, only: str_lower
 use m_ship
 use geometry_module ! , only: dbdistance, normalout, half
 use m_physcoef, only: backgroundwatertemperature
 use m_alloc
 use unstruc_channel_flow, only: network
 use m_1d_structures, only: initialize_structures_actual_params, t_structure
 use m_oned_functions, only: updateFreeboard, set_max_volume_for_1d_nodes, updateDepthOnGround, updateVolOnGround, updateTotalInflow1d2d, updateTotalInflowLat
 use m_waves
 use m_structures
 use m_longculverts
 use unstruc_channel_flow, only: useVolumeTables
 use m_VolumeTables
 use timers
 use m_setucxcuy_leastsquare, only: reconst2nd
 use mathconsts, only: sqrt2_hp

 implicit none

 ! locals
 integer :: k, L, k1, k2, n, jw, msam
 integer :: kb, kt, ki, LL

 double precision :: r,eer,r0,dep,Rossby,amp,csth, sqghi, ss, snth, x, y, uth
 double precision :: xx, yy, zz, ux, uy, pin, xli, atet, slope, cz, z00, cs
 double precision :: xx1, yy1, xx2, yy2, ux1, uy1, ux2, uy2, csl, snl
 double precision :: xzmin, xzmax, yzmin, yzmax, bobmin
 double precision :: fout, foutk, aa, blm, dis, dmu, var, rho1, hunsat, fac, zi, zido, ziup, saldo, salup
 double precision, external :: rho_Eckart
 double precision, allocatable:: u1_tmp(:)

 integer :: itest = 1, ispecials, ierr, mrst, kk, La, j, nq, ierror, N1, N2, Lb, Lt, nat, ntmp
 integer :: ihandle
 double precision :: psi, samp, ct, st, omeg, t, rr, rmx, x0, y0, dxx, dyy, blmn, blmx, dbl, bot, rms, ucmk, phi, dphi
 character(len=255) :: rstfile
 character(len=Idlen) :: fileName
 character(len=4)   :: EXT
 logical :: jawel, jawelrestart
 integer :: nstrucsg, L0, istru
 type(t_structure), pointer :: pstru
 integer, external :: flow_initexternalforcings
 double precision, external :: setrho

 double precision  :: xm, ym
 double precision  :: trshcorioi
 double precision  :: Ds
 double precision  :: hw,tw,csw,snw, uorbi,rkw,ustt,hh, upot, ukin, ueaa

 iresult = DFM_GENERICERROR

 if (ndx == 0) then
    iresult = DFM_MODELNOTINITIALIZED
    goto 888
 end if

 Lnmax = 0 ; Lnmin = 0
 ndmax = 0 ; ndmin = 0

 call inisferic()                                    ! also set coriolis :<
 if (icorio > 0) then
    call inifcori()
 endif
 if (Corioadamsbashfordfac > 0) then
    if (allocated(fvcoro) ) deallocate(fvcoro)
    allocate ( fvcoro(lnkx), stat = ierr ) ; fvcoro = 0d0
    call aerr('fvcoro(lnkx)', ierr, lnkx  )
 endif
 if (jsferic == 0) then
    if (jatidep > 0) then
       ! call qnerror('Tide generating potential only supported for sferical models ' ,' ',' ')
       jatidep = 0
    endif
    if (jaselfal > 0) then
       ! call qnerror('Self attraction and loading only supported for sferical models ' ,' ',' ')
       jaselfal = 0
    end if
 endif
 call inidensconstants()                             ! Some density parameters

if (ti_waq > 0d0 .and. max(limtypmom, limtypsa, limtypTM) <= 0) then
    call qnerror('DELWAQ requires at least one limiter (Numerical Parameters). DELWAQ output disabled for now.', ' ', ' ')
    ti_waq = 0d0
end if

 s1  = sini                                  ! initial values
 if (waterdepthini1D > 0) then
    do k = ndx2D+1, ndxi
       s1(k) = bl(k) + waterdepthini1D
    enddo
 endif

 u1  = uini
 if (jasal > 0) then
    sa1 = salini
 endif
 if (jatem > 0) then
    tem1 = temini
 endif

 ! spiral flow
 if (jasecflow > 0 ) then
    spirint = spirini
 endif

 if (jased > 0 .and. jased < 4) then
    do k = 1,ndkx
       do j = 1,mxgr
          sed(j,k) = sedini(j)
       enddo
    enddo
 endif

 if (jasal == 0 .and. jatem == 0 .and. jased == 0) then
    idensform = 0
 endif
 volerror = 0d0 ; squ = 0 ; sqi = 0

 call statisticsini()

 call setkbotktop(1)                                 ! prior to correctblforzlayerpoints, setting kbot

 call mess(LEVEL_INFO, 'Start initializing external forcings...')
 call timstrt('Initialise external forcings', handle_iniext)
 iresult = flow_initexternalforcings()               ! this is the general hook-up to wind and boundary conditions
 call timstop(handle_iniext)

 ! from hereon, the processes are in sync
 if (jampi == 1) then
    ! globally reduce the error
    call reduce_error(iresult)
 end if
 if (iresult /= DFM_NOERR) then
    call qnerror('Error occurred while initializing external forcings, please inspect the preceding lines in the diagnostic output for more details.',' ', ' ')
    if (jampi == 1) then
        call qnerror('Error occurs on one or more processes when initializing external forcings.',' ', ' ')
    end if
   goto 888
 end if
 call mess(LEVEL_INFO, 'Done initializing external forcings.')

 ! Set ihorvic related to horizontal viscosity
 ihorvic = 0
 if (vicouv > 0 .or. javiusp == 1 .or. Smagorinsky > 0 .or. Elder > 0) then
    ihorvic = 1
 endif

 ! If constituents have been added at this point, the sum-arrays in crs require redimensioning
 if (allocated(crs)) then
    call ReallocCrossSectionSums(crs)
 endif

 ! call setbobsongullies()

 if (isimplefixedweirs == 0) then
    call setbobs_fixedweirs()
 else
    call setfixedweirs()
 endif

 call delpol()

 if (Slopedrop2D > 0) then !todo, uitsluitende test maken
    do L  = lnx1D+1,lnxi
       k1 = ln(1,L) ; k2 = ln(2,L)
       if (iadv(L) /=0 .and. .not. (iadv(L) >= 21 .and. iadv(L) <=25) .and. dxi(L)*abs(bl(k1) - bl(k2)) > Slopedrop2D) then ! Not for fixed weirs itself.
           iadv(L) = 8
       endif
    enddo
 endif

 do L  = 1,lnxi
    if (iadv(L) /= 0) then
       ! only when advection calculation is on, change the advection type
       if (kcu(L) == 3) then         ! lateral overflow
          if (iadveccorr1D2D == 2) then ! Switch off advection on 1D2D lateral/embedded links
             iadv(L) = 0
          else
             iadv(L) = 8                ! Original lateral overflow
          endif
       else if (kcu(L) == 5 .or. kcu(L) == 7) then   ! pipe connection
          iadv(L) = 8
       endif
    endif
 enddo

 if (japure1D > 0) then
    call setiadvpure1D()
 endif

! check if at most one structure claims a flowlink
 call check_structures_and_fixed_weirs()

 ! First call to setexternalforcingsonboundaries, here only for the structure timeseries (prior to adjust_bobs_for_dams_and_structs())
 call setzminmax()                                 ! our side of preparation for 3D ec module
 call setsigmabnds()
 call flow_setexternalforcingsonboundaries(tstart_user, iresult)  ! set structure (and bnd) external forcings. Error handling later in 2nd call for bnds.
 call initialize_structures_actual_params(network%sts)            ! After structure time series, and prior to adjust_bobs, to use proper crest levels.


 call adjust_bobs_for_dams_and_structs()
 call setup_structures_and_weirs_list()

 ! Floodfill water levels based on sample file.
 if (len_trim(md_s1inifile) > 0) then
    call savesam()
    NS = 0
    call oldfil(msam, md_s1inifile)
    if (msam /= 0) then
        call reasam(msam, 0)
        if (jampi > 0) then
            call mess(LEVEL_WARN, 'Filling water level using [geometry] WaterLevIniFile in .mdu does not exchange information between partitions')
        endif
        call flow_initfloodfill()
    end if
    call restoresam()
 end if

 if (allocated(ibot)) then
    deallocate(ibot)                                    ! after meteoiniti of ibedlevtype
 end if

 if (ndx > 800000) then
    if (allocated (zk) ) then
       ! deallocate(zk)
    endif
 endif

 call setFrictionForLongculverts()

 do L = 1,lnx
    if (frcu(L) == dmiss) then
       if (L <= lnx1D) then
          if (kcu(L) == 3) then
             frcu(L)  = frcuni1d2d
          else if (kcu(L) == 5) then       ! 1D2D internal link pipe streetinlet
             ! Because frcunistreetinlet is not available in the mdu file, the friction type is
             ! always manning.
             frcu(L)  = frcunistreetinlet
             ifrcutp(L) = 1
          else if (kcu(L) == 7) then       ! 1D2D internal link pipe roofgutterpipe
             ! Because frcuniroofgutterpipe is not available in the mdu file, the friction type is
             ! always manning
             frcu(L)  = frcuniroofgutterpipe
             ifrcutp(L) = 1
          else
             frcu(L)  = frcuni1d
          endif
       else
          frcu(L) =  frcuni
       endif
    endif
    if (ifrcutp(L) == -999) then
        ifrcutp(L) = ifrctypuni
    endif
    if (frcu(L) > frcmax) then
        frcmax = frcu(L)
    endif
 enddo

 if (jafrculin == 1) then                           ! plus uniform value on not found
    do L = 1,lnx
       if (frculin(L) == dmiss) then
           frculin(L) =  frcunilin
       endif
    enddo
 endif

 if ( jaFrcInternalTides2D.eq.1 ) then
    do k=1,Ndx
       if ( frcInternalTides2D(k).eq.DMISS ) then
          frcInternalTides2D(k) = 0d0
       end if
    end do
 end if

 call setupwslopes()                                   ! set upwind slope pointers and weightfactors

 if (iuvfield > 0) call setvelocityfield()           ! only when testing

 ! remember initial water levels at the water level boundaries
 ! so that reading rst file won't influence it. This is used for restart a model with Riemann boundary conditions.
 do n=1,nbndz
    k2 = kbndz(2,n)
    zbndz0(n) = max(bl(k2), s1(k2)) ! NOTE: the s1=max(bl, s1) step can only be done later, so do it here as well.
 end do

 ! Load restart file (*_map.nc) assigned in the *.mdu file OR read a *.rst file
 jawel = .false.
 if (len_trim(md_restartfile) > 0 ) then
!    Find file extension based on first full stop symbol '.' at the back of the string.
     N1  = INDEX (md_restartfile,'.', .true.)
     N2  = len_trim(md_restartfile)
     EXT = ' '
     EXT = md_restartfile(N1:N2)
     ! Restart from *.rst:
     if ( index(md_restartfile, '.rst') > 0 .or. index(md_restartfile, '.RST') > 0) then
         INQUIRE(FILE = rstFILE, EXIST=JAWEL)
         IF (JAWEL) THEN
            call oldfil(mrst, rstfile)
            call rearst(mrst, jw)
            JAWEL = (jw == 1)
         endif
     else ! Restart from *_yyyymmdd_hhmmss_rst.nc or from *_map.nc
       call read_restart_from_map(md_restartfile, iresult)
       if (jased > 0 .and. stm_included) then
          call setbobs()
       end if

       if (jampi == 1) then
          ! globally reduce the error
          call reduce_error(iresult)
       end if
       if (iresult /= DFM_NOERR) then
          if (jampi == 1) then
              call qnerror('Error occurs on one or more processes when reading the restart file.',' ', ' ')
          else
              call qnerror('Error occurs when reading the restart file.',' ', ' ')
          end if

          goto 888
       else
          JAWEL = .true.
       end if

       ntmp = size(u1)
       allocate(u1_tmp(ntmp))
       u1_tmp = u1
       u1     = u0
       hs     = s0 - bl
      if (iperot == -1) then
         call reconst2nd ()
      endif
       call fill_onlyWetLinks()
       call setucxucyucxuucyunew() !reconstruct cell-center velocities
       u1     = u1_tmp
       deallocate(u1_tmp)
     end if
 end if

 jawelrestart = jawel
 if (jawel) jarestart = 1                                       ! in the module

 call flow_setstarttime()                                       ! the flow time0 and time1 are managed by flow
                                                                ! this is the only function that a user can use to influence the flow times
                                                                ! TSTART MAY BE OVERWRITTEN IN REARST

 if (jased > 0 .and. stm_included) then
    if (stmpar%morpar%morft < eps10) then
        !
        ! if the morphological start time is not set to some positive value due
        ! to restart from mapfile, then make sure that the morphological start
        ! time corresponds to the hydrodynamic start time. This includes TStart!
        !
        stmpar%morpar%morft  = tstart_user/86400d0
        stmpar%morpar%morft0 = stmpar%morpar%morft
    endif
 endif

 call setkbotktop(1)                                            ! set sigmabnds for ec

 if ( janudge.eq.1 ) call setzcs()
 call flow_setexternalforcings(tstart_user, .true., iresult)             ! set field oriented external forcings, flag that the call is from the initialization phase

 if (iresult /= DFM_NOERR) then
    goto 888
 end if

 ! Actual boundary forcing (now that initial water levels, etc. are also known):
 call flow_setexternalforcingsonboundaries(tstart_user, iresult)         ! set bnd   oriented external forcings
 if (jampi == 1) then
    ! globally reduce the error
    call reduce_error(iresult)
 end if
 if (iresult /= DFM_NOERR) then
    if (jampi == 1) then
        call qnerror('Error occurs on one or more processes when setting external forcings on boundaries.',' ', ' ')
    end if
    goto 888
 end if
 tim1bnd = tstart_user
 tim1fld = tstart_user

 if (jaoldrstfile==1) then ! If the restart file is of old version (which does not have waterlevel etc info on boundaries), then need to set.
    call sets01zbnd(0, 0)
 endif
 call sets01zbnd(1, 1)

 do n  = 1, nbndn                                  ! for normal velocity boundaries, also initialise velocity on link
    LL = kbndn(3,n)
    call getLbotLtop(LL,Lb,Lt)
    do L = Lb, Lt
       u1(L) = zbndn(n)
    enddo
 enddo

 do nq = 1,nqbnd                                    ! discharge boundaries
    nat = 0 ; bobmin = huge(1d0)
    do n   = L1qbnd(nq), L2qbnd(nq)                 ! initially only
       kb  = kbndu(1,n)
       k2  = kbndu(2,n)
       L   = kbndu(3,n)
       blm = min( bob(1,L), bob(2,L) )
       if (s1(k2) - blm > epshu) then
          nat = 1
       endif
       bobmin = min(bobmin,blm)
    enddo

!   boundary is dry: add 1 cm of water above lowest bed level
    if (nat == 0) then
       do n   = L1qbnd(nq), L2qbnd(nq)             ! initially only
          kb  = kbndu(1,n)
          k2  = kbndu(2,n)
          s1(k2) = max(s1(k2), bobmin + 0.01d0)
          s1(kb) = s1(k2)
       enddo
    else
       do n   = L1qbnd(nq), L2qbnd(nq)             ! initially only
          kb  = kbndu(1,n)
          k2  = kbndu(2,n)
          s1(kb) = s1(k2)
       enddo
    endif
 enddo


 Do L = Lnxi+1, Lnx   ! copy 1D bnd arrays to that of internal attached link
    if (kcu(L) == -1) then
        LL         = Lbnd1D(L)
        frcu(L)    = frcu(LL)
        ifrcutp(L) = ifrcutp(LL)
        IF (jaconveyance2D > 0) aifu(L)    = aifu(LL)
    endif
 enddo

 if (lnx > lnxi) then                               ! boundaries always implicit
    teta(lnxi+1:lnx) = 1d0
 endif

 if (nonlin1d == 2 .or. nonlin1d == 3 .or.nonlin2D == 2) then
    if (allocated(s1mini) ) deallocate(s1mini)
    allocate  ( s1mini(ndx) , stat= ierr)
    call aerr ('s1mini(ndx)', ierr, ndx ) ; s1mini = bl
    do L = 1,lnx1D
       k1 = ln(1,L) ; k2 = ln(2,L)
       if ( abs(PROF1D(3,L)) == 1 ) then
          teta(L)    = 1d0                          ! closed pipes always implicit
          s1mini(k1) = max( s1mini(k1), bl(k1) + prof1D(1,L) )
          s1mini(k2) = max( s1mini(k2), bl(k2) + prof1D(1,L) )
       else if ( abs(PROF1D(3,L)) == 2 .or. abs(PROF1D(3,L)) == 3 ) then
          if ( prof1D(2,L) .ne. dmiss ) then
              teta(L) = 1d0
          endif
       endif
    enddo
    do LL = lnxi+1, lnx
       if (kcu(LL) == -1 ) then
          k1 = ln(1,LL) ; k2 = ln(2,LL)
          L  = LBND1D(LL)
          if ( abs(PROF1D(3,L)) == 1 ) then
             teta(LL)   = 1d0                          ! closed pipes always implicit
             s1mini(k1) = max( s1mini(k1), bl(k1) + prof1D(1,L) )
             s1mini(k2) = max( s1mini(k2), bl(k2) + prof1D(1,L) )
          else if ( abs(PROF1D(3,L)) == 2 .or. abs(PROF1D(3,L)) == 3 ) then
             if ( prof1D(2,L) .ne. dmiss ) then
                 teta(L) = 1d0
             endif
          endif
       endif
    enddo

 endif

! Set teta for all structure links to 1.0 (implicit)
 nstrucsg = network%sts%count
 do istru = 1, nstrucsg
    pstru => network%sts%struct(istru)
    do L0 = 1, pstru%numlinks
       L = iabs(pstru%linknumbers(L0))
       teta(L) = 1d0
    enddo
 enddo

 if (.not. jawelrestart) then
    if (japatm > 0 .and. PavIni > 0) then
       if ( trshcorio.gt.0d0 ) then
          trshcorioi = 1d0/trshcorio
       end if
       do k = 1,ndxi
!          ss    = min( 1d0,  max(s1(k),bl(k)) - bl(k) )        ! reduced correction values at low depths
          Ds = - ( patm(k) - PavIni ) / (ag*rhomean)
          s1(k) = s1(k) + Ds
       end do
    endif

    if  ( jaselfal.gt.0 .and. jaSELFALcorrectWLwithIni.eq.1 ) then
       do k=1,Ndxi
           s1init(k) = max(s1(k), bl(k))
       end do
    end if

    u0 = u1
 endif

 s1  = max(bl, s1)
 s00 = s1

 nonlin = max(nonlin1D, nonlin2D)
 
 if (nonlin >= 2) then
    if (allocated(s1m) ) deallocate (s1m, a1m)
    allocate ( s1m(ndx), a1m(ndx) , STAT=ierr) ; s1m = s1
    call aerr('s1m(ndx), a1m(ndx)', ierr, ndx)
 endif

 if (nshiptxy > 0) then
    call setship()                               ! in flowinit
    if (kmx > 0 .and. jasal > 0) then
       inquire(file = 'verticalsalinityprofile.pli', exist=success )
       call setkbotktop(1)
       if (success) then
          call setinitialverticalprofile(sa1 , ndkx, 'verticalsalinityprofile.pli') ; success = .true.
       endif
    endif
 endif

 call setkbotktop(1)

 if (.not. jawelrestart) then
    s0 = s1
    hs = s1 - bl
 else ! If one restarts a simulation, then use s0 to compute hs
    s0 = max(s0,bl)
    hs = s0 - bl
 endif

 if ( jaselfal.gt.0 ) then
!  with hs available: recompute SAL potential
   call flow_settidepotential(tstart_user/60d0)
 end if


 if (jagrw > 0) then
    do k = 1,ndx

       if (hs(k) > epshs) then
           sgrw1(k) = bl(k)
       else
           if (allocated(h_unsat) ) then
               sgrw1(k) = bl(k) - h_unsat(k)
           else if (infiltrationmodel == 1) then
               sgrw1(k) = bl(k) - Hinterceptionlayer
           else if (h_unsatini > 0d0) then
               sgrw1(k) = bl(k) - h_unsatini
           else
               sgrw1(k) = sgrwini
           endif
       endif
       sgrw1(k) = min( bl(k), sgrw1(k) )
       bgrw(k)  = min( bl(k), bgrw (k) )
       hunsat   = bl(k) - sgrw1(k)
       fac      = min ( 1d0, max(0d0,  hunsat / h_transfer  )  )      ! 0 at bed, 1 at sgrw
       pgrw(k)  = sgrw1(k)*fac + s1(k)*(1d0-fac)

    enddo

    if (allocated(h_unsat) ) deallocate(h_unsat)
    sgrw0 = sgrw1

 endif

 if (infiltrationmodel ==  DFM_HYD_INFILT_CONST .or. infiltrationmodel == DFM_HYD_INFILT_DARCY) then  ! set infiltcap=0 for closed links only
    call realloc(kcsini, ndx, keepExisting=.false., fill = 0)
    do L = 1,lnx  ! only one connected open profile will open surface runoff
       n1 = ln(1,L) ; n2 = ln(2,L)
       if (L <= lnx1D) then
          if (prof1D(3,L) < 0) then ! closed profile
          else
              kcsini(n1) = 1 ; kcsini(n2) = 1
          endif
       else
          kcsini(n1) = 1 ; kcsini(n2) = 1
       endif
    enddo
    do n = 1,ndx
       infiltcap(n) = infiltcap(n)*kcsini(n)  ! 0 for all links closed
    enddo
    if (allocated(kcsini)) deallocate(kcsini)
 endif

 call sethu(1)
 if (kmx > 0) then ! temporary fix for sepr 3D
    do L = 1,lnx
       if (abs(kcu(L)) == 1) then
          call addlink1D(L,1)
          if (hu(L) > 0d0) then
             wu(L) = au(L) / hu(L)
          endif
       endif
    enddo
 endif

  ihandle = 0
 call timstrt('makeVolumeTables', ihandle)
 if (useVolumeTables) then
     filename = defaultFilename('volumeTables')
     call makeVolumeTables(filename)
 endif
 call timstop(ihandle)

 call volsur()          ! flowinit
 call a1vol1tot()
 vol0tot = vol1tot
 a0tot   = a1tot
 vol0    = vol1

 ! initial velocity in 3D (needs Lbot, Ltop)
 if ( inivel.eq.1 .and. kmx.gt.0 ) then
    do L=1,lnx
       call getLbotLtop(L,Lb,Lt)
       u1(Lb:Lt) = u1(L)
    end do
 end if

  if ((jawave==3 .or. jawave==6) .and. .not. flowWithoutWaves) then
    ! Normal situation: use wave info in FLOW
    hs = max(hs, 0d0)
    if (jawave == 6) then
      ! HSIG is read from SWAN NetCDF file. Convert to HRMS
      hwav = hwavcom / sqrt2_hp
    else
      hwav = hwavcom
    endif
    hwav = min(hwav, gammax*hs)
    !
    call wave_uorbrlabda()                      
    if( kmx == 0 ) then
       call wave_comp_stokes_velocities()
       call tauwave()
    end if
    call setwavfu()
    call setwavmubnd()
 end if

 if ((jawave==3 .or. jawave==6) .and. flowWithoutWaves) then
    ! Exceptional situation: use wave info not in FLOW, only in WAQ
    ! Only compute uorb
    ! Works both for 2D and 3D
    if (jawave == 6) then
      ! HSIG is read from SWAN NetCDF file. Convert to HRMS
      hwav = hwavcom / sqrt2_hp
    else
      hwav = hwavcom
    endif
    hwav = min(hwav, gammax*hs)
    call wave_uorbrlabda()                       ! hwav gets depth-limited here
 end if

 if (jawave==5 .and. .not. flowWithoutWaves) then
    hs = max(hs, 0d0)
    hwav = min(hwavcom, gammax*hs)
    call wave_uorbrlabda()
    if( kmx == 0 ) then
       do L=1,lnx
          k1=ln(1,L); k2=ln(2,L)
          hh = hu(L); hw=0.5d0*(hwav(k1)+hwav(k2));tw=.5d0*(twav(k1)+twav(k2))
          csw = 0.5*(cos(phiwav(k1)*dg2rd)+cos(phiwav(k2)*dg2rd))
          snw = 0.5*(sin(phiwav(k1)*dg2rd)+sin(phiwav(k2)*dg2rd))
          call tauwavehk(hw, tw, hh, uorbi, rkw, ustt)
          ustokes(L) = ustt*(csu(L)*csw + snu(L)*snw)
          vstokes(L) = ustt*(-snu(L)*csw + csu(L)*snw)
       enddo
       call tauwave()
    endif
 endif

 if (jasal > 0) then
    if (kmx > 0 .and. inisal2D >= 1 .and. jarestart.eq.0 ) then
    do kk = 1,ndx
       call getkbotktop(kk,kb,kt)
       if (inisal2D == 2) then
          do k = kb, kt
             if (kt == kb) then
                rr  = 1d0
             else
                rr  = dble(k-kb)/dble(kt-kb)
             endif
             sa1(k) = (1d0 - rr)*sa1(kb) + rr*satop(kk)
          enddo
       else if (inisal2D == 3) then          ! uniform below is specified
          do k = kb, kt
             zz = 0.5d0*( zws(k) + zws(k-1) )
             if (zz < uniformsalinitybelowz .and. sabot(kk) .ne. dmiss) then
                sa1(k) = sabot(kk)
             else
                sa1(k) = sa1(kk)
             endif
          enddo
       endif
       do k = kt+1, kb+kmxn(kk)-1
          sa1(k) = sa1(max(kt,kb))
       enddo
    enddo

    if ( allocated(satop) ) then
       deallocate (satop)
    endif
    if ( allocated(sabot) ) then
       deallocate (sabot)
    endif
 endif

    do k = 1,ndkx
       sa1(k) = max( 0d0,  sa1(k) )
    enddo

    if (Sal0abovezlev .ne. dmiss) then
    do kk = 1,ndx
       call getkbotktop(kk,kb,kt)
       do k = kb, kt
             if (zws(k) > Sal0abovezlev) then
                 sa1(k) = 0d0
             endif
       enddo
    enddo
 endif

    salmax = maxval(sa1)

 endif

 if (kmx > 0 .and. inised2D > 0 ) then
    do kk = 1,ndx
       if (sedh(kk) .ne. dmiss) then
          call getkbotktop(kk,kb,kt)
          do k = kb, kt
              sed(1:mxgr,k) = sedh(kk)
          enddo
       endif
    enddo
    deallocate(sedh)
    inised2D = 0
 endif

  ! When restart, initialize salinity, temperature, sed on waterlevel boundaries
 ! NOTE: keep this identical to how it's done at the end of transport()
 ! hk: and, make sure this is done prior to fill constituents
 if (jarestart > 0) then
    do LL = lnxi + 1, lnx                           ! copy on outflow
        call getLbotLtop(LL,Lb,Lt)
        if (Lt < Lb) then
            cycle
        endif
        do L = Lb, Lt
           if (q1(L) <= 0d0) then
              kb = ln(1,L) ; ki = ln(2,L)
              if (jasal > 0) then
                  sa1(kb)  = sa1(ki)
              endif
              if (jatem > 0) then
                  tem1(kb)  = tem1(ki)
              endif
              if (jased > 0) then
                 do j = 1,mxgr
                    sed(j,kb) = sed(j,ki)
                 enddo
              endif
           endif
        enddo
    enddo
    call restore_au_q1_3D_for_1st_history_record()
 endif

 if ( janudge.eq.1 ) then  ! and here last actions on sal/tem nudging, before we set rho
    call set_nudgerate()
    if ( jainiwithnudge > 0 ) then
       call set_saltem_nudge()
       if (jainiwithnudge == 2) then
           janudge = 0
           deallocate (nudge_tem, nudge_sal, nudge_rate , nudge_time)
       endif
    end if
 end if

 if (jasal > 0) then ! used to be in fill_constituents each step
         do k = 1,ndkx
       constituents(isalt,k) = max( 0d0,  sa1(k) )
    enddo
 endif

 if (itemp > 0) then ! used to be in fill_constituents each step
    do k = 1,ndkx
            constituents(itemp,k) = tem1(k)
         enddo
     endif

 if (jainirho == 1) then
    do kk = 1,ndx  ! initialise rho's
       call setrhokk(kk)
       call getkbotktop(kk,kb,kt)
       do k = kt+1 , kb + kmxn(kk) - 1
          if (stm_included) rhowat(k) = rhowat(kt)   ! UNST-5170
       enddo
    enddo
 endif

 if (allocated (rho0))  then
    rho0 = rho
 endif

 if (jaFlowNetChanged == 1 .or. nodtot /= ndx .or. lintot /= lnx) then
       call reducept(Ndx,Ndxi,Lnx)                              ! also alloc arrays for reduce
       if (icgsolver == 10) then
          call alloc_jacobi(ndx,lnx)
       endif
 end if


! BEGIN DEBUG
! if ( jampi.eq.1 ) then
!    call update_ghosts(ITYPE_U,Lnx,teta,ierr)
! endif
! END DEBUG

 if (kmx < 2) then                    ! in 2D, use 1
    if ( ja_timestep_auto.ne.-123 ) then
       ja_timestep_auto = 1
    else
       ja_timestep_auto = 0
    end if
 endif

 if ( jaimplicit.eq.1 ) then
    call inisolver_advec(ierror)
 end if

 call ini_filter(jafilter, filterorder, jacheckmonitor, ierror)

 if (jabarrieradvection == 3) then
    call setstruclink()
 endif

 if (japillar > 0) then
    call setpillars()
 endif

 ! for 1D only
 if (network%loaded .and. ndxi-ndx2d > 0) then
    if (jamapVolOnGround > 0) then
       call set_max_volume_for_1d_nodes() ! set maximal volume, it will be used to update the volume on ground level for the output
    end if
 end if

 call upotukinueaa(upot, ukin, ueaa)

 iresult = DFM_NOERR
 return

888 continue  ! Some error occurred, prevent further flow
 ndx = 0

 if (nonlin>=2) then
    s1m = bl
 endif

end function flow_flowinit

    
!> restore au and q1 for 3D case for the first write into a history file    
subroutine restore_au_q1_3D_for_1st_history_record()
use m_flow,                 only : q1, LBot, kmx, kmxL
use m_flowexternalforcings, only : fusav, rusav, ausav
use m_flowgeom,             only : lnx

implicit none

integer                       :: i_q1_v, i_q1_0
double precision, allocatable :: fu_temp(:,:), ru_temp(:,:), au_temp(:,:)

if ( kmx > 0 ) then
   fu_temp = fusav
   ru_temp = rusav
   au_temp = ausav
   call furusobekstructures() ! to have correct au values but it provides incorrect q1 values for structures
   fusav = fu_temp
   rusav = ru_temp
   ausav = au_temp
!  restore correct discharge values
   do i_q1_0 = 1, lnx
      q1(i_q1_0) = 0d0 
      do i_q1_v = Lbot(i_q1_0), Lbot(i_q1_0) - 1 + kmxL(i_q1_0)
         q1(i_q1_0) = q1(i_q1_0) + q1(i_q1_v)       ! depth integrated result
      end do
   end do
end if
    
end subroutine restore_au_q1_3D_for_1st_history_record