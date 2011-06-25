module mod_trisim
!!--copyright-------------------------------------------------------------------
! Copyright (c) 2007, WL | Delft Hydraulics. All rights reserved.
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
!    Function: Main routine for the 2d / 3d program
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!
contains
!


!
!==============================================================================
integer function trisim_init(numdom, nummap, context_id, fsm_flags, fsm_tracefile, runid_arg, gdp) result(retval)
!
!!--declarations----------------------------------------------------------------
!
    use precision
    use SyncRtcFlow
    use dfparall
    use timers
    use meteo_data
    !
    use m_openda_exchange_items, only : openda_buffer_initialize
    !
    ! global data declaration; compare with include 'globdat.igd'
    !
    use globaldata
    implicit none
    type(globDat),target :: gdp
    !
    integer                       , pointer :: lundia
    integer                       , pointer :: lunprt
    integer                       , pointer :: iphisi
    integer, dimension(:)         , pointer :: ipmap
    logical                       , pointer :: dredge
    logical                       , pointer :: struct
    integer                       , pointer :: numdomains
    integer                       , pointer :: nummappers
    character(6)                  , pointer :: prognm
    integer                       , pointer :: rtcmod
    include 'flow_steps_f.inc'
    include 'fsm.i'
!
! Parameters
!
    integer       , intent(in)  :: context_id
    integer       , intent(in)  :: fsm_flags
    integer       , intent(in)  :: numdom        ! Number of subdomains (0 means single domain)
                                                 ! as detected by hydra
    integer       , intent(in)  :: nummap        ! Number of mappers (one for each DD boundaries connected with this subdomain)
                                                 ! as detected by hydra
    character(*)  , intent(in)  :: fsm_tracefile
    character(256)              :: runid_arg

!
! Local variables
!
    integer        , pointer            :: itinit
    integer        , pointer            :: itstrt
    integer                             :: fsmstatus
    integer                             :: i
    integer                             :: ic           ! Length of character parameter CASE 
    integer                             :: icheck
    integer        , pointer            :: initia
    integer        , pointer            :: initi        ! Control parameter 
                                                        ! = 1 initialization
                                                        ! = 2 initialization and read restart information from the communication file
                                                        ! = 3 no initialization 
    integer        , pointer            :: it01         ! Reference date in yymmdd 
    integer        , pointer            :: it02         ! Reference time in hhmmss 
    integer        , pointer            :: itb          ! Start time of computational interval 
    integer        , pointer            :: ite          ! End time of computational interval 
    integer        , pointer            :: itima        ! Time to start simulation (N * tscale) according to DELFT3D conventions 
    integer        , pointer            :: itlen        ! Lenght of the tide cycle 
    integer                             :: lenid
    integer                             :: lunid
    integer                             :: luntri       ! Unit number for trigger file for TRISIM for running programs simultaniously 
    integer                             :: nhystp
    integer                  , external :: newlun
    integer                  , external :: fsmtrf
    logical        , pointer            :: alone        ! TRUE when flow runs stand-alone, FALSE when flow is part of morsys 
    logical                             :: ex
    logical                             :: init         ! Flag=TRUE when initialisation is re- quired (always the case if FLOW is used stand alone) 
    logical                             :: lexist
    logical        , pointer            :: mainys       ! Logical flag for FLOW is main porgram (TRUE) for writing output 
    logical                             :: opend        ! Help logical var. to determine whether each of the output files was opened 
    real(fp)       , pointer            :: tscale       ! Basic unit time 
    character(12)                       :: filmrs       ! File name for DELFT3D_MOR FLOW input file (MD-flow.xxx) 
    character(12)                       :: filsim       ! Name for trigger file for TRISIM for running programs simultaniously 
    character(256)                      :: case         ! Project identification (a non-blank character string presumed) 
    character(256) , pointer            :: comfil       ! Communication file name 
    character(256) , pointer            :: filmd        ! File name for MD FLOW file 
    character(256) , pointer            :: trifil       ! File name for FLOW NEFIS output files (tri"h/m"-"casl""labl".dat/def) 
    character(4)                        :: subsys       ! Sub-system definition of Delft3D here SUBSYS = 'flow' 
    character(5)                        :: filid
    character(5)   , pointer            :: versio       ! Version nr. of the current package 
    character(256) , pointer            :: runid
!
!! executable statements -------------------------------------------------------
!
    !
    ! Start simulation performance measurement
    !
    call StartComputation
    !
    ! Initialization using a semaphore
    ! Related vseminit is in tricom.f90
    !
    call pseminit
    !
    ! Initializes MPI
    !
    call dfinitmpi
    !
    ! initialize GDP structure
    ! allocate(gdp) MUST have been done in trisim/trisim_dll
    !
    ! openDA buffer
    !
    call openda_buffer_initialize
    !    
    call gdp_alloc(gdp)
    call initsafe(gdp)
    call timers_init(gdp)
    call timer_start(timer_total, gdp)
    call timer_start(timer_init, gdp)
    !
    ! esm/fsm initialization
    !
    fsmstatus = fsmini (context_id, fsm_flags)
    if (len (fsm_tracefile) > 0) then
       fsmstatus = fsmtrf (fsm_tracefile)
    endif
    !
    lundia       => gdp%gdinout%lundia
    lunprt       => gdp%gdinout%lunprt
    iphisi       => gdp%gdinttim%iphisi
    ipmap        => gdp%gdinttim%ipmap
    itinit       => gdp%gdinttim%itinit
    itstrt       => gdp%gdinttim%itstrt    
    dredge       => gdp%gdprocs%dredge
    struct       => gdp%gdprocs%struct
    numdomains   => gdp%gdprognm%numdomains
    nummappers   => gdp%gdprognm%nummappers
    prognm       => gdp%gdprognm%prognm
    rtcmod       => gdp%gdrtc%rtcmod
    initia       => gdp%gdtricom%initia
    initi        => gdp%gdtricom%initi 
    it01         => gdp%gdtricom%it01
    it02         => gdp%gdtricom%it02
    itb          => gdp%gdtricom%itb
    ite          => gdp%gdtricom%ite
    itima        => gdp%gdtricom%itima
    itlen        => gdp%gdtricom%itlen
    alone        => gdp%gdtricom%alone
    mainys       => gdp%gdtricom%mainys
    tscale       => gdp%gdtricom%tscale
    comfil       => gdp%gdtricom%comfil
    runid        => gdp%gdtricom%runid
    trifil       => gdp%gdtricom%trifil
    versio       => gdp%gdtricom%versio
    filmd        => gdp%gdtricom%filmd
    !
    ! Initialize local parameters, including IPHISI and IPMAP(1)
    ! in case program crashes the test below can be performed anyway
    !
    runid    = runid_arg
    
    init     = .true.
    filmrs   = ' '
    subsys   = 'flow'
    alone    = .true.
    !
    iphisi   = 0
    ipmap(1) = -1
    !
    rtcmod   = noRTC
    icheck   = 0
    !
    ! Store numdom (counted by Hydra) in numdomains (in GDP-structure)
    ! For single domain cases, Hydra does not count at all and numdom is zero.
    !
    numdomains = max(1,numdom)
    !
    ! Store nummap (counted by Hydra) in nummapperss (in GDP-structure)
    !
    nummappers = nummap
    if (runid==' ') then
       !
       ! First try to read runid from file called RUNID
       ! This simplifies fluidmud synchronisation
       !
       filid = 'runid'
       inquire (file = filid, exist = ex)
       if (ex) then
          lunid = newlun(gdp)
          open (lunid, file = filid, form = 'formatted', status = 'old')
          read (lunid, '(a)') runid
          close (lunid)
       else
          runid = ' '
       endif
    else
       !
       ! Remove (possible) trailing '\0' from c-code
       !
       do i = 1, len(runid)
          if (ichar(runid(i:i)) == 0 .or. ichar(runid(i:i)) == 10) runid(i:i) = ' '
       enddo
    endif
    runid = adjustl(runid)
    !
    ! Platform dependent initialization
    !
    call pldep
    !
    ! Run TDATOM
    !
    if (.not.parll .or. inode == master) then
       call tdatmain(runid, alone, subsys, filmrs, icheck, gdp) 
    endif
    call dfbroadc(icheck, 1, dfint,gdp)
    call dfsync(gdp)
    if (icheck /= 0 ) then
       write (*, '(a)') 'ABORT: error returned by tdatmain'
       lundia = 0
       call d3stop(1,gdp)
    endif
    !
    ! Set program name (after running tdatom)
    !
    prognm = 'TRISIM'
    !
    ! Determine by trigger-file if RTC is running as well
    !
    luntri = newlun(gdp)
    filsim = 'TMP_SYNC.RUN'
    inquire (file = filsim, exist = lexist)
    if (lexist) then
       open (luntri, file = filsim, form = 'unformatted', status = 'unknown')
       read (luntri) icheck
       close (luntri)
       !
       ! Check 'RUNRTC' by telephone
       !
       if (icheck==786782) then
          rtcmod = dataFromRTCToFLOW
       else
          write (*, '(a)') 'Trigger-file TMP_SYNC.RUN not made by TDATOM'
          call d3stop(1         ,gdp       )
       endif
    endif
    !
    ! Initialize sub-system for Delft3D-FLOW
    !
    call defsub(subsys    ,gdp       )
    !
    ! Start FLOW simulation program
    !
    call noextspaces(runid     ,lenid     )
    if (init) then
       !
       ! Read  dimensions of arrays and declare array pointers
       !
       call tripoi(runid, filmrs, versio, filmd, &
                 & alone, gdp)
       if (gdp%errorcode > 0) then
          if (rtcmod == dataFromRTCToFLOW) then
             call timer_start(timer_wait, gdp)
             call syncflowrtc_quit
             call timer_stop(timer_wait, gdp)
             rtcmod = noRTC
          endif
          !
          ! Error: try to close and return -1
          !
          call d3stop(1, gdp)
          retval = trisim_close(gdp)
          retval = -1
          return
          !
       endif
    endif
    !
    ! Initialize time frame parameters for stand alone program
    !
    initi  = 1
    !
    it01   = 0
    it02   = 0
    !
    itb    = 1
    ite    = -1
    itlen  = 0
    tscale = 1.0
    !
    itima  = 0
    !
    mainys = .true.
    itima  = 0
    !
    ! Initialize communication file name
    ! NOTE: case may never be only blanks
    !
    case = runid
    call noextspaces(case      ,ic        )
    !
    comfil(1:4) = 'com-'
    comfil(5:)  = case(1:ic)
    !
    trifil(1:5) = 'trix-'
    trifil(6:)  = case(1:ic)
    !
    ! Insert node number in file names
    !
    if ( parll ) then
       write(comfil(5+ic:5+ic+4),'(a,i3.3)') '-',inode
    endif
    !
    ! The following is necessary; originally, tricom was called with parameter initi
    ! while its n-th argument was initia. Confusing. (VT)
    !
    initia = initi
    !
    ! Start FLOW simulation
    !
    call tricom_init(gdp)
    !
    ! Error status of tricom_init is returned via gdp%errorcode
    !
    if (gdp%errorcode /= 0) then
        retVal = -1
    else
        retVal = 0
    endif
    !
    ! for OpenDA purposes, itinit and itfinish replace itstrt and itstop
    !
    itinit = itstrt
    !    
end function trisim_init
!
!
!
!----------------------------------------------------------------------
integer function trisim_step(gdp) result(retval)
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    RetVal = 0
    !
    ! the subroutine called 'tricom_verify' form the BAW-version contains
    ! part VII and VIII of the initialisation. That subroutine is not needed here and
    ! part VII and VIII can be found at the end of tricom_init. (VT)
    call tricom_step(gdp)
    if (gdp%errorcode /= 0) then
        retVal = -1
    else
        retVal = 0
    endif
    !
end function trisim_step
!
!
!
!-----------------------------------------------------------------------
integer function trisim_finish(gdp) result(retVal)
    use globaldata
    !    
    implicit none
    !    
    type(globdat),target :: gdp
    integer            , pointer :: lundia
 
    lundia       => gdp%gdinout%lundia 
    retval = 0

    call tricom_finish(gdp)

    write(lundia,*)
    write(lundia,'(a)') '*** Simulation finished *******************************************************'
    write(lundia,*)

    if (gdp%errorcode>0) then
    endif
    !
    ! Write diagnostics and close file
    !
    retval = trisim_close(gdp)
    !
end function trisim_finish
!
!
!
!-----------------------------------------------------------------------
integer function trisim_close(gdp) result (retval)
    use timers
    use meteo_data
    use globaldata
    !    
    implicit none
    !    
    type(globdat),target :: gdp
    !  
    character(256)                , pointer :: filmd        ! File name for MD FLOW file 
    integer                       , pointer :: lundia
    integer                       , pointer :: lunprt
    integer                       , pointer :: iphisi    
    integer                       , pointer :: numdomains
    integer                       , pointer :: nummappers
    integer        , dimension(:) , pointer :: ipmap    
    character(256)                , pointer :: runid
    character(5)                  , pointer :: versio
    !
    include 'flow_steps_f.inc'
    include 'fsm.i'
    !
    integer                     :: nhystp
    logical                     :: opend        ! Help logical var. to determine whether each of the output files was opened 
    logical                     :: init         ! Flag=TRUE when initialisation is re- quired (always the case if FLOW is used stand alone) 
    !
    ! body
    !
    runid        => gdp%gdtricom%runid
    versio       => gdp%gdtricom%versio
    filmd        => gdp%gdtricom%filmd
    lundia       => gdp%gdinout%lundia
    lunprt       => gdp%gdinout%lunprt
    iphisi       => gdp%gdinttim%iphisi    
    ipmap        => gdp%gdinttim%ipmap
    nummappers   => gdp%gdprognm%nummappers
    !
    ! Write diagnostics and close file
    !
    init = .true.
    !
    call timer_stop(timer_total, gdp)
    call timers_finish(gdp)
    if (init) then
       call triend(runid, gdp)
       if (gdp%errorcode > 0) then
          !
          ! User interaction is removed
          !
          !
          ! To avoid statemachine problems:
          !
          nhystp = nxtstp(d3dflow_init, gdp)
       endif
    endif
    !
    ! Close intermediate files TMP//runid//.*
    !
    call delfil(runid     ,filmd     ,gdp       )
    !
    ! Close diagnostic file, print file and dredge file
    ! No prints requested (LUNPRT still unit 8 see SYSINI) delete file
    !
    inquire (lundia, opened = opend)
    if (opend) close (lundia)
    if (iphisi>0 .or. ipmap(1)>=0) then
       !
       ! Only then the tri-prt file can actually be present
       !
       inquire (lunprt, opened = opend)
       if (opend) then
          close (lunprt)
       endif
    endif
    !
    call deallocmeteo(gdp%runid)
    !
    ! Finish using a semaphore
    ! Related psemnefis is in tricom.f90
    ! This used to be a vsemfinish
    !
    call vsemfinish
    !
    ! Tell gaws (Global ADI Wang Solver) and mapper we're done.
    !
    if (nummappers >= 1) then
        call gwsslv(-1)
        nhystp = nxtstp(d3dflow_finish, gdp)
    endif
    !
    ! Stop simulation performance measurement
    !
    call EndComputation
    !
    retval = 0  
    !
end function trisim_close
!
!
!
!==============================================================================
integer function trisim_initialise_single_step(gdp) result(retval)
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    integer        , pointer :: itinit
    integer        , pointer :: itstrt
    integer        , pointer :: itstop
    integer        , pointer :: itfinish
    !
    ! body
    !
    itinit       => gdp%gdinttim%itinit
    itstrt       => gdp%gdinttim%itstrt
    itstop       => gdp%gdinttim%itstop
    itfinish     => gdp%gdinttim%itfinish
    retval = 0
    itinit = itstrt
    itstop = itstrt + 1
    if (itstop > itfinish) then
      retVal = -1
    endif
end function trisim_initialise_single_step
!
!
!
!==============================================================================
integer function trisim_prepare_next_step(gdp) result(retval)
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    integer           , pointer :: itinit
    integer           , pointer :: itstrt
    integer           , pointer :: itstop
    integer           , pointer :: itfinish
    !
    ! body
    !
    itinit       => gdp%gdinttim%itinit
    itstrt       => gdp%gdinttim%itstrt
    itstop       => gdp%gdinttim%itstop
    itfinish     => gdp%gdinttim%itfinish
    retval = 0
    itstrt = itstop
    itstop = itstrt + 1
end function trisim_prepare_next_step
!
!
!
!==============================================================================
integer function trisim_check_step(gdp) result(retval)
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp

    integer                             , pointer :: itstop
    integer                             , pointer :: itfinish
    !
    ! body
    !
    itstop       => gdp%gdinttim%itstop
    itfinish     => gdp%gdinttim%itfinish
    retval = 0
    if (itstop > itfinish) then
      retVal = -1
    endif
end function trisim_check_step
!
!
!
!---------------------------------------------------------------
end module mod_trisim
