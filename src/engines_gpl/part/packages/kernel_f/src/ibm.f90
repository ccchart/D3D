module ibm_mod
!
!  data definition module(s)
!
use precision_part               ! single/double precision
use timers
!
!  module procedure(s)
!
!use stop_exit_mod           ! explicit interface
use larvm2_mod              ! explicit interface
use tidal_state_mod         ! explicit interface
!
implicit none

contains
      subroutine ibm    ( lunrep   , itime    , idelt    , nmax     , mmax     ,    &
                          layt     , nosegl   , nolay    , mnmaxk   , lgrid    ,    &
                          lgrid2   , lgrid3   , nopart   , npwndw   , nosubs   ,    &
                          npart    , mpart    , kpart    , xpart    , ypart    ,    &
                          zpart    , wpart    , iptime   , wsettl   , locdep   ,    &
                          nocons   , const    , conc     , xa       , ya       ,    &
                          angle    , vol1     , vol2     , volume   ,flow     , salin1   ,    &
                          temper1  , v_swim   , d_swim   , itstrtp  , vel1     ,    &
                          vel2     )

      ! function  : calculates individual based model specifics

      ! arguments :

      integer(ip), intent(in)    :: lunrep              ! report file
      integer(ip), intent(in)    :: itime               ! time in seconds
      integer(ip), intent(in)    :: idelt               ! time step size in seconds
      integer(ip), intent(in)    :: nmax                ! first grid dimension
      integer(ip), intent(in)    :: mmax                ! second grid dimension
      integer(ip), intent(in)    :: layt                ! number of layers of hydr. database
      integer(ip), intent(in)    :: nosegl              ! number segments per layer
      integer(ip), intent(in)    :: nolay               ! number of layers in calculation
      integer(ip), intent(in)    :: mnmaxk              ! total number of active grid cells
      real   (sp), intent(in   ) :: volume( * )         ! grid cell volume (full grid)
      integer(ip), pointer       :: lgrid ( : , : )     ! grid with active grid numbers, negatives for open boundaries
      integer(ip), pointer       :: lgrid2( : , : )     ! total grid
      integer(ip), pointer       :: lgrid3( : , : )     ! original grid (conc array)
      integer(ip), intent(in)    :: nopart              ! total number of particles
      integer(ip), intent(in)    :: npwndw              ! first active particle
      integer(ip), intent(in)    :: nosubs              ! number of substances per particle
      integer(ip), pointer       :: npart ( : )         ! first  grid index of the particles
      integer(ip), pointer       :: mpart ( : )         ! second grid index of the particles
      integer(ip), pointer       :: kpart ( : )         ! third grid index of the particles
      real   (sp), pointer       :: xpart ( : )         ! x-value (0.0-1.0) first  direction within grid cell
      real   (sp), pointer       :: ypart ( : )         ! y-value (0.0-1.0) second direction within grid cell
      real   (sp), pointer       :: zpart ( : )         ! z-value (0.0-1.0) third  direction within grid cell
      real   (sp), pointer       :: wpart ( : , :)      ! weight factors of the subs per particle
      integer(ip), pointer       :: iptime( : )         ! particle age in seconds
      real   (sp), pointer       :: wsettl( : )         ! settling per particle
      real   (sp), pointer       :: locdep( : , : )     ! depth per layer
      integer(ip), intent(in)    :: nocons              ! number of constants
      real   (sp), pointer       :: const ( : )         ! user-defined constants
      real   (sp), pointer       :: conc  ( : , : )     ! concentration array in transport grid
      real   (sp), pointer       :: xa    ( : )         ! x-coordiante in real world
      real   (sp), pointer       :: ya    ( : )         ! y-coordinate in real world
      real   (sp), pointer       :: angle ( : )         ! angle with horizontal
      real   (sp), pointer       :: vol1  ( : )         ! volume begin hydr step
      real   (sp), pointer       :: vol2  ( : )         ! volume end hydr step
      real   (sp), pointer       :: flow  ( : )         ! all flows
      real   (sp), pointer       :: salin1( : )         ! salinity segment numbering
      real   (sp), pointer       :: temper1( : )        ! temperature segment numbering
      real   (sp), pointer       :: v_swim( : )         ! horizontal swimming velocity m/s
      real   (sp), pointer       :: d_swim( : )         ! horizontal swimming direction (degree)
      integer(ip), intent(in)    :: itstrtp             ! start time
      real   (sp), pointer       :: vel1  ( : )         ! velocity begin hydr step
      real   (sp), pointer       :: vel2  ( : )         ! velocity end hydr step

      ! from input (const)

      integer, save              :: nstage              ! number of stages
      integer, save              :: iday = 0            ! for csv output
      integer, save              :: ncum = 0            ! for csv output
      character*256              :: filcsv
      integer, save              :: luncsv = 631
      real(sp), allocatable,save :: astage(:)           ! a coefficient in stage development (-)
      real(sp), allocatable,save :: bstage(:)           ! b coefficient in stage development (-)
      real(sp), allocatable,save :: length(:)           ! Length to grow within the stage (mm)
      integer , allocatable,save :: btype(:)            ! behaviour type
      real(sp), allocatable,save :: mort1(:)            ! base mortality begin stage
      real(sp), allocatable,save :: mort2(:)            ! base mortality end stage
      real(sp), allocatable,save :: tcmort(:)           ! temperature coefficient mortality
      real(sp), allocatable,save :: buoy1(:)            ! buoyancy begin stage
      real(sp), allocatable,save :: buoy2(:)            ! buoyancy end stage
      real(sp), allocatable,save :: vzact1(:)           ! active vertical velocity begin stage
      real(sp), allocatable,save :: vzact2(:)           ! active vertical velocity end stage
      real(sp), allocatable,save :: ztop1(:)            ! depth top layer begin stage
      real(sp), allocatable,save :: ztop2(:)            ! depth top layer end stage
      real(sp), allocatable,save :: zbot1(:)            ! depth bot layer begin stage
      real(sp), allocatable,save :: zbot2(:)            ! depth bot layer end stage
      real(sp), allocatable,save :: phase(:)            ! phase in diurnal behaviour
      real(sp), allocatable,save :: hbtype(:)           ! horizontal behaviour type
      real(sp), allocatable,save :: ampli1(:)            ! amplitude of signal
      real(sp), allocatable,save :: ampli2(:)            ! amplitude of signal
      real(sp), allocatable,save :: vswim1(:)           ! swim velocity at begin of stage
      real(sp), allocatable,save :: vswim2(:)           ! swim velocity at end of stage
      integer, save              :: istage_nursery      ! stage from where larvae can settle in nursery area
      integer, save              :: layer_release       ! layer where eggs are released
      integer, save              :: it_start_m2         ! start time m2 output
      integer, save              :: it_stop_m2          ! stop time m2 output
      integer, save              :: idt_m2              ! timestep m2 output
      integer, save              :: idt_csv             ! timestep csv output, 0 = not, times as m2
      integer, save              :: iniday              ! release of initial condition, not used here
      integer, save              :: tide_opt            ! option in tidal state determination

      ! local :

      real(sp), allocatable,save :: phase_diurn(:)      ! phase in diurnal behaviour
      real   (sp)                :: ampli               ! amplitude in diurnal behaviour
      integer(ip)                :: ipart               ! particle index
      integer(ip)                :: ilay                ! layer index (for m2 output)
      real   (sp)                :: stage               ! stage development
      real   (sp)                :: fstage              ! fraction of current stage
      real   (sp)                :: stage_dur           ! length development (herring)
      real   (sp)                :: fstage_dur          ! fraction of length stage (herring)
      integer(ip)                :: istage              ! integer stage development
      integer(ip)                :: istage_t            ! index stage development
      integer(ip)                :: istag0              ! previous stage
      real   (sp)                :: stime               ! time since last stage change
      real   (sp)                :: stemp               ! average temperature since last stage change
      real   (sp)                :: dtemp               ! average temp day
      real   (sp)                :: ddepth              ! average depth day
      real   (sp)                :: sdepth              ! average depth since last stage change
      real   (sp)                :: a                   ! a coefficient in development (-)
      real   (sp)                :: b                   ! b coefficient in development (-)
      real   (sp)                :: temp                ! temperature
      real   (sp)                :: salinity            ! salinity
      real   (sp)                :: salprev             ! salinity previous
      logical                    :: eb_tide             ! eb_tide
      real   (sp)                :: delt                ! timestep in days
      real   (sp)                :: day                 ! time in days
      real   (sp)                :: duration            ! duration of current stage
!      real   (sp)                :: dur_20mm            ! duration of stage 2 herring up to 21mm
      real   (sp)                :: her_length          ! herring length in mm (stage 2)
      logical, save              :: is_nursery          ! if arived at nursery area
      real   (sp)                :: vz                  ! vz
      real   (sp)                :: vzact               ! vzact
      real   (sp)                :: buoy                ! buoy
      real   (sp)                :: mort                ! mortality
      real   (sp)                :: ztop                ! ztop
      real   (sp)                :: zbot                ! zbot
      integer                    :: m                   ! m
      integer                    :: n                   ! n
      integer                    :: k                   ! k
      integer                    :: kb                  ! k
      integer                    :: iseg                ! iseg
      integer                    :: isegl               ! isegl
      integer                    :: isegt               ! isegl
      integer                    :: isegt3               ! isegl
      integer                    :: isegb               ! isegb
      real   (sp)                :: z                   ! z
      real   (sp)                :: zdepth              ! z relative to water surface
      real   (sp)                :: zlevel              ! z relative to bottom
      real   (sp)                :: laydep              ! laydep
      real   (sp)                :: totdep              ! totdep
      real   (sp)                :: volseg              ! volseg
      real   (sp)                :: totvol              ! totvol for m2 output
      real   (sp)                :: totconc             ! total of particles over the watercolumn (for m2 output)
      logical, save              :: l_csv               ! csv output
      logical                    :: l_csv_now           ! csv output at this timestep
      logical,pointer, save      :: ebb_flow( : )       ! true if flow is ebb

      integer                    :: behaviour_type      ! actual behaviour type

      integer, parameter         :: behaviour_none     = 0 ! behaviour type none
      integer, parameter         :: behaviour_bottom   = 1 ! behaviour type bottom
      integer, parameter         :: behaviour_midlow   = 2 ! behaviour type midlow
      integer, parameter         :: behaviour_midtop   = 3 ! behaviour type midtop
      integer, parameter         :: behaviour_pelagic  = 4 ! behaviour type pelagic
      integer, parameter         :: behaviour_demersal = 5 ! behaviour type demersal
      integer, parameter         :: behaviour_diurnal  = 6 ! behaviour type diurnal
      integer, parameter         :: behaviour_herring  = 7 ! behaviour type herring
      integer, parameter         :: behaviour_stst_dem = 8 ! behaviour type stst_dem
      integer, parameter         :: behaviour_stst_pel = 9 ! behaviour type stst_pel
      integer, parameter         :: behaviour_float    = 10 ! behaviour type floating

      integer                    :: h_behaviour_type       ! horizontal behaviour

      integer,parameter          :: h_behaviour_none     = 0  ! no horizontal behaviour
      integer,parameter          :: h_behaviour_salinity = 1  ! horizontal behaviour towards lowest salinity

      real                       :: vswim                  ! swimming velocity
      real                       :: low_sal                ! lowest salinity
      integer                    :: n0
      integer                    :: n1
      integer                    :: n12
      integer                    :: n2
      integer                    :: n23
      integer                    :: n3
      integer                    :: n34
      integer                    :: n4
      integer                    :: n41
      integer                    :: n_low                  ! segment with lowest salinity
      real                       :: x_low                  ! x lowest salinity
      real                       :: y_low                  ! y lowest salinity
      real                       :: local_angle            ! angle towards lowest salinity in grid
      logical                    :: thd_n1                 ! thin dam towards n1
      logical                    :: thd_n2                 ! thin dam towards n2
      logical                    :: thd_n3                 ! thin dam towards n3
      logical                    :: thd_n4                 ! thin dam towards n4


      real   , parameter         :: pi = 3.141592654
      real   , parameter         :: twopi = pi*2.0
      integer, parameter         :: nfix               = 9 ! fixed number of constants
      integer, parameter         :: nvar               =20 ! variable number of constants per stage

      integer, save              :: ifirst = 1
      logical, save              :: l_larvae
      integer(4),save            :: ithndl = 0             ! handle to time this subroutine

      if ( itime .eq. itstrtp ) then

         if ( nocons .eq. 0 ) then
            write(lunrep,*) ' no constants, no larvae model activated'
            l_larvae = .false.
            return
         endif
         nstage         = nint(const(1))
         if ( nstage .le. 0 ) then
            write(lunrep,*) ' no stages, no larvae model activated'
            l_larvae = .false.
            return
         endif
         if ( nocons .ne. nfix+nstage*nvar ) then
            write(lunrep,*) ' no of constants inconsistent with stages, no larvae model activated'
            l_larvae = .false.
            return
         endif
         l_larvae       = .true.
         istage_nursery = nint(const(2))
         layer_release  = nint(const(3))
         it_start_m2    = nint(const(4)*86400.)
         it_stop_m2     = nint(const(5)*86400.)
         idt_m2         = nint(const(6))
         idt_csv        = nint(const(7))
         iniday         = nint(const(8))
         tide_opt       = nint(const(9))
         if ( idt_csv .le. 0 ) then
            l_csv = .false.
         else
            l_csv = .true.
         endif

         allocate(astage(nstage))
         allocate(bstage(nstage))
         allocate(length(nstage))
         allocate(btype(nstage))
         allocate(mort1(nstage))
         allocate(mort2(nstage))
         allocate(tcmort(nstage))
         allocate(buoy1(nstage))
         allocate(buoy2(nstage))
         allocate(vzact1(nstage))
         allocate(vzact2(nstage))
         allocate(ztop1(nstage))
         allocate(ztop2(nstage))
         allocate(zbot1(nstage))
         allocate(zbot2(nstage))
         allocate(phase(nstage))
         allocate(phase_diurn(nstage))
         allocate(hbtype(nstage))
         allocate(ampli1(nstage))
         allocate(ampli2(nstage))
         allocate(vswim1(nstage))
         allocate(vswim2(nstage))

         do istage = 1, nstage
            astage(istage) = const(nfix+(istage-1)*nvar+ 1)
            bstage(istage) = const(nfix+(istage-1)*nvar+ 2)
            btype(istage)  = const(nfix+(istage-1)*nvar+ 3)
            mort1(istage)  = const(nfix+(istage-1)*nvar+ 4)
            mort2(istage)  = const(nfix+(istage-1)*nvar+ 5)
            tcmort(istage) = const(nfix+(istage-1)*nvar+ 6)
            buoy1(istage)  =-const(nfix+(istage-1)*nvar+ 7)/86400.   !buoyancy with a minus, so positive buoyance is a negative settling units are m/d
            buoy2(istage)  =-const(nfix+(istage-1)*nvar+ 8)/86400.
            vzact1(istage) =-const(nfix+(istage-1)*nvar+ 9)/86400.
            vzact2(istage) =-const(nfix+(istage-1)*nvar+10)/86400.
            ztop1(istage)  = const(nfix+(istage-1)*nvar+11)
            ztop2(istage)  = const(nfix+(istage-1)*nvar+12)
            zbot1(istage)  = const(nfix+(istage-1)*nvar+13)
            zbot2(istage)  = const(nfix+(istage-1)*nvar+14)
            phase(istage)  = const(nfix+(istage-1)*nvar+15)
            ampli1(istage)  = const(nfix+(istage-1)*nvar+16)
            ampli2(istage)  = const(nfix+(istage-1)*nvar+17)
            hbtype(istage)  = const(nfix+(istage-1)*nvar+18)
            vswim1(istage)  = const(nfix+(istage-1)*nvar+19)
            vswim2(istage)  = const(nfix+(istage-1)*nvar+20)
         enddo
         length(1) = 3.0    !growth in mm in stage 1 for herring
         length(2) = 10.0   ! growth in mm in stage 2 for herring
         length(3) = 500.0
         allocate(ebb_flow(nosegl))

         ! set some stuff hard coded for the moment

         is_nursery     = .true.

      endif

      ! exit if no larvae model

      if ( .not. l_larvae ) return
      if ( timon ) call timstrt( "larvae", ithndl )

      ! global forcing

      delt     = idelt/86400.
      day      = itime/86400.
      do istage = 1, nstage
         phase_diurn(istage) = sin((day-phase(istage))*twopi)/2.+.5   !varies from 0-1
      enddo

      ! determine the tide for the complete grid

      call tidal_state( itime, itstrtp, nosegl, tide_opt, vol1, vol2, vel1, vel2, ebb_flow)

      ! output larvae per m2

      if ( itime .ge. it_start_m2 .and. itime .le. it_stop_m2 .and. mod(itime-it_start_m2,idt_m2) .lt. idelt ) then
         call larvm2 ( lunrep   , itime    , nosegl   , nolay    , nosubs   ,    &
                       conc     )
      endif

      ! every day output position to csv file

      l_csv_now = .false.
      if ( l_csv ) then
         if ( itime .ge. it_start_m2 .and. itime .le. it_stop_m2 .and. mod(itime-it_start_m2,idt_csv) .lt. idelt ) then
            l_csv_now = .true.
            iday = itime/idt_csv
            write(filcsv,'(a,i5.5,a)') 'position_day_',iday,'.csv'
            open(luncsv,file=filcsv)
            write(luncsv,*)'xa, ya, zdepth, sdepth, ddepth, stemp , her_length, stage, stime'
         endif
         ncum = ncum + 1
      endif
      dtemp  = 0.0
      ddepth = 0.0
      sdepth = 0.0
      her_length = 0.0
      do ipart = npwndw, nopart

         ! position and forcing

         m      = mpart(ipart)
         n      = npart(ipart)
         k      = min(kpart(ipart),nolay)
         kb     = kpart(ipart)
         z      = zpart(ipart)
         iseg   = lgrid3(n,m)  !or lgrid?
!         iseg   = lgrid(n,m) 
         if (iseg .gt. 0 ) then
         isegl  = iseg + (k-1)*nosegl
         isegb  = iseg + (nolay-1)*nosegl
         isegt3  = lgrid3(n,m)+(k-1)*nosegl
         isegt  = lgrid2(n,m)+(k-1)*nmax*mmax
         volseg = volume(isegt)

         laydep = locdep(lgrid2(n,m),k)
         if( k .ne. 1 ) laydep = laydep - locdep(lgrid2(n,m),k-1)
         totdep = locdep(lgrid2(n,m),nolay)
         zlevel = locdep(lgrid2(n,m),nolay) - locdep(lgrid2(n,m),k) + (1.-z)*laydep
         zdepth = locdep(lgrid2(n,m),k) - (1.-z)*laydep
         totconc =  locdep(lgrid2(n,m),k)*conc(nosubs+2,lgrid3(n,m))
         totvol  = volume (lgrid2(n,m))
         do ilay =2, nolay
           totconc = totconc +(locdep(lgrid2(n,m),ilay)-locdep(lgrid2(n,m),ilay-1))*conc(nosubs+2,lgrid3(n,m)+(ilay-1)*nosegl) !number of particles in the watercolumn
           totvol = totvol + volume (lgrid2(n,m)+(ilay-1)*nmax*mmax)
         enddo
         

         temp     = temper1(isegl)
         salinity = salin1(isegl)

         ! stage development

         stage  = wpart(2,ipart)
         stime  = wpart(3,ipart)
         stemp  = wpart(4,ipart)
!         if ( nosubs .ge. 4+nstage+1 ) dtemp  = wpart(4+nstage+1,ipart)
!         if ( nosubs .ge. 4+nstage+2 ) ddepth = wpart(4+nstage+2,ipart)
!         if ( nosubs .ge. 4+nstage+3 ) sdepth = wpart(4+nstage+3,ipart)
         istage = stage
         istage = max(istage,0)
         istag0 = istage

         stemp  = (stemp*stime + temp*delt)/(stime+delt)
 !        if ( nosubs .ge. 4+nstage+1 ) dtemp  = dtemp + temp
 !        if ( nosubs .ge. 4+nstage+2 ) ddepth = ddepth + zdepth
 !        if ( nosubs .ge. 4+nstage+3 ) sdepth = (sdepth*stime + zdepth*delt)/(stime+delt)
         stime  = stime+delt
         wpart(8,ipart)=0.0

         if ( istage .eq. 0 ) then !still at the bottom
            kpart(ipart) = nolay + 1
            if ( stime .gt. -(stage+delt) ) then
               stage  = 1.0
               istage = 1
               fstage = 0.0
               stime  = 0.0
               kpart(ipart) = layer_release
               zpart(ipart) = 0.5
               k            = layer_release
               kb           = layer_release
            endif
         else
            a = astage(istage)
            b = bstage(istage)
            behaviour_type = btype(istage)
            
            duration = exp(a+b*stemp) ! general duration curve
            if (behaviour_type.eq.7) then
              duration = length(istage)/(a * exp(b*stemp))
            endif
            
!            dur_20mm = 129.17*exp(-0.18*stemp) !in days, to be adapted for each stage
            if ( duration .lt. stime ) then
               istage = istage + 1
               istage = min(istage,5)
               stime  = 0.0
            endif
            fstage = stime/duration
            stage  = istage + fstage
            ! calculate actual herring length (after stage 0 length 7)
            select case (istage)
              case(1)
                her_length = 7.0 + fstage * length (istage)
     !           her_length = 12.5  ! for testing only
              case(2)
                her_length = 10.0 + fstage * length (istage)
     !           her_length = 12.5  ! for testing only
              case(3)
                her_length = 20.0 + fstage * length (istage)
     !           her_length = 12.5  ! for testing only
              case default
                her_length = 0.
            end select
  
              
!            fstage_dur = stime/dur_21mm ! for herring in stage 2
!            fstage_dur = iptime(ipart)/dur_20mm/86400. ! for herring in stage 2
!            stage_dur = istage + fstage_dur
!            stage_dur = fstage_dur
         endif
         wpart(2,ipart) = stage
         wpart(3,ipart) = stime
         wpart(4,ipart) = stemp

         ! mortality

         if ( istage .gt. 0 ) then
            mort = mort1(istage) + fstage*(mort2(istage)-mort1(istage))
            if ( abs(tcmort(istage)) .gt. 1.e-20 ) then
               mort = mort*exp(tcmort(istage)*temp)
            endif
            wpart(1,ipart) = wpart(1,ipart) - mort*wpart(1,ipart)*delt
         endif

         ! if there is substances per fraction set weight according to stage

         do istage_t = 1, nstage
            if ( istage_t+4 .le. nosubs ) then
               if ( istage_t .eq. istage ) then
                  wpart(istage_t+4,ipart) = wpart(1,ipart)
               else
                  wpart(istage_t+4,ipart) = 0.0
               endif
            endif
         enddo

         

         ! deactivate particles in nursery stage who have arrived in nursery area or if still 

         if ( istage .ge. istage_nursery .and. is_nursery ) then
            kpart(ipart) = nolay +1
         elseif ( istage .gt. 0 ) then

            ! release any particle from bottom

            if ( kb .gt. nolay ) then
               kpart(ipart) = nolay
               k            = nolay
               kb           = nolay
               zpart(ipart) = 0.5
            endif

            ! set layer and/or settling velocity according to stage and type of behaviour

            behaviour_type = btype(istage)
            wpart(8,ipart)=0.0
            select case ( behaviour_type )
               case ( behaviour_none )
                  wsettl(ipart) = 0.0
               case ( behaviour_bottom )
                  wsettl(ipart) = 0.0
                  kpart(ipart) = nolay + 1
                  zpart(ipart) = 0.5
               case ( behaviour_midlow )
                  wsettl(ipart) = 0.0
                  kpart(ipart) = nolay
                  zpart(ipart) = 0.5
               case ( behaviour_midtop )
                  wsettl(ipart) = 0.0
                  kpart(ipart) = 1
                  zpart(ipart) = 0.5
               case ( behaviour_pelagic )
                  wsettl(ipart) = buoy1(istage) + fstage*(buoy2(istage)-buoy1(istage))
               case ( behaviour_demersal )
                  ztop  = ztop1(istage)  + fstage*(ztop2(istage)-ztop1(istage))
                  zbot  = zbot1(istage)  + fstage*(zbot2(istage)-zbot1(istage))
                  buoy  = buoy1(istage)  + fstage*(buoy2(istage)-buoy1(istage))
                  vzact = vzact1(istage) + fstage*(vzact2(istage)-vzact1(istage))
                  if ( zlevel .gt. ztop ) then
                     vz = buoy
                  elseif ( zlevel .lt. zbot ) then
                     vz = vzact
                  else
                     vz = buoy + phase_diurn(istage)*(vzact-buoy)
                  endif
                  wsettl(ipart) = vz
               case ( behaviour_diurnal )
                  ztop  = ztop1(istage)  + fstage*(ztop2(istage)-ztop1(istage))
                  zbot  = zbot1(istage)  + fstage*(zbot2(istage)-zbot1(istage))
                  buoy  = buoy1(istage)  + fstage*(buoy2(istage)-buoy1(istage))
                  vzact = vzact1(istage) + fstage*(vzact2(istage)-vzact1(istage))
                  if ( zdepth .lt. ztop ) then
                     vz = buoy
                  elseif ( zlevel .lt. zbot ) then
                     vz = vzact
                  else
                     vz = buoy + phase_diurn(istage)*(vzact-buoy)
                  endif
                  wsettl(ipart) = vz
               case ( behaviour_herring )
                  ampli = ampli1(istage)+fstage*(ampli2(istage)-ampli1(istage))
                  ztop  = ztop1(istage)  + fstage*(ztop2(istage)-ztop1(istage))
                  zbot  = phase_diurn(istage)*ampli*2+(zbot1(istage)  + fstage*(zbot2(istage)-zbot1(istage))) !needs to be adapated for herring
                  buoy  = buoy1(istage)  + fstage*(buoy2(istage)-buoy1(istage))
                  vzact = vzact1(istage) + fstage*(vzact2(istage)-vzact1(istage))
                  if ( zdepth .lt. ztop ) then
                     vz = buoy
                  elseif ( zdepth .gt. ztop + zbot ) then
                     vz = vzact
                  else
                     vz = (buoy + vzact)/2.
                  endif
                  wsettl(ipart) = vz
!                  if (istage == 2 .and. fstage_dur <=1.0) then 
!                    her_length = 21.0 * max(0.0,(min(1.0,fstage_dur))) !estimated length; for testing use all stages (otherwise length underestimated due to particles in other stages
                  if ((laydep > 0) .and. (conc (nosubs+2,isegt3)>0)) then
!                  if (laydep > 0 ) then
                    wpart(8,ipart) = her_length * volseg /(totconc)  !herring length seg ! so that in the output it is the average length in that segment
                  else
                    wpart(8,ipart)=0.0
                  endif
                  
               case ( behaviour_stst_dem )
                  vzact = vzact1(istage) + fstage*(vzact2(istage)-vzact1(istage))
                  if ( ebb_flow(iseg) ) then

                     if ( k .ge. nolay ) then

                        ! settle on bed if arrived in lowest layer

                        wsettl(ipart) = 0.0
                        kpart(ipart) = nolay + 1
                        zpart(ipart) = 0.5

                     else

                        ! swim downwards

                        wsettl(ipart) = -vzact

                     endif

                  else

                     ! swim upwards

                     wsettl(ipart) = vzact

                  endif
               case ( behaviour_stst_pel )
                  vzact = vzact1(istage) + fstage*(vzact2(istage)-vzact1(istage))
                  if ( ebb_flow(iseg) ) then

                     ! swim downwards

                     wsettl(ipart) = -vzact

                  else

                     ! swim upwards

                     wsettl(ipart) = vzact

                  endif
               case ( behaviour_float )
               
                  !stay at the surface
                  
                  wsettl(ipart) = 0.0
                  kpart(ipart) = 1
                  zpart(ipart) = 0.001
                  
               case default
                  write(lunrep,*) ' error, larval behaviour type not defined'
                  call stop_exit(1)
            end select

            ! horizontal swimming

            h_behaviour_type = hbtype(istage)
            select case ( h_behaviour_type )
               case ( h_behaviour_none )
                  v_swim(ipart) = 0.0
                  d_swim(ipart) = 0.0
               case ( h_behaviour_salinity )
                  vswim = vswim1(istage) + fstage*(vswim2(istage)-vswim1(istage))
                  v_swim(ipart) = vswim

                  n0  = lgrid (n,m)
                  n1  = lgrid2(n-1,m)
                  n2  = lgrid2(n,m-1)
                  thd_n1 = ( flow(n1) .eq. 0.0 )
                  thd_n2 = ( flow(n2+mnmaxk) .eq. 0.0 )
                  thd_n3 = ( flow(n0) .eq. 0.0 )
                  thd_n4 = ( flow(n0+mnmaxk) .eq. 0.0 )

                  n0  = lgrid3(n,m)
                  n1  = lgrid3(n-1,m)
                  n12 = lgrid3(n-1,m-1)
                  n2  = lgrid3(n,m-1)
                  n23 = lgrid3(n+1,m-1)
                  n3  = lgrid3(n+1,m)
                  n34 = lgrid3(n+1,m+1)
                  n4  = lgrid3(n,m+1)
                  n41 = lgrid3(n-1,m+1)
                  low_sal = salin1(n0)
                  n_low = n0
                  if ( n1 .gt. 0 .and. .not. thd_n1 ) then
                     if ( salin1(n1) .lt. low_sal ) then
                        low_sal = salin1(n1)
                        n_low = n1
                        x_low = 0.5
                        y_low =-0.5
                     endif
                  endif
                  if ( n12 .gt. 0 .and. .not. ( thd_n1 .and. thd_n2 ) ) then
                     if ( salin1(n12) .lt. low_sal ) then
                        low_sal = salin1(n12)
                        n_low = n12
                        x_low =-0.5
                        y_low =-0.5
                     endif
                  endif
                  if ( n2 .gt. 0 .and. .not. thd_n2 ) then
                     if ( salin1(n2) .lt. low_sal ) then
                        low_sal = salin1(n2)
                        n_low = n2
                        x_low =-0.5
                        y_low = 0.5
                     endif
                  endif
                  if ( n23 .gt. 0 .and. .not. ( thd_n2 .and. thd_n3 ) ) then
                     if ( salin1(n23) .lt. low_sal ) then
                        low_sal = salin1(n23)
                        n_low = n23
                        x_low =-0.5
                        y_low = 1.5
                     endif
                  endif
                  if ( n3 .gt. 0 .and. .not. thd_n3 ) then
                     if ( salin1(n3) .lt. low_sal ) then
                        low_sal = salin1(n3)
                        n_low = n3
                        x_low = 0.5
                        y_low = 1.5
                     endif
                  endif
                  if ( n34 .gt. 0 .and. .not. ( thd_n3 .and. thd_n4 ) ) then
                     if ( salin1(n34) .lt. low_sal ) then
                        low_sal = salin1(n34)
                        n_low = n34
                        x_low = 1.5
                        y_low = 1.5
                     endif
                  endif
                  if ( n4 .gt. 0 .and. .not. thd_n4 ) then
                     if ( salin1(n4) .lt. low_sal ) then
                        low_sal = salin1(n4)
                        n_low = n4
                        x_low = 1.5
                        y_low = 0.5
                     endif
                  endif
                  if ( n41 .gt. 0 .and. .not. ( thd_n4 .and. thd_n1 ) ) then
                     if ( salin1(n41) .lt. low_sal ) then
                        low_sal = salin1(n41)
                        n_low = n41
                        x_low = 1.5
                        y_low =-0.5
                     endif
                  endif
                  if ( n_low .eq. n0 ) then
                     ! stay put, againts the flow with velocity of the flow


                     v_swim(ipart) = 0.0
                     local_angle   = 0.0

                  else

                     ! towards grid centre of n_low
                     a  = x_low - xpart(ipart)
                     b  = y_low - ypart(ipart)
                     local_angle = atan2(a,b)
                  endif

                  ! convert angle in local assumed rectangular grid towards angle in grid coordinate system, convert to degrees

                  local_angle = (local_angle - angle(n0))*360./twopi
                  if ( local_angle .lt. 0.0 ) local_angle = local_angle + 360.
                  if ( local_angle .gt. 360.0 ) local_angle = local_angle - 360.

                  d_swim(ipart) = local_angle

               case default
                  write(lunrep,*) ' error, horizontal behaviour type not defined'
                  call stop_exit(1)
            end select

         endif

         ! every day output position to csv file

         if ( l_csv_now ) then
            ddepth = ddepth/ncum
            dtemp  = dtemp/ncum
            write(luncsv,'(f10.2,'','',f10.2,7('','',f10.4))') xa(ipart),ya(ipart),zdepth,sdepth,ddepth,stemp,her_length,stage, stime
            ddepth = 0.0
            dtemp  = 0.0
         endif
!         if ( nosubs .ge. 4+nstage+1 ) wpart(4+nstage+1,ipart) = dtemp
!         if ( nosubs .ge. 4+nstage+2 ) wpart(4+nstage+2,ipart) = ddepth
!         if ( nosubs .ge. 4+nstage+3 ) wpart(4+nstage+3,ipart) = sdepth

         else
            ! inactive

         if ( l_csv_now ) then
            zdepth = -999.
            sdepth = -999.
            ddepth = -999.
            stemp  = -999.
            dtemp  = -999.
            write(luncsv,'(f10.2,'','',f10.2,5('','',f10.4),2('','',f10.1))') -999.,-999.,zdepth,sdepth,ddepth,stemp,-999.,-999.,-999.
         endif

         endif

         !

      enddo


      ! every day output position to csv file

      if ( l_csv_now ) then
         close(luncsv)
         ncum = 0
      endif

      continue

      if ( timon ) call timstop ( ithndl )
      return
      end subroutine
end module
