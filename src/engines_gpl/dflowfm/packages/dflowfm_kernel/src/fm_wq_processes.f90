!----- AGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2018-2019.!
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

! $Id$
! $HeadURL$
   subroutine fm_wq_processes_ini()
      use m_fm_wq_processes
      use m_alloc
      use unstruc_messages
      use m_flow, only: kmx, Lnkx
      use m_flowgeom, only: Ndxi, ba, Lnx, Lnxi, ln, lne2ln
      use m_flowparameters, only: jasal, jatem, jawave
      use m_flowexternalforcings
      use m_partitioninfo
      use unstruc_model
      use unstruc_files
      use m_flowtimes
      use m_wind, only: jawind, jarain
      
      implicit none
      
      character(20), allocatable :: syname_sub(:)       !< substance names from sub-file
      character(20), allocatable :: syunit_sub(:)       !< substance unit from sub-file
      character(20), allocatable :: coname_sub(:)       !< constant names from sub-file
      real         , allocatable :: covalue_sub(:)      !< values for contants from sub-file
      character(20), allocatable :: ouname_sub(:)       !< output names from sub-file
      character(80), allocatable :: oudesc_sub(:)       !< output decriptions from sub-file

      character(20), allocatable :: syname_eho(:)       !< substance names from sub-file
      character(20), allocatable :: syunit_eho(:)       !< substance names from sub-file
      character(20), allocatable :: coname_eho(:)       !< constant names from sub-file
      real         , allocatable :: covalue_eho(:)      !< values for contants from sub-file
      character(20), allocatable :: ouname_eho(:)       !< output names from sub-file
      character(80), allocatable :: oudesc_eho(:)       !< output decriptions from sub-file

      integer                    :: ierr_sub            !< error status
      integer                    :: ierr_eho            !< error status
      character(256)             :: cerr                !< error message

      !     for wq_processes_initialise
      character(len=256)       :: substance_file        !< substance file
      character(len=256)       :: his_output_file       !< extra history output file
      character(len=256)       :: proc_log_file         !< processes log file
      character(len=256)       :: proc_def_file         !< processes definition file
      character(len=256)       :: bloom_file            !< BLOOM algae spiecies paramter file
      
      type(procespropcoll)     :: statprocesdef   !< the statistical proces definition
      integer  ( 4), parameter :: nomult = 0      !< number of multiple substances
      integer  ( 4)            :: imultp(2,nomult)!< multiple substance administration
      type(t_dlwq_item)        :: constants       !< delwaq constants list
      integer                  :: noinfo          !< count of informative message
      integer                  :: nowarn          !< count of warnings
      integer                  :: ierr            !< error count

!     Other
      integer( 4)              :: lunlsp
      integer( 4)              :: ierr2, iresult
      integer( 4)              :: nerror
      integer( 4)              :: nrvarn
      integer( 4)              :: noout_sub
      integer( 4)              :: nosys_eho,notot_eho,nocons_eho
      integer( 4)              :: noout_eho
      integer( 4)              :: i, j, ip, ipar, ifun, isfun, ivar
      integer( 4)              :: iastat

      character*10,parameter   :: cbloom = 'd40blo'

      integer :: janew, itrac, iex
      integer :: kk, k, kb, kt, ktmax, kdum
      
      integer :: ierror
      
      logical :: Lsub, Lpdf, Lblm, Leho, opened, Lallocated, writebalance
      
      ierror = 1
      
      jawaqproc = 1
      ibflag = 0
      
      substance_file = md_subfile
      his_output_file = md_ehofile
      proc_def_file = md_pdffile
      bloom_file = md_blmfile
      proc_log_file = defaultfilename('lsp')
      
      open (newunit=lunlsp , file=proc_log_file, status='unknown', iostat=ierr)
      if (ierr.ne.0) then
         call mess(LEVEL_ERROR, 'Could not open processes log file: ', trim(proc_log_file))
      end if
      call setmlu(lunlsp)
      
!     check if sub- and pdffile have a value and exist
      inquire(file=substance_file,exist=Lsub)
      if ( .not.Lsub) then
         call mess(LEVEL_ERROR, 'Substance file does not exist: ', trim(substance_file))
      end if
      
      if (his_output_file.ne.' ') then
         inquire(file=his_output_file,exist=Leho)
         if ( .not.Leho) then
            call mess(LEVEL_ERROR, 'Extra histrory output file specified, but does not exist: ', trim(his_output_file))
         end if
      else
         Leho = .false.
      endif

      if (proc_def_file.ne.' ') then
         inquire(file=proc_def_file,exist=Lpdf)
         if ( .not.Lpdf) then
            call mess(LEVEL_ERROR, 'Process library file does not exist: ', trim(proc_def_file))
         end if
      else
         call mess(LEVEL_ERROR, 'No process library file specified. Use commandline argument --processlibrary "<path>/<name>"')
      endif   
      
      if (bloom_file.ne.' ') then
         inquire(file=bloom_file,exist=Lblm)
         if ( .not.Lblm) then
            call mess(LEVEL_ERROR, 'BLOOM species definition file specified, but does not exist: ', trim(bloom_file))
         end if
      else
         Lblm = .false.
      end if

!     water column definition
      if ( kmx.gt.0 ) then
         call getkbotktopmax(1,kbx,kt,kdum)
         call getkbotktopmax(Ndxi,kdum,kt,ktx)
         noseg = ktx-kbx+1 ! includes D-Flow FM dummy layer
         noq1 = 0
         noq2 = 0
         noq3 = noseg - 2*Ndxi + 1
         noq4 = 0
         
!        allocate vertical exchanges array
         call realloc(iexpnt, [4, noq3], keepExisting=.false., fill=0)
      
!        allocate array that indicates active cells (segments)
         call realloc(iknmrk, noseg, keepExisting=.false., fill=0)
         
!        allocate exchange to interface array
         call realloc(iex2k, noq3, keepExisting=.false., fill=0)
         
!        set vertical exchanges
         iex = 0
         do kk=1,Ndxi
            call getkbotktopmax(kk,kb,kt,ktmax)
            do k=ktmax,kb+1,-1
               iex = iex+1
               iexpnt(1,iex) = k - kbx+1
               iexpnt(2,iex) = k-1 - kbx+1
               
               iex2k(iex) = k-1
            end do
         end do
            
!        set array that indicates active cells (segments)
         iknmrk = 0
         do kk=1,Ndxi
            call getkbotktopmax(kk,kb,kt,ktmax)
            if ( kb.eq.ktmax ) then
               iknmrk(kb-kbx+1) = 1101
            else
               iknmrk(kb-kbx+1) = 1131
               do k=kb+1,ktmax-1
                  iknmrk(k-kbx+1) = 1121
               end do
               iknmrk(ktmax-kbx+1) = 1111
            end if
         end do
      else
         call mess(LEVEL_ERROR, 'WAQ: 2D not supported, set kmx=1')
      end if
   
! ======================
! Start initialising WAQ
! ======================

!     Read the substance file for the process defintion
!     Reset number of messages
      noinfo=0
      nowarn=0
      nerror=0
      ierr=0
      
      call mess(LEVEL_INFO, 'Opening substance file: ', trim(substance_file))
      Lallocated = .false.
      call rd_sub(Lallocated,substance_file,nosys,notot,nocons,noout_sub,syname_sub,syunit_sub,coname_sub, &
                  covalue_sub,ouname_sub,oudesc_sub,ierr_sub,cerr)
      if (ierr_sub.ne.0) call mess(LEVEL_ERROR, cerr)
      call realloc (syname_sub, notot, keepExisting=.false., fill=' ')
      call realloc (syunit_sub, notot, keepExisting=.false., fill=' ')
      call realloc (coname_sub, nocons, keepExisting=.false., fill=' ')
      call realloc (covalue_sub, nocons, keepExisting=.false., fill=0.0e0)
      call realloc (ouname_sub, noout_sub, keepExisting=.false., fill=' ')
      call realloc (oudesc_sub, noout_sub, keepExisting=.false., fill=' ')
      Lallocated = .true.
      call rd_sub(Lallocated,substance_file,nosys,notot,nocons,noout_sub,syname_sub,syunit_sub,coname_sub, &
                  covalue_sub,ouname_sub,oudesc_sub,ierr_sub,cerr)
      if (ierr_sub.ne.0) call mess(LEVEL_ERROR, cerr)
      noout_map = noout_sub

      if (Leho) then
         call mess(LEVEL_INFO, 'Opening extra history output file: ', trim(his_output_file))
         Lallocated = .false.
         call rd_sub(Lallocated,his_output_file,nosys_eho,notot_eho,nocons_eho,noout_eho,syname_eho,syunit_eho,coname_eho, &
                     covalue_eho,ouname_eho,oudesc_eho,ierr_eho,cerr)
         if (ierr_eho.ne.0) call mess(LEVEL_ERROR, cerr)
         if (nosys_eho==0 .and. notot_eho==0 .and. nocons_eho==0) then
            call realloc (ouname_eho, noout_eho, keepExisting=.false., fill=' ')
            call realloc (oudesc_eho, noout_eho, keepExisting=.false., fill=' ')
            Lallocated = .true.
            call rd_sub(Lallocated,his_output_file,nosys_eho,notot_eho,nocons_eho,noout_eho,syname_eho,syunit_eho,coname_eho, &
                        covalue_eho,ouname_eho,oudesc_eho,ierr_eho,cerr)
            if (ierr_eho.ne.0) call mess(LEVEL_ERROR, cerr)
         else
            ! Error: extra history outputfile contains other definitions that output alone!
         endif
      else
         ! Warning: extra history outputfile does not exist!
         noout_eho = 0
      end if

      
      noout = noout_map + noout_eho

!     The active substances should be initialised as 'constituents' in DFM.
!     Initial concentration (fields), boundary conditions and additional (waste) loads
!     should be specified in DFM
      call realloc(syname, notot)
      call realloc(syunit, notot)
      do i = 1, notot
          syname(i) = syname_sub(i)
          syunit(i) = syunit_sub(i)
      end do
      
      call realloc(amass, [notot, noseg], keepExisting=.false., fill=0.0d0)       !< mass array to be updated

!     Additional  data that comes from DFM should be added to the parameter/function/segment function list before the wq_processes_initialise call

!     No spatial parameters for now, they should come from DFM
      nopa = 0
      call realloc(paname, nopa)
      
!      Use functions to set 2D (or 0D variables) from DFM per column
      nofun = 0
      allocate(funame(nofun))
      call realloc(funame, nofun)
      
      call dfm_waq_initexternalforcings(ierr)
      if (ierr.ne.0) then
         call mess(LEVEL_ERROR, 'Error reading water quality processes external forcings from ext-file')
      endif

!     Use segment functions to set 2D (or 0D variables) from DFM per column (e.g. salinity or temperature)
      nosfun = 3
      call realloc(sfname, nosfun)
      sfname(1) = 'SURF                '
      sfname(2) = 'tau                 '
      sfname(3) = 'velocity            '
      
      i=3
      if ( jasal.eq.1 ) then
        nosfun = nosfun+1
        call realloc(sfname, nosfun, keepExisting=.true., fill='Salinity')
        isfsal = nosfun
      else
        isfsal = 0
      end if
      
      if ( jatem.gt.0 ) then
        nosfun = nosfun+1
        call realloc(sfname, nosfun, keepExisting=.true., fill='Temperature')
        isftem = nosfun
      else
        isftem = 0
      end if
      
      if ( jawind.eq.1 ) then
        nosfun = nosfun+1
        call realloc(sfname, nosfun, keepExisting=.true., fill='VWind')
        isfvwind = nosfun
      else
         isfvwind = 0
      end if
      
      if ( jawave.eq.1 .or. jawave.eq.2 ) then  ! copied from "flow_setexternalforcings", call to "tauwavefetch"
         nosfun = nosfun+1
         call realloc(sfname, nosfun, keepExisting=.true., fill='Fetch')
         isffetch = nosfun
      else
         isffetch = 0
      end if
      
      if ( jatem.gt.1 ) then  ! copied from "heatun"
         nosfun = nosfun+1
         call realloc(sfname, nosfun, keepExisting=.true., fill='RadSurf')
         isfradsurf = nosfun
      else
         isfradsurf = 0
      end if
      
      if ( jarain.gt.0 ) then
         nosfun = nosfun+1
         call realloc(sfname, nosfun, keepExisting=.true., fill='Rain')
         isfrain = nosfun
      else
         isfrain = 0
      end if

      noconm = nocons + 1000
      call realloc(coname, noconm)
      ierr2 = dlwq_init_item(constants)
      ierr2 = dlwq_resize(constants,noconm)
      j=0

!     Skip constants from the sub-file that are will be added by DFM
      do i = 1, nocons
         call zoekns(coname_sub(i),nopa,paname,20,ipar)
         call zoekns(coname_sub(i),nofun,funame,20,ifun)
         call zoekns(coname_sub(i),nosfun,sfname,20,isfun)
         if (ipar<0 .and. ifun<0 .and. isfun<0) then
            j = j + 1
            coname(j) = coname_sub(i)
            constants%ipnt(j) = i
            constants%name(j) = coname_sub(i)
            constants%constant(j) = covalue_sub(i)
         end if
      end do
      nocons = j
      constants%no_item = nocons

!     Set the required output data.
!     -> When outputs are not specified, they might end up in A(1), a 'black hole' location that is constantly overwritten
      allocate(outputs%names(noout))
      allocate(outputs%stdnames(noout))
      allocate(outputs%pointers(noout))
      allocate(outputs%units(noout))
      allocate(outputs%descrs(noout))
      
      call realloc(waqoutputs, [noout, noseg], keepExisting=.false., fill = 0.0d0)
      
      outputs%cursize  = noout
      do i = 1, noout_sub
          outputs%names(i) = ouname_sub(i)
          outputs%stdnames(i) = ' '
          outputs%units(i) = ' '
          outputs%descrs(i) = trim(oudesc_sub(i))//' ('//trim(ouname_sub(i))//') '
          outputs%pointers(i) = -1
      enddo
      do j = 1, noout_eho
          i = noout_map+j
          outputs%names(i) = ouname_eho(j)
          outputs%stdnames(i) = ' '
          outputs%units(i) = ' '
          outputs%descrs(i) = trim(oudesc_eho(j))//' ('//trim(ouname_eho(j))//') '
          outputs%pointers(i) = -1
      enddo

      deallocate (syname_sub)
      deallocate (coname_sub)
      deallocate (covalue_sub)
      deallocate (ouname_sub)

!     calculation timers need to be known for the statistical processes (start time/stop time)
      itfact = 86400              !< time scale factor processes
      itstrt = nint(tstart_user)
      itstop = nint(tstop_user)
      
!     for now, disable statistical processes
      statprocesdef%cursize = 0

!     Finally, evaluate the processes using the proces library
!     --------------------------------------------------------
      call wq_processes_initialise ( lunlsp, proc_def_file, bloom_file, statprocesdef, outputs, &
                                     nomult, imultp, constants, noinfo, nowarn, ierr)

!     proces fractional step multiplier is 1 for all
      prondt = 1

!     Allocate the work arrays for the pointers
      call realloc(ipmsa, nipmsa, keepExisting=.false., fill=0)
      call realloc(increm,nipmsa, keepExisting=.false., fill=0)
      
!     allocate flux and deriv arrays
      call realloc(flux, [nflux, noseg], keepExisting=.false., fill=0.0 )       !< Proces fluxes
      call realloc(deriv, [notot, noseg], keepExisting=.false., fill=0.0 )      !< Model derivatives (= stochi(notot ,noflux) * flux(noflux, noseg))
      call realloc(velonw, [nveln, noq3], keepExisting=.false., fill=0.0 )      !< New velocity array   

!     Determine size of a array from process system and noseg/noq3, and allocate it
      call wq_processes_pmsa_size( lunlsp, noseg, noq3, sizepmsa )
!     And actually allocate and zero the A array
      call realloc(pmsa, sizepmsa , keepExisting=.false., fill=0.0)

!     constants from the substance file
      ip = arrpoi(iicons)
      do i = 1, nocons
          pmsa(ip+i-1) = constants%constant(i)
      end do

!     defaults from the proces library
      ip = arrpoi(iidefa)
      do i = 1, nodef
          pmsa(ip+i-1) = defaul(i)
      end do

!     parameters from ext-file
!     spatially varying constants provided trough the ext-file that remain fixed during the run
!      -> fill using abcdabcdabcd partern
      if (nopa.gt.0) then
         ip = arrpoi(iiparm)
         do k=0,ktx-kbx
            do j = 1,nopa
               pmsa(ip) = painp(j,k+kbx)
               ip = ip + 1
            end do
         end do
      end if            

!      functions from DFM, one value for the whole system
!      -> can be used for timeseries that are updated by DFM, e.g. global irradiation
!      -> or can be used for 2D values when calculating only one water column
!      -> updated in waqfil

!     exchange areas
      ip = arrpoi(iiarea)
      iex = 0
      do kk=1,Ndxi
         call getkbotktopmax(kk,kb,kt,ktmax)
         do k=ktmax,kb+1,-1
            pmsa(ip+iex) = ba(kk)
            iex = iex+1
         end do
      end do

      call zoekns(cbloom,nproc,pronam,10,ipbloo)
      if (ipbloo.gt.0) then
          ioffbl = prvpnt(ipbloo)
          write ( lunlsp, * ) ' MESSAGE: Bloom fractional step switched on'
      else
          ipbloo = 0
          ioffbl = 0
      endif

!     add corresponding tracers, if not already defined by initial and/or boundary conditions
      transformcoef = 0.0_hp
      call realloc(isys2trac,notot,keepExisting=.false.,fill=0)
      do i=1,notot
         call add_bndtracer(trim(syname(i)), syunit(i), isys2trac(i), janew)
      end do
!      write(lunrep,'(/a)') '  var number var name               var array   var index'
!      write(lunrep,'(i,x,a20,2i)') ((i,varnam(i),vararr(i),varidx(i)),i=1,novar)
!     lookup output variables in the variables list (output%pointers is to complex)
      noout = outputs%cursize   
      call realloc(outvar,noout,keepExisting=.false.,fill=0)
      do j=1,noout
          call zoekns(outputs%names(j),novar,varnam,20,ivar)
          outvar(j) = ivar
      enddo  

      ierror = 0
   
!     error handling   
 1234 continue

      return
   end subroutine fm_wq_processes_ini

 !! @return Integer result status (0 if successful)
   subroutine dfm_waq_initexternalforcings(iresult)
   use unstruc_boundaries
   use m_alloc
   use m_flowexternalforcings
   use m_flowparameters
   use m_flowtimes
   use m_flowgeom
   use m_partitioninfo
   use unstruc_model
   use unstruc_messages
   use timespace
   use m_flow
   use m_meteo
   use m_ec_instance
   use dfm_error
   use processes_input
   use m_fm_wq_processes

   implicit none
   integer, intent (out)         :: iresult

   character(len=256)            :: filename, sourcemask
   integer                       :: kb, k, ja, method, kk, kt, lenqidnam, ipa, ifun, imba, imna
   character (len=NAMTRACLEN)    :: qidnam
   character (len=20)            :: waqinput
   integer                       :: minp0, npli, inside, filetype0, iad, needextramba, needextrambar
   double precision, allocatable :: viuh(:)            ! temporary variable
   integer, external             :: findname

   iresult = DFM_NOERR

   success = .true.    ! default if no valid providers are present in *.ext file (m_flowexternalforcings::success)

   call settimespacerefdat(refdat, julrefdat, Tzone, Timjan)

   ! initialise mass balance areas
   call realloc(mbadef, Ndkx, keepExisting=.false., fill =-999)
   call realloc(mbadefdomain, Ndkx, keepExisting=.false., fill =-999)

   if (mext /= 0) then 
      ja = 1

      do while (ja .eq. 1)                                ! read *.ext file
 
         call delpol()                                    ! remove a possibly existing polygon
         call readprovider(mext,qid,filename,filetype,method,operand,transformcoef,ja,sourcemask)
         if (ja == 1) then
            call mess(LEVEL_INFO, 'External Forcing or Initialising '''//trim(qid)//''' from file '''//trim(filename)//'''.')
            ! Initialize success to be .false.
            success = .false.

            qidnam = qid
            call get_waqinputname(qid, waqinput, qidnam)
            lenqidnam = len_trim(qidnam)
            if (filetype == 7 .and. method == 4) then
               method = 5                                   ! upward compatible fix
            endif
      
            if (qid(1:12) == 'waqparameter') then
               ipa = findname(nopa, paname, waqinput)
      
               if ( ipa.eq.0 ) then
                  nopa = nopa + 1
                  ipa = nopa
                  call realloc(paname, nopa, keepExisting=.true., fill=waqinput)
                  call realloc(painp, [nopa, Ndkx], keepExisting=.true., fill=0.0)
               end if
               call realloc(viuh, Ndkx, keepExisting=.false., fill=dmiss)
      
!              copy existing parameter values (if they existed) in temp array
               do kk=1,Ndxi
                  call getkbotktop(kk,kb,kt)
                  viuh(kk) = painp(ipa,kk)
                  do k=kb,kb+kmxn(kk)-1
                     viuh(k) = painp(ipa,k)
                  end do
               end do
      
!              will only fill 2D part of viuh          
               success = timespaceinitialfield(xz, yz, viuh, Ndx, filename, filetype, method, operand, transformcoef, 2)
      
               if (success) then
                  do kk = 1,Ndxi
                     if (viuh(kk) .ne. dmiss) then
                        painp(ipa,kk) = viuh(kk)
                        call getkbotktop(kk,kb,kt)
                        do k=kb,kb+kmxn(kk)-1
                           painp(ipa,k) = painp(ipa,kk)
                        end do
                     endif
                  enddo
               endif
               deallocate(viuh)
      
            else if (qid(1:11) == 'waqfunction') then
      
               ifun = findname(nofun, funame, waqinput)
      
               if ( ifun.eq.0 ) then
                  nofun = nofun + 1
                  ifun = ifun
                  call realloc(funame, nofun, keepExisting=.true., fill=waqinput)
                  call reallocP(funinp, [nofun, 1], keepExisting=.true., fill=0.0d0)
               end if
      
            else if (qid(1:18) == 'waqmassbalancearea') then
               imba = findname(nomba, mbaname, waqinput)
      
               if ( imba.eq.0 ) then
                  nomba = nomba + 1
                  imba = nomba
                  call realloc(mbaname,nomba,keepExisting=.true.,fill=waqinput)
               end if
               call realloc(viuh,Ndkx,keepExisting=.false.,Fill=dmiss)
      
!              will only fill 2D part of viuh          
               success = timespaceinitialfield(xz, yz, viuh, Ndx, filename, filetype, method, operand, transformcoef, 2)
      
               if (success) then
                  do kk=1,Ndxi
                     if (viuh(kk).ne.dmiss) then
                        if (mbadef(kk).ne. -999) then
                           ! warn that segment nn at xx, yy is nog mon area imba
                        endif
                        mbadef(kk) = imba
                        call getkbotktop(kk,kb,kt)
                        do k=kb,kb+kmxn(kk)-1
                           mbadef(k) = imba
                        end do
                     endif
                  end do
               endif
               deallocate(viuh)
      
            else if (qid(1:17) == 'waqmonitoringarea') then
               imna = findname(nomon, monname, waqinput)
      
               if ( imna.eq.0 ) then
                  nomon = nomon + 1
                  imna = nomon
                  call realloc(monname, nomon, keepExisting=.true., fill=waqinput)
                  call realloc(mondef, [nomon, Ndkx], keepExisting=.true., fill=2)
               end if
      
               call realloc(viuh, Ndkx, keepExisting=.false., fill=dmiss)
      
!              will only fill 2D part of viuh          
               success = timespaceinitialfield(xz, yz, viuh, Ndx, filename, filetype, method, operand, transformcoef, 2)
      
               if (success) then
                  do kk=1,Ndxi
                     if (viuh(kk).ne.dmiss) then
                        mondef(imna,kk) = 1
                        call getkbotktop(kk,kb,kt)
                        do k=kb,kb+kmxn(kk)-1
                           mondef(imna,k) = 1
                        end do
                     endif
                  end do
               endif
               deallocate(viuh)
            else
!              just accept any other keyword as success, they are evaluated again in unstruc.F90
               success = .true.
            endif
      
        endif
      
      enddo
 
   endif ! read mext file 

!  Check if there are any cells left that are not part of a mass balance area, and if we need an extra area.
   needextramba = 0
   do kk=1,Ndxi
      if (mbadef(kk).eq.-999) then
         needextramba = 1
         exit
      endif
   end do
   
   if (jampi.eq.1) then
!     check this among all domains (it could be that there are no remaing cels in this domain, while there are in other domains).
      call reduce_int_sum(needextramba, needextrambar)
      needextramba = needextrambar
   endif
   
   if(needextramba.ne.0) then
!     add the extra mass balance area, and assing the undefined cells to this area.
      nomba = nomba + 1
      call realloc(mbaname,nomba,keepExisting=.true.,fill="Remaining cells")
      imba = nomba
      do kk=1,Ndxi
         if (mbadef(kk).eq.-999) then
            mbadef(kk) = imba
            call getkbotktop(kk,kb,kt)
            do k=kb,kb+kmxn(kk)-1
               mbadef(k) = imba
            end do
         endif
      end do
   endif
   
   do kk=1,Ndxi
      if ( jampi.eq.1 ) then
!        do not include ghost cells
         if ( idomain(kk).ne.my_rank ) cycle
      end if
      mbadefdomain(kk) = mbadef(kk)
      call getkbotktop(kk,kb,kt)
      do k=kb,kb+kmxn(kk)-1
            mbadefdomain(k) = mbadef(k)
      end do
   end do

   if (loglevel_StdOut == LEVEL_DEBUG) then
      call ecInstancePrintState(ecInstancePtr,callback_msg,LEVEL_DEBUG)
   endif

   if (.not. success) then
      iresult = DFM_EXTFORCERROR
   end if
    
   if (mext /= 0) then
      rewind(mext) ! rewind ext file
   end if
   end subroutine dfm_waq_initexternalforcings


   
!> Convert qid (from .ext file) to waq input name (split in generic qidname and specific input name).
!! If the input qid is not waq input name, then the same qid is returned (and no waq input name)
   subroutine get_waqinputname(qid, inputname, qidname)
      implicit none
      
      character(len=*), intent(in)  :: qid       !< Original quantityid, e.g., 'waqfunctionradsurf'.
      character(len=*), intent(out) :: inputname !< The trimmed waq input name, e.g., 'fluor'.
      character(len=*), intent(out) :: qidname   !< The base input name for further use in external file analisys, e.g., 'tracerbnd'.
      
      inputname = ''
      if ( qid(1:11).eq.'waqfunction' ) then
         qidname = qid(1:11)
         if ( len_trim(qid).gt.11 ) then
            inputname = trim(qid(12:))
         end if
      else if (qid(1:12).eq.'waqparameter' ) then
         qidname = qid(1:12)
         if ( len_trim(qid).gt.12 ) then
            inputname = trim(qid(13:))
         end if
      else if (qid(1:18).eq.'waqmassbalancearea' ) then
         qidname = qid(1:18)
         if ( len_trim(qid).gt.18 ) then
            inputname = trim(qid(19:))
         end if
      else if (qid(1:17).eq.'waqmonitoringarea' ) then
         qidname = qid(1:17)
         if ( len_trim(qid).gt.17 ) then
            inputname = trim(qid(18:))
         end if
      end if
      
      return
   end subroutine get_waqinputname

   subroutine fm_wq_processes_step(dt,time)
      use m_fm_wq_processes
      use m_missing, only: dmiss
      use unstruc_model, only: md_flux_int
      implicit none
      
      double precision, intent(in) :: dt   !< timestep for waq in seconds
      double precision, intent(in) :: time !< time     for waq in seconds
      
      integer                      :: ipoiconc
      integer                      :: i, j, ip
                                   
      integer                      :: ipoivol, ipoisurf, ipoiarea
      integer                      :: ipoivelx, ipoidefa
      
      integer                      :: ivar, iarr, iv_idx
      integer                      :: iarknd, ip_arr, idim1, idim2
      integer                      :: incr

      integer                      :: idt, itime
                                   
      integer                      :: ierr
      
      double precision             :: dti
      
      integer                      :: ipvol, isys, k
      
      if ( jawaqproc.ne.1 ) return
      flux_int = md_flux_int
      
!     copy data from D-FlowFM to WAQ 
      call copy_data_from_fm_to_wq_processes(time)
      
      ipoiconc = arrpoi(iiconc)
      ipoivol  = arrpoi(iivol)
      ipoivelx = arrpoi(iivelx)
      ipoidefa = arrpoi(iidefa)
      ipoisurf = arrpoi(iisfun)
      ipoiarea = arrpoi(iiarea)
      
      idt   = int(dt)
      itime = int(time)
      
      call wq_processes_proces (notot , noseg , pmsa(ipoiconc), pmsa(ipoivol) , itime , idt   , deriv , ndmpar, &
                                nproc , nflux , ipmsa , prvnio, promnr, iflux , increm, flux  , flxdmp, stochi, &
                                ibflag, ipbloo, ioffbl, amass , nosys , itfact, intopt, iexpnt, iknmrk, noq1  , &
                                noq2  , noq3  , noq4  , pmsa(ipoiarea), ndspn , idpnew, dispnw, ndspx , dspx  , &
                                dsto  , nveln , ivpnw , velonw, nvelx , pmsa(ipoivelx), vsto  , mbadefdomain(kbx:ktx), &
                                pmsa(ipoidefa), prondt, prvvar, prvtyp, vararr, varidx, arrpoi, arrknd, arrdm1, &
                                arrdm2, novar , pmsa  , nomba , pronam, prvpnt, nodef , pmsa(ipoisurf), flux_int )

!     copy data from WAQ to D-FlowFM
      call copy_data_from_wq_processes_to_fm()
      return
   end subroutine fm_wq_processes_step
  
!
!  copy data from D-FlowFM to WAQ 
!
   subroutine copy_data_from_fm_to_wq_processes(time) 
      use m_flowgeom,       only: Ndxi, ba
      use m_flow,           only: vol1, sa1, tem1, ucx, ucy
      use m_fm_wq_processes           
      use m_transport,      only: itrac2const, constituents
      use m_sferic,         only: twopi
      use m_wind  
      use m_meteo
      use processes_input
      use m_waves,          only: fetch, nwf
      use unstruc_messages
      implicit none

      double precision, intent(in) :: time !< time     for waq in seconds

      double precision :: taucurc, czc
      double precision :: u10, dir
      double precision :: alfa1, alfa2
      double precision :: qsu
      
      integer          :: isys, iconst
      integer          :: ipoisurf, ipoitau, ipoivel
      integer          :: ipoivol, ipoiconc, ipoisal, ipoitem
      integer          :: ipoivwind, ipoifetch, ipoiradsurf, ipoirain
      integer          :: i, iex, ip, ifun
      integer          :: kk, k, kb, kt, ktmax, kwaq
      integer          :: L, nw1, nw2
                       
      integer          :: iknmrk_dry, iknmrk_wet
      logical, save    :: first = .true.
      
      if (nofun>0) then
         do ifun=1,nofun
            success = ec_gettimespacevalue(ecInstancePtr, item_waqfun(ifun), time)
         end do
         ip = arrpoi(iifunc)
         do ifun=1,nofun
            pmsa(ip+ifun-1) = funinp(ifun,1)
         end do
      end if

      ipoisurf = arrpoi(iisfun)           ! surface is first segment function
      ipoitau  = ipoisurf + noseg
      ipoivel  = ipoitau  + noseg
      
      ipoivol  = arrpoi(iivol)
      
      ipoiconc = arrpoi(iiconc)
      
      do kk=1,Ndxi
         call getkbotktop(kk,kb,kt)
         do k=kb,kt
            pmsa(ipoisurf + k-kbx) = ba(kk)
            pmsa(ipoivel  + k-kbx) = sqrt(ucx(k)**2 + ucy(k)**2)
         end do
         
         call gettau(kk,taucurc,czc)
         pmsa(ipoitau+kb-kbx) = taucurc
      end do
      
      do k=0,ktx-kbx
         pmsa(ipoivol + k) = vol1(k+kbx)
      end do
      
      if ( isfsal.gt.0 ) then
         ipoisal = arrpoi(iisfun) + (isfsal-1)*noseg
         do k=0,ktx-kbx
            pmsa(ipoisal + k) = sa1(k+kbx)
         end do
      end if
      
      if ( isftem.gt.0 ) then
         ipoitem = arrpoi(iisfun) + (isftem-1)*noseg
         do k=0,ktx-kbx
            pmsa(ipoitem + k) = tem1(k+kbx)
         end do
      end if
      
!     copy 2D arrays for wind velocity magnitude, fetch length, solar radiation and rain to 3D waq arrays, fill over whole column (safety)
      if ( isfvwind.gt.0 ) then
         ipoivwind = arrpoi(iisfun) + (isfvwind-1)*noseg
      end if
      
      if ( isffetch.gt.0 ) then   ! note: no fetch without wind
         ipoifetch = arrpoi(iisfun) + (isffetch-1)*noseg
      end if
      
      if ( isfradsurf.gt.0 ) then
         ipoiradsurf = arrpoi(iisfun) + (isfradsurf-1)*noseg
      end if
      
      if ( isfrain.gt.0 ) then
         ipoirain = arrpoi(iisfun) + (isfrain-1)*noseg
      end if
         
         
         
      do kk=1,Ndxi
         call getkbotktop(kk,kb,kt)
         
         if ( isfvwind.gt.0 ) then
!           copied from "tauwavefetch"
            call getlink1(k,L)   ! see comment in "tauwavefetch"
            u10 = sqrt( wx(L)*wx(L) + wy(L)*wy(L) )
         
            pmsa(ipoivwind + kb-kbx : ipoivwind + kt-kbx) = u10
         end if
         
         if ( isffetch.gt.0 ) then
!           copied from "tauwavefetch"
            if (u10.lt.1)  then
               dir = atan2(wy(L),wx(L))
               if ( dir.lt.0d0 ) dir = dir + twopi
               dir = dir/twopi
               if ( dir>=1d0 ) dir = 0d0
               
               nw1    = dir*(nwf-1) + 1
               nw2    = nw1 + 1
               alfa2  = (nwf-1)*( dir - dble(nw1-1) / dble(nwf-1) )
               alfa1  = 1d0 - alfa2
               
               pmsa(ipoifetch + kb-kbx : ipoifetch + kt-kbx) = alfa1*fetch(nw1,kk) + alfa2*fetch(nw2,kk)
            else
               pmsa(ipoifetch + kb-kbx : ipoifetch + kt-kbx) = 0d0
            end if
         end if
         
         if ( isfradsurf.gt.0 ) then
            qsu = 0d0
            
!           TODO: check and copy from "heatun"
            
!            pmsa(ipoiradsurf + kb-kbx : ipoiradsurf + kt-kbx) = qsu
         end if
         
         if ( isfrain.gt.0 ) then
            
!           rain: mm/day            
!            pmsa(ipoirain + kb-kbx : ipoirain + kt-kbx ) = rain(kk)
         end if
      end do
      
     
      
!     fill concentrations   
      do k=kbx,ktx
         do isys=1,notot
            iconst = isys2const(isys)
            pmsa(ipoiconc+(k-kbx)*(notot)+isys-1) = constituents(iconst,k)
         end do
      end do

!     fill masses (transported)
      do k=kbx,ktx
         if (vol1(k).gt.waq_vol_dry_thr) then
            do isys=1,nosys
               iconst = isys2const(isys)
               amass(isys,k-kbx+1) = constituents(iconst,k)*vol1(k)
            end do
         else
            do isys=1,nosys
               iconst = isys2const(isys)
               amass(isys,k-kbx+1) = 0.0d0
            end do
         endif
      end do
      
!     fill masses (not transported, only first time)
      if (notot>nosys.and.first) then
         first = .false.
         do kk=1,Ndxi
            call getkbotktopmax(kk,kb,kt,ktmax)
            do isys=nosys+1,notot
               iconst = isys2const(isys)
               amass(isys,kb-kbx+1) = constituents(iconst,kb)*ba(kk)
            end do
         end do
      end if

!     set dry/wet indicator
      do kk=1,Ndxi
         call getkbotktopmax(kk,kb,kt,ktmax)
         do k=kb,kt
            kwaq = k-kbx+1
            
            iknmrk_dry = int(iknmrk(kwaq)/10) * 10
            iknmrk_wet = iknmrk_dry + 1
            
            if ( vol1(k).gt.waq_vol_dry_thr ) then
               iknmrk(kwaq) = iknmrk_wet
            else
               iknmrk(kwaq) = iknmrk_dry
            end if
         end do
         
         do k=kt+1,ktmax
            kwaq = k-kbx+1
            
            iknmrk_dry = int(iknmrk(kwaq)/10) * 10
            iknmrk_wet = iknmrk_dry + 1
            
            if ( vol1(k).gt.waq_vol_dry_thr ) then
               iknmrk(kwaq) = iknmrk_wet
            else
               iknmrk(kwaq) = iknmrk_dry
            end if
         end do
      end do
      
      return
   end subroutine copy_data_from_fm_to_wq_processes

!
!  copy data from WAQ to D-FlowFM
!
   subroutine copy_data_from_wq_processes_to_fm() 
      use m_missing,        only: dmiss
      use m_flowgeom,       only: Ndxi, ba
      use m_flow,           only: vol1
      use m_fm_wq_processes           
      use m_transport,      only: itrac2const, constituents
      use m_sferic,         only: twopi
      use m_wind            
      use m_waves,          only: fetch, nwf
      use unstruc_messages
      implicit none
      
      double precision :: taucurc, czc
      double precision :: u10, dir
      double precision :: alfa1, alfa2
      double precision :: qsu
      
      integer          :: isys, iconst
      integer          :: ipoiconc
      integer          :: ivar, iarr, iv_idx
      integer          :: iarknd, ip_arr, idim1, idim2
      integer          :: incr
      integer          :: i, j, ip
      integer          :: kk, k, kb, kt, ktmax                 
      
!     fill concentrations (transported)
      do kk=1,Ndxi
         call getkbotktopmax(kk,kb,kt,ktmax)
         do k=kb,kt
            if (vol1(k).gt.waq_vol_dry_thr) then
               do isys=1,nosys
                  iconst = isys2const(isys)
                  constituents(iconst,k) = amass(isys,k-kbx+1) / vol1(k)
               end do
            end if
         end do
      end do
      
!     fill concentrations (not transported)
      if (notot>nosys) then
         do kk=1,Ndxi
            call getkbotktopmax(kk,kb,kt,ktmax)
            do isys=nosys+1,notot
               iconst = isys2const(isys)
               constituents(iconst,kb) = amass(isys,kb-kbx+1) / ba(kk)
            end do
         end do
      end if
      
! Ouputs to waq outputs array
      waqoutputs=dmiss
      noout = outputs%cursize   
      do j = 1, noout   
         ivar   = outvar(j)  ! which variable is it
         if (ivar > 0) then
            iarr   = vararr(ivar)         ! which array in pmsa
            iv_idx = varidx(ivar)         ! which index within the array
            iarknd = arrknd(iarr)         ! which type of array (increm is 0, dim1 or 1)
            ip_arr = arrpoi(iarr)         ! start point of the array in pmsa
            idim1  = arrdm1(iarr)         ! dimension in the 1e direction
            idim2  = arrdm2(iarr)         ! dimension in the 2e direction
            if ( iarknd .eq. 1 ) then
               ip = ip_arr + iv_idx - 1
               incr = 0
            elseif ( iarknd .eq. 2 ) then
               ip = ip_arr + iv_idx - 1
               incr = idim1
            elseif ( iarknd .eq. 3 ) then
               ip = ip_arr + (iv_idx-1)*idim1
               incr = 1
            endif
            do i = 1, noseg
                waqoutputs(j,i) = pmsa(ip)
                ip = ip + incr
            enddo
        endif
      enddo 
      return
   end subroutine copy_data_from_wq_processes_to_fm

   !> defaults for process library (WAQ)
   subroutine default_fm_wq_processes()
      use m_fm_wq_processes
      use unstruc_model
      implicit none
      
      jawaqproc = 0
      jamba = 0      
      md_subfile = ''
      md_ehofile = ''
      md_thetav_waq = 0d0
      md_dt_waqproc = 0d0
      
      return
   end subroutine default_fm_wq_processes
