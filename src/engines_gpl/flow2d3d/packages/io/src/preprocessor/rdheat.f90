subroutine rdheat(lunmd     ,lundia      ,error     ,nrrec     ,mdfrec      , &
                & noui      ,runid       ,filtem    ,fmttem    ,ktemp       , &
                & rttem     ,dt          ,itstrt    ,itstop    ,mxtemt      , &
                & ntemtm    ,ivapop      ,rhum      ,tdryb     ,qsolar      , &
                & tback     ,tair        ,vapres    ,cfclou    ,solrad_read , &
                & gdp       )
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
!    Function: - Reads the time dependent boundary data for the
!                heat model directly from the MD-file or indi-
!                rectly from the attribute file Filtmp.
!              - Tests the file or data consistency.
!              - Checks whether the file exists or the required
!                data is not empty.
!              - An essential assumption is that the data has to
!                be specified sequentially in time. This imply
!                that NT times Background temperature records
!                should exist in the file (NT unrestricted).
!              - Writes the data to an unformatted file
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use properties
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer , pointer :: itis
!
! Global variables
!
    integer                                      :: itstop      !  Description and declaration in inttim.igs
    integer                                      :: itstrt      !  Description and declaration in inttim.igs
    integer                        , intent(in)  :: ivapop      !  Description and declaration in heat.igs
    integer                        , intent(in)  :: ktemp       !  Description and declaration in tricom.igs
    integer                                      :: lundia      !  Description and declaration in inout.igs
    integer                                      :: lunmd       !  Description and declaration in inout.igs
    integer                        , intent(in)  :: mxtemt      !!  Maximum number of times for which heat model data is allowed in the MD-file
    integer                                      :: nrrec       !!  Pointer to the record number in the MD-file
    integer                        , intent(out) :: ntemtm      !!  Actual number of times for which heat model data is allowed in the MD-file
    logical                                      :: error       !!  Flag=TRUE if an error is encountered
    logical                        , intent(in)  :: noui        !!  Flag for reading from User Interface
    logical                        , intent(in)  :: solrad_read !  Description and declaration in heat.igs
    real(fp)                                     :: dt          !  Description and declaration in rjdim.f90
    real(fp)    , dimension(mxtemt), intent(out) :: cfclou      !  Description and declaration in heat.igs
    real(fp)    , dimension(mxtemt), intent(out) :: qsolar      !!  Sun Insolation (Cloud-less if KTEMP=1 or 4,
                                                                !!  QSUN in HEATU, or net if KTEMP=2, QRADIN in HEATU)
    real(fp)    , dimension(mxtemt), intent(out) :: rhum        !  Description and declaration in heat.igs
    real(fp)    , dimension(mxtemt)              :: rttem       !!  At most MXTEMT times for time varying heat model data
    real(fp)    , dimension(mxtemt), intent(out) :: tair        !  Description and declaration in heat.igs
    real(fp)    , dimension(mxtemt), intent(out) :: tback       !  Description and declaration in heat.igs
    real(fp)    , dimension(mxtemt), intent(out) :: tdryb       !  Description and declaration in heat.igs
    real(fp)    , dimension(mxtemt), intent(out) :: vapres      !  Description and declaration in heat.igs
    character(*)                                 :: filtem      !!  File name for the time varying heat model file
    character(*)                                 :: mdfrec      !!  Standard rec. length in MD-file (300)
    character(*)                                 :: runid       !!  Run identification code for the current simulation (used to determine
                                                                !!  the names of the in- /output files used by the system)
    character(2)                                 :: fmttem      !!  File format for the time varying heat model file
!
! Local variables
!
    integer                      :: iocond   ! IO status for reading 
    integer                      :: item     ! Help var. for times read 
    integer                      :: itold    ! Help var. to store last read time to test accending order 
    integer                      :: ittdep   ! Help var. for the time read (now defined as multiples of DT, but in future it may take any value) 
    integer                      :: l        ! Help var. 
    integer                      :: lenc     ! Help var. (length of var. cvar to be looked for in the MD-file) 
    integer                      :: lf
    integer                      :: lkw      ! Length (in characters) of keyword 
    integer                      :: lrid     ! Length of character string runid 
    integer                      :: lunout   ! Unit number for unformatted FLOW help file between TDATOM and TRISIM  
    integer                      :: n        ! Help var. 
    integer                      :: newlun
    integer                      :: nlook    ! Help var.: nr. of data to look for in the MD-file 
    integer                      :: nrval    ! Number of values to read from file 
    integer                      :: ntrec    ! Help. var to keep track of NRREC 
    logical                      :: dtn
    logical                      :: ex       ! Flag to test if file exists 
    logical                      :: rec1st   ! Flag set to TRUE if the record read is the first record 
    real(fp)                     :: rdummy
    real(fp)      , dimension(5) :: rval     ! Help array (real) where the data, recently read from the MD-file, are stored temporarily 
    character(11)                :: fmtdef   ! Default file format (usually = blank) 
    character(11)                :: fmttmp   ! Help variable for file format 
    character(12)                :: fildef   ! Default file name (usually = blank) 
    character(256)               :: filout   ! Help variable for file name 
    character(6)                 :: keyw     ! Name of record to look for in the MD-file (usually KEYWRD or RECNAM)
    character(300)               :: message
!
!! executable statements -------------------------------------------------------
!
    itis  => gdp%gdrdpara%itis
    !
    ! initialize global parameters
    !
    item   = 1
    itold  = -1
    rec1st = .true.
    nlook  = 1
    fildef = ' '
    filout = ' '
    fmtdef = 'FRformatted'
    fmttmp = ' '
    !
    lunout = 8
    !
    ! Locate 'Filtem' record for time varying heat model data in separate input file containing time-varying data
    ! Keyword does not have to exist, temperature might be specified space-varying via the meteo-module
    !
    filtem = fildef
    call prop_get_string(gdp%mdfile_ptr, '*', 'Filtmp', filtem)
    !
    ! Time varying heat module data in file?
    !
    if (filtem /= fildef) then
       !
       ! It is assumed the file is free-formatted
       !
       fmttem = 'FR'
       fmttmp = 'formatted'
       error  = .false.
       !
       ! If not UI then:
       ! Check filename "filtem" <> TMP file or
       ! "filtem" = TMP file and access is unformatted
       ! Define length of RUNID
       ! open output file (ONLY VERSION 2.48 or lower) +
       ! Set name for the constituents, Regenerated locally
       !
       if (noui) then
          call noextspaces(runid     ,lrid      )
          filout = 'TMP_' // runid(:lrid) // '.tem'
          !
          ! define length of file name
          !
          call noextspaces(filtem, lf)
          !
          ! open unformatted tem-file
          !
          lunout = newlun(gdp)
          inquire (file = filout(:8 + lrid), exist = ex)
          if (ex) then
             open (lunout, file = filout(:8 + lrid))
             close (lunout, status = 'delete')
          endif
          open (lunout, file = filout(:8 + lrid), form = 'unformatted',      &
              & status = 'unknown')
          !
          write (message, '(2a)') 'Reading Heat module file ', filtem(:lf)
          call prterr(lundia, 'G051', trim(message))
          nrval = 1
          if (ktemp /= 3) then
             nrval = 3
             if (ivapop == 1 .or. solrad_read) then
                nrval = 4
             endif
             
          endif
          call rdtdf(lundia    ,lunout    ,error     ,filtem    ,fmttmp    , &
                   & nrval     ,rval      ,dt        ,itstrt    ,itstop    , &
                   & gdp       )
       endif
    else
       !
       ! time varying heat module data in a .tem file
       !
       ! Remove file TMP_<runid>.tem if it exists
       !
       call noextspaces(runid, lrid)
       filout = 'TMP_' // runid(:lrid) // '.tem'
       !
       ! Open file
       !
       lunout = newlun(gdp)
       inquire (file = filout(:8 + lrid), exist = ex)
       if (ex) then
          open (lunout, file = filout(:8 + lrid))
          close (lunout, status = 'delete')
       endif
       !
       ! Check whether heat parameters are handled by the (new) meteo module
       !
       filtem = fildef
       ex     = .false.
                     call prop_get_string(gdp%mdfile_ptr, '*', 'Fwndgr', filtem)
       if (filtem /= fildef) ex = .true.
       if (.not. ex) call prop_get_string(gdp%mdfile_ptr, '*', 'Fwndgt', filtem)
       if (filtem /= fildef) ex = .true.
       if (.not. ex) call prop_get_string(gdp%mdfile_ptr, '*', 'Fwndgc', filtem)
       if (filtem /= fildef) ex = .true.
       if (.not. ex) call prop_get_string(gdp%mdfile_ptr, '*', 'Filwr' , filtem)
       if (filtem /= fildef) ex = .true.
       if (.not. ex) call prop_get_string(gdp%mdfile_ptr, '*', 'Filwt' , filtem)
       if (filtem /= fildef) ex = .true.
       if (.not. ex) call prop_get_string(gdp%mdfile_ptr, '*', 'Filwc' , filtem)
       if (ex) then
          !
          ! Reading/checking heat data is done by the meteo module,
          ! when starting the calculation
          ! Skip this subroutine without error
          !
           error  = .false.
           lunout = 8
           goto 9999
       endif
       !
       ! -->
       !
  110  continue
       !
       ! If an NaN is read -> error
       !
       do n = 1, nrval
           if (isnan(rval(n))) then
              write(message,'(2a)') 'NaN in ', trim(filtem)
              call prterr(lundia    ,'P004'    ,message      )
              !
              error = .true.
              goto 9999
           endif
       enddo
       !
       ! NOTE : in the future one should be able to interpolate across dt
       !
       rttem(item) = rval(1)
       !
       ittdep = nint(rttem(item)/dt)
       if (dtn(ittdep, rttem(item), dt)) then
          if (noui) error = .true.
          call prterr(lundia    ,'U044'    ,'Tstmp'   )
       !
       endif
       !
       ! test times and define minimum time
       !
       if (rec1st) then
          if (ittdep > itstrt) then
             call prterr(lundia    ,'U041'    ,'First time Tstmp  >' )
             if (noui) error = .true.
          endif
          rec1st = .false.
       endif
       !
       if (ittdep <= itold) then
          call prterr(lundia    ,'U060'    ,'Tstmp'   )
          if (noui) error = .true.
       endif
       !
       ! define relative humidity, dry bulb temperature and solar
       ! insolation (cloudiness) data in array for KTEMP=1
       !
       if (ktemp == 1) then
          rhum(item)   = rval(2)
          tdryb(item)  = rval(3)
          qsolar(item) = rval(4)
       endif
       !
       ! define relative humidity, dry bulb temperature and solar
       ! insolation (net) data in array for KTEMP=2
       !
       if (ktemp == 2) then
          rhum(item)   = rval(2)
          tdryb(item)  = rval(3)
          qsolar(item) = rval(4)
       endif
       !
       ! define background temperature data in array for KTEMP=3
       !
       if (ktemp == 3) then
          tback(item) = rval(2)
       endif
       !
       ! define relative humidity, air temperature, solar insolation
       ! (net) and (if IVAPOP = 1) vapour pressure data for KTEMP=4
       !
       if (ktemp == 4) then
          rhum(item)   = rval(2)
          tair(item)   = rval(3)
          qsolar(item) = rval(4)
          if (ivapop == 1) then
             vapres(item) = rval(5)
          endif
       endif
       !
       ! define relative humidity, air temperature and fraction of
       ! sky coverd with clouds data in array for KTEMP=5
       !
       if (ktemp == 5) then
          rhum(item) = rval(2)
          tair(item) = rval(3)
          cfclou(item) = rval(4)
          if (solrad_read) then
             qsolar(item) = rval(5)
          endif
       endif
       !
       ! writing to LUNOUT only if NOUI = .true.
       !
       if (noui) then
          write (lunout) (rval(l), l = 1, nlook)
       endif
       !
       ! check if item exceeds maximum value, then item will be
       ! reset, for NOUI = .true. this will never appear
       !
       if (.not.noui) item = item + 1
       !
       if (item > mxtemt) then
          call prterr(lundia    ,'U157'    ,' '       )
          !
          item = mxtemt
          goto 500
       endif
       !
       ! next time to read
       !
       itold = ittdep
       !
       ! <--
       !
       goto 110
       !
       ! stop reading
       !
       !
       ! define actual number of times for time varying heat data
       ! and define maximum time
       !
  500  continue
       if (itold/= - 1) then
          if (itold<itstop) then
             write(message,'(3a)') 'Last time in file ', trim(filtem), ' <' 
             call prterr(lundia    ,'U042'    ,message)
             error = .true.
             goto 9999
          endif
       endif
       ntemtm = item - 1
    endif
    !
    ! close files
    !
 9999 continue
    if (filtem == fildef) then
       error = .true.
       write (message, '(a)') 'No heat module file(s) found'
       call prterr(lundia, 'P004', message)
    endif
    if (noui .and. lunout /= 8) then
       if (error) then
          close (lunout, status = 'delete')
       else
          close (lunout)
       endif
    endif
end subroutine rdheat
