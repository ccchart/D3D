!> process command line option
module m_commandline_option
implicit none
   integer, parameter  :: MAXOPTLEN = 255          !< maximum option string length
   integer, parameter  :: MAXKEYLEN = 255          !< maximum key string length
   integer, parameter  :: MAXKEYS   = 16           !< maximum number of key-value pairs
   integer, parameter  :: MAXSTRLEN = 255          !< maximum value string length
   integer, parameter  :: IMISS = -999999          !< unset value

!  input files from command line
   integer, parameter                             :: lenfile  = 255         !< maximum filename length
   integer, parameter                             :: maxnumfiles = 10       !< maximum number of files
   integer                                        :: numfiles               !< number of files to be loaded
   character(len=lenfile), dimension(maxnumfiles) :: inputfiles             !< files to be loaded
   character(len=lenfile)                         :: iarg_outfile = ' '     !< Output filename for several commandline/batch-mode operations (not related to model runs).
   integer                                        :: iarg_autostart         !< autostart/autstartstop or not set (-1)
   integer                                        :: iarg_usecaching = -1   !< use cache file or not or not set (-1)

contains
!> read, from a string, command line option with key-value pair(s) of the form
!>  <->[-[...]]<option>[:key1=val1[:key2=val2[...]]]
!> key=val are string-integer pairs
   subroutine read_commandline_option(str, Soption, Nkeys, Skeys, ivals, svals)
      implicit none

      character(len=*),                             intent(in)  :: str        !< string
      character(len=MAXOPTLEN),                     intent(out) :: Soption    !< option
      integer,                                      intent(out) :: Nkeys      !< number of keys
      character(len=MAXKEYLEN), dimension(MAXKEYS), intent(out) :: Skeys      !< keys
      integer,                  dimension(MAXKEYS), intent(out) :: ivals      !< values
      character(len=MAXSTRLEN), dimension(MAXKEYS), intent(out) :: svals

      integer                                                   :: ipoint, ibegin, iend, iequal, Lenstr

      Soption = ''
      Nkeys = 0
      Skeys = ''
      ivals = IMISS

      Lenstr = len_trim(str)

   !  ignore leading dashes -/--/... (but dashes inside option names/values are allowed)
      ipoint=1
      do while (ipoint <= Lenstr)
         if (str(ipoint:ipoint) /= '-') exit
         ipoint = ipoint+1
      end do

   !  proceed if argument is option only
      if ( ipoint.eq.1 ) return

   !  read option
      ibegin = ipoint
      iend   = index(str(ipoint:Lenstr), ':')+ibegin-2
      if ( iend.lt.ibegin ) iend = Lenstr
      Soption=trim(adjustL(str(ibegin:iend)))
      ipoint = iend+1

   !  read key-value pairs
      do while ( ipoint.lt.Lenstr )
   !     find begin of this substring, that succeeds the last colon
         ibegin = ipoint+1

   !     find end of this substring, that preceeds the next colon or end of string
         iend = index(str(ibegin:Lenstr), ':')+ibegin-2
         if ( iend.le.ibegin ) iend = Lenstr

   !     find equal sign
         iequal = index(str(ibegin:iend), '=')+ibegin-1
         if ( iequal.lt.ibegin ) then
   !        no equal sign: no value
            iequal = iend+1
         end if
         if ( iequal.gt.ibegin ) then
            Nkeys=Nkeys+1
   !        read key
            Skeys(Nkeys) = trim(adjustL(str(ibegin:iequal-1)))
            if ( iequal.lt.iend ) then
   !           read value as string
               read(str((iequal+1):iend), *, err=666) svals(Nkeys)
   !           read value as integer
               read(str((iequal+1):iend), *, err=666) ivals(Nkeys)
   666         continue
            end if
         end if

         ipoint = iend+1
      end do

      return
   end subroutine read_commandline_option

end module
