module partini_mod
!
!  module declarations
!
!  data definition module(s)
!
use precision_part ! single/double precision
use timers
!
!  module procedure(s)
!
! use partmem
use get_key_mod
use grid_search_mod
use pinpok_mod
use wait_mod
!
implicit none ! force explicit typing
!
contains
      subroutine partini( lgrid, lgrid2, nmax, mmax, xcor, ycor, nopart, nosubs, ini_file, wpart  , xpart , &
                          ypart , zpart , npart   , mpart  , kpart , &
                          iptime, lunpr )
!
!     programmer : jan van beek
!     function   : set up of initial condition for a specific number of particles
!     date       : sep 2011
!
!
!     method     : for each partical the position and weigth of the first subtance is read and assigned
!
      integer(ip)                       :: nopart,nosubs
      integer(ip)                       :: lunpr
      character(len=*)                  :: ini_file

      integer(ip), pointer, dimension(:)          :: iptime
      integer(ip), pointer, dimension(:)          :: npart, mpart, kpart
      real   (sp), pointer, dimension(:)          :: xpart, ypart, zpart
      real   (sp), pointer, dimension(:,:)        :: wpart
      integer  ( ip), intent(in   ) :: nmax                    !< first dimension matrix
      integer  ( ip), intent(in   ) :: mmax                    !< second dimension matrix
      integer  ( ip), intent(in   ) :: lgrid (nmax,mmax)       !< active grid matrix
      integer  ( ip), intent(in   ) :: lgrid2(nmax,mmax)       !< total grid matrix
      real     ( rp), intent(in   ) :: xcor  (nmax*mmax)
      real     ( rp), intent(in   ) :: ycor  (nmax*mmax)


      integer(ip)                       :: lun_ini=50
      integer(ip)                       :: ios
      integer(ip)                       :: i
      integer(ip)                       :: np
      integer(ip)                       :: nopart_ini
      integer(ip)                       :: nosubs_ini
      integer(ip)                       :: nosubs_ext
      real                              :: rdummy
!
!     local scalars
!
      real   (sp)                       :: xxcel, yycel

      integer(ip) :: len_file
      integer(ip)                       :: nnpart, mmpart
      integer(ip)                       :: ier, nerr
      !
!     required, otherwise under linux the built-in
!     random generator will be used, rather than the
!     part generator.
!
      integer(4) ithndl              ! handle to time this subroutine
      data       ithndl / 0 /
      if ( timon ) call timstrt( "partini", ithndl )

      len_file          =  len_trim(ini_file)

!      open(lun_ini,file=ini_file,form='binary',status='old',iostat=ios)
      open(lun_ini,file=ini_file,status='old',iostat=ios)
      if (ios /= 0) go to 900
      read(lun_ini,*,end=920,err=930) nopart_ini,nosubs_ini
      if ( nosubs_ini .gt. nosubs ) then
         write ( lunpr,*) ' warning, extra nosubs in initial condition not used'
         nosubs_ext = nosubs_ini - nosubs
         nosubs_ini = nosubs
      else
         nosubs_ext = 0
      endif
      write (*,*)' Reading initial particle distribution from file'
      do np = 1, nopart_ini
         do i = 1, nosubs
           wpart(i,np) = 0.0
         enddo
         
         read (lun_ini,*,end=920,err=930) kpart(np),xpart(np),ypart(np), &
                                        (wpart(i,np),i=1,nosubs_ini),(rdummy,i=1,nosubs_ext)
         iptime(np)  = 0
! now conversion to mnk and local coordinates
! transform world coordinates (xx,yy) into cell-coordinates
! (xxcel,yycel)
!
         call part07 (lgrid  , lgrid2 , nmax   , mmax   , xcor  ,       &
                      ycor   , xpart(np)  , ypart(np)  , nnpart , mmpart,       &
                      xxcel  , yycel  , ier )
            if ( ier == 0 ) then
              nerr = 0
!
!            nopart = nopart + 1 !count total number spread particles

              xpart(np)     = xxcel    !  0 <= xxcel <= 1
              ypart(np)     = yycel    !  0 <= yycel <= 1
            
            
              npart(np)     = nnpart
              mpart(np)     = mmpart
!
!          locate particles at water surface
!          (top of surface layer (z=0.0 ; k=1))
!
              zpart(np)     = 0.5   ! top of layer

!            wpart(isub,nopart) = avgmass*conc2(np+1)  !mass of the particle is scaled with the factor in the dataset (to allow for concentration variations)
              iptime(np)     = 0     ! age of particle
            else
                npart(np)=1
                mpart(np)=1
              nerr = nerr + 1
              if (nerr .gt. 10000) go to 930
            end if

      enddo
      nopart = nopart_ini

      if ( timon ) call timstop ( ithndl )
      return
!     error handling

  900 write(*,'(//a,a)')       ' Error: problem with ini-file ',ini_file(:len_file)
      write(*,'(a)')           ' Could not open/find ini-file ??'
      call wait
      write(lunpr,'(//a,a)')   ' Error: problem with ini-file ',ini_file(:len_file)
      write(lunpr,'(a,a)')     ' Could not open/find ini-file ??'
      stop  ' Part aborted'

  920 write(*,'(//a,a)')       ' Error: problem with ini-file ',ini_file(:len_file)
      write(*,'(//a,a)')       ' End-of-file found on ini-file '
      call wait
      write(lunpr,'(//a,a)')   ' Error: problem with ini-file ',ini_file(:len_file)
      write(lunpr,'(//a,a)')   ' End-of-file found on ini-file '
      stop  ' Part aborted'

  930 write(*,'(//a,a)')       ' Error: problem with ini-file ',ini_file(:len_file)
      write(*,'(//a,a)')       ' Error while reading ini-file'
      call wait
      write(lunpr,'(//a,a)')   ' Error: problem with ini-file ',ini_file(:len_file)
      write(lunpr,'(//a,a)')   ' Error while reading ini-file'
      stop  ' Part aborted'

      end subroutine partini

end module
