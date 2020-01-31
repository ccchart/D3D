      real function t_forcing(itime,iseg)

!     read and evaluates temperatures forcing from file

      implicit none

!     Arguments

      integer                      :: itime
      integer                      :: iseg

      save

!     Locals

      integer                      :: ifirst = 1
      integer                      :: notime
      integer                      :: noseg
      integer                      :: it, it1, it2
      integer, allocatable         :: t_time(:)
      real   , allocatable         :: temperature(:)
      character(len=256)           :: t_file
      integer                      :: io_err

      if ( ifirst .eq. 1 ) then
         ifirst = 0
         open(849,file='t_forcing.dat')
         read(849,*) notime
         if ( notime .eq. -2 ) then
            ! binary temperature file
            read(849,*) noseg
            read(849,*) t_file
            allocate(temperature(noseg))
            open(848,file=t_file,form='binary')
            read(848) it1,temperature
            read(848) it2
         else
            allocate(t_time(notime))
            allocate(temperature(notime))
            do it = 1, notime
               read(849,*) t_time(it),temperature(it)
            enddo
            it = 1
         endif
      endif

      if ( notime .eq. -2 ) then
         do
            if ( itime .ge. it2 ) then
               it1 = it2
               read(848) temperature
               read(848,iostat=io_err) it2
            else
               exit
            endif
         enddo
         t_forcing = temperature(iseg)
      else
         if ( it .ne. notime ) then
            do
               if ( itime .lt. t_time(it+1) ) exit
               it = it + 1
            enddo
         endif
         t_forcing = temperature(it)
      endif

      return
      end
