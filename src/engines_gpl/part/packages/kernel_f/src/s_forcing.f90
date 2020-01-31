      real function s_forcing(itime,iseg,s_prev)

!     read and evaluates salinity forcing from file

      implicit none

!     Arguments

      integer                      :: itime
      integer                      :: iseg
      real                         :: s_prev

      save

!     Locals

      integer                      :: ifirst = 1
      integer                      :: notime
      integer                      :: noseg
      integer                      :: it, it0, it1, it2
      integer, allocatable         :: t_time(:)
      real   , allocatable         :: salinity(:)
      real   , allocatable         :: salprev(:)
      character(len=256)           :: t_file
      integer                      :: io_err

      if ( ifirst .eq. 1 ) then
         ifirst = 0
         open(847,file='s_forcing.dat')
         read(847,*) notime
         if ( notime .eq. -2 ) then
            ! binary salinity file
            read(847,*) noseg
            read(847,*) t_file
            allocate(salinity(noseg))
            allocate(salprev(noseg))
            open(846,file=t_file,form='binary')
            read(846) it1,salinity
            read(846) it2
         else
            allocate(t_time(notime))
            allocate(salinity(notime))
            do it = 1, notime
               read(847,*) t_time(it),salinity(it)
            enddo
            salprev = salinity
            it = 1
         endif
      endif

      if ( notime .eq. -2 ) then
         do
            if ( itime .ge. it2 ) then
               it1 = it2
               salprev = salinity
               read(846) salinity
               read(846,iostat=io_err) it2
            else
               exit
            endif
         enddo
         s_prev    = salprev(iseg)
         s_forcing = salinity(iseg)
      else
         if ( it .ne. notime ) then
            do
               if ( itime .lt. t_time(it+1) ) exit
               it = it + 1
            enddo
         endif
         it0 = max(1,it-1)
         s_prev    = salinity(it0)
         s_forcing = salinity(it)
      endif

      return
      end
