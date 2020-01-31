module wait_mod
      use timers
implicit none
contains
      subroutine wait
!
!     purpose: replacement for F90 PAUSE statement
!
      integer(4) ithndl              ! handle to time this subroutine
      data       ithndl / 0 /
      if ( timon ) call timstrt( "wait", ithndl )
      write(*,'(a)') ' Press <Enter> to continue ...'
      read (*,'(a)',advance='no')
!
      if ( timon ) call timstop ( ithndl )
      return
      end subroutine
end module
