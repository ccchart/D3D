 subroutine slnabs(n,sx1,sy1)
 implicit none
 integer          :: n
 double precision :: sx1,sx2,sy1,sy2
 call shipcoor(n,sx1,sy1,sx2,sy2)
 call lnabs(sx2,sy2)
 end subroutine slnabs
