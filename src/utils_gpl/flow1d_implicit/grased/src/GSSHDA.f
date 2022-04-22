      subroutine gsshda (ngrid  ,nfrac  ,nlayer ,ptrla  ,pexla ,deff )
c***********************************************************************
c CVS log information:
c
c $Id$
c
c History:
c $Log: gsshda.F,v $
c Revision 1.2  1995/09/27  10:12:54  kuipe_j
c Maintenance
c
c
c***********************************************************************
c
c     Graded Sediment SHift DAta to old time step
c
c     Declaration of parameters
c
      integer    ngrid ,nfrac      ,nlayer
      real       ptrla (ngrid,nfrac,2)   ,pexla (ngrid,nfrac,2),
     &           deff  (ngrid,2)
c
c     Declaration of local parameters
c
      integer    igr   ,jf
      do 20 jf=1,nfrac
         do 10 igr = 1,ngrid
             ptrla(igr,jf,1) = ptrla(igr,jf,2)
   10    continue
   20 continue

      if (nlayer .eq. 2) then
         do 40 jf=1,nfrac
            do 30 igr = 1,ngrid
                pexla(igr,jf,1) = pexla(igr,jf,2)
   30       continue
   40    continue
      endif

      do 50 igr = 1,ngrid
         deff (igr,1) = deff (igr,2)
   50 continue

      end
