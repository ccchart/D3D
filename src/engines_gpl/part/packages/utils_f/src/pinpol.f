      subroutine pinpol (x0    , y0    , n     , x     , y     , inside)
c
c----------------------------------------------------------------------
c     date                          : 10-11-1987
c----------------------------------------------------------------------
c
c     userid                        :
c     programmer                    : j.mooiman
c     owner                         :
c     project                       :
c
c     language                      : f77/standard
c     system                        :
c     revisions                     :
c     backup                        :
c     resident library              :
c
c-----------------------------------------------------------------------
c     description :
c     this routine determines whether a point is inside a concave/
c     convex polygon or not.
c     a vertical line is drawn through the point in question;
c     if it crosses the polygon an odd number of times then the
c     point is inside the polygon.
c     (only the upper-half of the vertical-line is tested).
c     reference(s): k van smeden, delft hydraulics
c------------------------------------------------------------------
c
c     variables:
c     det            real          1  orientation of triangle
c     x0,y0   in/    real          1  coords of point in question
c     x ,y    in/    real          n  coords of vertices
c     xt,yt          real       nmax  transformed coords of vertices
c     inside    /out logical          true ;point on edge or at vertex
c                                     true ;point inside polygon
c                                     false;point outside polygon
c     iquadr         integer    nmax  quatre administration
c     n       in/    integer       1  number of vertices
c     ncross         integer       1  counter of intersection
c     nqua           integer       1  added quatres
c
c------------------------------------------------------------------
c
c     quatre administration
c
c     iquadr = 1;  if  xt >= 0  and  yt >= 0
c     iquadr = 2;  if  xt >  0  and  yt <  0
c     iquadr = 4;  if  xt <  0  and  yt >  0
c     iquadr = 5;  if  xt <= 0  and  yt <= 0
c
c------------------------------------------------------------------
c
c     i/o other than by argument list :
c     used logical unit numbers       :
c     referenced routines             :
c     subprogram is referenced by     :
c
c------------------------------------------------------------------
c
      parameter (nmax = 10)
c
      real      x0    , y0    ,
     *          x     ( n     ), y     ( n     ),
     *          xt    ( nmax  ), yt    ( nmax  )
c
      integer   i     , j     ,
     *          nqua  , ncross,
     *          iquadr( nmax  )
c
      logical   inside
c
      if( n .gt. nmax ) stop ' nmax too small, pinpol'
c
c initialisation
c
      inside = .false.
      ncross = 0
c
c translate all vertex coordinates,
c associate vertex with quatre
c
      do 100 i = 1, n
        xt (i) = x (i) - x0
        yt (i) = y (i) - y0
        if( xt (i) .gt. 0.0 )then
          if( yt (i) .lt. 0.0 )then
            iquadr (i) = 2
          else
            iquadr (i) = 1
          endif
        else if( xt (i) .lt. 0.0 )then
          if( yt (i) .gt. 0.0 )then
            iquadr (i) = 4
          else
            iquadr (i) = 5
          endif
        else
          if( yt(i) .gt. 0.0 )then
            iquadr (i) = 1
          else if( yt (i) .lt. 0.0 )then
            iquadr (i) = 5
          else
            inside = .true.
            goto 999
          endif
        endif
  100 continue
c
c compute intersections with positive y-ax
c
      do 200 i = 1, n
        j    = 1 + mod (i,n)
        nqua = iquadr (i) + iquadr (j)
        if( nqua .eq. 5 )then
          ncross = ncross + 1
        else if( nqua .eq. 6 )then
          det = xt (i) * yt (j) - xt (j) * yt (i)
          if( det .gt. 0.0 )then
            if( iquadr (i) .lt. iquadr (j) ) ncross = ncross + 1
          else if (det .lt. 0.0) then
            if( iquadr (i) .gt. iquadr (j) ) ncross = ncross + 1
          else
            inside = .true.
            goto 999
          endif
        endif
  200 continue
c
      if( mod (ncross,2) .ne. 0 ) inside = .true.
c
  999 return
      end
