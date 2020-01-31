      SUBROUTINE FINDNM (XLOC  , YLOC  , X     , Y     , NMAX , MMAX  ,
     *                   NC    , MC    , INSIDE, NLOC  , MLOC )
c***********************************************************************
c
c delft hydraulics laboratory        department estuaries and se
c
c depint       .001  august 10 1990          herman kernkamp
c
c last update        june   06 1997          johan boon
c                                            inverse nmx, mmax,uitkleding
c last update        june   06 1992          theo van der kaaij
c                    depth points instead of zeta points
c                    x,y 2-dimensional arrays
c                    no "smart-search"
c
c function           : given a pair of x,y coordinates, find the
c                      m- and n-coordinates
c entry              : subroutine call
c exit               : calling module
c messages           : none
c logical units      : none
c subroutines called : pinpol
c machine dependent  : none
c
c parameters         :
c
c name    type     i o t length    description
c ----    ----     - - - ------    -----------
c
c xp,yp  real      *     1         point to search
c x,y    real      *     mmax,nmax coordinates depth points
c mmax   integer   *     1         array dimension x,y
c nmax   integer   *     1         array dimension x,y
c mc     integer   *     1         actual number of grid rows
c nc     integer   *     1         actual number of grid columns
c inside boolean     *   1         true if found
c m,n    integer     *   1         trisula m,n coordinates
c rm,rn  real        *   1         mantisse of trisula coordinat
c loctyp char*7    *     1         'station', for STATION
c                                  'drogue ', for DROGUE TRACK
c***********************************************************************
c
      integer       mmax  , nmax  , mc    , nc    , m     , n     ,
     *              i     , j
c
c
      REAL          X     ( NMAX  , MMAX  ),
     *              Y     ( NMAX  , MMAX  )
c
      real          xx    ( 5     )       , yy    ( 5     )
c
      logical       inside
c
c-----------------------------------------------------------------------
c
      eps    = 1.0 e-10
      inside = .false.
c
      do 10 i = 2, NC
         do 10 j = 2, MC
c
            xx (1) = x (i-1,j-1)
            xx (2) = x (i  ,j-1)
            xx (3) = x (i  ,j  )
            xx (4) = x (i-1,j  )
            xx (5) = x (i-1,j-1)
            yy (1) = y (i-1,j-1)
            yy (2) = y (i  ,j-1)
            yy (3) = y (i  ,j  )
            yy (4) = y (i-1,j  )
            yy (5) = y (i-1,j-1)
c
            IF( ABS(XX (1)) .GT. EPS .AND. ABS(XX (2)) .GT. EPS .AND.
     *          ABS(XX (3)) .GT. EPS .AND. ABS(XX (4)) .GT. EPS )THEN
c
c------------- search grid cell
c
               CALL PINPOL (XLOC  , YLOC  , 5    , XX   , YY   , INSIDE)
c
               IF( INSIDE )THEN
                   NLOC = I
                   MLOC = J
                   GOTO 20
               ENDIF
            ENDIF
c
   10 continue
   20 continue
c
      return
      end
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
