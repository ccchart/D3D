!!  Copyright (C)  Stichting Deltares, 2012-2015.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

module pmov_boom_mod
!
!  module declarations
!
!  data definition module(s)
!
use precision_part      ! single/double precision
      use timers
!
implicit none      ! force explicit typing
!
contains
      subroutine pmov_boom(xl, yl, n, x, y, xn,yn)

!     programmer : frank kleissen
!     function   : moves particle to new position of moving boom
!     date       : aug 2015
!
      real(sp), dimension(:) :: x, y
!
!     local scalars
!
      integer     :: i    , i1    , i2  , inside ,  n , np
      real        :: amiss, rechts, rl  , rm
      real(rp)     :: x1   , x2    , xl  , y1  , y2  , yl, xn, yn, xdif, ydif
      
      integer(4) ithndl              ! handle to time this subroutine
      
      ! locals
      real    :: distp1,distp2,distp4,distlo,distln,distpx, distpy, dist14,dist23 
      real    :: distmin  = 5.0
      
      data       ithndl / 0 /
      if ( timon ) call timstrt( "pmov_boom", ithndl )
      !first determine the projected distance on the original line, for both directions
      ! element 1,2 of the x, y vector
      distp1 = sqrt((xl-x(1))*(xl-x(1))+(yl-y(1))*(yl-y(1))) ! x direction
      distp2 = sqrt((xl-x(2))*(xl-x(2))+(yl-y(2))*(yl-y(2)))
      distp4 = sqrt((xl-x(4))*(xl-x(4))+(yl-y(4))*(yl-y(4)))
      distlo = sqrt((x(1)-x(2))*(x(1)-x(2))+(y(1)-y(2))*(y(1)-y(2)))
      distln = sqrt((x(3)-x(4))*(x(3)-x(4))+(y(3)-y(4))*(y(3)-y(4)))
      dist14 = sqrt((x(1)-x(4))*(x(1)-x(4))+(y(1)-y(4))*(y(1)-y(4)))
      dist23 = sqrt((x(3)-x(2))*(x(3)-x(2))+(y(3)-y(2))*(y(3)-y(2)))
      
      if (distlo.gt.0) then
        distpx = (distp1*distp1+distlo*distlo-distp2*distp2)/2.0/distlo
      ! this distance can be negative or greater than the original line length
      ! but then move it to the corner coordinates of the nwe line and move it a little inwards
      ! new coordinates on the new line
      ! note that element 1 of the old line moves to element 4 of the new line
      ! the 1.001 is to ensure that the new point is just outside the polygon
      if (distpx.lt.0) then
        xn = x(4)
        yn = y(4)
      elseif (distpx.gt.distlo) then
        xn = x(3)
        yn = y(3)
      else
        xn = x(4) + (x(3)-x(4))*distpx/distlo
        yn = y(4) + (y(3)-y(4))*distpx/distlo
      endif
      
        !difference x and y of new and old point
        xdif=xn-xl
        ydif=yn-yl
        ! extra movement
        yn = yn + (y(4)-y(1))*0.5
        xn = xn + (x(4)-x(1))*0.5
        ! nieuwe positie moet buiten de polygoon vallen, dus een aantal meters verder leggen.
      else
        write(*,*)' Distance is zero, no movement of particle'
      endif
        
!
99999 if ( timon ) call timstop ( ithndl )
      return
      end subroutine pmov_boom
end module
