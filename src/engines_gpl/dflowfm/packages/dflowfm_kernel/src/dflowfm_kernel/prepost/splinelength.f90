!> approximate spline length
double precision function splinelength(num, xspl, yspl)

   use geometry_module, only: dbdistance
   use m_missing, only: dmiss
   use m_sferic, only: jsferic, jasfer3D

   implicit none

   integer,                          intent(in) :: num          !< number of spline control points
   double precision, dimension(num), intent(in) :: xspl, yspl   !< coordinates of slpine control points

   double precision, dimension(num)             :: xspl2, yspl2 !  second order derivates of spline coordinates

   double precision                             :: xL, yL, xR, yR, tL, tR, dt

   integer                                      :: i,j

   integer, parameter                           :: NSAM = 100   ! sample factor

   call splinxy(xspl, yspl, xspl2, yspl2, num)

   tR = 0d0
   call splintxy(xspl,yspl,xspl2,yspl2,num,tR,xR,yR)

   dt = 1d0/dble(NSAM)

   splinelength = 0d0
   do i=1,num-1
      tR = dble(i-1)
      do j=1,NSAM
         tL = tR
         xL = xR
         yL = yR
         tR = tR + dt
         call splintxy(xspl,yspl,xspl2,yspl2,num,tR,xR,yR)
         splinelength = splinelength + dbdistance(xL,yL,xR,yR, jsferic, jasfer3D, dmiss)
      end do
   end do

   return
end function splinelength
