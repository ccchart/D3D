   module m_Dambreak
   !----- AGPL --------------------------------------------------------------------
   !
   !  Copyright (C)  Stichting Deltares, 2017-2018.
   !
   !  This program is free software: you can redistribute it and/or modify
   !  it under the terms of the GNU Affero General Public License as
   !  published by the Free Software Foundation version 3.
   !
   !  This program is distributed in the hope that it will be useful,
   !  but WITHOUT ANY WARRANTY; without even the implied warranty of
   !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   !  GNU Affero General Public License for more details.
   !
   !  You should have received a copy of the GNU Affero General Public License
   !  along with this program.  If not, see <http://www.gnu.org/licenses/>.
   !
   !  contact: delft3d.support@deltares.nl
   !  Stichting Deltares
   !  P.O. Box 177
   !  2600 MH Delft, The Netherlands
   !
   !  All indications and logos of, and references to, "Delft3D" and "Deltares"
   !  are registered trademarks of Stichting Deltares, and remain the property of
   !  Stichting Deltares. All rights reserved.
   !
   !-------------------------------------------------------------------------------
   !  $Id$
   !  $HeadURL$
   !-------------------------------------------------------------------------------

   use m_GlobalParameters
   use m_struc_helper

   implicit none

   public prepareComputeDambreak
   public setCoefficents

   type, public :: t_dambreak
      double precision :: start_location_x
      double precision :: start_location_y
      integer          :: algorithm
      double precision :: crestlevelini
      double precision :: breachwidthini
      double precision :: crestlevelmin
      double precision :: timetobreachtomaximumdepth
      double precision :: dischargecoeff
      double precision :: f1
      double precision :: f2
      double precision :: ucrit
      double precision :: t0
      integer          :: hasTable
      integer          :: materialtype

      ! State variables, not to be read
      integer          :: phase
      double precision :: width
      double precision :: crl
      double precision :: aCoeff
      double precision :: bCoeff
      double precision :: maximumAllowedWidth

   end type

   private

   contains

   subroutine prepareComputeDambreak(dambreak, s1m1, s1m2, u0, time0, time1, dt, maximumWidth)


   type(t_dambreak), pointer, intent(inout) :: dambreak
   double precision, intent(in)             :: s1m1
   double precision, intent(in)             :: s1m2
   double precision, intent(in)             :: u0
   double precision, intent(in)             :: time0
   double precision, intent(in)             :: time1
   double precision, intent(in)             :: dt
   double precision, intent(in)             :: maximumWidth

   !locals
   double precision :: smax
   double precision :: smin
   double precision :: hmx
   double precision :: hmn
   double precision :: deltaLevel
   double precision :: deltaWidth
   double precision :: timeFromBreaching
   double precision :: breachWidth
   double precision :: actualMaximumWidth

   timeFromBreaching = time1 - dambreak%t0

   ! breaching not started
   if (timeFromBreaching < 0d0) return

   !vdKnaap(2000) formula: to do: implement table 
   if(dambreak%algorithm == 1) then     

      ! The linear part
      if (timeFromBreaching < dambreak%timetobreachtomaximumdepth ) then
         dambreak%crl    = dambreak%crestlevelini - timeFromBreaching / dambreak%timetobreachtomaximumdepth * (dambreak%crestlevelini - dambreak%crestlevelmin)
         breachWidth     = dambreak%breachwidthini
      else
      ! The logarithmic part, timeFromBreaching in seconds 
         breachWidth = dambreak%aCoeff * dlog(timeFromBreaching/dambreak%bCoeff)
      endif
      
      ! breach width must increase monotonically 
      if (breachWidth > dambreak%width ) then
         dambreak%width = breachWidth
      endif
      
      ! in vdKnaap(2000) the maximum allowed branch width is limited (see sobek manual and setCoefficents subroutine below) 
      actualMaximumWidth = min(dambreak%maximumAllowedWidth, maximumWidth)
      
   ! Verheij-vdKnaap(2002) formula
   else if (dambreak%algorithm == 2) then  

      ! phase 1: lowering
      if (timeFromBreaching < dambreak%timetobreachtomaximumdepth) then
         dambreak%crl    = dambreak%crestlevelini - timeFromBreaching / dambreak%timetobreachtomaximumdepth * (dambreak%crestlevelini - dambreak%crestlevelmin)
         dambreak%width  = dambreak%breachwidthini
         dambreak%phase  = 1
      ! phase 2: widening
      else
         dambreak%crl = dambreak%crestlevelmin
         smax = max(s1m1, s1m2)
         smin = min(s1m1, s1m2)
         hmx = max(0d0,smax - dambreak%crl)
         hmn = max(0d0,smin - dambreak%crl)
         deltaLevel = gravity*(hmx - hmn)**1.5d0

         if (dambreak%width < maximumWidth .and. (.not.isnan(u0)) .and. dabs(u0) > dambreak%ucrit) then
            dambreak%width = dambreak%width  + (dambreak%f1*dambreak%f2/dlog(10D0)) * &
                             (deltaLevel/(dambreak%ucrit*dambreak%ucrit)) * &
                             (1.0/(1.0 + (dambreak%f2*gravity*timeFromBreaching/(3600.0d0*dambreak%ucrit)))) * &
                             (dt/3600.0d0)
         endif
      endif
      
      ! in Verheij-vdKnaap(2002), there is no limitation for the maximum allowed branch width
      actualMaximumWidth =  maximumWidth
      
   endif

   !width cannot exceed the width of the snapped polyline
   if(dambreak%width >= actualMaximumWidth) then
      dambreak%width = actualMaximumWidth
   endif

   end subroutine prepareComputeDambreak

   
   subroutine setCoefficents(dambreak)

   type(t_dambreak), pointer, intent(inout) :: dambreak

   if (dambreak%algorithm == 1) then
      ! clay
      if (dambreak%materialtype == 1) then 
         dambreak%aCoeff = 20
         dambreak%bCoeff = 288
         dambreak%maximumAllowedWidth = 75  !meters
      ! sand
      else if(dambreak%materialtype == 2) then 
         dambreak%aCoeff = 67
         dambreak%bCoeff = 522
         dambreak%maximumAllowedWidth = 200 !meters
      endif
   endif

   end subroutine setCoefficents

   end