module m_Orifice
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

   implicit none

   public ComputeOrifice

   type, public :: t_orifice
      double precision         :: crestlevel
      double precision         :: crestwidth
      double precision         :: contrcoeff
      double precision         :: latcontrcoeff
      integer                  :: allowedflowdir
      double precision         :: openlevel
      logical                  :: uselimitflowpos
      double precision         :: limitflowpos
      logical                  :: uselimitflowneg
      double precision         :: limitflowneg
   end type

contains

subroutine ComputeOrifice(orifice, fum, rum, aum, dadsm, kfum, s1m1, s1m2, qm, q0m,   &
                 & u1m, u0m, dxm, dt)
!!--copyright-------------------------------------------------------------------
! Copyright (c) 2003, Deltares. All rights reserved.
!!--disclaimer------------------------------------------------------------------
! This code is part of the Delft3D software system. Deltares has
! developed c.q. manufactured this code to its best ability and according to the
! state of the art. Nevertheless, there is no express or implied warranty as to
! this software whether tangible or intangible. In particular, there is no
! express or implied warranty as to the fitness for a particular purpose of this
! software, whether tangible or intangible. The intellectual property rights
! related to this software code remain with Deltares at all times.
! For details on the licensing agreement, we refer to the Delft3D software
! license and any modifications to this license, if applicable. These documents
! are available upon request.
!!--version information---------------------------------------------------------
! $Author$
! $Date$
! $Revision$
!!--description-----------------------------------------------------------------
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    !=======================================================================
    !                       Deltares
    !                One-Two Dimensional Modelling System
    !                           S O B E K
    !
    ! Subsystem:          Flow Module
    !
    ! Programmer:         G. Stelling
    !
    ! Module:             orifice (ORIFICE)
    !
    ! Module description: coefficients for momentum equation in wet orifice point
    !
    !
    !     update information
    !     person                    date
    !
    !
    !
    !     Include Pluvius data space
    !
    use m_GlobalParameters
    use m_struc_helper
    implicit none
!
! Global variables
!
    type(t_orifice), pointer        :: orifice
    integer, intent(inout)          :: kfum
    double precision, intent(inout) :: aum
    double precision, intent(out)   :: dadsm
    double precision, intent(out)   :: fum
    double precision, intent(out)   :: q0m
    double precision, intent(out)   :: qm
    double precision, intent(out)   :: rum
    double precision, intent(in)    :: u0m
    double precision, intent(inout) :: u1m
    double precision, intent(in)    :: s1m1
    double precision, intent(in)    :: s1m2
    double precision, intent(in)    :: dxm
    double precision, intent(in)    :: dt
!
!
! Local variables
!
    integer                        :: allowedflowdir
    double precision               :: cu
    double precision               :: dabs
    double precision               :: dsqrt
    double precision               :: fr
    double precision               :: qtm
    double precision               :: rhsc
    double precision               :: scf
    double precision               :: scr
    double precision               :: smax
    double precision               :: smin
    double precision               :: sop
    double precision               :: swi
    double precision               :: uweir
    double precision               :: dxdt
    double precision               :: cmu

!
!
!! executable statements -------------------------------------------------------
!
    !
    !
    !

    scr = orifice%crestlevel
    swi = orifice%crestwidth
    scf = orifice%contrcoeff
    sop = orifice%openlevel
    allowedflowdir = orifice%allowedflowdir
    cmu = orifice%latcontrcoeff
    smax = max(s1m1, s1m2)
    smin = min(s1m1, s1m2)
    !
    if ((allowedflowdir==3) .or.                        &
        (s1m1 - s1m2>=0.0 .and. allowedflowdir==2) .or. &
        (s1m1 - s1m2<=0.0 .and. allowedflowdir==1) .or. &
        (sop - scr<=.0)) then
       kfum = 0
    !           ARS 9698 .lt. is change to .le.
    elseif (smax - scr<thresholdDry) then
       kfum = 0
    elseif ((smax - scr>thresholdFlood)) then
       kfum = 1
    else
    endif
    if (kfum==0) then
       fum = 0.0
       rum = 0.0
       u1m = 0.0
       qm = 0.0
       q0m = 0.0
       return
    endif

    ! ARS 7367, 7966 wrong value of dynstruc set here
    !     dynstruc=1.0  dynstruc is set by ini-file

    rhsc = 0.0

    !     made similar condition as in wtcrsu
    if (2./3.*(smax - scr)<sop - scr) then
       !        weir flow; 2/3 ( h_1 - z_s) < d_g = s_op - z_s
       if (smax - scr<=1.5*(smin - scr)) then
          !          submerged weir flow;  h_1 - z_s <= 3/2 (h_2 - z_s)
          cu = cmu**2*2.0d0*gravity/(StructureDynamicsFactor*dxm)
          !          ARS 4681 improved wetted area computation
          !          ARS 3479 wetted area orifice limited to opening
          aum = min(max(smax - u0m**2/(2.0d0*gravity) - scr, 2.0d0/3.0d0*(smax - scr))*swi,   &
              & scf*(sop - scr)*swi)
          uweir = cmu*dsqrt(gravity*2.0d0*(smax - smin))
          fr = dabs(uweir)/(StructureDynamicsFactor*dxm)
          rhsc = 0.0
          dadsm = swi
       else
          !          free weir flow;  h_1 - z_s > 3/2 (h_2 - z_s)
          aum = (2.0d0/3.0d0)*(smax - scr)*swi
          cu = cmu**2*gravity/(1.5d0*(StructureDynamicsFactor*dxm))
          uweir = cmu*dsqrt(2.0d0/3.0d0*gravity*(smax - scr))
          fr = dabs(uweir)/(StructureDynamicsFactor*dxm)
          if (s1m2>s1m1) then
             rhsc = -cu*(s1m1 - scr)
          else
             rhsc = cu*(s1m2 - scr)
          endif
          dadsm = swi
       endif
    else
       !        orifice flow; 2/3 ( h_1 - z_s) >= d_g = s_op - z_s
       if (smin>sop) then
          !          submerged orifice flow;  h_2 > z_s + d_g = z_s + s_op - z_s = s_op
          cu = cmu**2*2.0d0*gravity/(StructureDynamicsFactor*dxm)
          uweir = cmu*dsqrt(gravity*2.0d0*(smax - smin))
          fr = dabs(uweir)/(StructureDynamicsFactor*dxm)
          rhsc = 0.0
       else
          !          free orifice floww;  h_2 <= z_s + d_g = z_s + s_op - z_s = s_op
          cu = cmu**2*2.0d0*gravity/(StructureDynamicsFactor*dxm)
          uweir = cmu*dsqrt(2.0d0*gravity*(smax - (scr + scf*(sop - scr))))
          fr = uweir/(StructureDynamicsFactor*dxm)
          if (s1m2>s1m1) then
             rhsc = -cu*(smin - (scr + scf*(sop - scr)))
          else
             rhsc = cu*(smin - (scr + scf*(sop - scr)))
          endif
       endif
       aum = scf*(sop - scr)*swi
       dadsm = swi
    endif
    dxdt = 1.0/dt
    
    call furu_iter(fum, rum, s1m2, s1m1, u1m, u0m, q0m, aum, fr, cu, rhsc, dxdt)

    qtm = aum*u1m
    !
    !     check for restriction on flow
    if ((orifice%uselimitflowpos) .and. (qtm>0)) then
       if (dabs(qtm)>orifice%limitflowpos ) then
          fum = 0.0
          rum = orifice%limitflowpos/max(aum, 1.0D-4)
          u1m = rum
          qm = orifice%limitflowpos
       endif
    elseif ((orifice%uselimitflowneg) .and. (qtm<0)) then
       if (dabs(qtm)>orifice%limitflowneg) then
          fum = 0.0
          rum = -1.0*orifice%limitflowneg/max(aum, 1.0D-4)
          u1m = rum
          qm = -1.0*orifice%limitflowneg
       endif
    else
    endif
end subroutine ComputeOrifice


end module m_Orifice
