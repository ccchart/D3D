subroutine tranb7(utot      ,d50       ,d90       ,h         ,par       , &
                & sbot        ,ssus      ,vonkar    ,mudfrac   ,i2d3d     , &
                & chezy       ,sag       ,f_VR84    ,utotVARIES,hVARIES   , &
                & rmuc        ,fc        ,wsCOMP    ,kmax)
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2017.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
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
!!--description-----------------------------------------------------------------
!
! computes sediment transport according to
! van rijn (1984)
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !use glob_bankPLIC ,only: ratio_ca_c2d,simpleVR84,prescVR93refHEIGHT,prescVR93settl
    implicit none
!
! Call variables
!
    real(fp)               , intent(in)  :: d50     ! grain size diameter (first specified diameter)
    real(fp)               , intent(in)  :: d90     ! grain size diameter (first specified diameter)
    real(fp)               , intent(in)  :: h       ! water depth
    real(fp)               , intent(in)  :: mudfrac ! fraction of mud
    real(fp)               , intent(out) :: sbot    ! bed load transport
    real(fp)               , intent(out) :: ssus    ! suspended sediment transport
    real(fp)               , intent(in)  :: utot    ! flow velocity
    real(fp)               , intent(in)  :: utotVARIES   ! space and time varying flow velocity 
    real(fp)               , intent(in)  :: hVARIES      ! space and time varying depth
    real(fp), dimension(30), intent(in)  :: par     ! sediment parameter list
    real(fp)               , intent(in)  :: vonkar
    real(fp)               , intent(in)  :: chezy
    real(fp)               , intent(out) :: f_VR84
    real(fp)               , intent(in)  :: sag
    real(fp)               , intent(in)  :: wsCOMP(0:kmax)
    integer                , intent(in)  :: i2d3d
    integer                , intent(in)  :: kmax
!
! Local variables
!
    real(fp)       :: z0cur
    real(fp)       :: a
    real(fp)       :: ah
    real(fp)       :: alf1
    real(fp)       :: beta   ! lowest level of integration interval over vertical
    real(fp)       :: betam  !  power coefficient for adjusting critical bed shear stress for sand-mud interaction
    real(fp)       :: ca
    real(fp)       :: del
    real(fp)       :: dstar
    real(fp)       :: fc
    real(fp)       :: ff     ! coriolis coefficient
    real(fp)       :: ag     ! gravity acceleration
    real(fp)       :: psi
    real(fp)       :: rhosol ! density of sediment
    real(fp)       :: rhowat ! density of water
    real(fp)       :: rksc
    real(fp)       :: rmuc
    real(fp)       :: rnu    ! laminar viscosity of water
    real(fp)       :: t      ! dimensionless relative shear stress
    real(fp)       :: tbc
    real(fp)       :: tbce
    real(fp)       :: tbcr
    real(fp)       :: thetcr
    real(fp)       :: ustar
    real(fp)       :: ws     ! settling velocity
    real(fp)       :: zc
    real(fp), external :: shld
!
!! executable statements -------------------------------------------------------
!
    sbot = 0.0
    ssus = 0.0
    !
    ag = par(1)
    rhowat = par(2)
    rhosol = par(3)
    del  = par(4)
    rnu = par(5)
    alf1 = par(11)
    rksc = par(13)
    ws = par(14)
    betam = par(15)
    !
    if (h/rksc<1.33 .or. utot<1.E-3) then
       return
    endif
    !
    !if (prescVR93refHEIGHT) THEN ! I overwrite the reference height rksc (read from the tra file) with the van rijn 1993 ref height (expressions taken from subroutine tram1)
    !   z0cur = h/(exp(1._fp)*(exp(vonkar*chezy/sag) - 1.0_fp)) !it has to be recomputed so h is on the axis when called with simpleVR84==1
    !   rksc = max(30._fp*z0cur, 0.01_fp*h)
    !   rksc = min(rksc, 0.2_fp * h)
    !endif
    !if (prescVR93settl) THEN ! I overwrite the reference settling (read from the tra file) with the van rijn   settling velocity computed in fallve
    !   ws = wsCOMP(kmax)
    !endif
    a = rksc
   dstar = d50*(del*ag/rnu/rnu)**(1._fp/3._fp)
    !
    !if (simpleVR84==2) then
    !   tbce = rhowat*ag*utot**2/Chezy**2       
    !   fc = 99999999999 ! not used
    !   thetcr = shld(dstar)
    !   tbcr = (rhosol - rhowat)*ag*d50*thetcr
    !   t = (tbce - tbcr)/tbcr
    !   if (t<.000001) t = .000001
    !   ca = .015_fp*alf1*d50/a*t**1.5_fp/dstar**.3_fp
    !   !if (ratio_ca_c2d>0) then
    !   !   ff = 1._fp/ratio_ca_c2d
    !   !else
    !      ustar = sqrt(tbce/rhowat)
    !      zc = 0.
    !      beta = 1. + 2.*(ws/ustar)**2
    !      beta = min(beta, 1.5_fp)
    !      psi = 2.5*(ws/ustar)**0.8*(ca/0.65)**0.4
    !      if (ustar>0.) zc = ws/vonkar/ustar/beta + psi
    !      if (zc>20.) zc = 20.
    !      ah = a/h
    !      ff = 0.
    !      if (abs(zc - 1.2)>1.E-4) then
    !         ff = (ah**zc - ah**1.2)/(1. - ah)**zc/(1.2 - zc)
    !      else
    !         ff = -(ah/(1. - ah))**1.2*log(ah)
    !      endif
    !   !endif
    !else
       rmuc = (log10(12.*h/rksc)/log10(12.*h/3./d90))**2
       fc = .24*(log10(12.*h/rksc))**( - 2)
       tbc = .125*rhowat*fc*utot**2
       tbce = rmuc*tbc
       thetcr = shld(dstar)
       tbcr = (rhosol - rhowat)*ag*d50*thetcr*(1.0+mudfrac)**betam
       t = (tbce - tbcr)/tbcr
       !
       if (t<.000001) t = .000001
       ca = .015*alf1*d50/a*t**1.5/dstar**.3 !at the axis if simpleVR84==1
       !
       ustar = sqrt(.125*fc)*utot
       zc = 0.
       beta = 1. + 2.*(ws/ustar)**2
       beta = min(beta, 1.5_fp)
       psi = 2.5*(ws/ustar)**0.8*(ca/0.65)**0.4
       if (ustar>0.) zc = ws/vonkar/ustar/beta + psi
       if (zc>20.) zc = 20.
       ah = a/h
       ff = 0.
       if (abs(zc - 1.2)>1.E-4) then
          ff = (ah**zc - ah**1.2)/(1. - ah)**zc/(1.2 - zc)
       else
          ff = -(ah/(1. - ah))**1.2*log(ah)
       endif
       !if (simpleVR84==1) then !compute the local Ceq with the time and spacial varying utotVARIES
       !   tbc = .125*rhowat*fc*utotVARIES**2
       !   tbce = rmuc*tbc
       !   t = (tbce - tbcr)/tbcr
       !   !
       !   if (t<.000001) t = .000001
       !   ca = .015*alf1*d50/a*t**1.5/dstar**.3 !at the axis if simpleVR84==1
       !endif
    !endif
    ssus = ff*utotVARIES*hVARIES*ca !ff*utotVARIES*hVARIES*ca
    f_VR84 = ff ! store it in order to pass it out.
    !
    if (t<3.) then
       sbot = 0.053*(del)**0.5*sqrt(ag)*d50**1.5*dstar**( - 0.3)*t**2.1
    else
       sbot = 0.100*(del)**0.5*sqrt(ag)*d50**1.5*dstar**( - 0.3)*t**1.5
    endif
end subroutine tranb7
