subroutine calseddf1993(ustarc    ,ws        ,tp        ,delr      ,dstar     , &
                      & uwb       ,hrms      ,h1        ,seddif    ,kmax      , &
                      & sig       ,thick     ,dicww     ,tauwav    ,tauc      , &
                      & ltur      ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
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
!!--description-----------------------------------------------------------------
!
! Compute sediment diffusion coefficient
! Van Rijn (1993)
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    real(fp)               , pointer :: eps
    real(fp)               , pointer :: bed
    real(fp)               , pointer :: vonkar
    logical                , pointer :: const
    logical                , pointer :: wave
    logical                , pointer :: sedim
!
! Global variables
!
    integer                    , intent(in)  :: kmax   !  Description and declaration in iidim.f90
    integer                    , intent(in)  :: ltur   !  Description and declaration in iidim.f90
    real(fp)                   , intent(in)  :: delr
    real(fp)                   , intent(in)  :: dstar
    real(fp)                   , intent(in)  :: h1
    real(fp)                   , intent(in)  :: hrms   !  Description and declaration in rjdim.f90
    real(fp)                   , intent(in)  :: tauc
    real(fp)                   , intent(in)  :: tauwav
    real(fp)                   , intent(in)  :: tp     !  Description and declaration in rjdim.f90
    real(fp)                   , intent(in)  :: ustarc
    real(fp)                   , intent(in)  :: uwb
    real(fp), dimension(0:kmax), intent(in)  :: dicww  !  Description and declaration in rjdim.f90
    real(fp), dimension(0:kmax), intent(out) :: seddif !  Description and declaration in rjdim.f90
    real(fp), dimension(0:kmax), intent(in)  :: ws     !  Description and declaration in rjdim.f90
    real(fp), dimension(kmax)  , intent(in)  :: sig    !  Description and declaration in rjdim.f90
    real(fp), dimension(kmax)  , intent(in)  :: thick  !  Description and declaration in rjdim.f90
!
! Local variables
!
    integer                        :: k
    real(fp)                       :: beta
    real(fp)                       :: betaef
    real(fp)                       :: deltas
    real(fp)                       :: epsbed
    real(fp)                       :: epsmax
    real(fp)                       :: epsmxc
    real(fp)                       :: lci
    real(fp), dimension(0:kmax)    :: epscur
    real(fp), dimension(0:kmax)    :: epswav
!
!! executable statements -------------------------------------------------------
!
    eps                 => gdp%gdconst%eps
    bed                 => gdp%gdmorpar%bed
    vonkar              => gdp%gdphysco%vonkar
    const               => gdp%gdprocs%const
    wave                => gdp%gdprocs%wave
    sedim               => gdp%gdprocs%sedim
    !
    ! **** CALCULATE VERTICAL SEDIMENT DIFFUSION COEFFICIENT ****
    !
    if (ustarc>eps) then
       !
       ! Beta factor assumed constant over the depth, using ws at
       ! water surface (approx. clear water). Beta limited to
       ! range 1 - 1.5 following Van Rijn (1993)
       !
       beta = 1.0 + 2.0*(ws(1)/ustarc)**2.0
       beta = max(1.0_fp, beta)
       beta = min(1.5_fp, beta)
    else
       beta = 1.0
    endif
    !
    if ((ltur==0) .or. (ltur==1)) then
       !
       ! if algebraic or K-L turbulence model
       ! if waves are active then calculate sediment mixing according to Van Rijn
       ! use Van Rijn's parabolic-linear mixing distribution for current-related
       ! mixing
       !
       epsmxc = vonkar*beta*ustarc
       if (tp>0. .and. wave) then
          !
          ! calculate sediment mixing due to waves following Van Rijn 1993
          !
          deltas = 3.0*delr
          deltas = min(0.2_fp, deltas)
          deltas = max(0.05_fp, deltas)
          epsbed = 0.004*dstar*deltas*uwb
          epsmax = 0.035*h1*1.4142*hrms/tp
       else
          deltas = 0.05_fp
          epsbed = 0.0
          epsmax = 0.0
       endif
       !
       ! set vertical sediment mixing values for waves and currents at water surface
       !
       epswav(0) = epsmax
       epscur(0) = 0.25*epsmxc*h1
       seddif(0) = sqrt(epscur(0)**2 + epswav(0)**2)
       !
       ! loop through vertical
       !
       do k = 1, kmax
          !
          ! calculate level of lower cell interface
          !
          lci = (1.0 + sig(k) - thick(k)/2.0)*h1
          if (lci>=(0.5*h1)) then
             epswav(k) = epsmax
             epscur(k) = 0.25*epsmxc*h1
          elseif (lci>deltas) then
             epswav(k) = epsbed + (epsmax - epsbed)*(lci - deltas)              &
                       & /(0.5*h1 - deltas)
             epscur(k) = epsmxc*lci*(1.0 - lci/h1)
          else
             epswav(k) = epsbed
             epscur(k) = epsmxc*lci*(1.0 - lci/h1)
          endif
          !
          ! set vertical sediment mixing values for waves and currents
          !
          seddif(k) = sqrt(epscur(k)**2 + epswav(k)**2)
       enddo
    else
       !
       ! set vertical sediment mixing values from K-epsilon turbulence model
       ! note beta factor should only be applied to current-related mixing
       ! this is rather rough method to estimate the proportion of the mixing
       ! that is due to current.
       !
       if (tauwav + tauc>eps) then
          betaef = (1.0 + (beta - 1.0)*tauc/(tauwav + tauc))
       else
          betaef = beta
       endif
       do k = 0, kmax
          seddif(k) = dicww(k)*betaef
       enddo
    endif
end subroutine calseddf1993
