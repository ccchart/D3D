subroutine trab14(kode      ,ntrsi     ,utot      ,d50       ,chezy     , &
                & par       ,hidexp    ,sbot      ,ssus      )
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
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
! computes sediment transport according to
! generalized version of Ashida and Michiue
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
!
! Global variables
!
    integer                , intent(in)  :: kode   ! indicates active grid point
                                                   !    = 0 active point
                                                   !    < 0 point not used now
    integer                , intent(out) :: ntrsi  ! indicator for output representation
                                                   !    of sediment transport.
                                                   !    =1 :  magnitude of transport
                                                   !    =2 :  components of transport
    real(fp)               , intent(in)  :: chezy  ! chezy value
    real(fp)               , intent(in)  :: d50    ! Grain size specified as d50
    real(fp)               , intent(in)  :: hidexp ! hiding & exposure factor
    real(fp)               , intent(out) :: sbot   ! bed load transport
    real(fp)               , intent(out) :: ssus   ! suspended sediment transport
    real(fp)               , intent(in)  :: utot   ! flow velocity
    real(fp), dimension(30), intent(in)  :: par    ! sediment parameter list
!
! Local variables
!
    real(fp) :: a
    real(fp) :: ag
    real(fp) :: delta  ! relative density of sediment particle
    real(fp) :: m
    real(fp) :: p
    real(fp) :: q
    real(fp) :: sag
    real(fp) :: sgd
    real(fp) :: ssgd3
    real(fp) :: t
    real(fp) :: tc
    real(fp) :: tct
    real(fp) :: ustar
!
!! executable statements -------------------------------------------------------
!
    ag    = par(1)
    !        rhosol = par( 2)
    !        rhow   = par( 3)
    delta = par(4)      ! (rhosol - rhow) / rhow
    a     = par(11)     ! acal: tuning constant of sediment transport
    tc    = par(12)
    m     = par(13)
    p     = par(14)
    q     = par(15)
    sag   = sqrt(ag)
    sgd   = delta*ag*d50
    ssgd3 = sqrt(sgd*d50*d50)
    ntrsi = 1
    sbot  = 0.0
    ssus  = 0.0
    if (kode /= - 1) then
       if (chezy >= 1.e-6) then
          if (utot >= 1.e-6) then
             ustar = sag*utot/chezy
             t     = ustar**2/sgd
             tct   = hidexp*tc/t
             if (tct < 1.0) then
                sbot = a * ssgd3 * t**m * (1.0 - tct)**p * (1.0 - sqrt(tct))**q
             endif
          endif
       endif
    endif
end subroutine trab14
