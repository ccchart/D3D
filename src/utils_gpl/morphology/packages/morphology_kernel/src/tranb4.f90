subroutine tranb4(utot      ,d         ,c         ,npar      ,par       , &
                & hidexp    ,sbot      ,ssus      )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2021.                                
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
! computes sediment transport according to
! general formula
! -
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
!
! Arguments
!
    integer                  , intent(in)    :: npar
    real(fp)                 , intent(in)    :: c
    real(fp)                 , intent(in)    :: d
    real(fp)                 , intent(in)    :: hidexp !< hiding & exposure factor
    real(fp), dimension(npar), intent(in)    :: par
    real(fp)                 , intent(in)    :: utot
    !
    real(fp)                 , intent(out)   :: sbot
    real(fp)                 , intent(out)   :: ssus
!
! Local variables
!
    real(fp) :: acal
    real(fp) :: ag    ! gravity acceleration
    real(fp) :: b     ! correction coefficient shear stress
    real(fp) :: cc
    real(fp) :: delta ! velocity (es/ew)  relative density of sediment particle
    real(fp) :: f     ! real help array
    real(fp) :: rmu
    real(fp) :: th
    real(fp) :: thcr
!
!! executable statements -------------------------------------------------------
!
    sbot  = 0.0
    ssus  = 0.0
    !
    ag    = par( 1)
    delta = par( 4)
    acal  = par(11)
    b     = par(12)
    cc    = par(13)
    rmu   = par(14)
    thcr  = par(15)
    !
    if ((c<1.0e-6) .or. (utot<1.0e-6)) then
       return
    endif
    th = (utot/c)**2/(delta*d)
    f  = rmu*th - hidexp*thcr
    if (f>1.0e-8) sbot = acal*d**1.5*sqrt(ag*delta)*th**b*f**cc
    ssus = 0.0
end subroutine tranb4
