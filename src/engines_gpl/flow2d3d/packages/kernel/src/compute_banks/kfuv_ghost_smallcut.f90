subroutine kfuv_ghost_smallcut(kfu,kfv,nlb,nub,mlb,mub,constant, gdp)
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
!  Ndryact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: My_intersec.f90 
!  $HeadURL:
!!--description-----------------------------------------------------------------
!
!   Function:  deactivate small partially cut edges before uzd to turn off update there, and reactivate it before sud in order to compute
!              fluxes
!
!   Author: Alberto Canestrelli
!
!!--declarations----------------------------------------------------------------
!
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    integer                 , pointer :: totGHOSTu1
    integer                 , pointer :: totGHOSTv1
    integer, dimension(:)   , pointer :: nGPu1
    integer, dimension(:)   , pointer :: mGPu1
    integer, dimension(:)   , pointer :: nGPv1
    integer, dimension(:)   , pointer :: mGPv1
    real(fp), dimension(:,:), pointer :: aguu
    real(fp), dimension(:,:), pointer :: agvv
!
! global variables
!
    integer, dimension(nlb:nub,mlb:mub)                                 , intent(inOUT) :: kfu 
    integer, dimension(nlb:nub,mlb:mub)                                 , intent(inOUT) :: kfv
    integer                                                             , intent(in)    :: nlb
    integer                                                             , intent(in)    :: nub
    integer                                                             , intent(in)    :: mlb
    integer                                                             , intent(in)    :: mub
    integer                                                             , intent(in)    :: constant
!
!
! local variables
!
  integer                    :: I,mGP,nGP
!
! executable statements -------------------------------------------------------
!
    totGHOSTu1 => gdp%gdimbound%totGHOSTu1
    totGHOSTv1 => gdp%gdimbound%totGHOSTv1
    nGPu1      => gdp%gdimbound%nGPu1
    mGPu1      => gdp%gdimbound%mGPu1
    nGPv1      => gdp%gdimbound%nGPv1
    mGPv1      => gdp%gdimbound%mGPv1
    aguu       => gdp%gdimbound%aguu
    agvv       => gdp%gdimbound%agvv
!
! set u-velocity mask to constant
!
    do i = 1,totGHOSTu1   
       mGP = mGPu1(i)
       nGP = nGPu1(i) 
       if  (comparereal(aguu(nGP,mGP),0.5_fp).le.0.and.comparereal(aguu(nGP,mGP),0._fp).gt.0)  then   !i dont want to change it if aguu=0
         kfu(nGP,mGP) = constant 
       endif
    enddo
!
! set v-velocity mask to constant
!
    do i = 1,totGHOSTv1   
       mGP = mGPv1(i)
       nGP = nGPv1(i) 
       if  (comparereal(agvv(nGP,mGP),0.5_fp).le.0.and.comparereal(agvv(nGP,mGP),0._fp).gt.0)  then   !i dont want to change it if agvv=0
         kfv(nGP,mGP) = constant 
       endif
    enddo
!
    return 
end subroutine kfuv_ghost_smallcut
