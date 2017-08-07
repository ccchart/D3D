subroutine periodic_xGyG_V(xGV,yGV,nlb,nub,mlb,mub,kmax, gdp)
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
!  $Id: My_intersec.f90 
!  $HeadURL:
!!--description-----------------------------------------------------------------
!
!   Function:    
!
!  Author: Alberto Canestrelli
!
!!--declarations----------------------------------------------------------------
!
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    integer              , pointer :: nrPER
    logical              , pointer :: twoCELLSperiod
    integer, dimension(:), pointer :: mPH_ext
    integer, dimension(:), pointer :: nPH_ext
    integer, dimension(:), pointer :: mPQ_ext
    integer, dimension(:), pointer :: nPQ_ext
    integer, dimension(:), pointer :: mPH_int
    integer, dimension(:), pointer :: mPQ_int
    integer, dimension(:), pointer :: mPH_intint
    integer, dimension(:), pointer :: mPQ_intint
    integer, dimension(:), pointer :: nPH_int
    integer, dimension(:), pointer :: nPQ_int
    integer, dimension(:), pointer :: nPH_intint
    integer, dimension(:), pointer :: nPQ_intint
    integer, dimension(:), pointer :: mPH_extext
    integer, dimension(:), pointer :: mPQ_extext
    integer, dimension(:), pointer :: nPH_extext
    integer, dimension(:), pointer :: nPQ_extext
    real(fp)             , pointer :: LchanPERprojX
    real(fp)             , pointer :: LchanPERprojY
!
! global variables
!
  real(fp)  , dimension(nlb:nub,mlb:mub)            ,intent(inout)   :: xGV  
  real(fp)  , dimension(nlb:nub,mlb:mub)            ,intent(inout)   :: yGV
  integer, intent(in)        :: nlb
  integer, intent(in)        :: nub
  integer, intent(in)        :: mlb
  integer, intent(in)        :: mub
  integer, intent(in)        :: kmax
!  integer, intent(in)        :: icx
!
! local variables
!
  integer                    :: k 
!
! executable statements -------------------------------------------------------
!                         
    nrPER          => gdp%gdimbound%nrPER
    twoCELLSperiod => gdp%gdimbound%twoCELLSperiod
    mPH_ext        => gdp%gdimbound%mPH_ext
    nPH_ext        => gdp%gdimbound%nPH_ext
    mPQ_ext        => gdp%gdimbound%mPQ_ext
    nPQ_ext        => gdp%gdimbound%nPQ_ext
    mPH_int        => gdp%gdimbound%mPH_int
    mPQ_int        => gdp%gdimbound%mPQ_int
    mPH_intint     => gdp%gdimbound%mPH_intint
    mPQ_intint     => gdp%gdimbound%mPQ_intint
    nPH_int        => gdp%gdimbound%nPH_int
    nPQ_int        => gdp%gdimbound%nPQ_int
    nPH_intint     => gdp%gdimbound%nPH_intint
    nPQ_intint     => gdp%gdimbound%nPQ_intint
    mPH_extext     => gdp%gdimbound%mPH_extext
    mPH_ext        => gdp%gdimbound%mPH_ext
    mPQ_extext     => gdp%gdimbound%mPQ_extext
    nPH_extext     => gdp%gdimbound%nPH_extext
    nPQ_extext     => gdp%gdimbound%nPQ_extext
    LchanPERprojX  => gdp%gdimbound%LchanPERprojX
    LchanPERprojY  => gdp%gdimbound%LchanPERprojY
      do k=0,nrPER+1  
         xGV(nPH_ext(k),mPH_ext(k)) = xGV(nPQ_int(k),mPQ_int(k)) + LchanPERprojX 
         xGV(nPQ_ext(k),mPQ_ext(k)) = xGV(nPH_int(k),mPH_int(k)) - LchanPERprojX   
         yGV(nPH_ext(k),mPH_ext(k)) = yGV(nPQ_int(k),mPQ_int(k)) + LchanPERprojY 
         yGV(nPQ_ext(k),mPQ_ext(k)) = yGV(nPH_int(k),mPH_int(k)) - LchanPERprojY   
         if (twoCELLSperiod) THEN
            xGV(nPH_extext(k),mPH_extext(k)) =  xGV(nPQ_intint(k),mPQ_intint(k)) + LchanPERprojX 
            xGV(nPQ_extext(k),mPQ_extext(k)) =  xGV(nPH_intint(k),mPH_intint(k)) - LchanPERprojX 
            yGV(nPH_extext(k),mPH_extext(k)) =  yGV(nPQ_intint(k),mPQ_intint(k)) + LchanPERprojY
            yGV(nPQ_extext(k),mPQ_extext(k)) =  yGV(nPH_intint(k),mPH_intint(k)) - LchanPERprojY
         endif
      enddo 
RETURN
end subroutine periodic_xGyG_V
