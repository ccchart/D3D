!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id$
! $HeadURL$


module m_f1dimp_data
   !
   ! Contains variables and parameters related to flow 1d implicit module
   !
   use precision
   
   implicit none
   
   public f1dimppar_type
   
   type f1dimppar_type

      ! 
      !variables in <flwpar>
      !
      
      logical                          :: lconv                  !< Flag if simulation should continue without conververgence if max.number of iterations is reached. 
      
      integer                          :: flitmx                 !< maximum number of iterationsteps
      integer                          :: iterbc                 !< Maximum iterations to be performed in BICGST
      
      !don't use <real(fp)>. In <FLOWIT> it is defined as <real>
      real                             :: g                      !< Value for gravity acceleration
      real                             :: psi                    !< Spatial weigth factor in Preissmann scheme
      real                             :: theta                  !< Temporal weigth factor in Preissmann scheme
      real                             :: epsh                   !< Convergence criterium for water levels
      real                             :: epsq                   !< Convergence criterium for discharges(absolute)
      real                             :: rhow                   !< Density of water
      real                             :: omega                  !< Under relaxation coefficient
      real                             :: epsqrl                 !< Convergence criterium for discharges(relative)
      real                             :: lambda                 !< Extra resistance in general structure
      real                             :: relstr                 !< Under relaxation factor for structures
      real                             :: dhstru                 !< dh used for numerical differentation
      real                             :: cflpse                 !< (initial) pseudo Courant number
      real                             :: overlp                 !< sumerdike transition height
      real                             :: omcfl                  !< parameter for computing next value of dx/dt_pseudo
      real                             :: dhtyp                  !< parameter for computing next value of dx/dt_pseudo
      real                             :: exrstp                 !< Extra resistance in momentum equation
      
      double precision                 :: resid                  !< Allowable convergence measure for BICGST
      
      !
      !input variables to <SOFLOW>
      !
      
      logical                          :: steady                 !<  steady state flag
      
      integer                          :: istep                  !<  Current time step number (at t n+1 ) 
      !integer, dimension(2)            :: itim                   !<   Actual time level (at t n+1 ) expressed in date and time. Format:
      !                                                               !� itim(1) = YYYYMMDD (year,month,day)
      !                                                               !� itim(2) = HHMMSSHH (hour,minute,second, hundredth of a second)
      
      double precision                :: time                   !<  Actual time level (at t n+1 ) in seconds.
      double precision                :: dtf                    !<  Flow time step in seconds.     

      !
      !dimensions
      !
      
      integer                          :: ngrid                 !< Number of cells in network.
      integer                          :: ngridm                 !< Maximum number of cells in a branch
      integer                          :: nbran                 !<  Maximum number of connected branches to one node.
      integer                          :: maxlev                 !<  Maximum+1 number of nlev(1:ngrid).
      
      !
      !network
      !
      integer, allocatable, dimension(:,:) :: branch                !< branch : P, double(4, nbran ). Branch information.     
                                                                  !- branch(1,i) : integer. Node number n1 at begin of branch i.
                                                                  !- branch(2,i) : integer. Node number n2 at end of branch i.
                                                                  !- branch(3,i) : integer. Grid point i1 at begin of branch i.
                                                                  !- branch(4,i) : integer. Grid point i2 at end of branch i.
      
      !character(256)                   :: flbdfh                 !< File specifying Bedform-height
      !real(fp)      , dimension(:)    , pointer :: bedformD50    !< 50-percentile of sediment diameters (if no sediment simulated)

   end type f1dimppar_type
   
end module m_f1dimp_data