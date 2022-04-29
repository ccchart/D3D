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

!> Performs a single computational timestep, calling <SOFLOW> of Sobek-RE
    
subroutine SOFLOW_wrap()                
use m_f1dimp
implicit none

!
!pointer
!

logical                          , pointer :: lconv                   

integer                          , pointer :: flitmx                 
integer                          , pointer :: iterbc                 

real                             , pointer :: g
!real(fp)                         , pointer :: g                      
real(fp)                         , pointer :: psi                    
real(fp)                         , pointer :: theta                  
real(fp)                         , pointer :: epsh                   
real(fp)                         , pointer :: epsq                   
real(fp)                         , pointer :: rhow                   
real(fp)                         , pointer :: omega                  
real(fp)                         , pointer :: epsqrl                 
real(fp)                         , pointer :: lambda                 
real(fp)                         , pointer :: relstr                 
real(fp)                         , pointer :: dhstru                 
real(fp)                         , pointer :: cflpse                 
real(fp)                         , pointer :: resid                  
real(fp)                         , pointer :: overlp                 
real(fp)                         , pointer :: omcfl                  
real(fp)                         , pointer :: dhtyp                  
real(fp)                         , pointer :: exrstp                 
      


!real(fp)      , dimension(:)         , pointer :: sedd50
!
!f1dimp variables
!

!#BEGIN# MOVE TO INITIALIZATION PARAMETERS

f1dimppar%g=9.81d0
f1dimppar%psi=0.5d0
f1dimppar%theta=1.0d0 !I think it is rewriten if steady flow is chosen anyhow
f1dimppar%epsh=1.0d-10 
f1dimppar%epsq=1.0d-10 
f1dimppar%rhow=1000.0d0 
f1dimppar%flitmx=10
f1dimppar%omega=0.5d0 !check sensible value
f1dimppar%epsqrl=1.0d-10
f1dimppar%lambda=0 
f1dimppar%relstr=1.0d0
f1dimppar%dhstru=1.0d-5
f1dimppar%cflpse=1000.0d0
f1dimppar%iterbc=10
f1dimppar%resid=1.0d-8 !check sensible value
f1dimppar%overlp=0 !change to summerdiketransitionheight
f1dimppar%lconv=1 !the input is converted to logical by calling soipar? setting to true we can break the simulation in FM code to handle the messages
f1dimppar%omcfl=0.9d0 !default in SRE
f1dimppar%dhtyp=0.1d0 !default in SRE
f1dimppar%exrstp=0.0d0 !default in SRE

f1dimppar%istep=1
f1dimppar%itim(1)=20000101
f1dimppar%itim(2)=00000000
f1dimppar%steady=.true.
f1dimppar%time=1.0d0
f1dimppar%dtf=1.0d0

!#END# MOVE TO INITIALIZATION PARAMETERS

g => f1dimppar%g

call SOFLOW(g)
    
        
end subroutine SOFLOW_wrap