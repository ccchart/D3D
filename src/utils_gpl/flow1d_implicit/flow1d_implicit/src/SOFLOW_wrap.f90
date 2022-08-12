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

logical                                  , pointer :: lconv                   
logical                                  , pointer :: steady    
                                         
integer                                  , pointer :: flitmx                 
integer                                  , pointer :: iterbc                 
integer                                  , pointer :: ngrid   
integer                                  , pointer :: ngridm   
integer                                  , pointer :: nbran   
integer                                  , pointer :: maxlev
integer                                  , pointer :: nnode
integer                                  , pointer :: nhstat
integer                                  , pointer :: nqstat
integer                                  , pointer :: maxtab
integer                                  , pointer :: ntabm
integer                                  , pointer :: nbrnod

integer, dimension(:)                    , pointer :: nlev
integer, dimension(:)                    , pointer :: numnod

integer, dimension(:,:)                  , pointer :: branch
integer, dimension(:,:)                  , pointer :: bfrict
integer, dimension(:,:)                  , pointer :: hbdpar
integer, dimension(:,:)                  , pointer :: qbdpar
integer, dimension(:,:)                  , pointer :: ntab
integer, dimension(:,:)                  , pointer :: node
integer, dimension(:,:)                  , pointer :: nodnod

real                                     , pointer :: g
real                                     , pointer :: psi                    
real                                     , pointer :: theta                  
real                                     , pointer :: epsh                   
real                                     , pointer :: epsq                   
real                                     , pointer :: rhow                   
real                                     , pointer :: omega                  
real                                     , pointer :: epsqrl                 
real                                     , pointer :: lambda                 
real                                     , pointer :: relstr                 
real                                     , pointer :: dhstru                 
real                                     , pointer :: cflpse                               
real                                     , pointer :: overlp                 
real                                     , pointer :: omcfl                  
real                                     , pointer :: dhtyp                  
real                                     , pointer :: exrstp     

real, dimension(:)                       , pointer :: table
real, dimension(:)                       , pointer :: x

real, dimension(:,:)                     , pointer :: bfricp
real, dimension(:,:)                     , pointer :: wft
real, dimension(:,:)                     , pointer :: aft
real, dimension(:,:)                     , pointer :: wtt
real, dimension(:,:)                     , pointer :: att
real, dimension(:,:)                     , pointer :: of


double precision                         , pointer :: time
double precision                         , pointer :: dtf
double precision                         , pointer :: resid

double precision, dimension(:,:)         , pointer :: hpack
double precision, dimension(:,:)         , pointer :: qpack
double precision, dimension(:,:)         , pointer :: hlev
      
!
!f1dimp variables
!

!#BEGIN# MOVE TO CONVERSION ROUTINE FOR EVERY TIME STEP

!<SOFLOW> input
f1dimppar%istep=1 
!<itim> only used for writing to file when error
!f1dimppar%itim(1)=20000101 
!f1dimppar%itim(2)=00000000
f1dimppar%time=1.0d0
f1dimppar%dtf=1.0d0

!deal with 
!hpack
!qpack

!#END# MOVE TO CONVERSION ROUTINE FOR EVERY TIME STEP

!<flwpar> variables
g      => f1dimppar%g
psi    => f1dimppar%psi
theta  => f1dimppar%theta
epsh   => f1dimppar%epsh
epsq   => f1dimppar%epsq
rhow   => f1dimppar%rhow
omega  => f1dimppar%omega
epsqrl => f1dimppar%epsqrl
lambda => f1dimppar%lambda
relstr => f1dimppar%relstr
dhstru => f1dimppar%dhstru
cflpse => f1dimppar%cflpse
resid  => f1dimppar%resid
overlp => f1dimppar%overlp
omcfl  => f1dimppar%omcfl
dhtyp  => f1dimppar%dhtyp
exrstp => f1dimppar%exrstp
flitmx => f1dimppar%flitmx

!<SOFLOW> variables
time   => f1dimppar%time
dtf    => f1dimppar%dtf
steady => f1dimppar%steady

!dimensions
ngrid  => f1dimppar%ngrid
ngridm => f1dimppar%ngridm
nbran  => f1dimppar%nbran
maxlev => f1dimppar%maxlev
nnode  => f1dimppar%nnode 
nhstat => f1dimppar%nhstat 
nqstat => f1dimppar%nqstat
maxtab => f1dimppar%maxtab
ntabm  => f1dimppar%ntabm
nbrnod => f1dimppar%nbrnod
nlev   => f1dimppar%nlev

!dependent on branch
branch => f1dimppar%branch
bfrict => f1dimppar%bfrict

!dependent on gridpoints 
bfricp => f1dimppar%bfricp
hpack  => f1dimppar%hpack
qpack  => f1dimppar%qpack
x      => f1dimppar%x

!cross-sectional shape
wft  => f1dimppar%wft 
aft  => f1dimppar%aft 
wtt  => f1dimppar%wtt 
att  => f1dimppar%att 
of   => f1dimppar%of  
hlev => f1dimppar%hlev

!boundary conditions
hbdpar => f1dimppar%hbdpar
qbdpar => f1dimppar%qbdpar

!tables
table  => f1dimppar%table
ntab   => f1dimppar%ntab

!dependent on node
node   => f1dimppar%node
numnod => f1dimppar%numnod
nodnod => f1dimppar%nodnod

call SOFLOW( &
!<flwpar> input
        &   g      , psi    , theta  , epsh   , epsq   , &
        &   rhow   , omega  , epsqrl , lambda , relstr , &
        &   dhstru , cflpse , resid  , overlp , omcfl  , &
        &   dhtyp  , exrstp , flitmx                   , &
!<SOFLOW> input
        &   time   , dtf    , steady                   , &
!dimensions 
        &   ngrid  , ngridm , nbran  , maxlev , nnode  , &
        &   nhstat , nqstat , maxtab , ntabm  , nbrnod , &
        &   nlev                                       , &
!dependent on branch
        &   branch , bfrict                            , &
!dependent on gridpoints 
        &   bfricp , hpack  , qpack  ,x                , & 
!cross-sectional shape
        &   wft    , aft    ,wtt     ,att     , of     , & 
        &   hlev                                       , &
!boundary conditions
        &   hbdpar , qbdpar                            , &
!tables
        &   table  , ntab                              , &
!dependent on node
        &   node   , numnod ,nodnod                    ,  &
!close
        &)
    
        
end subroutine SOFLOW_wrap