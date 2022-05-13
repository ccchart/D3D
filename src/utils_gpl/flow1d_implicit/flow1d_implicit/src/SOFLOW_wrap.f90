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

!#BEGIN# MOVE TO INITIALIZATION PARAMETERS

!<flwpar>
f1dimppar%g=9.81d0 !read from FM
f1dimppar%psi=0.5d0
f1dimppar%theta=1.0d0 !I think it is rewriten if steady flow is chosen anyhow
f1dimppar%epsh=1.0d-10 
f1dimppar%epsq=1.0d-10 
f1dimppar%rhow=1000.0d0 !read from FM
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

!<SOFLOW> input
f1dimppar%steady=.true.

!dimensions
f1dimppar%ngrid=10 !read from FM
f1dimppar%ngridm=10 !for one branch it is fine if it is the same as <ngrid>. Otherwise compute.
f1dimppar%nbran=1 !Fine if there is only one branch. Otherwise compute.
f1dimppar%maxlev=20 !Properly compute based on cross-section data. 
f1dimppar%nnode=2 !Properly compute based on cross-section data. 
f1dimppar%nhstat=1 !Properly compute based on cross-section data. 
f1dimppar%nqstat=1 !Properly compute based on cross-section data. 
f1dimppar%ntabm=10 !Properly compute based on cross-section data. 
f1dimppar%maxtab=2 !Properly compute based on cross-section data. 
f1dimppar%nbrnod=3

!dependent on branch
allocate(f1dimppar%branch(4,nbran)) !deal with allocate and deallocate properly
allocate(f1dimppar%bfrict(3,nbran)) 

!dependent on gridpoints 
allocate(f1dimppar%bfricp(6,f1dimppar%ngrid)) !deal with allocate and deallocate properly
allocate(f1dimppar%hpack(3,f1dimppar%ngrid)) 
allocate(f1dimppar%qpack(3,f1dimppar%ngrid)) 
allocate(f1dimppar%x(f1dimppar%ngrid)) 
allocate(f1dimppar%nlev(f1dimppar%ngrid)) 

!cross-sectional information
allocate(f1dimppar%wft(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%aft(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%wtt(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%att(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%of(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%hlev(f1dimppar%ngrid,f1dimppar%maxlev)) 
 
!boundary conditions
allocate(f1dimppar%hbdpar(3,f1dimppar%nhstat)) !deal with allocate and deallocate properly
allocate(f1dimppar%qbdpar(3,f1dimppar%nqstat)) 

!tables
allocate(f1dimppar%table(f1dimppar%ntabm)) 
allocate(f1dimppar%ntab(4,f1dimppar%maxtab)) 

!nodes
allocate(f1dimppar%node(4,f1dimppar%nnode))
allocate(f1dimppar%numnod(f1dimppar%nnode))

!#END# MOVE TO INITIALIZATION PARAMETERS

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
node   => f1dimppar%ntab
numnod => f1dimppar%numnod

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
        &   bfricp , hpack  , qpack                    , & 
!cross-sectional shape
        &   wft    , aft    ,wtt     ,att     , of     , & 
        &   hlev                                       , &
!boundary conditions
        &   hbdpar , qbdpar                            , &
!tables
        &   table  , ntab                              , &
!dependent on node
        &   node   , numnod                              &
!close
        &)
    
        
end subroutine SOFLOW_wrap