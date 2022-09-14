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
    
subroutine SOFLOW_wrap(ndx,lnx,s0,umag,au,wu,s1,u1)   

use m_f1dimp
!use m_flow, only: s1, u1 !this subroutine is not in flow_kernel and cannot use this module

implicit none

!
!input
!

integer, intent(in) :: ndx
integer, intent(in) :: lnx

double precision, dimension(ndx), intent(in) :: s0
double precision, dimension(ndx), intent(in) :: umag
double precision, dimension(lnx), intent(in) :: au
double precision, dimension(lnx), intent(in) :: wu

!
!output
!

double precision, dimension(ndx), intent(out) :: s1
double precision, dimension(lnx), intent(out) :: u1

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
     

!local
integer                              :: k, k2

real, dimension(:,:), allocatable :: waoft

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

!time dependent
allocate(waoft(ngrid,6))

do k=1,ngrid
    !I don't think <k> below is correct. It should be an index mapping gridpoint and internal gridpoint
    !FM1DIMP2DO: needs to be interpolated from link to cell centre
    !FM1DIMP2DO: needs to be separated between flow and total
    !FM1DIMP2DO: wetted perimeter get from results
    !check right order in <FLNORM> and not in documentation. 
    waoft(k,1)=wu(k) !wf = actual flow width 
    waoft(k,2)=wu(k) !wt = actual total width
    waoft(k,3)=au(k) !af = actual flow area
    waoft(k,4)=au(k) !at = actual total area n
    waoft(k,5)=au(k) !at = actual total area n+1
    waoft(k,6)=au(k)/wu(k) !o = actual wetted perimeter
enddo

!FM1DIMP2D: <au> needs to be interpolated from link to cell centre.
!FM1DIMP2D: check what happens with several branches. We have to get rid of ghosts. 
do k=1,3
    hpack(:,k)=s0(1:ngrid)
    do k2=1,ngrid
        qpack(k2,k)=umag(k2)*au(k2) !incorrect, see above. 
    enddo
    !qpack(:,k)=100 !check if hampers iteration
enddo

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
        &   bfricp , hpack  , qpack  ,x       , waoft  , & 
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
    
     
!Interpolate at links
!FM1DIMP2D: compute new flow area <au> to find right velocity before interpolating
!FM1DIMP2D: check what happens with several branches. We have to get rid of ghosts.
s1=hpack(:,3);
u1=qpack(:,3)/au(1:ngrid-1)

end subroutine SOFLOW_wrap