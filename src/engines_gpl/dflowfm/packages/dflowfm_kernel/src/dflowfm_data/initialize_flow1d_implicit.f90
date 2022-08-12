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

!> Initializes variables of flow1d implicit solver
subroutine initialize_flow1d_implicit(iresult)

!use m_flowparameters
use m_f1dimp
use m_physcoef
use m_flowgeom, only: ndx, ndxi
use unstruc_channel_flow, only: network
use m_flowexternalforcings
use unstruc_messages
use m_flow, only: s0, ucmag, hs

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
      
!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: k
integer :: k2
integer :: idx_crs

!!
!! CALC
!!

iresult=0

!<flwpar>
f1dimppar%g=ag 
f1dimppar%psi=0.5d0
f1dimppar%theta=1.0d0 !I think it is rewriten if steady flow is chosen anyhow
f1dimppar%epsh=1.0d-10 
f1dimppar%epsq=1.0d-10 
f1dimppar%rhow=rhomean 
f1dimppar%flitmx=10
f1dimppar%omega=0.5d0 !check sensible value
f1dimppar%epsqrl=1.0d-10
f1dimppar%lambda=0 
f1dimppar%relstr=1.0d0
f1dimppar%dhstru=1.0d-5
f1dimppar%cflpse=1000.0d0
f1dimppar%iterbc=10
f1dimppar%resid=1.0d-8 !check sensible value
f1dimppar%overlp=0 !change to summerdiketransitionheight -> check if <network%csdefinitions%cs(1)%summerdike> allocated?
f1dimppar%lconv=1 !the input is converted to logical by calling soipar? setting to true we can break the simulation in FM code to handle the messages
f1dimppar%omcfl=0.9d0 !default in SRE
f1dimppar%dhtyp=0.1d0 !default in SRE
f1dimppar%exrstp=0.0d0 !default in SRE

!<SOFLOW> input
f1dimppar%steady=.true.

!dimensions
f1dimppar%ngrid=network%numk !total number of mesh nodes (water level points)
f1dimppar%nbran=network%brs%count !Fine if there is only one branch. Otherwise compute.

f1dimppar%ngridm=0
do k=1,f1dimppar%nbran
    f1dimppar%ngridm=f1dimppar%ngridm+network%BRS%BRANCH(k)%GRIDPOINTSCOUNT
end do

f1dimppar%maxlev=0 
do k=1,network%CSDEFINITIONS%COUNT 
    f1dimppar%maxlev=max(f1dimppar%maxlev,network%CSDEFINITIONS%CS(1)%LEVELSCOUNT)
end do

f1dimppar%nnode=network%nds%count 
f1dimppar%nhstat=nzbnd 
f1dimppar%nqstat=nqbnd 
f1dimppar%ntabm=8 !Properly compute. I cannot find the location of the read BC tables!
f1dimppar%maxtab=ndx - ndxi !<we have as many tables as open boundaries
f1dimppar%nbrnod=network%NDS%MAXNUMBEROFCONNECTIONS

!!dependent on branch
allocate(f1dimppar%branch(4,f1dimppar%nbran)) 
allocate(f1dimppar%bfrict(3,f1dimppar%nbran))
do k=1,f1dimppar%nbran
    !branch
    f1dimppar%branch(1,k)=network%BRS%BRANCH(k)%NODEINDEX(1)
    f1dimppar%branch(2,k)=network%BRS%BRANCH(k)%NODEINDEX(2)
    f1dimppar%branch(3,k)=network%BRS%BRANCH(k)%FROMNODE%GRIDNUMBER
    f1dimppar%branch(4,k)=network%BRS%BRANCH(k)%TONODE%GRIDNUMBER
    
    !bfrict
    do k2=1,3 !< main channel, floodplain 1, floodplain 2
        select case (network%RGS%ROUGH(k2)%FRICTIONTYPE) !< where is the information per branch? add when several branches!
            case (0)
                f1dimppar%bfrict(k2,k)=1
            case default
                write (msgbuf, '(a)') 'Only constant Chezy friction is supported at the moment.'
                call err_flush()
                iresult=1
        end select
    enddo
!
!bfrict(1,i)=cfrchc (1) : Chézy constant
!bfrict(1,i)=cfrchq (2) : Chézy function of discharge
!bfrict(1,i)=cfrchh (3) : Chézy function of water level
!bfrict(1,i)=cfrman (4) : Manning constant
!bfrict(1,i)=cfrskn (5) : Strickler 1 constant ( k n )
!bfrict(1,i)=cfrsks (6) : Strickler 2 constant ( k s )
!bfrict(1,i)=cfrnik (7) : Nikuradze constant
!bfrict(1,i)=cfreng (8) : Engelund predicto

end do
 

!dependent on gridpoints 
allocate(f1dimppar%bfricp(6,f1dimppar%ngrid)) !needs the part with FP1, FP2
allocate(f1dimppar%x(f1dimppar%ngrid))
allocate(f1dimppar%nlev(f1dimppar%ngrid)) 
allocate(f1dimppar%hpack(3,f1dimppar%ngrid)) 
allocate(f1dimppar%qpack(3,f1dimppar%ngrid))

    !cross-sectional information (gridpoint,level)
allocate(f1dimppar%wft(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%aft(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%wtt(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%att(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%of(f1dimppar%ngrid,f1dimppar%maxlev)) 
allocate(f1dimppar%hlev(f1dimppar%ngrid,f1dimppar%maxlev))

do k=1,f1dimppar%ngrid
    idx_crs=network%ADM%gpnt2cross(k)%C1 !< index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node 		
    
    !bfrictp
    f1dimppar%bfricp(1,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    f1dimppar%bfricp(2,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    !deal properly with the values below when friction per section varies. 
    f1dimppar%bfricp(3,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    f1dimppar%bfricp(4,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    f1dimppar%bfricp(5,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    f1dimppar%bfricp(6,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
!bfricp(6,ngrid)   I  Bed friction parameters:
!                     (1,i) = Parameter for positive flow direction
!                             in main section (depending on friction
!                             type):
!                             =     Chezy constant value
!                             =     Table pointer (Q or h table)
!                             =     Nikuradse parameter kn for Ni-
!                                   kuradse formula
!                             =     Manning parameter nm for Manning
!                                   formula
!                             =     Strickler coefficient ks for
!                                   Strickler formula
!                     (2,i) = Parameter for negative flow direction
!                             in main section (depending on friction
!                             type) Same definitions as bfricp(1,i).
!                     (3,i) = Parameter for positive flow direction
!                             in sub sec 1 (depending on friction
!                             type):
!                             =     Chezy constant value
!                             =     Nikuradse parameter kn for Niku-
!                                   radse formula
!                             =     Manning parameter nm for Manning
!                                   formula
!                             =     Strickler coefficient ks for
!                                   Strickler formula
!                     (4,i) = Parameter for negative flow direction
!                             in sub sec 1 (depending on friction
!                             type) Same definition as bfricp (3,i):
!                     (5,i) = Parameter for positive flow direction
!                             in sub sec 2 (depending on friction
!                             type) Same definition as bfricp (3,i).
!                     (6,i) = Parameter for negative flow direction
!                             in sub sec 2 (depending on friction
!                             type) Same definition as bfricp (3,i).
    
    !x
    f1dimppar%x(k)=network%CRS%CROSS(k)%CHAINAGE
    
    !nlev
    f1dimppar%nlev(k)=network%CRS%CROSS(101)%TABDEF%LEVELSCOUNT
    do k2=1,f1dimppar%nlev(k)
        f1dimppar%wft(k,k2)=network%CRS%CROSS(k)%TABDEF%FLOWWIDTH(k2) 
        f1dimppar%aft(k,k2)=network%CRS%CROSS(k)%TABDEF%FLOWAREA(k2)  
        f1dimppar%wtt(k,k2)=network%CRS%CROSS(k)%TABDEF%TOTALWIDTH(k2)   
        f1dimppar%att(k,k2)=network%CRS%CROSS(k)%TABDEF%TOTALAREA(k2)    
        f1dimppar%of(k,k2)=network%CRS%CROSS(k)%TABDEF%WETPERIMETER(k2)     
        f1dimppar%hlev(k,k2)=network%CRS%CROSS(k)%TABDEF%HEIGHT(k2)
    end do !k2

    !dependent variables
    do k2=1,3 !< time step before, intermediate, after
        !I don't think <k> below is correct. It should be an index mapping gridpoint and internal gridpoint
        f1dimppar%hpack(k2,k)=s0(k)
        f1dimppar%qpack(k2,k)=hs(k)*ucmag(k)
    end do !k2
end do !k

! 
!boundary conditions
!   h
allocate(f1dimppar%hbdpar(3,f1dimppar%nhstat)) 
do k=1,f1dimppar%nhstat
    f1dimppar%hbdpar(1,k)=kbndz(2,1) !< first s1 point on the inside of the domain
    f1dimppar%hbdpar(2,k)=1
    f1dimppar%hbdpar(3,k)=k
end do
!   q
allocate(f1dimppar%qbdpar(3,f1dimppar%nqstat)) 
do k=1,f1dimppar%nqstat
    f1dimppar%qbdpar(1,k)=kbndu(2,1) !< first s1 point on the inside of the domain
    f1dimppar%qbdpar(2,k)=1
    f1dimppar%qbdpar(3,k)=f1dimppar%nhstat+k !< table number after the ones of <hbdpar>
end do

     !hbdpar(3,nhstat)  Hydrodynamic conditions for H-stations:
     !    (1,i) = Location [grid point] for H-station.
     !    (2,i) = Type of condition
     !            cbftim (1) : h = f(t)
     !            cbfqoh (2) : h = h(Q)
     !            cbfour (3) : h = fourier
     !            cbtidl (4) : h = tidal components
     !    (3,i) = Table number for f(t), h(Q), fourier
     !            or tidal components table.
      

!tables
!I cannot find the location of BC. Hardcoded for now:
allocate(f1dimppar%table(f1dimppar%ntabm)) 
f1dimppar%table=(/ 0,10000,0,86400,100,100,1,2 /)

!                     call INTTAB (ntab(1,itab), ntab(4,itab),
!     +                            table(ntab(3,itab)),
!     +                            table(ntab(2,itab)),
!     +                            dble(watlev),qh    )
     
allocate(f1dimppar%ntab(4,f1dimppar%maxtab)) 
!   h
f1dimppar%ntab(1,1)=2
f1dimppar%ntab(2,1)=1
f1dimppar%ntab(3,1)=1
f1dimppar%ntab(4,1)=0
!   q
f1dimppar%ntab(1,2)=2
f1dimppar%ntab(2,2)=3
f1dimppar%ntab(3,2)=3
f1dimppar%ntab(4,2)=0

!nodes
allocate(f1dimppar%node(4,f1dimppar%nnode))
allocate(f1dimppar%numnod(f1dimppar%nnode))
!everything should be in the loop, but there are some issues...
!upstream
k=1
f1dimppar%node(1,k)=3
f1dimppar%node(3,k)=1
f1dimppar%node(4,k)=1 
!downstream
k=2
f1dimppar%node(1,k)=2
f1dimppar%node(3,k)=2
f1dimppar%node(4,k)=2
do k=1,f1dimppar%nnode
    !f1dimppar%node(1,k)= 
    !for some reason at this stage <network%NDS%NODE(2)%NODETYPE> is 0 for all cases. 
                                                         !! - -2    boundary node
                                                         !! - -1    not set
                                                         !! -  0    node with one reach connected
                                                         !! -  1    connection node with more than one reach connected
                                                         !! -  2    water level boundary
                                                         !! -  3    Discharge boundary 
                                                         !! -  4    Discharge boundary as tabulated function of water level
                                                         !! -  5    Embedded node
    f1dimppar%node(2,k)=network%NDS%NODE(k)%GRIDNUMBER
    
    f1dimppar%numnod(k)=network%NDS%NODE(1)%NUMBEROFCONNECTIONS
end do
!

end subroutine initialize_flow1d_implicit

