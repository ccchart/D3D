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
use m_flowgeom, only: ndx, ndxi, wu, teta, lnx, lnx1D, lnx1Db, ln, lnxi
use unstruc_channel_flow, only: network
use m_flowexternalforcings !FM1dIMP2DO: do I need it?
use unstruc_messages
use m_flow, only: s0, u1, au !<ucmag> is velocity at cell centres, but we initialize <u1>
use m_sediment, only: stmpar, jased, stm_included
use m_fm_erosed, only: link1sign2
!use m_meteo !boundary conditions -> eventually done every time step

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
integer                                  , pointer :: table_length

integer, dimension(:)                    , pointer :: nlev
integer, dimension(:)                    , pointer :: numnod
integer, dimension(:)                    , pointer :: grd_sre_fm

integer, dimension(:,:)                  , pointer :: grd_fmL_sre

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
real, dimension(:,:)                     , pointer :: waoft

double precision                         , pointer :: time
double precision                         , pointer :: dtf
double precision                         , pointer :: resid

double precision, dimension(:,:)         , pointer :: hpack
double precision, dimension(:,:)         , pointer :: qpack
double precision, dimension(:,:)         , pointer :: hlev

!debug
integer, pointer :: fm1dimp_debug_k1

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: k
integer :: k2
integer :: idx_crs
integer :: n1, n2, nint, nout
integer :: table_number
integer :: idx_fr, idx_to
integer :: idx_i, idx_f, idx_fm, nl, L
!integer :: nlink !I don't think I need it global

integer, allocatable, dimension(:)   :: kcol

!integer, allocatable, dimension(:,:) :: nodnod_aux

real :: swaoft


!point
!pointer cannot be before the array is allocated, here only non-allocatable arrays

table_length => f1dimppar%table_length
maxtab       => f1dimppar%maxtab
nnode        => f1dimppar%nnode
ntabm        => f1dimppar%ntabm
nbran        => f1dimppar%nbran
ngrid        => f1dimppar%ngrid
nbrnod       => f1dimppar%nbrnod
maxlev       => f1dimppar%maxlev
ngridm       => f1dimppar%ngridm
nhstat       => f1dimppar%nhstat
nqstat       => f1dimppar%nqstat

!!
!! CALC
!!

iresult=0

!<flwpar>
f1dimppar%g=ag 
f1dimppar%psi=0.5d0
f1dimppar%theta=1.0d0 !I think it is rewriten if steady flow is chosen anyhow
f1dimppar%epsh=1.0d-6 !FM1DIMP2DO: make input 
f1dimppar%epsq=1.0d-6 !FM1DIMP2DO: make input 
f1dimppar%rhow=rhomean 
f1dimppar%flitmx=100 !FM1DIMP2DO: make input 
f1dimppar%omega=0.5d0 !check sensible value
f1dimppar%epsqrl=1.0d-10
f1dimppar%lambda=0 
f1dimppar%relstr=1.0d0
f1dimppar%dhstru=1.0d-5
f1dimppar%cflpse=1000.0d0
f1dimppar%iterbc=100
f1dimppar%resid=1.0d-8 !check sensible value
f1dimppar%overlp=0 !change to summerdiketransitionheight -> check if <network%csdefinitions%cs(1)%summerdike> allocated?
f1dimppar%lconv=1 !the input is converted to logical by calling soipar? setting to true we can break the simulation in FM code to handle the messages
f1dimppar%omcfl=0.9d0 !default in SRE
f1dimppar%dhtyp=0.1d0 !default in SRE
f1dimppar%exrstp=0.0d0 !default in SRE

!<SOFLOW> input
f1dimppar%steady=.true.

!dimensions
nbran=network%brs%count 

ngrid=0
ngridm=0
!nlink
do k=1,nbran
    ngrid=ngrid+network%BRS%BRANCH(k)%GRIDPOINTSCOUNT
    ngridm=max(ngridm,network%BRS%BRANCH(k)%GRIDPOINTSCOUNT)
    !nlink=nlink+network%BRS%BRANCH(k)%UPOINTSCOUNT
enddo

!construct branches

if (allocated(f1dimppar%grd_sre_fm)) then
    deallocate(f1dimppar%grd_sre_fm)
endif
allocate(f1dimppar%grd_sre_fm(ngrid)) 
grd_sre_fm => f1dimppar%grd_sre_fm

if (allocated(f1dimppar%grd_fmL_sre)) then
    deallocate(f1dimppar%grd_fmL_sre)
endif
allocate(f1dimppar%grd_fmL_sre(lnx1D,2)) 
grd_fmL_sre => f1dimppar%grd_fmL_sre

if (allocated(f1dimppar%branch)) then
    deallocate(f1dimppar%branch)
endif
allocate(f1dimppar%branch(4,nbran)) 
branch => f1dimppar%branch

if (allocated(f1dimppar%x)) then
    deallocate(f1dimppar%x)
endif
allocate(f1dimppar%x(ngrid))
x => f1dimppar%x

idx_i=1
do k=1,nbran
    !update index final
    idx_f=idx_i+network%BRS%BRANCH(k)%GRIDPOINTSCOUNT-1
    
    grd_sre_fm(idx_i:idx_f)=network%BRS%BRANCH(k)%GRD
    x(idx_i:idx_f)=network%BRS%BRANCH(k)%GRIDPOINTSCHAINAGES !chainage

    nl=network%BRS%BRANCH(k)%UPOINTSCOUNT
    do k2=1,nl
        L=network%BRS%BRANCH(k)%LIN(k2)
        grd_fmL_sre(L,:)=(/ idx_i+k2-1, idx_i+k2 /)
        
        !search for the GRD with <n1>? 
    	n1 = ln(1,L) 
        n2 = ln(2,L)	
        if (.not. ((grd_sre_fm(grd_fmL_sre(L,1)) .eq. n1) .or. (grd_sre_fm(grd_fmL_sre(L,1)) .eq. n2))) then
           write (msgbuf, '(a)') 'Links and nodes do not match.'
           call err_flush()
           iresult=1
        endif
    enddo
        
    !branch
    branch(1,k)=network%BRS%BRANCH(k)%NODEINDEX(1)
    branch(2,k)=network%BRS%BRANCH(k)%NODEINDEX(2)
    branch(3,k)=idx_i
    branch(4,k)=idx_f
    
    !update index initial
    idx_i=idx_f+1
enddo

!ngrid=network%numk !total number of mesh nodes (internal water level points)


maxlev=0 
do k=1,network%CSDEFINITIONS%COUNT 
    maxlev=max(maxlev,network%CSDEFINITIONS%CS(1)%LEVELSCOUNT)
end do

nnode=network%nds%count 
nhstat=nzbnd 
nqstat=nqbnd 
maxtab=ndx - ndxi !<we have as many tables as open boundaries
!if (comparereal(nzbnd+nqbnd,ndx-ndxi,1d-10)/=0) then !FM1DIMP2DO: why does the compiler complain when using <comparereal>?
if ((nzbnd+nqbnd).ne.(ndx-ndxi)) then
    write (msgbuf, '(a)') 'Number of open boundaries is different than number of water level + discharge boundaries'
    call err_flush()
    iresult=1    
endif
table_length=2 !length of each table. All have only 2 times.
ntabm=maxtab*table_length*2 !last 2 is for <time> and <values> 
nbrnod=network%NDS%MAXNUMBEROFCONNECTIONS

if (allocated(f1dimppar%bfrict)) then
    deallocate(f1dimppar%bfrict)
endif
allocate(f1dimppar%bfrict(3,nbran))
bfrict => f1dimppar%bfrict

do k=1,nbran
    
    !bfrict
    do k2=1,3 !< main channel, floodplain 1, floodplain 2
        select case (network%RGS%ROUGH(k2)%FRICTIONTYPE) !< where is the information per branch? add when several branches!
            case (0)
                bfrict(k2,k)=1
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
if (allocated(f1dimppar%bfricp)) then
    deallocate(f1dimppar%bfricp)
endif
allocate(f1dimppar%bfricp(6,ngrid)) !needs the part with FP1, FP2
bfricp => f1dimppar%bfricp

if (allocated(f1dimppar%nlev)) then
    deallocate(f1dimppar%nlev)
endif
allocate(f1dimppar%nlev(ngrid)) 

if (allocated(f1dimppar%hpack)) then
    deallocate(f1dimppar%hpack)
endif
allocate(f1dimppar%hpack(ngrid,3)) 

if (allocated(f1dimppar%qpack)) then
    deallocate(f1dimppar%qpack)
endif
allocate(f1dimppar%qpack(ngrid,3))

if (allocated(f1dimppar%waoft)) then
    deallocate(f1dimppar%waoft)
endif
allocate(f1dimppar%waoft(ngrid,18))
swaoft=size(f1dimppar%waoft,dim=2)

    !cross-sectional information (gridpoint,level)
if (allocated(f1dimppar%wft)) then
    deallocate(f1dimppar%wft)
endif
allocate(f1dimppar%wft(ngrid,maxlev)) 

if (allocated(f1dimppar%aft)) then
    deallocate(f1dimppar%aft)
endif
allocate(f1dimppar%aft(ngrid,maxlev)) 

if (allocated(f1dimppar%wtt)) then
    deallocate(f1dimppar%wtt)
endif
allocate(f1dimppar%wtt(ngrid,maxlev)) 

if (allocated(f1dimppar%att)) then
    deallocate(f1dimppar%att)
endif
allocate(f1dimppar%att(ngrid,maxlev)) 

if (allocated(f1dimppar%of)) then
    deallocate(f1dimppar%of)
endif
allocate(f1dimppar%of(ngrid,maxlev)) 

if (allocated(f1dimppar%hlev)) then
    deallocate(f1dimppar%hlev)
endif
allocate(f1dimppar%hlev(ngrid,maxlev))

!call fm1dimp_update_network(iresult) !update of the flow variables (change every time step)

do k=1,ngrid
    
    idx_fm=grd_sre_fm(k) !index of the global grid point in fm for the global gridpoint <k> in SRE
    
    !idx_crs
    !index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node, but 
    !due to precision, <network%ADM%gpnt2cross%F> may not be exactly 1. 
    if (network%ADM%gpnt2cross(idx_fm)%F>0.5) then
        if (network%ADM%gpnt2cross(idx_fm)%F<1-1.0e-6) then
           write (msgbuf, '(a)') 'It seems there is not 1 cross-section per mesh node.'
           call err_flush()
           iresult=1
        endif
        idx_crs=network%ADM%gpnt2cross(idx_fm)%C2     
    endif
    if (network%ADM%gpnt2cross(idx_fm)%F<0.5) then
        if (network%ADM%gpnt2cross(idx_fm)%F>1.0e-6) then
           write (msgbuf, '(a)') 'It seems there is not 1 cross-section per mesh node.'
           call err_flush()
           iresult=1
        endif
        idx_crs=network%ADM%gpnt2cross(idx_fm)%C1    
    endif
    
    !bfrictp
    bfricp(1,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(2,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    !deal properly with the values below when friction per section varies. 
    bfricp(3,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(4,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    bfricp(5,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(6,k)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
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
    !f1dimppar%x(k)=network%CRS%CROSS(idx_crs)%CHAINAGE !we cannot take the chaninage from the cross-section because when there are several branches, the first cross-section may come from another branch than the one we are dealing.
    
    !done in <fm1dimp_update_network>, which is called before the time step
    !
    !!nlev
    !f1dimppar%nlev(k)=network%CRS%CROSS(idx_crs)%TABDEF%LEVELSCOUNT
    !do k2=1,f1dimppar%nlev(k)
    !    f1dimppar%wft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWWIDTH(k2) 
    !    f1dimppar%aft(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%FLOWAREA(k2)  
    !    f1dimppar%wtt(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALWIDTH(k2) 
    !    !FM1DIMP2DO: deal with (at least error) case of rectangular cross-section with only two elevation points. 
    !    !In this case, the area for the low point is 0. This causes a huge interpolation error in SRE. 
    !    f1dimppar%att(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%TOTALAREA(k2)    
    !    f1dimppar%of(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%WETPERIMETER(k2)     
    !    f1dimppar%hlev(k,k2)=network%CRS%CROSS(idx_crs)%TABDEF%HEIGHT(k2)
    !end do !k2

    !dependent variables and waoft dealt with time step initialization.
    
    !!dependent variables
    !!FM1DIMP2DO: not sure this is needed here. Done in <SOFLOW_wrap>? -> Better here and keep <hpack> and <qpack> as SRE computing variables
    !!
    !do k2=1,3 !< time step before, intermediate, after
    !    !I don't think <k> below is correct. It should be an index mapping gridpoint and internal gridpoint
    !    f1dimppar%hpack(k,k2)=s0(k)
    !    !FM1DIMP2DO: for the given water level we have to compute the flow area and multiply by the velocity <u1>
    !    f1dimppar%qpack(k,k2)=100 
    !end do !k2
    
    !!waoft
    !!I don't think <k> below is correct. It should be an index mapping gridpoint and internal gridpoint
    !!FM1DIMP2DO: needs to be interpolated from link to cell centre
    !!FM1DIMP2DO: needs to be separated between flow and total
    !!FM1DIMP2DO: wetted perimeter get from results
    !!check right order in <FLNORM> and not in documentation. 
    !f1dimppar%waoft(k,1)=real(wu(k)) !wf = actual flow width 
    !f1dimppar%waoft(k,2)=real(wu(k)) !wt = actual total width
    !f1dimppar%waoft(k,3)=real(au(k)) !af = actual flow area
    !f1dimppar%waoft(k,4)=real(au(k)) !at = actual total area n
    !f1dimppar%waoft(k,5)=real(au(k)) !at = actual total area n+1
    !f1dimppar%waoft(k,6)=real(au(k)/wu(k)) !o = actual wetted perimeter
    !do k2=7,swaoft
    !    f1dimppar%waoft(k,k2)=0
    !enddo

end do !k
    
! 
!boundary conditions
!

if (allocated(f1dimppar%ntab)) then
    deallocate(f1dimppar%ntab)
endif     
allocate(f1dimppar%ntab(4,maxtab)) 
ntab => f1dimppar%ntab

table_number=0 !counter position in which the BC is saved in the table

!   h
if (allocated(f1dimppar%hbdpar)) then
    deallocate(f1dimppar%hbdpar)
endif
allocate(f1dimppar%hbdpar(3,nhstat)) 
hbdpar       => f1dimppar%hbdpar

do k=1,nhstat
    table_number=table_number+1
    
    hbdpar(1,k)=kbndz(2,k) !< first s1 point on the inside of the domain
    hbdpar(2,k)=1
    hbdpar(3,k)=table_number
    
    !FM1DIMP2DO: make this a subroutine? called in every loop for every BC
    ntab(1,table_number)=table_length !length of table 
    ntab(2,table_number)=table_number*table_length-1 !start address X
    ntab(3,table_number)=maxtab*table_length+table_number*table_length-1 !start address Y
    ntab(4,table_number)=0 !access method (0=continuous interpolation)
end do
!   q
if (allocated(f1dimppar%qbdpar)) then
    deallocate(f1dimppar%qbdpar)
endif
allocate(f1dimppar%qbdpar(3,nqstat)) 
qbdpar       => f1dimppar%qbdpar

do k=1,nqstat
    table_number=table_number+1
    
    qbdpar(1,k)=kbndu(2,k) !< first s1 point on the inside of the domain
    qbdpar(2,k)=1
    qbdpar(3,k)=table_number !< table number after the ones of <hbdpar>
    
    ntab(1,table_number)=table_length !length of table 
    ntab(2,table_number)=table_number*table_length-1 !start address X
    ntab(3,table_number)=maxtab*table_length+table_number*table_length-1 !start address Y
    ntab(4,table_number)=0 !access method (0=continuous interpolation)
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
if (allocated(f1dimppar%table)) then
    deallocate(f1dimppar%table)
endif
allocate(f1dimppar%table(ntabm)) 
table => f1dimppar%table

!nodes
if (allocated(f1dimppar%node)) then
    deallocate(f1dimppar%node)
endif 
allocate(f1dimppar%node(4,nnode))
node => f1dimppar%node

if (allocated(f1dimppar%numnod)) then
    deallocate(f1dimppar%numnod)
endif 
allocate(f1dimppar%numnod(nnode))
numnod => f1dimppar%numnod

!FM1DIMP2DO: REMOVE
!!everything should be in the loop, but there are some issues...
!!upstream
!k=1
!f1dimppar%node(1,k)=3
!f1dimppar%node(3,k)=1 !station number of the ones that are Q-stations. 
!f1dimppar%node(4,k)=1 
!!downstream
!k=2
!f1dimppar%node(1,k)=2
!f1dimppar%node(3,k)=1 !station number of the ones that are H-stations.
!f1dimppar%node(4,k)=1

do k=1,nnode
    !f1dimppar%node(1,k)= 
    !for some reason at this stage <network%NDS%NODE%NODETYPE> is only either 0 or 1
    !hence, if type 0 it means BC and I have to search which one is it
    !network%NDS%NODE%NODETYPE
    !! - -2    boundary node
    !! - -1    not set
    !! -  0    node with one reach connected
    !! -  1    connection node with more than one reach connected
    !! -  2    water level boundary
    !! -  3    Discharge boundary 
    !! -  4    Discharge boundary as tabulated function of water level
    !! -  5    Embedded node
    if (network%NDS%NODE(k)%NODETYPE .eq. 0) then !BC noce
        
        do k2=1,nhstat !search in hbdpar
            if (hbdpar(1,k2) .eq. network%NDS%NODE(k)%GRIDNUMBER) then
                node(1,k)=2 !H boundary
                node(2,k)=hbdpar(1,k2) !gridpoint
                node(3,k)=k2
                node(4,k)=1 !not sure what should it be
            endif
        enddo !nhstat
        
        if (node(1,k) .eq. 0) then !it is not hbdpar, we search in qbdpar
        
            do k2=1,nqstat !search in hbdpar
                if (qbdpar(1,k2) .eq. network%NDS%NODE(k)%GRIDNUMBER) then
                    node(1,k)=3 !Q boundary
                    node(2,k)=qbdpar(1,k2) !gridpoint
                    node(3,k)=k2
                    node(4,k)=1 !not sure what should it be
                endif
            enddo
        
        endif
        
        if (node(1,k) .eq. 0) then !it is not hbdpar nor qbdpar => error
            write (msgbuf, '(a)') 'There is a node which is neither internal nor H nor Q boundary.'
            call err_flush()
            iresult=1
        endif
        
    elseif (network%NDS%NODE(k)%NODETYPE .eq. 1) then !internal node
        node(1,k)=1
        !node(2:end) is undefined if internal node
    else
        write (msgbuf, '(a)') 'The type of node is not what I expected.'
        call err_flush()
        iresult=1
    endif
    
    numnod(k)=network%NDS%NODE(k)%NUMBEROFCONNECTIONS
    
end do

!nodes

if (allocated(f1dimppar%nodnod)) then
    deallocate(f1dimppar%nodnod)
endif 
allocate(f1dimppar%nodnod(nnode,nbrnod+1))
nodnod => f1dimppar%nodnod

!<kcol> saves the index to write in the row of <nonnod_aux>
if (allocated(kcol)) then
    deallocate(kcol)
endif 
allocate(kcol(nnode))
kcol=2 !first index is filled with its own node

!nodnod_ncol=20 !preallocating varibales
!!<nodnod_aux> has the values 
!if (allocated(nodnod_aux)) then
!    deallocate(nodnod_aux)
!endif 
!allocate(nodnod_aux(nnode,nodnod_ncol)) !we have to be sure it is less than <nodnod_ncol>
!nodnod_aux=0 !we have to be sure ther is no 0 node

!filling first index with its own node
do k=1,nnode
    nodnod(k,1)=k
enddo

do k=1,nbran
    
    idx_fr=network%BRS%BRANCH(k)%NODEINDEX(1)
    idx_to=network%BRS%BRANCH(k)%NODEINDEX(2)
    
    nodnod(idx_fr,kcol(idx_fr))=idx_to
    kcol(idx_fr)=kcol(idx_fr)+1
    
    nodnod(idx_to,kcol(idx_to))=idx_fr
    kcol(idx_to)=kcol(idx_to)+1
    
enddo

!FM1DIMP2DO: remove debug
f1dimppar%fm1dimp_debug_k1=1

!in steady computation <theta> in SRE is set to 1. We set it here to properly compute <u1> <q1>.
!maybe better to make a case
!do k=1,lnx
!    teta(k)=1
!enddo

!we use the pure1d morpho implementation, only data on x!
!I am not sure that implementation is correct though. 
if (jased > 0 .and. stm_included) then !passing if no morphpdynamics
    stmpar%morpar%mornum%pure1d=1


!the most downstream link points inside and we have to consider this for the flux.
call init_1dinfo() !<initialize_flow1d_implicit> is called before <init_1dinfo>. We have to call to call it here and it will not be called again because it will be allocated. 
!FM1DIMP2DO: I don't know why it fails compiling in case I ask to allocate here. It does not have the allocatable attribute, but neither it does <link1sign> 
!if (allocated(link1sign2)) then
!    deallocate(link1sign2)
!endif 
allocate(link1sign2(lnx)) 
do k=1,lnx1d !internal links
    link1sign2(k)=1
enddo
do k=lnxi+1,lnx1Db !boundary links
    !FM1DIMP2DO: we could create a variable with this mapping to prevent computation of max, min and x(nint) every timestep
    n1 = ln(1,k) 
    n2 = ln(2,k)
    nint=min(n1,n2) !from the two cells that this link connects, the minimum is internal, and hence we have data
    nout=max(n1,n2) !from the two cells that this link connects, the maximum is extrernal, and it is the one in which we have to set the water level
    if (f1dimppar%x(nint).eq.0) then !upstream
        link1sign2(k)=1 
    else ! downstream
        link1sign2(k)=-1
    endif
enddo

endif

!because the <height> is used in the cross-sections of SRE, <shift> cannot be used in cross-section

end subroutine initialize_flow1d_implicit

