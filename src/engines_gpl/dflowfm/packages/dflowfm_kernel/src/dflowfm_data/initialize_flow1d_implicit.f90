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
use m_alloc
use m_physcoef
use m_flowgeom, only: ndx, ndxi, wu, teta, lnx, lnx1D, lnx1Db, ln, lnxi, nd, kcs, tnode, wcl, dx, kcu
use unstruc_channel_flow, only: network
use m_flowexternalforcings !FM1DIMP2DO: do I need it?
use unstruc_messages
use m_flow, only: s0, s1, u1, au, hu, u_to_umain, frcu_mor, frcu, ifrcutp, ustb, qa !<ucmag> is velocity at cell centres, but we initialize <u1>
use m_sediment, only: stmpar, jased, stm_included
use m_fm_erosed, only: link1sign2, ndx_mor, ucyq_mor, hs_mor, ucxq_mor, kfsed, nd_mor, uuu, vvv, umod, zumod !ucx_mor, ucy_mor, 
use m_oned_functions, only: gridpoint2cross
use m_waves, only: taubxu

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
integer, dimension(:)                    , pointer :: grd_fm_sre
integer, dimension(:)                    , pointer :: idx_cs !FM1DIMP2DO: not a good name. Rename to <grd_sre_cs>
integer, dimension(:)                    , pointer :: lin
integer, dimension(:)                    , pointer :: grd
integer, dimension(:)                    , pointer :: kcs_sre
      
integer, dimension(:,:)                  , pointer :: grd_fmL_sre
integer, dimension(:,:)                  , pointer :: grd_fmLb_sre

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

type(tnode)    , allocatable :: nd_o(:) !type defined in <m_f1dimp_data>
!type(tnode_sre), pointer     :: nd_mor(:) !type defined in <m_f1dimp_data>

!debug
integer, pointer :: fm1dimp_debug_k1

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: kbr, knod, k1, k2, kbe, klnx, kndx, ksre, kn, kl, kd !FM1DIMP2DO: make the variables names consistent
integer :: c_lnx, c_ndx !counters
integer :: idx_crs, idx_sre, idx_fm !indices
integer :: n1, n2, nint, nout, pointscount, jpos
integer :: table_number
integer :: idx_fr, idx_to
integer :: idx_i, idx_f, nl, L, idx_fm_r, idx_fm_l, idx_l1, idx_l2, idx_sre_p, idx_sre_c
integer :: j
integer :: lnx_mor 

!move to function
integer :: idx_aux
integer :: min_1, min_2

!integer :: nlink !I don't think I need it global

integer, allocatable, dimension(:)   :: kcol
integer, allocatable, dimension(:)   :: grd_ghost_link_closest
integer, allocatable, dimension(:)   :: node_fm_processed
integer, allocatable, dimension(:)   :: grd_fmmv_fmsv !from FM multi-valued to FM single-valued
integer, allocatable, dimension(:,:) :: ln_fm

real :: swaoft

double precision :: wu_int, au_int

double precision, allocatable, dimension(:) :: frcu_mor_fm
double precision, allocatable, dimension(:) :: ifrcutp_fm

double precision, allocatable, dimension(:,:) :: wcl_fm

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
idx_cs       => f1dimppar%idx_cs

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
do kbr=1,nbran
    ngrid=ngrid+network%BRS%BRANCH(kbr)%GRIDPOINTSCOUNT
    ngridm=max(ngridm,network%BRS%BRANCH(kbr)%GRIDPOINTSCOUNT)
    !nlink=nlink+network%BRS%BRANCH(k)%UPOINTSCOUNT
enddo

!construct branches

if (allocated(f1dimppar%grd_sre_fm)) then
    deallocate(f1dimppar%grd_sre_fm)
endif
allocate(f1dimppar%grd_sre_fm(ngrid)) 
grd_sre_fm => f1dimppar%grd_sre_fm

if (allocated(f1dimppar%grd_fm_sre)) then
    deallocate(f1dimppar%grd_fm_sre)
endif
allocate(f1dimppar%grd_fm_sre(ndx+network%NDS%COUNT)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes.
grd_fm_sre => f1dimppar%grd_fm_sre

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

if (allocated(f1dimppar%idx_cs)) then
    deallocate(f1dimppar%idx_cs)
endif
allocate(f1dimppar%idx_cs(ngrid))
idx_cs => f1dimppar%idx_cs

if (allocated(f1dimppar%hpack)) then
    deallocate(f1dimppar%hpack)
endif
allocate(f1dimppar%hpack(ngrid,3)) 
hpack => f1dimppar%hpack

if (allocated(f1dimppar%qpack)) then
    deallocate(f1dimppar%qpack)
endif
allocate(f1dimppar%qpack(ngrid,3))
qpack => f1dimppar%qpack

if (allocated(f1dimppar%waoft)) then
    deallocate(f1dimppar%waoft)
endif
allocate(f1dimppar%waoft(ngrid,18))
waoft => f1dimppar%waoft
swaoft=size(f1dimppar%waoft,dim=2)

if (allocated(nd_mor)) then
    deallocate(nd_mor)
endif
allocate(nd_mor(ngrid))
!nd_mor => nd_mor_sre

if (allocated(f1dimppar%kcs_sre)) then
    deallocate(f1dimppar%kcs_sre)
endif
allocate(f1dimppar%kcs_sre(ngrid))
kcs_sre => f1dimppar%kcs_sre 
kcs_sre=1
    
!FM1DIMP2DO: add check on allocation
allocate(node_fm_processed(ndx+network%NDS%COUNT)) !more than we need
allocate(grd_fmmv_fmsv(ndx+network%NDS%COUNT)) !more than we need

!allocate(nd_mor(ndx+network%NDS%COUNT)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes.
!we cannot make a pointer to it because it has the same variable name
!do k=1,ndx
!    nd_mor(k)%lnx=nd(k)%lnx
!    nd_mor(k)%ln=nd(k)%ln
!enddo
nd_o=nd

if (allocated(grd_ghost_link_closest)) then
    deallocate(grd_ghost_link_closest)
endif
allocate(grd_ghost_link_closest(lnx+network%NDS%maxnumberofconnections*network%NDS%COUNT)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes. 

do kl=1,lnx
    grd_ghost_link_closest(kl)=kl
enddo

!if (allocated(node_processed)) then
!    deallocate(node_processed)
!endif
!allocate(node_processed(ndxi))
!node_processed=0


!
!BEGIN (LOB)
!
!Loop On Branches

idx_i=1
idx_sre=0
c_lnx=lnx
c_ndx=ndx
do kbr=1,nbran
    !update index final
    idx_f=idx_i+network%BRS%BRANCH(kbr)%GRIDPOINTSCOUNT-1
    
    grd_sre_fm(idx_i:idx_f)=network%BRS%BRANCH(kbr)%GRD
    x(idx_i:idx_f)=network%BRS%BRANCH(kbr)%GRIDPOINTSCHAINAGES !chainage

    nl=network%BRS%BRANCH(kbr)%UPOINTSCOUNT !only internal
    do kl=1,nl
        L=network%BRS%BRANCH(kbr)%LIN(kl)
        grd_fmL_sre(L,:)=(/ idx_i+kl-1, idx_i+kl /)
        
        !search for the GRD with <n1>? 
    	n1 = ln(1,L) 
        n2 = ln(2,L)	
        if (.not. ((grd_sre_fm(grd_fmL_sre(L,1)) .eq. n1) .or. (grd_sre_fm(grd_fmL_sre(L,1)) .eq. n2))) then
           write (msgbuf, '(a)') 'Links and nodes do not match.'
           call err_flush()
           iresult=1
        endif
        if (.not. ((grd_sre_fm(grd_fmL_sre(L,2)) .eq. n1) .or. (grd_sre_fm(grd_fmL_sre(L,2)) .eq. n2))) then
           write (msgbuf, '(a)') 'Links and nodes do not match.'
           call err_flush()
           iresult=1
        endif
    enddo !kl
    
    pointscount=network%BRS%BRANCH(kbr)%GRIDPOINTSCOUNT !FM1DIMP2DO: also make pointer?
    lin      => network%brs%branch(kbr)%lin
    grd      => network%brs%branch(kbr)%grd
    
    do kn=1,pointscount
        idx_sre=idx_sre+1
        idx_fm=grd(kn) 

        grd_fm_sre(idx_fm)=idx_sre
                
        !cross-section
        !FM1DIMP2DO: This part of the code is part of <set_cross_sections_to_gridpoints>, could be modularized.
        !->start 01
        if (kn==1 .or. kn==pointscount) then
            
           ! search for correct location
           if (kn==1) then 
              L = lin(1)
           else
              L = lin(pointscount-1)
           endif
           do kl = 1,nd(idx_fm)%lnx
              if (L == iabs(nd(idx_fm)%ln(kl))) then
                 jpos = kl
              endif
           enddo !kl
           
           !add ghost link
           if (nd(idx_fm)%lnx>2) then !bifurcation
               if (node_fm_processed(idx_fm).eq.0) then !not yet processed. We only add the links once per bifurcation. 
                   node_fm_processed(idx_fm)=1 !set to processed
                   
                   !!save the closest link associated to the new link
                   !if (allocated(Lv))then
                   !    deallocate(Lv)
                   !endif
                   !allocate(Lv(nd(idx_fm)%lnx))
                   !Lv=nd(idx_fm)%ln 
                   
                   !set ghost link as the one connected to bifurcation node
                   do kl=1,nd(idx_fm)%lnx !loop on number links attached to the bifurcation mesh node
                       c_lnx=c_lnx+1 !update link number to ghost link
                       grd_ghost_link_closest(c_lnx)=abs(nd_o(idx_fm)%ln(kl))
                       nd(idx_fm)%ln(kl)=c_lnx !set new link number to bifurcation
                   enddo !kl
                   
                   !nd_mor(idx_sre)%lnx=2
                   !if (allocated(nd_mor(idx_sre)%ln)) then
                   !   deallocate(nd_mor(idx_sre)%ln)
                   !endif
                   !allocate(nd_mor(idx_sre)%ln(2))
                   !if (i==1) then 
                   !   nd_mor(idx_sre)%ln(1)=c_lnx
                   !   nd_mor(idx_sre)%ln(2)=nd(idx_fm)%ln(jpos)
                   !else
                   !   nd_mor(idx_sre)%ln(1)=nd(idx_fm)%ln(jpos)
                   !   nd_mor(idx_sre)%ln(2)=-c_lnx
                   !endif
               
               endif !(node_fm_processed(idx_f).eq.0)
               
               
               if (jpos.ne.1) then !if it is equal to 1, it has been stored in the non-additional locations (`grd_fm_sre(idx_fm)=idx_sre`)
                  c_ndx=c_ndx+1
                  grd_fm_sre(c_ndx)=idx_sre
                  grd_fmmv_fmsv(c_ndx)=idx_fm
               endif !(jpos.ne.1)
           endif !(nd(idx_fm)%lnx>2)
           
        else
           jpos = 1
           
           !fill <nd>
           nd_mor(idx_sre)%lnx=nd(idx_fm)%lnx
           nd_mor(idx_sre)%ln=nd(idx_fm)%ln
           
           grd_fmmv_fmsv(idx_fm)=idx_fm
        endif  
        !-> end 01
        
        idx_cs(idx_sre)=gridpoint2cross(idx_fm)%cross(jpos) !cross-section index associated to the FM gridpoint per branch
        !if there is not a unique cross-section per gridpoint per branch, <ic=-999>. It is not needed to check
        !this here because it is already checked in <flow_sedmorinit>, which is called before <initialize_flow1d_implicit>
        
        !icd=network%crs%cross(ic) !cross-section associated to the FM gridpoint per branch
        
    enddo !kn

    !deal with values at begin and end of the branch
    !copying the value after and before, respectively
    !as there is information at the link, which is closer to the
    !end and beginning of the SRE node we are filling than 
    !the previous (or later) SRE node, it would be more accurate
    !to fill using the link info rather than the SRE info. 
    do kbe=1,2 !upstream and downstram
        if (kbe.eq.1) then !begin of branch
            idx_sre_p=idx_sre-pointscount+1 !paste
            idx_sre_c=idx_sre-pointscount+2 !copy
        else !end of branch
            idx_sre_p=idx_sre !paste
            idx_sre_c=idx_sre-1 !copy
        endif 
    
        do k2=1,3 !< time step in SRE [before, intermediate, after]
            !discharge
            qpack(idx_sre_p,k2)=qpack(idx_sre_c,k2) 
        enddo
        
        !waoft
        do k2=1,swaoft
            waoft(idx_sre_p,k2)=waoft(idx_sre_c,k2)
        enddo
    enddo !kbe
    
    !branch    
    branch(1,kbr)=network%BRS%BRANCH(kbr)%NODEINDEX(1)
    branch(2,kbr)=network%BRS%BRANCH(kbr)%NODEINDEX(2)
    branch(3,kbr)=idx_i
    branch(4,kbr)=idx_f
    
    !update index initial
    idx_i=idx_f+1
enddo !branch
lnx_mor=c_lnx !store new number of links (considering ghost links)
ndx_mor=c_ndx !store new number of flow nodes (considering multivaluedness)

!
!END (LOB)
!

!BEGIN (FAAL)
!
!Fill Arrays that need Additional Link
!
!FM1DIMP2DO: If friction varies with time, <frcu_mor> is updated. The subroutine that
!does that must be modified to also adapt the friction in the ghost links.

!frcu_mor_fm=frcu_mor !copy to temporary array
!if (allocated(frcu_mor)) then
!    deallocate(frcu_mor)
!endif
!allocate(frcu_mor(lnx_mor)) 
!frcu_mor_fm=frcu_mor

!ifrcutp_fm=ifrcutp !copy to temporary array
!if (allocated(ifrcutp)) then
!    deallocate(ifrcutp)
!endif
!allocate(ifrcutp(lnx_mor)) 

    !allocate
!call realloc(frcu_mor,lnx_mor)
!call realloc(ifrcutp,lnx_mor)

!FM1DIMP2DO: how do I <realloc> with more than 1 dimension?
!call realloc(wcl,lnx_mor) 
!wcl_fm=wcl !copy to temporary array
!if (allocated(wcl)) then
!    deallocate(wcl)
!endif
!allocate(wcl(2,lnx_mor)) 
!
!    !copy data from links existing in FM
!do klnx=1,lnx
!!    frcu_mor(klnx)=frcu_mor_fm(klnx)
!!    ifrcutp(klnx)=ifrcutp_fm(klnx)
!    
!    do L=1,2
!        wcl(L,klnx)=wcl_fm(L,klnx)
!    enddo !L
!enddo !klnx

    !add data of closest link to new links
!FM1DIMP2DO: we should copy data from the SRE flow node associated to the ghost link, rather than the closest link.
!FM1DIMP2DO: Are we sure that the width and other parameters are updated correctly with time?
!do klnx=lnx+1,lnx_mor
!    frcu_mor(klnx)=frcu_mor_fm(grd_ghost_link_closest(klnx))
!    ifrcutp(klnx)=ifrcutp_fm(grd_ghost_link_closest(klnx))
!    
!    do L=1,2
!        wcl(L,klnx)=1-wcl(L,grd_ghost_link_closest(klnx)) !it should be equal to 0.5
!    enddo !L
!enddo !klnx
 
!
!END (FAAL)
!

!
!BEGIN (FFMA)
!
!Fill FM Arrays

!FM1DIMP2DO: Make asubroutine for this in which the input is the variable (<s0>), the map (<grd_fmmv_fmsv>), and the indices (<ndx>, <ndx_mor>)

!nodes
    !allocate
call reallocate_fill(s0,grd_fmmv_fmsv,ndx,ndx_mor) 
call reallocate_fill(s1,grd_fmmv_fmsv,ndx,ndx_mor) 

!FM1DIMP2DO: how do I <realloc> with more than 1 dimension?


!call realloc(s0,ndx_mor)
!call realloc(s1,ndx_mor)
!
!    !fill
!do kndx=ndx+1,ndx_mor
!
!enddo

!links
    !allocate
call reallocate_fill(u1,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(au,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(hu,grd_ghost_link_closest,lnx,lnx_mor)
call reallocate_fill(dx,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(wu,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(frcu_mor,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(frcu,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(qa,grd_ghost_link_closest,lnx,lnx_mor) 

call reallocate_fill_int(ifrcutp,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill_int(kcu,grd_ghost_link_closest,lnx,lnx_mor)

!call realloc(u1,lnx_mor)
!call realloc(au,lnx_mor)
!call realloc(wu,lnx_mor)
!call realloc(frcu_mor,lnx_mor)
!call realloc(ifrcutp,lnx_mor)
!
!FM1DIMP2DO: how do I <realloc> with more than 1 dimension?
!call realloc(wcl,lnx_mor) 
wcl_fm=wcl !copy to temporary array
if (allocated(wcl)) then
    deallocate(wcl)
endif
allocate(wcl(2,lnx_mor)) 

ln_fm=ln
if (allocated(ln)) then
    deallocate(ln)
endif
allocate(ln(2,lnx_mor))

    !copy data from links existing in FM
do kl=1,lnx
!    frcu_mor(klnx)=frcu_mor_fm(klnx)
!    ifrcutp(klnx)=ifrcutp_fm(klnx)
    
    do kd=1,2
        wcl(kd,kl)=wcl_fm(kd,kl)
        ln(kd,kl)=ln_fm(kd,kl)
    enddo !kl
enddo !klnx

    !fill
do kl=lnx+1,lnx_mor
    !u1(kl)=u1(grd_ghost_link_closest(kl))
    !au(kl)=au(grd_ghost_link_closest(kl))
    !wu(kl)=wu(grd_ghost_link_closest(kl))
    !frcu_mor(kl)=frcu_mor(grd_ghost_link_closest(kl))
    !ifrcutp(kl)=ifrcutp(grd_ghost_link_closest(kl))
    
    do kd=1,2
        wcl(kd,kl)=1-wcl(kd,grd_ghost_link_closest(kl)) !it should be equal to 0.5
        ln(kd,kl)=ln(kd,grd_ghost_link_closest(kl))
    enddo !L
enddo

!
!END (FFMA)
!

!BEGIN (FIC)
!
!Fill Initial Condition
!
!Data must be available already at ghost links and nodes
do ksre=1,ngrid
    idx_sre=ksre 
    idx_fm=grd_sre_fm(idx_sre)
    
    !links connected to a given fm grid node
    idx_l1=abs(nd(idx_fm)%ln(1))
    idx_l2=abs(nd(idx_fm)%ln(2))
        
        !initial condition
        do k2=1,3 !< time step in SRE [before, intermediate, after]
            !water level
            hpack(idx_sre,k2)=s1(idx_fm)
            !discharge
            qpack(idx_sre,k2)=0.5*(au(idx_l1)*u1(idx_l1)+au(idx_l2)*u1(idx_l2))
        enddo !k2
        
        !waoft
        wu_int=0.5*(wu(idx_l1)+wu(idx_l2))
        au_int=0.5*(au(idx_l1)+au(idx_l2))
        
        !FM1DIMP2DO: needs to be separated between flow and total
        !check right order in <FLNORM> and not in documentation. 

        waoft(idx_sre,1)=real(wu_int) !wf = actual flow width 
        waoft(idx_sre,2)=real(wu_int) !wt = actual total width
        waoft(idx_sre,3)=real(au_int) !af = actual flow area
        waoft(idx_sre,4)=real(au_int) !at = actual total area n
        waoft(idx_sre,5)=real(au_int) !at = actual total area n+1
        waoft(idx_sre,6)=real(au_int/wu_int) !o = actual wetted perimeter
        do k2=7,swaoft
            waoft(idx_sre,k2)=0
        enddo !k2
enddo !ksre

!END (FIC)       
        
if (allocated(f1dimppar%grd_fmLb_sre)) then
    deallocate(f1dimppar%grd_fmLb_sre)
endif
allocate(f1dimppar%grd_fmLb_sre(lnx1Db-lnxi,2))
grd_fmLb_sre => f1dimppar%grd_fmLb_sre

idx_fm=0
do L=lnxi+1,lnx1Db !boundary links
    idx_fm=idx_fm+1
    n1 = ln(1,L) 
    n2 = ln(2,L)
    nint=min(n1,n2) !from the two cells that this link connects, the minimum is internal, and hence we have data
    nout=max(n1,n2) !from the two cells that this link connects, the maximum is extrernal, and it is the one in which we have to set the water level
    
    !FM1DIMP2DO: move to function or search for smarter way
    !grd_fmLb_sre(k,1)=findloc(grd_sre_fm,nint) !sre index with <nint> FM value !not working fine due to type of array I guess. 
    idx_aux=1
    min_1=abs(grd_sre_fm(1)-nint)
    do k2=2,size(grd_sre_fm)
        min_2=abs(grd_sre_fm(k2)-nint)
        if (min_2 < min_1) then
            min_1=min_2
            idx_aux=k2
        endif
    enddo
    
    grd_fmLb_sre(idx_fm,1)=idx_aux !SRE index of the boundary cell
    grd_fmLb_sre(idx_fm,2)=nout !FM index of the ghost cell centre associated to link <L>
    
    !mask grid
    kcs_sre(idx_aux)=-1 !FM1DIMP2DO: I am not sure I need this or I better deal with directions in <fm_erosed> and here just set to 1 but the right dimensions.
enddo

!ngrid=network%numk !total number of mesh nodes (internal water level points)


maxlev=0 
do k1=1,network%CSDEFINITIONS%COUNT 
    maxlev=max(maxlev,network%CSDEFINITIONS%CS(1)%LEVELSCOUNT)
enddo

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

do kbr=1,nbran
    
    !bfrict
    do k2=1,3 !< main channel, floodplain 1, floodplain 2
        select case (network%RGS%ROUGH(k2)%FRICTIONTYPE) !< where is the information per branch? add when several branches!
            case (0)
                bfrict(k2,kbr)=1
            case default
                write (msgbuf, '(a)') 'Only constant Chezy friction is supported at the moment.'
                call err_flush()
                iresult=1
        end select
    enddo
!
!bfrict(1,i)=cfrchc (1) : Ch�zy constant
!bfrict(1,i)=cfrchq (2) : Ch�zy function of discharge
!bfrict(1,i)=cfrchh (3) : Ch�zy function of water level
!bfrict(1,i)=cfrman (4) : Manning constant
!bfrict(1,i)=cfrskn (5) : Strickler 1 constant ( k n )
!bfrict(1,i)=cfrsks (6) : Strickler 2 constant ( k s )
!bfrict(1,i)=cfrnik (7) : Nikuradze constant
!bfrict(1,i)=cfreng (8) : Engelund predicto

enddo !kbr
 

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

if (allocated(f1dimppar%bedlevel)) then
    deallocate(f1dimppar%bedlevel)
endif
allocate(f1dimppar%bedlevel(ngrid))

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

do ksre=1,ngrid
    
    idx_fm=grd_sre_fm(ksre) !index of the global grid point in fm for the global gridpoint <k> in SRE
    
    !idx_crs
    !index of the cross-section at grid-node <k>. Should be the same as C2 as there is a cross-section per node, but 
    !due to precision, <network%ADM%gpnt2cross%F> may not be exactly 1. 
    !if (network%ADM%gpnt2cross(idx_fm)%F>0.5) then
    !    if (network%ADM%gpnt2cross(idx_fm)%F<1-1.0e-6) then
    !       write (msgbuf, '(a)') 'It seems there is not 1 cross-section per mesh node.'
    !       call err_flush()
    !       iresult=1
    !    endif
    !    idx_crs=network%ADM%gpnt2cross(idx_fm)%C2     
    !endif
    !if (network%ADM%gpnt2cross(idx_fm)%F<0.5) then
    !    if (network%ADM%gpnt2cross(idx_fm)%F>1.0e-6) then
    !       write (msgbuf, '(a)') 'It seems there is not 1 cross-section per mesh node.'
    !       call err_flush()
    !       iresult=1
    !    endif
    !    idx_crs=network%ADM%gpnt2cross(idx_fm)%C1    
    !endif
    
    
    !bfrictp
    idx_crs=idx_cs(ksre)
    
    bfricp(1,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(2,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    !deal properly with the values below when friction per section varies. 
    bfricp(3,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(4,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
    bfricp(5,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUEPOS(1)
    bfricp(6,ksre)=network%CRS%CROSS(idx_crs)%FRICTIONVALUENEG(1)
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
    
    !<hpack>, <qpack>, and <waoft> are saved as <fm1dimp> variables but need to be initialized here. 
    
    !dependent variables
    !
    !do k2=1,3 !< time step [before, intermediate, after]
    !    f1dimppar%hpack(k,k2)=s0(idx_fm)
    !    if (nd(idx_fm)%lnx>2) then
    !        f1dimppar%qpack(k,k2)=au(idx_fm)*u1(idx_fm) 
    !    else
    !        f1dimppar%qpack(k,k2)=au(idx_fm)*u1(idx_fm) 
    !    endif
    !end do !k2

enddo !ksre
    
! 
!boundary conditions
!

!<table> will contain 4 elements per BC:
!   -X1 (time)
!   -X2 (time)
!   -Y1 (variable, e.g., water level)
!   -Y2 (variable, e.g., water level)
!the order in <table> is:
!   -H-boundaries 
!   -Q-boundaires

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

do k1=1,nhstat
    table_number=table_number+1
    
    hbdpar(1,k1)=kbndz(2,k1) !< first s1 point on the inside of the domain
    hbdpar(2,k1)=1
    hbdpar(3,k1)=table_number
    
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

do k1=1,nqstat
    table_number=table_number+1
    
    qbdpar(1,k1)=kbndu(2,k1) !< first s1 point on the inside of the domain
    qbdpar(2,k1)=1
    qbdpar(3,k1)=table_number !< table number after the ones of <hbdpar>
    
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
node = -999 !we use this value to check that it has not been filled.

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

do knod=1,nnode
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
    if (network%NDS%NODE(knod)%NODETYPE .eq. 0) then !BC noce
        
        do k2=1,nhstat !search in hbdpar
            if (hbdpar(1,k2) .eq. network%NDS%NODE(knod)%GRIDNUMBER) then
                node(1,knod)=2 !H boundary
                node(2,knod)=hbdpar(1,k2) !gridpoint
                node(3,knod)=k2
                node(4,knod)=1 !not sure what should it be
            endif
        enddo !nhstat
        
        if (node(1,knod) .eq. -999) then !it is not hbdpar, we search in qbdpar
        
            do k2=1,nqstat !search in hbdpar
                if (qbdpar(1,k2) .eq. network%NDS%NODE(knod)%GRIDNUMBER) then
                    node(1,knod)=3 !Q boundary
                    node(2,knod)=qbdpar(1,k2) !gridpoint
                    node(3,knod)=k2
                    node(4,knod)=1 !not sure what should it be
                endif
            enddo
        
        endif
        
        if (node(1,knod) .eq. 0) then !it is not hbdpar nor qbdpar => error
            write (msgbuf, '(a)') 'There is a node which is neither internal nor H nor Q boundary.'
            call err_flush()
            iresult=1
        endif
        
    elseif (network%NDS%NODE(knod)%NODETYPE .eq. 1) then !internal node
        node(1,knod)=1
        !node(2:end) is undefined if internal node
    else
        write (msgbuf, '(a)') 'The type of node is not what I expected.'
        call err_flush()
        iresult=1
    endif
    
    numnod(knod)=network%NDS%NODE(knod)%NUMBEROFCONNECTIONS+1
    
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
do knod=1,nnode
    nodnod(knod,1)=knod
enddo

do kbr=1,nbran
    
    idx_fr=network%BRS%BRANCH(kbr)%NODEINDEX(1)
    idx_to=network%BRS%BRANCH(kbr)%NODEINDEX(2)
    
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
    
    !call realloc(s0,ndx_mor)
    !call realloc(s1,ndx_mor)
    !call realloc(ucxq_mor,ndx_mor)
    !call realloc(ucyq_mor,ndx_mor)
    !call realloc(hs_mor,ndx_mor)
    !call realloc(u_to_umain,ndx_mor)
    !call realloc(uuu,ndx_mor)
    !call realloc(vvv,ndx_mor)
    !call realloc(umod,ndx_mor)
    !call realloc(zumod,ndx_mor)
    !call realloc(kfsed,ndx_mor)
    
    !FM1DIMP2DO: not needed here. It will be updated every time step in <
    !do k=ndx+1,ndx_mor
    !    s1(k)=hpack(grd_fm_sre(k),3)
    !enddo
    
    
    !these variables cannot be allocated in <inipointers_erosed> because at that point we have not initialized 1dimp
    !if (allocated(ucxq_mor)) then
    !    deallocate(ucxq_mor,ucyq_mor,hs_mor) !ucx_mor,ucy_mor
    !endif
    !allocate(ucxq_mor(1:ndx_mor), ucyq_mor(1:ndx_mor), hs_mor(1:ndx_mor), stat=iresult)   !ucx_mor(1:ndx_mor), ucy_mor(1:ndx_mor)
    !ucxq_mor = 0d0; ucyq_mor = 0d0; hs_mor = 0d0 
    !!ucx_mor = 0d0; ucy_mor = 0d0
    !
    !if (allocated(u_to_umain)) then
    !    deallocate(u_to_umain)
    !endif
    !allocate(u_to_umain(1:ndx_mor))
    !
    !allocate(uuu(1:ndx_mor),vvv(1:ndx_mor),umod(1:ndx_mor),zumod(1:ndx_mor))
    !uuu=0.0_fp; vvv=0.0_fp; umod=0.0_fp; zumod=0.0_fp
    
    !if (allocated(kfsed)) then
    !    deallocate(kfsed)
    !endif
    !allocate(kfsed(1:ndx_mor))
    
    
!the most downstream link points inside and we have to consider this for the flux.
call init_1dinfo() !<initialize_flow1d_implicit> is called before <init_1dinfo>. We have to call to call it here and it will not be called again because it will be allocated. 
!FM1DIMP2DO: I don't know why it fails compiling in case I ask to allocate here. It does not have the allocatable attribute, but neither it does <link1sign> 
!if (allocated(link1sign2)) then
!    deallocate(link1sign2)
!endif 
allocate(link1sign2(lnx)) 
do kl=1,lnx1d !internal links
    link1sign2(kl)=1
enddo
do kl=lnxi+1,lnx1Db !boundary links
    !FM1DIMP2DO: we could create a variable with this mapping to prevent computation of max, min and x(nint) every timestep
    n1 = ln(1,kl) 
    n2 = ln(2,kl)
    nint=min(n1,n2) !from the two cells that this link connects, the minimum is internal, and hence we have data
    nout=max(n1,n2) !from the two cells that this link connects, the maximum is extrernal, and it is the one in which we have to set the water level
    if (f1dimppar%x(nint).eq.0) then !upstream
        link1sign2(kl)=1 
    else ! downstream
        link1sign2(kl)=-1
    endif
enddo

endif

!because the <height> is used in the cross-sections of SRE, <shift> cannot be used in cross-section

    contains

    !FM1DIMP2DO: Move to module
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END MAIN SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!
!BEGIN reallocate_fill
!
    
subroutine reallocate_fill(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

double precision, dimension(:), allocatable, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call realloc(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill
    
!
!END reallocate_fill_int
!

!
!BEGIN reallocate_fill_int
!

!FM1DIMP2DO: There must be a better way to do it. This is just the same but for integers
    
subroutine reallocate_fill_int(val,idx_mask,idxi,idxf)

use m_alloc

implicit none

integer, dimension(:), allocatable, intent(inout) :: val
integer, dimension(idxf), intent(in) :: idx_mask
integer, intent(in) :: idxi, idxf

integer :: k


call realloc(val,idxf)

do k=idxi+1,idxf
    val(k)=val(idx_mask(k))
enddo

end subroutine reallocate_fill_int
    
!
!END reallocate_fill
!

    end subroutine initialize_flow1d_implicit