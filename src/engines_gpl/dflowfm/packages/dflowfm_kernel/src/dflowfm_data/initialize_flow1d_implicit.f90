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
use m_flowgeom, only: ndx, ndxi, wu, teta, lnx, lnx1D, lnx1Db, ln, lnxi, nd, kcs, tnode, wcl, dx, kcu, acl, snu, csu, wu_mor, bai_mor, bl, griddim, dxi, wcx1, wcx2, wcy1, wcy2, ba
use unstruc_channel_flow, only: network
use m_flowexternalforcings !FM1DIMP2DO: do I need it?
use unstruc_messages
use m_flow, only: s0, s1, u1, au, hu, u_to_umain, frcu_mor, frcu, ifrcutp, ustb, qa, kmx, ndkx, ndkx_mor, z0urou
use m_sediment, only: stmpar, jased, stm_included, sedtra, vismol, kcsmor
use m_initsedtra, only: initsedtra
use m_fm_erosed, only: link1, link1sign, link1sign2, ndx_mor, lnx_mor, lnxi_mor, ndxi_mor, ucyq_mor, hs_mor, ucxq_mor, kfsed, nd_mor, uuu, vvv, umod, zumod, e_dzdn, e_sbcn, lsedtot, e_sbn, dbodsd, dzbdt, pmcrit, frac, ln_mor
use m_oned_functions, only: gridpoint2cross, t_gridp2cs
use m_waves, only: taubxu
use morphology_data_module, only: allocsedtra
use m_turbulence, only: rhowat
use m_xbeach_data, only: ktb
use m_bedform, only: bfmpar

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
integer, dimension(:)                    , pointer :: grd_sre_cs 
integer, dimension(:)                    , pointer :: grd_ghost_link_closest
integer, dimension(:)                    , pointer :: grd_fmmv_fmsv
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

type(tnode)    , allocatable :: nd_o(:) !Copy of <nd> for reworking <nd>
!type(tnode)    , pointer     :: nd_mor(:) !Modified <nd> for <bott3d>
type(t_gridp2cs), dimension(:), allocatable :: gridpoint2cross_o

!debug
integer, pointer :: fm1dimp_debug_k1

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: kbr, knod, k1, k2, kbe, klnx, ksre, kn, kl, kd, ksed !FM1DIMP2DO: make the variables names consistent
integer :: ndx_max, lnx_max
integer :: c_lnx, c_ndx !counters
integer :: idx_crs, idx_sre, idx_fm !indices
integer :: n1, n2, nint, nout, pointscount, jpos
integer :: table_number
integer :: idx_fr, idx_to
integer :: idx_i, idx_f, nl, L, L2, idx_fm_r, idx_fm_l, idx_l1, idx_l2, idx_sre_p, idx_sre_c, idx_n
integer :: j

integer, dimension(1) :: idx_findloc
!integer :: lnx_mor 

!move to function
integer :: idx_aux
integer :: min_1, min_2

!integer :: nlink !I don't think I need it global

integer, allocatable, dimension(:)   :: kcol
!integer, allocatable, dimension(:)   :: grd_ghost_link_closest
!integer, allocatable, dimension(:)   :: node_fm_processed
!integer, allocatable, dimension(:)   :: grd_fmmv_fmsv !from FM multi-valued to FM single-valued
!integer, allocatable, dimension(:,:) :: ln_o

real :: swaoft

double precision :: wu_int, au_int

!double precision, allocatable, dimension(:) :: frcu_mor_fm
!double precision, allocatable, dimension(:) :: ifrcutp_fm

double precision, allocatable, dimension(:,:) :: wcl_fm
!double precision, allocatable, dimension(:,:) :: e_sbcn_fm
!double precision, allocatable, dimension(:,:) :: e_sbn_fm
double precision, allocatable, dimension(:,:) :: bodsed_o
double precision, allocatable, dimension(:,:) :: frac_o

!!
!! POINT
!!

!pointer cannot be before the array is allocated, here only non-allocatable arrays

table_length           => f1dimppar%table_length
maxtab                 => f1dimppar%maxtab
nnode                  => f1dimppar%nnode
ntabm                  => f1dimppar%ntabm
nbran                  => f1dimppar%nbran
ngrid                  => f1dimppar%ngrid
nbrnod                 => f1dimppar%nbrnod
maxlev                 => f1dimppar%maxlev
ngridm                 => f1dimppar%ngridm
nhstat                 => f1dimppar%nhstat
nqstat                 => f1dimppar%nqstat
grd_sre_cs             => f1dimppar%grd_sre_cs
grd_ghost_link_closest => f1dimppar%grd_ghost_link_closest
grd_fmmv_fmsv          => f1dimppar%grd_fmmv_fmsv

!!
!! CALC
!!

f1dimp_initialized=.true.
iresult=0

!<flwpar>
f1dimppar%g=ag 

f1dimppar%rhow=rhomean 


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

ndx_max=ndx+network%NDS%COUNT !maximum number of multivalued flownodes
lnx_max=lnx+network%NDS%maxnumberofconnections*network%NDS%COUNT !maximum number of links considering added ghost links

!construct branches

if (allocated(f1dimppar%grd_sre_fm)) then
    deallocate(f1dimppar%grd_sre_fm)
endif
allocate(f1dimppar%grd_sre_fm(ngrid)) 
grd_sre_fm => f1dimppar%grd_sre_fm

if (allocated(f1dimppar%grd_fm_sre)) then
    deallocate(f1dimppar%grd_fm_sre)
endif
allocate(f1dimppar%grd_fm_sre(ndx_max)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes.
grd_fm_sre => f1dimppar%grd_fm_sre
grd_fm_sre=0

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

if (allocated(f1dimppar%grd_sre_cs)) then
    deallocate(f1dimppar%grd_sre_cs)
endif
allocate(f1dimppar%grd_sre_cs(ngrid))
grd_sre_cs => f1dimppar%grd_sre_cs

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
allocate(nd_mor(ndx_max)) !more than we need
do kd=1,ndx
    nd_mor(kd)=nd(kd)
enddo

if (allocated(f1dimppar%kcs_sre)) then
    deallocate(f1dimppar%kcs_sre)
endif
allocate(f1dimppar%kcs_sre(ngrid))
kcs_sre => f1dimppar%kcs_sre 
kcs_sre=1
    
if (allocated(f1dimppar%grd_fmmv_fmsv)) then
    deallocate(f1dimppar%grd_fmmv_fmsv)
endif
allocate(f1dimppar%grd_fmmv_fmsv(ndx_max)) !more than we need
grd_fmmv_fmsv => f1dimppar%grd_fmmv_fmsv
!allocate every node with itself
do kd=1,ndx_max
    grd_fmmv_fmsv(kd)=kd
enddo

!allocate(nd_mor(ndx_max)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes.
!we cannot make a pointer to it because it has the same variable name
!do k=1,ndx
!    nd_mor(k)%lnx=nd(k)%lnx
!    nd_mor(k)%ln=nd(k)%ln
!enddo
nd_o=nd

!ln_o=ln
if (allocated(ln_mor)) then
    deallocate(ln_mor)
endif
allocate(ln_mor(2,lnx_max))
do kl=1,lnx
    do kd=1,2
        !ln(kd,kl)=ln_o(kd,kl)
        ln_mor(kd,kl)=ln(kd,kl)
    enddo
enddo


if (allocated(f1dimppar%grd_ghost_link_closest)) then
    deallocate(f1dimppar%grd_ghost_link_closest)
endif
allocate(f1dimppar%grd_ghost_link_closest(lnx_max)) !we allocate more than we need. The maximum number of bifurcations and confluences is less than the number of nodes. 
grd_ghost_link_closest => f1dimppar%grd_ghost_link_closest
do kl=1,lnx
    grd_ghost_link_closest(kl)=kl
enddo

!FM1DIMP2DO: I am now adapting the input for using the morphodynamic implementation of Pure 1D. However, 
!I amnot sure it is the best. THis should be revisited with Bert :). 
stmpar%morpar%mornum%pure1d=1
call init_1dinfo() !<initialize_flow1d_implicit> is called before <init_1dinfo>. We have to call it here and it will not be called again because it will be allocated. 

allocate(link1sign(lnx_max))
link1sign=1

allocate(link1sign2(lnx_max))
link1sign2=0
!All internal links have direction 1
do kl=1,lnxi
    link1sign2(kl)=1
enddo
!do kl=lnxi+1,lnx
!    link1sign2(kl)=1
!enddo

!if (allocated(node_processed)) then
!    deallocate(node_processed)
!endif
!allocate(node_processed(ndxi))
!node_processed=0

!copy to <gridpoint2cross_o>
if (allocated(gridpoint2cross_o)) then
    deallocate(gridpoint2cross_o)
endif
allocate(gridpoint2cross_o(ndx_max))
do kd=1,ndxi
    gridpoint2cross_o(kd)=gridpoint2cross(kd) 
    !if (gridpoint2cross_o(kd)%num_cross_sections>1) then
    !     write (msgbuf, '(a)') 'There is more than 1 cross-section per node (this error should have been captured before).'
    !     call err_flush()
    !     iresult=1
    !endif
enddo
!allocate
if (allocated(gridpoint2cross)) then
    deallocate(gridpoint2cross)
endif
allocate(gridpoint2cross(ndx_max))
!internal cross-sections are the same as they were (1 CS per flownode).
do kd=1,ndxi
    gridpoint2cross(kd)=gridpoint2cross_o(kd) 
enddo
!at ghost-boundary flownodes we set the number of CS to 0 to prevent looping on them (there is no CS)
do kd=ndxi+1,ndx
    gridpoint2cross(kd)%num_cross_sections=0 !This prevents it is looped in <fm_update_crosssections>
enddo


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
        
        !FM1DIMP2DO: Do we need this?
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

        !cross-section
        if (kn==1 .or. kn==pointscount) then
           
           !FM1DIMP2DO: This part of the code is part of <set_cross_sections_to_gridpoints>, could be modularized.
           if (kn==1) then 
              L = lin(1)
           else
              L = lin(pointscount-1)
           endif
           do kl = 1,nd_o(idx_fm)%lnx
              if (L == iabs(nd_o(idx_fm)%ln(kl))) then
                 jpos = kl
              endif
           enddo !kl
               
           !add ghost link
           if (nd(idx_fm)%lnx>2) then !bifurcation
               
               !link
                c_lnx=c_lnx+1 !update link number to ghost link
                
                grd_ghost_link_closest(c_lnx)=abs(nd_o(idx_fm)%ln(jpos))
                
                !In <nd> we keep the junction node <idx_fm> connected to several branches via the ghost link, as <nd> is used for the nodal point relation.
                !nd(idx_fm)%ln(jpos)=c_lnx !set ghost link as the one connected to junction flownode
                
                !FM1DIMP2DO: The link sign is set to positive for all cases and the sign is given in <link1sign2>. This should be revisited
                !when revisiting <fm_upwbed>. -> NO, set sign in <nd>
                !Impact on nodal point relation. I kept the signs there. 
                
                !FM1DIMP2DO: It seems all links point toward a junction in the standard scheme (based on nodal point relation). Is it true?
                if (kn==1) then 
                   link1sign2(c_lnx)=1 !link direction for morphodynamics
                   nd(idx_fm)%ln(jpos)=-c_lnx !set ghost link as the one connected to junction flownode
                else
                   link1sign2(c_lnx)=1
                   nd(idx_fm)%ln(jpos)=c_lnx !set ghost link as the one connected to junction flownode
                endif

                !node
                c_ndx=c_ndx+1 !update node number to multivalued node

                !save the flownode closest to the junction node in the branch under consideration
                !and
                !change the flownode connected to the first link connected to the junction flownode along the branch under consideration to the new flownode
                n1=ln(1,grd_ghost_link_closest(c_lnx)) !flownode 1 associated to new link
                n2=ln(2,grd_ghost_link_closest(c_lnx)) !flownode 2 associated to new link
                !either <n1> or <n2> is the junction node <idx_fm>. We take the other one. 
                if (idx_fm.eq.n1) then
                    grd_fmmv_fmsv(c_ndx)=n2
                    ln_mor(1,grd_ghost_link_closest(c_lnx))=c_ndx
                else
                    grd_fmmv_fmsv(c_ndx)=n1
                    ln_mor(2,grd_ghost_link_closest(c_lnx))=c_ndx
                endif
                
                nd_mor(c_ndx)%lnx=2 !in <nd_mor> only two links are connected to each node. For ghost nodes these are:
                allocate(nd_mor(c_ndx)%ln(2))
                if (kn==1) then 
                   nd_mor(c_ndx)%ln(1)=c_lnx !new ghost link
                   nd_mor(c_ndx)%ln(2)=-grd_ghost_link_closest(c_lnx) !existing link
                   ln_mor(1,c_lnx)=c_ndx
                   ln_mor(2,c_lnx)=grd_fmmv_fmsv(c_ndx)
                else
                   nd_mor(c_ndx)%ln(1)=grd_ghost_link_closest(c_lnx) !existing link
                   nd_mor(c_ndx)%ln(2)=-c_lnx !new ghost link
                   ln_mor(2,c_lnx)=c_ndx
                   ln_mor(1,c_lnx)=grd_fmmv_fmsv(c_ndx)
                endif

                grd_fm_sre(c_ndx)=idx_sre
                
                !node <idx_fm> (at the junction) does not play any role anymore in <nd_mor>. Still, 
                !we save here the index of one of the SRE points associated
                !to it for the sake of writing a value for output. 
                grd_fm_sre(idx_fm)=idx_sre
                
                !add CS at multivalued-ghost flownode
                gridpoint2cross(c_ndx)%num_cross_sections=1
                allocate(gridpoint2cross(c_ndx)%cross(gridpoint2cross(c_ndx)%num_cross_sections))
                gridpoint2cross(c_ndx)%cross(1)=gridpoint2cross(idx_fm)%cross(jpos)
                

           
                !remove CS at junction flownode
                gridpoint2cross(idx_fm)%num_cross_sections=0 !This prevents it is looped in <fm_update_crosssections>
                !gridpoint2cross(idx_fm)%cross(jpos)=-999 !This prevents it is passed in <fm_update_crosssections> -> NO. -999 causes error when parsing the number of CS per node. 
                
           else !not a bifurcation (i.e., boundary)
                !FM1DIMP2DO: This part could be condensed. The same is called above and below but just changed the index. 
                !copy values from <nd_o>
                !nd_mor(idx_fm)%lnx=nd_o(idx_fm)%lnx
                !nd_mor(idx_fm)%ln=nd_o(idx_fm)%ln
                
                grd_fmmv_fmsv(idx_fm)=idx_fm !the closest value is itself
                grd_fm_sre(idx_fm)=idx_sre 
                
                !relate ghost flownode also to <idx_sre>
                idx_l1=abs(nd_o(idx_fm)%ln(1))
                idx_l2=abs(nd_o(idx_fm)%ln(2))
                !there are only two links
                L=max(idx_l1,idx_l2) !the one which is external (the largest of the two) points to the ghost flownode
                L2=min(idx_l1,idx_l2) !the one which is internal (the smallest of the two) points to the internal cell
                n1=ln(1,L)
                n2=ln(2,L)
                idx_n=max(n1,n2) !the maximum flownode is the ghost one
                grd_fm_sre(idx_n)=idx_sre
                
                !link direction for morphodynamics
                link1sign2(L2)=1
                if (kn==1) then 
                   link1sign2(L)=1
                else
                   link1sign2(L)=-1
                endif
                
                !FM1DIMP2DO: I wonder whether we need this or we can use the adapted <gridpoint2cross> in which there is a cross-section for 1:ndx_mor 
                !grd_sre_cs(idx_sre)=gridpoint2cross(idx_fm)%cross(jpos) !cross-section index associated to the FM gridpoint per branch
           endif !(nd(idx_fm)%lnx>2)
           
        else !internal point of a branch, not beginning or end. 
           jpos = 1
           
           !copy values from <nd_o>
           !nd_mor(idx_fm)%lnx=nd_o(idx_fm)%lnx
           !nd_mor(idx_fm)%ln=nd_o(idx_fm)%ln
           
           grd_fmmv_fmsv(idx_fm)=idx_fm !the closest value is itself
           grd_fm_sre(idx_fm)=idx_sre 
           
           !link1sign2()=1 !link direction for morphodynamics
           
           !FM1DIMP2DO: I wonder whether we need this or we can use the adapted <gridpoint2cross> in which there is a cross-section for 1:ndx_mor 
           !grd_sre_cs(idx_sre)=gridpoint2cross(idx_fm)%cross(jpos) !cross-section index associated to the FM gridpoint per branch
        endif  
        

                !FM1DIMP2DO: I wonder whether we need this or we can use the adapted <gridpoint2cross> in which there is a cross-section for 1:ndx_mor 
                grd_sre_cs(idx_sre)=gridpoint2cross(idx_fm)%cross(jpos) !cross-section index associated to the FM gridpoint per branch
                
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

!new dimensions
lnx_mor=c_lnx !store new number of links (considering ghost links)
lnxi_mor=lnx_mor !there are no ghosts in SRE
ndx_mor=c_ndx !store new number of flow nodes (considering multivaluedness)
ndxi_mor=ndx_mor !there are no ghosts in SRE
ndkx_mor=ndx_mor





!ndkx=ndx_mor !used to preallocate <ucxq_mor> and similar. !Cannot be changed because it is used in output data. The only solution is to specifically reallocate these variables. 

!
!END (LOB)
!


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

!frac_o=frac !needs to be copied before <allocsedtra>

allocate(link1(ndx_mor)) !we allocate with
link1=0
do L = 1,lnx_mor
    k1 = ln_mor(1,L)
    !k2 = ln(2,L)
    link1(k1) = L
enddo

!FM1DIMP2DO: When not having a morhodynamic simulation, morpho variables are not initialized. The best
!would be to <return> in <reallocate_fill>

!nodes
    !allocate
call reallocate_fill(s0     ,grd_fmmv_fmsv,ndx,ndx_mor) 
call reallocate_fill(s1     ,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill(bai_mor,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill(ba     ,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill(rhowat ,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill(ktb    ,grd_fmmv_fmsv,ndx,ndx_mor) !FM1DIMP2DO: It could be better to allocate the wave part with <ndx_mor>. To limit the mess it is now done here. 
call reallocate_fill(bl     ,grd_fmmv_fmsv,ndx,ndx_mor)
!FM1DIMP2DO
!The value of <bl> here is not correct, as it copies the value from the closest node. 
!However, this is not a big issue, as this value is only used for checking whether flow
!depth is above or below threshold. The bed level value for computing flow comes from
!the cross-sections, which are treated independently. Nevertheless, we could here fill
!the right value from cross-sections. 



!call reallocate_fill_pointer(ucxq_mor   ,grd_fmmv_fmsv,ndx,ndx_mor)
!call reallocate_fill_pointer(ucyq_mor   ,grd_fmmv_fmsv,ndx,ndx_mor)
!call reallocate_fill_pointer(hs_mor     ,grd_fmmv_fmsv,ndx,ndx_mor)
!call reallocate_fill_pointer(dzbdt      ,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill_pointer(bfmpar%rksr,grd_fmmv_fmsv,ndx,ndx_mor)
!call reallocate_fill_pointer(pmcrit     ,grd_fmmv_fmsv,ndx,ndx_mor)

!call reallocate_fill_int_pointer(kfsed  ,grd_fmmv_fmsv,ndx,ndx_mor)

!multidimensional nodes

    !copy arrays to temporary array
!dbodsd_fm=dbodsd 
!if (allocated(dbodsd)) then
!    deallocate(dbodsd)
!endif
!allocate(dbodsd(lsedtot,ndx_mor))

    !copy pointers to temporary array
bodsed_o=stmpar%morlyr%state%bodsed
allocate(stmpar%morlyr%state%bodsed(lsedtot,ndx_mor))

    !copy data from nodes existing in FM
do kn=1,ndx
    !arrays sediment
    do ksed=1,lsedtot
        stmpar%morlyr%state%bodsed(ksed,kn)=bodsed_o(ksed,kn)
    enddo !ksed
enddo !kn

    !fill
do kn=ndx+1,ndx_mor    
    !arrays sediment
    do ksed=1,lsedtot
        stmpar%morlyr%state%bodsed(ksed,kn)=bodsed_o(ksed,grd_fmmv_fmsv(kn))
    enddo !ksed
enddo !kl

!links
    !allocate
call reallocate_fill(u1      ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(au      ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(hu      ,grd_ghost_link_closest,lnx,lnx_mor)
call reallocate_fill(dx      ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(dxi     ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(wu      ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(frcu_mor,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(frcu    ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(qa      ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(acl     ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(snu     ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(csu     ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(wu_mor  ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(z0urou  ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(taubxu  ,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill(wcx1    ,grd_ghost_link_closest,lnx,lnx_mor)
call reallocate_fill(wcx2    ,grd_ghost_link_closest,lnx,lnx_mor)
call reallocate_fill(wcy1    ,grd_ghost_link_closest,lnx,lnx_mor)
call reallocate_fill(wcy2    ,grd_ghost_link_closest,lnx,lnx_mor)

call reallocate_fill_int(ifrcutp,grd_ghost_link_closest,lnx,lnx_mor) 
call reallocate_fill_int(kcu,grd_ghost_link_closest,lnx,lnx_mor)

!call reallocate_fill_pointer(e_dzdn  ,grd_ghost_link_closest,lnx,lnx_mor) 

!multidimensional links

    !copy arrays to temporary array
!wcl_fm=wcl 
if (allocated(wcl)) then
    deallocate(wcl)
endif
allocate(wcl(2,lnx_mor)) 
do kl=1,lnx_mor
    do kd=1,2
        wcl(kd,kl)=0.5d0 !each link corresponds to 50% of the flownode. 
    enddo
enddo
do kl=lnxi+1,lnx
    wcl(1,kl)=1.0d0
enddo

!ln_fm=ln
!if (allocated(ln)) then
!    deallocate(ln)
!endif
!allocate(ln(2,lnx_mor))

    !copy pointers to temporary array
!e_sbcn_fm=e_sbcn 
!allocate(e_sbcn(lnx_mor,lsedtot))
!
!e_sbn_fm=e_sbn
!allocate(e_sbn(lnx_mor,lsedtot))

    !copy data from links existing in FM
!do kl=1,lnx
!    !arrays left-right
!    do kd=1,2
!        wcl(kd,kl)=wcl_fm(kd,kl)
!        !ln(kd,kl)=ln_fm(kd,kl)
!    enddo !kd
!    !arrays sediment
!    do ksed=1,lsedtot
!        !e_sbcn(kl,ksed)=e_sbcn_fm(kl,ksed)
!        !e_sbn (kl,ksed)=e_sbn_fm (kl,ksed)
!    enddo !ksed
!enddo !kl
!
!    !fill
!do kl=lnx+1,lnx_mor    
!    !arrays left-right
!    do kd=1,2
!        wcl(kd,kl)=1-wcl(kd,grd_ghost_link_closest(kl)) !it should be equal to 0.5
!        !ln(kd,kl)=ln(kd,grd_ghost_link_closest(kl))
!    enddo !kd
!    !arrays sediment
!    do ksed=1,lsedtot
!        !e_sbcn(kl,ksed)=e_sbcn(grd_ghost_link_closest(kl),ksed)
!        !e_sbn (kl,ksed)=e_sbn (grd_ghost_link_closest(kl),ksed)
!    enddo !ksed
!enddo !kl


!morphodynamics initialization is done before fm1dimp initialization. Hence, we have to reallocate here using <ndx_mor> and <lnx_mor>
!It must be before the calls to <reallocate_~> because some variables (e.g., <ucxq_mor>) are set to the wrong size (i.e., <ndkx>) in <allocsedtra>
!We have to copy <bodsed> before initialization.
!griddim%nmub=ndx_mor !This is used for allocating in <rdsed>. We change it here and then back to the original for not messing up with other calls. -> This does not work because the code is only passed if `if (.not. associated(sedpar%sedd50)) then`
call flow_sedmorinit()
!griddim%nmub=ndx 
!call allocsedtra(sedtra, stmpar%morpar%moroutput, max(kmx,1), stmpar%lsedsus, stmpar%lsedtot, 1, ndx_mor, 1, lnx_mor, stmpar%morpar%nxx, stmpar%morpar%moroutput%nstatqnt)
!call inipointers_erosed()
!call initsedtra(sedtra, stmpar%sedpar, stmpar%trapar, stmpar%morpar, stmpar%morlyr, rhomean, ag, vismol, 1, ndx_mor, ndx_mor, stmpar%lsedsus, stmpar%lsedtot)

!We could do the same trick and call <lnx_mor> in <flow_waveinit>, but some variables have been moved to another module after JR merge. Hence, we reallocate in this routine. 
!call flow_waveinit()

!needs to be after <flow_sedmorinit>, where it is allocated. 
!FM1DIMP2DO: Ideally, these variables are allocated in <flow_sedmorinit>
call reallocate_fill_pointer(pmcrit                   ,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill_pointer(stmpar%morlyr%state%dpsed,grd_fmmv_fmsv,ndx,ndx_mor)
call reallocate_fill_int(kcsmor  ,grd_fmmv_fmsv,ndx,ndx_mor)

ucyq_mor=0d0 !set to 0 once rather than every time step. Somewhere in the code is changed. I have to set it every time step. 

stmpar%morlyr%settings%nmub=ndx_mor

!
!END (FFMA)
!

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
!bfrict(1,i)=cfrchc (1) : Chézy constant
!bfrict(1,i)=cfrchq (2) : Chézy function of discharge
!bfrict(1,i)=cfrchh (3) : Chézy function of water level
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
    idx_crs=grd_sre_cs(ksre)
    
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
    !stmpar%morpar%mornum%pure1d=1
    
    !the most downstream link points inside and we have to consider this for the flux.
    !call init_1dinfo() !<initialize_flow1d_implicit> is called before <init_1dinfo>. We have to call it here and it will not be called again because it will be allocated. 
    !FM1DIMP2DO: I don't know why it fails compiling in case I ask to allocate here. It does not have the allocatable attribute, but neither it does <link1sign> 
    !if (allocated(link1sign2)) then
    !    deallocate(link1sign2)
    !endif 
    !allocate(link1sign(lnx_mor))
    !allocate(link1sign2(lnx_mor)) 
    !do kl=1,lnx_mor 
    !    link1sign(kl)=1
    !    link1sign2(kl)=1
    !    
    !enddo
    !do kl=lnxi+1,lnx1Db !boundary links
    !    n1 = ln(1,kl) 
    !    n2 = ln(2,kl)
    !    nint=min(n1,n2) !from the two cells that this link connects, the minimum is internal, and hence we have data
    !    nout=max(n1,n2) !from the two cells that this link connects, the maximum is extrernal, and it is the one in which we have to set the water level
    !    if (f1dimppar%x(nint).eq.0) then !upstream
    !        link1sign2(kl)=1 
    !    else ! downstream
    !        link1sign2(kl)=-1
    !    endif
    !enddo
    !FM1DIMP2DO: Ghost links will be rechanged in the nodal point relation. Not sure we need them here. 
    !do kn=1,ndx_mor 
    !    n1=grd_fm_sre(nd_mor(kn)%ln(1))
    !    n2=grd_fm_sre(nd_mor(kn)%ln(2))
    !    if (f1dimppar%x(nint).eq.0) then !upstream
    !        link1sign2(kl)=1 
    !    else ! downstream
    !        link1sign2(kl)=-1
    !    endif
    !        
    !    !link1sign2(kl)
    !enddo 
    
endif !jased

!because the <height> is used in the cross-sections of SRE, <shift> cannot be used in cross-section

!
!BEGIN (CHK)
!
!CHecKs

!all cross-sections must be mapped to an SRE gridpoint
do kd=1,ngrid
    idx_findloc=findloc(grd_sre_cs,kd)
    if (idx_findloc(1).eq.0) then
        write (msgbuf, '(a)') 'Not all SRE nodes are related to a cross-section.'
        call err_flush()
        iresult=1
    endif
enddo

!
!END (CHK)
!


    contains

    !FM1DIMP2DO: Move to module
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END MAIN SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine initialize_flow1d_implicit