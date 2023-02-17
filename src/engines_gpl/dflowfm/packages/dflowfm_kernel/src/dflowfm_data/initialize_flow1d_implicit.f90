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
use m_initialize_flow1d_implicit
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
integer                                  , pointer :: juer
integer                                  , pointer :: nlyr

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
!double precision, dimension(:,:)         , pointer :: bodsed
!double precision, dimension(:,:)         , pointer :: thlyr

!double precision, dimension(:,:,:)       , pointer :: msed

type(tnode)    , allocatable :: nd_o(:) !Copy of <nd> for reworking <nd>
!type(tnode)    , pointer     :: nd_mor(:) !Modified <nd> for <bott3d>
type(t_gridp2cs), dimension(:), allocatable :: gridpoint2cross_o

!debug
integer, pointer :: fm1dimp_debug_k1

!output
integer, intent(out) :: iresult !< Error status, DFM_NOERR==0 if succesful.

!local
integer :: kbr, knod, k1, k2, kbe, klnx, ksre, kn, kl, kd, ksed, klyr !FM1DIMP2DO: make the variables names consistent
integer :: ndx_max, lnx_max
integer :: c_lnx, c_ndx !counters
integer :: idx_crs, idx_sre, idx_fm !indices
integer :: n1, n2, nint, nout, pointscount, jpos
integer :: table_number
integer :: idx_fr, idx_to
integer :: idx_i, idx_f, nl, L, L2, idx_l1, idx_l2, idx_sre_p, idx_sre_c, idx_n
integer :: j
integer :: stat

character(len=512) :: msg

integer, dimension(1) :: idx_findloc
integer, dimension(:), allocatable :: grd_fm_sre2
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

!double precision, allocatable, dimension(:,:) :: wcl_fm
!double precision, allocatable, dimension(:,:) :: e_sbcn_fm
!double precision, allocatable, dimension(:,:) :: e_sbn_fm
double precision, allocatable, dimension(:,:) :: bodsed_o
!double precision, allocatable, dimension(:,:) :: frac_o
double precision, allocatable, dimension(:,:) :: thlyr_o
double precision, allocatable, dimension(:,:) :: sedshort_o
double precision, allocatable, dimension(:,:) :: svfrac_o
double precision, allocatable, dimension(:,:) :: preload_o

double precision, allocatable, dimension(:,:,:) :: msed_o

!!
!! POINT
!!

!pointer cannot be before the array is allocated, here only non-allocatable arrays

!f1dimppar
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
juer                   => f1dimppar%juer
bfrict => f1dimppar%bfrict
ntab => f1dimppar%ntab
qbdpar       => f1dimppar%qbdpar
hbdpar       => f1dimppar%hbdpar
table => f1dimppar%table
node => f1dimppar%node
numnod => f1dimppar%numnod
nodnod => f1dimppar%nodnod

!stmpar
if (jased > 0 .and. stm_included) then !passing if no morphpdynamics
nlyr                   => stmpar%morlyr%SETTINGS%NLYR
endif

!!
!! CALC
!!

call inifm1dimp_ini(iresult) !INItialize arrays
call inifm1dimp_lob(iresult) !Loop On Branches
call inifm1dimp_faap(iresult)!Fill Arrays that need Additional Point
call inifm1dimp_fic(iresult) !Fill Initial Condition
call inifm1dimp_fbrp(iresult)!Fill Branch PRoperties 
call inifm1dimp_fbc(iresult) !Fill Boundary Conditions
call inifm1dimp_fnod(iresult)!Fill NODes
call inifm1dimp_chk(iresult) !CHecK

end subroutine initialize_flow1d_implicit