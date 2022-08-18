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

 subroutine alloclinktocenterweights()                 ! set center related linkxy weights

 use m_flow
 use m_netw
 use m_flowgeom
 use m_sferic
 use m_longculverts
 implicit none

 double precision       :: wud, wuL1, wuL2, wuk, cs, sn
 integer                :: k, L, ierr, n, kk, n12, lnxmax
 integer                :: k1, k2, k3, k4, nn, LL, jaclosedcorner
 integer                :: ilongc, L1dlink

 double precision       :: xloc, yloc, beta, aa1, wcw, alf

 double precision, external :: lin2nodx, lin2nody

 if ( allocated (wcx1) )  deallocate(wcx1)
 if ( allocated (wcy1) )  deallocate(wcy1)
 if ( allocated (wcx2) )  deallocate(wcx2)
 if ( allocated (wcy2) )  deallocate(wcy2)
 if ( allocated (wcxy) )  deallocate(wcxy)
 if ( allocated (wcL) )   deallocate(wcL)
 if ( allocated (wc) )    deallocate(wc)

 allocate ( wcx1(lnx) , stat  = ierr) ; wcx1 = 0
 call aerr('wcx1(lnx)', ierr, lnx)
 allocate ( wcy1(lnx) , stat  = ierr) ; wcy1 = 0
 call aerr('wcy1(lnx)', ierr, lnx)
 allocate ( wcx2(lnx) , stat  = ierr) ; wcx2 = 0
 call aerr('wcx2(lnx)', ierr, lnx)
 allocate ( wcy2(lnx) , stat  = ierr) ; wcy2 = 0
 call aerr('wcy2(lnx)', ierr, lnx)
 allocate ( wcxy (2,ndx) , stat  = ierr) ; wcxy  = 0
 call aerr('wcxy (2,ndx)', ierr, 2*ndx)
 allocate ( wcL  (2,Lnx) , stat  = ierr) ; wcL   = 0
 call aerr('wcL  (2,Lnx)', ierr, 2*Lnx)
 allocate ( wc     (ndx) , stat  = ierr) ; wc    = 0
 call aerr('wc     (ndx)', ierr, ndx)
 
 lnxmax = 0
 do n   = 1, mxwalls                                        ! wall contribution to scalar linktocenterweights
    k1  = walls(1,n)
    lnxmax = max(lnxmax,  nd(k1)%lnx)
    call realloc(wwLL, lnxmax, keepExisting = .false.)
 enddo

 end subroutine alloclinktocenterweights
