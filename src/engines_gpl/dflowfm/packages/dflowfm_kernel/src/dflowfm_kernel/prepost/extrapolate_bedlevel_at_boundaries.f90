!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2023.                                
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

! 
! 

subroutine extrapolate_bedlevel_at_boundaries()   
 !use timespace_data
 !use timespace
 !use unstruc_model
 use m_flowgeom
 use m_flow
 use m_netw !  only : xk, yk, zk
 !use m_missing
 !use system_utils, only: split_filename
 !use unstruc_files, only: resolvePath
 !use string_module, only: strcmpi
 !use unstruc_inifields, only: readIniFieldProvider, checkIniFieldFileVersion
 !use dfm_error
 !use unstruc_netcdf
 use m_flowparameters, only: jadpsopt 
 use m_flowexternalforcings, only: zbndu, zbndq, kbndz
 
 implicit none

 logical, external :: timespaceinitialfield_mpi

 logical :: jawel
 logical :: bl_set_from_zkuni = .false.
 integer              :: mxyb, ja, ja1, ja2, method, iprimpos
 integer              :: k, L, k1, k2, mx

 double precision :: dxz, dyz, bedslopex, bedslopey
 
 if (ibedlevtyp == 1) then
    ! To improve: bed levels at boundary to be set from net file, instead of mirroring
    do L = Lnxi + 1, Lnx
        k1 = ln(1,L)
        k2 = ln(2,L)
        if (jaextrbl.eq.1) then 
            !add a check that izbndpos must be 0 (i.e., D3D4 style). 
            !
            !we first copy and then add a correction to `bl(k1)`
            dxz=xz(k1)-xz(k2)
            dyz=yz(k1)-yz(k2)
            !m_fm_erosed::e_dzdn check if this can be somehow used
            bedslopex=-0.007d0/9.81d0 !AD-HOC make input or read from some existing variable.
            bedslopey=0d0 !AD-HOC make input or read from some existing variable.
            if (L.ne.101) then !downstream boundary (H-boundary) mega ad-hoc 1D case!!!!
                !we want that at an H-boundary the bed level of the ghost cell is
                !the correct one, as if the sloping bed would continue. In this 
                !way, when computing the mean of the bed levels at cell centres,
                !the elevation at the dowsntream link will be the correct one and
                !`hu` will be correct if the water level at the cell centre of the 
                !ghost cell (i.e., using `izbndpos=0` is `\Delta x` down from the 
                !correct value. Note that it is `\Delta x` and not `\Delta x/2`.
                bl(k1)=bl(k2)+bedslopex*dxz+bedslopey*dyz !u/s Q
                !bl(k1)=bl(k2)-bedslopex*dxz+bedslopey*dyz !u/s etaw
            else !upstream boundary (u-boundary)
                !we want that at a u-boundary the bed level of the ghost cell is 
                !the same as that of the second-inside cell. In this way, when 
                !computing the mean of the bed levels at cell centres, the elevation
                !at the upstream link will be the same as that of the second link. 
                !As the water level of the first inside cell is (1) downwinded for
                !computing `hu` at the second link and (2) copied to the upstream 
                !ghost cell, we will obtain the same correct `hu` at the first link.
                !
                !For a case of a sloping bed in the positive x direction and a u-boundary
                !imposed upstream, `dxz` will be negative upstream, hence obtaining the
                !elevation of the second inside cell. 
                !bl(k1)=bl(k2)-bedslopex*dxz+bedslopey*dyz !d/s etaw
                bl(k1)=bl(k2)+bedslopex*dxz+bedslopey*dyz !d/s Q
            endif
        endif !jaextrbl
    enddo !L
 endif !ibedlevtyp==1

 end subroutine extrapolate_bedlevel_at_boundaries 
