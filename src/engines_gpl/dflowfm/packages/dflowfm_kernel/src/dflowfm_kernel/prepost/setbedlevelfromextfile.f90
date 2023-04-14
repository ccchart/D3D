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

subroutine setbedlevelfromextfile()    ! setbedlevels()  ! check presence of old cell centre bottom level file
 use timespace_data
 use timespace
 use unstruc_model
 use m_flowgeom
 use m_flow
 use m_netw !  only : xk, yk, zk
 use m_missing
 use system_utils, only: split_filename
 use unstruc_files, only: resolvePath
 use string_module, only: strcmpi
 use unstruc_inifields, only: readIniFieldProvider, checkIniFieldFileVersion
 use dfm_error
 use unstruc_netcdf
 use m_flowparameters, only: jadpsopt 
 use m_flowexternalforcings, only: zbndu, zbndq, kbndz
 
 implicit none

 logical, external :: timespaceinitialfield_mpi

 logical :: jawel
 logical :: bl_set_from_zkuni = .false.
 integer              :: mxyb, ja, ja1, ja2, method, iprimpos
 integer              :: k, L, k1, k2, mx
 integer, allocatable :: kcc(:), kc1D(:), kc2D(:)
 integer              :: ibathyfiletype
 integer              :: kc_size_store

 character(len=256) :: filename
 character(len=64)  :: varname
! character(len=1)   :: operand
! double precision   :: transformcoef(25) !< Transform coefficients a+b*x

 type(tree_data),  pointer       :: inifield_ptr        !< tree of inifield-file's [Initial] or [Parameter] blocks
 type(tree_data),  pointer       :: node_ptr
 integer                         :: istat
 integer                         :: num_items_in_file
 integer, parameter              :: ini_key_len   = 32
 integer, parameter              :: ini_value_len = 256
 character(len=ini_key_len)      :: groupname
 character(len=255)              :: fnam
 character(len=255)              :: basedir
 integer :: i, iLocType

 double precision :: dxz, dyz, bedslopex, bedslopey
 
 kc_size_store = 0
 inifield_ptr => null()

 ! Attempt to read cell centred bed levels directly from net file:
 call setbedlevelfromnetfile()
 call mess(LEVEL_INFO, 'setbedlevelfromextfile: Using bedlevel as specified in net-file.')

 ! ibedlevtyp determines from which source data location the bed levels are used to derive bobs and bl.
 ! These types need to be mapped to one of three possible primitive locations (center/edge/corner).
 select case (ibedlevtyp)
 case (1)       ! position = waterlevelpoint, cell centre
    iprimpos = 2 ; mx = max(numk, ndx)
 case (2)       ! position = velocitypoint, cellfacemid
    iprimpos = 1 ; mx = max(numk, lnx)
 case (3,4,5,6) ! position = netnode, cell corner
    iprimpos = 3 ; mx = numk
 end select

 if (mext /= 0 .or. len_trim(md_inifieldfile) > 0) then
    ! 0.a Prepare masks for 1D/2D distinctions
    kc_size_store = size(kc)
    allocate(kcc(mx),kc1d(mx),kc2d(mx)) ; kcc = 1; kc1D = 0 ; kc2D = 0
    call realloc(kc, mx, keepExisting = .false., fill = 0)

    do L = 1, numL1D
       if (kn(3,L) == 1 .or. kn(3,L) == 6) then ! TODO: AvD: why not also type 3/4/5/7?
           k1 = kn(1,L) ; k2 = kn(2,L)
           if (nmk(k1) > 1) kc1D(k1) = 1
           if (nmk(k2) > 1) kc1D(k2) = 1
       endif
    enddo

    if (iprimpos == 3) then
       do L = 1, numL
          if (kn(3,L) == 2) then
              k1 = kn(1,L) ; k2 = kn(2,L)
              kc2D(k1) = 1
              kc2D(k2) = 1
          endif
       enddo
    else if (iprimpos == 1) then
       kc2D(lnx1d+1:lnxi) = 1
    else if (iprimpos == 2) then
       kc2D(1:ndx2D) = 1
    endif

    ja = 0
    ja1= 0
    ja2= 0
    ! 0.b Prepare loop across old ext file:
    if (mext /= 0) then
       rewind(mext)
       ja1 = 1
    end if

    ! 0.c Prepare loop across new initial field file:
    if (len_trim(md_inifieldfile) > 0) then
       call tree_create(trim(md_inifieldfile), inifield_ptr)
       call prop_file('ini',trim(md_inifieldfile),inifield_ptr,istat)
       call split_filename(md_inifieldfile, basedir, fnam)
       istat = checkIniFieldFileVersion(md_inifieldfile, inifield_ptr)
       if (istat /= DFM_NOERR) then
          num_items_in_file = 0
       end if
       if (associated(inifield_ptr%child_nodes)) then
           num_items_in_file = size(inifield_ptr%child_nodes)
       endif
       if (num_items_in_file > 0) then
          i  = 1
          ja2 = 1
       end if
    end if

    ! Trick: loop across the 2 supported file types (*.ext and *.ini), most inner do-loop code is the same for both.
bft:do ibathyfiletype=1,2
    if (ibathyfiletype == 1) then
       call split_filename(md_extfile,      basedir, fnam) ! Remember base dir of *.ext file, to resolve all refenced files below w.r.t. that base dir.
       if (ja1 == 1) then
          ja = 1
       end if
    else if (ibathyfiletype == 2) then
       call split_filename(md_inifieldfile, basedir, fnam) ! Remember base dir of *.ini file, to resolve all refenced files below w.r.t. that base dir.
       if (ja2 == 1) then
          ja = 1
       end if
    end if

    do while (ja == 1)
       if (ibathyfiletype == 1) then       ! read *.ext file
          call delpol()
          call readprovider(mext,qid,filename,filetype,method,operand,transformcoef,ja,varname)
       else if (ibathyfiletype == 2) then  ! read *.ini file
          if (i > num_items_in_file) then
             ja = 0
             exit
          end if
          node_ptr => inifield_ptr%child_nodes(i)%node_ptr
          call readIniFieldProvider(md_inifieldfile, node_ptr,groupname,qid,filename,filetype,method,iLocType,operand,transformcoef,ja,varname) !,smask, maxSearchRadius)
          i = i + 1
          if (.not. strcmpi(groupname, 'Initial')) then
             cycle
          end if
       end if

       ! Initialize bedlevel based on the read provider info
       if (ja == 1) then
          call resolvePath(filename, basedir, filename)
          if (index(qid,'bedlevel') > 0 .and. ibathyfiletype == 1 .and. len_trim(md_inifieldfile) > 0) then
             ! Don't support bedlevel in *.ext file when there is ALSO a *.ini file.
             call mess(LEVEL_WARN, 'Bed level info should be defined in file '''//trim(md_inifieldfile)//'''. Quantity '//trim(qid)//' ignored in external forcing file '''//trim(md_extfile)//'''.')
             cycle bft ! Try ini field file next
          end if
          success = .true.
          if (strcmpi(qid, 'bedlevel1D') .or. (strcmpi(qid, 'bedlevel') .and. ibathyfiletype == 2 .and. iLocType == ILATTP_1D)) then
             call mess(LEVEL_INFO, 'setbedlevelfromextfile: Setting 1D bedlevel from file '''//trim(filename)//'''.')
             kc(1:mx) = kc1D
             success = timespaceinitialfield_mpi(xk, yk, zk, numk, filename, filetype, method, operand, transformcoef, 3, kc) ! see meteo module
          else if (strcmpi(qid,'bedlevel', 8)) then
             if ((strcmpi(qid, 'bedlevel') .and. ibathyfiletype == 1) .or. (strcmpi(qid, 'bedlevel') .and. ibathyfiletype == 2 .and. iLocType == ILATTP_ALL))  then
                call mess(LEVEL_INFO, 'setbedlevelfromextfile: Setting both 1D and 2D bedlevel from file '''//trim(filename)//'''.')
                kc(1:mx) = kcc
             else if (strcmpi(qid, 'bedlevel2D') .or. (strcmpi(qid, 'bedlevel') .and. ibathyfiletype == 2 .and. iLocType == ILATTP_2D)) then
                call mess(LEVEL_INFO, 'setbedlevelfromextfile: Setting 2D bedlevel from file '''//trim(filename)//'''.')
                kc(1:mx) = kc2D
             endif

             if (ibedlevtyp == 3) then
                success = timespaceinitialfield_mpi(xk, yk, zk, numk, filename, filetype, method, operand, transformcoef, iprimpos, kc) ! see meteo module
             else if (ibedlevtyp == 2) then
                success = timespaceinitialfield_mpi(xu, yu, blu, lnx, filename, filetype, method, operand, transformcoef, iprimpos, kc) ! see meteo module
             else if (ibedlevtyp == 1) then
                success = timespaceinitialfield_mpi(xz, yz, bl, ndx, filename, filetype, method, operand, transformcoef, iprimpos, kc) ! see meteo module
             endif
          endif
          if ( .not. success) then
             call mess(LEVEL_FATAL, "Error reading " // trim(qid) // " from " // trim(filename) // ".")
          end if
       endif

    end do ! ja==1 provider loop
    end do bft ! ibathyfiletype=1,2

    ! Clean up *.ext file
    if (mext /= 0) then
       rewind (mext)
    endif   

    ! Clean up *.ini file.
    call tree_destroy(inifield_ptr)

    ! Interpreted values for debugging.
    if ( md_jasavenet == 1 ) then
!      save network
       select case (ibedlevtyp)
          case (3,4,5,6) ! primitime position = netnode, cell corner
             call unc_write_net('DFM_interpreted_network_'//trim(md_ident)//'_net.nc')
       end select
    end if

    deallocate(kcc,kc1d,kc2d)

 endif



 if (ibedlevtyp == 1) then
    do k = 1, ndxi
        if (bl(k) == dmiss) then
           bl(k) = zkuni
           bl_set_from_zkuni = .true.
        endif
    enddo
    if (bl_set_from_zkuni) then
        call mess(LEVEL_INFO, 'setbedlevelfromextfile: Unspecified bedlevels replaced using value from BedlevUni.')
    endif

    ! To improve: bed levels at boundary to be set from net file, instead of mirroring
    do L = Lnxi + 1, Lnx
        k1 = ln(1,L)
        k2 = ln(2,L)
        bl(k1) = bl(k2) !original
        if (jaextrbl.eq.1) then 
            !add a check that izbndpos must be 0 (i.e., D3D4 style). 
            !
            !we first copy and then add a correction to `bl(k1)`
            dxz=xz(k1)-xz(k2)
            dyz=yz(k1)-yz(k2)
            !m_fm_erosed::e_dzdn check if this can be somehow used
            bedslopex=-0.007d0/9.81d0 !AD-HOC make input or read from some existing variable.
            bedslopey=0d0 !AD-HOC make input or read from some existing variable.
            if (L.ne.101) then !H-boundary mega ad-hoc 1D case!!!!
                !we want that at an H-boundary the bed level of the ghost cell is
                !the correct one, as if the sloping bed would continue. In this 
                !way, when computing the mean of the bed levels at cell centres,
                !the elevation at the dowsntream link will be the correct one and
                !`hu` will be correct if the water level at the cell centre of the 
                !ghost cell (i.e., using `izbndpos=0` is `\Delta x` down from the 
                !correct value. Note that it is `\Delta x` and not `\Delta x/2`.
                bl(k1)=bl(k2)+bedslopex*dxz+bedslopey*dyz !u/s Q
                !bl(k1)=bl(k2)-bedslopex*dxz+bedslopey*dyz !u/s etaw
            else !u-boundary
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
                !
                !When the u-boundary is at the downstream end (check on negative velocity)
                !then we simply have to extrapolate the bed level to the correct value.
                bl(k1)=bl(k2)+bedslopex*dxz+bedslopey*dyz !d/s Q
            endif
        endif
        
    enddo
    call mess(LEVEL_INFO, 'setbedlevelfromextfile: Mirroring input bedlevels at open boundaries.')

 endif

 if ( kc_size_store.gt.0 ) then
    call realloc(kc, kc_size_store, keepExisting = .false., fill = 0)
 endif


 end subroutine setbedlevelfromextfile ! setbottomlevels
