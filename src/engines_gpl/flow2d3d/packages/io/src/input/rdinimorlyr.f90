subroutine rdinimorlyr(filcomp   ,lyrfrac   ,thlyr     ,cdryb     , &
                     & lsedtot   ,mmax      ,nlyr      ,nmax      , &
                     & nmaxus    ,nmmax     ,lundia    ,kcs       , &
                     & icx       ,icy       ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!!--description-----------------------------------------------------------------
!
! Reads attribute file for 3D morphology computation
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use bedcomposition_module
    use properties
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
!
! Global variables
!
    integer                                                , intent(in)  :: icx
    integer                                                , intent(in)  :: icy
    integer                                                , intent(in)  :: lsedtot !  Description and declaration in iidim.f90
    integer                                                              :: lundia  !  Description and declaration in inout.igs
    integer                                                , intent(in)  :: mmax    !  Description and declaration in iidim.f90
    integer                                                , intent(in)  :: nlyr    !  Description and declaration in iidim.f90
    integer                                                , intent(in)  :: nmax    !  Description and declaration in iidim.f90
    integer                                                , intent(in)  :: nmaxus  !  Description and declaration in iidim.f90
    integer                                                              :: nmmax   !  Description and declaration in iidim.f90
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)             , intent(in)  :: kcs     !  Description and declaration in iidim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub,nlyr,lsedtot), intent(out) :: lyrfrac
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub,nlyr)        , intent(out) :: thlyr
    real(fp), dimension(lsedtot)                           , intent(in)  :: cdryb   !  Description and declaration in rjdim.f90
    character(*)                                                         :: filcomp
!
! Local variables
!
    integer                               :: i
    integer                               :: ilyr
    integer                               :: istat
    integer                               :: l
    integer                               :: length
    integer                               :: nm
    integer                               :: nm2
    integer                               :: nmlb
    integer                               :: nmub
    real(fp)                              :: cdrybavg
    real(fp)                              :: fraction
    real(fp)                              :: sedbed
    real(fp)                              :: sedthick
    real(fp)                              :: rmissval
    real(fp)                              :: totfrac
    real(fp), dimension(:,:), allocatable :: rtemp
    logical                               :: anyfrac
    logical                               :: anysedbed
    logical                               :: err
    logical                               :: ex
    logical                               :: error
    character(10)                         :: lstr
    character(10)                         :: versionstring
    character(11)                         :: fmttmp   ! Format file ('formatted  ') 
    character(80)                         :: parname
    character(80)                         :: layertype
    character(256)                        :: filename
    character(300)                        :: message
    type(tree_data), pointer              :: mor_ptr
    type(tree_data), pointer              :: layer_ptr
!
!! executable statements -------------------------------------------------------
!
    rmissval      = -999.0_fp
    fmttmp        = 'formatted'
    nmlb          = gdp%d%nmlb
    nmub          = gdp%d%nmub
    error         = .false.
    !
    ! Create Initial Morphology branch in input tree
    !
    call tree_create_node( gdp%input_tree, 'Initial Morphology', mor_ptr )
    call tree_put_data( mor_ptr, transfer(trim(filcomp),node_value), 'STRING' )
    !
    ! Put mor-file in input tree
    !
    call prop_file('ini', trim(filcomp), mor_ptr, istat)
    if (istat /= 0) then
       select case (istat)
       case(1)
          call prterr(lundia, 'G004', trim(filcomp))
       case(3)
          call prterr(lundia, 'G006', trim(filcomp))
       case default
          call prterr(lundia, 'G007', trim(filcomp))
       endselect
       call d3stop(1, gdp)
    endif
    !
    ! Check version number of mor input file
    !
    versionstring = 'n.a.'
    call prop_get_string(mor_ptr, 'BedCompositionFileInformation', 'FileVersion', versionstring)
    if (trim(versionstring) == '01.00') then
       !
       ! reset fractions to zero
       !
       lyrfrac = 0.0_fp
       !
       ! allocate temporary array
       !
       ilyr = 0
       allocate(rtemp(nmlb:nmub,lsedtot), stat = istat)
       if (istat /= 0) then
          call prterr(lundia, 'U021', 'RdIniMorLyr: memory alloc error')
          call d3stop(1, gdp)
       endif
       !
       do i = 1, size(mor_ptr%child_nodes)
          !
          ! Does sed_ptr contain a child with name 'Layer' (converted to lower case)?
          !
          layer_ptr => mor_ptr%child_nodes(i)%node_ptr
          parname = tree_get_name( layer_ptr )
          call small(parname, len(parname))
          if ( trim(parname) /= 'layer') cycle
          !
          ! Increment ilyr, but do not exceed nlyr
          !
          ilyr = min(nlyr, ilyr+1)
          !
          ! Initialize/reset the temporary array
          !
          do l = 1, lsedtot
             do nm = 1, nmmax
                rtemp(nm, l) = 0.0_fp
             enddo
          enddo
          !
          ! Layer group found, scan it for the layer composition
          !
          layertype = ' '
          call prop_get_string(layer_ptr, '*', 'Type', layertype)
          call small(layertype, len(layertype))
          if (layertype == ' ') then
             !
             ! no Type field found
             !
             write (message,'(a,i2,2a)') 'No type specified for layer ', ilyr, &
                                       & ' in file ', trim(filcomp)
             call prterr(lundia, 'U021', trim(message))
             call d3stop(1, gdp)          
          elseif (layertype == 'mass fraction' .or. &
                & layertype == 'volume fraction') then
             !
             ! mass or volume fraction and layer thickness specified
             !
             parname  = 'Thick'
             filename = ' '
             call prop_get_string(layer_ptr, '*', parname, filename)
             !
             ! Intel 7.0 crashes on an inquire statement when file = ' '
             !
             if (filename == ' ') filename = 'dummyname'
             inquire (file = filename, exist = ex)
             if (.not. ex) then
                !
                ! Constant thickness
                !
                sedbed = rmissval
                call prop_get(layer_ptr, '*', parname, sedbed)
                if (comparereal(sedbed,rmissval) == 0) then
                   write (message,'(a,i2,2a)')  &
                       & 'Missing Thick keyword for layer ', ilyr, &
                       & ' in file ', trim(filcomp)
                   call prterr(lundia, 'U021', trim(message))
                   call d3stop(1, gdp)          
                endif
                do nm = 1, nmmax
                   thlyr(nm, ilyr) = sedbed
                enddo
             else
                !
                ! Spatially varying thickness
                !
                call depfil(lundia    ,error     ,filename, fmttmp    , &
                          & mmax      ,nmax      ,nmaxus  , &
                          & thlyr(nmlb,ilyr)     ,gdp       )
                if (error) then
                   write (message,'(3a,i2,2a)')  &
                       & 'errorreading thickness from ', trim(filename), &
                       & ' for layer ', ilyr, ' in file ', trim(filcomp)
                   call prterr(lundia, 'U021', trim(message))
                   call d3stop(1, gdp)          
                endif
             endif
             !
             anyfrac = .false.
             totfrac = 0.0_fp
             do l = 1, lsedtot
                !
                ! Scan file for fractions
                !
                write(lstr,'(i10)') l
                length = 10
                call noextspaces(lstr, length)
                !
                ! Keyword SedBed<i> may not be used when layertype is fraction
                !
                parname  = 'SedBed' // trim(lstr)
                filename = ' '
                call prop_get_string(layer_ptr, '*', parname, filename)
                if (filename /= ' ') then
                   write (message,'(7a,i2,2a)')  &
                       & 'Use Fraction' ,trim(lstr), ' instead of SedBed', &
                       & trim(lstr), ' for ', trim(layertype), ' layer ', &
                       & ilyr, ' in file ', trim(filcomp)
                   call prterr(lundia, 'U021', trim(message))
                   call d3stop(1, gdp)          
                endif
                !
                parname  = 'Fraction' // trim(lstr)
                filename = ' '
                call prop_get_string(layer_ptr, '*', parname, filename)
                !
                ! Intel 7.0 crashes on an inquire statement when file = ' '
                !
                if (filename == ' ') filename = 'dummyname'
                inquire (file = filename, exist = ex)
                if (.not. ex) then
                   !
                   ! Constant fraction
                   !
                   fraction = 0.0_fp
                   call prop_get(layer_ptr, '*', parname, fraction)
                   if (fraction > 0.0_fp) anyfrac = .true.
                   do nm = 1, nmmax
                      rtemp(nm, l) = fraction
                   enddo
                else
                   !
                   ! Spatially varying fraction
                   !
                   anyfrac = .true.
                   call depfil(lundia    ,error     ,filename, fmttmp    , &
                             & mmax      ,nmax      ,nmaxus  , &
                             & rtemp(nmlb,l)        ,gdp       )
                   if (error) then
                      write (message,'(a,i2,3a,i2,2a)')  &
                          & 'errorreading fraction ', l, 'from ', &
                          & trim(filename), ' for layer ', ilyr, ' in file ', &
                          & trim(filcomp)
                      call prterr(lundia, 'U021', trim(message))
                      call d3stop(1, gdp)          
                   endif
                endif
             enddo
             !
             ! Check if we have found any information that makes sense.
             !
             if (.not. anyfrac) then
                write (message,'(a,i2,2a)')  &
                    & 'No Fraction keywords found for layer ', ilyr, &
                    & ' in file ', trim(filcomp)
                call prterr(lundia, 'U021', trim(message))
                call d3stop(1, gdp)          
             endif
             !
             ! Check validity of input data.
             !
             do nm = 1, nmmax
                if (kcs(nm) == 1) then
                   !
                   ! At an internal point the composition of all layers is important.
                   ! Check the values carefully before continuing.
                   !
                   totfrac = 0.0_fp
                   do l = 1, lsedtot
                      if (rtemp(nm, l) < 0.0_fp) then
                         write (message,'(2a,i2,a,i2,3a,i0)')  &
                             & 'Negative ', trim(layertype), l, ' in layer ', &
                             & ilyr, ' in file ', trim(filcomp), ' at nm=', nm
                         call prterr(lundia, 'U021', trim(message))
                         call d3stop(1, gdp)          
                      elseif (rtemp(nm, l) > 1.0_fp) then
                         write (message,'(a,i2,a,i2,3a,i0)')  &
                             & trim(layertype), l, ' bigger than 1 in layer ', &
                             & ilyr, ' in file ', trim(filcomp), ' at nm=', nm
                         call prterr(lundia, 'U021', trim(message))
                         call d3stop(1, gdp)
                      endif
                      totfrac = totfrac + rtemp(nm,l)
                   enddo
                   if (abs(totfrac-1.0_fp) > 1e-4_fp) then
                      write (message,'(3a,i2,3a,i0)')  &
                          & 'Sum of ', trim(layertype), ' not equal to 1 in layer ', &
                          & ilyr, ' in file ', trim(filcomp), ' at nm=', nm
                      call prterr(lundia, 'U021', trim(message))
                      call d3stop(1, gdp)          
                   else
                      totfrac = 0.0_fp
                      do l = 1, lsedtot-1
                         totfrac = totfrac + rtemp(nm, l)
                      enddo
                      rtemp(nm, lsedtot) = 1.0_fp - totfrac
                   endif
                elseif (kcs(nm) == 2 .and. ilyr == 1) then
                   !
                   ! At an open boundary only the composition of the transport layer
                   ! is important. If it is not valid, mark the data as dummy data:
                   ! the data will be overwritten with data coming from the neighbouring
                   ! internal point.
                   !
                   totfrac = 0.0_fp
                   err     = .false.
                   do l = 1, lsedtot
                      if (rtemp(nm, l) < 0.0_fp .or. rtemp(nm, l) > 1.0_fp) err=.true.
                      totfrac = totfrac + rtemp(nm,l)
                   enddo
                   if (comparereal(totfrac,1.0_fp) /= 0) err=.true.
                   if (err) then
                      !
                      ! dummy
                      !
                      rtemp(nm, 1) = -1.0_fp
                   endif
                else
                   !
                   ! Point/layer that will never be used: don't care about the
                   ! values. Just replace whatever was read by something valid.
                   !
                   do l = 1, lsedtot
                      rtemp(nm, l) = 0.0_fp
                   enddo
                   thlyr(nm, ilyr) = 0.0_fp
                endif
             enddo
             !
             ! Copy RTEMP data to open boundary points that have not
             ! yet been assigned valid data.
             !
             do nm = 1, nmmax
                if (kcs(nm) == 2 .and. rtemp(nm,1) < 0.0_fp) then
                   if (kcs(nm-icx) == 1) then
                      ! ndm
                      nm2 = nm - icx
                   elseif (kcs(nm+icx) == 1) then
                      ! num
                      nm2 = nm + icx
                   elseif (kcs(nm-icy) == 1) then
                      ! nmd
                      nm2 = nm - icy
                   else
                      ! nmu
                      nm2 = nm + icy
                   endif
                   do l = 1, lsedtot
                      rtemp(nm, l) = rtemp(nm2, l)
                   enddo
                   thlyr(nm, ilyr) = thlyr(nm2, ilyr)
                endif
             enddo
             !
             ! convert mass fractions into volume fractions
             !
             if (layertype == 'mass fraction') then
                do nm = 1, nmmax
                   cdrybavg = 0.0_fp
                   do l = 1, lsedtot
                      cdrybavg = cdrybavg + rtemp(nm, l)/cdryb(l)
                   enddo
                   if (cdrybavg > 0.0_fp) then
                      cdrybavg = max(cdrybavg, 1.0e-8_fp)
                      cdrybavg = 1.0_fp / cdrybavg
                      do l = 1, lsedtot
                         rtemp(nm, l) = rtemp(nm, l) * cdrybavg / cdryb(l)
                      enddo
                   endif
                enddo
             endif
             !
             ! add thicknesses in lyrfrac
             !
             do l = 1, lsedtot
                do nm = 1, nmmax
                   lyrfrac(nm, ilyr, l) = lyrfrac(nm, ilyr, l) + rtemp(nm,l)*thlyr(nm, ilyr)
                enddo
             enddo
             !
          elseif (layertype == 'sediment mass' .or. &
                & layertype == 'sediment thickness') then
             !
             ! sediment mass as specified in sediment input file
             !
             anysedbed = .false.
             do l = 1, lsedtot
                !
                ! Scan file for sediment masses
                !
                write(lstr,'(i10)') l
                length = 10
                call noextspaces(lstr, length)
                !
                ! Keyword Fraction<i> may not be used when layertype is sediment
                !
                parname  = 'Fraction' // trim(lstr)
                filename = ' '
                call prop_get_string(layer_ptr, '*', parname, filename)
                if (filename /= ' ') then
                   write (message,'(6a,i2,2a)')  &
                       & 'Use SedBed', trim(lstr), ' instead of Fraction', &
                       & trim(lstr), ' for ', trim(layertype), ' layer ', &
                       & ilyr, ' in file ', trim(filcomp)
                   call prterr(lundia, 'U021', trim(message))
                   call d3stop(1, gdp)
                endif
                !
                parname  = 'SedBed' // trim(lstr)
                filename = ' '
                call prop_get_string(layer_ptr, '*', parname, filename)
                !
                ! Intel 7.0 crashes on an inquire statement when file = ' '
                !
                if (filename == ' ') filename = 'dummyname'
                inquire (file = filename, exist = ex)
                if (.not. ex) then
                   !
                   ! Constant thickness or mass
                   !
                   sedbed = 0.0_fp
                   call prop_get(layer_ptr, '*', parname, sedbed)
                   if (sedbed > 0.0_fp) anysedbed = .true.
                   do nm = 1, nmmax
                      rtemp(nm, l) = sedbed
                   enddo
                else
                   !
                   ! Spatially varying thickness or mass
                   !
                   anysedbed = .true.
                   call depfil(lundia    ,error     ,filename, fmttmp    , &
                             & mmax      ,nmax      ,nmaxus  , &
                             & rtemp(nmlb,l)        ,gdp       )
                   if (error) then
                      write (message,'(5a,i2,2a)')  &
                          & 'errorreading ', layertype, '  from ', trim(filename), &
                          & ' for layer ', ilyr, ' in file ', trim(filcomp)
                      call prterr(lundia, 'U021', trim(message))
                      call d3stop(1, gdp)          
                   endif
                endif
             enddo
             !
             ! Check if we have found any information that makes sense.
             !
             if (.not. anysedbed) then
                write (message,'(a,i2,2a)')  &
                    & 'No SedBed keywords found for layer ' ,ilyr, &
                    & ' in file ' ,trim(filcomp)
                call prterr(lundia, 'U021', trim(message))
                call d3stop(1, gdp)          
             endif
             !
             ! Check validity of input data.
             !
             do nm = 1, nmmax
                if (kcs(nm) == 1) then
                   !
                   ! At an internal point the composition of all layers is important.
                   ! Check the values carefully before continuing.
                   !
                   do l = 1, lsedtot
                      if (rtemp(nm, l) < 0.0_fp) then
                         write (message,'(2a,i2,a,i2,3a,i0)')  &
                             & 'Negative ', trim(layertype), l, ' in layer ', &
                             & ilyr, ' in file ', trim(filcomp), ' at nm=', nm
                         call prterr(lundia, 'U021', trim(message))
                         call d3stop(1, gdp)          
                      endif
                   enddo
                elseif (kcs(nm) == 2 .and. ilyr == 1) then
                   !
                   ! At an open boundary only the composition of the transport layer
                   ! is important. If it is not valid, mark the data as dummy data:
                   ! the data will be overwritten with data coming from the neighbouring
                   ! internal point.
                   !
                   err = .false.
                   do l = 1, lsedtot
                      if (rtemp(nm, l) < 0.0_fp) err=.true.
                   enddo
                   if (err) then
                      !
                      ! dummy
                      !
                      rtemp(nm, 1) = -1.0_fp
                   endif
                else
                   !
                   ! Point/layer that will never be used: don't care about the
                   ! values. Just replace whatever was read by something valid.
                   !
                   do l = 1, lsedtot
                      rtemp(nm, l) = 0.0_fp
                   enddo
                endif
             enddo
             !
             ! Copy RTEMP data to open boundary points that have not
             ! yet been assigned valid data.
             !
             do nm = 1, nmmax
                if (kcs(nm) == 2 .and. rtemp(nm,1) < 0.0_fp) then
                   if (kcs(nm-icx) == 1) then
                      ! ndm
                      nm2 = nm - icx
                   elseif (kcs(nm+icx) == 1) then
                      ! num
                      nm2 = nm + icx
                   elseif (kcs(nm-icy) == 1) then
                      ! nmd
                      nm2 = nm - icy
                   else
                      ! nmu
                      nm2 = nm + icy
                   endif
                   do l = 1, lsedtot
                      rtemp(nm, l) = rtemp(nm2, l)
                   enddo
                endif
             enddo
             !
             ! convert sediment mass to sediment thickness
             !
             if (layertype == 'sediment mass') then
                do l = 1, lsedtot
                   do nm = 1, nmmax
                      rtemp(nm, l) = rtemp(nm, l) / cdryb(l)
                   enddo
                enddo
             endif
             !
             ! add thicknesses in lyrfrac
             !
             do l = 1, lsedtot
                do nm = 1, nmmax
                   lyrfrac(nm, ilyr, l) = lyrfrac(nm, ilyr, l) + rtemp(nm, l)
                enddo
             enddo
          else
             write (message,'(3a,i2,2a)') 'Unknown layer type ''', &
              & trim(layertype), ''' specified for layer ', ilyr, &
              & ' in file ', trim(filcomp)
             call prterr(lundia, 'U021', trim(message))
             call d3stop(1, gdp)          
          endif
          !
       enddo
       !
       deallocate(rtemp, stat = istat)
       !
       ! Process lyrfrac (sediment thicknesses) to thlyr and lyrfrac
       !
       do ilyr = 1, nlyr
          do nm = 1, nmmax
             sedthick = 0.0_fp
             do l = 1, lsedtot
                sedthick = sedthick + lyrfrac(nm, ilyr, l)
             enddo
             sedthick        = max(sedthick, 0.0_fp)
             thlyr(nm, ilyr) = sedthick
             if (sedthick > 0.0_fp) then
                do l = 1, lsedtot
                   lyrfrac(nm, ilyr, l) = lyrfrac(nm, ilyr, l) / sedthick
                enddo
             else
                do l = 1, lsedtot-1
                   lyrfrac(nm, ilyr, l) = 0.0_fp
                enddo
                lyrfrac(nm, ilyr, lsedtot) = 1.0_fp
             endif
          enddo
       enddo
       !
    else
       write (message,'(2a)') 'Invalid file version of ', trim(filcomp)
       call prterr(lundia, 'U021', trim(message))
       call d3stop(1, gdp)          
    endif
    !
end subroutine rdinimorlyr
