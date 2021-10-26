module m_readBoundaries
!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2019.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  This program is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
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
!  $Id$
!  $HeadURL$
!-------------------------------------------------------------------------------

   use m_GlobalParameters
   use MessageHandling
   use m_network
   use m_boundaryConditions

   use properties
   use m_tables
   use m_hash_search
   use string_module

   use m_ec_module

   implicit none

   private

   type(tEcInstance), pointer :: ec => NULL()

   ! needed to communicate with ec
   integer                  , dimension(:), allocatable :: ec_target_items_ids
   character(len=maxNameLen), dimension(:), allocatable :: ec_loc_names
   character(len=maxNameLen), dimension(:), allocatable :: ec_quant_names
   logical                  , dimension(:), allocatable :: ec_item_is_lateral

   ! meteo item in ec
   integer, public            :: ec_itemId_air_temperature
   integer, public            :: ec_itemId_cloudiness
   integer, public            :: ec_itemId_humidity
   integer, public            :: ec_itemId_winddirection
   integer, public            :: ec_itemId_windvelocity
   integer, public            :: ec_itemId_boundaries = ec_undef_int
   integer, public            :: ec_itemId_laterals = ec_undef_int

   integer, public, parameter :: ec_lat_vartype_disch  = 1  ! todo: move to m_laterals, and use where applicable
   integer, public, parameter :: ec_lat_vartype_sal    = 2
   integer, public, parameter :: ec_lat_vartype_temp   = 3
   integer, public, parameter :: ec_lat_num_vartypes   = 3
   
   logical                  , dimension(:)  , allocatable :: ec_item_has_been_set_externally
   integer, public          , dimension(:,:), allocatable :: ec_bnd_2_ec_index
   integer, public          , dimension(:,:), allocatable :: ec_lat_2_ec_index
   
   public ec_target_items_ids, ec_loc_names, ec_quant_names, ec_item_is_lateral, ec_item_has_been_set_externally
   public readBoundaryLocations, readBoundaryConditions, getBoundaryValue, updateBoundaryValue, closeBoundaryConditionFiles
   public getBoundaryValuePtr, getLateralValuePtr
   public disableScalarUpdate
   public write_laterals_and_boundaries

   contains

   subroutine readBoundaryLocations(network, boundaryLocationFile)

      implicit none
      
      type(t_network), target, intent(inout) :: network
      character*(*), intent(in)              :: boundaryLocationFile

      logical                                       :: success
      type(tree_data), pointer                      :: md_ptr 
      integer                                       :: istat
      integer                                       :: numstr
      integer                                       :: i
      integer                                       :: icon
      integer                                       :: ityp

      character(len=IdLen)                          :: nodeID
      integer                                       :: nodeIdx
      integer                                       :: boundaryType
      integer                                       :: icount
      double precision                              :: returnTime = 0.0d0
      
      type(t_boundarySet), pointer                  :: boundaries        !< set containing boundary conditions
      type(t_nodeSet), pointer                      :: nds               !< set containing node data

      boundaries => network%boundaries
      nds        => network%nds

      call tree_create(trim(boundaryLocationFile), md_ptr, maxlenpar)
      call prop_file('ini',trim(boundaryLocationFile),md_ptr, istat)
      
      numstr = 0
      if (associated(md_ptr%child_nodes)) then
         numstr = size(md_ptr%child_nodes)
      end if

      do i = 1, numstr
         
         if (tree_get_name(md_ptr%child_nodes(i)%node_ptr) .eq. 'boundary') then
            
            ! Read Data
            
            call prop_get_string(md_ptr%child_nodes(i)%node_ptr, 'boundary', 'nodeid', nodeID, success)
            if (success) call prop_get_integer(md_ptr%child_nodes(i)%node_ptr, 'boundary', 'type', boundaryType, success)
            if (.not. success) then
               call SetMessage(LEVEL_FATAL, 'Error Reading Boundary '''//trim(nodeID)//'''')
            endif
            
            nodeIdx = hashsearch(nds%hashlist, nodeID)
            if (NodeIdx <= 0) then
               call SetMessage(LEVEL_FATAL, 'Error Reading Boundary '''//trim(nodeID)//''': Node Not Found')
            endif

            if (.not. doReadCache) then  ! When Reading Cache This Already has been Checked and Set.

               if (nds%count > 0) then
                  if ((boundaryType == 1 .or. boundaryType == 2)  .and. &
                      (nds%node(nodeIdx)%nodetype == nt_LevelBoun .or. nds%node(nodeIdx)%nodetype == nt_DischBoun)) then

                     ! two boundaries on one node
                     call setMessage(LEVEL_FATAL, 'Error Reading Boundary: On node '//trim(nds%node(nodeIdx)%id) //' more than one boundary condition is defined.')

                  elseif (boundaryType == 2 .and. nds%node(nodeIdx)%nodetype /= nt_EndNode) then
                     call setMessage(LEVEL_FATAL, 'Error Reading Boundary: Discharge Boundary Node '//trim(nds%node(nodeIdx)%id) //' is not an end node.')
                     return
                  endif
               endif

               if (nds%count > 0) then
                  if (boundaryType==1) then
                     ! water level boundary
                     nds%node(nodeIdx)%nodetype = nt_LevelBoun
                  elseif (boundaryType==2) then
                     ! discharge boundary
                     nds%node(nodeIdx)%nodetype = nt_DischBoun
                  endif

               endif
               
            endif

            boundaries%tp(boundaryType)%Count = boundaries%tp(boundaryType)%Count + 1
            icount = boundaries%tp(boundaryType)%Count
            if (icount > boundaries%tp(boundaryType)%Size) then
               call realloc(boundaries%tp(boundaryType), boundaryType)
            endif

            boundaries%tp(boundaryType)%bd(icount)%nodeID        = nodeID
            boundaries%tp(boundaryType)%bd(icount)%node          = nodeIdx
            boundaries%tp(boundaryType)%bd(icount)%salinityIndex = -1
            boundaries%tp(boundaryType)%bd(icount)%iopt          = 1
            boundaries%tp(boundaryType)%bd(icount)%length        = 2
            boundaries%tp(boundaryType)%bd(icount)%boundaryValue = 0d0


            do icon = 1, transportPars%constituents_count
               
               ! For each constituent add a boundary condition
               ityp = transportPars%co_h(icon)%boundary_index
               boundaries%tp(ityp)%Count = boundaries%tp(ityp)%Count + 1
               if (boundaries%tp(ityp)%Count > boundaries%tp(ityp)%Size) then
                  call realloc(boundaries%tp(ityp), ityp)
               endif
         
               ! Get Thatcher-Harleman Return Time
               call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'boundary', 'thatcher-harlemancoeff', returnTime, success)
               if (success) then 
                  if (returnTime < 0.0d0) then
                     call setMessage(LEVEL_ERROR, 'Negative Thatcher-Harleman Coefficient at Boundary Node '//trim(nds%node(nodeIdx)%id))
                     returnTime = 0.0d0
                  endif
               else
                  returnTime = 0.0d0
               endif
               
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%nodeID        = nodeID
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%node          = nodeIdx
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%salinityIndex = -1
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%iopt          = 2     ! Thatcher-Harleman
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%length        = 2
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%returnTime    = returnTime
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%boundaryValue = 0d0
               
               boundaries%tp(boundaryType)%bd(icount)%salinityIndex = boundaries%tp(ityp)%Count
               
               boundaries%tp(ityp)%bd(boundaries%tp(ityp)%Count)%returnTime = returnTime
         
            enddo
         endif

      end do
      
      boundaries%count = 0
      do i = 1, NUM_BOUN_TYPE
         boundaries%count = boundaries%count + boundaries%tp(i)%Count
         call fill_hashtable(network%boundaries%tp(i))
      enddo
      
      call tree_destroy(md_ptr)
      
   end subroutine readBoundaryLocations


function getBoundaryValuePtr(network, quantityId, locationId) result (valuePtr)
   use m_boundaryConditions
   double precision, pointer               :: valuePtr
   type(t_network) , intent(inout), target :: network
   character(len=*), intent(in)            :: quantityID
   character(len=*), intent(in)            :: locationID
   integer :: locIndex, varType
   valuePtr => null()
   select case (quantityId)
      case ('wind_speed')
         varType = B_WINDVEL
      case ('wind_from_direction')
         varType = B_WINDDIR
      case('water_level')
         varType = H_BOUN
      case('water_discharge')
         varType = Q_BOUN
      case('water_salinity')
         varType = S_BOUN
      case('water_temperature')
         varType = T_BOUN
      case default
         call SetMessage(LEVEL_ERROR, 'Quantity type ' // trim(quantityId) // ' not yet implemented in set_external_value')
         return
   end select
   if (trim(locationID) == 'model_wide') then
      locIndex = 1
   else
      locIndex = hashsearch(network%boundaries%tp(varType)%hashlist, locationID)
   endif
   valuePtr => network%boundaries%tp(varType)%bd(locIndex)%boundaryValue
end function getBoundaryValuePtr


function getLateralValuePtr(network, quantityId, locationId) result (valuePtr)
   use m_boundaryConditions
   double precision, pointer               :: valuePtr
   type(t_network) , intent(inout), target :: network
   character(len=*), intent(in)            :: quantityID
   character(len=*), intent(in)            :: locationID
   integer :: locIndex, varType
   locIndex = hashsearch(network%lts%hashlist, locationID)
   valuePtr => null()
   select case(quantityID)
      case('water_discharge')
         varType = ec_lat_vartype_disch
         valuePtr => network%lts%lat(locIndex)%discharge
      case('water_salinity')
         varType = ec_lat_vartype_sal
         valuePtr => network%lts%lat(locIndex)%concentration(transportPars%salt_index)
      case('water_temperature')
         varType = ec_lat_vartype_temp
         valuePtr => network%lts%lat(locIndex)%concentration(transportPars%temp_index)
      case default
         call SetMessage(LEVEL_ERROR, 'Quantity type ' // trim(quantityId) // ' not yet implemented in set_external_value')
         return
   end select         
end function getLateralValuePtr
   

subroutine readBoundaryConditions(network, boundaryConditionsFile)

   use m_ec_bccollect
   use m_globalParameters

   implicit none
      
   type(t_network) ,  target, intent(inout) :: network                 ! the model's network
   character(len=*),          intent(in)    :: boundaryConditionsFile  ! file containing the lateral conditions

   logical              :: success        ! return value of function call
   integer              :: istat          ! status after function call
   integer              :: i, inod, ityp  ! loop counter, node index, boundary type
   integer              :: total_items    ! total number of items expected in bc-file (boundaries, laterals)

   character(len=IdLen) :: locationID     ! boundary location identifier (node name)
   character(len=IdLen) :: quantityID     ! water level or discharge
   
   integer              :: bc_count       ! counting bc-elements other than boundaries and laterals 
   integer              :: ec_src_item    ! external forcing source item id in ec-module
   integer              :: ec_tgt_item    ! external forcing target item id in ec-module (i.e. post temporal interpolation)

   logical                         :: is_qh_bound         ! Q(h) boundary?
   real(hp), dimension(:), pointer :: h_values, q_values  ! Q(h) table
   integer                         :: numValues           ! #rows in Q(h) table
   type(t_tp), pointer             :: p_bnds              ! pointer to boundaries (q-bnds, h-bnds, etc.)
   integer     :: numpars
   double precision, pointer       :: valueptr            ! dummy storage pointer to boundary value passed to ec
   integer                         :: connectionId        ! temporarily stored ID of newly created ec-module connection instance

   inquire(file=boundaryConditionsFile, exist=success)

   if (.not. success) then
      call setmessage(level_error, trim(boundaryConditionsFile) // ' not found.')
      return
   endif

   if (associated(ec)) then
      success = ecInstanceFree(ec)
      ec => null()
      if (.not. success) then
         call setmessage(LEVEL_FATAL, 'Reading Boundaries: Could not free/deallocate ec-instance')
      endif
   endif
   
   if (.not. ecInstanceCreate(ec))then
      call setmessage(LEVEL_FATAL, 'Reading Boundaries: Could not create ec-instance')
   endif

   bc_count = ecCollectBCBlocks(ec, boundaryConditionsFile, istat)

   if (bc_count /= 0 .and. istat /= 0)then ! boundaries/laterals found, but error in reading them
      call setmessage(LEVEL_ERROR, 'Reading Boundaries: could not read boundary blocks from file '// trim(boundaryConditionsFile))
      call setmessage(LEVEL_ERROR, dumpECMessageStack(LEVEL_ERROR, setmessage))
      return
   endif
   
   istat = 0
   total_items = network%lts%count + network%boundaries%count
   
   if (bc_count < total_items) then
      call SetMessage(LEVEL_ERROR, 'Reading Boundaries: Not enough boundaries in file ' // trim(boundaryConditionsFile))
	else
      total_items = bc_count   
   endif

   call realloc(ec_target_items_ids, total_items, stat=istat)
   if (istat == 0) call realloc(ec_loc_names, total_items, stat=istat)
   if (istat == 0) call realloc(ec_quant_names, total_items, stat=istat)
   if (istat == 0) call realloc(ec_item_is_lateral, total_items, stat=istat)
   if (istat == 0) call realloc(ec_item_is_lateral, total_items, stat=istat)
   if (istat == 0) call realloc(ec_item_has_been_set_externally, total_items, stat=istat)
   if (istat == 0) call realloc(ec_bnd_2_ec_index, NUM_BOUN_TYPE, network%boundaries%count, stat=istat)
   if (istat == 0) call realloc(ec_lat_2_ec_index, ec_lat_num_vartypes, network%lts%count, stat=istat)
   if (istat .ne. 0) then
      call SetMessage(LEVEL_FATAL, 'Reading Boundaries: Error Allocating Arrays')
   endif
   
   ec_target_items_ids = -1
   ec_loc_names = ' '
   ec_quant_names = ' '
   ec_item_is_lateral = .false.
   ec_item_has_been_set_externally = .false.
   ec_bnd_2_ec_index = -1
   ec_lat_2_ec_index = -1

   bc_count = 0
   ! laterals
   ec_itemId_laterals = ecInstanceCreateItem(ec)
   if (.not.ecSetItemRole(ec, ec_itemId_laterals, itemType_target)) then
       return
   endif
   do i = 1, network%lts%count
      locationID = network%lts%lat(i)%id
      numpars = 1 + transportPars%constituents_count
      
      do ityp = 1, numpars
         if (ityp==1) then
            quantityID = 'water_discharge'
         else
            quantityID = trim(transportPars%co_h(ityp-1)%name)
         endif
         
         ec_src_item = ecFindItemByQuantityLocation(ec, locationID, quantityID, isLateral = .true.)
         if (ec_src_item < 0)then
            ! use default concentration for this constituent
            continue
         else
            is_qh_bound = ecItemGetQHtable(ec, ec_src_item, h_values, q_values, success)
            if (is_qh_bound) then
               if (.not.success) then
                  call setmessage(LEVEL_FATAL, 'Error getting lateral discharge qh-table for ' &
                                               // trim(locationID) // '.' // trim(quantityID) )
               else
                  call SetQlatQh(network%lts%lat(i), h_values, q_values, size(h_values))
               endif
            else
               valueptr => getLateralValuePtr(network, quantityId, locationId) 
               ec_tgt_item = ecCreateTimeInterpolatedItem(ec, ec_src_item, scalarPtr=valueptr)
               ec_lat_2_ec_index(ityp,i) = ec_tgt_item
               connectionId = ecCreateConnection(ec)
               if (.not.ecAddConnectionSourceItem(ec, connectionId, ec_tgt_item)) return
               if (.not.ecAddConnectionTargetItem(ec, connectionId, ec_itemId_laterals)) return
               if (.not.ecAddItemConnection(ec, ec_itemId_laterals, connectionId)) return
            endif
         endif
      enddo
   enddo

   ! boundaries
   ec_itemId_boundaries = ecInstanceCreateItem(ec)
   if (.not.ecSetItemRole(ec, ec_itemId_boundaries, itemType_target)) then
       return
   endif
   do ityp = 1, NUM_BOUN_TYPE
      p_bnds => network%boundaries%tp(ityp)
      do i = 1, p_bnds%count
         inod = p_bnds%bd(i)%node
         locationID = network%nds%node(inod)%id
         if (ityp==H_BOUN) then
            quantityID = 'water_level'
         else if (ityp==Q_BOUN) then
            quantityID = 'water_discharge'
         else if (ityp==S_BOUN) then
            quantityID = 'water_salinity'
         else if (ityp==T_BOUN) then
            quantityID = 'water_temperature'
         else
            call setmessage(LEVEL_ERROR, 'Boundary type not supported')
            cycle
         endif
         ec_src_item = ecFindItemByQuantityLocation(ec, locationID, quantityID)
         if (ec_src_item < 0)then
            call setmessage(LEVEL_ERROR, 'Could not find ' // trim(locationID) // '.' // trim(quantityID) // ' in file ' // trim(boundaryConditionsFile))
         else
            is_qh_bound = .false.
            if ((ityp==Q_BOUN) .or. (ityp==H_BOUN)) then
               is_qh_bound = ecItemGetQHtable(ec, ec_src_item, h_values, q_values, success)
               if (is_qh_bound) then
                  if (.not.success) then
                     call setmessage(LEVEL_FATAL, 'Error getting boundary condition qh-table for ' &
                                            // trim(locationID) // '.' // trim(quantityID) )
                  endif
               endif
            endif
            if (is_qh_bound) then
               numValues = size(h_values)
               p_bnds%bd(i)%length = numValues
               p_bnds%bd(i)%iopt   = 2
               if (ityp==h_boun) then
                  call setTable(p_bnds%bd(i)%table, 0, q_values, h_values, numValues)
               else
                  call setTable(p_bnds%bd(i)%table, 0, h_values, q_values, numValues)
               endif
               p_bnds%bd(i)%table%interpoltype = 0
               
            else
               valueptr => getBoundaryValuePtr(network, quantityId, locationId) 
               ec_tgt_item = ecCreateTimeInterpolatedItem(ec, ec_src_item, scalarPtr=valueptr)
               ec_bnd_2_ec_index(ityp,i) = ec_tgt_item
               connectionId = ecCreateConnection(ec)
               if (.not.ecAddConnectionSourceItem(ec, connectionId, ec_tgt_item)) return
               if (.not.ecAddConnectionTargetItem(ec, connectionId, ec_itemId_boundaries)) return
               if (.not.ecAddItemConnection(ec, ec_itemId_boundaries, connectionId)) return
            endif
         endif
         !TODO invullen: function type for interpolation\n
         ! - div(interpoltype, 10) == 1: periodical function (nvt)
         ! - mod(interpoltype, 10) == 0: linear function
         ! - mod(interpoltype, 10) == 1: block function
         ! - mod(interpoltype, 10) == 2: linear function in degrees
      enddo
   enddo

   ! Meteo
   ec_src_item = ecFindItemByQuantityLocation(ec, 'model_wide', 'wind_speed')
   if (ec_src_item > 0) then
      ec_itemId_windvelocity  = ecCreateTimeInterpolatedItem(ec, ec_src_item)
   else
      ec_itemId_windvelocity  = 0
   endif
   
   ec_src_item = ecFindItemByQuantityLocation(ec, 'model_wide', 'wind_from_direction')
   if (ec_src_item > 0) then
      ec_itemId_winddirection  = ecCreateTimeInterpolatedItem(ec, ec_src_item)
   else
      ec_itemId_winddirection  = 0
   endif

   ec_src_item = ecFindItemByQuantityLocation(ec, 'model_wide', 'humidity')
   if (ec_src_item > 0) then
      ec_itemId_humidity  = ecCreateTimeInterpolatedItem(ec, ec_src_item)
   else
      ec_itemId_humidity  = 0
   endif

   ec_src_item = ecFindItemByQuantityLocation(ec, 'model_wide', 'air_temperature')
   if (ec_src_item > 0) then
      ec_itemId_air_temperature  = ecCreateTimeInterpolatedItem(ec, ec_src_item)
   else
      ec_itemId_air_temperature  = 0
   endif

   ec_src_item = ecFindItemByQuantityLocation(ec, 'model_wide', 'cloudiness')
   if (ec_src_item > 0) then
      ec_itemId_cloudiness  = ecCreateTimeInterpolatedItem(ec, ec_src_item)
   else
      ec_itemId_cloudiness  = 0
   endif
   
   ! Prevent crash in Model-API
   ! TODO: Remove this when not necessary anymore
   allocate(network%boundaries%tp(B_WINDVEL)%bd(1))
   network%boundaries%tp(B_WINDVEL)%bd(1)%boundaryValue = 0d0
   allocate(network%boundaries%tp(B_WINDDIR)%bd(1))
   network%boundaries%tp(B_WINDDIR)%bd(1)%boundaryValue = 0d0

   ! Dump EC item hierarchy
!  call ecInstancePrintState(ec,callback_msg,LEVEL_DEBUG)
   call ecInstancePrintState(ec,callback_msg,LEVEL_INFO)

end subroutine readBoundaryConditions

subroutine callback_msg(lvl,msg)
   integer, intent(in)              :: lvl
    character(len=*), intent(in)    :: msg 
    call setmessage(lvl, trim(msg))
end subroutine

function getBoundaryValue(ec_target_item, timeAsMJD) result(value_from_ec)
   double precision               :: value_from_ec
   integer         , intent(in)   :: ec_target_item
   double precision, intent(in)   :: timeAsMJD
   double precision, dimension(1) :: array_values_from_ec
   
   value_from_ec = 0.0
   if (.not. ecGetValues(ec, ec_target_item, timeAsMJD, array_values_from_ec) ) then
      call SetMessage(LEVEL_FATAL, 'Error ec_target_item value from EC file')
   else
      value_from_ec = array_values_from_ec(1)
   endif

end function getBoundaryValue

subroutine updateBoundaryValue(ec_target_item, timeAsMJD) 
   integer         , intent(in)   :: ec_target_item
   double precision, intent(in)   :: timeAsMJD
   
   if (.not. ecGetValues(ec, ec_target_item, timeAsMJD) ) then
      call SetMessage(LEVEL_FATAL, 'Error ec_target_item value from EC file')
   endif

end subroutine updateBoundaryValue

subroutine disableScalarUpdate(ec_target_item) 
   integer         , intent(in)   :: ec_target_item
   type(tEcItem), pointer         :: itemPtr
   
   itemPtr => ecSupportFindItem(ec, ec_target_item)
   if (associated(itemPtr)) then
      if (associated(itemPtr%targetFieldPtr)) then
         itemPtr%targetFieldPtr%scalarptr => null()
      endif
   endif
end subroutine disableScalarUpdate

subroutine closeBoundaryConditionFiles()
   logical :: success
   if (associated(ec)) then
      success = ecInstanceFree(ec)
      ec => null()
   endif
end subroutine closeBoundaryConditionFiles

subroutine write_laterals_and_boundaries(network, timeAsMJD, f_lat, f_bnd)
   type(t_network), target, intent(inout) :: network
   double precision, intent(in)   :: timeAsMJD
   character(len=*), intent(in) :: f_lat         ! file to write the lateral values
   character(len=*), intent(in) :: f_bnd         ! file to write the boundaries values

   integer :: devlat     ! device to write the lateral values
   integer :: devbnd     ! device to write the boundary values
   open(file=f_lat, newunit=devlat, access='append')
   open(file=f_bnd, newunit=devbnd, access='append')
   call write_laterals_and_boundaries_stream(network, timeAsMJD, devlat, devbnd)
   close(devbnd)
   close(devlat)
end subroutine write_laterals_and_boundaries

subroutine write_laterals_and_boundaries_stream(network, timeAsMJD, devlat, devbnd)
   type(t_network), target, intent(inout) :: network
   double precision, intent(in)   :: timeAsMJD
   integer, intent(in) :: devlat     ! device to write the lateral values
   integer, intent(in) :: devbnd     ! device to write the boundary values
   integer :: nType, nLoc, locIndex, varType
   ! write laterals for a single timestep
   nLoc = size(network%lts%lat)
   write(devlat,'(f15.5,a$)') timeAsMJD, '|'
   do locIndex = 1, nLoc
      associate (lateral => network%lts%lat(locIndex))
         write(devlat,'(i5.5,a,3e12.4,a$)') locIndex, ':', &
                                         lateral%discharge, &
                                         lateral%concentration(transportPars%temp_index), &
                                         lateral%concentration(transportPars%salt_index), '|'
      end associate ! lateral
   enddo
   write(devlat,*)
   
   ! write boundaries for a sngle timestep
   write(devbnd,'(f15.5,a$)') timeAsMJD, '|'
   nType = size(network%boundaries%tp)
   do varType = 1, nType  
      write(devbnd,'(i5.5,a$)') varType, ':'
      nLoc = size(network%boundaries%tp(varType)%bd)
      do locIndex = 1, nLoc
         write(devbnd,'(e12.4,a$)') network%boundaries%tp(varType)%bd(locIndex)%boundaryValue, ' '
      enddo
      write(devbnd,'(a$)') '|'
   enddo
   write(devbnd,*)
end subroutine write_laterals_and_boundaries_stream
 

end module m_readBoundaries
