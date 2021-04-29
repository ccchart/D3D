module m_1d_networkreader
!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
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

   use MessageHandling
   use properties
   use m_hash_search
   use m_hash_list
   use m_network
   
   implicit none

   private
   
   public NetworkReader
   public NetworkUgridReader
   public read_1d_ugrid
   public read_network_cache
   public write_network_cache
   public construct_network_from_meshgeom

   contains
    
   subroutine NetworkReader(network, networkFile)
   
      use m_hash_search
      
      implicit none
   
      type(t_network), target, intent(inout) :: network
      character(len=*), intent(in) :: networkFile
      
      type(tree_data), pointer  :: md_ptr
      integer :: istat
      integer :: numstr
      integer :: i

      call tree_create(trim(networkFile), md_ptr, maxlenpar)
      call prop_file('ini',trim(networkFile), md_ptr, istat)

      numstr = 0
      if (associated(md_ptr%child_nodes)) then
         numstr = size(md_ptr%child_nodes)
      end if

      ! Get the Nodes First
      do i = 1, numstr
         
         if (tree_get_name(md_ptr%child_nodes(i)%node_ptr) .eq. 'node') then
           call readNode(network%nds, md_ptr%child_nodes(i)%node_ptr)
         endif

      enddo
      call fill_hashtable(network%nds)
      
      ! Get the Branches
      do i = 1, numstr
         
         if (tree_get_name(md_ptr%child_nodes(i)%node_ptr) .eq. 'branch') then
           call readBranch(network%brs, network%nds, md_ptr%child_nodes(i)%node_ptr)
         endif

      enddo
      
      call adminBranchOrders(network%brs)
      call fill_hashtable(network%brs)
      
      call tree_destroy(md_ptr)
      
   end subroutine NetworkReader
   
   subroutine NetworkUgridReader(network, networkUgridFile)
   
      use io_netcdf
      use io_ugrid
      use m_hash_search
      use gridgeom
      use meshdata
      
      implicit none
   
      type(t_network), target, intent(inout) :: network
      character(len=*), intent(in)           :: networkUgridFile
      
      integer                   :: ierr
      integer                   :: ioncid
      
      ! Open UGRID-File
      ierr = ionc_open(networkUgridFile, NF90_NOWRITE, ioncid)
      if (ierr .ne. 0) then
         call SetMessage(LEVEL_FATAL, 'Error Opening UGRID-File: '''//trim(networkUgridFile)//'''')
      endif

      ! Do the actual read
      call read_1d_ugrid(network, ioncid)
      
      ! Close UGRID-File
      ierr = ionc_close(ioncid)
      if (ierr .ne. 0) then
         call SetMessage(LEVEL_FATAL, 'Error Closing UGRID-File: '''//trim(networkUgridFile)//'''')
      endif
       
   end subroutine NetworkUgridReader  


   !> Constructs a flow1d t_network datastructure, based on meshgeom read from a 1D UGRID file.
   !! meshgeom is only used for reading a UGRID NetCDF file, whereas network is used during
   !! a model computation.
   !! The network data structure assumes to have grid points on all start and end points of branches,
   !! but the input meshgeom often will only have one unique grid point on a connection node.
   !! In that case, parameter nodesOnBranchVertices allows to automatically create duplicate start/end points.
   integer function construct_network_from_meshgeom(network, meshgeom, branchids, branchlongnames, nodeids, nodelongnames, & !1d network character variables
      gpsID, gpsIDLongnames, network1dname, mesh1dname, nodesOnBranchVertices, jampi, my_rank) result(ierr)

   use gridgeom
   use meshdata
   use m_hash_search
   use odugrid
   use io_netcdf, only : ionc_get_node_coordinates

   !in variables
   type(t_network),                             intent(inout) :: network
   type(t_ug_meshgeom),                         intent(in   ) :: meshgeom
   character(len=*), allocatable, dimension(:), intent(in   ) :: branchids
   character(len=*), allocatable, dimension(:), intent(in   ) :: branchlongnames
   character(len=*), allocatable, dimension(:), intent(in   ) :: nodeids
   character(len=*), allocatable, dimension(:), intent(in   ) :: nodelongnames
   character(len=*), allocatable, dimension(:), intent(inout) :: gpsID
   character(len=*), allocatable, dimension(:), intent(inout) :: gpsIDLongnames
   character(len=*),                            intent(in   ) :: network1dname
   character(len=*),                            intent(in   ) :: mesh1dname
   integer,                                     intent(in   ) :: nodesOnBranchVertices !< Whether or not (1/0) the input meshgeom itself already contains duplicate points on each connection node between multiple branches.
                                                                                       !! If not (0), additional grid points will be created.
   integer,                                     intent(in   ), optional :: jampi       !< running in parallel mode (1) or not (0)
   integer,                                     intent(in   ), optional :: my_rank     !< my rank in parallel mode, for (debugging) output

   !locals
   integer, allocatable, dimension(:)               :: gpFirst
   integer, allocatable, dimension(:)               :: gpLast
   double precision, allocatable, dimension(:)      :: gpsX
   double precision, allocatable, dimension(:)      :: gpsY
   integer                                          :: ibran, inode, jsferic
   integer                                          :: gridPointsCount
   double precision, allocatable, dimension(:)      :: localOffsets
   double precision, allocatable, dimension(:)      :: localOffsetsSorted
   integer, allocatable, dimension(:)               :: localSortedIndexses
   double precision, allocatable, dimension(:)      :: localGpsX
   double precision, allocatable, dimension(:)      :: localGpsY
   character(len=len(gpsId)), allocatable, dimension(:)  :: localGpsID
   character(len=len(gpsId)), allocatable, dimension(:)  :: idMeshNodesInNetworkNodes
   integer                                          :: firstNode, lastNode
   double precision, parameter                      :: snapping_tolerance = 1e-10
   double precision                                 :: distance, meanLength
   integer                                          :: jampi_, my_rank_
   integer                                          :: maxGridPointCount
   logical, allocatable                             :: active_nodes(:)
   integer, allocatable                             :: active_branches(:)

   ierr = -1

   jampi_ = 0
   if (present(jampi)) jampi_ = jampi
   my_rank_ = -1
   if (present(my_rank)) my_rank_ = my_rank

   ! check data are present and correct
   if (meshgeom%numnode .eq. -1) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%numnode')
   endif
   if (meshgeom%nbranches .eq. -1) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nbranches')
   endif
   if (.not.associated(meshgeom%nodex)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nodex')
   endif
   if (.not.associated(meshgeom%nodey)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nodey')
   endif
   if (.not.allocated(nodeids)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in nodeids')
   endif
   if (.not.allocated(nodelongnames)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in nodelongnames')
   endif

   if (.not.allocated(branchids)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nbranchids')
   endif
   if (.not.associated(meshgeom%nedge_nodes)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nedge_nodes')
   endif
   if (.not.associated(meshgeom%nbranchorder)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nbranchorder')
   endif
   if (.not.associated(meshgeom%nodex)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nodex')
   endif
   if (.not.associated(meshgeom%nodey)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nodey')
   endif
   if (.not.associated(meshgeom%nodeoffsets)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in meshgeom%nodeoffsets')
   endif
   if (.not.allocated(gpsID)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in gpsID')
   endif
   if (.not.allocated(gpsIDLongnames)) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error in gpsIDLongnames')
   endif

   ! Store Network Node Data ('connection nodes') into Data Structures.
   call storeNodes(network%nds, meshgeom%nnodes, meshgeom%nnodex, meshgeom%nnodey, nodeids, nodelongnames)

   ! Calculate mesh1d x,y coordinates (computational grid points), based on UGRID branchId-based notation.
   allocate(gpsX(meshgeom%numnode), stat = ierr)
   if (ierr == 0) allocate(gpsY(meshgeom%numnode), stat = ierr)
   if (ierr .ne. 0) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Allocating Memory for Grid Point Coordinates')
   endif
   
   if (meshgeom%epsg==4326) then ! TODO: UNST-2510: LC: %epsg is never set (unless via c_meshgeom?), so this code does not work for lat/lon 1D networks.
      jsferic = 1
   else
      jsferic = 0
   endif
   if (jampi_ == 0) then
      ierr = ggeo_get_xy_coordinates(meshgeom%nodebranchidx, meshgeom%nodeoffsets, meshgeom%ngeopointx, meshgeom%ngeopointy, &
         meshgeom%nbranchgeometrynodes, meshgeom%nbranchlengths, jsferic, gpsX, gpsY)
   else
      ierr = ionc_get_node_coordinates(1, 1, gpsX, gpsY)
   end if
   if (ierr /= 0) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Getting Mesh Coordinates From UGrid Data')
   endif

   ! Get the starting and ending mesh1d grid point indexes for each network branch.
   ibran = 0
   allocate(gpFirst(meshgeom%nbranches), stat = ierr)
   if (ierr == 0) allocate(gpLast(meshgeom%nbranches), stat = ierr)
   if (ierr /= 0) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Allocating Memory for Branches')
   endif

   ierr = ggeo_get_start_end_nodes_of_branches(meshgeom%nodebranchidx, gpFirst, gpLast)
   if (jampi_ /= 0) then
      call active_network_nodes_branches(meshgeom, active_nodes, active_branches)
   else
      allocate(active_branches(meshgeom%nbranches))
      active_branches = 1
   end if
   if (ierr /= 0) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Getting first and last nodes of the network branches')
   endif

   ! Fill the array storing the mesh1d node ids for each network node.
   if(nodesOnBranchVertices==0) then
      allocate(idMeshNodesInNetworkNodes(meshgeom%nnodes))
      idMeshNodesInNetworkNodes = ' '
      do ibran = 1, meshgeom%nbranches
         firstNode = gpFirst(ibran)
         lastNode  = gpLast(ibran)
         ! if no mesh points in the branch, cycle
         if(firstNode < 0 .or. lastNode < 0) then
            cycle
         endif
         do inode = firstNode, lastNode
            if(meshgeom%nodeoffsets(inode)<snapping_tolerance) then
               idMeshNodesInNetworkNodes(meshgeom%nedge_nodes(1,ibran))(1:len_trim(gpsID(inode))) = gpsID(inode)(1:len_trim(gpsID(inode)))
            endif
            if(abs(meshgeom%nodeoffsets(inode)-meshgeom%nbranchlengths(ibran))<snapping_tolerance) then
               idMeshNodesInNetworkNodes(meshgeom%nedge_nodes(2,ibran))(1:len_trim(gpsID(inode))) = gpsID(inode)(1:len_trim(gpsID(inode)))
            endif
         enddo
      enddo
   endif

   ! allocate local arrays
   maxGridPointCount = 3 + maxval(gpLast - gpFirst)
   allocate(localOffsets(maxGridPointCount))
   allocate(localOffsetsSorted(maxGridPointCount))
   allocate(localSortedIndexses(maxGridPointCount))
   allocate(localGpsX(maxGridPointCount))
   allocate(localGpsY(maxGridPointCount))
   allocate(localGpsID(maxGridPointCount))

   ! Store the branches + computational grid points on them into Data Structures.
   do ibran = 1, meshgeom%nbranches

      firstNode = gpFirst(ibran)
      lastNode  = gpLast(ibran)
      ! if no mesh points in the branch, cycle
      if(firstNode < 0 .or. lastNode < 0) then
         ! end node and begin node are missing and no internal gridpoints on branch.
         localOffsets = 0d0
         gridpointscount = 0
         ! set dummy local offset to half the branch length
         localGpsX   = 0d0
         localGpsY   = 0d0
         localGpsID  = ''
      else
         localOffsets = 0d0
         gridPointsCount                 = lastNode - firstNode + 1
         localOffsets(1:gridPointsCount) = meshgeom%nodeoffsets(firstNode:lastNode)
         localGpsX(1:gridPointsCount)    = gpsX(firstNode:lastNode)
         localGpsY(1:gridPointsCount)    = gpsY(firstNode:lastNode)
         localGpsID(1:gridPointsCount)   = gpsID(firstNode:lastNode)
      endif

      if (nodesOnBranchVertices==0 .and. (jampi_ == 0 .or. active_branches(ibran) == 1)) then
         if(localOffsets(1)>snapping_tolerance .or. gridpointsCount == 0) then
            !start point missing
            call add_point(.true., localOffsets, localGpsX, localGpsY, localGpsID, idMeshNodesInNetworkNodes, gridPointsCount, ibran, meshgeom)
         endif
         ! TODO: consider using a relative tolerance
         if(abs(localOffsets(gridPointsCount)-meshgeom%nbranchlengths(ibran))> snapping_tolerance .or. gridpointsCount == 1) then
            !end point missing
            call add_point(.false., localOffsets, localGpsX, localGpsY, localGpsID, idMeshNodesInNetworkNodes, gridPointsCount, ibran, meshgeom)
         endif
      else if(nodesOnBranchVertices==0 .and. jampi_ == 1 .and. active_branches(ibran) == -1) then
         if (firstNode > 0 .and. lastNode > 0) then
            if (gridPointsCount > 1) then
               meanLength = (localOffsets(gridPointsCount) - localOffsets(1)) / dble(gridPointsCount - 1)
            else
               meanLength = meshgeom%nbranchlengths(ibran)
            end if
            distance = localOffsets(1)
            if (distance > snapping_tolerance .and. distance < 2d0 * meanLength) then
               !start point missing
               call add_point(.true., localOffsets, localGpsX, localGpsY, localGpsID, idMeshNodesInNetworkNodes, gridPointsCount, ibran, meshgeom)
            endif
            distance = abs(localOffsets(gridPointsCount)-meshgeom%nbranchlengths(ibran))
            if (distance > snapping_tolerance .and. distance < 2d0 * meanLength) then
               !end point missing
               call add_point(.false., localOffsets, localGpsX, localGpsY, localGpsID, idMeshNodesInNetworkNodes, gridPointsCount, ibran, meshgeom)
            endif
         endif
      endif

      call storeBranch(network%brs, network%nds, branchids(ibran), nodeids(meshgeom%nedge_nodes(1,ibran)), &
         nodeids(meshgeom%nedge_nodes(2,ibran)), meshgeom%nbranchorder(ibran), gridPointsCount, localGpsX(1:gridPointsCount), &
         localGpsY(1:gridPointsCount),localOffsets(1:gridPointsCount), localGpsID(1:gridPointsCount), active_branches(ibran), my_rank_)
   enddo

   call adminBranchOrders(network%brs)
   call fill_hashtable(network%brs)

   network%loaded = .true.

   !free local memory
   deallocate(gpsX)
   deallocate(gpsY)

   end function construct_network_from_meshgeom

   !> determine which network_nodes and network_branches are active in the current domain.
   !! TODO: move to partitioner or part of t_ug_meshgeom
   subroutine active_network_nodes_branches(meshgeom, active_nodes, active_branches)
   use meshdata, only : t_ug_meshgeom
   type(t_ug_meshgeom), intent(in   ) :: meshgeom          !< struct for mesh geometry
   logical, allocatable, intent(out) :: active_nodes(:)    !< true if node is in current domain
   integer, allocatable, intent(out) :: active_branches(:) !< 0 if outside current domain; 1 if completely in current domain; -1 if partly in current domain

   integer :: ibran, istart, istop
   logical, allocatable :: active_b(:)  ! array for simple check if branch is at least partly in current domain

   if (allocated(active_nodes))    deallocate(active_nodes)
   if (allocated(active_branches)) deallocate(active_branches)

   allocate(active_nodes(meshgeom%nnodes))
   allocate(active_branches(meshgeom%nbranches))

   allocate(active_b(meshgeom%nbranches))
   do ibran = 1,  meshgeom%nbranches
      active_b(ibran) = (any(meshgeom%edgebranchidx(:) == ibran))
   enddo

   active_nodes(:) = .false.
   do ibran = 1,  meshgeom%nbranches
      if (.not. active_b(ibran)) then
         active_branches(ibran) = 0  ! outside current domain
         cycle
      end if
      istart = meshgeom%nedge_nodes(1, ibran)
      istop  = meshgeom%nedge_nodes(2, ibran)
      if (.not. active_nodes(istart)) then
         ! behind if to minimize calls to coordinate_check
         active_nodes(istart) = coordinate_check(istart)
      end if
      if (.not. active_nodes(istop)) then
         active_nodes(istop)  = coordinate_check(istop)
      end if

      if (active_nodes(istart) .and. active_nodes(istop)) then
         active_branches(ibran) = 1   ! completely in current domain
      else if (active_b(ibran)) then
         active_branches(ibran) = -1  ! partly in current domain
      else
         active_branches(ibran) = 0   ! outside current domain
      end if
   end do

   contains

   function coordinate_check(icmp) result (found)
      integer, intent(in) :: icmp
      logical             :: found

      integer             :: i

      found = .false.
      do i = 1, meshgeom%numnode
         if (meshgeom%nnodex(icmp) == meshgeom%nodex(i) .and. meshgeom%nnodey(icmp) == meshgeom%nodey(i)) then
            found = .true.
            exit
         endif
      end do
   end function coordinate_check

   end subroutine active_network_nodes_branches

   !> helper function to add a point at the start or end of a branch
   subroutine add_point(atStart, localOffsets, localGpsX, localGpsY, localGpsID, idMeshNodesInNetworkNodes, gridPointsCount, ibran, meshgeom)
   use meshdata, only : t_ug_meshgeom
   logical,                                         intent(in   ) :: atStart     !< add point at start (true) or end (false) of branch
   double precision,     allocatable, dimension(:), intent(inout) :: localGpsX
   double precision,     allocatable, dimension(:), intent(inout) :: localGpsY
   double precision,     allocatable, dimension(:), intent(inout) :: localOffsets
   character(len=*),     allocatable, dimension(:), intent(inout) :: localGpsID
   character(len=*),     allocatable, dimension(:), intent(in   ) :: idMeshNodesInNetworkNodes
   integer,                                         intent(inout) :: gridPointsCount
   integer,                                         intent(in   ) :: ibran
   type(t_ug_meshgeom),                             intent(in   ) :: meshgeom

   if (atStart) then
      localOffsets(1:gridPointsCount+1)=(/ 0.0d0, localOffsets(1:gridPointsCount) /)
      localGpsX(1:gridPointsCount+1)=(/ meshgeom%nnodex(meshgeom%nedge_nodes(1,ibran)), localGpsX(1:gridPointsCount) /)
      localGpsY(1:gridPointsCount+1)=(/ meshgeom%nnodey(meshgeom%nedge_nodes(1,ibran)), localGpsY(1:gridPointsCount) /)
      localGpsID(1:gridPointsCount+1)=(/ idMeshNodesInNetworkNodes(meshgeom%nedge_nodes(1,ibran)), localGpsID(1:gridPointsCount) /)
   else
      localOffsets(1:gridPointsCount+1)=(/ localOffsets(1:gridPointsCount), meshgeom%nbranchlengths(ibran) /)
      localGpsX(1:gridPointsCount+1)=(/ localGpsX(1:gridPointsCount), meshgeom%nnodex(meshgeom%nedge_nodes(2,ibran)) /)
      localGpsY(1:gridPointsCount+1)=(/ localGpsY(1:gridPointsCount), meshgeom%nnodey(meshgeom%nedge_nodes(2,ibran)) /)
      localGpsID(1:gridPointsCount+1)=(/ localGpsID(1:gridPointsCount), idMeshNodesInNetworkNodes(meshgeom%nedge_nodes(2,ibran)) /)
   end if
   gridPointsCount = gridPointsCount + 1
   end subroutine add_point

   subroutine read_1d_ugrid(network, ioncid, dflowfm)

   use io_netcdf
   use io_ugrid
   use m_hash_search
   use gridgeom
   use meshdata

   implicit none

   type(t_network), target, intent(inout) :: network
   integer, intent(in)                    :: ioncid
   logical, optional, intent(inout)       :: dflowfm

   integer                   :: ierr
   integer                   :: numMesh
   integer                   :: meshIndex
   integer                   :: networkIndex
   integer                   :: startIndex
   
   character(len=ug_idsLen), allocatable, dimension(:)              :: branchids
   character(len=ug_idsLongNamesLen), allocatable, dimension(:)     :: branchlongnames
   character(len=ug_idsLen), allocatable, dimension(:)              :: nodeids
   character(len=ug_idsLongNamesLen), allocatable, dimension(:)     :: nodelongnames
   character(len=ug_idsLen), allocatable, dimension(:)              :: gpsID
   character(len=ug_idsLongNamesLen), allocatable, dimension(:)     :: gpsIDLongnames
   character(len=ug_idsLongNamesLen)                                :: network1dname
   character(len=ug_idsLongNamesLen)                                :: mesh1dname
   integer, parameter                                               :: nodesOnBranchVertices = 0


   !< Structure where all mesh is stored. Internal arrays are all pointers
   type(t_ug_meshgeom) :: meshgeom

   ! Make indexes 1-based
   startIndex = 1


   ierr = ionc_get_mesh_count(ioncid, numMesh)
   if (ierr .ne. 0) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Reading Number of Meshes')
   endif
   if (numMesh < 1) then
      call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Data Missing')
   endif

   ! Get Index of 1D-Mesh
   ierr = ionc_get_1d_mesh_id_ugrid(ioncid, meshIndex)
   if (ierr .ne. 0 .or. meshIndex <= 0) then
      if (present(dflowfm)) then
         call SetMessage(LEVEL_INFO, 'Network UGRID-File: No 1D-Mesh Present, Skipped 1D')
         network%loaded = .false.
         return
      else
         call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Reading Mesh ID')
      endif
   endif

   ! Get Index of 1D-Network
   ierr = ionc_get_1d_network_id_ugrid(ioncid, networkIndex)
   if (ierr .ne. 0 .or. networkIndex <= 0) then
      if (present(dflowfm)) then
         call SetMessage(LEVEL_INFO, 'Network UGRID-File: No 1D-Network Present, Skipped 1D')
         network%loaded = .false.
         return
      else
         call SetMessage(LEVEL_FATAL, 'Network UGRID-File: Error Reading Mesh ID')
      endif
   endif

   ! Get all data in one call: mesh and network
   ierr =  ionc_get_meshgeom(ioncid, meshIndex, networkIndex, meshgeom, startIndex, .true., &
      branchids, branchlongnames, nodeids, nodelongnames, & !1d network character variables
      gpsID, gpsIDLongnames, network1dname, mesh1dname)     !1d grid character variables

   ! Fill the flow1d m_network::network data structure based on the meshgeom from file.
   ierr = construct_network_from_meshgeom(network, meshgeom, branchids, branchlongnames, nodeids, nodelongnames, & 
      gpsID, gpsIDLongnames, network1dname, mesh1dname, nodesOnBranchVertices)

   
   !deallocate memory
   ierr = t_ug_meshgeom_destructor(meshgeom)
   deallocate(branchids)
   deallocate(branchlongnames)
   deallocate(nodeids)
   deallocate(nodelongnames)
   deallocate(gpsID)
   deallocate(gpsIDLongnames)
   
   end subroutine read_1d_ugrid

   
   subroutine adminBranchOrders(brs)
      type (t_branchset), intent(inout) :: brs
      
      integer i, ibr, j
      type(t_branch), pointer :: pbr, pbr2
      
      do ibr = 1, brs%count
         pbr => brs%branch(ibr)
         if (pbr%orderNumber > 0) then
            do j = 1, 2
               if (pbr%nextBranch(j) < 0) then
                  ! find neighbouring branch
                  do i = ibr+1, brs%count
                     pbr2 => brs%branch(i)
                     if (pbr%ordernumber == pbr2%ordernumber) then
                        if (pbr%nodeIndex(j) == pbr2%nodeIndex(1)) then
                           ! found one
                           pbr%nextBranch(j) = i
                           pbr2%nextBranch(1)= ibr
                           ! finished
                           cycle
                        elseif (pbr%nodeIndex(j) == pbr2%nodeIndex(2)) then
                           ! found one
                           pbr%nextBranch(j) = i
                           pbr2%nextBranch(2)= ibr
                           ! finished
                           cycle
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      
   end subroutine adminBranchOrders

   subroutine readNode(nds, md_ptr)
   
      implicit none
   
      type(t_nodeSet), target, intent(inout) :: nds
      type(tree_data), pointer, intent(in)   :: md_ptr    

      character(len=IdLen)                   :: nodeId
      character(len=IdLen)                   :: nodeName
      double precision                       :: x
      double precision                       :: y
      logical                                :: success
      
      call  prop_get_string(md_ptr, 'node', 'id', nodeId, success)
      if (.not. success) then
         call SetMessage(LEVEL_FATAL, 'Error Reading Node ID')
      endif
      
      call prop_get_string(md_ptr, 'node', 'name', nodeName, success)

      call prop_get_double(md_ptr, 'node', 'x', x, success)
      if (success) call prop_get_double(md_ptr, 'node', 'y', y, success)
      
      if (.not. success) then
         call SetMessage(LEVEL_FATAL, 'Error Reading Node '''//trim(nodeId)//'''')
      endif

      nds%Count = nds%Count+1
      if (nds%Count > nds%Size) then
         call realloc(nds)
      endif
      
      nds%node(nds%Count)%id       = nodeId
      nds%node(nds%Count)%name     = nodeName
      nds%node(nds%Count)%index    = nds%count
      nds%node(nds%Count)%nodetype = nt_NotSet
      nds%node(nds%Count)%numberOfConnections = 0
      nds%node(nds%Count)%x = x
      nds%node(nds%Count)%y = y
      
   end subroutine readNode
   
   subroutine storeNodes(nds, nNodes, nodesX, nodesY, nodeids, nodelongnames)
   
      implicit none
   
      type(t_nodeSet), target, intent(inout)             :: nds
      integer, intent(in)                                :: nNodes
      double precision, dimension(nNodes), intent(in)    :: nodesX 
      double precision, dimension(nNodes), intent(in)    :: nodesY
      character(len=*), dimension(nNodes), intent(in)    :: nodeids
      character(len=*), dimension(nNodes) , intent(in)   :: nodelongnames
      
      integer                                :: iNode
      
      
      do iNode = 1, nNodes
      
         nds%Count = nds%Count + 1
         if (nds%Count > nds%Size) then
            call realloc(nds)
         endif
      
         nds%node(nds%Count)%id                  = nodeids(iNode)
         nds%node(nds%Count)%name                = nodelongnames(iNode)
         nds%node(nds%Count)%index               = nds%count
         nds%node(nds%Count)%nodetype            = nt_NotSet
         nds%node(nds%Count)%numberOfConnections = 0
         nds%node(nds%Count)%x                   = nodesX(iNode)
         nds%node(nds%Count)%y                   = nodesY(iNode)
         
      enddo
      
      call fill_hashtable(nds)
      
   end subroutine storeNodes
   
   subroutine readBranch(brs, nds, md_ptr)
   
      use m_branch
      
      implicit none

      type(t_branchSet), target, intent(inout) :: brs
      type(t_nodeSet), target, intent(inout)   :: nds
      type(tree_data), pointer, intent(in)     :: md_ptr
      
      ! Local Variables
      integer                                  :: ibr
      type(t_branch), pointer                  :: pbr
      type(t_node), pointer                    :: node
      logical                                  :: success
      integer                                  :: istat
      integer                                  :: ibegNode
      integer                                  :: iendNode
      integer                                  :: orderNumber
      integer                                  :: gridPointsCount
      integer                                  :: uPointsCount
      integer                                  :: igr
      integer                                  :: gridIndex
      integer                                  :: j
      integer                                  :: ip1
      integer                                  :: ip2
      character(len=IdLen)                     :: branchId
      character(len=IdLen)                     :: begNodeId
      character(len=IdLen)                     :: endNodeId
      character(len=IdLen)                     :: Chainage
      
      double precision, allocatable, dimension(:)     :: gpX
      double precision, allocatable, dimension(:)     :: gpY
      double precision, allocatable, dimension(:)     :: gpchainages
      character(len=IdLen), allocatable, dimension(:) :: gpID
      
      brs%Count = brs%Count + 1
      ibr = brs%Count
      if (brs%Count > brs%Size) then
         call realloc(brs)
      endif
      
      pbr =>brs%branch(brs%Count)
      
      call  prop_get_string(md_ptr, 'branch', 'id', branchId, success)
      if (.not. success) then
         call SetMessage(LEVEL_FATAL, 'Error Reading Branch ID')
      endif

      call  prop_get_string(md_ptr, 'branch', 'fromnode', begNodeId, success)
      if (success) call  prop_get_string(md_ptr, 'branch', 'tonode', endNodeId, success)
      if (success) call  prop_get_integer(md_ptr, 'branch', 'order', ordernumber, success)
      if (success) call  prop_get_integer(md_ptr, 'branch', 'gridpointscount', gridPointsCount, success)

      if (.not. success) then
         call SetMessage(LEVEL_FATAL, 'Error Reading Branch '''//trim(branchId)//'''')
      endif
      
      ibegNode = hashsearch(nds%hashlist, begNodeId)
      if (ibegNode <= 0) then
         write(msgbuf, '(4a)') trim(branchId), ': fromNode ''', trim(begNodeId), ''' does not exist'
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif
      
      iendNode = hashsearch(nds%hashlist, endNodeId)
      if (iendNode <= 0) then
         write(msgbuf ,'(4a)') trim(branchId), ': toNode ''', trim(endNodeId), ''' does not exist'
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif

      if (ibegNode == iendNode) then
         write(msgbuf, '(5a)') trim(branchId), ': fromNode ''', trim(begNodeId) , ''' is identical to toNode ''', trim(endNodeId)//''''
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif
      
      pbr%id                           = branchID
      pbr%index                        = ibr
      pbr%FromNode                     => nds%node(ibegNode)
      pbr%FromNode%numberOfConnections = pbr%FromNode%numberOfConnections + 1
      pbr%ToNode                       => nds%node(iendNode)
      pbr%ToNode%numberOfConnections   = pbr%toNode%numberOfConnections + 1
      pbr%orderNumber                  = orderNumber
      pbr%nextBranch                   = -1
      pbr%nodeIndex(1)                 = ibegNode
      pbr%nodeIndex(2)                 = iendNode
      ! The Gridpoints
      call realloc(gpX, gridPointsCount, stat=istat)
      if (istat == 0) call realloc(gpY, gridPointsCount, stat=istat)
      if (istat == 0) call realloc(gpchainages, gridPointsCount, stat=istat)
      if (istat .ne. 0) then
         call SetMessage(LEVEL_FATAL, 'Reading Branch: Error allocating Gridpoint Arrays')
      endif
      
      call prop_get_doubles(md_ptr, 'branch', 'gridPointX', gpX, gridPointsCount, success)
      if (success) call prop_get_doubles(md_ptr, 'branch', 'gridPointY', gpY, gridPointsCount, success)
      if (success) call prop_get_doubles(md_ptr, 'branch', 'gridPointOffsets', gpchainages, gridPointsCount, success)
      if (.not. success) then
         call SetMessage(LEVEL_FATAL, 'Error Reading Grid Data from Branch '''//trim(branchId)//'''')
      endif
      
      ! Check chainages of Grid Points
      do igr = 1, gridPointsCount - 1
         if (gpchainages(igr) + minSectionLength > gpchainages(igr + 1)) then
            ! Two adjacent gridpoints too close
            write (msgbuf, '(a, a, a, g11.4, a, g11.4)' ) 'Two grid points on branch ''', trim(branchid), ''' are too close at chainage ',               &
                                    gpchainages(igr), ' and ', gpchainages(igr + 1)
            call SetMessage(LEVEL_WARN, msgbuf)
         endif 
      enddo

      pbr%gridPointsCount = gridPointsCount
      uPointsCount        = pbr%gridPointsCount - 1
      pbr%uPointsCount    = uPointsCount
      
      call realloc(pbr%gridPointschainages, pbr%gridPointsCount)
      call realloc(pbr%uPointschainages, pbr%uPointsCount)
      call realloc(pbr%Xs, pbr%gridPointsCount)
      call realloc(pbr%Ys, pbr%gridPointsCount)

      ip1 = brs%gridPointsCount + 1
      brs%gridPointsCount = brs%gridPointsCount + gridPointsCount
      ip2 = brs%gridPointsCount
      pbr%StartPoint        = ip1
      pbr%gridPointschainages = gpchainages
      pbr%uPointschainages    = (pbr%gridPointschainages(1:uPointsCount) + pbr%gridPointschainages(2:uPointsCount+1) ) / 2.0d0
      pbr%length            = gpchainages(gridPointsCount)
      pbr%Xs                = gpX
      pbr%Ys                = gpY
      pbr%fromNode%x        = gpX(1)
      pbr%fromNode%y        = gpY(1)
      pbr%toNode%x          = gpX(gridPointsCount)
      pbr%toNode%y          = gpY(gridPointsCount)
      pbr%iTrench           = 0

      do j = 1, 2
         if (j==1) then
            node => pbr%fromNode
            gridIndex = ip1
         else
            node => pbr%toNode
            gridIndex = ip2
         endif
         if (node%nodeType == nt_NotSet) then
            ! probably end node (until proved otherwise
            node%nodeType = nt_endNode
         node%gridNumber = gridIndex
         elseif (node%nodeType == nt_endNode) then
            ! Already one branch connected, so not an endNode
            node%nodeType = nt_LinkNode
            node%gridNumber = 0
         endif
         if (node%numberOfConnections > nds%maxNumberOfConnections) then
            nds%maxNumberOfConnections = node%numberOfConnections
         endif
      enddo
      
      ! Get and Set grid point IDs
      call realloc(gpID, gridPointsCount, stat=istat)
      if (istat .ne. 0) then
         call SetMessage(LEVEL_FATAL, 'Reading Branch: Error allocating Gridpoint Array for IDs')
      endif
      gpID = ' '

      call prop_get_strings(md_ptr, 'branch', 'gridPointIds', gridPointsCount, gpID, success)
      call realloc(pbr%gridPointIDs, gridPointsCount)

      if (success) then
         pbr%gridPointIDs = gpID
      else
         ! Create grid point IDs
         pbr%gridPointIDs = ' '
         do igr = 1, gridPointsCount
            write(Chainage, '(F35.3)') pbr%gridPointschainages(igr)
            pbr%gridPointIDs(igr) = trim(pbr%id)//'_'//trim(adjustl(Chainage))
         enddo
         
      endif
      
      ! Clear Grid Point Arrays
      if (allocated(gpX)) deallocate(gpX, stat=istat)
      if (istat == 0 .and. allocated(gpY)) deallocate(gpY, stat=istat)
      if (istat == 0 .and. allocated(gpchainages)) deallocate(gpchainages, stat=istat)
      if (istat == 0 .and. allocated(gpID)) deallocate(gpID, stat=istat)
      if (istat .ne. 0) then
         call SetMessage(LEVEL_FATAL, 'Reading Branch: Error Deallocating Gridpoint Arrays')
      endif
      
   end subroutine readBranch

   subroutine storeBranch(brs, nds, branchId, begNodeId, endNodeId, ordernumber, gridPointsCount, gpX, gpY, gpchainages, gpID, active_branch, my_rank)
   
      use m_branch
      
      implicit none

      type(t_branchSet), target, intent(inout)       :: brs
      type(t_nodeSet), target, intent(inout)         :: nds
      character(len=*), intent(in)                   :: branchId
      character(len=*), intent(in)                   :: begNodeId
      character(len=*), intent(in)                   :: endNodeId
      integer, intent(in)                            :: orderNumber

      integer, intent(in)                                          :: gridPointsCount
      double precision, dimension(gridPointsCount), intent(in)     :: gpX
      double precision, dimension(gridPointsCount), intent(in)     :: gpY
      double precision, dimension(gridPointsCount), intent(in)     :: gpchainages
      character(len=*), dimension(gridPointsCount), intent(in)     :: gpID
      integer                                     , intent(in)     :: my_rank
      integer                                     , intent(in)     :: active_branch

      ! Local Variables
      integer                                  :: ibr
      type(t_branch), pointer                  :: pbr
      type(t_node), pointer                    :: node
      integer                                  :: ibegNode
      integer                                  :: iendNode
      integer                                  :: uPointsCount
      integer                                  :: igr
      integer                                  :: gridIndex
      integer                                  :: j
      integer                                  :: ip1
      integer                                  :: ip2
      character(len=IdLen)                     :: Chainage
      
      brs%Count = brs%Count + 1
      ibr = brs%Count
      if (brs%Count > brs%Size) then
         call realloc(brs)
      endif
      
      pbr =>brs%branch(brs%Count)
      
      ibegNode = hashsearch(nds%hashlist, begNodeId)
      if (ibegNode <= 0) then
         write(msgbuf, '(4a)') trim(branchId), ': fromNode ''', trim(begNodeId), ''' does not exist'
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif
      
      iendNode = hashsearch(nds%hashlist, endNodeId)
      if (iendNode <= 0) then
         write(msgbuf ,'(4a)') trim(branchId), ': toNode ''', trim(endNodeId), ''' does not exist'
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif

      if (ibegNode == iendNode) then
         write(msgbuf, '(5a)') trim(branchId), ': fromNode ''', trim(begNodeId) , ''' is identical to toNode ''', trim(endNodeId)//''''
         call SetMessage(LEVEL_FATAL, msgbuf)
      endif
      
      pbr%id                           = branchID
      pbr%index                        = ibr
      pbr%FromNode                     => nds%node(ibegNode)
      pbr%FromNode%numberOfConnections = pbr%FromNode%numberOfConnections + 1
      pbr%ToNode                       => nds%node(iendNode)
      pbr%ToNode%numberOfConnections   = pbr%toNode%numberOfConnections + 1
      pbr%orderNumber                  = orderNumber
      pbr%nextBranch                   = -1
      pbr%nodeIndex(1)                 = ibegNode
      pbr%nodeIndex(2)                 = iendNode
      pbr%active_branch                = active_branch

      ! The Gridpoints

      ! Check chainages of Grid Points
      do igr = 1, gridPointsCount - 1
         if (gpchainages(igr) + minSectionLength > gpchainages(igr + 1)) then
            ! Two adjacent gridpoints too close
            write (msgbuf, '(3a, g11.4, a, g11.4)' ) 'Two grid points on branch ''', trim(branchid), ''' are too close at chainage ',               &
                                    gpchainages(igr), ' and ', gpchainages(igr + 1)
            if (my_rank >= 0) then
               write(msgbuf,'(a,i4,2a)') 'Rank = ', my_rank, ': ', trim(msgbuf)
            end if
            call SetMessage(LEVEL_WARN, msgbuf)
         endif 
      enddo

      pbr%gridPointsCount = gridPointsCount
      uPointsCount        = max(0, pbr%gridPointsCount - 1)
      pbr%uPointsCount    = uPointsCount
      
      call realloc(pbr%gridPointschainages, pbr%gridPointsCount)
      call realloc(pbr%uPointschainages, pbr%uPointsCount)
      call realloc(pbr%Xs, pbr%gridPointsCount)
      call realloc(pbr%Ys, pbr%gridPointsCount)

      ip1 = brs%gridPointsCount + 1
      brs%gridPointsCount = brs%gridPointsCount + gridPointsCount
      ip2 = brs%gridPointsCount
      pbr%StartPoint         = ip1
      pbr%gridPointschainages = gpchainages
      pbr%uPointschainages    = (pbr%gridPointschainages(1:uPointsCount) + pbr%gridPointschainages(2:uPointsCount+1) ) / 2.0d0
      pbr%length            = gpchainages(gridPointsCount)
      pbr%Xs                = gpX
      pbr%Ys                = gpY
      pbr%fromNode%x        = gpX(1)
      pbr%fromNode%y        = gpY(1)
      pbr%toNode%x          = gpX(gridPointsCount)
      pbr%toNode%y          = gpY(gridPointsCount)
      pbr%iTrench           = 0

      do j = 1, 2
         if (j==1) then
            node => pbr%fromNode
            gridIndex = ip1
         else
            node => pbr%toNode
            gridIndex = ip2
         endif
         if (node%nodeType == nt_NotSet) then
            ! probably end node (until proved otherwise)
            node%nodeType = nt_endNode
            if (gridPointsCount > 0) then
               node%gridNumber = gridIndex
            else
               node%gridNumber = -1
            endif
         elseif (node%nodeType == nt_endNode) then
            ! Already one branch connected, so not an endNode
            node%nodeType = nt_LinkNode
            node%gridNumber = 0
         endif
         if (node%numberOfConnections > nds%maxNumberOfConnections) then
            nds%maxNumberOfConnections = node%numberOfConnections
         endif
      enddo
      
      ! Set grid point IDs
      call realloc(pbr%gridPointIDs, gridPointsCount)

      do igr = 1, gridPointsCount
      
         if (gpID(igr) .eq. ' ') then
            write(Chainage, '(F35.3)') pbr%gridPointschainages(igr)
            pbr%gridPointIDs(igr) = trim(pbr%id)//'_'//trim(adjustl(Chainage))
         else
            pbr%gridPointIDs(igr) = gpID(igr)
         endif
         
      enddo
      
   end subroutine storeBranch

   subroutine read_network_cache(ibin, network)
   
      type(t_network), intent(inout)  :: network
      integer, intent(in)             :: ibin

      call read_node_cache(ibin, network)
      
      call read_branch_cache(ibin, network)
   
   end subroutine read_network_cache
   
   subroutine write_network_cache(ibin, network)
   
      type(t_network), intent(in)     :: network
      integer, intent(in)             :: ibin

      call write_node_cache(ibin, network%nds)
      
      call write_branch_cache(ibin, network%brs)
   
   end subroutine write_network_cache
   
   subroutine read_node_cache(ibin, network)
   
      type(t_network), intent(inout)    :: network
      integer, intent(in)               :: ibin
      
      type(t_node), pointer             :: pnod
      integer                           :: i

      read(ibin) network%nds%count
      
      network%nds%growsby = network%nds%count + 2
      call realloc(network%nds)

      read(ibin) network%nds%maxNumberOfConnections
      
      network%nds%LevelBoundaryCount = 0
      network%nds%DisBoundaryCount   = 0
      network%nds%bndCount           = 0

      do i = 1, network%nds%Count
      
         pnod => network%nds%node(i)
         
         read(ibin) pnod%id
         read(ibin) pnod%name
         read(ibin) pnod%index
         read(ibin) pnod%nodeType

         read(ibin) pnod%x
         read(ibin) pnod%y
         read(ibin) pnod%gridNumber
         
         read(ibin) pnod%numberOfConnections
      
      enddo
 
      call read_hash_list_cache(ibin, network%nds%hashlist)

   end subroutine read_node_cache
   
   subroutine write_node_cache(ibin, nds)
   
      type(t_nodeSet), intent(in)     :: nds
      integer, intent(in)             :: ibin

      type(t_node), pointer           :: pnod
      integer                         :: i

      write(ibin) nds%Count
      
      write(ibin) nds%maxNumberOfConnections

      do i = 1, nds%Count
      
         pnod => nds%node(i)
         
         write(ibin) pnod%id
         write(ibin) pnod%name
         write(ibin) pnod%index
         write(ibin) pnod%nodeType

         write(ibin) pnod%x
         write(ibin) pnod%y
         write(ibin) pnod%gridNumber
         
         write(ibin) pnod%numberOfConnections
      
      enddo
 
      call write_hash_list_cache(ibin, nds%hashlist)
      
   end subroutine write_node_cache
   
   subroutine read_branch_cache(ibin, network)
   
      type(t_network), intent(inout)   :: network
      integer, intent(in)              :: ibin
   
      type(t_branch), pointer          :: pbrn
      integer                          :: i
      integer                          :: j

      read(ibin) network%brs%Count
      
      network%brs%growsby = network%brs%count + 2
      call realloc(network%brs)

      read(ibin) network%brs%gridpointsCount

      do i = 1, network%brs%Count
      
         pbrn => network%brs%branch(i)
         
         read(ibin) pbrn%id
         read(ibin) pbrn%index
         read(ibin) pbrn%name
         read(ibin) pbrn%length
         read(ibin) pbrn%orderNumber
         
         read(ibin) pbrn%iTrench
         read(ibin) pbrn%flapGate
         
         read(ibin) (pbrn%nextBranch(j), j = 1, 2)
         
         read(ibin) (pbrn%nodeIndex(j), j = 1, 2)
         pbrn%FromNode => network%nds%node(pbrn%nodeIndex(1))
         pbrn%ToNode => network%nds%node(pbrn%nodeIndex(2))

         read(ibin) pbrn%gridPointsCount
         
         allocate(pbrn%gridPointschainages(pbrn%gridPointsCount))
         allocate(pbrn%gridPointIDs(pbrn%gridPointsCount))
         allocate(pbrn%Xs(pbrn%gridPointsCount))
         allocate(pbrn%Ys(pbrn%gridPointsCount))
         
         read(ibin) (pbrn%gridPointschainages(j), j = 1, pbrn%gridPointsCount)   
         read(ibin) (pbrn%gridPointIDs(j), j = 1, pbrn%gridPointsCount)   
         read(ibin) (pbrn%Xs(j), j = 1, pbrn%gridPointsCount)   
         read(ibin) (pbrn%Ys(j), j = 1, pbrn%gridPointsCount)   
      
         read(ibin) pbrn%uPointsCount
         
         allocate(pbrn%uPointschainages(pbrn%uPointsCount))

         read(ibin) (pbrn%uPointschainages(j), j = 1, pbrn%uPointsCount)   

         read(ibin) pbrn%StartPoint

      enddo
 
      call read_hash_list_cache(ibin, network%brs%hashlist)

   end subroutine read_branch_cache
   
   subroutine write_branch_cache(ibin, brs)
   
      type(t_branchSet), intent(in)   :: brs
      integer, intent(in)             :: ibin

      type(t_branch), pointer         :: pbrn
      integer                         :: i
      integer                         :: j

      write(ibin) brs%Count
      
      write(ibin) brs%gridpointsCount

      do i = 1, brs%Count
      
         pbrn => brs%branch(i)
         
         write(ibin) pbrn%id
         write(ibin) pbrn%index
         write(ibin) pbrn%name
         write(ibin) pbrn%length
         write(ibin) pbrn%orderNumber
         
         write(ibin) pbrn%iTrench
         write(ibin) pbrn%flapGate

         write(ibin) (pbrn%nextBranch(j), j = 1, 2)
         write(ibin) (pbrn%nodeIndex(j), j = 1, 2)

         write(ibin) pbrn%gridPointsCount
         write(ibin) (pbrn%gridPointschainages(j), j = 1, pbrn%gridPointsCount)   
         write(ibin) (pbrn%gridPointIDs(j), j = 1, pbrn%gridPointsCount)   
         write(ibin) (pbrn%Xs(j), j = 1, pbrn%gridPointsCount)   
         write(ibin) (pbrn%Ys(j), j = 1, pbrn%gridPointsCount)   
      
         write(ibin) pbrn%uPointsCount
         write(ibin) (pbrn%uPointschainages(j), j = 1, pbrn%uPointsCount)   

         write(ibin) pbrn%StartPoint

      enddo
 
      call write_hash_list_cache(ibin, brs%hashlist)

   end subroutine write_branch_cache
   
   
end module m_1d_networkreader
