module m_readRetentions
!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2018.                                
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
   use m_network
   use m_Storage
   use m_tables

   use properties
   use m_hash_search

   implicit none

   private

   public readRetentions

   contains

   subroutine readRetentions(network, retentionFile)

      implicit none
      
      type(t_network), intent(inout) :: network
      character*(*), intent(in)      :: retentionFile

      logical                                       :: success
      type(tree_data), pointer                      :: md_ptr 
      integer                                       :: istat
      integer                                       :: numstr
      integer                                       :: i

      character(len=IdLen)                          :: retentionID
      character(len=IdLen)                          :: branchID
      character(len=IdLen)                          :: storageType
      logical                                       :: useTable
      
      double precision                              :: Chainage
      integer                                       :: branchIdx
      integer                                       :: gridPoint
      type(t_storage), pointer                      :: pSto
      
      integer                                       :: nLevels
      integer                                       :: interPolate
      double precision, allocatable, dimension(:)   :: storLevels
      double precision, allocatable, dimension(:)   :: storAreas

      call tree_create(trim(retentionFile), md_ptr, maxlenpar)
      call prop_file('ini',trim(retentionFile),md_ptr, istat)
      
      numstr = 0
      if (associated(md_ptr%child_nodes)) then
         numstr = size(md_ptr%child_nodes)
      end if

      do i = 1, numstr
         
         if (tree_get_name(md_ptr%child_nodes(i)%node_ptr) .eq. 'retention') then
            
            ! Read Data
            
            call prop_get_string(md_ptr%child_nodes(i)%node_ptr, 'retention', 'id', retentionID, success)
            if (success) call prop_get_string(md_ptr%child_nodes(i)%node_ptr, 'retention', 'branchid', branchID, success)
            if (success) call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'retention', 'chainage', Chainage, success)
            if (.not. success) then
               call SetMessage(LEVEL_FATAL, 'Error Reading Retention '''//trim(retentionID)//'''')
            endif

            call prop_get_string(md_ptr%child_nodes(i)%node_ptr, 'retention', 'storagetype', storageType, success)
            if (.not. success) storageType = 'Reservoir'
            
            branchIdx = hashsearch(network%brs%hashlist, branchID)
            if (branchIdx <= 0) Then
               call SetMessage(LEVEL_ERROR, 'Error Reading Retention '''//trim(retentionID)//''': Branch: '''//trim(branchID)//''' not Found')
               exit
            endif

            network%storS%Count = network%storS%Count + 1
            if (network%storS%Count > network%storS%Size) then
               call realloc(network%storS)
            endif
      
            pSto => network%storS%stor(network%storS%Count)
            nullify(pSto%storageArea)
            nullify(pSto%streetArea)

            ! Bcause of the complicated data structure of SOBEK storage in 'connection nodes'
            ! must be separated from the ordinary gridpoints
            gridPoint = getCalcPoint(network%brs, branchIdx, Chainage)
            if (gridPoint == network%brs%branch(branchIdx)%points(1) ) then
               gridPoint = -network%brs%branch(branchIdx)%fromNode%index
               psto%node_index = network%brs%branch(branchIdx)%fromNode%index
               psto%local_grid_index = -1
               branchIdx             = -1
            elseif (gridPoint == network%brs%branch(branchIdx)%points(2)) then
               gridPoint = -network%brs%branch(branchIdx)%toNode%index
               psto%node_index = network%brs%branch(branchIdx)%toNode%index
               branchIdx             = -1
               psto%local_grid_index = -1
            else
               psto%local_grid_index = gridPoint - network%brs%branch(branchIdx)%points(1) + 1 
            endif
            
            pSto%id        = retentionID
            pSto%gridPoint = gridPoint
            network%storS%mapping(gridPoint) = network%storS%Count
            psto%branch_index     = branchIdx

            if (storageType == 'Closed') then
               pSto%storageType = nt_Closed
            elseif (storageType == 'Loss') then
               pSto%storageType = nt_Loss
            else
               pSto%storageType = nt_Reservoir
            endif

            call prop_get_integer(md_ptr%child_nodes(i)%node_ptr, 'retention', 'numlevels', nLevels, success)
            if (success) then
               useTable = .true.
            else
               nLevels  = 2
               useTable = .false.
            endif
            
            ! Allocate Arrays
            call realloc(storLevels, nLevels, stat=istat)
            if (istat == 0) call realloc(storAreas, nLevels, stat=istat)
            if (istat .ne. 0) then
               call SetMessage(LEVEL_FATAL, 'Reading Retentions: Error Allocating Arrays')
            endif
            
            if (useTable) then

               call prop_get_doubles(md_ptr%child_nodes(i)%node_ptr, 'retention', 'levels', storLevels, nLevels, success)
               if (success)call prop_get_doubles(md_ptr%child_nodes(i)%node_ptr, 'retention', 'storageArea', storAreas, nLevels, success)
               if (.not. success) then
                  call SetMessage(LEVEL_FATAL, 'Reading Retentions: Error in Level/Storage Data')
               endif
      
               call prop_get_integer(md_ptr%child_nodes(i)%node_ptr, 'retention', 'interpolate', interPolate, success)
               if (.not. success) interPolate = 0
               
            else
               
               call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'retention', 'bedlevel', storLevels(1), success)
               if (success) call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'retention', 'area', storAreas(1), success)
               if (.not. success) then
                  call SetMessage(LEVEL_FATAL, 'Reading Retentions: Error in Level/Storage Data')
               endif
               
               call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'retention', 'streetlevel', storLevels(2), success)
               if (.not. success) storLevels(2) = storLevels(1) + 1.0d0
               
               call prop_get_double(md_ptr%child_nodes(i)%node_ptr, 'retention', 'streetlevelarea', storAreas(2), success)
               if (.not. success) storLevels(2) = storLevels(1)
               
               interPolate = 1
               
            endif

            call setTable(pSto%storageArea, interPolate, storLevels, storAreas, nLevels)

            ! Clear Arrays
            istat = 0
            if (allocated(storLevels)) deallocate(storLevels, stat=istat)
            if (istat == 0 .and. allocated(storAreas)) deallocate(storAreas, stat=istat)
            if (istat .ne. 0) then
               call SetMessage(LEVEL_FATAL, 'Reading Retentions: Error Deallocating Arrays')
            endif
      
         endif
      
      end do
      
      call fill_hashtable(network%storS)
      
      call tree_destroy(md_ptr)

   end subroutine readRetentions

end module m_readRetentions