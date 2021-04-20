!> Returns the index of a structure in the controllable value arrays.
!! Structure is identified by strtypename, e.g. 'pumps', and structure name, e.g., 'Pump01'.
!! Returned index can be used to directly address variables like, m_flowexternalforcings::qpump, zgate, etc.
subroutine getStructureIndex(strtypename, strname, index, is_in_network)
! NOTE: this will only return the GUI-used structures (i.e., the new gates and weirs via general structure, not the old ext-based damlevel and gateloweredgelevel).
! TODO: longer-term all structure sets run via channel_flow and t_structureset, cleanup this function then.
   use m_flowexternalforcings
   use m_hash_search, only: hashsearch
   use unstruc_channel_flow, only: network
   use m_longculverts

   implicit none
   character(len=*), intent(in   ) :: strtypename   !< the type of the structure: 'pumps', 'weirs', 'gates', ...
   character(len=*), intent(in   ) :: strname       !< Id/name of the requested structure, e.g. 'Pump01'
   integer,          intent(  out) :: index         !< Returned index of the found structure in its controllable value arrays. -1 when not found.
   logical,          intent(  out) :: is_in_network !< Whether or not the found structure is inside the network%sts set, or in FM global structure set. No meaning when structure not found.

   integer :: i, nstr, icgen
   integer, pointer :: cgen_mapping(:)
   index = -1
   is_in_network = .false.

   if (network%sts%count > 0) then
      ! TODO: when we allow non-unique ids between different structure types, select proper hashlist.
      index = hashsearch(network%sts%hashlist_structure, trim(strname))
      if (index > 0) then
         is_in_network = .true.
      end if
      return
   else
      ! Retry on the 2D structures in code below
      continue
   end if


   if (trim(strtypename) == 'pumps') then
      do i=1,npumpsg
         if (trim(pump_ids(i)) == trim(strname)) then
            if (L2pumpsg(i) - L1pumpsg(i) >= 0) then
               ! Only return this pump index if pump is active in flowgeom (i.e., at least 1 flow link associated)
               index = i
               exit
            end if
         end if
      end do
   else if (trim(strtypename) == 'sourcesinks') then
      do i=1,numsrc
         if (trim(srcname(i)) == trim(strname)) then
            index = i
            exit
         end if
      end do
   else if (trim(strtypename) == 'dambreak') then
      do i=1,ndambreaksg
         if (trim(dambreak_ids(i)) == trim(strname)) then
            if (L2dambreaksg(i) - L1dambreaksg(i) >= 0) then
               ! Only return this dambreak index if dambreak is active in flowgeom (i.e., at least 1 flow link associated)
               index = i
               exit
            end if
         end if
      end do
   else if (trim(strtypename) == 'longculverts') then
      do i=1,nlongculvertsg
         if (trim(longculverts(i)%id) == trim(strname)) then
            index = i
            exit
         end if
      end do
   else  ! generalstructure-based types
      select case(strtypename)
      case('weirs')
         cgen_mapping => weir2cgen
         nstr = nweirgen
      case('gates')
         cgen_mapping => gate2cgen
         nstr = ngategen
      case('generalstructures')
         cgen_mapping => genstru2cgen
         nstr = ngenstru
      case default
         nstr = 0
      end select

      do i=1,nstr
         icgen = cgen_mapping(i)
         if (trim(cgen_ids(icgen)) == trim(strname)) then
            if (L2cgensg(icgen) - L1cgensg(icgen) >= 0) then
               ! Only return this structure index if structure is active in flowgeom (i.e., at least 1 flow link associated)
               index = icgen
               exit
            end if
         end if
      end do
   end if

end subroutine getStructureIndex
