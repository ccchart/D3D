module m_network
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
   use networkTypes
   use m_crossSections
   use m_1d_structures
   use m_roughness
   use m_trenches
   use m_ExtraResistance
   use m_df1d_transport
   use m_ObservCrossSections
   use m_compound
   
   implicit none

   public realloc
   public dealloc
   public admin_network
   public initialize_1dadmin
   public getFrictionValue
 
   interface realloc
      module procedure realloc_1dadmin
   end interface realloc

   interface dealloc
      module procedure deallocNetwork
      module procedure dealloc_1dadmin
   end interface dealloc
   

   ! !TODO JN: zorg voor allocatie en initialisatie. en vullen van lin2ibr en lin2point uit adm%lin. -1 is missing value e.g. for a 2d link, length LINALL
   
   type, public :: t_administration_1d
      integer, allocatable          :: lin2str(:)                          !< indirection list, containing structure numbers for flowlinks.
                                                                           !< These structure numbers refer to the elements of network%sts%struct.
      logical                       :: hasStructures = .false.             !< Flag, true if one or more structures in model
      integer, allocatable          :: lin2ibr(:)                          !< indirection list, containing branch number on which the flow link is positioned 
      integer, allocatable          :: lin2point(:)                        !< indirection list, containing relative index of link on branch adm%lin2ibr(l)
                                                                           !> indirection list, containing local link index for 1d arrays.
      integer, allocatable          :: lin2local(:)
      integer, allocatable          :: lin2grid(:)
      type(t_chainage2cross), pointer :: line2cross(:) => null()             !< list containing cross section indices per u-chainage
      type(t_chainage2cross), pointer :: gpnt2cross(:) => null()             !< list containing cross section indices per gridpoint-chainage
      logical, allocatable          :: hysteresis_for_summerdike(:,:)      !< array indicating for hysteresis in summerdikes

      double precision, allocatable :: au_1d(:)
      double precision, allocatable :: conv_1d(:)
      double precision, allocatable :: minwidth1d(:)
   
      double precision, allocatable :: minconv(:)                      ! minimal 1d conveyance in link: C*P*sqrt(A/P)
   end type

   type, public   :: t_network
      integer                                   :: gridpointsCount         !< total number of gridpoints in network NDS%count - NDS%bndCount
      integer                                   :: l1dall                  !< total number of links (internal, boundary and compound links)
      integer                                   :: l1d                     !< total number of links (internal and boundary)
      integer                                   :: numk                    !< total number of links (internal and boundary)
      integer                                   :: numl                    !< total number of links (internal and boundary)
      logical                                   :: sferic                  !< flag indicating whether the used coordinate system is sferical or metric
      type(t_administration_1d)                 :: adm                     !< network administration
      type(t_nodeSet)                           :: nds                     !< set of nodes
      type(t_branchSet)                         :: brs                     !< set of branches
      type(t_CrossSectionSet)                   :: crs                     !< set of Cross-Sections
      type(t_StructureSet)                      :: sts                     !< structure list
      type(t_CompoundSet)                       :: cmps                    !< Administration compound structures
      type(t_RoughnessSet)                      :: rgs                     !< set containing roughness sections
      type(t_ObservationPointSet)               :: obs                     !< set of observation points
      type(t_lateralSet)                        :: lts                     !< set of lateral discharges
      type(t_storageSet)                        :: storS                   !< set containing storage in gridpoints
      type(t_trenchSet)                         :: trs                     !< set containing trenches (infiltration sewer tubes)
      type(t_CSDefinitionSet)                   :: CSDefinitions
      type(t_spatial_dataSet)                   :: spData
      type(t_boundarySet)                       :: boundaries
      type(t_transportSet)                      :: trans
      type(t_ObservCrossSectionSet)             :: observcrs               !< set of observation Cross-Sections 
      type(t_forcingList)                       :: forcinglist             !< Work list of read-in (structure) forcing data, to be initialized by calling kernel later.
      logical                                   :: loaded      = .false.
      logical                                   :: initialized = .false.
   end type
   
contains

   subroutine realloc_1dadmin(adm, all_links_count, oned_links_count)

      type(t_administration_1d)  :: adm
      integer, intent(in)        :: all_links_count
      integer, intent(in)        :: oned_links_count
      
      if (.not. allocated(adm%lin2str)) allocate(adm%lin2str(all_links_count))  
      if (.not. allocated(adm%lin2ibr)) allocate(adm%lin2ibr(all_links_count))   
      if (.not. allocated(adm%lin2point)) allocate(adm%lin2point(all_links_count)) 
      if (.not. allocated(adm%lin2local)) allocate(adm%lin2local(all_links_count)) 
      if (.not. allocated(adm%lin2grid)) allocate(adm%lin2grid(all_links_count)) 
      if (.not. associated(adm%line2cross)) allocate(adm%line2cross(all_links_count))
      if (.not. associated(adm%gpnt2cross)) allocate(adm%gpnt2cross(all_links_count))
      if (.not. allocated(adm%hysteresis_for_summerdike)) allocate(adm%hysteresis_for_summerdike(2,all_links_count))
      adm%hysteresis_for_summerdike = .true.
      if (.not. allocated(adm%au_1d)) allocate(adm%au_1d(oned_links_count))
      if (.not. allocated(adm%conv_1d)) allocate(adm%conv_1d(oned_links_count))
      if (.not. allocated(adm%minwidth1d)) allocate(adm%minwidth1d(oned_links_count))
      if (.not. allocated(adm%minconv)) allocate(adm%minconv(oned_links_count))   
      
   end subroutine realloc_1dadmin

   subroutine dealloc_1dadmin(adm)
      type(t_administration_1d)     :: adm

      if (allocated(adm%lin2str))      deallocate(adm%lin2str)
      if (allocated(adm%lin2ibr))      deallocate(adm%lin2ibr)
      if (allocated(adm%lin2point))    deallocate(adm%lin2point)
      if (allocated(adm%lin2local))    deallocate(adm%lin2local)
      if (associated(adm%line2cross))  deallocate(adm%line2cross)
      if (allocated(adm%lin2grid))    deallocate(adm%lin2grid)
      if (associated(adm%gpnt2cross))  deallocate(adm%gpnt2cross)
      if (allocated(adm%hysteresis_for_summerdike)) deallocate(adm%hysteresis_for_summerdike)
      if (allocated(adm%au_1d))        deallocate(adm%au_1d)
      if (allocated(adm%conv_1d))      deallocate(adm%conv_1d)
      if (allocated(adm%minwidth1d))   deallocate(adm%minwidth1d)
      if (allocated(adm%minconv))      deallocate(adm%minconv)

   end subroutine dealloc_1dadmin


   subroutine deallocNetwork(network)
      ! Modules
      use messagehandling
   
      implicit none
      ! Input/output parameters
      type(t_network), intent(inout) :: network
   
      ! Local variables

      ! Program code
      call dealloc(network%adm)
      call dealloc(network%nds)
      call dealloc(network%brs)
      call dealloc(network%crs)
      call dealloc(network%sts)
      call dealloc(network%cmps)
      call dealloc(network%rgs)
      call dealloc(network%obs)
      call dealloc(network%lts)
      call dealloc(network%storS)
      call dealloc(network%trs)
      call dealloc(network%CSDefinitions)
      call dealloc(network%spData)
      call dealloc(network%boundaries)
      call dealloc(network%trans)
      call dealloc(network%observcrs)
      network%loaded = .false.
   
   end subroutine deallocNetwork

   subroutine admin_network(network, ngrid, nlink)
      use m_node
      use m_branch
      use m_boundaryConditions
   
      type(t_network), intent(inout) :: network
      integer, intent(inout) :: ngrid
      integer, intent(inout) :: nlink
   
      integer ibr
      integer nnode
      integer ityp
      integer i
      integer icon, ibnd, typ
      integer nod
      integer, allocatable, dimension(:) :: itype
      integer, allocatable, dimension(:) :: iboun
      type(t_boundary), pointer :: bd, bd_co
   
      type(t_branch), pointer :: pbr
   
      call admin_nodes(network%nds, ngrid)
      call admin_branch(network%brs, ngrid, nlink)
      call admin_transport(network%trans, network%brs)
      network%gridpointsCount = ngrid
      network%l1d    = nlink
      network%l1dall = nlink

      nnode = network%nds%count
      allocate(itype(nnode), iboun(nnode))

      itype = 0
      iboun = 0
      do ityp = 1, 2
         do i = 1, network%boundaries%tp(ityp)%count
            bd =>network%boundaries%tp(ityp)%bd(i)
            itype(bd%node) = ityp
            iboun(bd%node) = i
         enddo
      enddo
   
      do ibr = 1, network%brs%count
         pbr => network%brs%branch(ibr)
         if (pbr%FromNode%nodeType == nt_LevelBoun .or. &
             pbr%FromNode%nodeType == nt_DischBoun) then
            nod = pbr%FromNode%index
            bd => network%boundaries%tp(itype(nod))%bd(iboun(nod))
            bd%igrid = pbr%FromNode%gridNumber
            bd%branch = ibr
            bd%intern = pbr%grd(2)
            bd%direction = 1
            bd%linknumber = pbr%lin(1)

            do icon = 1, transportPars%constituents_count
               typ = transportPars%co_h(icon)%boundary_index
               do ibnd = 1, network%boundaries%tp(typ)%Count
                  bd_co => network%boundaries%tp(typ)%bd(ibnd)
                  if (bd_co%node == nod) then
                     bd_co%igrid      = bd%igrid     
                     bd_co%intern     = bd%intern    
                     bd_co%branch     = bd%branch    
                     bd_co%direction  = bd%direction 
                     bd_co%linknumber = bd%linknumber
                     exit
                  endif
               enddo
            enddo
         endif
             
               
          if (pbr%ToNode%nodeType == nt_LevelBoun .or. &
             pbr%ToNode%nodeType == nt_DischBoun) then
            nod = pbr%ToNode%index
            bd => network%boundaries%tp(itype(nod))%bd(iboun(nod))
            bd%igrid      = pbr%ToNode%gridNumber
            bd%intern     = pbr%grd(pbr%gridPointsCount-1)
            bd%branch     = ibr
            bd%direction  = -1
            bd%linknumber = pbr%lin(pbr%uPointsCount)

            do icon = 1, transportPars%constituents_count
               typ = transportPars%co_h(icon)%boundary_index
               do ibnd = 1, network%boundaries%tp(typ)%Count
                  bd_co => network%boundaries%tp(typ)%bd(ibnd)
                  if (bd_co%node == nod) then
                     bd_co%igrid      = bd%igrid     
                     bd_co%intern     = bd%intern    
                     bd_co%branch     = bd%branch    
                     bd_co%direction  = bd%direction 
                     bd_co%linknumber = bd%linknumber
                     exit
                  endif
               enddo
            enddo
         endif   
      
         continue
          
      enddo
   
      !call initializeLaterals(network%brs, network%lts)

      deallocate(itype, iboun)
   
      !do ityp = 1, 3
      !   do i = 1, network%boundaries%tp(ityp)%count
      !      bd =>network%boundaries%tp(ityp)%bd(i)
      !      inod = network%boundaries%tp(ityp)%bd(i)%node
      !      if (ityp==H_BOUN) node(1,inod) = nt_LevelBoun
      !      if (ityp==Q_BOUN) node(1,inod) = nt_DischBoun
      !      node(3,inod) = i
      !      if (node(3,inod) == 0) then
      !         node(4,inod) = 0
      !      else
      !         node(4,inod) = inod
      !      endif
      !   enddo
      !enddo
      !do i = 1, network%boundaries%tp(3)%count
      !   inod = network%boundaries%tp(ityp)%bd(i)%node
      !   node(4,inod) = i
      !enddo


   end subroutine admin_network


   subroutine initialize_1dadmin(network, linall)
   
      use m_CrossSections
      use m_GlobalParameters
      
      type(t_network), intent(inout) :: network
      integer, intent(in)            :: linall
      
      integer :: ilnk
      integer :: igpt
      integer :: ll
      integer :: istru
      integer :: ibran
      integer :: m
      integer :: icrs1
      integer :: icrs2
      
      double precision                   :: f
      double precision                   :: dpu1
      double precision                   :: dpu2
      double precision                   :: chainage1
      double precision                   :: chainage2
      double precision                   :: chainageu
      double precision                   :: chainageg
      double precision                   :: chezy
      double precision                   :: as
      double precision                   :: wetperimeter
      type(t_administration_1d)          :: adm
      type(t_branch), pointer            :: pbran
      type(t_structure), pointer         :: pstru

      integer, allocatable, dimension(:) :: crossOrder
      integer, allocatable, dimension(:) :: lastAtBran
      integer                            :: icrsBeg
      integer                            :: icrsEnd
      double precision                   :: xBeg
      double precision                   :: xEnd
      integer                            :: i
      logical                            :: interpolDone
      logical                            :: initError = .false.

      character(20)                      :: chainageString
      
      call realloc(network%adm, linall, network%gridpointsCount)
      
      adm = network%adm
      
      adm%lin2str = -huge(1)
      if (network%sts%Count > 0) then
         adm%hasStructures = .true.
      else
         adm%hasStructures = .false.
      endif

      adm%lin2ibr   = -huge(1)
      adm%lin2point = -huge(1)  
      adm%lin2local = -huge(1)
      adm%lin2grid  = -huge(1)
      
      do ibran = 1, network%brs%Count
         pbran => network%brs%branch(ibran)
         do m = 1, pbran%uPointsCount
            if (pbran%lin(m) > 0) then
               adm%lin2ibr(pbran%lin(m)) = ibran
               adm%lin2point(pbran%lin(m)) = m
               adm%lin2local(pbran%lin(m)) = m
               adm%lin2grid(pbran%lin(m))  = pbran%grd(m)
            endif
         enddo
      enddo
      
      if (network%crs%Count > 0) then

         call useBranchOrders(network%crs, network%brs)
        
         allocate(crossOrder(network%crs%Count))
         allocate(lastAtBran(network%brs%Count))
                  
         call crossSectionsSort(network%crs, network%brs, crossOrder, lastAtBran)
         
         ! Fill in the Index for checking later
         do icrs1 = 1, network%crs%Count
            network%crs%cross(icrs1)%crossIndx = icrs1
         enddo
         
         ! Cross-Section indices for links
         adm%minconv    = 0.0d0
         adm%minwidth1d = 0.0d0
         
         do ibran = 1, network%brs%Count

            pbran   => network%brs%branch(ibran)

            if (ibran .eq. 1) then
               icrsBeg = 1
            else
               icrsBeg =lastAtBran(ibran - 1) + 1
            endif
            icrsEnd = lastAtBran(ibran)
            
            
            if (icrsbeg > icrsend) then
               
               call setmessage(LEVEL_WARN, 'No cross sections found on branch '//trim(pbran%id)//'. Using default rectangular cross section')
               do i = 1, pbran%uPointsCount
                  ilnk = pbran%lin(i)
                  adm%line2cross(ilnk)%c1 = -1
                  adm%line2cross(ilnk)%c2 = -1
                  adm%line2cross(ilnk)%f  = 1.0d0
                  adm%line2cross(ilnk)%distance  = 0d0
               enddo
               
               cycle
            endif
            
            icrs1 = icrsBeg
            icrs2 = icrsBeg

            xBeg = network%crs%cross(crossOrder(icrsBeg))%chainage
            xEnd = network%crs%cross(crossOrder(icrsEnd))%chainage

            do m = 1, pbran%uPointsCount

               chainageu = pbran%uPointschainages(m)
               ilnk = pbran%lin(m)
               
               if (icrsBeg == icrsEnd) then
                  
                  ! Just one Cross-Section
                  adm%line2cross(ilnk)%c1 = crossOrder(icrsBeg)
                  adm%line2cross(ilnk)%c2 = crossOrder(icrsBeg)
                  adm%line2cross(ilnk)%f  = 1.0d0
                  adm%line2cross(ilnk)%distance  = 0d0
                  interpolDone            = .true.
                  
               elseif (chainageu <= xBeg) then
                  
                  ! Before First Cross-Section
                  adm%line2cross(ilnk)%c1 = crossOrder(icrsBeg)
                  adm%line2cross(ilnk)%c2 = crossOrder(icrsBeg)
                  adm%line2cross(ilnk)%f  = 1.0d0
                  adm%line2cross(ilnk)%distance  = 0d0
                  interpolDone            = .true.
                  
               elseif (chainageu >= xEnd) then
                  
                  ! After Last Cross-Section
                  adm%line2cross(ilnk)%c1 = crossOrder(icrsEnd)
                  adm%line2cross(ilnk)%c2 = crossOrder(icrsEnd)
                  adm%line2cross(ilnk)%f  = 1.0d0
                  adm%line2cross(ilnk)%distance  = 0d0
                  interpolDone            = .true.
                  
               else
                  
                  chainage1 = network%crs%cross(crossOrder(icrs1))%chainage
                  chainage2 = network%crs%cross(crossOrder(icrs2))%chainage
                  adm%line2cross(ilnk)%distance  = chainage2 - chainage1
                  
                  if (.not. ((chainage1 <= chainageu) .and. (chainage2 >= chainageu))) then
                     
                     do i = icrs1, icrsEnd
                        if (network%crs%cross(crossOrder(i))%chainage >= chainageu) then
                           chainage2 = network%crs%cross(crossOrder(i))%chainage
                           icrs2 = i
                           exit
                        endif
                     enddo
                     
                     do i = icrsEnd, icrsBeg, -1
                        if (network%crs%cross(crossOrder(i))%chainage <= chainageu) then
                           chainage1 = network%crs%cross(crossOrder(i))%chainage
                           icrs1 = i
                           exit
                        endif
                     enddo
                     
                  endif
                  
                  interpolDone = .false.

               endif
                  
               if (ilnk > 0) then
                  
                  if (ibran == network%crs%cross(crossOrder(icrs2))%branchid) then
                     
                     if (.not. interpolDone) then
                        if (icrs1 == icrs2) then 
                            adm%line2cross(ilnk)%c1 = crossOrder(icrs1)
                            adm%line2cross(ilnk)%c2 = crossOrder(icrs2)
                            f = 1.0d0
                            adm%line2cross(ilnk)%f = f
                        else    
                            adm%line2cross(ilnk)%c1 = crossOrder(icrs1)
                            adm%line2cross(ilnk)%c2 = crossOrder(icrs2)
                            f  = (chainageu - chainage1) / (chainage2 - chainage1)
                            f = max(f, 0.0d0) 
                            f = min(f, 1.0d0) 
                            adm%line2cross(ilnk)%f = f
                        endif 
                     endif
                     
                  else
                     adm%line2cross(ilnk)%c1 = crossOrder(icrs1)
                     adm%line2cross(ilnk)%c2 = crossOrder(icrs1)
                     adm%line2cross(ilnk)%f  = 1.0d0
                  endif
                  
                  chezy = 0.0d0
                  ll = adm%lin2local(ilnk)
                  f = adm%line2cross(ilnk)%f
                  dpu1 = -network%crs%cross(adm%line2cross(ilnk)%c1)%bedLevel
                  dpu2 = -network%crs%cross(adm%line2cross(ilnk)%c2)%bedLevel
                  call GetCSParsFlow(adm%line2cross(ilnk), network%crs%cross, thresholdDry, as, wetPerimeter, adm%minwidth1d(ll))
               endif
                  
            enddo

         enddo
         
         ! Cross-Section indices for gridpoints
         do ibran = 1, network%brs%Count

            pbran   => network%brs%branch(ibran)

            if (ibran .eq. 1) then
               icrsBeg = 1
            else
               icrsBeg =lastAtBran(ibran - 1) + 1
            endif
            icrsEnd = lastAtBran(ibran)
            
            if (icrsBeg > icrsEnd) then
               ! branch without cross sections
               do m = 1, 2
                  igpt = pbran%grd(m)
                  adm%gpnt2cross(igpt)%c1 = 0
                  adm%gpnt2cross(igpt)%c2 = 0
                  adm%gpnt2cross(igpt)%f  = 1.0d0
               enddo
               cycle   
            endif
            
            icrs1 = icrsBeg
            icrs2 = icrsBeg

            xBeg = network%crs%cross(crossOrder(icrsBeg))%chainage
            xEnd = network%crs%cross(crossOrder(icrsEnd))%chainage
            
            do m = 1, pbran%gridPointsCount

               chainageg = pbran%gridPointschainages(m)
               igpt = pbran%grd(m)
               
               if (icrsBeg == icrsEnd) then
                  
                  ! Just one Cross-Section
                  adm%gpnt2cross(igpt)%c1 = crossOrder(icrsBeg)
                  adm%gpnt2cross(igpt)%c2 = crossOrder(icrsBeg)
                  adm%gpnt2cross(igpt)%f  = 1.0d0
                  interpolDone            = .true.   
                  
               elseif (chainageg <= xBeg) then
                  
                  ! Before First Cross-Section
                  adm%gpnt2cross(igpt)%c1 = crossOrder(icrsBeg)
                  adm%gpnt2cross(igpt)%c2 = crossOrder(icrsBeg)
                  adm%gpnt2cross(igpt)%f  = 1.0d0
                  interpolDone            = .true.   
                  
               elseif (chainageg >= xEnd) then
                  
                  ! After Last Cross-Section
                  adm%gpnt2cross(igpt)%c1 = crossOrder(icrsEnd)
                  adm%gpnt2cross(igpt)%c2 = crossOrder(icrsEnd)
                  adm%gpnt2cross(igpt)%f  = 1.0d0
                  interpolDone            = .true.   
                  
               else
                  
                  chainage1 = network%crs%cross(crossOrder(icrs1))%chainage
                  chainage2 = network%crs%cross(crossOrder(icrs2))%chainage
                  
                  if (.not. ((chainage1 <= chainageg) .and. (chainage2 >= chainageg))) then
                     
                     do i = icrs1, icrsEnd
                        if (network%crs%cross(crossOrder(i))%chainage >= chainageg) then
                           chainage2 = network%crs%cross(crossOrder(i))%chainage
                           icrs2 = i
                           exit
                        endif
                     enddo
                     
                     do i = icrsEnd, icrsBeg, -1
                        if (network%crs%cross(crossOrder(i))%chainage <= chainageg) then
                           chainage1 = network%crs%cross(crossOrder(i))%chainage
                           icrs1 = i
                           exit
                        endif
                     enddo
                     
                  endif
                  
                  interpolDone = .false.   

               endif
                  
               ! Interpolation data for Grid Point
               if (igpt > 0) then
                  if (ibran == network%crs%cross(crossOrder(icrs2))%branchid) then
                     
                     if (.not. interpolDone) then
                        adm%gpnt2cross(igpt)%c1 = crossOrder(icrs1)
                        adm%gpnt2cross(igpt)%c2 = crossOrder(icrs2)
                        if (icrs1 == icrs2) then 
                           f = 1.0d0
                        else    
                           if (chainage1 == chainage2) then 
                               write(chainageString, '(F10.3)') chainage1 
                               call setmessage(LEVEL_ERROR, 'Mulitple cross sections defined at same chainage ('// trim(chainageString) //') on branch '//trim(pbran%id))
                               initError = .true.
                           endif
                           f = (chainageg - chainage1) / (chainage2 - chainage1)
                        endif    
                        f = max(f, 0.0d0) 
                        f = min(f, 1.0d0) 
                        adm%gpnt2cross(igpt)%f = f
                     endif
                     
                  else
                     adm%gpnt2cross(igpt)%c1 = crossOrder(icrs1)
                     adm%gpnt2cross(igpt)%c2 = crossOrder(icrs1)
                     adm%gpnt2cross(igpt)%f  = 1.0d0
                  endif
               endif

            enddo

         enddo
         
         deallocate(crossOrder)
         deallocate(lastAtBran)
      
      endif
      
      network%adm = adm
      
      if (initError) then 
          call setmessage(LEVEL_FATAL, 'Error initialising network')
      endif 
      
   end subroutine initialize_1dadmin

   subroutine reassign_pointers(network)

      use m_culvert
   
      type (t_network) :: network
   
      !type (t_branch), pointer :: pbr
      type(t_CrossSection), pointer :: pcs
      type(t_structure), pointer :: pstru
      type(t_culvert), pointer :: pcul
      integer i
   
      !do i = 1, network%brs%count
      !   pbr => network%brs%branch(i)
      !   pbr 
      !enddo
   
      do i = 1, network%crs%Count
         pcs => network%crs%cross(i)
         pcs%tabdef => network%CSDefinitions%cs(pcs%iTabDef)
      enddo
   
      do i = 1, network%sts%Count
         pstru => network%sts%struct(i)
         select case (pstru%type) 
         case (ST_CULVERT) 
            pcul => pstru%culvert
            pcul%pcross => network%crs%cross(pcul%crosssectionnr)
         end select
      enddo
   
   end subroutine reassign_pointers

   subroutine set_network_pointers(network)
      ! Modules
   
      implicit none
      ! Input/output parameters
      type(t_network), intent(inout) :: network
   
      ! Local variables
      type(t_CrossSection), pointer :: pcross
      integer :: istru
      integer :: icrs
   
      ! Program code
      do istru = 1, network%sts%count
         if (network%sts%struct(istru)%type==ST_CULVERT) then
            icrs = network%sts%struct(istru)%culvert%crosssectionnr
            pcross => network%crs%cross(icrs)
            network%sts%struct(istru)%culvert%pcross => pcross
         endif
      enddo
   
   end subroutine set_network_pointers

   !> In this subroutine arrays crossorder and lastAtBran are filled \n
   !! crossorder contains the cross section indices, where the cross sections are ordered 
   !! in ascending branchid and subsequently in chainage.\n
   !! crossorder(lastAtBran(ibr-1)+1) .. crossorder(lastAtBran(ibr)) contain the cross
   !! sections on branch ibr, in ascending branch chainage.
   !! sections on branch ibr, in ascending branch chainage.
   subroutine crossSectionsSort(crs, brs, crossOrder, lastAtBran)

      ! Ordering crs's on branches and x on branches
   
      implicit none
      
      ! Global Variables
      type(t_CrossSectionSet), intent(in)       :: crs
      type(t_branchSet), intent(in)             :: brs
      integer, intent(out),dimension(:)         :: crossOrder
      integer, intent(out), dimension(:)        :: lastAtBran
   
      ! Local Variables
      logical                                   :: changed
      integer                                   :: icrsn
      integer                                   :: ihulp
      integer                                   :: ilast
      integer                                   :: ibran
      integer                                   :: ifirst
      integer                                   :: crossCount
      integer                                   :: i

      do icrsn = 1, crs%count
          crossOrder(icrsn) = icrsn
      enddo

       ! Sort available crsn's in crossOrder on branch number
      crossCount = crs%count
      do i = crs%count, 1, -1
         
         if (.not. crs%cross(i)%branchid > 0) then
            
            ihulp                  = crossOrder(i)
            crossOrder(i)          = crossOrder(crossCount)
            crossOrder(crossCount) = ihulp
            crossCount             = crossCount - 1
            
         endif
         
      enddo

      changed = .true.
      
      do while (changed)
         
         changed = .false.
         
         do icrsn = 1, crossCount - 1
            
            if (crs%cross(crossOrder(icrsn))%branchid .gt. &
                crs%cross(crossOrder(icrsn + 1))%branchid) then
               
               ihulp                 = crossOrder(icrsn)
               crossOrder(icrsn)     = crossOrder(icrsn + 1)
               crossOrder(icrsn + 1) = ihulp
               changed               = .true.
               
            endif
                
         enddo
         
      enddo

!     fill temporal administration array lastAtBran
!     last index in crossOrder for branch(i)

      do ibran=1, brs%count !nbran
         lastAtBran(ibran) = 0
      enddo

      if (crossCount .gt. 0) then
         
         ilast = crs%cross(crossOrder(1))%branchid
         
         do icrsn = 2, crossCount
            if (ilast .ne. crs%cross(crossOrder(icrsn))%branchid ) then
              lastAtBran(ilast) = icrsn - 1
           endif
           ilast = crs%cross(crossOrder(icrsn))%branchid
         enddo
         
         lastAtBran(ilast) = crossCount
         
      endif

!     sort available crsn's in crossOrder on distance in branch

      do ibran = 1, brs%count !nbran

         ! Find range for sorting procedure per branch

         if ( ibran .eq. 1 ) then
            ifirst = 1
         else
            if (lastAtBran(ibran - 1) == 0) then
               lastAtBran(ibran-1) = lastAtBran(max(1,ibran-2))
            endif
            ifirst = lastAtBran(ibran - 1) + 1
         endif
         
         ilast = lastAtBran(ibran)

         ! Sorting procedure per branch

         changed = .true.
          
         do while (changed)
             
            changed = .false.
              
            do icrsn = ifirst, ilast - 1
                 
               if (crs%cross(crossOrder(icrsn))%chainage .gt.  &
                   crs%cross(crossOrder(icrsn + 1))%chainage) then
                  
                  ihulp               = crossOrder(icrsn)
                  crossOrder(icrsn)   = crossOrder(icrsn+1)
                  crossOrder(icrsn+1) = ihulp
                  changed             = .true.
                  
               endif
                     
            enddo
            
         enddo

      enddo

   end subroutine crossSectionsSort
   
   function getRetentionId(network, gridpoint) result(id)
   
      use m_Storage
      use m_node
   
      character(len=80)              :: id
      type(t_network), intent(in)    :: network
      integer, intent(in)            :: gridpoint
      
      integer :: i
      type(t_storage), pointer, dimension(:) :: stor
      
      stor => network%stors%stor
      
      do i = 1, network%stors%count
         if (stor(i)%gridPoint == gridpoint) then
            id = stor(i)%id
            return
         endif
      enddo

      id = 'No retention area assigned to node '//getnodeid(network%nds, gridpoint)
      
   end function getRetentionId


!> Get friction value (Chezy) for a specific point in the network
!> This function is moved from Roughness.f90 to Network.f90 to avoid circular references
double precision function getFrictionValue(rgs, spdata, cross, ibranch, section, igrid, h, q, u, r, d, chainage)

use m_tables
use m_tablematrices
!!--description-----------------------------------------------------------------
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    !=======================================================================
    !                       Deltares
    !                One-Two Dimensional Modelling System
    !                           S O B E K
    !
    ! Subsystem:          Flow Module
    !
    ! Programmer:
    !
    ! Function:           getFrictionValue, replacement of old FLCHZT (FLow CHeZy Friction coeff)
    !
    ! Module description: Chezy coefficient is computed for a certain gridpoint
    !
    !
    !
    !     update information
    !     person                    date
    !     Kuipers                   5-9-2001
    !     Van Putten                11-8-2011
    !
    !     Use stored table counters
    !
    !
    !
    !
    !     Declaration of Parameters:
    !

    implicit none
!
! Global variables
!
    type t_parr
       type(t_CSType), pointer        :: cross
    end type
    

    type(t_RoughnessSet), intent(in)        :: rgs         !< Roughness data
    type(t_spatial_dataSet), intent(in)     :: spData      !< spatial data
    type(t_CSType), pointer, intent(in)     :: cross
    integer, intent(in)                     :: igrid       !< gridpoint index
    integer, intent(in)                     :: ibranch     !< branch index
    integer, intent(in)                     :: section     !< section index (0=main, 1=Flood plane 1, 2=Flood plane 2)
    double precision, intent(in)            :: d           !< water depth
    double precision, intent(in)            :: h           !< water level
    double precision, intent(in)            :: q           !< discharge
    double precision, intent(in)            :: r           !< hydraulic radius
    double precision, intent(in)            :: u           !< velocity
    double precision, intent(in)            :: chainage    !< chainage (location on branch)
    
!
!
! Local variables
!
    integer                         :: isec, i
    double precision                :: cpar
    double precision                :: cz
    double precision                :: dep
    double precision                :: ys
    double precision                :: rad
    type(t_Roughness), pointer      :: rgh 
    type(t_spatial_data), pointer   :: values
    integer, dimension(:), pointer  :: rgh_type
    integer, dimension(:), pointer  :: fun_type

    !     Explanation:
    !     -----------
    !
    !     1. Each Chezy formula, apart from Engelund bed friction, is defined
    !        by 1 constant parameter. This constant is stored in bfricp.
    !        An exception is the Engelund bed friction defined by 10 parameters.
    !     2. For the Engelund bed friction the specific parameters are stored
    !        in the array engpar.
    !
    !
    !     Prevention against zero hydraulic radius and depth
    !
    rad = max(r, 1.d-6)
    dep = max(d, 1.d-6)
    !
    if (associated(cross)) then
       ! section refers to the roughness section *within* the cross section. Actual friction section index then comes from lookup.
       if (cross%frictionSectionsCount > 0) then
         isec = cross%frictionsectionIndex(min(cross%frictionSectionsCount, section))
       else
          isec = section
       endif
    else
       ! No cross section definition: section directly refers to a friction section index.
       isec = section
    endif
    
    if (rgs%version == 1) then
      rgh => rgs%rough(isec)
      if (q >= 0d0 .or. .not. associated(rgh%rgh_type_neg)) then
          values    => spData%quant(rgh%spd_pos_idx)
          rgh_type  => rgh%rgh_type_pos 
          fun_type  => rgh%fun_type_pos 
       else 
          values    => spData%quant(rgh%spd_neg_idx)
          rgh_type  => rgh%rgh_type_neg 
          fun_type  => rgh%fun_type_neg 
       endif   
       if (fun_type(ibranch) == R_FunctionDischarge) then
          cpar = interpolate(values%tables%tb(values%tblIndex(igrid))%table,  dabs(q))
       !
       !        Roughness function of water level depending on flow direction
       !
       elseif (fun_type(ibranch) == R_FunctionLevel) then
          cpar = interpolate(values%tables%tb(values%tblIndex(igrid))%table,  h)
       !
       !        Roughness constant depending on flow direction
       !
       else
          cpar = values%values(igrid)
       endif
       getFrictionValue = GetChezy(rgh_type(ibranch), cpar, rad, dep, u)    !
    !     Formulation = .not. Engelund
    !
    else ! Version 2 roughness
       do i = 1, 2
          if (isec < 0) then
             ! Current cross section does *not* refer to a friction section index, but has defined direct roughness type+coefficient.
             cz = GetChezy(cross%frictionType(section), cross%frictionValue(section), rad, dep, u)
          else
             rgh => rgs%rough(isec)
             if (rgh%useGlobalFriction)then
                cz = GetChezy(rgh%frictionType, rgh%frictionValue, rad, dep, u)
             else
                ! For now, direction independent, always *_pos values.
                rgh_type  => rgh%rgh_type_pos 
                fun_type  => rgh%fun_type_pos 
                if (rgh_type(ibranch) == -1)  then ! This branch has no own roughness definition, use global.
                   cz = GetChezy(rgh%frictionType, rgh%frictionValue, rad, dep, u)
                else
                   if (fun_type(ibranch) == R_FunctionDischarge) then
                      ys = abs(q)
                   elseif (fun_type(ibranch) == R_FunctionLevel) then
                      ys = h
                   else
                      ys = 0d0
                   endif
       
                   cpar = interpolate(rgh%table(ibranch), chainage, ys)
                   cz = GetChezy(rgh_type(ibranch), cpar, rad, dep, u)
                endif
             endif
          endif
       enddo
       getFrictionValue = cz
    endif
    
    
end function getFrictionValue

end module m_network
