!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2019.                                
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

module m_oned_functions
   use m_missing, only: dmiss
   implicit none
   private

   public set_1d_roughnesses
   public set_1d_indices_in_network
   public save_1d_nrd_vars_in_stm
   public setbobs_1d
   public gridpoint2cross
   public computePump_all_links
   public convert_cross_to_prof
   public set_ground_level_for_1d_nodes
   public set_max_volume_for_1d_nodes
   public updateFreeboard
   public updateTimeWetOnGround
   public updateDepthOnGround
   public updateVolOnGround

   type, public :: t_gridp2cs
      integer :: num_cross_sections
      integer, allocatable, dimension(:) :: cross
   end type

   type(t_gridp2cs), allocatable, dimension(:) :: gridpoint2cross

   contains

   !> IFRCUTP and FRCu are filled, using 1D roughness values from Network structure 
   subroutine set_1d_roughnesses()
      use m_flowgeom
      use m_flow, only: frcu, ifrcutp, frcu_mor
      use unstruc_channel_flow

      implicit none

      ! FRCU and FRCU_MOR should only be used after SETAU - VOL12D. 
      ! Therefore initialise these arrays with a negative value.
      if (network%loaded) then
         where (kcu(1:lnx1d) == 1)
            frcu(1:lnx1d) = -10d0
            ifrcutp(1:lnx1d) = 0
            frcu_mor(1:lnx1d) = -10d0
         end where
      endif

   end subroutine set_1d_roughnesses

   !> Sets the flowgeom link and node numbers of the computational grid 
   !! into the 1D network structure for branches, storage nodes, 
   !! cross sections and structures, etc. 
   subroutine set_1d_indices_in_network()
      use m_sediment
      use m_flowgeom
      use m_flow
      use m_cross_helper
      use m_flowparameters
      use unstruc_channel_flow
      
      
      implicit none
      
      default_width = wu1DUNI
      
      if (network%loaded) then
         ! nonlinear computation is required for 1d flow
         if (nonlin1D == 0) then
            nonLin1D = 1
         elseif (nonlin1D == 2) then
            CSCalculationOption = CS_TYPE_PLUS
         endif
         
         nonlin = max(nonlin, nonlin1D)
      endif
      
      if (.not. network%initialized) then
         call set_linknumbers_in_branches()
         call set_node_numbers_for_storage_nodes()
         call set_structure_grid_numbers()
      
         if (jased > 0 .and. stm_included) then
            ! 
            call set_cross_sections_to_gridpoints()
         endif
         call set_structure_indices()
         
         network%initialized = .true.         
      endif
      
   end subroutine set_1d_indices_in_network

   !> set the flowgeom linknumbers and node numbers in the branches
   subroutine set_linknumbers_in_branches()
   
      use m_globalParameters
      use unstruc_channel_flow
      use m_flowgeom
      use m_sediment
      use messageHandling

      implicit none

      integer :: L
      integer :: ibr
      integer :: nbr, upointscount, pointscount
      integer :: storageCount
      integer :: i, j, jpos, linkcount
      integer :: k1, k2, igrid
      integer :: c1, c2
      integer :: storage_count
      type(t_branch), pointer                 :: pbr
      type(t_storage), pointer                :: pstor
      integer, dimension(:), pointer          :: lin
      integer, dimension(:), pointer          :: grd
      double precision, dimension(:), pointer :: chainage
      type(t_chainage2cross), pointer           :: gpnt2cross(:)                   !< list containing cross section indices per u-location
      type (t_CrossSection), pointer          :: cross1, cross2

      nbr = network%brs%count
      do ibr = 1, nbr
         pbr => network%brs%branch(ibr)
         call realloc(pbr%lin, pbr%uPointsCount)
         call realloc(pbr%grd, pbr%gridPointsCount)
         lin => pbr%lin
         grd => pbr%grd
         L = lin(1)
         k1  =  iabs(ln(1,L))
         pbr%FromNode%gridNumber = k1
         upointscount = pbr%uPointsCount
         do i = 1, uPointsCount
            L = lin(i)
            k1 = iabs(ln(1,L))
            grd(i) = k1
         enddo
         k2 = ln(2,iabs(lin(upointscount)))
         pbr%tonode%gridnumber = k2
         grd(upointscount+1) = k2
      enddo
   end subroutine set_linknumbers_in_branches

   !> Set the node numbers from flowgeom for the storage nodes
   subroutine set_node_numbers_for_storage_nodes()
   
      use unstruc_channel_flow
      use m_flowgeom
      use m_sediment
      use messageHandling
      use m_GlobalParameters, only: INDTP_ALL

      implicit none

      integer :: L
      integer :: ibr
      integer :: nbr, upointscount, pointscount
      integer :: storageCount
      integer :: i, j, jpos, linkcount
      integer :: k1, k2, igrid
      integer :: c1, c2
      integer :: storage_count
      type(t_branch), pointer                 :: pbr
      type(t_storage), pointer                :: pstor
      integer, allocatable                    :: ixy2stor(:), k_tmp(:)
      double precision, allocatable           :: x_tmp(:), y_tmp(:)
      character(len=IdLen), allocatable       :: name_tmp(:)
      integer                                 :: nxy, countxy, jakdtree

      
      countxy = network%storS%Count_xy
      if (countxy > 0) then
         call realloc(ixy2stor,    countxy, keepExisting=.false.)
         call realloc(k_tmp,       countxy, keepExisting=.false.)
         call realloc(x_tmp,       countxy, keepExisting=.false.)
         call realloc(y_tmp,       countxy, keepExisting=.false.)
         call realloc(name_tmp,    countxy, keepExisting=.false.)
      end if
      
      nxy = 0
      do i = 1, network%storS%count
         pstor => network%storS%stor(i)
         if (pstor%node_index < 0) then
            nxy = nxy + 1
            ixy2stor(nxy) = i
            x_tmp(nxy)    = pstor%x
            y_tmp(nxy)    = pstor%y
            name_tmp(nxy) = pstor%id
         else
            pStor%gridPoint = network%nds%node(pstor%node_index)%gridNumber
         end if
      end do
      
      if (nxy > 0) then ! find flow nodes for storage nodes that are defined by x-, y-coordinates
         jakdtree = 1
         call find_flownode(nxy, x_tmp(1:nxy), y_tmp(1:nxy), name_tmp(1:nxy), k_tmp(1:nxy), jakdtree, 0, INDTP_1D)
         do i = 1, nxy
            if (k_tmp(i) > 0) then
               pstor => network%storS%stor(ixy2stor(i))
               pstor%gridPoint = k_tmp(i)
            else
               call SetMessage(LEVEL_ERROR, 'Error when snapping storage node '''//trim(name_tmp(i))//''' to a flow node.')
            end if
         end do
      
         if (allocated(k_tmp))    deallocate(k_tmp)
         if (allocated(x_tmp))    deallocate(x_tmp)
         if (allocated(y_tmp))    deallocate(y_tmp)
         if (allocated(ixy2stor)) deallocate(ixy2stor)
         if (allocated(name_tmp)) deallocate(name_tmp)
      end if
      
      
   end subroutine set_node_numbers_for_storage_nodes
   
   subroutine set_structure_grid_numbers()
      use unstruc_channel_flow
      use m_flowgeom
      use m_flowexternalforcings
      use m_inquire_flowgeom

      implicit none
    
      integer :: istru, local_index, nstru, ierr
      type(t_structure), pointer :: pstru
      type(t_branch), pointer :: pbranch

      
      nstru = network%sts%count
      if (nstru>0) then
         call realloc(L1strucsg, nstru)
         call realloc(L2strucsg, nstru)
      endif
      
      
   end subroutine set_structure_grid_numbers
   
   !> For sediment transport on each node a cross section is required
   !! Fills gridpoint2cross with for each gridpoint a cross section index. \n
   !! Note: On connection nodes we have multiple cross sections (one for each 
   !!       incoming or outgoing branch (link). \n
   !!       A connection node is located at the beginning or end of the branch.
   subroutine set_cross_sections_to_gridpoints()
   
      use unstruc_channel_flow
      use m_flowgeom
      use m_sediment
      use messageHandling

      implicit none

      integer :: L
      integer :: ibr
      integer :: nbr, upointscount, pointscount
      integer :: storageCount
      integer :: i, j, jpos, linkcount
      integer :: k1, k2, igrid
      integer :: c1, c2
      integer :: storage_count
      double precision :: d1, d2, dh
      type(t_branch), pointer                 :: pbr
      type(t_storage), pointer                :: pstor
      integer, dimension(:), pointer          :: lin
      integer, dimension(:), pointer          :: grd
      double precision, dimension(:), pointer :: chainage
      type(t_chainage2cross), pointer           :: gpnt2cross(:)                   !< list containing cross section indices per u-location
      type (t_CrossSection), pointer          :: cross1, cross2


      ! cross sections (in case of sediment transport every gridpoint requires a unique
      ! cross section)
      if (jased > 0 .and. stm_included) then
         if (allocated(gridpoint2cross)) deallocate(gridpoint2cross)
         allocate(gridpoint2cross(ndxi))
         gpnt2cross => network%adm%gpnt2cross
         do i = 1, ndxi
            gridpoint2cross(i)%num_cross_sections = 0
         enddo

         ! allocate space for local cross section numbers on connection nodes (multiple cross sections)
         do i = 1, network%nds%count
            k1 = network%nds%node(i)%gridNumber
            linkcount = nd(k1)%lnx
            if (allocated(gridpoint2cross(k1)%cross)) deallocate(gridpoint2cross(k1)%cross)
            allocate(gridpoint2cross(k1)%cross(linkcount))
            gridpoint2cross(k1)%num_cross_sections = linkcount
            gridpoint2cross(k1)%cross = -999
         enddo
         
         igrid = 0
         nbr = network%brs%count
         do ibr = 1, nbr
            pbr => network%brs%branch(ibr)
            lin => pbr%lin
            grd => pbr%grd
            chainage => pbr%gridPointschainages
            pointscount = pbr%gridPointsCount
            do i = 1, pointscount
               igrid = igrid+1
               k1 = grd(i)
               if (i==1 .or. i==pointscount) then
                  ! search for correct location
                  ! this entry (gridpoint2cross(k1)) is already allocated
                  if (i==1) then 
                     L = lin(1)
                     dh = (chainage(i+1)-chainage(i))/2d0
                  else
                     L = lin(pointscount-1)
                     dh = (chainage(i)-chainage(i-1))/2d0
                  endif
                  do j = 1,nd(k1)%lnx
                     if (L == iabs(nd(k1)%ln(j))) then
                        jpos = j
                     endif
                  enddo
               else
                  ! Internal gridpoint on branch, only 1 cross section attached
                  if (allocated(gridpoint2cross(k1)%cross)) deallocate(gridpoint2cross(k1)%cross)
                  allocate(gridpoint2cross(k1)%cross(1))
                  gridpoint2cross(k1)%num_cross_sections = 1
                  jpos = 1
                  dh = min(chainage(i)-chainage(i-1),chainage(i+1)-chainage(i))/2d0
               endif
               c1 = gpnt2cross(igrid)%c1
               c2 = gpnt2cross(igrid)%c2
               d1 = abs(network%crs%cross(c1)%chainage - chainage(i))
               d2 = abs(network%crs%cross(c2)%chainage - chainage(i))
               ! cross1%branchid and cross2%branchid should correspond to ibr
               if (d1 < dh) then
                  gridpoint2cross(k1)%cross(jpos) = c1
               elseif (d2 < dh) then
                  gridpoint2cross(k1)%cross(jpos) = c2
               else
                  gridpoint2cross(k1)%cross(jpos) = -999
               endif
            enddo
         enddo
      endif
   end subroutine set_cross_sections_to_gridpoints
      
   ! function to store variables related to the nodal relation variables
   subroutine save_1d_nrd_vars_in_stm
      use m_branch
      use m_node
      use m_sediment, only: sedtra, stmpar
      use unstruc_channel_flow
      use morphology_data_module, only : t_nodefraction, t_noderelation
      use string_module

      implicit none

      integer :: inod, ibr, iFrac, iNodeRel, directionLink = 0
      type(t_branch), pointer :: pbr
      type(t_node)  , pointer :: pnod
      type(t_nodefraction), pointer          :: pFrac
      type(t_noderelation),pointer           :: pNodRel

      if (network%loaded) then
          do iFrac = 1, stmpar%nrd%nFractions
              pFrac => stmpar%nrd%nodefractions(iFrac)
              do iNodeRel = 1, pFrac%nNodeRelations
                  pNodRel => pFrac%noderelations(iNodeRel)
                  do ibr = 1, network%brs%Count
                      pbr => network%brs%branch(ibr)
                      if (pNodRel%node == pbr%fromNode%id) then
                          pNodRel%nodeIdx = pbr%fromNode%gridnumber
                      endif
                      if (pNodRel%node == pbr%toNode%id) then
                          pNodRel%nodeIdx = pbr%toNode%gridnumber
                      endif
                      if (pNodRel%BranchIn == pbr%id) then
                          if (pNodRel%node == pbr%fromNode%id) then
                              pNodRel%BranchInLn = pbr%lin(1)                    ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchInLn = pbr%lin(pbr%uPointsCount)     ! (positive = at end of branch)
                          endif
                      endif
                      if (pNodRel%BranchOut1 == pbr%id) then
                          if (pNodRel%node == pbr%fromNode%id) then
                              pNodRel%BranchOut1Ln = pbr%lin(1)                  ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchOut1Ln = pbr%lin(pbr%uPointsCount)   ! (positive = at end of branch)
                          endif
                      endif
                      if (pNodRel%BranchOut2 == pbr%id) then
                          if (pNodRel%node == pbr%fromNode%id) then
                              pNodRel%BranchOut2Ln = pbr%lin(1)                  ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchOut2Ln = pbr%lin(pbr%uPointsCount)   ! (positive = at end of branch)
                          endif
                      endif
                  enddo
              enddo
          enddo
      endif

   end subroutine save_1d_nrd_vars_in_stm

   !> 
   subroutine set_structure_indices()
   end subroutine set_structure_indices

   subroutine setbobs_1d()
   
   use m_network
   use m_flowgeom
   use m_flowtimes
   use messagehandling
   use unstruc_messages
   use unstruc_channel_flow
   use m_1d_structures
   use m_cross_helper
   use network_data
   
   implicit none
   
   integer :: i
   integer :: L, L0
   integer :: n1
   integer :: n2
   integer :: nstor
   integer :: nnode
   integer :: nstruc
   double precision :: crest_level
   type(t_structure), pointer :: pstruc
   type(t_storage),   pointer :: pstor
   
   do i = ndx2D+1, ndxi
      bl(i) = huge(1d0)
   enddo
   
   
   nstor = network%storS%count
   do i = 1, nstor
      pstor => network%storS%stor(i)
      n1 = pstor%gridPoint
      bl(n1) = min(bl(n1), pstor%storageArea%x(1))
   enddo
      
   do L = 1, lnx1D
      if (kcu(L) ==1) then
         bob(:,L)  = getbobs(network, L)
         bob0(:,L) = bob(:,L)
         n1  = ln(1,L)
         n2 = ln(2,L)                    ! flow ref
         bl(n1) = min(bl(n1), bob(1,L))
         bl(n2) = min(bl(n2), bob(2,L))
      endif
   enddo

   ! In case of compound structures bob is set to the lowest crest level
   ! Pumping stations get the bob of the channel
   ! First step is to initialise all structure bobs to huge
   nstruc = network%sts%count
   do i = 1, nstruc
      pstruc => network%sts%struct(i)
      do L0 = 1, pstruc%numlinks
         L = pstruc%linknumbers(L0)
         bob(:,L) = huge(1d0)
      enddo
   enddo
   
   nstruc = network%sts%count
   do i = 1, nstruc
      pstruc => network%sts%struct(i)
      crest_level = get_crest_level(pstruc)
      do L0 = 1, pstruc%numlinks
         L = pstruc%linknumbers(L0)
         if (crest_level < huge(1d0)) then
            bob(1,L) = min(bob(1,L), crest_level)
            bob(2,L) = min(bob(2,L), crest_level)
         else
            ! pumping station
            bob(1,L) = min(bob(1,L), bob0(1,L))
            bob(2,L) = min(bob(2,L), bob0(2,L))
         endif
      enddo
   enddo
   
   if (time_user<= tstart_user) then
      ! check if all manholes are lower than or equal to the invert level of all incoming pipes
      nstor = network%storS%count
      do i = 1, nstor
         pstor => network%storS%stor(i)
         n1 = pstor%gridPoint
         if (bl(n1) < pstor%storageArea%x(1)) then
            call setmessage(LEVEL_WARN, 'At node '//trim(network%nds%node(i)%id)//' the bedlevel is below the bedlevel of the assigned storage area.')
            write(msgbuf, '(''The bedlevel (due to invert levels of incoming channels/pipes) = '', g14.2, '' and the bottom level of the storage area is '', g14.2)') &
                        bl(n1), pstor%storageArea%x(1)
            call setmessage(-LEVEL_WARN, msgbuf)
         
         endif
      
      enddo
   endif
   

   do i = ndx2D+1, ndxi
      if (bl(i) > 0.5d0*huge(1d0)) then
         write(msgbuf, '(a,i0,a)') 'Bedlevel is missing on calculation flow node ', i, '. No nearby cross sections nor storage nodes.'
         call warn_flush()
         bl(i) = zkuni
      end if
   enddo

   ! look for missing bobs
   do L = 1, lnx1d
      if (bob(1,L) > 0.5d0*huge(1d0)) then
         bob(1,L)  = bl(ln(1,L))
         bob0(1,L) = bob(1,L)
      endif
      if (bob(2,L) > 0.5d0*huge(1d0)) then
         bob(2,L)  = bl(ln(2,L))
         bob0(2,L) = bob(2,L)
      endif
   enddo
   
   do L = lnxi+1, lnx1Db
       ! mirror 1d bed level points at boundary 
       n1 = ln(1,L)
       n2 = ln(2,L)
       bl(n1) = bl(n2)
       bob(1,L)  = bl(n1)
       bob(2,L)  = bl(n2)       
       bob0(1,L) = bl(n1)
       bob0(2,L) = bl(n2)       
   enddo    
   
   end subroutine setbobs_1d

   !> Compute FU and RU coefficients for each flow link that is part of 
   !! the pump. Values are stored in struct%fu(:), etc. and *also* set
   !! in m_flow::fu(:), etc.
   subroutine computePump_all_links(struct)
      use m_1d_structures
      use m_pump
      use m_flowtimes
      use m_flowgeom 
      use m_flow
      
      type(t_structure), intent(inout) :: struct !< The parent structure of the pump (which also contains the flow link information).
                
      double precision     :: s1k1
      double precision     :: s1k2
      double precision     :: qp
      double precision     :: ap
      double precision     :: vp1, vp2, vp
      integer              :: L   
      integer              :: L0   
      integer              :: k1   
      integer              :: k2
      integer              :: dir
      integer              :: n
      
      ! First compute average waterlevels on suction side and delivery side of the pump
      s1k1 = 0d0
      s1k2 = 0d0
      ap = 0d0
      vp1 = 0d0
      vp2 = 0d0
      vp = 0d0
      qp = 0d0
      do L0 = 1, struct%numlinks
         L = struct%linknumbers(L0)
         ! Note: Link L may have negative sign if flow link is opposite pump's orientation
         ! (pump spatial orientation is polyline+righthand rule, or network branch direction).
         ! Note 2: do not account for pumping direction here, that is done in prepareComputePump.
         dir = sign(1, L) ! only includes flow link w.r.t. structure spatial orientation.
         L = iabs(L)
         if ( dir > 0) then         
            k1 = ln(1,L)
            k2 = ln(2,L)
         else
            k1 = ln(2,L)
            k2 = ln(1,L)
         endif
         
         if (hs(k1) > 1d-2) then
            ! NOTE: pump area-weighting across links is uniform for all links (au=1).
            au(L) = 1d0
            ap    = ap + au(L)
            vp1    = vp1 + vol1(k1)
            vp2    = vp2 + vol1(k2)
            s1k1 = s1k1 + au(L)*s1(k1)
            s1k2 = s1k2 + au(L)*s1(k2)
         endif
      enddo
      
      ! With these average waterlevels, evaluate the pump discharge.
      if (ap > 0d0) then
         s1k1 = s1k1/ap
         s1k2 = s1k2/ap
         call PrepareComputePump(struct%pump, s1k1, s1k2)
         qp    = struct%pump%discharge ! Already in our local structure spatial orientation.
         
         ! Choose available volume on suction side.
         if (qp > 0d0) then
            vp = vp1
         else
            vp = vp2
         end if
      endif          

      ! Finally, redistribute the requested pump discharge across all flow links.
      if (qp == 0d0 .or. ap == 0 .or. vp == 0d0) then
         ! Pump is off
         struct%fu = 0d0
         struct%ru = 0d0
         struct%au = 0d0
      else

         ! Limit the pump discharge in case the volume in the cells at the suction side is limited.
          if (abs(qp) > 0.9d0*vp/dts) then
            qp = sign(qp,0.9d0*vp/dts)
         endif
         
         do L0  = 1, struct%numlinks
            L = struct%linknumbers(L0)
            dir = int(sign(1d0, L*qp)) ! Includes both pumping direction and flow link w.r.t. structure spatial orientation.
            L = iabs(L)
            if ( dir > 0) then         
               k1 = ln(1,L)
            else
               k1 = ln(2,L)
            endif
         
            if (hs(k1) > 1d-2) then
               struct%fu(L0) =  0d0
               struct%ru(L0) =  qp/ap
               struct%au(L0) =  ap
            else 
               struct%fu(L0) = 0d0
               struct%ru(L0) = 0d0
               struct%au(L0) = 0d0
            endif
         enddo
      endif
      
      do L0  = 1, struct%numlinks
         L = iabs(struct%linknumbers(L0))
         fu(L) = struct%fu(L0)
         ru(L) = struct%ru(L0)
         au(L) = struct%au(L0)
      enddo

   end subroutine computePump_all_links


   !> Converts the currently active cross sections in network%crs to
   !! old-format profloc/profdef input.
   subroutine convert_cross_to_prof(basename)
   use unstruc_channel_flow
   use m_inquire_flowgeom
   use m_samples
   use m_polygon
   use m_missing
   use m_flowgeom
   
   character(len=*), intent(in) :: basename !< Basename for the profdef/loc output files.

   integer :: ic, it, nprof, LF, ierr, mloc, mdef, mxyz, nyz, numxyztype
   integer, allocatable :: cs2prof(:)
   type(t_CrossSection), pointer :: pcrs
   type(t_CSType), pointer :: pcs
   character(len=255) :: proflocfile, profdeffile, profdefxyzfile

   call savesam()
   call savepol()

   proflocfile = trim(basename)//'_profloc.xyz'
   profdeffile = trim(basename)//'_profdef.txt'
   profdefxyzfile = trim(basename)//'_profdefxyz.pliz'

   !defs%CS(pcross%iTabDef)
   call newfil(mdef, profdeffile)
   call newfil(mxyz, profdefxyzfile)

   allocate(cs2prof(network%CSDefinitions%Count))

   nprof = 0
   npl = 0
   numxyztype = 0
   call realloc(nampli, network%CSDefinitions%Count, fill=' ')
   do it=1,network%CSDefinitions%Count
      pcs => network%CSDefinitions%CS(it)
      select case (pcs%crossType)
      case(CS_YZ_PROF)
         nprof = nprof+1
         cs2prof(it) = nprof

         numxyztype = numxyztype+1
         ! First write definition one-liner
         write (mdef, '(a,i0,a)') 'PROFNR=', nprof, '     TYPE=201' ! FRCTP=1 FRCCF=.035

         ! Then write xyz definition pliz
         write (nampli(numxyztype), '(a,i0)') 'PROFNR=', nprof
         nyz = pcs%levelsCount
         xpl(npl+1:npl+nyz) = 0d0
         ypl(npl+1:npl+nyz) = pcs%y(1:nyz)
         zpl(npl+1:npl+nyz) = pcs%z(1:nyz)
         npl = npl+nyz+1
         xpl(npl) = dmiss; ypl(npl) = dmiss; zpl(npl) = dmiss ! Separator between pli/csdef
      case default
         call QNERROR('Error in convert_cross_to_prof(), profile type not supported:',CSTypeName(pcs%crossType),'')
      end select
   end do

   call wripol(mxyz)


   call doclose(mxyz)
   call doclose(mdef)

   call increasesam(network%crs%Count)
   NS = 0
   do ic=1,network%crs%Count
      pcrs => network%crs%cross(ic)
      
      ierr = findlink(pcrs%branchid, pcrs%chainage, Lf)
      if (Lf > 0) then
         NS = NS+1
         XS(NS) = xu(Lf)
         YS(NS) = yu(Lf)
         ZS(NS) = cs2prof(pcrs%iTabDef)
      end if
   end do
   
   call newfil(mloc, proflocfile)
   call wrisam(mloc)
   call doclose(mloc)
   
   !call restoresam()
   !call restorepol()
   
   end subroutine convert_cross_to_prof
   
   !> Set ground level for 1d nodes.
   !! Ground level should not be confused with bed level.
   !! It is defined as:
   !! * street level for storage nodes that have street storage (storageType is reservior or closed),
   !! * the highest nearby cross section level ("embankment") for other nodes,
   !! * dmiss, i.e. not applicable, if no cross section is defined at the node.
   subroutine set_ground_level_for_1d_nodes(network)
   use m_flowgeom, only: groundLevel, ndxi, ndx2d
   use m_Storage
   use m_CrossSections
   use m_network
   implicit none
   type(t_network), intent(inout), target :: network
   type(t_storage), pointer               :: pSto
   type(t_administration_1d), pointer     :: adm
   integer                                :: i, cc1, cc2

   groundlevel(:) = dmiss

   ! set for all 1D nodes, the ground level equals to the highest cross section "embankment" value
   do i = 1, ndxi-ndx2d
      adm => network%adm
      cc1 = adm%gpnt2cross(i)%c1
      cc2 = adm%gpnt2cross(i)%c2
      if (cc1 > 0 .and. cc2 > 0) then ! if there are defined cross sections
         groundLevel(i) = getHighest1dLevel(network%crs%cross(cc1), network%crs%cross(cc2), adm%gpnt2cross(i)%f)
      end if
   end do

   ! set for storage nodes that have prescribed street level, i.e. storageType is reservior and closed
   do i = 1, network%storS%Count
      pSto => network%storS%stor(i)
      if (pSto%useStreetStorage .and. (.not. pSto%useTable)) then
         groundLevel(pSto%gridPoint) = pSto%streetArea%x(1)
      end if
   end do

   end subroutine set_ground_level_for_1d_nodes
   
   !> Set maximal volume for 1d nodes
   subroutine set_max_volume_for_1d_nodes()
   use m_flowgeom, only: groundLevel, volMax, ndx, ndxi, ndx2d
   use m_flow,     only: s1, vol1, a1, vol1_f, a1m
   implicit none
   double precision, allocatable :: s1_tmp(:), vol1_tmp(:), a1_tmp(:), vol1_ftmp(:), a1m_tmp(:)
   integer                       :: ndx1d
   
   ! 1. copy curernt s1, vol1, vol1_f, a1 and a1m to a temporary array
   if (allocated(s1_tmp)) then
      deallocate(s1_tmp)
   end if
   allocate(s1_tmp(ndx))
   s1_tmp = s1
   
   if (allocated(vol1_tmp)) then
      deallocate(vol1_tmp)
   end if
   allocate(vol1_tmp(ndx))
   vol1_tmp = vol1
   
   if (allocated(vol1_ftmp)) then
      deallocate(vol1_ftmp)
   end if
   allocate(vol1_ftmp(ndx))
   vol1_ftmp = vol1_f
   
   if (allocated(a1_tmp)) then
      deallocate(a1_tmp)
   end if
   allocate(a1_tmp(ndx))
   a1_tmp = a1
   
   if (allocated(a1m_tmp)) then
     deallocate(a1m_tmp)
   end if
   allocate(a1m_tmp(ndx))
   a1m_tmp = a1

   ! 2. set s1 to be the ground level
   ndx1d = ndxi-ndx2d
   s1(ndx2d+1:ndxi) = groundLevel(1:ndx1d)
   vol1   = 0d0
   vol1_f = 0d0
   a1     = 0d0
   a1m    = 0d0
   

   ! 3. compute the maximal volume
   call vol12d(0)
   volMax(1:ndx1d) = vol1(ndx2d+1:ndxi)

   ! 4. set s1, vol1, vol1_f, a1, a1m back
   s1 = s1_tmp
   vol1 = vol1_tmp
   vol1_f = vol1_ftmp
   a1 = a1_tmp
   a1m = a1m_tmp
   
   if (allocated(s1_tmp))   deallocate(s1_tmp)
   if (allocated(vol1_tmp)) deallocate(vol1_tmp)
   if (allocated(a1_tmp))   deallocate(a1_tmp)
   if (allocated(vol1_ftmp))deallocate(vol1_ftmp)
   if (allocated(a1m_tmp))  deallocate(a1m_tmp)
   end subroutine set_max_volume_for_1d_nodes


   !> Update freeboard for each 1d node.
   !! Freeboard is the vertical distance between the ground level, i.e., not bed level, and the water surface.
   !! It can be negative value if the water level is above the ground level
   !! It has maximal value 0d0 if it has storageType "closed", or if the relevant cross sections are closed.
   subroutine updateFreeboard(network)
   use m_flow, only: freeboard, s1
   use m_flowgeom, only: ndxi, ndx2d, groundLevel
   use m_network
   implicit none
   type(t_network), intent(inout), target :: network
   type(t_administration_1d), pointer     :: adm
   integer :: i, ii, cc1, cc2

   do i = ndx2d+1, ndxi
      ii = i- ndx2d
      if (groundLevel(ii) .ne. dmiss) then ! if ground level is applicable
         adm => network%adm
         cc1 = adm%gpnt2cross(ii)%c1
         cc2 = adm%gpnt2cross(ii)%c2
         if (network%crs%cross(cc1)%closed .and. network%crs%cross(cc2)%closed) then
            freeboard(i) = max(0d0, groundLevel(ii) - s1(i))
         else
            freeboard(i) = groundLevel(ii) - s1(i)
         end if
      end if
   end do
   
   ! if the storage node has storageType "closed", then freeboard's maximal value is 0d0
   if (network%storS%Count_closed > 0) then
      call updateValueForClosedNodes(ndxi, freeboard, 0d0, network, 1)
   end if

   end subroutine updateFreeboard


   !> Compute the cumulative time when water is above ground level.
   !! hsOnGround has to be updated before calling this subroutine
   subroutine updateTimeWetOnGround(dts)
   use m_flowparameters, only: epshs
   use m_flowtimes, only: time_wetground
   use m_flow, only: hsOnGround
   use m_flowgeom,only: ndxi, ndx2d
   implicit none
   double precision, intent(in) :: dts !< computational time step
   integer                      :: i
   
   do i = ndx2d+1, ndxi
      if (hsOnGround(i) .ne. dmiss .and. hsOnGround(i) >= epshs) then
         time_wetground(i) = time_wetground(i) + dts
      end if
   end do

   end subroutine updateTimeWetOnGround


   !> Update waterdepth above ground level for each 1d node.
   !! This is the vertical distance between the water surface and the ground level, i.e. waterLevel minus groundLevel.
   !! It has minimal value 0d0
   !! It equals to dmiss if the node has storageType closed, or if the relavent cross sections are closed, or no cross section is defined.
   subroutine updateDepthOnGround(network)
   use m_flow, only: hsOnGround, s1
   use m_network
   use m_flowgeom, only: ndxi, ndx2d, groundLevel
   implicit none
   type(t_network), intent(inout), target :: network
   type(t_administration_1d), pointer     :: adm
   integer                                :: i, ii, cc1, cc2
   
   hsOnGround = dmiss
   do i = ndx2d+1, ndxi
      ii = i-ndx2d
      adm => network%adm
      cc1 = adm%gpnt2cross(ii)%c1
      cc2 = adm%gpnt2cross(ii)%c2
      if (network%crs%cross(cc1)%closed .and. network%crs%cross(cc2)%closed) then
         cycle
      else
         if (groundLevel(ii) .ne. dmiss) then ! if groudLevel is applicable
            hsOnGround(i) = max(0d0, s1(i) - groundLevel(ii))
         end if
      end if
   end do

   ! if the storage node has storageType "closed", then hsOnGround equals to dmiss
   if (network%storS%Count_closed > 0) then
      call updateValueForClosedNodes(ndxi, hsOnGround, dmiss, network, 0)
   end if
   
   end subroutine updateDepthOnGround


   !> Update volume above ground level for each 1d node.
   !! It has minimal value 0d0
   !! It equals to dmiss if the node has storageType closed, or if the relavent cross sections are closed, or no cross section is defined.
   subroutine updateVolOnGround(network)
   use m_flow,     only: vol1, volOnGround
   use m_flowgeom, only: volMax, ndxi, ndx2d, groundLevel
   use m_network
   implicit none
   type(t_network), intent(inout), target :: network
   type(t_administration_1d), pointer     :: adm
   integer:: i, ii, cc1, cc2
   
   volOnGround = dmiss
   do i = ndx2d+1, ndxi
      ii = i-ndx2d
      adm => network%adm
      cc1 = adm%gpnt2cross(ii)%c1
      cc2 = adm%gpnt2cross(ii)%c2
      if (network%crs%cross(cc1)%closed .and. network%crs%cross(cc2)%closed) then
         cycle
      else
         if (groundLevel(ii) .ne. dmiss) then ! if ground level is applicable
            volOnGround(i) = max(0d0, vol1(i) - volMax(ii))
         end if
      end if
   end do
   
   ! if the storage node has storageType "closed", then volOnGround equals to dmiss
   if (network%storS%Count_closed > 0) then
      call updateValueForClosedNodes(ndxi, volOnGround, dmiss, network, 0)
   end if

   end subroutine updateVolOnGround


   !< Update the given array on storage nodes of sotrageType "closed", given the value that is to be set on these nodes 
   !! If updates freeboard array, then set as the maximal value among the given value and the existing value
   subroutine updateValueForClosedNodes(n, array, value, network, jaFreeboard)
   use m_network
   use m_storage
   implicit none
   integer,                        intent(in   ) :: n
   double precision, dimension(n), intent(inout) :: array
   double precision,               intent(in   ) :: value
   type(t_network),                intent(inout) :: network
   integer,                        intent(in   ) :: jaFreeboard
   type(t_storage), pointer                      :: pSto
   integer                                       :: i, ii
   
   do i = 1, network%storS%Count
      pSto => network%storS%stor(i)
      if (pSto%useStreetStorage .and. (.not. pSto%useTable) .and. pSto%storageType == nt_Closed) then
         ii = pSto%gridPoint
         if (jaFreeboard == 0) then
            array(ii) = value
         else
            array(ii) = max(value, array(ii))
         end if
      end if
   end do

   end subroutine updateValueForClosedNodes

end module m_oned_functions
