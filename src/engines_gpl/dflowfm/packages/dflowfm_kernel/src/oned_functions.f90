module m_oned_functions

   implicit none
   private

   public set_1d_roughnesses
   public set_1d_indices_in_network
   public save_1d_nrd_vars_in_stm
   public setbobs_1d
   public gridpoint2cross

   type, public :: t_gridp2cs
      integer :: num_cross_sections
      integer, allocatable, dimension(:) :: cross
   end type

   type(t_gridp2cs), allocatable, dimension(:) :: gridpoint2cross

   contains

   !> IFRCUTP and FRCu are filled, using 1D roughness values from Network structure 
   subroutine set_1d_roughnesses()
      use m_flow, only: frcu, ifrcutp
      use unstruc_channel_flow
      use m_spatial_data
      use m_branch
      use m_hash_search

      implicit none

      integer :: L, i, k
      integer :: ibr
      integer :: cross
      integer :: irough

      type(t_branch), pointer                 :: pbr
      double precision, dimension(:), pointer :: cpar
      integer,          dimension(:), pointer :: rgh_type
      integer,          dimension(:), pointer :: fun_type
      integer,          dimension(9)          :: rgh_mapping
      type(t_CrossSection), dimension(:), pointer :: crs
      
      if (network%brs%Count > 0) then
         ! RGH_TYPE is similar to IFRCUTP, only with different type numbers
         ! Dflow1D also supports water level or discharge dependent roughness parameters (FUN_TYPE )
         rgh_mapping = -1
         rgh_mapping(R_Chezy         ) = 0
         rgh_mapping(R_Manning       ) = 1
         rgh_mapping(R_WhiteColebrook) = 3

         crs => network%crs%cross

         do ibr = 1, network%brs%Count
            pbr => network%brs%branch(ibr)
            L = pbr%lin(1)
            cross = network%adm%line2cross(l)%c1
            
            if (cross < 0) then
               !use default
               do i = 1, pbr%uPointsCount
                  L = pbr%lin(i)
                  ifrcutp(L) = 0
                  frcu(L) = 60d0
               enddo
            else
               iRough = hashsearch(network%rgs%hashlist, crs(cross)%frictionSectionID(1))
               rgh_type => network%rgs%rough(irough)%rgh_type_pos
               fun_type => network%rgs%rough(irough)%fun_type_pos
               cpar     => network%spData%quant(network%rgs%rough(irough)%spd_pos_idx)%values
               do i = 1, pbr%uPointsCount
                  L = pbr%lin(i)
                  k = pbr%points(1) -1 + i
                  ifrcutp(L) = rgh_mapping(rgh_type(ibr))
                  if (ifrcutp(L) >=0) then
                     ! R_FunctionConstant, R_FunctionDischarge, R_FunctionLevel are computed as Chezy, frcu(L) is set in getprof_1D ! Why is there another chezy computaton, getcz exists?
                     ifrcutp(L) = 0
                  else
                     call setmessage(LEVEL_FATAL, '1D roughness type on branch '// trim(pbr%name) //', '//pbr%id//' is not available in D-FlowFM')
                     ifrcutp(L) = 0
                     frcu(L)    = 45d0
                  endif
               enddo
            endif
         enddo
      endif

   end subroutine set_1d_roughnesses

   !> Sets the flowgeom link and node numbers of the computational grid 
   !! into the 1D network structure for branches, retention nodes, 
   !! cross sections and structures, etc. 
   subroutine set_1d_indices_in_network()
      use m_sediment
      use m_flowgeom
      use m_cross_helper
      use m_flowparameters
      use unstruc_channel_flow
      
      implicit none
      
      default_width = wu1DUNI
      
      if (network%brs%count > 0) then
         ! nonlinear computation is required for 1d flow
         if (nonlin1D == 0) nonLin1D = 1
         nonlin = max(nonlin, nonlin1D)
      endif
      
      call set_linknumbers_in_branches()
      call set_retention_grid_numbers()
      call set_structure_grid_numbers()
      
      if (jased > 0 .and. stm_included) then
         ! 
         call set_cross_sections_to_gridpoints()
      endif
      call set_structure_indices()
      
   end subroutine set_1d_indices_in_network

   !> set the flowgeom linknumbers and node numbers in the branches
   subroutine set_linknumbers_in_branches()
   
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
      double precision, dimension(:), pointer :: offset
      type(t_offset2cross), pointer           :: gpnt2cross(:)                   !< list containing cross section indices per u-location
      type (t_CrossSection), pointer          :: cross1, cross2

      nbr = network%brs%count
      do ibr = 1, nbr
         pbr => network%brs%branch(ibr)
         lin => pbr%lin
         grd => pbr%grd
         L = lin(1)
         k1  =  ln(1,L)
         pbr%FromNode%gridNumber = k1
         upointscount = pbr%uPointsCount
         do i = 1, uPointsCount
            L = lin(i)
            k1 = ln(1,L)
            grd(i) = k1
         enddo
         k2 = ln(2,lin(upointscount))
         pbr%tonode%gridnumber = k2
         grd(upointscount+1) = k2
      enddo
   end subroutine set_linknumbers_in_branches

   !> Set the node numbers from flowgeom in the retention structure
   subroutine set_retention_grid_numbers()
   
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
      double precision, dimension(:), pointer :: offset
      type(t_offset2cross), pointer           :: gpnt2cross(:)                   !< list containing cross section indices per u-location
      type (t_CrossSection), pointer          :: cross1, cross2


      storageCount = network%storS%count
      do i = 1, storageCount
         pstor => network%storS%stor(i)
         if (pstor%branch_index <= 0) then
            pstor%gridPoint = network%nds%node(pstor%node_index)%gridNumber
         else
            pbr => network%brs%branch(pstor%branch_index)
            pstor%gridPoint = pbr%grd(pstor%local_grid_index)
         endif
      enddo
   end subroutine set_retention_grid_numbers
   
   subroutine set_structure_grid_numbers()
      use unstruc_channel_flow
      use m_flowexternalforcings

      implicit none
    
      integer :: istru, local_index, nstru
      type(t_structure), pointer :: pstru
      type(t_branch), pointer :: pbranch
      
      nstru = network%sts%count
      if (nstru>0) then
         call realloc(L1strucsg, nstru)
         call realloc(L2strucsg, nstru)
      endif
      
      do istru = 1, nstru
         pstru   => network%sts%struct(istru)
         pbranch => network%brs%branch(pstru%ibran)
         local_index = pstru%left_calc_point - pbranch%Points(1)+1
         pstru%link_number      = pbranch%lin(local_index)
         pstru%left_calc_point  = pbranch%grd(local_index)
         pstru%right_calc_point = pbranch%grd(local_index+1)
         L1strucsg(istru) = istru
         L2strucsg(istru) = istru
      enddo
      
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
      double precision, dimension(:), pointer :: offset
      type(t_offset2cross), pointer           :: gpnt2cross(:)                   !< list containing cross section indices per u-location
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
         enddo
         
         igrid = 0
         nbr = network%brs%count
         do ibr = 1, nbr
            pbr => network%brs%branch(ibr)
            lin => pbr%lin
            grd => pbr%grd
            offset => pbr%gridPointsOffsets
            pointscount = pbr%gridPointsCount
            do i = 1, pointscount
               igrid = igrid+1
               k1 = grd(i)
               if (i==1 .or. i==pointscount) then
                  ! search for correct location
                  ! this entry (gridpoint2cross(k1)) is already allocated
                  if (i==1) then 
                     L = lin(1)
                     dh = (offset(i+1)-offset(i))/2d0
                  else
                     L = lin(pointscount-1)
                     dh = (offset(i)-offset(i-1))/2d0
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
                  dh = min(offset(i)-offset(i-1),offset(i+1)-offset(i))/2d0
               endif
               c1 = gpnt2cross(igrid)%c1
               c2 = gpnt2cross(igrid)%c2
               d1 = abs(network%crs%cross(c1)%location - offset(i))
               d2 = abs(network%crs%cross(c2)%location - offset(i))
               ! cross1%branchid and cross2%branchid should correspond to ibr
               if (d1 < dh) then
                  gridpoint2cross(k1)%cross(jpos) = c1
               elseif (d2 < dh) then
                  gridpoint2cross(k1)%cross(jpos) = c2
               else
                  gridpoint2cross(k1)%cross(jpos) = 0
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

      if (network%brs%Count > 0) then
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
                              pNodRel%BranchInLn = pbr%lin(pbr%upoints(1))    ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchInLn = pbr%lin(pbr%upoints(2))     ! (positive = at end of branch)
                          endif
                      endif
                      if (pNodRel%BranchOut1 == pbr%id) then
                          if (pNodRel%node == pbr%fromNode%id) then
                              pNodRel%BranchOut1Ln = pbr%lin(pbr%upoints(1))  ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchOut1Ln = pbr%lin(pbr%upoints(2))   ! (positive = at end of branch)
                          endif
                      endif
                      if (pNodRel%BranchOut2 == pbr%id) then
                          if (pNodRel%node == pbr%fromNode%id) then
                              pNodRel%BranchOut2Ln = pbr%lin(pbr%upoints(1))  ! (negative = at start of branch)
                          elseif (pNodRel%node == pbr%toNode%id) then
                              pNodRel%BranchOut2Ln = pbr%lin(pbr%upoints(2))   ! (positive = at end of branch)
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
   use messagehandling
   use unstruc_messages
   use unstruc_channel_flow
   use m_structure
   use m_cross_helper
   use network_data
   
   implicit none
   
   integer :: i
   integer :: L
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
         bob(:,L) = getbobs(network, L)
         n1  = ln(1,L)
         n2 = ln(2,L)                    ! flow ref
         if (bob(1,L) < bl(n1)) then
            bl(n1) = bob(1,L)
            write(msgbuf, '(f8.4)') bob(1,L)
            msgbuf = 'Bed level of retention area: '//trim(getRetentionId(network, n1))//'. is above adjoining invertlevel of pipe (= '//trim(msgbuf)//')'
            call warn_flush()
         endif
         if (bob(2,L) < bl(n2)) then
            bl(n2) = bob(2,L)
            write(msgbuf, '(f8.4)') bob(2,L)
            msgbuf = 'Bed level of retention area: '//trim(getRetentionId(network, n2))//'. is above adjoining invertlevel of pipe (= '//trim(msgbuf)//')'
            call warn_flush()
         endif
         
         bl(n2) = min(bl(n2), bob(2,L))
      endif
   enddo
      
   nstruc = network%sts%count
   do i = 1, nstruc
      pstruc => network%sts%struct(i)
      crest_level = get_crest_level(pstruc)
      if (crest_level < 0.5*huge(1d0) ) then
         L = pstruc%link_number
         bob(1,L) = crest_level
         bob(2,L) = crest_level
         !n1  = ln(1,L)
         !n2 = ln(2,L)                    ! flow ref
         !bl(n1) = min(bl(n1), bob(1,L))
         !bl(n2) = min(bl(n2), bob(2,L))
      endif
   enddo
   
   ! look for missing bobs
   do L = 1, lnx1d
      if (bob(1,L) > 0.5d0*huge(1d0)) then
         bob(1,L) = bl(ln(1,L))
      endif
      if (bob(2,L) > 0.5d0*huge(1d0)) then
         bob(2,L) = bl(ln(2,L))
      endif
   enddo
   
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

   ! check for missing storage on nodes
   nnode = networK%nds%count
   do i = 1, nnode
      n1 = network%nds%node(i)%gridNumber
      if (bl(n1) > 0.5d0*huge(1d0)) then
         call setmessage(LEVEL_ERROR, 'Storage is missing on node '//trim(network%nds%node(i)%id))
      endif
      
      if (getMaxErrorLevel() > LEVEL_ERROR) then
         call setmessage(LEVEL_FATAL, 'Error(s) occurred in initialisation.')
      endif
   enddo
         
   do i = ndx2D+1, ndxi
      if (bl(i) > 0.5d0*huge(1d0)) then
         bl(i) = zkuni
      endif
   enddo

   do L = lnxi+1, lnx1Db
       ! mirror 1d bed level points at boundary 
       n1 = ln(1,L)
       n2 = ln(2,L)
       bl(n1) = bl(n2)
       bob(1,L) = bl(n1)
       bob(2,L) = bl(n2)       
   enddo    
   
   end subroutine setbobs_1d

end module m_oned_functions
