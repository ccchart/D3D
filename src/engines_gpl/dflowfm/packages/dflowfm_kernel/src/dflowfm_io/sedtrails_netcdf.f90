module m_sedtrails_netcdf
   use m_sedtrails_data
   
   implicit none
   
   contains
   
   !> Increase the number of net links
   SUBROUTINE sedtrails_increasenetwork(K0,L0, also_dxe)
   use m_alloc
   use m_missing, only : xymis, dmiss

   implicit none
   integer,           intent(in) :: K0       !< New number of net nodes.
   integer,           intent(in) :: L0       !< New number of net links
   logical, optional, intent(in) :: also_dxe !< Also allocate the optional dxe array for edge lengths. Default: .false.

   integer :: ierr
   integer :: k
   integer :: knxx
   logical :: also_dxe_

   if (present(also_dxe)) then
      also_dxe_ = also_dxe
   else
      also_dxe_ = .false.
   end if

   if (also_dxe_) then
      ! Directly ensure dxe allocation, since it may not have been allocated before (as it is optional).
      call realloc(dxe, LMAX, keepExisting = .true., fill = dmiss)
   end if

   if (K0 < KMAX .and. L0 < LMAX) RETURN

   CALL sedtrails_SAVENET()

   IF (KMAX <= K0) THEN
      KMAX = K0 + 100000   ! 2 KAN WEG
      IF (allocated(nod)) then
         do k = 1,size(nod)
            if ( allocated(nod(k)%lin) ) deallocate (nod(k)%lin)
         enddo
         deallocate(nod, xk, yk, zk, kc, nmk)
      end if
      ALLOCATE ( NOD(KMAX) , STAT = IERR)
      CALL AERR('NOD(KMAX)', IERR, KMAX )
      ALLOCATE ( XK (KMAX), YK (KMAX), ZK (KMAX), KC (KMAX), NMK (KMAX) , STAT=IERR   )
      CALL AERR('XK (KMAX), YK (KMAX), ZK (KMAX), KC (KMAX), NMK (KMAX)', IERR, 6*KMAX)

      DO K = 1,KMAX
         IF ((allocated(NMK0)).and.(K .LE. SIZE(NMK0))) THEN
            KNXX = MAX(NMK0(K),KNX)
         ELSE
            KNXX = KNX
         ENDIF
         ALLOCATE(NOD(K)%LIN(KNXX) , STAT=IERR) ;
         NOD(K)%LIN = 0
      ENDDO

      NMK = 0 ; KC = 1 ; XK = XYMIS ; YK = XYMIS ; ZK = dmiss
   ENDIF

   IF (LMAX <= L0) THEN
      LMAX = L0 + 3*100000
      IF (SIZE(LC) > 0 .and. allocated(kn)) THEN
         DEALLOCATE(KN ,LC , RLIN)
      ENDIF
      ALLOCATE (KN (3,LMAX), LC (LMAX), STAT=IERR) ; KN = 0 ; LC = 0 ! TODO: AvD: catch memory error
      ALLOCATE (RLIN (LMAX), STAT=IERR)
      
      if (also_dxe_) then
         if (allocated(dxe)) deallocate(dxe)
         allocate(dxe(LMAX))
         dxe = dmiss
      end if
   ENDIF

   CALL sedtrails_RESTORE()

END SUBROUTINE 
   
   !> Restore variables with backup data
SUBROUTINE sedtrails_RESTORE()
   use m_sedtrails_data
   implicit none
   integer :: k, KX, LS, LS0, LX, NODSIZ, IERR

   !IF ( NUMK0.EQ.0 ) RETURN
   if (.not. allocated(xk0) .or. .not. allocated(kn0) .or. .not. allocated(nod0)) return

   KX = size(XK0) ! restore everything present (in case numk/numk0 has not yet been increased)

   XK (1:KX)  = XK0 (1:KX)
   YK (1:KX)  = YK0 (1:KX)
   ZK (1:KX)  = ZK0 (1:KX)

   NMK(1:KX)  = NMK0(1:KX)
   KC (1:KX)  = KC0 (1:KX)

   LX = size(LC0) ! restore everything present (in case numl/numl0 has not yet been increased)

   KN(:,1:LX) = KN0(:,1:LX)
   LC(  1:LX) = LC0(  1:LX)

   ! Only restore optional dxe when it is there already
   if (allocated(dxe)) then
      dxe(1:LX) = dxe0(1:LX)
   end if

   NODSIZ = SIZE(NOD)

   DO K = 1,KX
      IF (ALLOCATED(NOD0(K)%LIN) ) THEN
         LS0 = SIZE(NOD0(K)%LIN )  ! LS0 = NMK0(K)
      ELSE
         LS0 = 0
      ENDIF

      IF (LS0 .GE. 1) THEN
         ! IF (.NOT. ASSOCIATED(NOD(K)%LIN) ) THEN
         IF (ALLOCATED(NOD(K)%LIN) ) THEN
            LS = SIZE(NOD(K)%LIN )
         ELSE
            LS = 0
         ENDIF
         IF (LS .LT. LS0) THEN
            IF (LS .GE. 1) DEALLOCATE (NOD(K)%LIN )
            ALLOCATE   (NOD(K)%LIN(LS0), STAT = IERR )
            NOD(K)%LIN = 0
         ENDIF
         NOD(K)%LIN(1:LS0) = NOD0(K)%LIN(1:LS0)
      ENDIF
   ENDDO

   NUMK = NUMK0
   NUML = NUML0

   !  need findcells
   netstat = NETSTAT_CELLS_DIRTY
   RETURN
   END SUBROUTINE   

!> Reads the net data from a NetCDF file.
!! Processing is done elsewhere.
subroutine sedtrails_unc_read_net_ugrid(filename, numk_keep, numl_keep, numk_read, numl_read, ierr)
   use m_sedtrails_data
   use m_save_ugrid_state
   use io_netcdf
   use netcdf
   use netcdf_utils, only: ncu_get_att
   use m_sferic
   use m_missing
   use unstruc_messages
   use MessageHandling
   use dfm_error
   use m_alloc

   character(len=*), intent(in)    :: filename           !< Name of NetCDF file.
   integer,          intent(inout) :: numk_keep          !< Number of netnodes to keep in existing net (0 to replace all).
   integer,          intent(inout) :: numl_keep          !< Number of netlinks to keep in existing net (0 to replace all).
   integer,          intent(out)   :: numk_read          !< Number of new netnodes read from file.
   integer,          intent(out)   :: numl_read          !< Number of new netlinks read from file.
   integer,          intent(out)   :: ierr               !< Return status (NetCDF operations)

   integer                         :: ioncid, iconvtype, start_index, networkIndex
   integer                         :: im, nmesh, i, L, numk_last, numl_last
   integer                         :: ncid, id_netnodez
   integer, allocatable            :: kn12(:,:), kn3(:) ! Placeholder arrays for the edge_nodes and edge_types
   double precision                :: convversion, zk_fillvalue, altsign
   type(t_ug_meshgeom) :: meshgeom

   character(len=:), allocatable             :: tmpstring
   integer :: n1, n2, ibr_n1, ibr_n2, ibr
   double precision :: off1, off2
   integer :: numerr, threshold_abort_current
   logical :: need_edgelengths, includeArrays

   numk_read = 0
   numl_read = 0
   start_index = 1
   numk_last = numk_keep
   numl_last = numl_keep
   includeArrays = .true.
   networkIndex = 0
   numerr = 0
   
   allocate(character(len=0) :: tmpstring)
   ierr = ionc_open(filename, NF90_NOWRITE, ioncid, iconvtype, convversion)

   if (ierr /= ionc_noerr .or. iconvtype /= IONC_CONV_UGRID .or. convversion < 1.0) then ! NOTE: no check on conventions version number (yet?)
      ! No valid UGRID, not a problem, call site will fall back to trying old format.
      call mess(LEVEL_ERROR,  'sedtrails_unc_read_net_ugrid:: net file '''//trim(filename)//''' is not UGRID. Save network file with cell info.')
      ierr = DFM_EFILEFORMAT
      goto 999
   end if
 
   ierr = ionc_get_ncid(ioncid, ncid)
   tmpstring = ''
   ierr = ncu_get_att(ncid, nf90_global, 'Conventions', tmpstring)
   if (ierr == NF90_ENOTATT) then
      call mess(LEVEL_DEBUG,  'sedtrails_unc_read_net_ugrid::No NetCDF Conventions found. Defaulting to current format (>= "CF-1.8 UGRID-1.0 Deltares-0.10") for '''//trim(filename)//'''.')
   end if
   deallocate(tmpstring)

   ierr = ionc_get_coordinate_reference_system(ioncid, crs)
   ! ierr = ionc_get_crs(ioncid, crs) ! TODO: make this API routine.
   ! TODO: also get the %crs item from the ionc dataset, store it in unstruc, AND, use that one in unc_write_flowgeom_ugrid.
   if (ierr /= ionc_noerr) then
   call mess(LEVEL_WARN,  'sedtrails_unc_read_net_ugrid::ionc_get_coordinate_system: No epsg_code found in UGRID net file '''//trim(filename)//'''.')
   goto 999
   end if
   select case (crs%epsg_code)
   case (4326) ! WGS84
      jsferic  = 1
   case default
      jsferic  = 0
      jasfer3D = 0
   end select
   !
   ! Prepare for multiple (partial) meshes
   !
   ierr = ionc_get_mesh_count(ioncid, nmesh)
   if (ierr /= ionc_noerr) then
      call mess(LEVEL_WARN,  'sedtrails_unc_read_net_ugrid::: No grids found in UGRID net file '''//trim(filename)//'''.')
      goto 999
   end if
   
   !------------------------------------------------------------!
   ! meshes
   !------------------------------------------------------------!
   do im = 1, nmesh
      
      ierr = ionc_get_meshgeom(ioncid, im, networkIndex, meshgeom)
      
      if (meshgeom%dim == 2) then
         !Else 2d/3d mesh
         if (meshgeom%numnode < 0 .or. meshgeom%numface < 0) then
            cycle
         end if
         ierr = ionc_get_meshgeom(ioncid, im, networkIndex, meshgeom, start_index, includeArrays) 
         mesh2dname = meshgeom%meshname
         !ierr = ionc_get_face_coordinates(ioncid, im, xface, yface)
         !call read_mesh2d_face_z(ioncid, im, meshgeom%numface)
      else
         ! Only support 1D network and 2D grid
         write(msgbuf, '(a,i0,a,i0,a)') 'sedtrails_unc_read_net_ugrid: unsupported topology dimension ', meshgeom%dim, &
            ' in file '''//trim(filename)//' for mesh #', im, '.'
         call warn_flush()
         cycle
      end if
      
      !do_edgelengths = .false.
      need_edgelengths = .false. ! Either from a previous meshgeom, or now for the first time.
      
      !increasenetw 
      call sedtrails_increasenetwork(numk_last + meshgeom%numnode, numl_last + meshgeom%numedge, also_dxe=need_edgelengths) ! increases XK, YK, KN, optionally dxe
      if (meshgeom%dim == 2) then  ! always true, see above
      ! 2D, or 1D without network topology
         ierr = ionc_get_node_coordinates(ioncid, im, XK(numk_last+1:numk_last + meshgeom%numnode), YK(numk_last+1:numk_last + meshgeom%numnode)) ! TODO: LC: this duplicates with the above get_meshgeom with includearrays=.true.
      endif

      if (ierr /= ionc_noerr) then
         write (msgbuf, '(a,i0,a)') 'unc_read_net_ugrid: Could not read x/y node coordinates from mesh #', im, ' in UGRID net file '''//trim(filename)//'''.'
         call warn_flush()
         goto 999
      end if

      ierr = ionc_get_ncid(ioncid, ncid)
      if (ierr /= ionc_noerr) then
         write (msgbuf, '(a,i0,a)') 'unc_read_net_ugrid: Could not get direct access to UGRID NetCDF net file '''//trim(filename)//'''.'
         call warn_flush()
         goto 999
      end if

      ! zk values on nodes, not needed
      ZK(numk_last+1:numk_last+meshgeom%numnode) = dmiss

      !
      ! 3. Net links. Just append the edges from the mesh(es) as netlinks, later setnodadm() at call site will group them by 1D and 2D.
      !
      if (allocated(kn12)) deallocate(kn12)
      allocate(kn12(2, meshgeom%numedge))

      if (allocated(kn3))  deallocate(kn3)
      allocate(kn3(meshgeom%numedge))

      if (meshgeom%dim.ne.1) then 
         ! TODO: LC: these have already been read into meshgeom%edge_nodes, so maybe just copy it here?
         ierr = ionc_get_edge_nodes(ioncid, im, kn12, 1) !unstruct requires 1 based indexes
      else
         kn12 = meshgeom%edge_nodes
         ierr = ionc_noerr
      endif

      if (ierr /= ionc_noerr) then
         write (msgbuf, '(a,i0,a)') 'sedtrails_unc_read_net_ugrid: Could not read edge-node connectivity from mesh #', im, ' in UGRID net file '''//trim(filename)//'''.'
         call warn_flush()
         goto 999
      end if
      
      kn3(:) = 2

      do L=1,meshgeom%numedge
         ! Append the netlink table, and also increment netnode numbers in netlink array to ensure unique ids.
         kn(1:2,numl_last+L) = numk_last + kn12(:,L)
         kn(3,  numl_last+L) = kn3(L)
      enddo
      
      numk_read = numk_read + meshgeom%numnode 
      numk_last = numk_last + meshgeom%numnode  

      numl_read = numl_read + meshgeom%numedge
      numl_last = numl_last + meshgeom%numedge

   end do

   ! Success
888 continue    
   ierr = ionc_close(ioncid)
   ierr = dfm_noerr
   return

999 continue
   ! Some error occurred (error code previously set)
   ! Try to close+cleanup the data set anyway.
   i = ionc_close(ioncid) ! Don't overwrite actual ierr.

end subroutine sedtrails_unc_read_net_ugrid   
   
   
!> Reads the net data from a NetCDF file.
!! Processing is done elsewhere.
subroutine sedtrails_unc_read_net(filename, numk_keep, numl_keep, numk_read, numl_read, ierr)
    use m_sedtrails_data
    use m_sferic
    use m_missing
    use dfm_error
    use gridoperations
    use netcdf_utils, only: ncu_get_att, ncu_get_var_attset

    character(len=*), intent(in)     :: filename  !< Name of NetCDF file.
    integer,          intent(inout)  :: numk_keep !< Number of netnodes to keep in existing net.
    integer,          intent(inout)  :: numl_keep !< Number of netlinks to keep in existing net.
    integer,          intent(out)    :: numk_read !< Number of new netnodes read from file.
    integer,          intent(out)    :: numl_read !< Number of new netlinks read from file.
    integer,          intent(out)    :: ierr      !< Return status (NetCDF operations)

    logical :: stringsequalinsens
    
    character(len=:), allocatable :: coordsyscheck
    integer, dimension(:),   allocatable :: kn3read
    integer, dimension(:),   allocatable :: kn1read
    integer, dimension(:),   allocatable :: kn2read
    
    
    integer :: inetfile, &
               id_netnodedim, id_netlinkdim, &             !< Dimensions
               id_netnodex, id_netnodey, id_netnodez, &    ! Node variables
               id_netlink, id_netlinktype, &                !< Link variables
               id_crsvar

    integer :: L
    double precision :: zk_fillvalue

    call mess(LEVEL_INFO,'sedtrails_unc_read_net::Reading net data...')
    !
    ! Try and read as new UGRID NetCDF format
    !
    call sedtrails_unc_read_net_ugrid(filename, numk_keep, numl_keep, numk_read, numl_read, ierr)
    if (ierr == dfm_noerr) then
       ! UGRID successfully read, we're done.
       return
    else
       ! No UGRID, but just try to use the 'old' format now.
       call mess(LEVEL_ERROR,'sedtrails_unc_read_net::Could not read '''//trim(filename)//'''')
       return
    end if

end subroutine


   SUBROUTINE sedtrails_SAVENET()
   use m_sedtrails_data
   implicit none
   integer :: ierr
   integer :: k, KX, LS, LS0, LX, NN

   if (.not. allocated(xk) .or. .not. allocated(kn) .or. .not. allocated(nod)) return

   KX = KMAX ! backup everything present (in case numk has not yet been increased) ! KX = NUMK
   IF (ALLOCATED(nod0)) THEN
      DO K= 1, SIZE(NOD0)
         if ( allocated(nod0(k)%lin) ) DEALLOCATE(NOD0(K)%LIN)
      ENDDO
      DEALLOCATE(NOD0)
   ENDIF
   ALLOCATE ( NOD0(KX) , stat = ierr )
   !CALL AERR('NOD0(KX)', IERR, KX)

   if (allocated(xk0)) deallocate(xk0,yk0,zk0)
   allocate ( XK0(KX), YK0(KX), ZK0(KX) , STAT=IERR)
   !call aerr('XK0(KX), YK0(KX), ZK0(KX)', IERR, 3*kx)

   if (allocated (KC0) ) deallocate ( KC0 )
   ALLOCATE( KC0(KX), STAT=IERR)

   if (allocated (nmk0) ) deallocate ( NMK0 )
   ALLOCATE( NMK0(KX), STAT=IERR)

   XK0 (1:KX) = XK (1:kx)
   YK0 (1:KX) = YK (1:kx)
   ZK0 (1:KX) = ZK (1:kx)
   KC0( 1:KX) = KC (1:kx)
   NMK0(1:KX) = NMK(1:kx)

   IF (ALLOCATED(LC0)) DEALLOCATE(KN0 ,LC0)
   LX = LMAX ! backup everything present (in case numl has not yet been increased) ! LX = NUML
   ALLOCATE (KN0(3,LX), LC0(LX), STAT=IERR)

   KN0(:,1:LX) = KN(:,1:LX)
   LC0(  1:LX) = LC(  1:LX)

   ! Only save optional dxe when it is there already
   if (allocated(dxe)) then
      if (allocated(dxe0)) deallocate(dxe0)
      allocate(dxe0(LX), STAT=IERR)
      dxe0(1:LX) = dxe(1:LX)
   end if


   DO K   = 1,KX
      LS  = NMK(K) ! SIZE(NOD (K)%LIN )
      IF (LS .GE. 1) THEN
         !       IF (.NOT. ASSOCIATED(NOD0(K)%LIN) ) THEN
         IF (.NOT. ALLOCATED(NOD0(K)%LIN) ) THEN
            LS0 = 0
         ELSE
            LS0 = SIZE(NOD0(K)%LIN )
         ENDIF
         IF (LS0 .LT. LS) THEN
            IF (LS0 .GE. 1 .and. allocated(NOD0(K)%LIN)) DEALLOCATE (NOD0(K)%LIN )
            ALLOCATE   (NOD0(K)%LIN(LS) ) ; NOD0(K)%LIN = 0
         ENDIF
         NOD0(K)%LIN(1:LS) = NOD(K)%LIN(1:LS)
      ENDIF
   ENDDO

   NUMK0 = NUMK
   NUML0 = NUML
   RETURN
   END SUBROUTINE
  
   
end module