module m_sedtrails_network
   use m_sedtrails_data
   use m_alloc
   
   implicit none
   
   contains
   
   subroutine default_sedtrails_geom()
   ! JRE Clean this after throwing out all 1D stuff
      
       ! Remaining of variables is handled in reset_flowgeom()
       call reset_sedtrails_geom()
   end subroutine default_sedtrails_geom
   
   subroutine reset_sedtrails_geom()
       ! node (s) related : dim=numk
       numk=0
   end subroutine
   
   subroutine sedtrails_loadNetwork(filename, istat, jadoorladen)
       use m_sedtrails_netcdf, only : sedtrails_unc_read_net
       use unstruc_messages
       use m_missing
   
       implicit none
   
   
       character(*), intent(in)  :: filename !< Name of file to be read (in current directory or with full path).
       integer,      intent(out) :: istat    !< Return status (0=success)
       integer,      intent(in)  :: jadoorladen
       character(len=255)        :: data_file_1d
   
       integer                   :: iDumk
       integer                   :: iDuml
      
       integer                   :: minp,  K0, L0, L, NUMKN, NUMLN
       logical                   :: jawel
      
       inquire(file = filename, exist=jawel)
       if (.not. jawel) then
           call mess(LEVEL_WARN,'sedtrails_loadNetwork::Could not open '''//trim(filename)//'''')
           return
       end if
   
       IF (JADOORLADEN == 0) THEN
           K0 = 0
       ELSE
           K0 = NUMK
       ENDIF
   
       ! New NetCDF net file
       call sedtrails_unc_read_net(filename, K0, L0, NUMKN, NUMLN, istat)
      
       if (istat == 0) then
           NUMK = K0 + NUMKN
       else
          call qnerror('sedtrails_loadNetwork::Error while loading network from '''//trim(filename)//''', please inspect the preceding diagnostic output.', ' ',  ' ')
       endif
   end subroutine sedtrails_loadNetwork  
   
   !> Increase the number of net links
   SUBROUTINE sedtrails_increasenetwork(K0,L0)
   use m_alloc
   use m_missing, only : xymis, dmiss

   implicit none
   integer,           intent(in) :: K0       !< New number of net nodes.
   integer,           intent(in) :: L0       !< New number of net links

   integer :: ierr
   integer :: k
   integer :: knxx

   if (K0 < KMAX) RETURN

   CALL sedtrails_SAVENET()

   IF (KMAX <= K0) THEN
      KMAX = K0 + 100000   ! 2 KAN WEG
      IF (allocated(xk)) then
         deallocate(xk, yk)
      end if
      ALLOCATE ( XK (KMAX), YK (KMAX), STAT=IERR   )
      CALL AERR('XK (KMAX), YK (KMAX)', IERR, 2*KMAX)

      XK = XYMIS ; YK = XYMIS
   ENDIF

   CALL sedtrails_RESTORE()

END SUBROUTINE 
   
   !> Restore variables with backup data
SUBROUTINE sedtrails_RESTORE()
   use m_sedtrails_data
   implicit none
   integer :: k, KX, LS, LS0, LX, NODSIZ, IERR

   if (.not. allocated(xk0)) return

   KX = size(XK0) ! restore everything present (in case numk/numk0 has not yet been increased)

   XK (1:KX)  = XK0 (1:KX)
   YK (1:KX)  = YK0 (1:KX)

   NUMK = NUMK0

   RETURN
END SUBROUTINE   

 SUBROUTINE sedtrails_SAVENET()
   use m_sedtrails_data
   implicit none
   integer :: ierr
   integer :: k, KX, LS, LS0, LX, NN

   if (.not. allocated(xk)) return

   KX = KMAX ! backup everything present (in case numk has not yet been increased) ! KX = NUMK

   if (allocated(xk0)) deallocate(xk0,yk0,zk0)
   allocate ( XK0(KX), YK0(KX), ZK0(KX) , STAT=IERR)


   XK0 (1:KX) = XK (1:kx)
   YK0 (1:KX) = YK (1:kx)
   ZK0 (1:KX) = ZK (1:kx)

   NUMK0 = NUMK

   RETURN
 END SUBROUTINE
 
 ! Set mask to determine which sedtrails XK,YK points lie on present grid
 subroutine sedtrails_get_grid_on_network()
    use m_sedtrails_data
    use m_polygon
    use m_tpoly
    use m_partitioninfo
    use network_data, only: netstat, NETSTAT_OK
    use geometry_module, only: get_startend
    use m_missing
    use m_flowexternalforcings, only: transformcoef
    use m_flowgeom, only: xz, yz,ndx
    use m_ec_triangle, only: jagetwf, indxx, wfxx
    use m_ec_basic_interpolation, only: triintfast
    use m_sferic
    
    implicit none
    
    integer                               :: k
    integer                               :: ierr, netstat_store
    integer                               :: ipoint, ipoly, numpoints
    integer                               :: istart, iend
    integer                               :: nv, iorient,iinterior,inside
    integer                               :: jakdtree, jdla
    integer, allocatable                  :: indices(:)
    integer, allocatable                  :: sedtrails_idom(:)
    double precision, allocatable         :: dumin(:), dumout(:)
    type(tpoly),dimension(:), allocatable :: pli

    ! Detect grid enclosure for this partition/overlapping part of grids
    call savepol()
    if (jampi>0) then 
       netstat_store = netstat
       netstat = NETSTAT_OK
       call generate_partition_pol_from_idomain(ierr, myrank=my_rank)
       netstat = netstat_store
    else
       call copynetboundstopol(0, 0, 1, 0)
    endif
    call realloc(iistart, maxpoly, keepExisting=.false.)
    call realloc(iiend, maxpoly, keepExisting=.false.)
    ipoint = 1
    ipoly = 0
    numpoints = 0
    do while ( ipoint <= NPL)
       ipoly = ipoly+1
       if (ipoly > maxpoly) then
          maxpoly = ceiling(maxpoly*1.1)
          call realloc(iistart, maxpoly, keepExisting=.true.)
          call realloc(iiend, maxpoly, keepExisting=.true.)
       end if

      ! get polygon start and end pointer respectively
      call get_startend(NPL-ipoint+1,xpl(ipoint:NPL),ypl(ipoint:NPL), istart, iend, dmiss)
      istart = istart+ipoint-1
      iend   = iend  +ipoint-1

      if ( istart.ge.iend .or. iend.gt.NPL ) exit ! done
      
      iistart(ipoly) = istart
      iiend(ipoly)   = iend
      numpoints = numpoints + (iend-istart+1)

!     advance pointer
      ipoint = iend+2
    end do
    npoly = ipoly
    
    !
    ! Allocate poly index
    if (.not. allocated(sedtrails_idom)) then
       allocate(sedtrails_idom(1:numk))
       sedtrails_idom=0
    endif   
    !
    ! Convert to tpolies
    call pol_to_tpoly(npoly, pli, .false.)
    !
    ! Get sedtrials points inside domain
    inside=-1
    jins=1
    do k=1,numk
       call dbpinpol_tpolies(pli, xk(k), yk(k), inside)   ! takes into account inner pols
       if (inside==1) then
          sedtrails_idom(k)=1   
       endif   
    enddo
    !
    ! Get own nodes
    indices=find_nodes_idom_int(sedtrails_idom, 1)
    numk=size(indices)
    !
    ! Reallocate nodes arrays and copy values
    call realloc(xk1,size(xk),keepExisting=.false.,fill=0d0)  
    call realloc(yk1,size(xk),keepExisting=.false.,fill=0d0) 
    xk1=xk
    yk1=yk
    call realloc(xk,numk,keepExisting=.false.,fill=0d0)  
    call realloc(yk,numk,keepExisting=.false.,fill=0d0)   
    xk=xk1(indices)
    yk=yk1(indices)
    !
    ! Generate interpolation weights from flowgeom cell centres
    ! Use dummy interpolation
    ! Save in module variables st_ind, st_wf
    jagetwf = 1
    jakdtree = 1
    allocate ( indxx(3,numk), wfxx(3,numk), dumin(ndx),dumout(numk) )
    dumin=dmiss
    call TRIINTfast(xz,yz,dumin,ndx,1,xk,yk,dumout,numk,JDLA,jakdtree, jsferic, 1, jins, dmiss, jasfer3D, &
                    (/0d0/),(/0d0/),(/0d0/),transformcoef)
    !
    call realloc(st_ind,(/3,numk/), keepExisting=.false.,fill=0)
    call realloc(st_wf,(/3,numk/), keepExisting=.false.,fill=0d0)
    do k=1, numk
       st_ind(3,k)=indxx(3,k)
       st_wf(3,k)=wfxx(3,k)
    enddo   
    !
    call dealloc_tpoly(pli)
    call restorepol()

 end subroutine
 
 FUNCTION find_nodes_idom_int(array, min) RESULT(indices)
    use precision
    integer, INTENT(IN)  :: array(:)
    integer, INTENT(IN)  :: min
    INTEGER, ALLOCATABLE :: indices(:)
    INTEGER :: ii
    indices = PACK([(ii,ii=1,SIZE(array))], array == min)
 END FUNCTION find_nodes_idom_int

end module