module m_sedtrails_network
   use m_sedtrails_data
   use m_alloc
   
   implicit none
   
   contains
   
   subroutine default_sedtrails_geom()
      
       ! Remaining of variables is handled in reset_flowgeom()
       call reset_sedtrails_geom()
   end subroutine default_sedtrails_geom
   
   subroutine reset_sedtrails_geom()
       ! node (s) related : dim=numk
       numk=0
   end subroutine
      
   !> Increase the number of net links
   SUBROUTINE sedtrails_increasenetwork(K0)
   use m_alloc
   use m_missing, only : xymis, dmiss

   implicit none
   integer,           intent(in) :: K0       !< New number of net nodes.

   integer :: ierr
   integer :: k
   integer :: knxx

   if (K0 < KMAX) RETURN

   CALL sedtrails_SAVENET()

   IF (KMAX <= K0) THEN
      KMAX = K0 + 100000
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
   integer :: k, KX
   
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
    use m_partitioninfo, only: my_rank, jampi,generate_partition_pol_from_idomain
    use network_data, only: netstat, NETSTAT_OK
    use geometry_module, only: get_startend, dbdistance
    use m_missing
    use m_flowexternalforcings, only: transformcoef
    use m_flowgeom, only: xz, yz,ndx,bl
    use m_ec_triangle, only: jagetwf, indxx, wfxx
    use m_ec_basic_interpolation, only: triinterp2
    use m_sferic
    
    implicit none
    
    integer                               :: k
    integer                               :: ierr, netstat_store
    integer                               :: ipoint, ipoly, numpoints
    integer                               :: istart, iend
    integer                               :: nv, iorient,iinterior,inside
    integer                               :: jakdtree, jdla, ip1
    double precision                      :: dmaxsize
    
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
    call dealloc_tpoly(pli)
    call restorepol()
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
    deallocate(xk1, yk1)
    !
    ! Generate interpolation weights from flowgeom cell centres
    ! Use dummy interpolation
    ! Save in module variables st_ind, st_wf
    jagetwf = 1
    jakdtree = 1
    jdla=1
    call realloc(indxx,(/ 3,numk /),keepExisting=.false., fill=0)
    call realloc(wfxx,(/ 3,numk /),keepExisting=.false., fill=0d0)
    call realloc(dumout,numk,keepExisting=.false., fill=dmiss)
 
    transformcoef(6)=1.1d0
    CALL triinterp2(xk, yk, dumout, numk, jdla, &
            xz, yz, bl, ndx, dmiss, jsferic, jins, jasfer3D, NPL, 0, 0, XPL, YPL, ZPL, transformcoef)
    !
    call realloc(st_ind,(/3,numk/), keepExisting=.false.,fill=0)
    call realloc(st_wf,(/3,numk/), keepExisting=.false.,fill=0d0)
    do k=1, numk
       st_ind(:,k)=indxx(:,k)
       st_wf(:,k)=wfxx(:,k)
    enddo   
    !
    ! And now that we have the correct number of nodes:
    if (jampi>0) then
       call realloc(idomain,numk,keepExisting=.false.,fill=my_rank)
    endif   
    !
    deallocate (indxx, wfxx, dumout)

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