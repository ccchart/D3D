module m_sedtrails_data
   use coordinate_reference_system
   implicit none
   
   character(len=255)    :: sedtrails_analysis
   !
   ! flow geometry from md_sedtrailsfile
   ! 
   ! Copy of cell_geometry
   integer, target                       :: ndx2d      !< [-] Number of 2D flow cells (= NUMP). {"rank": 0}
   integer, target                       :: ndx        !< [-] Number of flow nodes (internal + boundary). {"rank": 0}
   integer, target                       :: ndxi       !< [-] Number of flow nodes (internal + boundary). {"rank": 0}
   integer, target                       :: ndx1db         !< [-] Number of flow nodes incl. 1D bnds (internal 2D+1D + 1D bnd). {"rank": 0}
   double precision, allocatable, target :: xz (:)     !< [m/degrees_east] waterlevel point / cell centre, x-coordinate (m) {"location": "face", "shape": ["ndx"]}
   double precision, allocatable         :: xz0(:)     !< backup of xz
   double precision, allocatable, target :: yz (:)     !< [m/degrees_north] waterlevel point / cell centre, y-coordinate (m) {"location": "face", "shape": ["ndx"]}
   double precision, allocatable         :: yz0(:)     !< backup of yz
   double precision, allocatable, target :: ba (:)     !< [m2] bottom area, if < 0 use table in node type {"location": "face", "shape": ["ndx"]}
   !
   ! copy of flowgeometry
   ! 1:ndx2D, ndx2D+1:ndxi, ndxi+1:ndx1Db, ndx1Db+1:ndx
   ! ^ 2D int ^ 1D int      ^ 1D bnd       ^ 2D bnd ^ total
   type tnode                                          !< node administration
     integer                         :: lnx            !< max nr of links attached to this node
     integer, allocatable            :: ln (:)         !< linknrs attached to this node, >0: to this flownode, <0: from this flownode
     
     integer, allocatable            :: nod(:)         !< Mapping to net nodes
     double precision, allocatable   :: x  (:)         !< for now, this is only for quick/aligned plotting, the corners of a cell
     double precision, allocatable   :: y  (:)         !< for now, this is only for quick/aligned plotting, the corners of a cell
     integer                         :: nwx            !< nr of walls attached
     integer, allocatable            :: nw (:)         !< wallnrs attached to this node
   end type tnode
   ! 
   type (tnode),     allocatable         :: nd(:)          !< (ndx) flow node administration
   integer,          allocatable         :: kcs(:)         !< node code permanent
   integer,          allocatable         :: kcsini(:)      !< node code during initialization, e.g., for initialwaterlevel1d/2d
   integer,          allocatable, target :: kfs(:)         !< [-] node code flooding {"shape": ["ndx"]}
   double precision, allocatable         :: bai(:)         !< inv bottom area (m2), if < 0 use table in node type  
   double precision, allocatable, target :: bl(:)          !< [m] bottom level (m) (positive upward) {"location": "face", "shape": ["ndx"]}
   !
   integer, target                       :: lnxi           !< [-] nr of flow links (internal, 1D+2D    ). {"rank": 0}
   integer, target                       :: lnx1D          !< [-] nr of flow links (internal, 1D+2D    ). {"rank": 0}
   integer, target                       :: lnx1Db         !< [-] nr of flow links including 1D bnds (internal, 1D+2D, boundary: only 1D. 2D bnd behind it). {"rank": 0}
   integer, target                       :: lnx            !< [-] nr of flow links (internal + boundary). First we have 1D links, next 2D links, next boundary links (first 1D, then 2D). {"rank": 0}
   integer,          allocatable, target :: ln    (:,:)    !< [-] 1D link (2,*) node   administration, 1=nd1,  2=nd2   linker en rechter celnr {"shape": [2, "lnkx"]}
   integer,          allocatable, target :: lncn  (:,:)    !< [-] 2D link (2,*) corner administration, 1=nod1, 2=nod2  linker en rechter netnr {"shape": [2, "lnkx"]}
   double precision, allocatable, target :: dx    (:)      !< [m] link length (m) {"location": "edge", "shape": ["lnx"]}
   double precision, allocatable         :: dxi   (:)      !< inverse dx
   double precision, allocatable, target :: wu(:)          !< [m] link initial width (m), if < 0 pointer to convtab {"location": "edge", "shape": ["lnx"]}
   integer,          allocatable         :: ln2lne(:)      !< flowlink to netlink nr dim = lnx
   integer,          allocatable         :: lne2ln(:)      !< netlink to flowlink nr dim = numL
   !
   type tcorn                                          !< corner administration
     integer                         :: lnx            !< max nr of links attached to this corner
     integer, allocatable            :: ln (:)         !< linknrs attached to this corner
   end type tcorn                                      !< corner administration
   type(tcorn)     , allocatable     :: cn  (:)        !< cell cornerpoints, (in counting order of nod)
   
   ! To do: clean
   double precision                  :: bamin          !< minimum 2D cell area
   double precision                  :: bamin1D        !< minimum cell area 1d nodes
   double precision                  :: dxmin=1d-3     !< minimum link length 1D (m)
   double precision                  :: dxmin1D        !< minimum link length 1D (m)
   double precision                  :: dxmin2D        !< minimum link length 2D (m)
   double precision                  :: dxwuimin2D     !< smallest fraction dx/wu , may increase dx if > 0
   double precision                  :: wu1DUNI        !< uniform 1D profile width
   double precision                  :: hh1DUNI        !< uniform 1D profile height
   
   ! netnode/flownode  related, dim = mxban
   double precision, allocatable     :: banf  (:)     !< horizontal netnode/flownode area (m2)
   double precision, allocatable     :: ban  (:)      !< horizontal netnode          area (m2)
   integer         , allocatable     :: nban  (:,:)   !< base area pointers to banf, 1,* = netnode number, 2,* = flow node number, 3,* = link number, 4,* = 2nd link number
   integer                           :: mxban         !< max dim of ban
   
   type tnod
      integer, ALLOCATABLE          :: lin(:)          !< Link nrs (==index in kn(:,L)).
   end type tnod
   type (tnod), allocatable              :: nod (:)         !< (numk) Net node connectivity.
   type (tnod), allocatable              :: nod0(:)         !< Backup for nod.
   double precision, allocatable, target :: xk(:)           !< [-] Net node x coordinate {"shape": ["numk"]}
   double precision, allocatable, target :: yk(:)           !< [-] Net node y coordinate {"shape": ["numk"]}
   double precision, allocatable, target :: zk(:)           !< [-] Net node z coordinate {"shape": ["numk"]}
   double precision, allocatable         :: XK0(:), YK0(:), ZK0(:) !< Backup for xk, etc.
   double precision, allocatable         :: XK1(:), YK1(:), ZK1(:) !< Work array for xk, etc.   
   
   type tface
     integer                         :: n               !< nr of nodes
     integer, allocatable            :: nod(:)          !< node nrs
     integer, allocatable            :: lin(:)          !< link nrs, kn(1 of 2,netcell(n)%lin(1)) =  netcell(n)%nod(1)
   end type tface                    
   type (tface), allocatable         :: netcell(:)      !< (nump1d2d) 1D&2D net cells (nodes and links)
   integer,  allocatable             :: cellmask(:)     !< (nump) Mask array for net cells
   
   ! Net link related :
   integer,  allocatable, target    :: kn(:,:)         !< [-] Net links: kn(1,:)=from-idx, kn(2,:)=to-idx, kn(3,:)=net link type (0/1/2/3/4) {"shape": [3, "numl"]}
   integer,  allocatable            :: KN0(:,:)        !< Backup for kn.
   integer,  allocatable            :: LC(:)           !< (numl) Mask array for net links.
   integer,  allocatable            :: LC0(:)          !< Backup for lc.
   real   , allocatable             :: RLIN(:)         !< (numl) Placeholder for link values to be displayed.
   double precision, allocatable    :: xe(:), ye(:)    !< (numl) Edge (link) center coordinates.
   double precision, allocatable    :: dxe(:)          !< (numl) Edge (link) actual length. OPTIONAL. When unallocated, we default to Euclidean distance between the netnodes xk,yk.
   double precision, allocatable    :: dxe0(:)         !< Backup for dxe.
   
   ! Edge (and cell) related :      ! there are more edges than flow links .....
   integer, allocatable             :: lne(:,:)        !< (2,numl) Edge administration 1=nd1 , 2=nd2, rythm of kn
                                                       !! flow nodes between/next to which this net link lies.
   integer, allocatable             :: lne0(:,:)       ! backup of lne
   integer, allocatable             :: LNN(:)          !< (numl) Nr. of cells in which link participates (ubound for non-dummy values in lne(:,L))
   integer, allocatable             :: LNN0(:)
   integer                          :: NUMK0
   integer, target                  :: numk            !< [-] Nr. of net nodes. {"shape": []}
   integer                          :: NUML0, NUML     !< Total nr. of net links. In link arrays: 1D: 1:NUML1D, 2D: NUML1D+1:NUML
   integer                          :: NUML1D          !< Nr. of 1D net links.
   integer                          :: NUMP0, NUMP     !< Nr. of 2d netcells.
   integer                          :: nump1d2d        !< nr. of 1D and 2D netcells (2D netcells come first)
   integer                          :: nump1d2d0       !< nr. of 1D and 2D netcells (2D netcells come first)
   integer                          :: KN3TYP = 2      !< Default netlink type (1D/2D).
   integer                          :: jconn = 0
   
   integer,  allocatable            :: NMK (:)         !< (numk) Nr. of neighbouring netnodes for each netnode (ubound for nod(k)%lin).
   integer,  allocatable            :: KC  (:)         !< (numk) Mask array for net nodes.
   integer,  allocatable            :: NMK0(:)         !< Backup for nmk.
   integer,  allocatable            :: KC0 (:)         !< Backup for kc.
   integer,  allocatable            :: NB  (:)         !< (numk) Node codes (corner/boundary/internal classification)   

   integer                          :: LNUMK = 0, LNUML = 0   
   integer                          :: linmin = 0, linmax = 0, nodmin= 0, nodmax = 0, netcelmax = 0, netcelmin=0
   integer                          :: makeorthocenters = 0              !< shift from circumcentre to orthocentre (acts as a maxiter)
   !  network administration status
   integer, parameter               :: NETSTAT_OK    = 0 !< Network administration up-to-date
   integer, parameter               :: NETSTAT_CELLS_DIRTY = 1 !< Network administration needs findcells call.
   integer                          :: netstat = NETSTAT_CELLS_DIRTY
 
   ! keep circumcenters before orthogonalization in case of quadtree meshes
   integer                          :: keepcircumcenters = 0    !< keep circumcenter (1) or not (0)
 
   !  netlink permutation by setnodadm
   integer, dimension(:), allocatable :: Lperm    !< permuation of netlinks by setnodadm, dim(numL): from current to old link numbers
   integer, dimension(:), allocatable :: Lperminv !< Inverse permutation of netlinks by setnodadm, dim(numl): from old to current link numbers
   !  netnode permutation by setnodadm
   integer, dimension(:), allocatable :: nodePermutation   !< permutation of netnodes by setnodadm, dim(numk)     
   
   type(t_crs), target :: crs !< crs read from net file, to be written to flowgeom. TODO: AvD: temp, move this global CRS into ug_meshgeom (now a bit difficult with old and new file format)
   
   integer, allocatable             :: linkcross(:,:)    !< (2,nlinkcross) Helparray, pairs of net links that cross each other.
   integer                          :: nlinkcross        !< Nr. of crossing net links detected.
   integer, allocatable             :: linkbadqual(:)    !< (nlinkbadortho+nlinktoosmall) Net link nrs with a bad flow link orthogonality (>cosphiutrsh) or too short (<removesmalllinktrsh*...).
                                                         !! 1:nlinkbadortho for badortho links
                                                         !! nlinkbadortho+1:nlinkbadortho+nlinktoosmall for too short flow links.
   integer                          :: nlinkbadortho = 0 !< Nr. of net links with bad orthogonality detected.
   integer                          :: nlinktoosmall = 0 !< Nr. of net links with too small flow links across them.
   
   integer                          :: MMAX_old = 3, NMAX_old = 3
   integer                          :: KMAX, LMAX, KNX, MXB
   
   integer                          :: jathindams = 0                    !< For quick check whether any of the kn(3,:)==0
   
   double precision                 :: TRIANGLEMINANGLE =  5d0 ! MINIMUM ANGLE IN CREATED TRIANGLES  IF MINANGLE > MAXANGLE: NO CHECK
   
   contains
   
   subroutine sedtrails_resetdata()
   
      implicit none
      
      sedtrails_analysis = ''
      ! node (s) related : dim=ndx
      ndx2D   = 0      ! nr of 2d FLOW CELLS = NUMP
      ndxi    = 0      ! max nr of internal flowcells  (internal = 2D + 1D )
      ndx     = 0      ! nr of flow nodes (internal + boundary)
      
      ! link (u) related : dim = lnx
      lnxi    = 0      ! nr of flow links (internal           )
      lnx     = 0      ! nr of flow links (internal + boundary)
      
      MMAX_old = 3
      NMAX_old = 3
      KMAX = 0
      LMAX = 0
      KNX   = 0
      MXB  = 0

   end subroutine
   
end module   