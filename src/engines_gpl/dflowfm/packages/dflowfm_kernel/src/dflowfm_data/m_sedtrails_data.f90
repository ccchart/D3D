module m_sedtrails_data
   use coordinate_reference_system
   implicit none
   
   character(len=255)                    :: sedtrails_analysis
   !
   ! flow geometry from md_sedtrailsfile
   ! 
   double precision, allocatable, target :: xk(:)           !< [-] Net node x coordinate {"shape": ["numk"]}
   double precision, allocatable, target :: yk(:)           !< [-] Net node y coordinate {"shape": ["numk"]}
   double precision, allocatable, target :: zk(:)           !< [-] Net node z coordinate {"shape": ["numk"]}
   double precision, allocatable         :: XK0(:), YK0(:), ZK0(:) !< Backup for xk, etc.
   double precision, allocatable         :: XK1(:), YK1(:), ZK1(:) !< Work array for xk, etc.   
   

   integer                               :: NUMK0
   integer, target                       :: numk            !< [-] Nr. of net nodes. {"shape": []}
   integer                               :: KMAX
   type(t_crs), target                   :: crs             !< crs read from net file, to be written to flowgeom. TODO: AvD: temp, move this global CRS into ug_meshgeom (now a bit difficult with old and new file format)
   
   integer, allocatable                  :: st_ind(:,:)     !< indexes and weight factors for interpolation cell centre->sedtrails grid
   double precision, allocatable         :: st_wf(:,:)      !< (3,:)
   
   integer, allocatable                  :: idomain(:) 
   integer, allocatable                  :: iglobal_s(:) 
   integer, allocatable                  :: iwork(:) 
      
   contains
   
   subroutine sedtrails_resetdata()
   
      implicit none
      
      sedtrails_analysis = 'all'

   end subroutine
   
end module   