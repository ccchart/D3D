module m_sedtrails_netcdf
   use m_sedtrails_data
   use unstruc_netcdf
   
   implicit none
   
   type t_unc_sedtrailsids
      integer                  :: ncid = 0 !< NetCDF data set id (typically NetCDF file pointer)
   end type t_unc_sedtrailsids
   
   contains
   
subroutine sedtrails_write_stats(tim)
   use m_flowparameters, only: eps10
   use m_flowtimes, only: ti_st, ti_sts, ti_ste, tstop_user, time_st   
   use precision_basics
   use m_sedtrails_stats
   
   implicit none
   
   double precision, intent(in)      :: tim
   integer                           :: ierr
   
   ierr = 1
   if (ti_st > 0) then
      if (comparereal(tim, time_st, eps10) >= 0) then
         call sedtrails_write_nc(time_st)
         call reset_sedtrails_stats()
         if (ti_st > 0) then
             time_st = max(ti_sts + (floor((tim-ti_sts)/ti_st)+1)*ti_st,ti_sts)
         else
             time_st = tstop_user
         endif
         if (comparereal(time_st, ti_ste, eps10) == 1) then
             time_st = tstop_user
         endif
      endif
   end if
   
   ierr = 0
   
1234 continue
   return
end subroutine

subroutine sedtrails_write_nc(tim)
    use netcdf
    use m_flowtimes, only: it_st
    use unstruc_model
    use unstruc_files , only: defaultFilename
    use Messagehandling
    implicit none

    double precision, intent(in) :: tim

    type(t_unc_sedtrailsids), save :: stids
    integer                        :: ierr
    character(len=256)             :: filnam

    if (stids%ncid /= 0 .and. it_st==0) then
       ierr = unc_close(stids%ncid)
       stids%ncid = 0
    end if


    if (stids%ncid == 0) then
       filnam = defaultFilename('sedtrails')
       ierr = unc_create(filnam , 0, stids%ncid)
       if (ierr /= nf90_noerr) then
           call mess(LEVEL_WARN, 'sedtrails_write_nc :: Could not create sedtrails output file.')
           stids%ncid = 0
       end if
    endif

    if (stids%ncid .ne. 0) then
       call unc_write_sedtrails_filepointer(stids%ncid,tim)  
    endif

    if (unc_noforcedflush == 0) then
       ierr = nf90_sync(stids%ncid) ! Flush file
    end if

end subroutine sedtrails_write_nc   
   
!> Reads the net data from a NetCDF file.
!! Processing is done elsewhere.
subroutine sedtrails_unc_read_net_ugrid(filename, numk_keep,  numk_read, ierr)
   use m_sedtrails_data
   use m_sedtrails_network, only: sedtrails_increasenetwork
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
   integer,          intent(out)   :: numk_read          !< Number of new netnodes read from file.
   integer,          intent(out)   :: ierr               !< Return status (NetCDF operations)

   integer                         :: ioncid, iconvtype, start_index, networkIndex
   integer                         :: im, nmesh, i, numk_last
   integer                         :: ncid
   double precision                :: convversion
   type(t_ug_meshgeom)             :: meshgeom

   character(len=:), allocatable             :: tmpstring
   logical :: need_edgelengths, includeArrays

   numk_read = 0
   start_index = 1
   numk_last = numk_keep
   includeArrays = .true.
   networkIndex = 0
   
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
   if (ierr /= ionc_noerr) then
      call mess(LEVEL_WARN,  'sedtrails_unc_read_net_ugrid::ionc_get_coordinate_system: No epsg_code found in UGRID net file '''//trim(filename)//'''.')
      goto 999
   end if
   !
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
      else
         ! Only support 2D grid
         write(msgbuf, '(a,i0,a,i0,a)') 'sedtrails_unc_read_net_ugrid: unsupported topology dimension ', meshgeom%dim, &
            ' in file '''//trim(filename)//' for mesh #', im, '.'
         call warn_flush()
         cycle
      end if
 
      need_edgelengths = .false. ! Either from a previous meshgeom, or now for the first time.
      
      !increasenetw 
      call sedtrails_increasenetwork(numk_last + meshgeom%numnode) ! increases XK, YK, KN, optionally dxe
      if (meshgeom%dim == 2) then  ! always true, see above
         ! 2D
         ierr = ionc_get_node_coordinates(ioncid, im, XK(numk_last+1:numk_last + meshgeom%numnode), YK(numk_last+1:numk_last + meshgeom%numnode))
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
      
      numk_read = numk_read + meshgeom%numnode 
      numk_last = numk_last + meshgeom%numnode  

   end do
   
   call realloc(zk,numk_read,keepExisting=.false.,fill=0d0)

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
subroutine sedtrails_unc_read_net(filename, numk_keep, numk_read, ierr)
    use m_sedtrails_data
    use m_sferic
    use m_missing
    use dfm_error
    use netcdf_utils, only: ncu_get_att, ncu_get_var_attset

    character(len=*), intent(in)     :: filename  !< Name of NetCDF file.
    integer,          intent(inout)  :: numk_keep !< Number of netnodes to keep in existing net.
    integer,          intent(out)    :: numk_read !< Number of new netnodes read from file.
    integer,          intent(out)    :: ierr      !< Return status (NetCDF operations)

    call mess(LEVEL_INFO,'sedtrails_unc_read_net::Reading net data...')
    !
    ! Try and read as new UGRID NetCDF format
    !
    call sedtrails_unc_read_net_ugrid(filename, numk_keep, numk_read, ierr)
    if (ierr == dfm_noerr) then
       ! UGRID successfully read, we're done.
       return
    else
       ! No UGRID, error.
       call mess(LEVEL_ERROR,'sedtrails_unc_read_net::Could not read '''//trim(filename)//'''')
       return
    end if

end subroutine
   
   
!> Writes the unstructured sedtrails geometry to an already opened netCDF dataset.
subroutine sedtrails_unc_write_flowgeom_filepointer(igeomfile)
    use m_sedtrails_data
    use m_sferic
    use m_missing
    use netcdf
    use m_partitioninfo, only: jampi
    use m_flowparameters, only: jafullgridoutput
    
    implicit none
    
    integer, intent(in)             :: igeomfile
    integer                         :: ndxndxi   
    integer :: ierr
    integer ::  id_flowelemdim,& 
                id_flowelemxcc, id_flowelemycc, &
                id_flowelemdomain, id_flowelemglobalnr

    integer :: jaInDefine

    jaInDefine = 0

    if (numk <= 0) then
        call mess(LEVEL_WARN, 'sedtrails_unc_write_flowgeom_filepointer :: No flow elements in model, will not write flow geometry.')
        return
    end if

    ndxndxi=numk

    ! Put dataset in define mode (possibly again) to add dimensions and variables.
    ierr = nf90_redef(igeomfile)
    if (ierr == nf90_eindefine) jaInDefine = 1 ! Was still in define mode.
    if (ierr /= nf90_noerr .and. ierr /= nf90_eindefine) then
        call mess(LEVEL_ERROR, 'sedtrails_unc_write_flowgeom_filepointer::Could not put header in sedtrails geometry file.')
        call check_error(ierr)
        return
    end if

    ierr = nf90_def_dim(igeomfile, 'nNodes',  ndxndxi, id_flowelemdim)

    ! Net nodes
    ierr = nf90_def_var(igeomfile, 'net_xcc', nf90_double, id_flowelemdim, id_flowelemxcc)
    ierr = nf90_def_var(igeomfile, 'net_ycc', nf90_double, id_flowelemdim, id_flowelemycc)
    
    ierr = unc_addcoordatts(igeomfile, id_flowelemxcc, id_flowelemycc, jsferic)
    ierr = nf90_put_att(igeomfile, id_flowelemxcc, 'long_name'    , 'x-coordinate of sedtrails grid corner')
    ierr = nf90_put_att(igeomfile, id_flowelemycc, 'long_name'    , 'y-coordinate of sedtrails grid corner')

    ! Coordinate/grid mapping
    ierr = unc_addcoordmapping(igeomfile, jsferic)
    
    !   domain numbers
    if ( jampi.eq.1 ) then
       ierr = nf90_def_var(igeomfile, 'FlowElemDomain', nf90_int, id_flowelemdim, id_flowelemdomain)
       ierr = nf90_put_att(igeomfile, id_flowelemdomain, 'long_name',   'domain number of sedtrails grid corner')
       ierr = nf90_def_var(igeomfile, 'FlowElemGlobalNr', nf90_int, id_flowelemdim, id_flowelemglobalnr)
       ierr = nf90_put_att(igeomfile, id_flowelemglobalnr, 'long_name',   'global sedtrails corner numbering')
    end if

    ierr = nf90_enddef(igeomfile)
    ! End of writing time-independent flow net data.

    ! -- Start data writing (time-independent data) ------------
    ! Flow cell cc coordinates
    ierr = nf90_put_var(igeomfile, id_flowelemxcc, xk(1:ndxndxi))
    ierr = nf90_put_var(igeomfile, id_flowelemycc, yk(1:ndxndxi))
    
    ! domain numbers
    if ( jampi.eq.1 ) then  
       ! flow cell domain number   
       ierr = nf90_put_var(igeomfile, id_flowelemdomain, idomain(1:ndxndxi) )  
       ierr = nf90_put_var(igeomfile, id_flowelemglobalnr, iglobal_s(1:ndxndxi))
    end if

    ! Leave the dataset in the same mode as we got it.
    if (jaInDefine == 1) then
        ierr = nf90_redef(igeomfile)
    end if

end subroutine sedtrails_unc_write_flowgeom_filepointer


subroutine unc_write_sedtrails_filepointer(imapfile,tim)
   use m_sedtrails_stats
   use m_alloc
   use m_flowtimes, only: Tudunitstr
   use m_fm_erosed
   use m_sediment, only: stm_included
   use m_missing, only: dmiss
   use m_flowtimes, only: it_st
   use m_flow, only: hs, ucx, ucy
   use m_flowgeom, only: bl
   
   implicit none
   
   integer, intent(in)                    :: imapfile
   double precision, intent(in)           :: tim
   
   ! Locals
   integer                                :: ndxndxi
   integer                                :: iid, k, l, ii
   integer                                :: ierr, itim
   integer, save                          :: ndim
   integer, dimension(2)                  :: idims
   integer, dimension(2), save            :: id_timedim, id_time, id_timestep, id_sbx, id_sby, id_ssx, id_ssy, id_ssc, &
                                             id_sedtotdim,id_flowelemdim, id_ucx, id_ucy, id_bl, id_hs, id_taus, id_tausmax
   double precision, allocatable          :: work(:,:)
   integer         , allocatable          :: nodes(:)
   
   ! Define variables and write time-invariant data
   if (numk <= 0) then
      call mess(LEVEL_WARN, 'unc_write_sedtrails_filepointer :: No flow elements in sedtrails grid, will not write geometry.')
      return
   end if

   iid = 1
   ndxndxi = numk
    
   ! Use nr of dimensions in netCDF file a quick check whether vardefs were written
   ! before in previous calls.
   ndim = 0
   ierr = nf90_inquire(imapfile, nDimensions=ndim)

   ! Only write net and flow geometry data the first time, or for a separate map file.
   if (ndim == 0) then
       call sedtrails_unc_write_flowgeom_filepointer(imapfile)      ! Write standard net data as well
       ierr = nf90_inq_dimid(imapfile, 'nNodes', id_flowelemdim(iid))

       ! Time
       ierr = nf90_def_dim(imapfile, 'time', nf90_unlimited, id_timedim(iid))
       call check_error(ierr, 'def time dim')
       ierr = nf90_def_var(imapfile, 'time', nf90_double, id_timedim(iid),  id_time(iid))
       ierr = nf90_put_att(imapfile, id_time(iid),  'units'        , trim(Tudunitstr))
       ierr = nf90_put_att(imapfile, id_time(iid),  'standard_name', 'time')
             
       ! Size of latest timestep
       ierr = nf90_def_var(imapfile, 'timestep', nf90_double, id_timedim(iid),  id_timestep(iid))
       ierr = nf90_put_att(imapfile, id_timestep(iid),  'units'        , 'seconds')
       ierr = nf90_put_att(imapfile, id_timestep(iid),  'standard_name', 'timestep')
      
       idims(1) = id_flowelemdim(iid)
       idims(2) = id_timedim(iid)

       ! Variables that will always be written
       call definencvar(imapfile,id_bl(iid)   ,nf90_double,idims,2, 'bedlevel', 'bed level', 'm', 'net_xcc net_ycc')
       
       if (trim(sedtrails_analysis)=='flowvelocity' .or. trim(sedtrails_analysis)=='all') then
          call definencvar(imapfile,id_ucx(iid),nf90_double,idims,2, 'sea_water_x_velocity', 'depth-averaged velocity on grid corner, x-component', 'm s-1', 'net_xcc net_ycc')
          call definencvar(imapfile,id_ucy(iid),nf90_double,idims,2, 'sea_water_y_velocity', 'depth-averaged velocity on grid corner, y-component', 'm s-1', 'net_xcc net_ycc')
       endif

       if ((trim(sedtrails_analysis)=='transport' .or. trim(sedtrails_analysis)=='all') .and. stm_included) then
          ! Fractions
          ierr = nf90_def_dim(imapfile, 'nSedTot', lsedtot, id_sedtotdim(iid))
          !
          call definencvar(imapfile,id_sbx(iid),nf90_double,(/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /),3, 'bedload_x_comp', 'bedload transport on grid corner, x-component', 'kg m-1 s-1', 'net_xcc net_ycc')
          call definencvar(imapfile,id_sby(iid),nf90_double,(/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /),3, 'bedload_y_comp', 'bedload transport on grid corner, y-component', 'kg m-1 s-1', 'net_xcc net_ycc')
          call definencvar(imapfile,id_ssx(iid),nf90_double,(/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /),3, 'susload_x_comp', 'suspended transport on grid corner, x-component', 'kg m-1 s-1', 'net_xcc net_ycc')
          call definencvar(imapfile,id_ssy(iid),nf90_double,(/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /),3, 'susload_y_comp', 'suspended transport on grid corner, y-component', 'kg m-1 s-1', 'net_xcc net_ycc')
       endif    
       
       if ((trim(sedtrails_analysis)=='soulsby' .or. trim(sedtrails_analysis)=='all') .and. stm_included) then
          ! Fractions
          ierr = nf90_inq_dimid(imapfile, 'nSedTot', id_sedtotdim(iid))
          if (ierr>0) then
             ierr = nf90_def_dim(imapfile, 'nSedTot', lsedtot, id_sedtotdim(iid)) ! only once
          endif
          !
          call definencvar(imapfile,id_hs(iid)     ,nf90_double,idims,2, 'waterdepth', 'water depth', 'm', 'net_xcc net_ycc')
          call definencvar(imapfile,id_taus(iid)   ,nf90_double,idims,2, 'mean_bss_magnitude', 'mean bed shear stress magnitude', 'Pa', 'net_xcc net_ycc')
          call definencvar(imapfile,id_tausmax(iid),nf90_double,idims,2, 'max_bss_magnitude' , 'max bed shear stress magnitude', 'Pa', 'net_xcc net_ycc')
          call definencvar(imapfile,id_ssc(iid)    ,nf90_double,(/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /),3, 'suspended_sed_conc', 'depth-averaged suspended sediment concentration', 'kg m-3', 'net_xcc net_ycc')
       endif   
       
       ierr = nf90_enddef(imapfile)
   endif
   
   ! Interpolate and write variables
   it_st   = it_st+1
   itim    = it_st ! Increment time dimension index  
    
   ! Time
   ierr = nf90_put_var(imapfile, id_time(iid), tim, (/ itim /))
   ierr = nf90_put_var(imapfile, id_timestep(iid), is_dtint, (/ itim /))
   
   ! Analysis:
   call realloc(work,(/ numk, lsedtot/), keepExisting=.false., fill=0d0)
   
   ! Bottom level
   do k=1, numk
      nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
      work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_BL,st_ind(nodes,k), 1))
   enddo 
   work=work/is_dtint
   ierr = nf90_put_var(imapfile, id_bl(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi, 1 /))
      
   ! 'FLOWVELOCITY'
   if ((trim(sedtrails_analysis)=='flowvelocity' .or. trim(sedtrails_analysis)=='all')) then
      do k=1, numk
         nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
         work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_UCX,st_ind(nodes,k), 1))
      enddo 
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_ucx(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi,1 /))
      !
      do k=1, numk
         nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
         work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_UCY,st_ind(nodes,k), 1))
      enddo 
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_ucy(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi, 1 /))  
   endif
 
   !'TRANSPORT'
   if ((trim(sedtrails_analysis)=='transport' .or. trim(sedtrails_analysis)=='all') .and. stm_included) then
      do l=1,lsedtot
         do k=1, numk
            nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
            work(k,l) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_SBX,st_ind(nodes,k), l))
         enddo 
      enddo
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_sbx(iid)  , work(1:numk,1:lsedtot), (/ 1, 1, itim /), (/ ndxndxi, lsedtot, 1 /))
      !
      do l=1,lsedtot
         do k=1, numk
            nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
            work(k,l) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_SBY,st_ind(nodes,k), l))
         enddo 
      enddo
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_sby(iid)  , work(1:numk,1:lsedtot), (/ 1, 1, itim /), (/ ndxndxi, lsedtot, 1 /))   
      !
      do l=1,lsedtot
         do k=1, numk
            nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
            work(k,l) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_SSX,st_ind(nodes,k), l))
         enddo 
      enddo
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_ssx(iid)  , work(1:numk,1:lsedtot), (/ 1, 1, itim /), (/ ndxndxi, lsedtot, 1 /))  
      !
      do l=1,lsedtot
         do k=1, numk
            nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
            work(k,l) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_SSY,st_ind(nodes,k), l))
         enddo 
      enddo
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_ssy(iid)  , work(1:numk,1:lsedtot), (/ 1, 1, itim /), (/ ndxndxi, lsedtot, 1 /))      
   endif
   
   !"SOULSBY"
   if ((trim(sedtrails_analysis)=='soulsby' .or. trim(sedtrails_analysis)=='all') .and. stm_included) then
      do k=1, numk
         nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
         work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_HS,st_ind(nodes,k), 1))
      enddo 
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_hs(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi, 1 /))
      !
      do k=1, numk
         nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
         work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_TAUS,st_ind(nodes,k), 1))         
      enddo 
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_taus(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi, 1 /))
      !
      do k=1, numk
         nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
         work(k,1) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_TAUSMAX,st_ind(nodes,k), 1))         
      enddo 
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_tausmax(iid)  , work(1:numk,1), (/ 1, itim /), (/ ndxndxi, 1 /))  
      !
      do l=1,lsedtot
         do k=1, numk
            nodes = PACK([(ii,ii=1,SIZE(st_ind(:,k)))], st_ind(:,k) > 0)
            work(k,l) = sum(st_wf(nodes,k)*is_sumvalsnd(IDX_SSC,st_ind(nodes,k), l))            
         enddo 
      enddo
      work=work/is_dtint
      ierr = nf90_put_var(imapfile, id_ssc(iid), work(1:numk,1:lsedtot), (/ 1, 1, itim /), (/ ndxndxi, lsedtot, 1 /))       
   endif   
   
end subroutine

subroutine sedtrails_loadNetwork(filename, istat, jadoorladen)
    use unstruc_messages
    use m_missing
    use m_alloc
    use m_partitioninfo, only: jampi
    use m_sedtrails_data

    implicit none


    character(*), intent(in)  :: filename !< Name of file to be read (in current directory or with full path).
    integer,      intent(out) :: istat    !< Return status (0=success)
    integer,      intent(in)  :: jadoorladen
   
    integer                   :: minp,  K0,  NUMKN
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
    call sedtrails_unc_read_net(filename, K0, NUMKN, istat)
   
    if (istat == 0) then
        NUMK = K0 + NUMKN
        
        ! Fill global domain mode numbers, before reducing the sedtrails grid to the submodel we're in
        if (jampi>0) then
           call realloc(iglobal_s,numk,keepExisting=.false., fill=0)
           iglobal_s=(/ 1:numk /)
        endif   
    else
       call qnerror('sedtrails_loadNetwork::Error while loading network from '''//trim(filename)//''', please inspect the preceding diagnostic output.', ' ',  ' ')
    endif
end subroutine sedtrails_loadNetwork  
 
end module