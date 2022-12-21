!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
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

!> Write history data in NetCDF format.
subroutine unc_write_his(hisids, tim, ihisdim)            ! wrihis
    use Timers
    use m_flowtimes
    use m_flow
    use m_flowgeom
    use network_data, only: xk, yk
    use m_observations
    use m_monitoring_crosssections
    use m_monitoring_runupgauges
    use m_missing
    use netcdf
    use netcdf_utils
    !use coordinate_reference_system, only: transform_and_put_latlon_coordinates
    use unstruc_files, only: defaultFilename
    use unstruc_netcdf, only: unc_create, unc_close, unc_addcoordatts, unc_addcoordmapping, unc_def_var_nonspatial, definencvar, unc_meta_add_user_defined
    use unstruc_netcdf, only: UNC_LOC_S3D, UNC_LOC_WU, UNC_LOC_W
    use unstruc_netcdf, only: unc_writeopts, unc_noforcedflush, UG_WRITE_LATLON
    use unstruc_netcdf, only: unc_add_time_coverage
    use unstruc_netcdf, only: unc_check_dimension, UNC_DIM_1D, UNC_DIM_2D, UNC_DIM_3D, UNC_DIM_ALL
    use unstruc_netcdf, only: t_unc_hisids, unc_put_3d_or_depthaveraged_valobs
    use unstruc_messages
    use m_sferic, only: jsferic
    use m_partitioninfo
    use m_timer
    use m_sediment
    use m_flowexternalforcings, only: numtracers, trnames
    use m_transport, only: NUMCONST_MDU, ITRA1, ITRAN, ISED1, ISEDN, const_names, const_units, NUMCONST, itemp, isalt
    use m_structures
    use m_particles, only: japart
    use m_fm_wq_processes
    use string_module
    use m_dad
    use m_filter, only: checkmonitor
    use m_alloc
    use unstruc_channel_flow, only: network
    use simple_geometry, only: sgeom_def_geometry_variables
    use m_1d_structures
    use m_structures
    use m_GlobalParameters
    use m_longculverts

    implicit none

    type(t_unc_hisids), intent(inout) :: hisids  !< Set of file and variable ids for this map-type file.
    double precision,   intent(in)    :: tim     !< Current time, should in fact be time1, since the data written is always s1, ucx, etc.
    integer,            intent(in)    :: ihisdim !< Dimensional type of his file to optionally select only 1D/2D or 3D output: sum of constants UNC_DIM_1D/2D/3D. Use UNC_DIM_ALL for regular his file.


    double precision, allocatable :: geom_x(:), geom_y(:)
    integer, allocatable          :: node_count(:)
    double precision, allocatable, save :: valobsT(:,:)
    integer :: maxlocT, maxvalT !< row+column count of valobsT

    integer                      :: IP, num, ntmp, n, nlyrs

    double precision, save       :: curtime_split = 0d0 ! Current time-partition that the file writer has open.
    integer                      :: ntot, k, i, j, jj, ierr, mnp, kk, idims(3),L, Lf, k3, k4, nNodeTot, nNodes, L0, k1, k2, nlinks
    character(len=255)           :: station_geom_container_name, crs_geom_container_name, weir_geom_container_name, orif_geom_container_name, &
                                    genstru_geom_container_name, uniweir_geom_container_name, culvert_geom_container_name, longculvert_geom_container_name, &
                                    gategen_geom_container_name, pump_geom_container_name, bridge_geom_container_name, src_geom_container_name, &
                                    lat_geom_container_name
    double precision             :: cof0

    integer                      :: strlen_netcdf  ! string length definition for (station) names on history file
    character(len=255)           :: filename
    character(len=25)            :: transpunit
    character(len=1024)          :: statcoordstring
    integer                      :: igen, istru
    integer                      :: ndims
    character(len=255)           :: tmpstr
    integer                      :: jawrizc = 0
    integer                      :: jawrizw = 0
    double precision             :: w1, pumplensum, pumplenmid, pumpxmid, pumpymid
    double precision             :: rhol
    double precision, allocatable:: toutput1(:), toutputx(:,:), toutputy(:,:), toutput3(:,:,:)
    double precision, allocatable:: toutput_cum, toutput_cur
    integer                      :: lsed
    logical                      :: add_latlon

    integer                      :: janormalhis   !< Normal his, also including space-independent quantities
    integer                      :: ja3d          !< Include 3D quantities in this his file

    ! Store 3D depth-varying variables in file?
    if (unc_check_dimension(ihisdim, UNC_DIM_3D)) then
       ja3d = 1
    else
       ja3d = 0
    end if

    ! Store 1D/2D and depth-averaged variables in file?
    if (unc_check_dimension(ihisdim, UNC_DIM_1D) .or. unc_check_dimension(ihisdim, UNC_DIM_2D)) then
       janormalhis = 1
    else
       janormalhis = 0
    end if

    if (jahiszcor > 0) then
       jawrizc = 1
       jawrizw = 1
    endif

    if (timon) call timstrt ( "unc_write_his", handle_extra(54))

    ! Another time-partitioned file needs to start, reset iteration count (and file).
    if (ti_split > 0d0 .and. curtime_split /= time_split0) then
        hisids%idx_curtime        = 0
        curtime_split = time_split0
    end if

    ! Close/reset any previous hisfile.
    if (hisids%ncid/=0) then  ! reset stord ncid to zero if file not open
       ierr = nf90_inquire(hisids%ncid, ndims)
       if (ierr/=0) hisids%ncid = 0
    end if

    if (hisids%ncid > 0 .and. hisids%idx_curtime == 0) then
        ierr = unc_close(hisids%ncid)
        hisids%ncid = 0
    end if

    ! When no crs/obs present, return immediately.
    if (numobs+nummovobs <= 0 .and. ncrs <= 0 .and. jahisbal <= 0 .and. jahiscgen <= 0 .and. nrug <= 0) then
        if (hisids%ncid == 0) then
            call mess(LEVEL_WARN, 'No observations nor cross sections defined. Will not produce a history file.')
        end if
        hisids%ncid = -1 ! -1 stands for: no file open, no obs/crs defined.
        return
    end if

    ! Only add auto-tranformed lat/lon coordinates if model is Cartesian and user has requested extra latlon output.
#ifdef HAVE_PROJ
    add_latlon = jsferic == 0 .and. iand(unc_writeopts, UG_WRITE_LATLON) == UG_WRITE_LATLON
#else
    add_latlon = .false.
#endif

    if (hisids%ncid == 0) then
        if (timon) call timstrt ( "unc_write_his INIT/DEF", handle_extra(61))

        call realloc(hisids%id_tra, ITRAN-ITRA1+1, keepExisting = .false.)
        call realloc(hisids%id_const, NUMCONST_MDU, keepExisting = .false.)
        call realloc(hisids%id_const_cum, NUMCONST_MDU, keepExisting = .false.)

        call realloc(hisids%id_voltot, MAX_IDX, keepExisting = .false.)

        ! Possibly a different model, so make valobs transpose at correct size again.
        maxlocT = max(size(valobs, 2), npumpsg, network%sts%numPumps, ngatesg, ncdamsg, ncgensg, ngategen, &
                      nweirgen, network%sts%numWeirs, ngenstru,  network%sts%numGeneralStructures, &
                      ndambreak, network%sts%numOrifices, network%sts%numBridges, network%sts%numculverts, &
                      network%sts%numuniweirs, network%cmps%count, nlongculvertsg)
        maxvalT = max(size(valobs, 1), NUMVALS_PUMP, NUMVALS_GATE, NUMVALS_CDAM, NUMVALS_CGEN, NUMVALS_GATEGEN, &
                      NUMVALS_WEIRGEN, NUMVALS_GENSTRU, &
                      NUMVALS_DAMBREAK, NUMVALS_ORIFGEN, NUMVALS_BRIDGE, NUMVALS_CULVERT, &
                      NUMVALS_UNIWEIR, NUMVALS_CMPSTRU, NUMVALS_LONGCULVERT)
        call realloc(valobsT, (/ maxlocT, maxvalT /), keepExisting = .false.)

        if (janormalhis > 0) then
           if (ti_split > 0d0) then
               filename = defaultFilename('his', timestamp=time_split0)
           else
               filename = defaultFilename('his')
           end if
        elseif (kmx > 0 .and. ja3d > 0) then
            filename = defaultFilename('his3d')
        end if

        ierr = unc_create(filename, 0, hisids%ncid, .false.)
        if (ierr /= nf90_noerr) then
            call mess(LEVEL_WARN, 'Could not create history file.')
        end if

        ierr = unc_meta_add_user_defined(hisids%ncid)

        ierr = unc_add_time_coverage(hisids%ncid, ti_hiss, ti_hise, ti_his)

        !if (unc_nounlimited > 0) then ! UNST-4764: His file has never shown good results with NcNoUnlimited option on.
        ierr = nf90_def_dim(hisids%ncid, 'time', nf90_unlimited, hisids%id_timedim)

        strlen_netcdf = idlen  !< Max string length of Ids.
        ierr = nf90_def_dim(hisids%ncid, 'name_len', strlen_netcdf, hisids%id_strlendim)

        if (ja3d > 0 .and. kmx > 0) then
           ierr = nf90_def_dim(hisids%ncid, 'laydim', kmx, hisids%id_laydim)
           ierr = nf90_def_dim(hisids%ncid, 'laydimw', kmx+1, hisids%id_laydimw)
        end if

        if (numobs+nummovobs > 0) then
            ierr = unc_addcoordmapping(hisids%ncid, jsferic)
            
            ierr = nf90_def_dim(hisids%ncid, 'stations', numobs+nummovobs, hisids%id_statdim)

            if (nummovobs > 0) then
               ierr = nf90_def_var(hisids%ncid, 'station_x_coordinate', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_statx) ! TODO: AvD: decide on UNST-1606 (trajectory_id vs. timeseries_id)
               ierr = nf90_def_var(hisids%ncid, 'station_y_coordinate', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_staty)
            else
               ierr = nf90_def_var(hisids%ncid, 'station_x_coordinate', nf90_double, hisids%id_statdim, hisids%id_statx)
               ierr = nf90_def_var(hisids%ncid, 'station_y_coordinate', nf90_double, hisids%id_statdim, hisids%id_staty)
            endif
            ierr = unc_addcoordatts(hisids%ncid, hisids%id_statx, hisids%id_staty, jsferic)
            ierr = nf90_put_att(hisids%ncid, hisids%id_statx, 'long_name', 'original x-coordinate of station (non-snapped)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_staty, 'long_name', 'original y-coordinate of station (non-snapped)')

            statcoordstring = 'station_x_coordinate station_y_coordinate station_name'
            if (add_latlon) then
               ierr = ncu_clone_vardef(hisids%ncid, hisids%ncid, hisids%id_statx, 'station_lon', hisids%id_statlon, &
                             'longitude', 'original lon-coordinate of station (non-snapped)', 'degrees_east')
               ierr = ncu_clone_vardef(hisids%ncid, hisids%ncid, hisids%id_staty, 'station_lat', hisids%id_statlat, &
                             'latitude', 'original lat-coordinate of station (non-snapped)', 'degrees_north')

               statcoordstring = trim(statcoordstring) // ' station_lon station_lat'
            end if

            ierr = nf90_def_var(hisids%ncid, 'station_id',         nf90_char,   (/ hisids%id_strlendim, hisids%id_statdim /), hisids%id_statid)
            ierr = nf90_put_att(hisids%ncid, hisids%id_statid,  'long_name'    , 'observation station identifier') ! REF

            ierr = nf90_def_var(hisids%ncid, 'station_name',         nf90_char,   (/ hisids%id_strlendim, hisids%id_statdim /), hisids%id_statname)
            ierr = nf90_put_att(hisids%ncid, hisids%id_statname,  'cf_role', 'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_statname,  'long_name'    , 'observation station name') ! REF

            ! Define geometry related variables
            station_geom_container_name = 'station_geom'
            nNodeTot = numobs+nummovobs
            ierr = sgeom_def_geometry_variables(hisids%ncid, station_geom_container_name, 'station', 'point', nNodeTot, hisids%id_statdim, &
               hisids%id_statgeom_node_count, hisids%id_statgeom_node_coordx, hisids%id_statgeom_node_coordy, add_latlon, hisids%id_statgeom_node_lon, hisids%id_statgeom_node_lat)


            if ( jahiswatlev > 0 .and. janormalhis > 0) then
               ierr = nf90_def_var(hisids%ncid, 'waterlevel', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_vars)
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, 'standard_name', 'sea_surface_height') ! sorry for inland water people
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, 'long_name', 'water level')
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_vars, '_FillValue', dmiss)
            endif

            if ( jahisbedlev > 0 .and. janormalhis > 0) then
               if( stm_included ) then
                  ierr = nf90_def_var(hisids%ncid, 'bedlevel', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varb)
               else
                  ierr = nf90_def_var(hisids%ncid, 'bedlevel', nf90_double, (/ hisids%id_statdim /), hisids%id_varb)
               endif
               ierr = nf90_put_att(hisids%ncid, hisids%id_varb, 'long_name', 'bottom level')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varb, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varb, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varb, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varb, '_FillValue', dmiss)
            endif

            idims(1) = hisids%id_statdim
            idims(2) = hisids%id_timedim
            if (jahiswatdep > 0 .and. janormalhis > 0) then
               call definencvar(hisids%ncid,hisids%id_hs, nf90_double, idims, 2, 'waterdepth'  , 'water depth', 'm', statcoordstring, station_geom_container_name)
            endif

            if( jahisvelvec > 0 ) then
               if ( kmx.gt.0) then
                  if (ja3d > 0) then
                     ierr = nf90_def_var(hisids%ncid, 'x_velocity', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varucx)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     ierr = nf90_def_var(hisids%ncid, 'y_velocity', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varucy)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     ierr = nf90_def_var(hisids%ncid, 'z_velocity', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varucz)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, 'standard_name', 'upward_sea_water_velocity')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, 'long_name', 'vertical/upward component of flow element center velocity vector') ! sorry for inland water people
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, 'units', 'm s-1')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, 'geometry', station_geom_container_name)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucz, '_FillValue', dmiss)
                     jawrizc = 1
                  end if
                  if (janormalhis > 0) then ! Depth avgd velocity always in full/normal his file, or in the non-3D his file
                     ierr = nf90_def_var(hisids%ncid, 'depth-averaged_x_velocity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varucxq)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, 'standard_name', 'sea_water_depth-averaged_x_velocity')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, 'long_name', 'flow element center depth-averaged velocity vector, x-component')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, 'units', 'm s-1')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, 'geometry', station_geom_container_name)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucxq, '_FillValue', dmiss)

                     ierr = nf90_def_var(hisids%ncid, 'depth-averaged_y_velocity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varucyq)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, 'standard_name', 'sea_water_depth-averaged_y_velocity')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, 'long_name', 'flow element center depth-averaged velocity vector, y-component')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, 'units', 'm s-1')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, 'geometry', station_geom_container_name)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varucyq, '_FillValue', dmiss)
                  end if
               else
                  ierr = nf90_def_var(hisids%ncid, 'x_velocity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varucx)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'coordinates', statcoordstring)
                  ierr = nf90_def_var(hisids%ncid, 'y_velocity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varucy)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'coordinates', statcoordstring)
               end if
               if ((kmx==0 .and. janormalhis > 0) .or. (kmx > 0 .and. ja3d > 0)) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'standard_name', 'sea_water_x_velocity')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'long_name', 'x-component of flow element center velocity vector')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'units', 'm s-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'standard_name', 'sea_water_y_velocity')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'long_name', 'y-component of flow element center velocity vector')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'units', 'm s-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, '_FillValue', dmiss)
               end if
            endif

            if ( kmx.gt.0 .and. ja3d > 0) then

               !idims(1) = hisids%id_laydim
               !idims(2) = hisids%id_statdim
               !idims(3) = hisids%id_timedim
               !call definencvar   (hisids%ncid, hisids%id_zcs, nf90_double, idims,3, 'zcoordinate_c' , 'vertical coordinate at center of flow element and layer'   , 'm',  'station_x_coordinate station_y_coordinate station_name zcoordinate_c')
               !ierr = nf90_put_att(hisids%ncid, hisids%id_zcs, 'positive' , 'up')

               idims(1) = hisids%id_laydimw
               idims(2) = hisids%id_statdim
               idims(3) = hisids%id_timedim
               !call definencvar   (hisids%ncid, hisids%id_zws, nf90_double, idims,3, 'zcoordinate_w' , 'vertical coordinate at center of flow element and at layer interface'   , 'm',  'station_x_coordinate station_y_coordinate station_name zcoordinate_w')
               !ierr = nf90_put_att(hisids%ncid, hisids%id_zws, 'positive' , 'up')

               if (iturbulencemodel >= 3 .and. jahistur > 0) then
                  call definencvar(hisids%ncid,hisids%id_turkin,nf90_double, idims,3, 'tke'   , 'turbulent kinetic energy'   , 'm2 s-2',  trim(statcoordstring) // ' zcoordinate_wu', station_geom_container_name, fillVal = dmiss)
                  jawrizw = 1
               endif
               if (iturbulencemodel > 1 .and. jahistur > 0 ) then
                  call definencvar(hisids%ncid,hisids%id_vicwwu,nf90_double, idims,3, 'vicww' , 'turbulent vertical eddy viscosity'    , 'm2 s-1' ,  trim(statcoordstring) // ' zcoordinate_wu', station_geom_container_name, fillVal = dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_turkin, 'standard_name', 'specific_turbulent_kinetic_energy_of_sea_water')
                  jawrizw = 1
               endif
               if (iturbulencemodel == 3 .and. jahistur > 0) then
                  call definencvar(hisids%ncid,hisids%id_tureps,nf90_double, idims,3, 'eps'   , 'turbulent energy dissipation', 'm2 s-3'  ,  trim(statcoordstring) // ' zcoordinate_wu', station_geom_container_name, fillVal = dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_tureps, 'standard_name', 'specific_turbulent_kinetic_energy_dissipation_in_sea_water')
                  jawrizw = 1
               else if (iturbulencemodel == 4 .and. jahistur > 0) then
                  call definencvar(hisids%ncid,hisids%id_tureps,nf90_double, idims,3, 'tau'   , 'turbulent time scale', 's-1'  ,  trim(statcoordstring) // ' zcoordinate_wu', station_geom_container_name, fillVal = dmiss)
                  jawrizw = 1
               endif

               if (jarichardsononoutput > 0) then
                  call definencvar(hisids%ncid,hisids%id_rich,nf90_double, idims,3, 'rich' , 'Richardson Nr'    , '  ' ,  trim(statcoordstring) // ' zcoordinate_wu', station_geom_container_name, fillVal = dmiss)
                  jawrizw = 1
               end if

            end if
            if ((kmx==0 .and. janormalhis > 0) .or. (kmx > 0 .and. ja3d > 0)) then
               if (jsferic == 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'standard_name', 'sea_water_x_velocity')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'standard_name', 'sea_water_y_velocity')
               else
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'standard_name', 'eastward_sea_water_velocity')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'standard_name', 'northward_sea_water_velocity')
               end if

               if (jaeulervel==0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'long_name', 'flow element center velocity vector, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'long_name', 'flow element center velocity vector, y-component')
               else
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'long_name', 'flow element center Eulerian velocity vector, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'long_name', 'flow element center Eulerian velocity vector, y-component')
               endif
               ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'units', 'm s-1')

               if (kmx > 0 .and. ja3d > 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
               else
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'coordinates', statcoordstring)
               end if

               ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varucx, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varucy, '_FillValue', dmiss)
            end if

            if (jasal > 0 .and. jahissal > 0) then
               if (kmx > 0) then
                  if (ja3d > 0) then
                     ierr = nf90_def_var(hisids%ncid, 'salinity', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varsal)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     jawrizc = 1
                  else
                     ierr = nf90_def_var(hisids%ncid, 'depth-averaged_salinity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varsal)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'coordinates', trim(statcoordstring))
                  endif
               else
                  ierr = nf90_def_var(hisids%ncid, 'salinity', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varsal)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'coordinates', statcoordstring)
               end if
               ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'units', '1e-3')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varsal, 'standard_name', 'sea_water_salinity')
            endif

            if ( jahisvelocity > 0 ) then
               if (kmx > 0) then
                  if (ja3d > 0) then
                     ierr = nf90_def_var(hisids%ncid, 'velocity_magnitude', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varumag)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                  else
                     ierr = nf90_def_var(hisids%ncid, 'depth-averaged_velocity_magnitude', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varumag)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'coordinates', trim(statcoordstring))
                  end if
               else 
                  ierr = nf90_def_var(hisids%ncid, 'velocity_magnitude', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varumag)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'coordinates', statcoordstring)
               endif
               if (jaeulervel==0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'long_name', 'velocity magnitude')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'standard_name', 'sea_water_speed')
               else
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'standard_name', 'sea_water_eulerian_speed')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'long_name', 'Eulerian velocity magnitude')
               end if
               ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varumag, '_FillValue', dmiss)
            endif

            if ( jahisdischarge > 0 ) then
               if (kmx > 0) then
                  if (ja3d > 0) then
                     ierr = nf90_def_var(hisids%ncid, 'discharge_magnitude', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varqmag)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'coordinates', trim(statcoordstring))
                  end if
               else
                  ierr = nf90_def_var(hisids%ncid, 'discharge_magnitude', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varqmag)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
               endif
               
               if (hisids%id_varqmag > 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'standard_name', 'water_volume_transport_in_river_channel') 
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'long_name', 'average discharge magnitude')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'units', 'm3 s-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varqmag, '_FillValue', dmiss)
               end if
            endif


!           JRE surfbeat
            if (jawave .eq. 4 .and. janormalhis > 0) then

               ierr = nf90_def_var(hisids%ncid, 'R',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_R)
               ierr = nf90_put_att(hisids%ncid, hisids%id_R,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_R,   'standard_name', 'sea_surface_bulk_roller_energy')                          ! not CF
               ierr = nf90_put_att(hisids%ncid, hisids%id_R,   'long_name'    , 'roller energy per square meter')
               ierr = nf90_put_att(hisids%ncid, hisids%id_R,   'units'        , 'J m-2')
               ierr = nf90_put_att(hisids%ncid, hisids%id_R, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_R, '_FillValue', dmiss)

            end if

            if (jawave > 0 .and. jahiswav > 0 .and. janormalhis > 0) then
               ierr = nf90_def_var(hisids%ncid, 'hwav',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_WH)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'standard_name', 'sea_surface_wave_significant_wave_height')     ! Default behaviour
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'long_name'    , 'Significant wave height')
               if (jahissigwav==0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'standard_name', 'sea_surface_wave_rms_height')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'long_name'    , 'Root mean square wave height based on wave energy')
               end  if
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH,   'units'        , 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WH, '_FillValue', dmiss)

               ierr = nf90_def_var(hisids%ncid, 'twav',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_WT)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT,   'standard_name', 'sea_surface_wave_period')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT,   'long_name'    , 'Wave period')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT,   'units'        , 's')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WT, '_FillValue', dmiss)

               ierr = nf90_def_var(hisids%ncid, 'phiwav',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_WD)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD,   'standard_name', 'sea_surface_wave_from_direction')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD,   'long_name'    , 'Wave from direction')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD,   'units'        , 'deg from N')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WD, '_FillValue', dmiss)

               ierr = nf90_def_var(hisids%ncid, 'rlabda',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_WL)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL,   'standard_name', 'sea_surface_wave_length')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL,   'long_name'    , 'Wave length')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL,   'units'        , 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WL, '_FillValue', dmiss)

               ierr = nf90_def_var(hisids%ncid, 'uorb',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_WU)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU,   'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU,   'standard_name', 'sea_surface_wave_orbital_velocity')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU,   'long_name'    , 'Orbital velocity')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU,   'units'        , 'm/s')
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_WU, '_FillValue', dmiss)

               if (kmx==0) then
                   ierr = nf90_def_var(hisids%ncid, 'ustokes',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_USTX)
                   ierr = nf90_put_att(hisids%ncid, hisids%id_USTX,   'coordinates'  , statcoordstring)
                   ierr = nf90_def_var(hisids%ncid, 'vstokes',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_USTY)
                   ierr = nf90_put_att(hisids%ncid, hisids%id_USTY,   'coordinates'  , statcoordstring)
               else
                   ierr = nf90_def_var(hisids%ncid, 'ustokes',  nf90_double, ((/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /)) , hisids%id_USTX)
                   ierr = nf90_put_att(hisids%ncid, hisids%id_USTX,   'coordinates'  , trim(statcoordstring) // ' zcoordinate_c')
                   ierr = nf90_def_var(hisids%ncid, 'vstokes',  nf90_double, ((/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /)) , hisids%id_USTY)
                   ierr = nf90_put_att(hisids%ncid, hisids%id_USTY,   'coordinates'  , trim(statcoordstring) // ' zcoordinate_c')
                   jawrizc = 1
               endif

               ierr = nf90_put_att(hisids%ncid, hisids%id_USTX,   'standard_name', 'sea_surface_wave_stokes_drift_x')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTX,   'long_name'    , 'Stokes drift, x-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTX,   'units'        , 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTX, '_FillValue', dmiss)

               ierr = nf90_put_att(hisids%ncid, hisids%id_USTY,   'standard_name', 'sea_surface_wave_stokes_drift_y')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTY,   'long_name'    , 'Stokes drift, y-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTY,   'units'        , 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_USTY, '_FillValue', dmiss)

            endif

            if (jahistaucurrent>0 .and. janormalhis > 0) then
               
               ierr = nf90_def_var(hisids%ncid, 'tausx',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_TAUX)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, 'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, 'standard_name', 'mean_bottom_shear_stress vector, x-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, 'long_name'    , 'Mean bottom shear stress vector, x-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, 'units'        , 'Pa')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUX, '_FillValue', dmiss)
               
               ierr = nf90_def_var(hisids%ncid, 'tausy',  nf90_double, ((/ hisids%id_statdim, hisids%id_timedim /)) , hisids%id_TAUY)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, 'coordinates'  , statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, 'standard_name', 'mean_bottom_shear_stress vector, y-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, 'long_name'    , 'Mean bottom shear stress vector, y-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, 'units'        , 'Pa')
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_TAUY, '_FillValue', dmiss)
            endif

            if (jatem > 0 .and. jahistem > 0) then
               if ( kmx.gt.0 ) then
                  if (ja3d > 0) then
                     ierr = nf90_def_var(hisids%ncid, 'temperature', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_vartem)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     jawrizc = 1
                  else
                     ierr = nf90_def_var(hisids%ncid, 'depth-averaged_temperature', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_vartem)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'coordinates', trim(statcoordstring))
                     jawrizc = 1
                  end if
               else
                  ierr = nf90_def_var(hisids%ncid, 'temperature', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_vartem)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'coordinates', statcoordstring)
               end if
               ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'units', 'degC')
               ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'geometry', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_vartem, 'standard_name', 'sea_water_temperature')

               if (jatem > 1 .and. jahisheatflux > 0 .and. janormalhis > 0) then ! here less verbose
                  idims(1) = hisids%id_statdim
                  idims(2) = hisids%id_timedim
                  call definencvar(hisids%ncid,hisids%id_wind   ,nf90_double,idims,2, 'wind'  , 'windspeed', 'm s-1', statcoordstring, station_geom_container_name)
                  call definencvar(hisids%ncid,hisids%id_tair   ,nf90_double,idims,2, 'Tair'  , 'air temperature', 'degC', statcoordstring, station_geom_container_name)
                  if (jatem == 5) then
                     call definencvar(hisids%ncid,hisids%id_rhum   ,nf90_double,idims,2, 'rhum'  , 'relative humidity', ' ',statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_clou   ,nf90_double,idims,2, 'clou'  , 'cloudiness', ' ', statcoordstring, station_geom_container_name)

                     call definencvar(hisids%ncid,hisids%id_qsun   ,nf90_double,idims,2, 'Qsun'  , 'solar influx', 'W m-2', statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_Qeva   ,nf90_double,idims,2, 'Qeva'  , 'evaporative heat flux', 'W m-2', statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_Qcon   ,nf90_double,idims,2, 'Qcon'  , 'sensible heat flux', 'W m-2', statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_Qlong  ,nf90_double,idims,2, 'Qlong' , 'long wave back radiation', 'W m-2', statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_Qfreva ,nf90_double,idims,2, 'Qfreva', 'free convection evaporative heat flux', 'W m-2', statcoordstring, station_geom_container_name)
                     call definencvar(hisids%ncid,hisids%id_Qfrcon ,nf90_double,idims,2, 'Qfrcon', 'free convection sensible heat flux', 'W m-2', statcoordstring, station_geom_container_name)
                  endif
                  if (jatem > 1) then
                     call definencvar(hisids%ncid,hisids%id_Qtot   ,nf90_double,idims,2, 'Qtot'  , 'total heat flux', 'W m-2', statcoordstring, station_geom_container_name)
                  end if
               endif
            endif

            if ((jasal > 0 .or. jatem > 0 .or. jased > 0) .and. jahisrho > 0) then
               if (kmx.gt.0) then
                  if (ja3d > 0) then
                     idims(1) = hisids%id_laydim
                     jawrizc = 1
                     call definencvar(hisids%ncid,hisids%id_varrhop,nf90_double, idims,3,'potential_density'       , 'potential_density'         , 'kg m-3' ,  trim(statcoordstring) // ' zcoordinate_c', station_geom_container_name, fillVal = dmiss)
                     if (idensform > 10) then 
                     call definencvar(hisids%ncid,hisids%id_varrho ,nf90_double, idims,3,'density'                 , 'density'                   , 'kg m-3' ,  trim(statcoordstring) // ' zcoordinate_c', station_geom_container_name, fillVal = dmiss)
                     endif
                     idims(1) = hisids%id_laydimw
                     call definencvar(hisids%ncid,hisids%id_bruv   ,nf90_double, idims,3,'Brunt_Vaisala_N2'  , 'Brunt_Vaisala_N2'  , '1/s2'   ,  trim(statcoordstring) // ' zcoordinate_w', station_geom_container_name, fillVal = dmiss)
                  end if
               else
                  call definencvar(hisids%ncid,hisids%id_varrhop,nf90_double, idims,2,'potential_density'       , 'potential_density'         , 'kg m-3' ,  trim(statcoordstring)                    , station_geom_container_name, fillVal = dmiss)
               end if  ! en misschien kan iemand de kmx > 0 repetitie ook in definencvar afvangen en // ' zcoordinate_c' en idims,3  automatisch maken voor kmx>0
            endif

            if (ITRA1 > 0 .and. jahisconst > 0) then
               do j=ITRA1,ITRAN
                  i = j-ITRA1+1  ! tracer nr

                  tmpstr = const_names(j)
                  ! Forbidden chars in NetCDF names: space, /, and more.
                  call replace_char(tmpstr,32,95)

                  if ( kmx > 0 ) then
                     if (ja3d > 0) then
                        ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_tra(i))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                        jawrizc = 1
                     else
                        ierr = nf90_def_var(hisids%ncid, 'depth-averaged_'//trim(tmpstr), nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_tra(i))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'coordinates', trim(statcoordstring))
                     end if
                  else
                     ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_tra(i))
                     ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'coordinates', statcoordstring)
                  end if
                  if (const_units(j).ne.' ') then
                     tmpstr = const_units(j)
                  else
                     tmpstr = '1e-3'
                  endif
                  ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'units', tmpstr)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), '_FillValue', dmiss)
                  if (kmx > 0 .and. ja3d > 0) then
                     ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'long_name', 'Depth-averaged '//trim(const_names(j))//' concentration')
                  else
                     ierr = nf90_put_att(hisids%ncid, hisids%id_tra(i), 'long_name', trim(const_names(j))//' concentration')
                  end if
               enddo
            endif

            if (numwqbots > 0) then
               call realloc(hisids%id_hwqb, numwqbots, keepExisting = .false.)
               do i=1,numwqbots
                  tmpstr = wqbotnames(i)
                  ! Forbidden chars in NetCDF names: space, /, and more.
                  call replace_char(tmpstr,32,95)

                  ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_hwqb(i))
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb(i), 'coordinates', statcoordstring)

                  tmpstr = wqbotunits(i)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb(i), 'units', tmpstr)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb(i), 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb(i), '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb(i), 'long_name', wqbotnames(i))
               enddo
               if (wqbot3D_output == 1 .and. ja3d > 0) then
                  call realloc(hisids%id_hwqb3d, numwqbots, keepExisting = .false.)
                  do i=1,numwqbots
                     tmpstr = wqbotnames(i)
                     ! Forbidden chars in NetCDF names: space, /, and more.
                     call replace_char(tmpstr,32,95)

                     ierr = nf90_def_var(hisids%ncid, trim(tmpstr)//'_3D', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_hwqb3d(i))
                     ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb3d(i), 'coordinates', statcoordstring)

                     tmpstr = wqbotunits(i)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb3d(i), 'units', tmpstr)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb3d(i), 'geometry', station_geom_container_name)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb3d(i), '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_hwqb3d(i), 'long_name', trim(wqbotnames(i))//' (3D)')
                  enddo
              endif
            endif

!          waq output
             if(jawaqproc > 0) then
               if (noout_user > 0) then
                  call realloc(hisids%id_hwq, noout_user, keepExisting = .false.)
                  do j=1,noout_user
                     tmpstr = ' '
                     write (tmpstr, "('water_quality_output_',I0)") j
                     if ( kmx > 0 ) then  !        3D
                        if (ja3d > 0) then
                           ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_hwq(j))
                           ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                           jawrizc = 1
                        end if
                     else
                        ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_hwq(j))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'coordinates', statcoordstring)
                     end if
                     if (hisids%id_hwq(j) > 0) then
                        tmpstr = trim(outputs%names(j))//' - '//trim(outputs%descrs(j))//' in flow element'
                        call replace_multiple_spaces_by_single_spaces(tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), '_FillValue', dmiss)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'long_name', trim(outputs%names(j)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'units', trim(outputs%units(j)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'description', tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(j), 'geometry', station_geom_container_name)
                     end if
                  enddo
               endif
               if (noout_statt > 0) then
                  call realloc(hisids%id_hwq, noout_user + noout_statt, keepExisting = .true.)
                  do j=1,noout_statt
                     jj = noout_user + j
                     tmpstr = ' '
                     write (tmpstr, "('water_quality_stat_',I0)") j
                     if ( kmx > 0 ) then  !        3D
                        if (ja3d > 0) then
                           ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_hwq(jj))
                           ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                           jawrizc = 1
                        end if
                     else
                        ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_hwq(jj))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'coordinates', statcoordstring)
                     end if
                     if (hisids%id_hwq(jj) > 0) then
                        tmpstr = trim(outputs%names(jj))//' - '//trim(outputs%descrs(jj))//' in flow element'
                        call replace_multiple_spaces_by_single_spaces(tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), '_FillValue', dmiss)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'long_name', trim(outputs%names(jj)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'units', trim(outputs%units(jj)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'description', tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'geometry', station_geom_container_name)
                     end if
                  enddo
               endif
               if (noout_state> 0) then
                  call realloc(hisids%id_hwq, noout, keepExisting = .true.)
                  do j=1,noout_state
                     jj = noout_user + noout_statt + j
                     tmpstr = ' '
                     write (tmpstr, "('water_quality_stat_',I0)") noout_statt + j
                     if ( kmx > 0 ) then  !        3D
                        if (ja3d > 0) then
                           ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_laydim, hisids%id_statdim /), hisids%id_hwq(jj))
                           ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'coordinates', statcoordstring)
                           jawrizc = 1
                        end if
                     else
                        ierr = nf90_def_var(hisids%ncid, trim(tmpstr), nf90_double, (/ hisids%id_statdim /), hisids%id_hwq(jj))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'coordinates', statcoordstring)
                     end if
                     if (hisids%id_hwq(jj) > 0) then
                        tmpstr = trim(outputs%names(jj))//' - '//trim(outputs%descrs(jj))//' in flow element'
                        call replace_multiple_spaces_by_single_spaces(tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), '_FillValue', dmiss)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'long_name', trim(outputs%names(jj)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'units', trim(outputs%units(jj)))
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'description', tmpstr)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_hwq(jj), 'geometry', station_geom_container_name)
                     end if
                  enddo
               endif
            endif

            if (stm_included .and. ISED1 > 0 .and. jahissed > 0) then
               ! New implementation, sedsus fraction is additional dimension
               ierr = nf90_def_dim(hisids%ncid, 'nSedTot', stmpar%lsedtot, hisids%id_sedtotdim)
               ierr = nf90_def_dim(hisids%ncid, 'nSedSus', stmpar%lsedsus, hisids%id_sedsusdim)
               !
               ierr = nf90_def_var(hisids%ncid, 'sedfrac_name', nf90_char, (/ hisids%id_strlendim, hisids%id_sedtotdim /), hisids%id_frac_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_frac_name,'long_name', 'sediment fraction identifier')
               !
               if (kmx>0) then
                  if (ja3d > 0) then ! NOTE: UNST-1765: No depth-averaged values for old-style 3d sed yet.
                     ierr = nf90_def_var(hisids%ncid, 'sed', nf90_double, (/  hisids%id_laydim, hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_sf)
                     ierr = nf90_def_var(hisids%ncid, 'ws', nf90_double, (/  hisids%id_laydimw, hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_ws)
                     ierr = nf90_def_var(hisids%ncid, 'seddif', nf90_double, (/  hisids%id_laydimw, hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_seddif)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_seddif, 'long_name', 'Sediment vertical diffusion')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_seddif, 'units', 'm2 s-1')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_seddif, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_seddif, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_seddif, 'geometry', station_geom_container_name)
                     !
                     jawrizc = 1
                     jawrizw = 1
                  end if
               else
                  ierr = nf90_def_var(hisids%ncid, 'sed', nf90_double, (/  hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_sf)
                  ierr = nf90_def_var(hisids%ncid, 'ws', nf90_double, (/ hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_ws)
               endif
               !
               if (hisids%id_sf > 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sf, 'long_name', 'Sediment concentration')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sf, 'units', 'kg m-3')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sf, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sf, 'geometry', station_geom_container_name)
               endif
               !
               if (hisids%id_ws > 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_ws, 'long_name', 'Sediment settling velocity')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_ws, 'units', 'm s-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_ws, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_ws, 'geometry', station_geom_container_name)
               end if
               !
            endif
            !
            ! Sediment transports
            !
            if (jased > 0 .and. stm_included .and. jahissed > 0 .and. janormalhis > 0) then
               !
               if (stmpar%morpar%moroutput%taub) then
                  ierr = nf90_def_var(hisids%ncid, 'taub', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_taub)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_taub, 'long_name', 'Bed shear stress for morphology')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_taub, 'units', 'Pa')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_taub, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_taub, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_taub, 'geometry', station_geom_container_name)
               endif
               !
               select case(stmpar%morpar%moroutput%transptype)
                  case (0)
                     transpunit = 'kg s-1 m-1'
                  case (1)
                     transpunit = 'm3 s-1 m-1'
                  case (2)
                     transpunit = 'm3 s-1 m-1'
               end select
               !
               if (stmpar%morpar%moroutput%sbcuv) then
                  ierr = nf90_def_var(hisids%ncid, 'sbcx', nf90_double, (/  hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sbcx)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcx, 'long_name', 'Current related bedload transport, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcx, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcx, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcx, 'geometry', station_geom_container_name)
                  ierr = nf90_def_var(hisids%ncid, 'sbcy', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sbcy)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcy, 'long_name', 'Current related bedload transport, y-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcy, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcy, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbcy, 'geometry', station_geom_container_name)
               endif
               if (stmpar%morpar%moroutput%sbwuv .and. jawave>0 .and. .not. flowWithoutWaves) then
                  ierr = nf90_def_var(hisids%ncid, 'sbwx', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sbwx)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwx, 'long_name', 'Wave related bedload transport, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwx, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwx, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwx, 'geometry', station_geom_container_name)
                  ierr = nf90_def_var(hisids%ncid, 'sbwy', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sbwy)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwy, 'long_name', 'Wave related bedload transport, y-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwy, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwy, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sbwy, 'geometry', station_geom_container_name)
               endif
               if (stmpar%morpar%moroutput%sswuv .and. jawave>0 .and. .not. flowWithoutWaves) then
                  ierr = nf90_def_var(hisids%ncid, 'sswx', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sswx)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswx, 'long_name', 'Wave related suspended transport, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswx, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswx, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswx, 'geometry', station_geom_container_name)
                  ierr = nf90_def_var(hisids%ncid, 'sswy', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sswy)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswy, 'long_name', 'Wave related suspended transport, y-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswy, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswy, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sswy, 'geometry', station_geom_container_name)
               end if
               if (stmpar%morpar%moroutput%sscuv) then
                  ierr = nf90_def_var(hisids%ncid, 'sscx', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sscx)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscx, 'long_name', 'Current related suspended transport, x-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscx, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscx, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscx, 'geometry', station_geom_container_name)
                  ierr = nf90_def_var(hisids%ncid, 'sscy', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_sscy)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscy, 'long_name', 'Current related suspended transport, y-component')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscy, 'units', transpunit)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscy, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sscy, 'geometry', station_geom_container_name)
               endif
               !
               ! Source and sink terms
               !
               if (stmpar%morpar%moroutput%sourcesink) then
                  ierr = nf90_def_var(hisids%ncid, 'sourse', nf90_double, (/ hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_sour)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sour, 'long_name', 'Source term suspended sediment transport')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sour, 'units', 'kg m-3 s-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sour, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sour, 'geometry', station_geom_container_name)

                  ierr = nf90_def_var(hisids%ncid, 'sinkse', nf90_double, (/ hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_sink)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sink, 'long_name', 'Sink term suspended sediment transport')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sink, 'units', 's-1')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sink, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sink, 'geometry', station_geom_container_name)
               endif
               !
               ! Bed composition variables
               !
               select case (stmpar%morlyr%settings%iunderlyr)
                  case (1)
                     ierr = nf90_def_var(hisids%ncid, 'bodsed', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_bodsed)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_bodsed, 'long_name', 'Available sediment mass in the bed')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_bodsed, 'units', 'kg m-2')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_bodsed, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_bodsed, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_bodsed, 'geometry', station_geom_container_name)
                     !
                     ierr = nf90_def_var(hisids%ncid, 'dpsed', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_dpsed)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_dpsed, 'long_name', 'Sediment thickness in the bed')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_dpsed, 'units', 'm')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_dpsed, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_dpsed, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_dpsed, 'geometry', station_geom_container_name)
                  case (2)
                     ierr = nf90_def_dim(hisids%ncid, 'nBedLayers', stmpar%morlyr%settings%nlyr, hisids%id_nlyrdim)
                     !
                     ierr = nf90_def_var(hisids%ncid, 'msed', nf90_double, (/ hisids%id_nlyrdim, hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_msed)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_msed, 'long_name', 'Available sediment mass in a layer of the bed')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_msed, 'units', 'kg m-2')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_msed, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_msed, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_msed, 'geometry', station_geom_container_name)
                     !
                     ierr = nf90_def_var(hisids%ncid, 'thlyr', nf90_double, (/ hisids%id_nlyrdim, hisids%id_statdim, hisids%id_timedim /), hisids%id_thlyr)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_thlyr, 'long_name', 'Thickness of a layer of the bed')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_thlyr, 'units', 'm')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_thlyr, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_thlyr, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_thlyr, 'geometry', station_geom_container_name)
                     !
                     if (stmpar%morlyr%settings%iporosity>0) then
                        ierr = nf90_def_var(hisids%ncid, 'poros', nf90_double, (/ hisids%id_nlyrdim, hisids%id_statdim, hisids%id_timedim /), hisids%id_poros)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_poros, 'long_name', 'Porosity of a layer of the bed')
                        ierr = nf90_put_att(hisids%ncid, hisids%id_poros, 'units', '-')
                        ierr = nf90_put_att(hisids%ncid, hisids%id_poros, '_FillValue', dmiss)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_poros, 'coordinates', statcoordstring)
                        ierr = nf90_put_att(hisids%ncid, hisids%id_poros, 'geometry', station_geom_container_name)
                      endif
                     !
                     ierr = nf90_def_var(hisids%ncid, 'lyrfrac', nf90_double, (/ hisids%id_nlyrdim, hisids%id_statdim,  hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_lyrfrac)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_lyrfrac, 'long_name', 'Volume fraction in a layer of the bed')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_lyrfrac, 'units', 'm')
                     ierr = nf90_put_att(hisids%ncid, hisids%id_lyrfrac, '_FillValue', dmiss)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_lyrfrac, 'coordinates', statcoordstring)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_lyrfrac, 'geometry', station_geom_container_name)
               end select
		         !
		         if (stmpar%morpar%moroutput%frac) then
                  ierr = nf90_def_var(hisids%ncid, 'frac', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_frac)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_frac, 'long_name', 'Availability fraction in top layer')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_frac, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_frac, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_frac, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_frac, 'geometry', station_geom_container_name)
               endif
               !
               if (stmpar%morpar%moroutput%mudfrac) then
                  ierr = nf90_def_var(hisids%ncid, 'mudfrac', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_mudfrac)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mudfrac, 'long_name', 'Mud fraction in top layer')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mudfrac, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mudfrac, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mudfrac, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mudfrac, 'geometry', station_geom_container_name)
               endif
               !
               if (stmpar%morpar%moroutput%sandfrac) then
                  ierr = nf90_def_var(hisids%ncid, 'sandfrac', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_sandfrac)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sandfrac, 'long_name', 'Sand fraction in top layer')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sandfrac, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sandfrac, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sandfrac, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sandfrac, 'geometry', station_geom_container_name)
               endif
               !
               if (stmpar%morpar%moroutput%fixfac) then
                  ierr = nf90_def_var(hisids%ncid, 'fixfac', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_fixfac)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_fixfac, 'long_name', 'Reduction factor due to limited sediment thickness')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_fixfac, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_fixfac, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_fixfac, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_fixfac, 'geometry', station_geom_container_name)
               endif
               !
               if (stmpar%morpar%moroutput%hidexp) then
                  ierr = nf90_def_var(hisids%ncid, 'hidexp', nf90_double, (/ hisids%id_statdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_hidexp)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hidexp, 'long_name', 'Hiding and exposure factor')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hidexp, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hidexp, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hidexp, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_hidexp, 'geometry', station_geom_container_name)
               endif
               !
               if (stmpar%morpar%flufflyr%iflufflyr>0 .and. stmpar%lsedsus>0) then
                  ierr = nf90_def_var(hisids%ncid, 'mfluff', nf90_double, (/ hisids%id_statdim, hisids%id_sedsusdim, hisids%id_timedim /), hisids%id_mfluff)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mfluff, 'long_name', 'Sediment mass in fluff layer')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mfluff, 'units', '-')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mfluff, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mfluff, 'coordinates', statcoordstring)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_mfluff, 'geometry', station_geom_container_name)
               end if
            endif
            !
            if (jased > 0 .and. .not. stm_included .and. jahissed > 0) then
               if ( kmx.gt.0 ) then
                  if (ja3d > 0) then ! NOTE: UNST-1765: No depth-averaged values for old-style 3d sed yet.
                     ierr = nf90_def_var(hisids%ncid, 'sediment_concentration', nf90_double, (/ hisids%id_laydim, hisids%id_statdim, hisids%id_timedim /), hisids%id_varsed)
                     ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, 'coordinates', trim(statcoordstring) // ' zcoordinate_c')
                     jawrizc = 1
                  end if
               else
                  ierr = nf90_def_var(hisids%ncid, 'sediment_concentration', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varsed)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, 'coordinates', statcoordstring)
               end if
               if (hisids%id_varsed > 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, 'units', 'kg m-3')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, 'geometry', station_geom_container_name)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, '_FillValue', dmiss)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varsed, 'long_name', 'sediment_concentration')
               end if
            endif

            if (japatm > 0 .and. jahiswind > 0 .and. janormalhis > 0) then
               call definencvar(hisids%ncid,hisids%id_varpatm   ,nf90_double,(/ hisids%id_statdim, hisids%id_timedim /),2, 'patm'  , 'atmospheric pressure', 'N m-2', statcoordstring, station_geom_container_name)
            endif

            if (jawind > 0 .and. jahiswind > 0 .and. janormalhis > 0) then
               ierr = nf90_def_var(hisids%ncid, 'windx', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varwx)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'geometry', station_geom_container_name)
               ierr = nf90_def_var(hisids%ncid, 'windy', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varwy)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'geometry', station_geom_container_name)
               if (jsferic == 0) then
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'standard_name', 'x_wind')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'standard_name', 'y_wind')
               else
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'standard_name', 'eastward_wind')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'standard_name', 'northward_wind')
               end if
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwx, 'long_name', 'velocity of air on flow element center, x-component')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varwy, 'long_name', 'velocity of air on flow element center, y-component')
            endif

            if (jarain > 0 .and. jahisrain > 0 .and. janormalhis > 0) then
               ierr = nf90_def_var(hisids%ncid, 'rain', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_varrain)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, 'units', 'mm day-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, 'standard_name', 'lwe_precipitation_rate')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, 'long_name', 'precipitation depth per time unit')
               ierr = nf90_put_att(hisids%ncid, hisids%id_varrain, 'geometry', station_geom_container_name)
            endif

            if ((infiltrationmodel == DFM_HYD_INFILT_CONST .or. infiltrationmodel == DFM_HYD_INFILT_HORTON) .and. jahisinfilt > 0 .and. janormalhis > 0) then
               ierr = nf90_def_var(hisids%ncid, 'infiltration_cap', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_infiltcap)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, 'units', 'mm hr-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, 'coordinates', statcoordstring)
               !ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, 'standard_name', 'infiltration_rate)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, 'long_name', 'Infiltration capacity')
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltcap, 'geometry', station_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'infiltration_actual', nf90_double, (/ hisids%id_statdim, hisids%id_timedim /), hisids%id_infiltact)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltact, 'units', 'mm hr-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltact, '_FillValue', dmiss)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltact, 'coordinates', statcoordstring)
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltact, 'long_name', 'Actual infiltration rate')
               ierr = nf90_put_att(hisids%ncid, hisids%id_infiltact, 'geometry', station_geom_container_name)
            endif

            if (kmx.gt.0 .and. jawrizc == 1 .and. ja3d > 0) then
               idims(1) = hisids%id_laydim
               idims(2) = hisids%id_statdim
               idims(3) = hisids%id_timedim
               call definencvar   (hisids%ncid, hisids%id_zcs, nf90_double, idims,3, 'zcoordinate_c' , 'vertical coordinate at center of flow element and layer'   , 'm',  'station_x_coordinate station_y_coordinate station_name zcoordinate_c', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_zcs, 'positive' , 'up')
               ierr = nf90_put_att(hisids%ncid, hisids%id_zcs, '_FillValue' , dmiss)
            endif
            if (kmx.gt.0 .and. jawrizw == 1 .and. ja3d > 0) then
               idims(1) = hisids%id_laydimw
               idims(2) = hisids%id_statdim
               idims(3) = hisids%id_timedim
               call definencvar   (hisids%ncid, hisids%id_zws, nf90_double, idims,3, 'zcoordinate_w' , 'vertical coordinate at centre of flow element and at layer interface'   , 'm',  'station_x_coordinate station_y_coordinate station_name zcoordinate_w', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_zws, 'positive' , 'up')
               ierr = nf90_put_att(hisids%ncid, hisids%id_zws, '_FillValue' , dmiss)
               call definencvar   (hisids%ncid, hisids%id_zwu, nf90_double, idims,3, 'zcoordinate_wu' , 'vertical coordinate at edge of flow element and at layer interface'   , 'm',  'station_x_coordinate station_y_coordinate station_name zcoordinate_wu', station_geom_container_name)
               ierr = nf90_put_att(hisids%ncid, hisids%id_zwu, 'positive' , 'up')
               ierr = nf90_put_att(hisids%ncid, hisids%id_zwu, '_FillValue' , dmiss)
            endif
        end if

        if (ncrs > 0 .and. janormalhis > 0) then
            mnp = 0
            do i=1,ncrs
                mnp = max(mnp, crs(i)%path%np)
            end do
            ierr = nf90_def_dim(hisids%ncid, 'cross_section', ncrs, hisids%id_crsdim)
            ierr = nf90_def_dim(hisids%ncid, 'cross_section_name_len', strlen_netcdf, hisids%id_crslendim)
            ierr = nf90_def_dim(hisids%ncid, 'cross_section_pts', mnp+1, hisids%id_crsptsdim)

            !ierr = nf90_def_var(hisids%ncid, 'cross_section_x_coordinate', nf90_double, (/ hisids%id_crsptsdim, hisids%id_crsdim /), hisids%id_crsx)
            !ierr = nf90_def_var(hisids%ncid, 'cross_section_y_coordinate', nf90_double, (/ hisids%id_crsptsdim, hisids%id_crsdim /), hisids%id_crsy)
            ierr = nf90_def_var(hisids%ncid, 'cross_section_name',         nf90_char,   (/ hisids%id_crslendim, hisids%id_crsdim /), hisids%id_crsname)
            ierr = nf90_put_att(hisids%ncid, hisids%id_crsname,  'cf_role', 'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_crsname,  'long_name', 'cross section name'    )
            !ierr = unc_addcoordatts(hisids%ncid, hisids%id_crsx, hisids%id_crsy, jsferic)

            ! Define geometry related variables
            crs_geom_container_name = 'cross_section_geom'
            nNodeTot = nNodesCrs

            ierr = sgeom_def_geometry_variables(hisids%ncid, crs_geom_container_name, 'cross section', 'line', nNodeTot, hisids%id_crsdim, &
               hisids%id_crsgeom_node_count, hisids%id_crsgeom_node_coordx, hisids%id_crsgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'cross_section_discharge',     nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varQ)
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQ,    'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQ,    'coordinates', 'cross_section_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQ, 'geometry', crs_geom_container_name)
            ierr = nf90_def_var(hisids%ncid, 'cross_section_cumulative_discharge', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varQint)
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQint, 'units', 'm3')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQint, 'coordinates', 'cross_section_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varQint, 'geometry', crs_geom_container_name)
            !ierr = nf90_def_var(hisids%ncid, 'cross_section_discharge_avg', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varQavg)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varQavg, 'units', 'm3 s-1')
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varQavg, 'coordinates', 'cross_section_name')

            ierr = nf90_def_var(hisids%ncid, 'cross_section_area',     nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varAu)
            ierr = nf90_put_att(hisids%ncid, hisids%id_varAu,    'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varAu,    'coordinates', 'cross_section_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varAu, 'geometry', crs_geom_container_name)
            !ierr = nf90_def_var(hisids%ncid, 'cross_section_area_avg', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varAuavg)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varAuavg, 'units', 'm2')
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varAuavg, 'coordinates', 'cross_section_name')

            ierr = nf90_def_var(hisids%ncid, 'cross_section_velocity',     nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varu)
            ierr = nf90_put_att(hisids%ncid, hisids%id_varu,    'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varu,    'coordinates', 'cross_section_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varu, 'geometry', crs_geom_container_name)
            ! Disable writing cross_section_velocity_avg (see UNST-1148), because in a parallel run, it is impossible to compute
            ! summation of area (denominator) at each computational time step in a cheap way, i.e. without communication between
            ! partitions. @see subroutines: sumvalueOnCrossSections, updateValuesOnCrossSections
            !ierr = nf90_def_var(hisids%ncid, 'cross_section_velocity_avg', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_varuavg)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varuavg, 'units', 'm s-1')
            !ierr = nf90_put_att(hisids%ncid, hisids%id_varuavg, 'coordinates', 'cross_section_name')

               do num = 1,NUMCONST_MDU
                  tmpstr = const_names(num)
                  ! Forbidden chars in NetCDF names: space, /, and more.
                  call replace_char(tmpstr,32,95)
                  call replace_char(tmpstr,47,95)
                  ierr = nf90_def_var(hisids%ncid, 'cross_section_cumulative_'//trim(tmpstr), nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_const_cum(num))
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const_cum(num), 'long_name', 'cumulative flux (based on upwind flow cell) for '//trim(tmpstr)//'.')
         
                  ierr = nf90_def_var(hisids%ncid, 'cross_section_'//trim(tmpstr), nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_const(num))
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const(num), 'long_name', 'flux (based on upwind flow cell) for '//trim(tmpstr)//'.')
         
                  if (num >= ISED1 .and. num <= ISEDN) then    ! if the constituent is sediment
                     select case(stmpar%morpar%moroutput%transptype)
                     case (0)
                        tmpstr = 'kg'
                     case (1, 2)
                        tmpstr = 'm3'
                     end select
                  else
                     if (const_units(num) /= ' ') then
                        tmpstr = trim(const_units(num)) // ' m3'
                     else
                        tmpstr = '-'
                     endif
                  endif
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const_cum(num), 'units', tmpstr)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const_cum(num), 'coordinates', 'cross_section_name')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const_cum(num), 'geometry', crs_geom_container_name)
         
                  if (num >= ISED1 .and. num <= ISEDN) then    ! if the constituent is sediment
                     select case(stmpar%morpar%moroutput%transptype)
                     case (0)
                        tmpstr = 'kg/s'
                     case (1, 2)
                        tmpstr = 'm3/s'
                     end select
                  else
                     if (const_units(num) /= ' ') then
                        tmpstr = trim(const_units(num)) // ' m3/s'
                     else
                        tmpstr = '-'
                     endif
                  endif
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const(num), 'units', tmpstr)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const(num), 'coordinates', 'cross_section_name')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_const(num), 'geometry', crs_geom_container_name)
               enddo
            endif
         
            if( jased == 4 .and. stmpar%lsedtot > 0 ) then
               ierr = nf90_def_var(hisids%ncid, 'cross_section_bedload_sediment_transport', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_sedbtrans)
               ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtrans, 'long_name', 'cumulative bed load sediment transport')
               ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtrans, 'units', 'kg')
               ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtrans, 'coordinates', 'cross_section_name')
               ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtrans, 'geometry', crs_geom_container_name)
               if( stmpar%lsedsus > 0 ) then
                  ierr = nf90_def_var(hisids%ncid, 'cross_section_suspended_sediment_transport', nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_sedstrans)
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedstrans, 'long_name', 'cumulative suspended load sediment transport')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedstrans, 'units', 'kg')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedstrans, 'coordinates', 'cross_section_name')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedstrans, 'geometry', crs_geom_container_name)
               endif
               if (.not. allocated(hisids%id_sedbtransfrac)) then
                  allocate(hisids%id_sedbtransfrac(stmpar%lsedtot))
                  hisids%id_sedbtransfrac = 0
               endif
               do lsed = 1,stmpar%lsedtot    ! Making bedload on crosssections per fraction
                  ierr = nf90_def_var(hisids%ncid, 'cross_section_bedload_sediment_transport_'//trim(stmpar%sedpar%namsed(lsed)), nf90_double, (/ hisids%id_crsdim, hisids%id_timedim /), hisids%id_sedbtransfrac(lsed))
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtransfrac(lsed), 'long_name', 'cumulative bed load sediment transport per fraction')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtransfrac(lsed), 'units', 'kg')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtransfrac(lsed), 'coordinates', 'cross_section_name')
                  ierr = nf90_put_att(hisids%ncid, hisids%id_sedbtransfrac(lsed), 'geometry', crs_geom_container_name)
               enddo
         
            endif
        

        ! runup gauges
        if (nrug > 0 .and. janormalhis > 0) then
           ierr = nf90_def_dim(hisids%ncid, 'runupgauges', nrug, hisids%id_rugdim)

           ierr = nf90_def_var(hisids%ncid, 'rug_x_coordinate', nf90_double, (/ hisids%id_rugdim, hisids%id_timedim /), hisids%id_rugx)
           ierr = nf90_def_var(hisids%ncid, 'rug_y_coordinate', nf90_double, (/ hisids%id_rugdim, hisids%id_timedim /), hisids%id_rugy)

            ierr = unc_addcoordatts(hisids%ncid, hisids%id_rugx, hisids%id_rugy, jsferic)
            ierr = nf90_put_att(hisids%ncid, hisids%id_rugx, 'long_name', 'time-varying x-coordinate of shoreline position')
            ierr = nf90_put_att(hisids%ncid, hisids%id_rugy, 'long_name', 'time-varying y-coordinate of shoreline position')

            ierr = nf90_def_var(hisids%ncid, 'rug_id', nf90_char,   (/ hisids%id_strlendim, hisids%id_rugdim /), hisids%id_rugid)
            ierr = nf90_put_att(hisids%ncid, hisids%id_rugid,  'long_name'    , 'runup gauge identifier') ! REF

            ierr = nf90_def_var(hisids%ncid, 'rug_name', nf90_char,   (/ hisids%id_strlendim, hisids%id_rugdim /), hisids%id_rugname)
            ierr = nf90_put_att(hisids%ncid, hisids%id_rugname,  'cf_role', 'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_rugname,  'long_name'    , 'runup gauge name') ! REF

            ierr = nf90_def_var(hisids%ncid, 'runup_height', nf90_double, (/ hisids%id_rugdim, hisids%id_timedim /), hisids%id_varruh)
            ierr = nf90_put_att(hisids%ncid, hisids%id_varruh, 'standard_name', 'runup_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varruh, 'long_name', 'runup height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varruh, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varruh, 'coordinates', 'rug_x_coordinate rug_y_coordinate rug_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_varruh, '_FillValue', dmiss)
        endif

        if (jahissourcesink > 0 .and. numsrc > 0 .and. janormalhis > 0) then
           ierr = nf90_def_dim(hisids%ncid, 'source_sink', numsrc, hisids%id_srcdim)
           ierr = nf90_def_dim(hisids%ncid, 'source_sink_name_len', strlen_netcdf, hisids%id_srclendim)
           ierr = nf90_def_dim(hisids%ncid, 'source_sink_pts', msrc, hisids%id_srcptsdim)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_name', nf90_char, (/ hisids%id_srclendim, hisids%id_srcdim/), hisids%id_srcname)
           ierr = nf90_put_att(hisids%ncid, hisids%id_srcname,  'cf_role', 'timeseries_id')
           ierr = nf90_put_att(hisids%ncid, hisids%id_srcname,  'long_name', 'source and sink name'    )

           ! Define geometry related variables
            src_geom_container_name = 'source_sink_geom'
            nNodeTot = 0
            do i = 1, numsrc
               nNodes = 0
               k1 = ksrc(1,i)
               k2 = ksrc(4,i)
               if (k1 /= 0) then
                  nNodes = nNodes + 1
               end if
               if (k2 /= 0) then
                  nNodes = nNodes + 1
               end if
               nNodeTot = nNodeTot + nNodes
            end do

            ierr = sgeom_def_geometry_variables(hisids%ncid, src_geom_container_name, 'source_sink', 'line', nNodeTot, hisids%id_srcdim, &
               hisids%id_srcgeom_node_count, hisids%id_srcgeom_node_coordx, hisids%id_srcgeom_node_coordy)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_x_coordinate', nf90_double, (/ hisids%id_srcdim, hisids%id_srcptsdim  /), hisids%id_srcx)
           ierr = nf90_def_var(hisids%ncid, 'source_sink_y_coordinate', nf90_double, (/ hisids%id_srcdim, hisids%id_srcptsdim /), hisids%id_srcy)
           ierr = unc_addcoordatts(hisids%ncid, hisids%id_srcx, hisids%id_srcy, jsferic)
           ierr = nf90_put_att(hisids%ncid, hisids%id_srcx, '_FillValue', dmiss)
           ierr = nf90_put_att(hisids%ncid, hisids%id_srcy, '_FillValue', dmiss)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_prescribed_discharge', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_pred)
           ierr = nf90_put_att(hisids%ncid, hisids%id_pred,    'units', 'm3 s-1')
           ierr = nf90_put_att(hisids%ncid, hisids%id_pred,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_pred,    'geometry', src_geom_container_name)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_prescribed_salinity_increment', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_presa)
           ierr = nf90_put_att(hisids%ncid, hisids%id_presa,    'units', '1e-3')
           ierr = nf90_put_att(hisids%ncid, hisids%id_presa,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_presa,    'geometry', src_geom_container_name)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_prescribed_temperature_increment', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_pretm)
           ierr = nf90_put_att(hisids%ncid, hisids%id_pretm,    'units', 'degC')
           ierr = nf90_put_att(hisids%ncid, hisids%id_pretm,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_pretm,    'geometry', src_geom_container_name)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_current_discharge', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_qsrccur)
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrccur,    'units', 'm3 s-1')
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrccur,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrccur,    'geometry', src_geom_container_name)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_cumulative_volume', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_vsrccum)
           ierr = nf90_put_att(hisids%ncid, hisids%id_vsrccum,    'long_name', 'Cumulative volume from the start time until current time at each source/sink')
           ierr = nf90_put_att(hisids%ncid, hisids%id_vsrccum,    'units', 'm3')
           ierr = nf90_put_att(hisids%ncid, hisids%id_vsrccum,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_vsrccum,    'geometry', src_geom_container_name)

           ierr = nf90_def_var(hisids%ncid, 'source_sink_discharge_average', nf90_double, (/ hisids%id_srcdim, hisids%id_timedim /), hisids%id_qsrcavg)
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrcavg,    'long_name', 'Average discharge in the past his-file output-interval at each source/sink')
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrcavg,    'units', 'm3 s-1')
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrcavg,    'coordinates', 'source_sink_name')
           ierr = nf90_put_att(hisids%ncid, hisids%id_qsrcavg,    'geometry', src_geom_container_name)
        endif

        if (timon) call timstrt ( "unc_write_his DEF bal", handle_extra(59))
        if (jahisbal > 0 .and. janormalhis > 0) then
           do num = 1,MAX_IDX
              if ( num.eq.IDX_InternalTidesDissipation ) then
                 if ( jaFrcInternalTides2D.eq.1 ) then
                    ierr = nf90_def_var(hisids%ncid, trim(voltotname(num)), nf90_double, (/ hisids%id_timedim /), hisids%id_voltot(num))
                    ierr = nf90_put_att(hisids%ncid, hisids%id_voltot(num), 'units', 'TJ')
                 end if
              else if ( num.eq.IDX_GravInput ) then
                 if ( jatidep.eq.1 ) then
                    ierr = nf90_def_var(hisids%ncid, trim(voltotname(num)), nf90_double, (/ hisids%id_timedim /), hisids%id_voltot(num))
                    ierr = nf90_put_att(hisids%ncid, hisids%id_voltot(num), 'units', 'TJ')
                 end if
              else if ( num.eq.IDX_SALInput .or. num.eq.IDX_SALInput2 ) then
                 if ( jaselfal.gt.0 ) then
                    ierr = nf90_def_var(hisids%ncid, trim(voltotname(num)), nf90_double, (/ hisids%id_timedim /), hisids%id_voltot(num))
                    ierr = nf90_put_att(hisids%ncid, hisids%id_voltot(num), 'units', 'TJ')
                 end if
              else
                 ierr = nf90_def_var(hisids%ncid, 'water_balance_'//trim(voltotname(num)), nf90_double, (/ hisids%id_timedim /), hisids%id_voltot(num))
                 ierr = nf90_put_att(hisids%ncid, hisids%id_voltot(num), 'units', 'm3')
              end if
           enddo
        end if
        if (timon) call timstop (handle_extra(59))

        if (timon) call timstrt ( "unc_write_his DEF structures", handle_extra(60))
        if (jaoldstr == 1) then
           ntmp = ncgensg
        else
           ntmp = ngenstru
        end if
        if(jahiscgen > 0 .and. ntmp > 0 .and. janormalhis > 0) then

            ierr = nf90_def_dim(hisids%ncid, 'general_structures', ntmp, hisids%id_genstrudim)
            ierr = nf90_def_var(hisids%ncid, 'general_structure_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_genstrudim /), hisids%id_genstru_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_id,  'long_name', 'Id of general structure'    )

            ! Define geometry related variables
            genstru_geom_container_name = 'general_structure_geom'
            nNodeTot = 0
            if (network%sts%numGeneralStructures > 0) then ! new general structure
               nNodeTot = nNodesGenstru
            else ! old general structure
               do n = 1, ngenstru
                  i = genstru2cgen(n)
                  nlinks = L2cgensg(i) - L1cgensg(i) + 1
                  if (nlinks > 0) then
                     nNodes = nlinks + 1
                  else if (nlinks == 0) then
                     nNodes = 0
                  end if
                  nNodeTot = nNodeTot + nNodes
               end do
            end if

            ierr = sgeom_def_geometry_variables(hisids%ncid, genstru_geom_container_name, 'general_structure', 'line', nNodeTot, hisids%id_genstrudim, &
               hisids%id_genstrugeom_node_count, hisids%id_genstrugeom_node_coordx, hisids%id_genstrugeom_node_coordy)


            ierr = nf90_def_var(hisids%ncid, 'general_structure_discharge',     nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_dis)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis, 'standard_name', 'integral_of_discharge_wrt_time') ! TODO: introduce time windows in nc
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis, 'long_name', 'Total discharge through general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_crest_level', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_crestl)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestl, 'standard_name', 'cgen_crest_level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestl, 'long_name', 'Crest level of general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestl, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestl, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_gate_lower_edge_level', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_edgel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_edgel, 'long_name', 'Gate lower edge level of general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_edgel, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_edgel, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_edgel, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_gate_opening_width', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_openw)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openw, 'long_name', 'Gate opening width of general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openw, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openw, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_s1up',     nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1up, 'long_name', 'Water level upstream of general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1up, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1up, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_s1dn',     nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1dn, 'long_name', 'Water level downstream of general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1dn, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1dn, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_head', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_head, 'long_name', 'Head difference across general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_head, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_head, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_flow_area', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au, 'long_name', 'Flow area at general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au, 'geometry', genstru_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'general_structure_velocity', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_vel, 'long_name', 'Velocity through general structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_vel, 'coordinates', 'general_structure_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_vel, 'geometry', genstru_geom_container_name)

            if (network%sts%numGeneralStructures > 0) then ! write extra fields for new general structure
               ierr = nf90_def_var(hisids%ncid, 'general_structure_crest_width', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_crestw)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestw, 'long_name', 'Crest width of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestw, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestw, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_crestw, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_discharge_through_gate_opening', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_dis_gate_open)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_open, 'long_name', 'Discharge through gate opening of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_open, 'units', 'm3 s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_open, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_open, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_discharge_over_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_dis_gate_over)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_over, 'long_name', 'Discharge over gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_over, 'units', 'm3 s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_over, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_over, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_discharge_under_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_dis_gate_under)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_under, 'long_name', 'Discharge under gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_under, 'units', 'm3 s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_under, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_dis_gate_under, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_gate_opening_height', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_openh)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openh, 'long_name', 'Gate opening height of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openh, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openh, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_openh, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_gate_upper_edge_level', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_uppl)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_uppl, 'long_name', 'Gate upper edge level of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_uppl, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_uppl, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_uppl, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_velocity_through_gate_opening', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_velgateopen)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateopen, 'long_name', 'Velocity through gate opening of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateopen, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateopen, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateopen, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_velocity_over_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_velgateover)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateover, 'long_name', 'Velocity over gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateover, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateover, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateover, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_velocity_under_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_velgateunder)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateunder, 'long_name', 'Velocity under gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateunder, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateunder, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_velgateunder, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_flow_area_in_gate_opening', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_au_open)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_open, 'long_name', 'Flow area in gate opening of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_open, 'units', 'm2')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_open, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_open, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_flow_area_over_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_au_over)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_over, 'long_name', 'Flow area over gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_over, 'units', 'm2')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_over, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_over, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_flow_area_under_gate', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_au_under)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_under, 'long_name', 'Flow area under gate of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_under, 'units', 'm2')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_under, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_au_under, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_state', nf90_int, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_stat)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'long_name', 'Flow state at general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'units', '-')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'units', '-')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'flag_values', (/ 0, 1, 2, 3, 4 /))
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'flag_meanings', 'no_flow weir_free weir_submerged gate_free gate_submerged')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'valid_range', (/ 0, 4 /))
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, '_FillValue', int(dmiss))
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_stat, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_s1_on_crest', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_s1crest)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1crest, 'long_name', 'Water level on crest of general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1crest, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1crest, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_s1crest, 'geometry', genstru_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'general_structure_force_difference', nf90_double, (/ hisids%id_genstrudim, hisids%id_timedim /), hisids%id_genstru_forcedif)
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_forcedif, 'long_name', 'Force difference per unit at general structure')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_forcedif, 'units', 'N m-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_forcedif, 'coordinates', 'general_structure_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_genstru_forcedif, 'geometry', genstru_geom_container_name)
            end if
        endif

        if(jahispump > 0 .and. npumpsg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'pumps', npumpsg, hisids%id_pumpdim)
            ierr = nf90_def_var(hisids%ncid, 'pump_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_pumpdim /), hisids%id_pump_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_id,  'long_name', 'Id of pump'    )

            ! Define geometry related variables
            pump_geom_container_name = 'pump_geom'
            nNodeTot = 0
            if (network%sts%numPumps > 0) then ! newpump
               nNodeTot = nNodesPump
            else ! old pump
               do n = 1, npumpsg
                  nlinks = L2pumpsg(n) - L1pumpsg(n) + 1
                  if (nlinks > 0) then
                     nNodes = nlinks + 1
                  else if (nlinks == 0) then
                     nNodes = 0
                  end if
                  nNodeTot = nNodeTot + nNodes
               end do
            end if

            ierr = sgeom_def_geometry_variables(hisids%ncid, pump_geom_container_name, 'pump', 'line', nNodeTot, hisids%id_pumpdim, &
               hisids%id_pumpgeom_node_count, hisids%id_pumpgeom_node_coordx, hisids%id_pumpgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'pump_xmid', nf90_double, (/ hisids%id_pumpdim /), hisids%id_pump_xmid)
            ierr = nf90_def_var(hisids%ncid, 'pump_ymid', nf90_double, (/ hisids%id_pumpdim /), hisids%id_pump_ymid)
            ierr = unc_addcoordatts(hisids%ncid, hisids%id_pump_xmid, hisids%id_pump_ymid, jsferic)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_xmid, 'long_name', 'x-coordinate of representative mid point of pump location (snapped polyline)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_ymid, 'long_name', 'y-coordinate of representative mid point of pump location (snapped polyline)')

            ierr = nf90_def_var(hisids%ncid, 'pump_structure_discharge', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_dis, 'long_name', 'Discharge through pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_dis, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_dis, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_capacity', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_cap)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_cap, 'long_name', 'Capacity of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_cap, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_cap, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_cap, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_discharge_dir', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_disdir)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_disdir, 'long_name', 'Discharge of pump w.r.t. pump orientation')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_disdir, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_disdir, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_disdir, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_s1up',     nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_s1up)    ! Nabi
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1up, 'long_name', 'Water level upstream of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1up, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1up, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_s1dn',     nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1dn, 'long_name', 'Water level downstream of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1dn, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1dn, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_structure_head', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_struhead)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_struhead, 'long_name', 'Head difference across pump structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_struhead, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_struhead, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_struhead, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_actual_stage',    nf90_int, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_stage)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_stage, 'long_name', 'Actual stage of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_stage, 'units', '-')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_stage, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_stage, '_FillValue', int(dmiss))
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_stage, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_head', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_head, 'long_name', 'Head difference in pumping direction')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_head, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_head, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_reduction_factor',    nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_redufact)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_redufact, 'long_name', 'Reduction factor of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_redufact, 'units', '-')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_redufact, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_redufact, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_s1_delivery_side', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_s1del)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1del, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1del, 'long_name', 'Water level at delivery side of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1del, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1del, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1del, 'geometry', pump_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'pump_s1_suction_side', nf90_double, (/ hisids%id_pumpdim, hisids%id_timedim /), hisids%id_pump_s1suc)
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1suc, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1suc, 'long_name', 'Water level at suction side of pump')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1suc, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1suc, 'coordinates', 'pump_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_pump_s1suc, 'geometry', pump_geom_container_name)
        endif

        if(jahisgate > 0 .and. ngatesg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'gates', ngatesg, hisids%id_gatedim)
            ierr = nf90_def_var(hisids%ncid, 'gate_name',  nf90_char,   (/ hisids%id_strlendim, hisids%id_gatedim /), hisids%id_gatename)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gatename,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gatename,  'long_name', 'gate name'    )

            ierr = nf90_def_var(hisids%ncid, 'gate_discharge', nf90_double, (/ hisids%id_gatedim, hisids%id_timedim /), hisids%id_gate_dis)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_gate_dis, 'standard_name', 'gate_discharge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_dis, 'long_name', 'gate discharge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_dis, 'coordinates', 'gate_name')

            ierr = nf90_def_var(hisids%ncid, 'gate_lower_edge_level', nf90_double, (/ hisids%id_gatedim, hisids%id_timedim /), hisids%id_gate_edgel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_edgel, 'standard_name', 'gate_lower_edge_level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_edgel, 'long_name', 'gate lower edge level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_edgel, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_edgel, 'coordinates', 'gate_name')

            ierr = nf90_def_var(hisids%ncid, 'gate_s1up',     nf90_double, (/ hisids%id_gatedim, hisids%id_timedim /), hisids%id_gate_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1up, 'long_name', 'gate water level up')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1up, 'coordinates', 'gate_name')

            ierr = nf90_def_var(hisids%ncid, 'gate_s1dn',     nf90_double, (/ hisids%id_gatedim, hisids%id_timedim /), hisids%id_gate_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1dn, 'long_name', 'gate water level down')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gate_s1dn, 'coordinates', 'gate_name')
        endif

         if(jahisgate > 0 .and. ngategen > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'gategens', ngategen, hisids%id_gategendim)
            ierr = nf90_def_var(hisids%ncid, 'gategen_name',  nf90_char,   (/ hisids%id_strlendim, hisids%id_gategendim /), hisids%id_gategenname)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategenname,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategenname,  'long_name', 'gate name'    )

            ! Define geometry related variables
            gategen_geom_container_name = 'gategen_geom'
            nNodeTot = 0
            do n = 1, ngategen
               i = gate2cgen(n)
               nlinks = L2cgensg(i) - L1cgensg(i) + 1
               if (nlinks > 0) then
                  nNodes = nlinks + 1
               else if (nlinks == 0) then
                  nNodes = 0
               end if
               nNodeTot = nNodeTot + nNodes
            end do
            ierr = sgeom_def_geometry_variables(hisids%ncid, gategen_geom_container_name, 'gategen', 'line', nNodeTot, hisids%id_gategendim, &
               hisids%id_gategengeom_node_count, hisids%id_gategengeom_node_coordx, hisids%id_gategengeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'gategen_discharge',     nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_dis)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_dis, 'standard_name', 'integral_of_discharge_wrt_time') ! TODO: introduce time windows in nc
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_dis, 'long_name', 'gate discharge (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_dis, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_dis, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_crest_level', nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_sillh)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillh, 'long_name', 'gate crest level (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillh, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillh, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillh, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_crest_width', nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_sillw)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillw, 'long_name', 'gate crest width (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillw, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_sillw, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_gate_lower_edge_level', nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_edgel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_edgel, 'long_name', 'gate lower edge level (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_edgel, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_edgel, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_edgel, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_flow_through_height', nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_flowh)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_flowh, 'standard_name', 'gategen_flow_through_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_flowh, 'long_name', 'gate flow through height (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_flowh, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_flowh, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_flowh, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_gate_opening_width', nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_openw)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_openw, 'long_name', 'gate opening width (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_openw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_openw, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_openw, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_s1up',     nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1up, 'long_name', 'gate water level up (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1up, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1up, 'geometry', gategen_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'gategen_s1dn',     nf90_double, (/ hisids%id_gategendim, hisids%id_timedim /), hisids%id_gategen_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1dn, 'long_name', 'gate water level down (via general structure)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1dn, 'coordinates', 'gategen_name')
            ierr = nf90_put_att(hisids%ncid, hisids%id_gategen_s1dn, 'geometry', gategen_geom_container_name)
        endif

        if(jahiscdam > 0 .and. ncdamsg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'dams', ncdamsg, hisids%id_cdamdim)
            !ierr = nf90_def_dim(hisids%ncid, 'cdam_name_len', strlen_netcdf, hisids%id_cdamlendim)
            ierr = nf90_def_var(hisids%ncid, 'cdam_name',  nf90_char,   (/ hisids%id_strlendim, hisids%id_cdamdim /), hisids%id_cdamname)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdamname,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdamname,  'long_name', 'controllable dam name'    )

            ierr = nf90_def_var(hisids%ncid, 'cdam_discharge',     nf90_double, (/ hisids%id_cdamdim, hisids%id_timedim /), hisids%id_cdam_dis)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_dis, 'standard_name', 'cdam_discharge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_dis, 'long_name', 'controllable dam discharge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_dis, 'coordinates', 'cdam_name')

            ierr = nf90_def_var(hisids%ncid, 'cdam_crest_level',     nf90_double, (/ hisids%id_cdamdim, hisids%id_timedim /), hisids%id_cdam_crestl)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_crestl, 'standard_name', 'cdam_crest_level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_crestl, 'long_name', 'controllable dam crest level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_crestl, 'coordinates', 'cdam_name')

            ierr = nf90_def_var(hisids%ncid, 'cdam_s1up',     nf90_double, (/ hisids%id_cdamdim, hisids%id_timedim /), hisids%id_cdam_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1up, 'long_name', 'controllable dam water level up')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1up, 'coordinates', 'cdam_name')

            ierr = nf90_def_var(hisids%ncid, 'cdam_s1dn',     nf90_double, (/ hisids%id_cdamdim, hisids%id_timedim /), hisids%id_cdam_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1dn, 'long_name', 'controllable dam water level down')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cdam_s1dn, 'coordinates', 'cdam_name')
        endif

        if(jahisweir > 0 .and. nweirgen > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'weirgens', nweirgen, hisids%id_weirgendim)
            ierr = nf90_def_var(hisids%ncid, 'weirgen_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_weirgendim /), hisids%id_weirgen_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_id,  'long_name', 'Id of weir'    )

            ! Define geometry related variables
            weir_geom_container_name = 'weirgen_geom'
            nNodeTot = 0
            if (network%sts%numWeirs > 0) then ! new weir
               nNodeTot = nNodesWeir
            else ! old weir
               do n = 1, nweirgen
                  i = weir2cgen(n)
                  nlinks = L2cgensg(i) - L1cgensg(i) + 1
                  if (nlinks > 0) then
                     nNodes = nlinks + 1
                  else if (nlinks == 0) then
                     nNodes = 0
                  end if
                  nNodeTot = nNodeTot + nNodes
               end do
            end if

            ierr = sgeom_def_geometry_variables(hisids%ncid, weir_geom_container_name, 'weir', 'line', nNodeTot, hisids%id_weirgendim, &
               hisids%id_weirgeom_node_count, hisids%id_weirgeom_node_coordx, hisids%id_weirgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'weirgen_discharge',     nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_dis)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_dis, 'standard_name', 'integral_of_discharge_wrt_time') ! TODO: introduce time windows in nc
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_dis, 'long_name', 'Discharge through weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_dis, 'coordinates', 'weirgen_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_dis, 'geometry', weir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'weirgen_crest_level', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_crestl)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestl, 'standard_name', 'weirgen_crest_level')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestl, 'long_name', 'Crest level of weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestl, 'coordinates', 'weirgen_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestl, 'geometry', weir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'weirgen_crest_width', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_crestw)
            !ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestw, 'standard_name', 'weirgen_crest_width')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestw, 'long_name', 'Crest width of weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestw, 'coordinates', 'weirgen_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_crestw, 'geometry', weir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'weirgen_s1up',     nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1up, 'long_name', 'Water level upstream of weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1up, 'coordinates', 'weirgen_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1up, 'geometry', weir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'weirgen_s1dn',     nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1dn, 'long_name', 'Water level downstream of weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1dn, 'coordinates', 'weirgen_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1dn, 'geometry', weir_geom_container_name)

            if (network%sts%numWeirs > 0) then ! write extra files for new weirs
               ierr = nf90_def_var(hisids%ncid, 'weirgen_structure_head', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_head)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_head, 'long_name', 'Head difference across weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_head, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_head, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_head, 'geometry', weir_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'weirgen_velocity', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_vel)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_vel, 'long_name', 'Velocity through weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_vel, 'units', 'm s-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_vel, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_vel, 'geometry', weir_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'weirgen_flow_area', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_au)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_au, 'long_name', 'Flow area at weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_au, 'units', 'm2')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_au, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_au, 'geometry', weir_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'weirgen_state', nf90_int, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weir_stat)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'long_name', 'Flow state at weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'units', '-')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'flag_values', (/ 0, 1, 2 /))
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'flag_meanings', 'no_flow weir_free weir_submerged')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'valid_range', (/ 0, 2 /))
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, '_FillValue', int(dmiss))
               ierr = nf90_put_att(hisids%ncid, hisids%id_weir_stat, 'geometry', weir_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'weirgen_force_difference', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_forcedif)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_forcedif, 'long_name', 'Force difference per unit width at weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_forcedif, 'units', 'N m-1')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_forcedif, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_forcedif, 'geometry', weir_geom_container_name)

               ierr = nf90_def_var(hisids%ncid, 'weirgen_s1_on_crest', nf90_double, (/ hisids%id_weirgendim, hisids%id_timedim /), hisids%id_weirgen_s1crest)
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1crest, 'long_name', 'Water level on crest of weir')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1crest, 'units', 'm')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1crest, 'coordinates', 'weirgen_id')
               ierr = nf90_put_att(hisids%ncid, hisids%id_weirgen_s1crest, 'geometry', weir_geom_container_name)
            end if
        endif

        ! Orifice
        if(jahisorif > 0 .and. network%sts%numOrifices > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'orifice', network%sts%numOrifices, hisids%id_orifgendim)
            ierr = nf90_def_var(hisids%ncid, 'orifice_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_orifgendim /), hisids%id_orifgen_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_id,  'long_name', 'Id of orifice')

            ! Define geometry related variables
            orif_geom_container_name = 'orifice_geom'
            nNodeTot = nNodesOrif

            ierr = sgeom_def_geometry_variables(hisids%ncid, orif_geom_container_name, 'orifice', 'line', nNodeTot, hisids%id_orifgendim, &
               hisids%id_orifgeom_node_count, hisids%id_orifgeom_node_coordx, hisids%id_orifgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'orifice_discharge',     nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_dis, 'long_name', 'Discharge through orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_dis, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_dis, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_crest_level', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_crestl)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestl, 'long_name', 'Crest level of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestl, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestl, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_crest_width', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_crestw)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestw, 'long_name', 'Crest width of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestw, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_crestw, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_gate_lower_edge_level', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_edgel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_edgel, 'long_name', 'Gate lower edge level of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_edgel, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_edgel, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_edgel, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_s1up',     nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1up, 'long_name', 'Water level upstream of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1up, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1up, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_s1dn',     nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1dn, 'long_name', 'Water level downstream of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1dn, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1dn, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_gate_opening_height', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_openh)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_openh, 'long_name', 'Gate opening height of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_openh, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_openh, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_openh, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_head', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_head, 'long_name', 'Head difference across orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_head, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_head, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_flow_area', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_au, 'long_name', 'Flow area at orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_au, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_au, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_state', nf90_int, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_stat)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'long_name', 'Flow state at orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'units', '-')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'units', '-')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'flag_values', (/ 0, 1, 2, 3, 4 /))
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'flag_meanings', 'no_flow weir_free weir_submerged gate_free gate_submerged')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'valid_range', (/ 0, 4 /))
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, '_FillValue', int(dmiss))
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_stat, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_s1_on_crest', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_s1crest)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1crest, 'long_name', 'Water level on crest of orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1crest, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1crest, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_s1crest, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_velocity', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_vel, 'long_name', 'Velocity through orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_vel, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_vel, 'geometry', orif_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'orifice_force_difference', nf90_double, (/ hisids%id_orifgendim, hisids%id_timedim /), hisids%id_orifgen_forcedif)
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_forcedif, 'long_name', 'Force difference per unit width at orifice')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_forcedif, 'units', 'N m-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_forcedif, 'coordinates', 'orifice_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_orifgen_forcedif, 'geometry', orif_geom_container_name)
        endif

        ! Bridge
        if(jahisbridge > 0 .and. network%sts%numBridges > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'bridge', network%sts%numbridges, hisids%id_bridgedim)
            ierr = nf90_def_var(hisids%ncid, 'bridge_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_bridgedim /), hisids%id_bridge_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_id,  'long_name', 'Id of bridge')

            ! Define geometry related variables
            bridge_geom_container_name = 'bridge_geom'
            nNodeTot = nNodesBridge

            ierr = sgeom_def_geometry_variables(hisids%ncid, bridge_geom_container_name, 'bridge', 'line', nNodeTot, hisids%id_bridgedim, &
               hisids%id_bridgegeom_node_count, hisids%id_bridgegeom_node_coordx, hisids%id_bridgegeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'bridge_discharge',     nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_dis, 'long_name', 'Discharge through bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_dis, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_dis, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_s1up',     nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1up, 'long_name', 'Water level upstream of bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1up, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1up, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_s1dn',     nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1dn, 'long_name', 'Water level downstream of bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1dn, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_s1dn, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_head', nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_head, 'long_name', 'Head difference across bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_head, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_head, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_flow_area', nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_au, 'long_name', 'Flow area at bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_au, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_au, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_velocity', nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_vel, 'long_name', 'Velocity through bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_vel, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_vel, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_blup',  nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_blup)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_blup, 'standard_name', 'altitude')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_blup, 'long_name', 'Bed level at upstream of bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_blup, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_blup, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_blup, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_bldn',  nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_bldn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bldn, 'standard_name', 'altitude')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bldn, 'long_name', 'Bed level at downstream of bridge')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bldn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bldn, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bldn, 'geometry', bridge_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'bridge_bl_actual',  nf90_double, (/ hisids%id_bridgedim, hisids%id_timedim /), hisids%id_bridge_bl_act)
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bl_act, 'standard_name', 'altitude')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bl_act, 'long_name', 'Actual bed level of bridge (crest)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bl_act, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bl_act, 'coordinates', 'bridge_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_bridge_bl_act, 'geometry', bridge_geom_container_name)
        endif

        ! Culvert
        if(jahisculv > 0 .and. network%sts%numculverts > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'culvert', network%sts%numculverts, hisids%id_culvertdim)
            ierr = nf90_def_var(hisids%ncid, 'culvert_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_culvertdim /), hisids%id_culvert_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_id,  'long_name', 'Id of culvert'    )

            ! Define geometry related variables
            culvert_geom_container_name = 'culvert_geom'
            nNodeTot = nNodesCulv

            ierr = sgeom_def_geometry_variables(hisids%ncid, culvert_geom_container_name, 'culvert', 'line', nNodeTot, hisids%id_culvertdim, &
               hisids%id_culvertgeom_node_count, hisids%id_culvertgeom_node_coordx, hisids%id_culvertgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'culvert_discharge',     nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_dis, 'long_name', 'Discharge through culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_dis, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_dis, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_crest_level', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_crestl)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_crestl, 'long_name', 'Crest level of culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_crestl, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_crestl, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_gate_lower_edge_level', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_edgel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_edgel, 'long_name', 'Gate lower edge level of culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_edgel, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_edgel, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_edgel, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_s1up',     nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1up, 'long_name', 'Water level upstream of culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1up, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1up, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_s1dn',     nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1dn, 'long_name', 'Water level downstream of culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1dn, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_s1dn, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_gate_opening_height', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_openh)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_openh, 'long_name', 'Gate opening height of culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_openh, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_openh, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_openh, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_head', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_head, 'long_name', 'Head difference across culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_head, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_head, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_flow_area', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_au, 'long_name', 'Flow area in culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_au, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_au, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_state', nf90_int, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_stat)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'long_name', 'Flow state in culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'units', '-')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'flag_values', (/ 0, 1, 2/))
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'flag_meanings', 'no_flow culvert_free culvert_submerged')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'valid_range', (/ 0, 2 /))
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, '_FillValue', int(dmiss))
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_stat, 'geometry', culvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'culvert_velocity', nf90_double, (/ hisids%id_culvertdim, hisids%id_timedim /), hisids%id_culvert_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_vel, 'long_name', 'Velocity in culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_vel, 'coordinates', 'culvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_culvert_vel, 'geometry', culvert_geom_container_name)
        endif


        ! Dambreak
        if (jahisdambreak > 0 .and. ndambreaksg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'dambreaks', ndambreaksg, hisids%id_dambreakdim)
            ierr = nf90_def_var(hisids%ncid, 'dambreak_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_dambreakdim /), hisids%id_dambreak_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_id,  'long_name', 'Id of dambreak')

            ierr = nf90_def_var(hisids%ncid, 'dambreak_s1up', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1up, 'long_name', 'Water level upstream of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1up, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1up, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_s1dn',     nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1dn, 'long_name', 'Water level downstream of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1dn, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_s1dn, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_discharge', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_discharge)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_discharge, 'long_name', 'Instantaneous discharge through dambreaks')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_discharge, 'units', 'm3 s-1') !link_sum
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_discharge, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_discharge, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_cumulative_discharge', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_cumulative_discharge)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cumulative_discharge, 'long_name', 'Cumulative discharge through dambreaks')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cumulative_discharge, 'units', 'm3') !link_sum
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cumulative_discharge, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cumulative_discharge, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_breach_width_time_derivative', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_breach_width_time_derivative)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_breach_width_time_derivative, 'long_name', 'Breach width time derivative of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_breach_width_time_derivative, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_breach_width_time_derivative, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_breach_width_time_derivative, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_water_level_jump', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_water_level_jump)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_water_level_jump, 'long_name', 'Breach water level jump of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_water_level_jump, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_water_level_jump, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_water_level_jump, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_normal_velocity', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_normal_velocity)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_normal_velocity, 'long_name', 'Normal velocity through dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_normal_velocity, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_normal_velocity, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_normal_velocity, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_structure_head', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_head, 'long_name', 'Head difference across dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_head, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_head, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_flow_area', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_au, 'long_name', 'Flow area at dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_au, 'coordinates', 'dambreak_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_au, '_FillValue', dmiss)

            ierr = nf90_def_var(hisids%ncid, 'dambreak_crest_level', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_cresth)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cresth, 'long_name', 'Crest level of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cresth, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_cresth, 'coordinates', 'dambreak_id')

            ierr = nf90_def_var(hisids%ncid, 'dambreak_crest_width', nf90_double, (/ hisids%id_dambreakdim, hisids%id_timedim /), hisids%id_dambreak_crestw)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_crestw, 'long_name', 'Crest width of dambreak')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_crestw, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dambreak_crestw, 'coordinates', 'dambreak_id')
        endif

        ! Universal weir
        if(jahisuniweir > 0 .and. network%sts%numuniweirs > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'universalWeirs', network%sts%numuniweirs, hisids%id_uniweirdim)
            ierr = nf90_def_var(hisids%ncid, 'uniweir_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_uniweirdim /), hisids%id_uniweir_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_id,  'long_name', 'Id of universal weir'    )

            ! Define geometry related variables
            uniweir_geom_container_name = 'uniweir_geom'
            nNodeTot = nNodesUniweir

            ierr = sgeom_def_geometry_variables(hisids%ncid, uniweir_geom_container_name, 'uniweir', 'line', nNodeTot, hisids%id_uniweirdim, &
               hisids%id_uniweirgeom_node_count, hisids%id_uniweirgeom_node_coordx, hisids%id_uniweirgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_discharge',     nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_dis, 'long_name', 'Discharge through universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_dis, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_dis, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_crest_level', nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_crestl)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_crestl, 'long_name', 'Crest level of universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_crestl, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_crestl, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_crestl, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_s1up',     nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1up, 'long_name', 'Water level upstream of universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1up, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1up, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_s1dn',     nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1dn, 'long_name', 'Water level downstream of universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1dn, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_s1dn, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_head', nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_head, 'long_name', 'Head difference across universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_head, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_head, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_flow_area', nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_au, 'long_name', 'Flow area in universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_au, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_au, 'geometry', uniweir_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'uniweir_velocity', nf90_double, (/ hisids%id_uniweirdim, hisids%id_timedim /), hisids%id_uniweir_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_vel, 'long_name', 'Velocity through universal weir')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_vel, 'coordinates', 'uniweir_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_uniweir_vel, 'geometry', uniweir_geom_container_name)
        endif

        ! compound structure
        if(jahiscmpstru > 0 .and. network%cmps%count > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'compoundStructures', network%cmps%count, hisids%id_cmpstrudim)
            ierr = nf90_def_var(hisids%ncid, 'cmpstru_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_cmpstrudim /), hisids%id_cmpstru_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_id,  'long_name', 'Id of compound structure')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_discharge',     nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_dis, 'long_name', 'Discharge through compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_dis, 'coordinates', 'cmpstru_id')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_s1up',     nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1up, 'long_name', 'Water level upstream of compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1up, 'coordinates', 'cmpstru_id')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_s1dn',     nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1dn, 'long_name', 'Water level downstream of compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_s1dn, 'coordinates', 'cmpstru_id')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_head', nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_head, 'long_name', 'Head difference across compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_head, 'coordinates', 'cmpstru_id')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_flow_area', nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_au, 'long_name', 'Flow area in compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_au, 'coordinates', 'cmpstru_id')

            ierr = nf90_def_var(hisids%ncid, 'cmpstru_velocity', nf90_double, (/ hisids%id_cmpstrudim, hisids%id_timedim /), hisids%id_cmpstru_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_vel, 'long_name', 'Velocity through compound structure')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_cmpstru_vel, 'coordinates', 'cmpstru_id')
        endif

        ! Long culvert
        if(jahislongculv > 0 .and. nlongculvertsg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'longculvert', nlongculvertsg, hisids%id_longculvertdim)
            ierr = nf90_def_var(hisids%ncid, 'longculvert_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_longculvertdim /), hisids%id_longculvert_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_id,  'long_name', 'Id of long culvert')

            ! Define geometry related variables
            longculvert_geom_container_name = 'longculvert_geom'
            nNodeTot = 0
            nNodeTot = nNodesLongCulv

            ierr = sgeom_def_geometry_variables(hisids%ncid, longculvert_geom_container_name, 'longculvert', 'line', nNodeTot, hisids%id_longculvertdim, &
               hisids%id_longculvertgeom_node_count, hisids%id_longculvertgeom_node_coordx, hisids%id_longculvertgeom_node_coordy)
            !ierr = ncu_clone_vardef(hisids%ncid, hisids%ncid, hisids%id_longculvertgeom_node_coordx, 'longculvert_geom_lon', hisids%id_laydim, & !hisids%id_statlon, &
            !              'longitude', 'blabla1', 'degrees_east')

            ierr = nf90_def_var(hisids%ncid, 'longculvert_discharge', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_dis, 'long_name', 'Discharge through long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_dis, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_dis, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_dis, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_s1up', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_s1up)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1up, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1up, 'long_name', 'Water level upstream of long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1up, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1up, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1up, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_s1dn', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_s1dn)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1dn, 'standard_name', 'sea_surface_height')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1dn, 'long_name', 'Water level downstream of long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1dn, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1dn, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_s1dn, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_head', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_head)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_head, 'long_name', 'Head difference across long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_head, 'units', 'm')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_head, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_head, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_flow_area', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_au)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_au, 'long_name', 'Flow area in long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_au, 'units', 'm2')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_au, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_au, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_velocity', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_vel)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_vel, 'long_name', 'Velocity in long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_vel, 'units', 'm s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_vel, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_vel, 'geometry', longculvert_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'longculvert_valve_relative_opening', nf90_double, (/ hisids%id_longculvertdim, hisids%id_timedim /), hisids%id_longculvert_valveopen)
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_valveopen, 'long_name', 'Valve relative opening in long culvert')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_valveopen, 'units', '1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_valveopen, 'coordinates', 'longculvert_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_longculvert_valveopen, 'geometry', longculvert_geom_container_name)
        endif

        ! Lateral
        if(jahislateral > 0 .and. numlatsg > 0 .and. janormalhis > 0) then
            ierr = nf90_def_dim(hisids%ncid, 'lateral', numlatsg, hisids%id_latdim)
            ierr = nf90_def_var(hisids%ncid, 'lateral_id',  nf90_char,   (/ hisids%id_strlendim, hisids%id_latdim /), hisids%id_lat_id)
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_id,  'cf_role',   'timeseries_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_id,  'long_name', 'Id of lateral')

            ! Define geometry related variables
            lat_geom_container_name = 'lateral_geom'
            nNodeTot = nNodesLat
            ierr = sgeom_def_geometry_variables(hisids%ncid, lat_geom_container_name, 'lateral', 'point', nNodeTot, hisids%id_latdim, &
               hisids%id_latgeom_node_count, hisids%id_latgeom_node_coordx, hisids%id_latgeom_node_coordy)

            ierr = nf90_def_var(hisids%ncid, 'lateral_prescribed_discharge_instantaneous', nf90_double, (/ hisids%id_latdim, hisids%id_timedim /), hisids%id_lat_predis_inst)
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_inst, 'long_name', 'Prescribed discharge through lateral at current time step (instantaneous)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_inst, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_inst, 'coordinates', 'lateral_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_inst, 'geometry', lat_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'lateral_prescribed_discharge_average', nf90_double, (/ hisids%id_latdim, hisids%id_timedim /), hisids%id_lat_predis_ave)
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_ave, 'long_name', 'Prescribed discharge through lateral, average over the last history time interval')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_ave, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_ave, 'coordinates', 'lateral_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_predis_ave, 'geometry', lat_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'lateral_realized_discharge_instantaneous', nf90_double, (/ hisids%id_latdim, hisids%id_timedim /), hisids%id_lat_realdis_inst)
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_inst, 'long_name', 'Realized discharge through lateral at current time step (instantaneous)')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_inst, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_inst, 'coordinates', 'lateral_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_inst, 'geometry', lat_geom_container_name)

            ierr = nf90_def_var(hisids%ncid, 'lateral_realized_discharge_average', nf90_double, (/ hisids%id_latdim, hisids%id_timedim /), hisids%id_lat_realdis_ave)
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_ave, 'long_name', 'Realized discharge through lateral, average over the last history time interval')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_ave, 'units', 'm3 s-1')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_ave, 'coordinates', 'lateral_id')
            ierr = nf90_put_att(hisids%ncid, hisids%id_lat_realdis_ave, 'geometry', lat_geom_container_name)
        endif
        if (timon) call timstop (handle_extra(60))

        if(dad_included .and. janormalhis > 0) then  ! Output for dredging and dumping
            ierr = nf90_def_dim(hisids%ncid, 'ndredlink', dadpar%nalink, hisids%id_dredlinkdim)
            ierr = nf90_def_dim(hisids%ncid, 'ndred', dadpar%nadred+dadpar%nasupl, hisids%id_dreddim)
            ierr = nf90_def_dim(hisids%ncid, 'ndump', dadpar%nadump, hisids%id_dumpdim)

            ierr = nf90_def_var(hisids%ncid, 'dredge_area_name',         nf90_char,   (/ hisids%id_strlendim, hisids%id_dreddim /), hisids%id_dred_name)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dred_name,  'long_name'    , 'dredge area identifier')

            ierr = nf90_def_var(hisids%ncid, 'dump_area_name',         nf90_char,   (/ hisids%id_strlendim, hisids%id_dumpdim /), hisids%id_dump_name)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dump_name,  'long_name'    , 'dump area identifier')

            ierr = nf90_def_var(hisids%ncid, 'dred_link_discharge',     nf90_double, (/ hisids%id_dredlinkdim, hisids%id_sedtotdim, hisids%id_timedim /), hisids%id_dredlink_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dredlink_dis, 'long_name', 'Cumulative dredged material transported via links per fraction')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dredlink_dis, 'units', 'm3') !link_sum

            ierr = nf90_def_var(hisids%ncid, 'dred_discharge',     nf90_double, (/ hisids%id_dreddim, hisids%id_timedim /), hisids%id_dred_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dred_dis, 'long_name', 'Cumulative dredged material for dredge areas')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dred_dis, 'units', 'm3') !totvoldred

            ierr = nf90_def_var(hisids%ncid, 'dump_discharge',     nf90_double, (/ hisids%id_dumpdim, hisids%id_timedim /), hisids%id_dump_dis)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dump_dis, 'long_name', 'Cumulative dredged material for dump areas')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dump_dis, 'units', 'm3') !totvoldump

            ierr = nf90_def_var(hisids%ncid, 'dred_time_frac',     nf90_double, (/ hisids%id_dreddim, hisids%id_timedim /), hisids%id_dred_tfrac)
            ierr = nf90_put_att(hisids%ncid, hisids%id_dred_tfrac, 'long_name', 'Time fraction spent dredging')
            ierr = nf90_put_att(hisids%ncid, hisids%id_dred_tfrac, 'units', '-') !ndredged

            ierr = nf90_def_var(hisids%ncid, 'plough_time_frac',   nf90_double, (/ hisids%id_dreddim, hisids%id_timedim /), hisids%id_plough_tfrac)
            ierr = nf90_put_att(hisids%ncid, hisids%id_plough_tfrac, 'long_name', 'Time fraction spent ploughing')
            ierr = nf90_put_att(hisids%ncid, hisids%id_plough_tfrac, 'units', '-') !nploughed
        endif

        if ( jacheckmonitor.eq.1 .and. janormalhis > 0) then
           ierr = nf90_def_var(hisids%ncid, 'checkerboard_monitor', nf90_double, (/ hisids%id_laydim, hisids%id_timedim /), hisids%id_checkmon)
           ierr = nf90_put_att(hisids%ncid, hisids%id_checkmon, 'long_name', 'Checkerboard mode monitor')
           ierr = nf90_put_att(hisids%ncid, hisids%id_checkmon, 'unit', 'm s-1')

           ierr = nf90_def_var(hisids%ncid, 'num_timesteps', nf90_int, hisids%id_timedim, hisids%id_num_timesteps)
           ierr = nf90_def_var(hisids%ncid, 'comp_time', nf90_double, hisids%id_timedim, hisids%id_comp_time)
        end if

        ierr = nf90_def_var(hisids%ncid, 'time', nf90_double, hisids%id_timedim, hisids%id_time)
        ierr = nf90_put_att(hisids%ncid, hisids%id_time,  'units'        , trim(Tudunitstr))
        ierr = nf90_put_att(hisids%ncid, hisids%id_time,  'standard_name', 'time')

        ! Size of latest timestep
        ierr = unc_def_var_nonspatial(hisids%ncid, hisids%id_timestep, nf90_double, (/ hisids%id_timedim /), 'timestep', '',     'latest computational timestep size in each output interval', 's')

        if ( japart.gt.0 .and. janormalhis > 0) then
!          write partiles header to hisfile
           call unc_write_part_header(hisids%ncid,hisids%id_timedim,hisids%id_partdim,hisids%id_parttime,hisids%id_partx,hisids%id_party,hisids%id_partz)
        end if

        ierr = nf90_enddef(hisids%ncid)
        if (timon) call timstop (handle_extra(61))

        if (timon) call timstrt ('unc_write_his timeindep data', handle_extra(63))
        do i=1,numobs+nummovobs
!           ierr = nf90_put_var(hisids%ncid, hisids%id_statx,    xobs(i),         (/ i /))
!           ierr = nf90_put_var(hisids%ncid, hisids%id_staty,    yobs(i),         (/ i /))
           ierr = nf90_put_var(hisids%ncid, hisids%id_statid, trimexact(namobs(i), strlen_netcdf), (/ 1, i /)) ! Extra for OpenDA-wrapper
           ierr = nf90_put_var(hisids%ncid, hisids%id_statname, trimexact(namobs(i), strlen_netcdf), (/ 1, i /))
        end do

        if (ncrs > 0 .and. janormalhis > 0) then
            do i=1,ncrs
                !ierr = nf90_put_var(hisids%ncid, hisids%id_crsx,     crs(i)%path%xp(1:crs(i)%path%np), (/ 1, i /))
                !ierr = nf90_put_var(hisids%ncid, hisids%id_crsy,     crs(i)%path%yp(1:crs(i)%path%np), (/ 1, i /))
                ierr = nf90_put_var(hisids%ncid, hisids%id_crsname,  trimexact(crs(i)%name, strlen_netcdf),      (/ 1, i /))
            end do
            if (hisids%idx_curtime == 0) then
               ierr = nf90_put_var(hisids%ncid, hisids%id_crsgeom_node_coordx, geomXCrs,     start = (/ 1 /), count = (/ nNodesCrs /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_crsgeom_node_coordy, geomYCrs,     start = (/ 1 /), count = (/ nNodesCrs /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_crsgeom_node_count,  nodeCountCrs)
               if (allocated(geomXCrs))     deallocate(geomXCrs)
               if (allocated(geomYCrs))     deallocate(geomYCrs)
               if (allocated(nodeCountCrs)) deallocate(nodeCountCrs)
               end if
                  end if

        if (nrug>0 .and. janormalhis > 0) then
            do i=1,nrug
                ierr = nf90_put_var(hisids%ncid, hisids%id_rugname,  trimexact(rug(i)%name, strlen_netcdf), (/ 1, i /))
                ierr = nf90_put_var(hisids%ncid, hisids%id_rugid,    trimexact(rug(i)%name, strlen_netcdf), (/ 1, i /))
            end do
        endif

        if (jahiscgen > 0 .and. ntmp > 0 .and. janormalhis > 0) then
            do i=1,ntmp
               if (jaoldstr == 1) then
                  igen = i
               else
                  if (network%sts%numGeneralStructures > 0) then
                     istru = network%sts%generalStructureIndices(i)
                     ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),  (/ 1, i /))
                     cycle
                  else
                     igen = genstru2cgen(i)
                  end if
               end if

               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_id,  trimexact(cgen_ids(igen), strlen_netcdf), (/ 1, i /))
            end do
        end if

        if (jahisorif > 0 .and. network%sts%numOrifices > 0 .and. janormalhis > 0) then
           do i = 1, network%sts%numOrifices
              istru = network%sts%orificeIndices(i)
              ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jahisbridge > 0 .and. network%sts%numBridges > 0 .and. janormalhis > 0) then
           do i = 1, network%sts%numBridges
              istru = network%sts%bridgeIndices(i)
              ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jahisculv > 0 .and. network%sts%numCulverts > 0 .and. janormalhis > 0) then
           do i = 1, network%sts%numCulverts
              istru = network%sts%culvertIndices(i)
              ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jahisuniweir > 0 .and. network%sts%numuniweirs > 0 .and. janormalhis > 0) then
           do i = 1, network%sts%numuniweirs
              istru = network%sts%uniweirIndices(i)
              ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jahiscmpstru > 0 .and. network%cmps%count > 0 .and. janormalhis > 0) then
           do i = 1, network%cmps%count
              ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_id,  trimexact(network%cmps%compound(i)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jahislateral > 0 .and. numlatsg > 0 .and. janormalhis > 0) then
           do i = 1, numlatsg
              ierr = nf90_put_var(hisids%ncid, hisids%id_lat_id,  trimexact(lat_ids(i), strlen_netcdf), (/ 1, i /))
           end do
        end if

        if (jahispump > 0 .and. npumpsg > 0 .and. janormalhis > 0) then
            do i=1,npumpsg
               ierr = nf90_put_var(hisids%ncid, hisids%id_pump_id,  trimexact(pump_ids(i), strlen_netcdf),      (/ 1, i /))

               ! NOTE: the code below is now only active for pumps (DELFT3D-36341). Should be
               ! generalized for all structure locations that are polyline based.
               !
               ! Store one single representative x/y point for each pump in the time series file,
               ! because CF conventions require that for variables on discrete geometries.
               ! Computed at half the total length of the snapped flow links
               ! (so, it lies on an edge, not per se on the input pump pli)).
               pumplensum = 0d0
               do k = L1pumpsg(i), L2pumpsg(i)
                  Lf = abs(kpump(3,k))
                  pumplensum = pumplensum + wu(Lf)
               end do
               pumplenmid = pumplensum / 2

               ! Find the mid point on the snapped flow link path
               pumplensum = 0d0
               do k = L1pumpsg(i), L2pumpsg(i)
                  Lf = abs(kpump(3,k))
                  if (pumplensum + wu(Lf) >= pumplenmid) then
                     if (kcu(Lf) == 2) then ! 2D
                        if (kpump(3,k) > 0) then
                           k3 = lncn(1,Lf)
                           k4 = lncn(2,Lf)
                        else
                           k3 = lncn(2,Lf)
                           k4 = lncn(1,Lf)
                        end if
                        w1 = (pumplenmid-pumplensum)/wu(Lf)
                        pumpxmid = w1*xk(k3) + (1d0-w1)*xk(k4)
                        pumpymid = w1*yk(k3) + (1d0-w1)*yk(k4)
                     else                   ! 1D
                        pumpxmid = xu(Lf)
                        pumpymid = yu(Lf)
                     end if
                     exit ! mid point was found
                  else
                     pumplensum = pumplensum + wu(Lf)
                  end if
               end do

               ierr = nf90_put_var(hisids%ncid, hisids%id_pump_xmid,  pumpxmid,      (/ i /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_pump_ymid,  pumpymid,      (/ i /))
            end do
        end if
        if (jahisgate > 0 .and. ngatesg > 0 .and. janormalhis > 0) then
            do i=1,ngatesg
               ierr = nf90_put_var(hisids%ncid, hisids%id_gatename,  trimexact(gate_ids(i), strlen_netcdf),      (/ 1, i /))
            end do
        end if
        if (jahisgate > 0 .and. ngategen > 0 .and. janormalhis > 0) then
           do i=1,ngategen
              igen = gate2cgen(i)
              ierr = nf90_put_var(hisids%ncid, hisids%id_gategenname,  trimexact(cgen_ids(igen), strlen_netcdf),      (/ 1, i /))
           end do
        end if
        if (jahiscdam > 0 .and. ncdamsg > 0 .and. janormalhis > 0) then
            do i=1,ncdamsg
               ierr = nf90_put_var(hisids%ncid, hisids%id_cdamname,  trimexact(cdam_ids(i), strlen_netcdf),      (/ 1, i /))
            end do
        end if
        if (jahisweir > 0 .and. nweirgen > 0 .and. janormalhis > 0) then
           if (allocated(weir2cgen)) then
              do i=1,nweirgen
                 igen = weir2cgen(i)
                 ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_id,  trimexact(cgen_ids(igen), strlen_netcdf),      (/ 1, i /))
              end do
           else if (network%sts%numWeirs > 0) then
              do i=1,nweirgen
                 istru = network%sts%weirIndices(i)
                 ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_id,  trimexact(network%sts%struct(istru)%id, strlen_netcdf),      (/ 1, i /))
              end do
           end if
        end if

        if (jahisdambreak > 0 .and. ndambreak > 0 .and. janormalhis > 0) then
            do i = 1,ndambreaksg
               ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_id, trimexact(dambreak_ids(i), strlen_netcdf),(/ 1, i /))
            end do
        end if

        if (jahislongculv > 0 .and. nlongculvertsg > 0 .and. janormalhis > 0) then
           do i = 1, nlongculvertsg
              ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_id,  trimexact(longculverts(i)%id, strlen_netcdf),  (/ 1, i /))
           end do
        end if

        if (jased>0 .and. stm_included .and. jahissed>0) then
           do i=1,stmpar%lsedtot
              ierr = nf90_put_var(hisids%ncid, hisids%id_frac_name, trimexact(stmpar%sedpar%namsed(i), strlen_netcdf), (/ 1, i /))
           enddo
        end if

        if (dad_included .and. janormalhis > 0) then
           !
           !do i=1,stmpar%lsedtot
           !   ierr = nf90_put_var(hisids%ncid, hisids%id_frac_name, trimexact(stmpar%sedpar%namsed(i), strlen_netcdf), (/ 1, i /))
           !enddo
           !ierr = nf90_put_var(hisids%ncid, hisids%id_frac_name, 'subsoil sediment', (/ 1, stmpar%lsedtot+1 /))        ! rest category
           !
           do i=1,(dadpar%nadred+dadpar%nasupl)
              ierr = nf90_put_var(hisids%ncid, hisids%id_dred_name, trimexact(dadpar%dredge_areas(i), strlen_netcdf), (/ 1, i /))
           enddo
           !
           do i=1,dadpar%nadump
              ierr = nf90_put_var(hisids%ncid, hisids%id_dump_name, trimexact(dadpar%dump_areas(i), strlen_netcdf), (/ 1, i /))
           enddo
        endif
        if (timon) call timstop ( handle_extra(63))
    endif
    ! Increment output counters in m_flowtimes.
    !time_his = tim
    hisids%idx_curtime   = hisids%idx_curtime + 1

    if (timon) call timstrt ('unc_write_his time data', handle_extra(64))

    ierr = nf90_put_var(hisids%ncid, hisids%id_time, tim, (/ hisids%idx_curtime /))
    ierr = nf90_put_var(hisids%ncid, hisids%id_timestep, dts, (/ hisids%idx_curtime /))
    if (timon) call timstop ( handle_extra(64))

!   write particles to hisfile (for now)
    if ( japart.gt.0 .and. janormalhis > 0) then
       call unc_write_part(hisids%ncid,hisids%idx_curtime,hisids%id_parttime,hisids%id_partx,hisids%id_party,hisids%id_partz)
    end if

!   Observation points (fixed+moving)

    ntot = numobs + nummovobs
    valobsT(1:ntot, 1:IPNT_NUM) = transpose(valobs)
    if (ntot > 0) then
       if (janormalhis > 0) then
          if (jahiswatlev > 0) then
             ierr = nf90_put_var(hisids%ncid,    hisids%id_vars,   valobsT(:,IPNT_S1),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          if (jahiswatdep > 0) then
             ierr = nf90_put_var(hisids%ncid,    hisids%id_hs  ,   valobsT(:,IPNT_HS),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          if( stm_included ) then
             ierr = nf90_put_var(hisids%ncid,    hisids%id_varb,   valobsT(:,IPNT_BL),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          else
             ierr = nf90_put_var(hisids%ncid,    hisids%id_varb,   valobsT(:,IPNT_BL),    start = (/ 1 /) )
          endif

       end if
       ! write geometry variables at the first time of history output
       if (hisids%idx_curtime == 1) then
          call realloc(node_count, numobs)
          node_count = 1
          ierr = nf90_put_var(hisids%ncid,    hisids%id_statgeom_node_count, node_count)
          ierr = nf90_put_var(hisids%ncid,    hisids%id_statgeom_node_coordx,  xobs(:), start = (/ 1 /), count = (/ numobs /))
          ierr = nf90_put_var(hisids%ncid,    hisids%id_statgeom_node_coordy,  yobs(:), start = (/ 1 /), count = (/ numobs /))
#ifdef HAVE_PROJ
          if (add_latlon) then
!             call transform_and_put_latlon_coordinates(hisids%ncid, hisids%id_statgeom_node_lon, hisids%id_statgeom_node_lat, nccrs%proj_string, xobs, yobs)
          end if
#endif
       end if

       if ( nummovobs > 0 ) then
          ierr = nf90_put_var(hisids%ncid,    hisids%id_statx,  xobs(:),            start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid,    hisids%id_staty,  yobs(:),            start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
#ifdef HAVE_PROJ
          if (add_latlon) then
!             call transform_and_put_latlon_coordinates(hisids%ncid, hisids%id_statlon, hisids%id_statlat, nccrs%proj_string, xobs, yobs, start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
#endif
       else
          ierr = nf90_put_var(hisids%ncid,    hisids%id_statx,  xobs(:),            start = (/ 1 /), count = (/ ntot /))
          ierr = nf90_put_var(hisids%ncid,    hisids%id_staty,  yobs(:),            start = (/ 1 /), count = (/ ntot /))
#ifdef HAVE_PROJ
          if (add_latlon) then
!             call transform_and_put_latlon_coordinates(hisids%ncid, hisids%id_statlon, hisids%id_statlat, nccrs%proj_string, xobs, yobs)
          endif
#endif
       endif
    endif

    if (timon) call timstrt('unc_write_his obs data 1', handle_extra(56))
    if (japatm > 0 .and. janormalhis > 0) then
       ierr = nf90_put_var(hisids%ncid, hisids%id_varpatm, valobsT(:,IPNT_patm), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
    endif

    if (jawind > 0 .and. janormalhis > 0) then
       ierr = nf90_put_var(hisids%ncid, hisids%id_varwx,  valobsT(:,IPNT_wx),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_varwy,  valobsT(:,IPNT_wy),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
    endif

    if ((jarain > 0) .and. (jahisrain > 0) .and. janormalhis > 0) then
       ierr = nf90_put_var(hisids%ncid, hisids%id_varrain,  valobsT(:,IPNT_rain),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
    endif

    if ((infiltrationmodel == DFM_HYD_INFILT_CONST .or. infiltrationmodel == DFM_HYD_INFILT_HORTON) .and. jahisinfilt > 0 .and. janormalhis > 0) then
       ierr = nf90_put_var(hisids%ncid, hisids%id_infiltcap,  valobsT(:,IPNT_infiltcap),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_infiltact,  valobsT(:,IPNT_infiltact),    start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
    endif
    if (timon) call timstop(handle_extra(56))

    if (timon) call timstrt('unc_write_his obs/crs data 2', handle_extra(57))

    if (numobs+nummovobs > 0) then
!      separate depth-averaged velocities
       if (jahisvelvec > 0 .and. kmx > 0 .and. janormalhis > 0) then
          ierr = nf90_put_var(hisids%ncid,    hisids%id_varucxq, valobsT(:,IPNT_UCXQ),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /)) ! depth-averaged velocity
          ierr = nf90_put_var(hisids%ncid,    hisids%id_varucyq, valobsT(:,IPNT_UCYQ),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       end if
!      3D
       if ( kmx>0 .and. ja3d > 0) then ! 3D only:
          if (jahisvelvec > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varucx, IPNT_UCX, ja3d)
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varucy, IPNT_UCY, ja3d)
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varucz, IPNT_UCZ, ja3d)
          end if
       end if

       if (kmx>0) then ! 3D or depth-averaged versions:
          if (jasal > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varsal, IPNT_SA1, ja3d)
          end if
          if (jatem > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_vartem, IPNT_TEM1, ja3d)
          end if
          if (IVAL_TRA1 > 0) then
             do j = IVAL_TRA1,IVAL_TRAN   ! enumerators of tracers in valobs array (not the pointer)
             i = j - IVAL_TRA1 + 1
             ierr = nf90_put_var(hisids%ncid, hisids%id_tra(i), valobsT(:,IPNT_TRA1 + (i-1)*kmx+kk-1), start = (/ kk, 1, hisids%idx_curtime /), count = (/ 1, ntot, 1/))
             enddo
          end if
          if (jahisvelocity > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varumag, IPNT_UMAG, ja3d)
          end if
       end if

       if ( kmx>0 .and. ja3d > 0) then
          if( (jasal > 0 .or. jatem > 0 .or. jased > 0 ) .and. jahisrho > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varrhop, IPNT_RHOP, ja3d)
             if (idensform > 10) then
                ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varrho, IPNT_RHO, ja3d)
             endif
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_bruv, IPNT_BRUV, ja3d) ! TODO
          end if
          if (jased > 0 .and. .not. stm_included) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varsed, IPNT_SED, ja3d)
          end if
          if (jahisdischarge > 0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_varqmag, IPNT_QMAG, ja3d)
          end if
          ! TODO
          !if (IVAL_HWQ1 > 0) then
          !   do j = IVAL_HWQ1,IVAL_HWQN   ! enumerators of waq output in valobs array (not the pointer)
          !      i = j - IVAL_HWQ1 + 1
          !      if (i .le. noout_user + noout_statt) then
          !         ierr = nf90_put_var(hisids%ncid, hisids%id_hwq(i), valobsT(:,IPNT_HWQ1 + (i-1)*kmx+kk-1), start = (/ kk, 1, hisids%idx_curtime /), count = (/ 1, ntot, 1/))
          !      else if (comparereal(tim, ti_hise, eps10) == 0) then
          !         ierr = nf90_put_var(hisids%ncid, hisids%id_hwq(i), valobsT(:,IPNT_HWQ1 + (i-1)*kmx+kk-1), start = (/ kk, 1 /), count = (/ 1, ntot, 1/))
          !      endif
          !   enddo       
          !end if
          if (IVAL_WQB3D1 > 0) then
             do j = IVAL_WQB3D1,IVAL_WQB3DN   ! enumerators of 3d waqbot output in valobs array (not the pointer)
                i = j - IVAL_WQB3D1 + 1
                ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_hwqb3d(i), IPNT_WQB3D1 + (i-1)*kmx, ja3d)
             enddo
          end if
          if (IVAL_SF1 > 0) then ! TODO
             do kk = 1,kmx
                call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
                do j = IVAL_SF1,IVAL_SFN
                  i = j - IVAL_SF1 + 1
                  toutputx(:,i) = valobsT(:,IPNT_SF1 + (i-1)*(kmx+1)+kk-1)
                enddo
                ierr = nf90_put_var(hisids%ncid, hisids%id_sf, toutputx, start = (/ kk, 1, 1, hisids%idx_curtime /), count = (/ 1, ntot, stmpar%lsedsus, 1/))
             end do
          end if
          if (jawave>0) then
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_ustx, IPNT_UCXST, ja3d)
             ierr = unc_put_3d_or_depthaveraged_valobs(hisids, hisids%id_usty, IPNT_UCYST, ja3d)
          endif
       end if
       if (kmx == 0) then
!      2D
          if (jahisvelvec > 0) then
             ierr = nf90_put_var(hisids%ncid,    hisids%id_varucx, valobsT(:,IPNT_UCX),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
             ierr = nf90_put_var(hisids%ncid,    hisids%id_varucy, valobsT(:,IPNT_UCY),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          if (jahisvelocity > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_varumag, valobsT(:,IPNT_UMAG),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          if (jahisdischarge > 0) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_varqmag, valobsT(:,IPNT_QMAG),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          if (jasal > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_varsal, valobsT(:,IPNT_SA1),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          endif
          if (jatem > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_vartem, valobsT(:,IPNT_TEM1), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
             if( (jasal > 0 .or. jatem > 0 .or. jased > 0 )  .and. jahisrho > 0) then
                ierr = nf90_put_var(hisids%ncid, hisids%id_varrhop, valobsT(:,IPNT_RHOP) ,  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if

          if (IVAL_TRA1 > 0) then
             do j = IVAL_TRA1,IVAL_TRAN   ! enumerators of tracers in valobs array (not the pointer)
               i = j - IVAL_TRA1 + 1
               ierr = nf90_put_var(hisids%ncid, hisids%id_tra(i), valobsT(:,IPNT_TRA1 + i-1), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1/))
             end do
          end if
          if (IVAL_HWQ1 > 0) then
             do j = IVAL_HWQ1,IVAL_HWQN   ! enumerators of extra waq output in valobs array (not the pointer)
               i = j - IVAL_HWQ1 + 1
               ierr = nf90_put_var(hisids%ncid, hisids%id_hwq(i), valobsT(:,IPNT_HWQ1 + i-1), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1/))
             end do
          end if
          if (IVAL_WQB3D1 > 0) then
             do j = IVAL_WQB3D1,IVAL_WQB3DN   ! enumerators of waqbot variables in valobs array (not the pointer)
               i = j - IVAL_WQB3D1 + 1
               ierr = nf90_put_var(hisids%ncid, hisids%id_hwqb3d(i), valobsT(:,IPNT_WQB3D1 + i-1), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1/))
             end do
          end if
          if (IVAL_SF1 > 0) then
             call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
             do j = IVAL_SF1,IVAL_SFN
               i = j - IVAL_SF1 + 1
               toutputx(:,i) = valobsT(:,IPNT_SF1 + i-1)
             end do
             ierr = nf90_put_var(hisids%ncid, hisids%id_sf, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedsus, 1/))
          end if
          if (IVAL_WS1 > 0) then
             call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
             do j = IVAL_WS1,IVAL_WSN
               i = j - IVAL_WS1 + 1
               toutputx(:,i) = valobsT(:,IPNT_WS1 + i-1)
               ierr = nf90_put_var(hisids%ncid, hisids%id_ws, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedsus, 1/))
             enddo
          end if
          !
          if (jased > 0 .and. .not. stm_included .and. janormalhis > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_varsed, valobsT(:,IPNT_SED),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          end if
          !
          if (jawave>0 .and. janormalhis > 0) then
             ierr = nf90_put_var(hisids%ncid,    hisids%id_ustx, valobsT(:,IPNT_UCXST),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
             ierr = nf90_put_var(hisids%ncid,    hisids%id_usty, valobsT(:,IPNT_UCYST),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          endif
       endif

       if (jahistaucurrent>0 .and. janormalhis > 0) then
          ierr = nf90_put_var(hisids%ncid, hisids%id_TAUX,   valobsT(:,IPNT_TAUX), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_TAUY,   valobsT(:,IPNT_TAUY), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif

       if ( jawave.eq.4  .and. janormalhis > 0) then
          ierr = nf90_put_var(hisids%ncid, hisids%id_R,      valobsT(:,IPNT_WAVER), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       end if

       if (jawave>0 .and. janormalhis > 0) then
          ierr = nf90_put_var(hisids%ncid, hisids%id_WH,      valobsT(:,IPNT_WAVEH),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_WD,      valobsT(:,IPNT_WAVED),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_WL,      valobsT(:,IPNT_WAVEL),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_WT,      valobsT(:,IPNT_WAVET),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_WU,      valobsT(:,IPNT_WAVEU),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif

      !    waq bottom variables are always 2D
       if (IVAL_WQB1 > 0) then ! TODO
          do j = IVAL_WQB1,IVAL_WQBN   ! enumerators of tracers in valobs array (not the pointer)
             i = j - IVAL_WQB1 + 1
             ierr = nf90_put_var(hisids%ncid, hisids%id_hwqb(i), valobsT(:,IPNT_WQB1 + i-1), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1/))
          end do
       endif
    endif ! numobs + nummovobs > 0

    if (jatem > 1 .and. jahisheatflux > 0 .and. janormalhis > 0) then
       ierr = nf90_put_var(hisids%ncid,    hisids%id_Wind   , valobsT(:,IPNT_WIND),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))

       if ( jatem.gt.1 ) then   ! also heat modelling involved
          ierr = nf90_put_var(hisids%ncid, hisids%id_Tair   , valobsT(:,IPNT_TAIR),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       end if

       if (jatem == 5 .and. allocated(Rhum) .and. allocated(Clou) ) then
           ierr = nf90_put_var(hisids%ncid, hisids%id_Rhum   , valobsT(:,IPNT_RHUM),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
           ierr = nf90_put_var(hisids%ncid, hisids%id_Clou   , valobsT(:,IPNT_CLOU),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       end if

       if (jatem == 5 ) then
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qsun   , valobsT(:,IPNT_QSUN),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qeva   , valobsT(:,IPNT_QEVA),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qcon   , valobsT(:,IPNT_QCON),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qlong  , valobsT(:,IPNT_QLON),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qfreva , valobsT(:,IPNT_QFRE),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_Qfrcon , valobsT(:,IPNT_QFRC),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif

       ierr = nf90_put_var(hisids%ncid,    hisids%id_Qtot   , valobsT(:,IPNT_QTOT),  start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))

    end if ! jahisheatflux > 0! jatem > 0

    ! 3d layer interface quantities
    if (kmx > 0 .and. ja3d > 0) then
       do kk = 1, kmx+1
          ierr = nf90_put_var(hisids%ncid,    hisids%id_zws,    valobsT(:,IPNT_ZWS+kk-1),   start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          ierr = nf90_put_var(hisids%ncid,    hisids%id_zwu,    valobsT(:,IPNT_ZWU+kk-1),   start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          if (kk > 1) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_zcs,    valobsT(:,IPNT_ZCS+kk-2),   start = (/ kk-1,1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          endif
          if (iturbulencemodel >= 3 .and. jahistur > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_turkin, valobsT(:,IPNT_TKIN +kk-1), start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
             ierr = nf90_put_var(hisids%ncid, hisids%id_tureps, valobsT(:,IPNT_TEPS +kk-1), start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          end if
          if (iturbulencemodel > 1) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_vicwwu, valobsT(:,IPNT_VICWW+kk-1), start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          end if
          if (idensform > 0 .and. jaRichardsononoutput > 0) then
             ierr = nf90_put_var(hisids%ncid, hisids%id_rich,   valobsT(:,IPNT_RICH +kk-1), start = (/ kk,  1, hisids%idx_curtime /), count = (/ 1, ntot, 1 /))
          end if
          !
          if (IVAL_WS1 > 0) then
             call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
             do j = IVAL_WS1,IVAL_WSN
             i = j - IVAL_WS1 + 1
             toutputx(:,i) = valobsT(:,IPNT_WS1 + (i-1)*(kmx+1)+kk-1)
             enddo
             ierr = nf90_put_var(hisids%ncid, hisids%id_ws, toutputx, start = (/ kk, 1, 1, hisids%idx_curtime /), count = (/ 1, ntot, stmpar%lsedsus, 1/))
          end if
          !
          if (IVAL_SEDDIF1 > 0) then
             call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
             do j = IVAL_SEDDIF1,IVAL_SEDDIFN
             i = j - IVAL_SEDDIF1 + 1
             toutputx(:,i) = valobsT(:,IPNT_SEDDIF1 + (i-1)*(kmx+1)+kk-1)
             enddo
             ierr = nf90_put_var(hisids%ncid, hisids%id_seddif, toutputx, start = (/ kk, 1, 1, hisids%idx_curtime /), count = (/ 1, ntot, stmpar%lsedsus, 1/))
          end if
          !
       enddo
    endif
    !
    ! Bed composition variables
    if (jahissed>0 .and. jased>0 .and. stm_included .and. janormalhis > 0) then
       if (stmpar%morpar%moroutput%taub) then
          ierr = nf90_put_var(hisids%ncid, hisids%id_taub, valobsT(:,IPNT_TAUB), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif
       !
       select case (stmpar%morlyr%settings%iunderlyr)
          case (1)
             ! dpsed
             ierr = nf90_put_var(hisids%ncid, hisids%id_dpsed, valobsT(:,IPNT_DPSED),   start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
             ! bodsed
             call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
             do j = IVAL_BODSED1, IVAL_BODSEDN
               i = j - IVAL_BODSED1 + 1
               toutputx(:,i) = valobsT(:,IPNT_BODSED1 + i-1)
             end do
             ierr = nf90_put_var(hisids%ncid, hisids%id_bodsed, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1/))
          case (2)
             nlyrs = stmpar%morlyr%settings%nlyr
             ! msed
             call realloc(toutput3, (/nlyrs, ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
             toutput3 = dmiss
             do j = IVAL_MSED1,IVAL_MSEDN
                i = j - IVAL_MSED1 + 1
                do kk = 1,nlyrs
                   toutput3(kk,:,i) = valobsT(:,IPNT_MSED1 + (i-1)*(nlyrs)+kk-1)
                enddo
             enddo
             ierr = nf90_put_var(hisids%ncid, hisids%id_msed, toutput3, start = (/ 1, 1, 1, hisids%idx_curtime /), count = (/ nlyrs, ntot, stmpar%lsedtot, 1/))
             ! lyrfrac
             toutput3=dmiss
             do j = IVAL_LYRFRAC1,IVAL_LYRFRACN
                i = j - IVAL_LYRFRAC1 + 1
                do kk = 1,nlyrs
                   toutput3(kk,:,i) = valobsT(:,IPNT_LYRFRAC1 + (i-1)*(nlyrs)+kk-1)
                enddo
             enddo
             ierr = nf90_put_var(hisids%ncid, hisids%id_lyrfrac, toutput3, start = (/ 1, 1, 1, hisids%idx_curtime /), count = (/ nlyrs, ntot, stmpar%lsedtot, 1/))
             ! thlyr
             call realloc(toutputx, (/nlyrs, ntot /), keepExisting=.false., fill = dmiss)
             do kk = 1,nlyrs
                toutputx(kk,:) = valobsT(:,IPNT_THLYR+kk-1)
             enddo
             ierr = nf90_put_var(hisids%ncid, hisids%id_thlyr, toutputx,  start = (/ 1, 1, hisids%idx_curtime /), count = (/ nlyrs, ntot, 1 /))
             ! poros
             if (stmpar%morlyr%settings%iporosity>0) then
                do kk = 1,nlyrs
                   toutputx(kk,:) = valobsT(:,IPNT_POROS+kk-1)
                enddo
                ierr = nf90_put_var(hisids%ncid, hisids%id_poros, toutputx,  start = (/ 1, 1, hisids%idx_curtime /), count = (/ nlyrs, ntot, 1 /))
             endif
       end select
       !
       if (stmpar%morpar%moroutput%frac) then
          ! frac
          call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
          do j = IVAL_FRAC1,IVAL_FRACN
            i = j - IVAL_FRAC1 + 1
            toutputx(:,i) = valobsT(:,IPNT_FRAC1 + i-1)
          enddo
          ierr = nf90_put_var(hisids%ncid, hisids%id_frac, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1/))
       endif
       if (stmpar%morpar%moroutput%mudfrac) then
          ! mudfrac
          ierr = nf90_put_var(hisids%ncid, hisids%id_mudfrac, valobsT(:,IPNT_MUDFRAC), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif
       if (stmpar%morpar%moroutput%sandfrac) then
          ! sandfrac
          ierr = nf90_put_var(hisids%ncid, hisids%id_sandfrac, valobsT(:,IPNT_SANDFRAC), start = (/ 1, hisids%idx_curtime /), count = (/ ntot, 1 /))
       endif
       if (stmpar%morpar%moroutput%fixfac) then
          !fixfac
          call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
          do j = IVAL_FIXFAC1,IVAL_FIXFACN
            i = j - IVAL_FIXFAC1 + 1
            toutputx(:,i) = valobsT(:,IPNT_FIXFAC1 + i-1)
          enddo
          ierr = nf90_put_var(hisids%ncid, hisids%id_fixfac, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1/))
       endif
       if (stmpar%morpar%moroutput%hidexp) then
          ! hidexp
          call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
          do j = IVAL_HIDEXP1,IVAL_HIDEXPN
            i = j - IVAL_HIDEXP1 + 1
            toutputx(:,i) = valobsT(:,IPNT_HIDEXP1 + i-1)
          enddo
          ierr = nf90_put_var(hisids%ncid, hisids%id_hidexp, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1/))
       endif
       !
       if (stmpar%morpar%flufflyr%iflufflyr>0 .and. stmpar%lsedsus>0) then
          call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
          do j = IVAL_MFLUFF1,IVAL_MFLUFFN
            i = j - IVAL_MFLUFF1 + 1
            toutputx(:,i) = valobsT(:,IPNT_MFLUFF1 + i-1)
          enddo
          ierr = nf90_put_var(hisids%ncid, hisids%id_mfluff, toutputx, start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedsus, 1/))
       end if
    endif

    !
    ! Cross sections
    if (ncrs > 0 .and. janormalhis > 0) then
       do i=1,ncrs
          ! Discharges Q
          ierr = nf90_put_var(hisids%ncid, hisids%id_varQ,    crs(i)%sumvalcur(IPNT_Q1C), (/ i, hisids%idx_curtime /))
          ierr = nf90_put_var(hisids%ncid, hisids%id_varQint, crs(i)%sumvalcum(IPNT_Q1C), (/ i, hisids%idx_curtime /))
!          ierr = nf90_put_var(hisids%ncid, hisids%id_varQavg, crs(i)%sumvalavg(IPNT_Q1C), (/ i, hisids%idx_curtime /))

          ! Cross sectional areas A*u
          ierr = nf90_put_var(hisids%ncid, hisids%id_varAu,    crs(i)%sumvalcur(IPNT_AUC), (/ i, hisids%idx_curtime /))
!          ierr = nf90_put_var(hisids%ncid, hisids%id_varAuavg, crs(i)%sumvalavg(IPNT_AUC), (/ i, hisids%idx_curtime /))

          ! Average velocity Q/Au
          ierr = nf90_put_var(hisids%ncid, hisids%id_varu,    crs(i)%sumvalcur(IPNT_U1A), (/ i, hisids%idx_curtime /))
!          ierr = nf90_put_var(hisids%ncid, hisids%id_varuavg, crs(i)%sumvalavg(IPNT_U1A), (/ i, hisids%idx_curtime /))

             IP = IPNT_HUA
             do num = 1,NUMCONST_MDU
                IP = IP + 1
                if (num >= ISED1 .and. num <= ISEDN) then
                   l = sedtot2sedsus(num-ISED1+1)
                   select case(stmpar%morpar%moroutput%transptype)
                   case (0)
                      rhol = 1d0
                   case (1)
                      rhol = stmpar%sedpar%cdryb(l)
                   case (2)
                      rhol = stmpar%sedpar%rhosol(l)
                   end select
                   toutput_cum = crs(i)%sumvalcum(IP)/rhol
                   toutput_cur = crs(i)%sumvalcur(IP)/rhol
                else
                  toutput_cum = crs(i)%sumvalcum(IP)
                  toutput_cur = crs(i)%sumvalcur(IP)
                endif
                ierr = nf90_put_var(hisids%ncid, hisids%id_const_cum(num), toutput_cum, (/ i, hisids%idx_curtime /))
                ierr = nf90_put_var(hisids%ncid, hisids%id_const(num),     toutput_cur, (/ i, hisids%idx_curtime /))
             end do
     
          if( jased == 4 .and. stmpar%lsedtot > 0 ) then
             IP = IPNT_HUA + NUMCONST_MDU + 1
             ierr = nf90_put_var(hisids%ncid, hisids%id_sedbtrans, crs(i)%sumvalcum(IP), (/ i, hisids%idx_curtime /))
             if( stmpar%lsedsus > 0 ) then
                IP = IP + 1
                ierr = nf90_put_var(hisids%ncid, hisids%id_sedstrans, crs(i)%sumvalcum(IP), (/ i, hisids%idx_curtime /))
             endif
             do lsed = 1,stmpar%lsedtot    ! Making bedload on crosssections per fraction
                IP = IP + 1
                ierr = nf90_put_var(hisids%ncid, hisids%id_sedbtransfrac(lsed), crs(i)%sumvalcum(IP), (/ i, hisids%idx_curtime /))
             enddo
          endif
       end do
    end if
    if (timon) call timstop(handle_extra(57))

    if (timon) call timstrt('unc_write_his RUG', handle_extra(65))

    ! runup gauges
    if (nrug>0 .and. janormalhis > 0) then
       call realloc(geom_x,   nrug, fill=dmiss)
       call realloc(geom_y,   nrug, fill=dmiss)
       call realloc(toutput1, nrug, fill=dmiss)
       do i=1,nrug
          geom_x(i)   = rug(i)%maxx
          geom_y(i)   = rug(i)%maxy
          toutput1(i) = rug(i)%maxruh
       end do
       ierr = nf90_put_var(hisids%ncid, hisids%id_rugx,   geom_x(:),   start = (/ 1, hisids%idx_curtime /), count = (/ nrug, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_rugy,   geom_y(:),   start = (/ 1, hisids%idx_curtime /), count = (/ nrug, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_varruh, toutput1(:), start = (/ 1, hisids%idx_curtime /), count = (/ nrug, 1 /))
    endif

    if (timon) call timstop(handle_extra(65))

    if (timon) call timstrt ( "unc_write_his str write", handle_extra(62))
    if (jahissourcesink > 0 .and. numsrc > 0) then
      if (tim == tstart_user) then
        do i = 1, numsrc
          ierr = nf90_put_var(hisids%ncid, hisids%id_srcname, trimexact(srcname(i), strlen_netcdf), (/ 1, i/) )
          ierr = nf90_put_var(hisids%ncid, hisids%id_qsrccur, qstss((numconst+1)*(i-1)+1), (/ i, hisids%idx_curtime /))
        enddo
        ierr = nf90_put_var(hisids%ncid, hisids%id_srcx, xsrc)
        ierr = nf90_put_var(hisids%ncid, hisids%id_srcy, ysrc)
      else
         ierr = nf90_put_var(hisids%ncid, hisids%id_qsrccur, qsrc, (/ 1, hisids%idx_curtime /))
      endif

      ! write geometry variables at the first time of history output
       if (hisids%idx_curtime == 1) then
          j = 1
          call realloc(node_count, numsrc, fill = 0)
          call realloc(geom_x, 2)
          call realloc(geom_y, 2)
          do i = 1, numsrc
             k1 = ksrc(1,i)
             k2 = ksrc(4,i)
             nNodes = 0
             if (k1 > 0) then
                nNodes = nNodes + 1
                geom_x(nNodes) = xz(k1)
                geom_y(nNodes) = yz(k1)
             end if
             if (k2 > 0) then
                nNodes = nNodes + 1
                geom_x(nNodes) = xz(k2)
                geom_y(nNodes) = yz(k2)
             end if

             node_count(i) = nNodes
             if (nNodes > 0) then
                ierr = nf90_put_var(hisids%ncid, hisids%id_srcgeom_node_coordx,  geom_x(1:nNodes), start = (/ j /), count = (/ nNodes /))
                ierr = nf90_put_var(hisids%ncid, hisids%id_srcgeom_node_coordy,  geom_y(1:nNodes), start = (/ j /), count = (/ nNodes /))
             end if

             j = j + nNodes
          end do
          ierr = nf90_put_var(hisids%ncid, hisids%id_srcgeom_node_count, node_count)
       end if

      ierr = nf90_put_var(hisids%ncid, hisids%id_vsrccum, vsrccum, (/ 1, hisids%idx_curtime /))
      ierr = nf90_put_var(hisids%ncid, hisids%id_qsrcavg, qsrcavg, (/ 1, hisids%idx_curtime /))
      do i = 1, numsrc
        ierr = nf90_put_var(hisids%ncid, hisids%id_pred,  qstss((numconst+1)*(i-1) + 1), (/ i, hisids%idx_curtime /))
        j = 1
        if (isalt > 0) then
           j = j + 1
           ierr = nf90_put_var(hisids%ncid, hisids%id_presa, qstss((numconst+1)*(i-1) + j), (/ i, hisids%idx_curtime /))
        endif
        if (itemp > 0) then
           j = j + 1
           ierr = nf90_put_var(hisids%ncid, hisids%id_pretm, qstss((numconst+1)*(i-1) + j), (/ i, hisids%idx_curtime /))
        endif
      enddo
    endif

      if (jahiscgen > 0 .and. janormalhis > 0) then
         if (ncgensg > 0) then
            do i = 1,ncgensg
               igen = i
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis   , valcgen(2,i)   , (/ i, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_crestl, zcgen(3*igen-2), (/ i, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_edgel , zcgen(3*igen-1), (/ i, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_openw , zcgen(3*igen  ), (/ i, hisids%idx_curtime /)) ! TODO: AvD: this part seems not entirely correct, double check with block below and duplication with gategen, etc.
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1up  , valcgen(3,i)   , (/ i, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1dn  , valcgen(4,i)   , (/ i, hisids%idx_curtime /))
            enddo
         end if
         if (ngenstru > 0) then
            valobsT(1:ngenstru, 1:NUMVALS_GENSTRU) = transpose(valgenstru)
            !do i=1,ngenstru
            !   !igen = genstru2cgen(i)
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis   , valgenstru(2,i),  (/ i, hisids%idx_curtime /))
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_crestl, valgenstru(9,i),  (/ i, hisids%idx_curtime /)) ! changed
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_edgel , valgenstru(14,i), (/ i, hisids%idx_curtime /)) ! changed
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_openw , valgenstru(13,i), (/ i, hisids%idx_curtime /)) ! changed
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1up  , valgenstru(3,i),  (/ i, hisids%idx_curtime /))
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1dn  , valgenstru(4,i),  (/ i, hisids%idx_curtime /))
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_head,          valgenstru(5,i),  (/ i, hisids%idx_curtime /))
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au,            valgenstru(6,i),  (/ i, hisids%idx_curtime /))
            !   ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_vel,           valgenstru(7,i),  (/ i, hisids%idx_curtime /))
            !   if (network%sts%numGeneralStructures > 0) then
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1crest,       valgenstru(8,i),  (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_crestw,        valgenstru(10,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_stat,     int(valgenstru(11,i)), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_forcedif,      valgenstru(12,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_openh,         valgenstru(15,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_uppl,          valgenstru(16,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_open, valgenstru(17,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_over, valgenstru(18,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_under,valgenstru(19,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_open,       valgenstru(20,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_over,        valgenstru(21,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_under,      valgenstru(22,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateopen,   valgenstru(23,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateover,   valgenstru(24,i), (/ i, hisids%idx_curtime /))
            !      ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateunder,  valgenstru(25,i), (/ i, hisids%idx_curtime /))
            !   end if
            !enddo
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis   , valobsT(1:ngenstru,2),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_crestl, valobsT(1:ngenstru,9),  (/ 1, hisids%idx_curtime /)) ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_edgel , valobsT(1:ngenstru,14), (/ 1, hisids%idx_curtime /)) ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_openw , valobsT(1:ngenstru,13), (/ 1, hisids%idx_curtime /)) ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1up  , valobsT(1:ngenstru,3),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1dn  , valobsT(1:ngenstru,4),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_head,          valobsT(1:ngenstru,5),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au,            valobsT(1:ngenstru,6),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_vel,           valobsT(1:ngenstru,7),  (/ 1, hisids%idx_curtime /))
               if (network%sts%numGeneralStructures > 0) then
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_s1crest,       valobsT(1:ngenstru,8),  (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_crestw,        valobsT(1:ngenstru,10), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_stat,      int(valobsT(1:ngenstru,11)),(/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_forcedif,      valobsT(1:ngenstru,12), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_openh,         valobsT(1:ngenstru,15), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_uppl,          valobsT(1:ngenstru,16), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_open, valobsT(1:ngenstru,17), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_over, valobsT(1:ngenstru,18), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_dis_gate_under,valobsT(1:ngenstru,19), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_open,       valobsT(1:ngenstru,20), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_over,       valobsT(1:ngenstru,21), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_au_under,      valobsT(1:ngenstru,22), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateopen,   valobsT(1:ngenstru,23), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateover,   valobsT(1:ngenstru,24), (/ 1, hisids%idx_curtime /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_genstru_velgateunder,  valobsT(1:ngenstru,25), (/ 1, hisids%idx_curtime /))
               end if
            ! write geometry variables at the first time of history output
            if (hisids%idx_curtime == 1) then
               if (network%sts%numGeneralStructures > 0) then ! new general structure
                  ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_coordx, geomXGenstru,     start = (/ 1 /), count = (/ nNodesGenstru /))
                  ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_coordy, geomYGenstru,     start = (/ 1 /), count = (/ nNodesGenstru /))
                  ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_count,  nodeCountGenstru, start = (/ 1 /), count = (/ network%sts%numGeneralStructures /))
                  if (allocated(geomXGenstru))     deallocate(geomXGenstru) ! Deallocate the geometry arrays after writing them to his-file.
                  if (allocated(geomYGenstru))     deallocate(geomYGenstru)
                  if (allocated(nodeCountGenstru)) deallocate(nodeCountGenstru)
               else ! old general structure
                  j = 1
                  call realloc(node_count, ngenstru, fill = 0)
                  do n = 1, ngenstru
                     i = genstru2cgen(n)
                     nlinks = L2cgensg(i) - L1cgensg(i) + 1
                     if (nlinks > 0) then
                        nNodes = nlinks + 1
                     else if (nlinks == 0) then
                        nNodes = 0
                     end if
                     node_count(n) = nNodes

                     if (nNodes > 0) then
                        call get_geom_coordinates_of_structure_old(i, nNodes, geom_x, geom_y)
                        ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_coordx, geom_x(1:nNodes), start = (/ j /), count = (/ nNodes /))
                        ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_coordy, geom_y(1:nNodes), start = (/ j /), count = (/ nNodes /))
                        j = j + nNodes
                     end if
                  end do
                  ierr = nf90_put_var(hisids%ncid, hisids%id_genstrugeom_node_count, node_count, start = (/ 1 /), count = (/ ngenstru /))
               end if
            end if
         endif
      endif

      if (jahispump > 0 .and. npumpsg > 0 .and. janormalhis > 0) then
         valobsT(1:npumpsg, 1:NUMVALS_PUMP) = transpose(valpump)
         !do i=1,npumpsg
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_dis,     valpump(2,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1up,    valpump(3,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1dn,    valpump(4,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_struhead,valpump(5,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_cap,     valpump(6,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_disdir,  valpump(12,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_stage,int(valpump(7,i)),(/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_head,    valpump(8,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_redufact,valpump(9,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1del,   valpump(10,i),(/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1suc,   valpump(11,i),(/ i, hisids%idx_curtime /))
         !end do
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_dis,     valobsT(1:npumpsg,2), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1up,    valobsT(1:npumpsg,3), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1dn,    valobsT(1:npumpsg,4), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_struhead,valobsT(1:npumpsg,5), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_cap,     valobsT(1:npumpsg,6), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_disdir,  valobsT(1:npumpsg,12), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_stage,int(valobsT(1:npumpsg,7)),(/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_head,    valobsT(1:npumpsg,8), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_redufact,valobsT(1:npumpsg,9), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1del,   valobsT(1:npumpsg,10),(/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_pump_s1suc,   valobsT(1:npumpsg,11),(/ 1, hisids%idx_curtime /))
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            if (network%sts%numPumps > 0) then ! new pump
               ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_coordx, geomXPump,     start = (/ 1 /), count = (/ nNodesPump /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_coordy, geomYPump,     start = (/ 1 /), count = (/ nNodesPump /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_count,  nodeCountPump, start = (/ 1 /), count = (/ network%sts%numPumps /))
               if (allocated(geomXPump))     deallocate(geomXPump)
               if (allocated(geomYPump))     deallocate(geomYPump)
               if (allocated(nodeCountPump)) deallocate(nodeCountPump)
            else
               j = 1
               call realloc(node_count, npumpsg, fill = 0)
               do i = 1, npumpsg
                  nlinks = L2pumpsg(i) - L1pumpsg(i) + 1
                  if (nlinks > 0) then
                     nNodes = nlinks + 1
                  else if (nlinks == 0) then
                     nNodes = 0
                  end if
                  node_count(i) = nNodes

                  if (nNodes > 0) then
                     call realloc(geom_x, nNodes)
                     call realloc(geom_y, nNodes)

                     L0 = L1pumpsg(i)
                     L = abs(kpump(3,L0))
                     k1 = lncn(1,L)
                     k2 = lncn(2,L)
                     geom_x(1) = xk(k1)
                     geom_x(2) = xk(k2)
                     geom_y(1) = yk(k1)
                     geom_y(2) = yk(k2)
                     if (nlinks > 1) then
                        k = 3
                        do L0 = L1pumpsg(i)+1, L2pumpsg(i)
                           L = abs(kpump(3,L0))
                           k3 = lncn(2,L)
                           geom_x(k) = xk(k3)
                           geom_y(k) = yk(k3)
                           k = k+1
                        end do
                     end if

                     ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_coordx, geom_x(1:nNodes), start = (/ j /), count = (/ nNodes /))
                     ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_coordy, geom_y(1:nNodes), start = (/ j /), count = (/ nNodes /))
                     j = j + nNodes
                  end if
               end do
               ierr = nf90_put_var(hisids%ncid, hisids%id_pumpgeom_node_count, node_count, start = (/ 1 /), count = (/ npumpsg /))
            end if
         end if
      end if

      if (jahisorif > 0 .and. network%sts%numOrifices > 0 .and. janormalhis > 0) then
         do i=1,network%sts%numOrifices
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_dis   ,        valorifgen(2,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_s1up  ,        valorifgen(3,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_s1dn  ,        valorifgen(4,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_head,          valorifgen(5,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_au,            valorifgen(6,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_vel,           valorifgen(7,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_s1crest,       valorifgen(8,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_crestl,        valorifgen(9,i),   (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_crestw,        valorifgen(10,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_stat,      int(valorifgen(11,i)), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_forcedif,      valorifgen(12,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_edgel ,        valorifgen(14,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgen_openh,         valorifgen(15,i),  (/ i, hisids%idx_curtime /))
         enddo
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgeom_node_coordx, geomXOrif,     start = (/ 1 /), count = (/ nNodesOrif /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgeom_node_coordy, geomYOrif,     start = (/ 1 /), count = (/ nNodesOrif /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_orifgeom_node_count,  nodeCountOrif, start = (/ 1 /), count = (/ network%sts%numOrifices /))
            if (allocated(geomXOrif))     deallocate(geomXOrif)
            if (allocated(geomYOrif))     deallocate(geomYOrif)
            if (allocated(nodeCountOrif)) deallocate(nodeCountOrif)
               end if
         end if

      if (timon) call timstrt('unc_write_his bridge data', handle_extra(58))
      if (jahisbridge > 0 .and. network%sts%numBridges > 0 .and. janormalhis > 0) then
         !do i=1,network%sts%numBridges
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_dis,   valbridge(2,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_s1up,  valbridge(3,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_s1dn,  valbridge(4,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_head,  valbridge(5,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_au,    valbridge(6,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_vel,   valbridge(7,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_blup,  valbridge(8,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_bldn,  valbridge(9,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_bl_act,valbridge(10,i),(/ i, hisids%idx_curtime /))
         !enddo
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_dis,   valbridge(2, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_s1up,  valbridge(3, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_s1dn,  valbridge(4, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_head,  valbridge(5, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_au,    valbridge(6, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_vel,   valbridge(7, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_blup,  valbridge(8, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_bldn,  valbridge(9, 1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_bridge_bl_act,valbridge(10,1:network%sts%numBridges), (/ 1, hisids%idx_curtime /))
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            if (timon) call timstrt('Bridge geom', handle_extra(74))
               ierr = nf90_put_var(hisids%ncid, hisids%id_bridgegeom_node_coordx, geomXBridge,     start = (/ 1 /), count = (/ nNodesBridge /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_bridgegeom_node_coordy, geomYBridge,     start = (/ 1 /), count = (/ nNodesBridge /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_bridgegeom_node_count,  nodeCountBridge, start = (/ 1 /), count = (/ network%sts%numBridges /))
               if (allocated(geomXBridge))     deallocate(geomXBridge)
               if (allocated(geomYBridge))     deallocate(geomYBridge)
               if (allocated(nodeCountBridge)) deallocate(nodeCountBridge)
            if (timon) call timstop(handle_extra(74))
               end if
         end if
      if (timon) call timstop(handle_extra(58))

      if (jahisculv > 0 .and. network%sts%numCulverts > 0 .and. janormalhis > 0) then
         do i=1,network%sts%numCulverts
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_dis,    valculvert(2,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_s1up,   valculvert(3,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_s1dn,   valculvert(4,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_head,   valculvert(5,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_au,     valculvert(6,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_vel,    valculvert(7,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_crestl, valculvert(8,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_stat,  int(valculvert(9,i)),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_edgel , valculvert(10,i),     (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvert_openh,  valculvert(11,i),     (/ i, hisids%idx_curtime /))
         enddo
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvertgeom_node_coordx, geomXCulv,     start = (/ 1 /), count = (/ nNodesCulv /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvertgeom_node_coordy, geomYCulv,     start = (/ 1 /), count = (/ nNodesCulv /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_culvertgeom_node_count,  nodeCountCulv, start = (/ 1 /), count = (/ network%sts%numCulverts /))
            if (allocated(geomXCulv))     deallocate(geomXCulv)
            if (allocated(geomYCulv))     deallocate(geomYCulv)
            if (allocated(nodeCountCulv)) deallocate(nodeCountCulv)
               end if
         end if

      if (jahisuniweir > 0 .and. network%sts%numuniweirs > 0 .and. janormalhis > 0) then
         do i=1,network%sts%numuniweirs
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_dis,    valuniweir(2,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_s1up,   valuniweir(3,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_s1dn,   valuniweir(4,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_head,   valuniweir(5,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_au,     valuniweir(6,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_vel,    valuniweir(7,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweir_crestl, valuniweir(8,i),      (/ i, hisids%idx_curtime /))
         enddo
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweirgeom_node_coordx, geomXUniweir,     start = (/ 1 /), count = (/ nNodesUniweir /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweirgeom_node_coordy, geomYUniweir,     start = (/ 1 /), count = (/ nNodesUniweir /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_uniweirgeom_node_count,  nodeCountUniweir, start = (/ 1 /), count = (/ network%sts%numuniweirs /))
            if (allocated(geomXUniweir))     deallocate(geomXUniweir)
            if (allocated(geomYUniweir))     deallocate(geomYUniweir)
            if (allocated(nodeCountUniweir)) deallocate(nodeCountUniweir)
               end if
         end if

      if (jahiscmpstru > 0 .and. network%cmps%count > 0 .and. janormalhis > 0) then
         do i=1,network%cmps%count
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_dis,            valcmpstru(2,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_s1up,           valcmpstru(3,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_s1dn,           valcmpstru(4,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_head,           valcmpstru(5,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_au,             valcmpstru(6,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cmpstru_vel,            valcmpstru(7,i), (/ i, hisids%idx_curtime /))
         enddo
      end if

      if (jahislongculv > 0 .and. nlongculvertsg > 0 .and. janormalhis > 0) then
         do i=1,nlongculvertsg
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_dis,       vallongculvert(2,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_s1up,      vallongculvert(3,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_s1dn,      vallongculvert(4,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_head,      vallongculvert(5,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_au,        vallongculvert(6,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_vel,       vallongculvert(7,i),      (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvert_valveopen, vallongculvert(8,i),      (/ i, hisids%idx_curtime /))
         enddo
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvertgeom_node_coordx, geomXLongCulv,     start = (/ 1 /), count = (/ nNodesLongCulv /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvertgeom_node_coordy, geomYLongCulv,     start = (/ 1 /), count = (/ nNodesLongCulv /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_longculvertgeom_node_count,  nodeCountLongCulv, start = (/ 1 /), count = (/ nlongculvertsg /))
            if (allocated(geomXLongCulv))     deallocate(geomXLongCulv)
            if (allocated(geomYLongCulv))     deallocate(geomYLongCulv)
            if (allocated(nodeCountLongCulv)) deallocate(nodeCountLongCulv)
               end if
         end if

      if (jahislateral > 0 .and. numlatsg > 0 .and. janormalhis > 0) then
         ierr = nf90_put_var(hisids%ncid, hisids%id_lat_predis_inst,  qplat,       start = (/1,hisids%idx_curtime/), count = (/numlatsg,1/))
         ierr = nf90_put_var(hisids%ncid, hisids%id_lat_predis_ave,   qplatAve,    start = (/1,hisids%idx_curtime/), count = (/numlatsg,1/))
         ierr = nf90_put_var(hisids%ncid, hisids%id_lat_realdis_inst, qLatReal,    start = (/1,hisids%idx_curtime/), count = (/numlatsg,1/))
         ierr = nf90_put_var(hisids%ncid, hisids%id_lat_realdis_ave,  qLatRealAve, start = (/1,hisids%idx_curtime/), count = (/numlatsg,1/))
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            ierr = nf90_put_var(hisids%ncid, hisids%id_latgeom_node_coordx, geomXLat(1:nNodesLat), start = (/ 1 /), count = (/ nlatnd /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_latgeom_node_coordy, geomYLat(1:nNodesLat), start = (/ 1 /), count = (/ nlatnd /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_latgeom_node_count,  nodeCountLat)
         end if
      end if

      if (jahisgate > 0 .and. ngatesg > 0 .and. janormalhis > 0) then
         do i=1,ngatesg
            ierr = nf90_put_var(hisids%ncid, hisids%id_gate_dis  , valgate(2,i) , (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_gate_edgel, zgate(i)     , (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_gate_s1up , valgate(3,i) , (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_gate_s1dn , valgate(4,i) , (/ i, hisids%idx_curtime /))
         end do
      end if

      if (jahisgate > 0 .and. ngategen > 0 .and. janormalhis > 0) then
         do i=1,ngategen
            igen = gate2cgen(i)
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_dis  , valgategen(2,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_sillh, valgategen(9,i), (/ i, hisids%idx_curtime /))   ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_edgel, valgategen(8,i), (/ i, hisids%idx_curtime /))   ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_flowh, valgategen(5,i), (/ i, hisids%idx_curtime /))   ! TODO: AvD sillw
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_openw, valgategen(7,i), (/ i, hisids%idx_curtime /))   ! changed
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_s1up , valgategen(3,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategen_s1dn , valgategen(4,i), (/ i, hisids%idx_curtime /))
         end do
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            j = 1
            call realloc(node_count, ngategen, fill = 0)
            do n = 1, ngategen
               i = gate2cgen(n)
               nlinks = L2cgensg(i) - L1cgensg(i) + 1
               if (nlinks > 0) then
                  nNodes = nlinks + 1
               else if (nlinks == 0) then
                  nNodes = 0
               end if
               node_count(n) = nNodes

               if (nNodes > 0) then
                  call get_geom_coordinates_of_structure_old(i, nNodes, geom_x, geom_y)
                  ierr = nf90_put_var(hisids%ncid, hisids%id_gategengeom_node_coordx, geom_x(1:nNodes), start = (/ j /), count = (/ nNodes /))
                  ierr = nf90_put_var(hisids%ncid, hisids%id_gategengeom_node_coordy, geom_y(1:nNodes), start = (/ j /), count = (/ nNodes /))
                  j = j + nNodes
               end if
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_gategengeom_node_count, node_count, start = (/ 1 /), count = (/ ngategen /))
         end if
      end if

      if (jahiscdam > 0 .and. ncdamsg > 0 .and. janormalhis > 0) then
         do i = 1,ncdamsg
            ierr = nf90_put_var(hisids%ncid, hisids%id_cdam_dis   , valcdam(2,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cdam_crestl, zcdam(i)    , (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cdam_s1up  , valcdam(3,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_cdam_s1dn  , valcdam(4,i), (/ i, hisids%idx_curtime /))
         end do
      end if

      if (jahisweir > 0 .and. nweirgen > 0 .and. janormalhis > 0) then
         valobsT(1:nweirgen, 1:NUMVALS_WEIRGEN) = transpose(valweirgen)
         !do i = 1,nweirgen
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_dis   , valweirgen(2,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1up  , valweirgen(3,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1dn  , valweirgen(4,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_crestl, valweirgen(9,i), (/ i, hisids%idx_curtime /))
         !   ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_crestw, valweirgen(10,i),(/ i, hisids%idx_curtime /))
         !   if (network%sts%numWeirs > 0) then ! write extra files for new weirs
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_head  , valweirgen(5,i),  (/ i, hisids%idx_curtime /))
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_au    , valweirgen(6,i),  (/ i, hisids%idx_curtime /))
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_vel   , valweirgen(7,i),  (/ i, hisids%idx_curtime /))
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1crest,valweirgen(8,i),  (/ i, hisids%idx_curtime /))
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weir_stat, int(valweirgen(11,i)),(/ i, hisids%idx_curtime /))
         !      ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_forcedif,valweirgen(12,i),(/i, hisids%idx_curtime /))
         !   end if
         !end do
         ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_dis   , valobsT(1:nweirgen,2), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1up  , valobsT(1:nweirgen,3), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1dn  , valobsT(1:nweirgen,4), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_crestl, valobsT(1:nweirgen,9), (/ 1, hisids%idx_curtime /))
         ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_crestw, valobsT(1:nweirgen,10),(/ 1, hisids%idx_curtime /))
            if (network%sts%numWeirs > 0) then ! write extra files for new weirs
            ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_head  , valobsT(1:nweirgen,5),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_au    , valobsT(1:nweirgen,6),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_vel   , valobsT(1:nweirgen,7),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_s1crest,valobsT(1:nweirgen,8),  (/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_weir_stat, int(valobsT(1:nweirgen,11)),(/ 1, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_weirgen_forcedif,valobsT(1:nweirgen,12),(/1, hisids%idx_curtime /))
            end if
         ! write geometry variables at the first time of history output
         if (hisids%idx_curtime == 1) then
            if (network%sts%numWeirs > 0) then ! new weir
               ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_coordx, geomXWeir,     start = (/ 1 /), count = (/ nNodesWeir /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_coordy, geomYWeir,     start = (/ 1 /), count = (/ nNodesWeir /))
               ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_count,  nodeCountWeir, start = (/ 1 /), count = (/ network%sts%numWeirs /))
               if (allocated(geomXWeir))     deallocate(geomXWeir)
               if (allocated(geomYWeir))     deallocate(geomYWeir)
               if (allocated(nodeCountWeir)) deallocate(nodeCountWeir)
            else
               j = 1
               call realloc(node_count, nweirgen, fill = 0)
               do n = 1, nweirgen
                  i = weir2cgen(n)
                  nlinks = L2cgensg(i) - L1cgensg(i) + 1
                  if (nlinks > 0) then
                     nNodes = nlinks + 1
                  else if (nlinks == 0) then
                     nNodes = 0
                  end if
                  node_count(n) = nNodes

                  if (nNodes > 0) then
                     call get_geom_coordinates_of_structure_old(i, nNodes, geom_x, geom_y)
                     ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_coordx, geom_x(1:nNodes), start = (/ j /), count = (/ nNodes /))
                     ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_coordy, geom_y(1:nNodes), start = (/ j /), count = (/ nNodes /))
                     j = j + nNodes
                  end if
               end do
               ierr = nf90_put_var(hisids%ncid, hisids%id_weirgeom_node_count, node_count, start = (/ 1 /), count = (/ nweirgen /))
            end if
         end if
      end if

      if (jahisdambreak > 0 .and. ndambreak > 0 .and. janormalhis > 0) then
         do i = 1,ndambreaksg
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_discharge,                    valdambreak(2,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_s1up,                         valdambreak(3,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_s1dn,                         valdambreak(4,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_head,                         valdambreak(5,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_au,                           valdambreak(6,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_normal_velocity,              valdambreak(7,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_cresth,                       valdambreak(8,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_crestw,                       valdambreak(9,i),  (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_water_level_jump,             valdambreak(10,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_breach_width_time_derivative, valdambreak(11,i), (/ i, hisids%idx_curtime /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_dambreak_cumulative_discharge,         valdambreak(12,i), (/ i, hisids%idx_curtime /))
         end do
      end if
      if (timon) call timstop ( handle_extra(62))
      !
      if (timon) call timstrt('unc_write_his sed data', handle_extra(66))
      if (jased>0 .and. stm_included .and. jahissed>0 .and. stmpar%lsedtot>0 .and. janormalhis > 0) then
         if (stmpar%morpar%moroutput%sbcuv) then
            call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            call realloc(toutputy, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            do l = 1, stmpar%lsedtot
               select case(stmpar%morpar%moroutput%transptype)
               case (0)
                  rhol = 1d0
               case (1)
                  rhol = stmpar%sedpar%cdryb(l)
               case (2)
                  rhol = stmpar%sedpar%rhosol(l)
               end select
               toutputy(:,l) = valobsT(:,IPNT_SBCY1+l-1)/rhol
               toutputx(:,l) = valobsT(:,IPNT_SBCX1+l-1)/rhol
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_sbcx, toutputx  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_sbcy, toutputy  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
         endif
         !
         if (stmpar%morpar%moroutput%sscuv) then
            call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            call realloc(toutputy, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            do l = 1, stmpar%lsedtot
               select case(stmpar%morpar%moroutput%transptype)
               case (0)
                  rhol = 1d0
               case (1)
                  rhol = stmpar%sedpar%cdryb(l)
               case (2)
                  rhol = stmpar%sedpar%rhosol(l)
               end select
               toutputy(:,l) = valobsT(:,IPNT_SSCY1+l-1)/rhol
               toutputx(:,l) = valobsT(:,IPNT_SSCX1+l-1)/rhol
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_sscx, toutputx  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_sscy, toutputy  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
         endif
         !
         if (stmpar%morpar%moroutput%sbwuv .and. jawave>0 .and. .not. flowWithoutWaves) then
            call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            call realloc(toutputy, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            do l = 1, stmpar%lsedtot
               select case(stmpar%morpar%moroutput%transptype)
               case (0)
                  rhol = 1d0
               case (1)
                  rhol = stmpar%sedpar%cdryb(l)
               case (2)
                  rhol = stmpar%sedpar%rhosol(l)
               end select
               toutputy(:,l) = valobsT(:,IPNT_SBWY1+l-1)/rhol
               toutputx(:,l) = valobsT(:,IPNT_SBWX1+l-1)/rhol
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_sbwx, toutputx  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_sbwy, toutputy  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
         endif
         !
         if (stmpar%morpar%moroutput%sswuv .and. jawave>0 .and. .not. flowWithoutWaves) then
            call realloc(toutputx, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            call realloc(toutputy, (/ntot, stmpar%lsedtot /), keepExisting=.false., fill = dmiss)
            do l = 1, stmpar%lsedtot
               select case(stmpar%morpar%moroutput%transptype)
               case (0)
                  rhol = 1d0
               case (1)
                  rhol = stmpar%sedpar%cdryb(l)
               case (2)
                  rhol = stmpar%sedpar%rhosol(l)
               end select
               toutputy(:,l) = valobsT(:,IPNT_SSWY1+l-1)/rhol
               toutputx(:,l) = valobsT(:,IPNT_SSWX1+l-1)/rhol
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_sswx, toutputx  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_sswy, toutputy  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedtot, 1 /))
         endif
         !
         !
         if (stmpar%morpar%moroutput%sourcesink .and. IVAL_SOUR1>0) then
            call realloc(toutputx, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
            call realloc(toutputy, (/ntot, stmpar%lsedsus /), keepExisting=.false., fill = dmiss)
            do l = 1, stmpar%lsedsus
               toutputx(:,l) = valobsT(:,IPNT_SOUR1+l-1)
               toutputy(:,l) = valobsT(:,IPNT_SINK1+l-1)
            end do
            ierr = nf90_put_var(hisids%ncid, hisids%id_sour, toutputx  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedsus, 1 /))
            ierr = nf90_put_var(hisids%ncid, hisids%id_sink, toutputy  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ ntot, stmpar%lsedsus, 1 /))
         endif
      endif

    if ( dad_included .and. janormalhis > 0) then  ! Output for dredging and dumping
       ierr = nf90_put_var(hisids%ncid, hisids%id_dredlink_dis, dadpar%link_sum  , start = (/ 1, 1, hisids%idx_curtime /), count = (/ dadpar%nalink, stmpar%lsedtot, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_dred_dis    , dadpar%totvoldred, start = (/ 1, hisids%idx_curtime /), count = (/ dadpar%nadred+dadpar%nasupl, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_dump_dis    , dadpar%totvoldump, start = (/ 1, hisids%idx_curtime /), count = (/ dadpar%nadump, 1 /))

       cof0 = 1d0 ; if( time_his > 0d0 ) cof0 = time_his
       ierr = nf90_put_var(hisids%ncid, hisids%id_dred_tfrac  , dadpar%tim_dredged/cof0  , start = (/ 1, hisids%idx_curtime /), count = (/ dadpar%nadred+dadpar%nasupl, 1 /))
       ierr = nf90_put_var(hisids%ncid, hisids%id_plough_tfrac, dadpar%tim_ploughed/cof0 , start = (/ 1, hisids%idx_curtime /), count = (/ dadpar%nadred+dadpar%nasupl, 1 /))
    endif
    if (timon) call timstop(handle_extra(66))

    if (timon) call timstrt('unc_write_his IDX data', handle_extra(67))
    if (janormalhis > 0) then
       do num = 1,MAX_IDX
          if ( num.eq.IDX_InternalTidesDissipation ) then
             if ( jaFrcInternalTides2D.eq.1 ) then
                ierr = nf90_put_var(hisids%ncid, hisids%id_voltot(num), 1d-12*voltot(num),  start=(/ hisids%idx_curtime /))
             end if
          else if ( num.eq.IDX_GravInput ) then
             if ( jatidep.eq.1 .or. jaselfal.gt.0 ) then
                ierr = nf90_put_var(hisids%ncid, hisids%id_voltot(num), 1d-12*voltot(num),  start=(/ hisids%idx_curtime /))
             end if
          else if ( num.eq.IDX_SALInput .or. num.eq.IDX_SALInput2 ) then
             if ( jaselfal.gt.0 ) then
                ierr = nf90_put_var(hisids%ncid, hisids%id_voltot(num), 1d-12*voltot(num),  start=(/ hisids%idx_curtime /))
             end if
          else
             ierr = nf90_put_var(hisids%ncid, hisids%id_voltot(num), voltot(num),  start=(/ hisids%idx_curtime /))
          end if
       enddo
    end if
    if (timon) call timstop(handle_extra(67))

    if( jahisgate > 0 .and. ngatesg+ngategen > 0 .and. janormalhis > 0) then
       ! todo: remove all do loops
       ! ngatesg ! Old-fashioned gates 'gateloweredgelevel'
       ! Actual discharge:
       !ierr = nf90_put_var(hisids%ncid, hisids%id_gatedisch, gatedisch(1:ngatesg+ngategen),  start=(/ 1, hisids%idx_curtime /), count = (/ ngatesg+ngategen, 1 /))

       !'pump_discharge_pumpA'
       !'pump_discharge_pumpB' (1:ntimes)

       !'pump_names' (1:npumps)
       !'pump_discharge'  (1:ntimes, 1:npumps)

       ! Door lower edge level
      ! ierr = nf90_put_var(hisids%ncid, hisids%id_zgate,     work ... (1:ngatesg+ngategen),  start=(/ 1, hisids%idx_curtime /), count = (/ ngatesg+ngategen, 1 /))

       ! hisids%id_gatesill: not for old style gates, they are just at bed level, so leave empty value in the file on columns 1:ngatesg
!       ierr = nf90_put_var(hisids%ncid, hisids%id_gatesill, gatesill(ngatesg+1:ngatesg+ngategen),  start=(/ ngatesg+1, hisids%idx_curtime /), count = (/ ngategen, 1 /))

!       ierr = nf90_put_var(hisids%ncid, hisids%id_zgate(num)    , zgate(num)    ,  start=(/ hisids%idx_curtime /))

!       do num=1,ngategen ! New-style gate, via generalstructure
!          igen=gate2cgen(num)
!          ipos = ! Just add new style gates at the back of the old style gates
!          ierr = nf90_put_var(hisids%ncid, hisids%id_gatedisch,      gatedisch(ipos),  start=(/ hisids%idx_curtime /))
!          ierr = nf90_put_var(hisids%ncid, hisids%id_zgate(num)    , zgate(num)    ,  start=(/ hisids%idx_curtime /))
    endif
        if( jahiscdam > 0 .and. ncdamsg + nweirgen > 0 .and. janormalhis > 0) then
           ! see gates
       do num = 1,ncdamsg
!          ierr = nf90_put_var(hisids%ncid, hisids%id_cdamdisch(num), cdamdisch(num),  start=(/ hisids%idx_curtime /))
!          ierr = nf90_put_var(hisids%ncid, hisids%id_zcdam(num)    , zcdam(num)    ,  start=(/ hisids%idx_curtime /))
       enddo
    endif

    if ( jacheckmonitor.eq.1 .and. janormalhis > 0) then
      ierr = nf90_put_var(hisids%ncid, hisids%id_checkmon, checkmonitor, start=(/ 1, hisids%idx_curtime /))

      ierr = nf90_put_var(hisids%ncid, hisids%id_num_timesteps, int(dnt), start=(/ hisids%idx_curtime /))
      ierr = nf90_put_var(hisids%ncid, hisids%id_comp_time, tim_get_wallclock(handle_steps), start=(/ hisids%idx_curtime /))
    end if

    if (unc_noforcedflush == 0) then
    ierr = nf90_sync(hisids%ncid) ! Flush file
    end if

    if (timon) call timstop (handle_extra(54))

end subroutine unc_write_his
