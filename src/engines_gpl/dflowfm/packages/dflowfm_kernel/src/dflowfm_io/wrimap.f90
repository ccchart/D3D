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

!> Toplevel map file writer.
!! Takes care of old/ugrid/tecplot format.
!! Does not take care of MapInterval checks, that is in flow_externaloutput().
subroutine wrimap(tim, imapdim, is_last_map_time)
    use m_flow
    use m_flowtimes
    use m_observations
    use unstruc_netcdf
    use unstruc_model
    use unstruc_files , only: defaultFilename
    use m_dad, only: dad_included
    use m_fm_update_crosssections, only: fm_update_mor_width_mean_bedlevel
    use m_flowgeom, only: ndx2d, ndxi
    use Timers

    implicit none
    double precision, intent(in   ) :: tim     !< Current time for writing a new map snapshot
    integer,          intent(in   ) :: imapdim !< Dimensional type of map file, see UNC_DIM_1D/2D/3D/ALL for valid options.
    logical,          intent(in   ) :: is_last_map_time !< Whether or not this is the last snapshot written to file. Used to trigger writing of some final run diagnostics (only for UGRID).

    ! locals
    type(t_unc_mapids), pointer :: mapidsptr
    integer            :: ierr
    integer            :: i
    integer            :: len
    integer, save      :: mtecfil = 0
    character(len=256) :: filnam
    logical            :: unitused
    double precision, save :: curtime_split = 0d0 ! Current time-partition that the file writer has open.
    integer            :: ndx1d, ndims
    integer            :: jabndnd

    ! Another time-partitioned file needs to start, reset iteration count (and file) (only for old map and tecplot output).
    if (ti_split > 0d0 .and. curtime_split /= time_split0) then
        mapids%id_tsp%idx_curtime = 0
        it_map       = 0
        it_map_tec   = 0
        curtime_split = time_split0
    end if

    if (md_mapformat == IFORMAT_UGRID) then
       if (imapdim == UNC_DIM_3D) then ! Separate 3D map file
          mapidsptr => mapids3d
       elseif (unc_check_dimension(imapdim, UNC_DIM_1D) .or. unc_check_dimension(imapdim, UNC_DIM_2D)) then ! Regular map file: 1D/2D, or 1D/2D+3D combined
          mapidsptr => mapids
       else
          write (msgbuf, '(a,i0)') 'wrimap error: invalid imapdim=', imapdim
          call err_flush()
          return
       end if
    else
       ! Old map output, only use regular map files.
       mapidsptr => mapids
    end if

     if ( md_mapformat.eq.IFORMAT_NETCDF .or. md_mapformat.eq.IFORMAT_NETCDF_AND_TECPLOT .or. md_mapformat == IFORMAT_UGRID) then   !   NetCDF output
       call timstrt('wrimap inq ncid', handle_extra(81))
       if (mapidsptr%ncid /= 0 .and. ((md_unc_conv == UNC_CONV_UGRID .and. mapidsptr%id_tsp%idx_curtime == 0) .or. (md_unc_conv == UNC_CONV_CFOLD .and. it_map == 0))) then
           ierr = unc_close(mapidsptr%ncid)
           mapidsptr%ncid = 0
       end if
       call timstop(handle_extra(81))

       call timstrt('wrimap inq ndims', handle_extra(80))
       if (mapidsptr%ncid/=0) then  ! reset stored ncid to zero if file not open
		  ierr = nf90_inquire(mapidsptr%ncid, ndims)
		  if (ierr/=0) mapidsptr%ncid = 0
       end if
       call timstop(handle_extra(80))

       call timstrt('wrimap unc_create', handle_extra(82))
       if (mapidsptr%ncid == 0) then
           if (ti_split > 0d0) then
               filnam = defaultFilename('map', timestamp=time_split0)
           else
               ! Separate 3D output only supported without timesplitting.
               if (imapdim == UNC_DIM_3D) then
                  filnam = defaultFilename('map3d')
               else
                  filnam = defaultFilename('map')
               end if
           end if
           ierr = unc_create(filnam , 0, mapidsptr%ncid)
           if (ierr /= nf90_noerr) then
               call mess(LEVEL_WARN, 'Could not create map file '''//trim(filnam)//'''.')
               mapidsptr%ncid = 0
           end if
       endif
       call timstop(handle_extra(82))

       if (mapidsptr%ncid .ne. 0) then
          if (md_unc_conv == UNC_CONV_UGRID) then
             ndx1d = ndxi - ndx2d
             if (ndx1d > 0 .and. stm_included) then
                if (stmpar%morpar%moroutput%blave) then
                   call fm_update_mor_width_mean_bedlevel()
                endif
             endif
             jabndnd = 0
             if (jamapbnd > 0) jabndnd = 1
             call unc_write_map_filepointer_ugrid(mapidsptr,tim,jabndnd, imapdim = imapdim)  ! wrimap

             if (is_last_map_time) then
                ierr = unc_put_run_diagnostics(mapidsptr%ncid, mapidsptr%diag_ids)
             end if

          else
             call unc_write_map_filepointer(mapidsptr%ncid,tim)  ! wrimap
          endif
       endif

       call timstrt('wrimap nf90_sync', handle_extra(83))
       if (unc_noforcedflush == 0) then
          ierr = nf90_sync(mapidsptr%ncid) ! Flush file
       end if
       call timstop(handle_extra(83))
    end if

    if (imapdim /= UNC_DIM_3D .and. (md_mapformat.eq.IFORMAT_TECPLOT .or. md_mapformat.eq.IFORMAT_NETCDF_AND_TECPLOT)) then      ! TecPlot output
       !if (mtecfil /= 0 .and. it_map_tec == 0) then
       !   call doclose(mtecfil)
       !end if

       !if (it_map_tec == 0) then
       !     if (ti_split > 0d0) then
       !         filnam = defaultFilename('tec', timestamp=time_split0)
       !     else
       !         filnam = defaultFilename('tec')
       !     end if
       !   call newfil(mtecfil, filnam)
       !endif

       !call tecplot_out(mtecfil, tim, it_map_tec==0)

!      write grid in Tecplot format only once
       if ( it_map_tec.eq.0 ) then
          filnam = defaultFilename('net.plt')
          call wrinet_tecplot(filnam)   ! do not call findcells
       end if

!      write solution in Tecplot format
       filnam = defaultFilename('map.plt', timestamp=tim)
       call wrimap_tecplot(filnam)

       it_map_tec = it_map_tec+1
    end if
 end subroutine wrimap
