subroutine corinp_gen2(error, gdp)
!----- GPL ---------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2011-2016.
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation version 3.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  contact: delft3d.support@deltares.nl
!  Stichting Deltares
!  P.O. Box 177
!  2600 MH Delft, The Netherlands
!
!  All indications and logos of, and references to, "Delft3D" and "Deltares"
!  are registered trademarks of Stichting Deltares, and remain the property of
!  Stichting Deltares. All rights reserved.
!
!-------------------------------------------------------------------------------
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
!    Function: Reads input needed for the coupling of Corjet/Cortime/Cormix
!              with Delft3d-Flow
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
!
    use precision
    use properties
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    integer                      , pointer :: lundia
    integer                      , pointer :: no_dis
    real(fp)      ,dimension(:)  , pointer :: x_diff
    real(fp)      ,dimension(:)  , pointer :: y_diff
    real(fp)      ,dimension(:)  , pointer :: x_amb
    real(fp)      ,dimension(:)  , pointer :: y_amb
    real(fp)      ,dimension(:)  , pointer :: x_intake
    real(fp)      ,dimension(:)  , pointer :: y_intake
    real(fp)      ,dimension(:)  , pointer :: z_intake
    real(fp)      ,dimension(:)  , pointer :: q_diff
    real(fp)      ,dimension(:)  , pointer :: t0_diff
    real(fp)      ,dimension(:)  , pointer :: s0_diff
    real(fp)      ,dimension(:)  , pointer :: d0
    real(fp)      ,dimension(:)  , pointer :: h0
    real(fp)      ,dimension(:)  , pointer :: sigma0
    real(fp)      ,dimension(:)  , pointer :: theta0
    character(256),dimension(:)  , pointer :: base_path
    character(256),dimension(:,:), pointer :: basecase
    character(256)               , pointer :: filename
    type(tree_data)              , pointer :: cosumofile_ptr
!
! Global variables
!
    logical, intent(out)   :: error
!
! Local variables
!
    integer                :: luntmp
    integer, external      :: newlun
    integer                :: i
    integer                :: idis
    integer                :: istat
    integer                :: no_dis_read
    real(sp)               :: version
    character(1)           :: slash
    character(300)         :: cdummy
    character(300)         :: errmsg
    type(tree_data), pointer :: cosumoblock_ptr
    type(tree_data), pointer :: node_ptr
!
!! executable statements -------------------------------------------------------
!
    lundia         => gdp%gdinout%lundia
    no_dis         => gdp%gdnfl%no_dis
    x_diff         => gdp%gdnfl%x_diff
    y_diff         => gdp%gdnfl%y_diff
    x_amb          => gdp%gdnfl%x_amb
    y_amb          => gdp%gdnfl%y_amb
    x_intake       => gdp%gdnfl%x_intake
    y_intake       => gdp%gdnfl%y_intake
    z_intake       => gdp%gdnfl%z_intake
    q_diff         => gdp%gdnfl%q_diff
    t0_diff        => gdp%gdnfl%t0_diff
    s0_diff        => gdp%gdnfl%s0_diff
    d0             => gdp%gdnfl%d0
    h0             => gdp%gdnfl%h0
    sigma0         => gdp%gdnfl%sigma0
    theta0         => gdp%gdnfl%theta0
    base_path      => gdp%gdnfl%base_path
    basecase       => gdp%gdnfl%basecase
    filename       => gdp%gdnfl%infile
    !
    if (gdp%arch=='win32' .or. gdp%arch=='win64') then
       slash = '\'
    else
       slash = '/'
    endif
    !
    if (.not.associated(gdp%gdnfl%cosumofile_ptr)) then
       !
       ! Create Cosumo input tree
       !
       write(lundia,'(3a)') "Reading file '", trim(filename), "' ..."
       call tree_create( 'TransportFormula Input', cosumofile_ptr )
       call tree_put_data( cosumofile_ptr, transfer(trim(filename),node_value), 'STRING' )
       !
       ! Put file in input tree
       !
       call prop_file('xml',trim(filename),cosumofile_ptr,istat)
       if (istat /= 0) then
          select case (istat)
          case(1)
             errmsg = FILE_NOT_FOUND // trim(filename)
             call write_error(errmsg, unit=lundia)
          case(3)
             errmsg = PREMATURE_EOF // trim(filename)
             call write_error(errmsg, unit=lundia)
          case default
             errmsg = FILE_READ_ERROR // trim(filename)
             call write_error(errmsg, unit=lundia)
          endselect
          nullify(gdp%gdnfl%cosumofile_ptr)
          error = .true.
          return
       endif
       !
       ! Store the file data(-pointer) in GDP
       !
       gdp%gdnfl%cosumofile_ptr => cosumofile_ptr
    else
       !
       ! File already read: reuse stored pointer
       !
       cosumofile_ptr => gdp%gdnfl%cosumofile_ptr
    endif
    !
    call tree_get_node_by_name( cosumofile_ptr, 'COSUMO', cosumoblock_ptr )
    if (.not.associated(cosumoblock_ptr)) then
       write(lundia, '(a)') "ERROR: Tag '<COSUMO>' not found"
       error = .true.
       return
    endif
    version     = -999.9
    no_dis_read = 0
    do i=1, size(cosumoblock_ptr%child_nodes)
       if (tree_get_name(cosumoblock_ptr%child_nodes(i)%node_ptr) == "settings") then
          no_dis_read = no_dis_read + 1
       endif
    enddo
    if (no_dis_read /= no_dis) then
       write(lundia,'(a,i0)') "ERROR: Unexpected number of discharges read: ", no_dis_read
    endif
    call prop_get(cosumofile_ptr, 'COSUMO/Fileversion', version)
    if (comparereal(version, 0.3_sp)) then
       write(lundia,'(a,f5.2)') "ERROR: Unexpected FileVersion number read: ", version
    endif
    !
    idis = 0
    do i=1, size(cosumoblock_ptr%child_nodes)
       node_ptr => cosumoblock_ptr%child_nodes(i)%node_ptr
       if (tree_get_name(node_ptr) /= "settings") cycle
       idis = idis + 1
       !
       ! Read position diffusor
       !
       call prop_get(node_ptr, 'data/Xdiff', x_diff(idis))
       call prop_get(node_ptr, 'data/Ydiff', y_diff(idis))
       !
       ! Read position ambient conditions
       !
       call prop_get(node_ptr, 'data/Xambient', x_amb(idis))
       call prop_get(node_ptr, 'data/Yambient', y_amb(idis))
       !
       ! Read intake location
       !
       call prop_get(node_ptr, 'data/Xintake', x_intake(idis))
       call prop_get(node_ptr, 'data/Yintake', y_intake(idis))
       call prop_get(node_ptr, 'data/Zintake', z_intake(idis))
       !
       ! Read discharge characteristics
       !
       call prop_get(node_ptr, 'data/M3s', q_diff(idis))
       call prop_get(node_ptr, 'data/T0', t0_diff(idis))
       call prop_get(node_ptr, 'data/S0', s0_diff(idis))
       !
       ! Read remainder of cormix general input
       !
       call prop_get(node_ptr, 'data/D0', d0(idis))
       call prop_get(node_ptr, 'data/H0', h0(idis))
       call prop_get(node_ptr, 'data/Theta0', theta0(idis))
       call prop_get(node_ptr, 'data/Sigma0', sigma0(idis))
       call prop_get(node_ptr, 'comm/FF2NFdir', base_path(idis))
       call prop_get(node_ptr, 'comm/FFrundir', basecase(idis,1))
       if (trim(basecase(idis,1)) == "rundir") then
          call getcwd(cdummy)
          basecase(idis,1) = trim(cdummy)//slash
          write(*,'(3a)') "corinp_gen2: 'rundir' substituted by '", trim(basecase(idis,1)), "'"
       endif
    enddo
    !
    ! Delete file info
    ! This will force rereading the Cosumo file at the next Cosumo calculation
    !
    call tree_destroy(gdp%gdnfl%cosumofile_ptr)
end subroutine corinp_gen2
