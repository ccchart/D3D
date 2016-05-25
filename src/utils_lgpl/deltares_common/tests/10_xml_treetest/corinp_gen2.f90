program corinp_gen2
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
    integer                                :: no_dis
    integer       ,dimension(:)  , allocatable :: m_diff
    integer       ,dimension(:)  , allocatable :: n_diff
    integer       ,dimension(:,:), allocatable :: mn_amb
    real(sp)      ,dimension(:)  , allocatable :: q_diff
    real(hp)      ,dimension(:,:), allocatable :: s0_diff
    character(256),dimension(:)  , allocatable :: linkDir
!
! Global variables
!
    logical   :: error
    logical   :: exist
    logical   :: test
!
! Local variables
!
    integer                :: luntmp
    integer, external      :: newlun
    integer                :: i
    integer                :: idis
    integer                :: istat
    real(fp)               :: version
    character(256)         :: filename
    character(300)         :: errmsg
    type(tree_data), pointer :: cosumofile_ptr
    type(tree_data), pointer :: cosumoblock_ptr
    type(tree_data), pointer :: node_ptr
!
!! executable statements -------------------------------------------------------
!
    !
    ! Create Cosumo input tree
    !
    filename = "COSUMOsettings.xml"
    write(*,*) "Reading XML file '", trim(filename), "'"
    call tree_create( 'COSUMO Input', cosumofile_ptr )
    call tree_put_data( cosumofile_ptr, transfer(trim(filename),node_value), 'STRING:XML' )
    !
    ! Put file in input tree
    !
    call prop_file('xml',trim(filename),cosumofile_ptr,istat)
    if (istat /= 0) then
       select case (istat)
       case(1)
          errmsg = "FILE_NOT_FOUND" // trim(filename)
          write(*,*) trim(errmsg)
       case(3)
          errmsg = "PREMATURE_EOF" // trim(filename)
          write(*,*) trim(errmsg)
       case default
          errmsg = "FILE_READ_ERROR" // trim(filename)
          write(*,*) trim(errmsg)
       endselect
       error = .true.
       stop
    endif
    call tree_get_node_by_name( cosumofile_ptr, 'COSUMO', cosumoblock_ptr )
    if (.not.associated(cosumoblock_ptr)) then
       error = .true.
       stop
    endif
    test    = .false.
    version = -999.9
    no_dis  = 0
    do i=1, size(cosumoblock_ptr%child_nodes)
       if (tree_get_name(cosumoblock_ptr%child_nodes(i)%node_ptr) == "settings") then
          no_dis = no_dis + 1
       endif
    enddo
    allocate(m_diff(no_dis), stat=istat); m_diff = -999
    allocate(n_diff(no_dis), stat=istat); n_diff = -999
    allocate(mn_amb(no_dis,2), stat=istat); mn_amb = -999
    allocate(q_diff(no_dis), stat=istat); q_diff = -999.9
    allocate(s0_diff(no_dis,3), stat=istat); s0_diff = -999.9
    allocate(linkDir(no_dis), stat=istat); linkDir = ' '
    
    call prop_get(cosumofile_ptr, 'COSUMO/test', test)
    call prop_get(cosumofile_ptr, 'COSUMO/version', version)
   
    idis = 0
    do i=1, size(cosumoblock_ptr%child_nodes)
       node_ptr => cosumoblock_ptr%child_nodes(i)%node_ptr
       if (tree_get_name(node_ptr) /= "settings") cycle
       idis = idis + 1
       call prop_get(node_ptr, 'data/Mdiff', m_diff(idis))
       call prop_get(node_ptr, 'data/Ndiff', n_diff(idis))
       call prop_get(node_ptr, 'data/mnambient', mn_amb(idis,:), 2)
       call prop_get(node_ptr, 'data/M3S', q_diff(idis))
       call prop_get(node_ptr, 'data/s0', s0_diff(idis,:), 3)
       call prop_get(node_ptr, 'data/linkInpDir', linkDir(idis))
    enddo

    write(*,*) "Test: ", test
    write(*,'(a,f10.5)') "Version: ", version
    do i=1, no_dis
       write(*,'(a,i0)') "Diffusor: ", i
       write(*,'(a,i0,a,i0,a)') "(Mdiff,Ndiff): (", m_diff(i), ",", n_diff(i), ")"
       write(*,'(a,i0,a,i0,a)') "(Mamb ,Namb ): (", mn_amb(i,1), ",", mn_amb(i,2), ")"
       write(*,'(a,f10.5)') "Q: ", q_diff(i)
       write(*,'(a,3f10.5)') "s0: ", s0_diff(i,:)
       write(*,'(a,a)') "linkdir: ", trim(linkdir(i))
    enddo

    luntmp = newunit()
    filename = "COSUMOsettings_output.xml"
    inquire(file=trim(filename),exist=exist)
    if (exist) then
       open (luntmp, file=trim(filename), iostat=istat, status='old')
       close (luntmp, status='delete')
    endif
    open (luntmp,file="COSUMOsettings_output.xml",iostat=istat,status='new')
    call prop_write_xmlfile(luntmp, cosumofile_ptr, 0, istat)
    close (luntmp)
    !
    call tree_destroy(cosumofile_ptr)
end program corinp_gen2
