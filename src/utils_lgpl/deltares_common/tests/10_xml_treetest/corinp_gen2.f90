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
    integer       ,dimension(:)  , allocatable :: m_amb
    integer       ,dimension(:)  , allocatable :: n_amb
    integer       ,dimension(:)  , allocatable :: m_intake
    integer       ,dimension(:)  , allocatable :: n_intake
    integer       ,dimension(:)  , allocatable :: k_intake
    real(fp)      ,dimension(:)  , allocatable :: q_diff
    real(fp)      ,dimension(:)  , allocatable :: t0_diff
    real(fp)      ,dimension(:)  , allocatable :: s0_diff
    real(fp)      ,dimension(:)  , allocatable :: d0
    real(fp)      ,dimension(:)  , allocatable :: h0
    real(fp)      ,dimension(:)  , allocatable :: sigma0
    real(fp)      ,dimension(:)  , allocatable :: theta0
    character(256),dimension(:,:), allocatable :: basecase
!
! Global variables
!
    logical   :: error
    logical   :: exist
!
! Local variables
!
    integer                :: luntmp
    integer, external      :: newlun
    integer                :: i
    integer                :: idis
    integer                :: istat
    real(fp)               :: dummy
    character(1)           :: slash
    character(256)         :: filename
    character(300)         :: cdummy
    character(300)         :: errmsg
    type(tree_data), pointer :: cosumofile_ptr
    type(tree_data), pointer :: cosumoblock_ptr
    type(tree_data), pointer :: node_ptr
!
!! executable statements -------------------------------------------------------
!
       slash = '\'
    !
    ! Create Cosumo input tree
    !
    filename = "COSUMOsettings.xml"
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
    no_dis = size(cosumoblock_ptr%child_nodes)
    allocate(m_diff(no_dis), stat=istat); m_diff = -999
    allocate(n_diff(no_dis), stat=istat); n_diff = -999
    !m_diff         => gdp%gdnfl%m_diff
    !n_diff         => gdp%gdnfl%n_diff
    !m_amb          => gdp%gdnfl%m_amb
    !n_amb          => gdp%gdnfl%n_amb
    !m_intake       => gdp%gdnfl%m_intake
    !n_intake       => gdp%gdnfl%n_intake
    !k_intake       => gdp%gdnfl%k_intake
    !q_diff         => gdp%gdnfl%q_diff
    !t0_diff        => gdp%gdnfl%t0_diff
    !s0_diff        => gdp%gdnfl%s0_diff
    !d0             => gdp%gdnfl%d0
    !h0             => gdp%gdnfl%h0
    !sigma0         => gdp%gdnfl%sigma0
    !theta0         => gdp%gdnfl%theta0
    !basecase       => gdp%gdnfl%basecase
    do i=1, no_dis
       node_ptr => cosumoblock_ptr%child_nodes(i)%node_ptr
       call prop_get(node_ptr, 'data/Mdiff', m_diff(i))
       call prop_get(node_ptr, 'data/Ndiff', n_diff(i))
    enddo

    do i=1, no_dis
       write(*,'(a,i0)') "Diffusor: ", i
       write(*,'(a,i0,a,i0,a)') "(Mdiff,Ndiff): (", m_diff(i), ",", n_diff(i), ")"
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
    ! Reading of the corinp.dat file by FLOW
    ! Parallel: should this be done by all partitions or just by one?
    ! Is concurrent file access possible?
    !
    !return
    !luntmp = newlun(gdp)
    !open (luntmp,file='corinp.dat')
    !!
    !! Read dummy line
    !!
    !!
    !! For each diffuser
    !!
    !do idis = 1, no_dis
    !   !
    !   ! Read position diffusor
    !   !
    !   read (luntmp,*) m_diff(idis)
    !   read (luntmp,*) n_diff(idis)
    !   !
    !   ! Read position ambient conditions
    !   !
    !   read (luntmp,*) m_amb(idis)
    !   read (luntmp,*) n_amb(idis)
    !   !
    !   ! Read intake location
    !   !
    !   read (luntmp,*) m_intake(idis)
    !   read (luntmp,*) n_intake(idis)
    !   read (luntmp,*) k_intake(idis)
    !   !
    !   ! Read discharge characteristics
    !   !
    !   read (luntmp,*) q_diff(idis)
    !   read (luntmp,*) t0_diff(idis)
    !   read (luntmp,*) s0_diff(idis)
    !   !
    !   ! Read remainder of cormix general input
    !   !
    !   read (luntmp,*) d0(idis)
    !   read (luntmp,*) h0(idis)
    !   read (luntmp,*) theta0(idis)
    !   read (luntmp,*) sigma0(idis)
    !   read (luntmp,'(a)') basecase(idis,1)
    !   if (trim(basecase(idis,1)) == "rundir") then
    !      call getcwd(cdummy)
    !      basecase(idis,1) = trim(cdummy)//slash
    !      write(*,'(3a)') "corinp_gen2: 'rundir' substituted by '", trim(basecase(idis,1)), "'"
    !   endif
    !enddo
    !!
    !! Close the general cormix input file
    !!
    !close (luntmp)
    !
end program corinp_gen2
