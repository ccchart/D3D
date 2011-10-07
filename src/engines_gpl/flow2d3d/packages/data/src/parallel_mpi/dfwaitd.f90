subroutine dfwaitd ( field, work, worksize, ks, ke, request, tag, gdp )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
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
!   Updates field array of type single precision through exchanging halo values
!   between neighbouring subdomains
!
!!--pseudo code and references--------------------------------------------------
!
!   for all neighbouring subdomains do
!      get subdomain number, pointer and size
!      store data to be sent in array WORK
!      send array WORK
!
!   for all neighbouring subdomains do
!      get subdomain number, pointer and size
!      receive next array and store in WORK
!      store the received data
!
!   Jan Thorbecke
!   June 2009
!
!!--declarations----------------------------------------------------------------
    use precision
#if defined (DFMPI)
    use mpi
#endif
    use dfparall
    use globaldata
    !
    implicit none
    !
    type(globdat), target    :: gdp
!
! Global variables
!
    integer                                         , intent(in)    :: ke           ! last index in vertical direction
    integer                                         , intent(in)    :: ks           ! first index in vertical direction
    integer                                         , intent(in)    :: tag          ! unique tag
    integer                                                         :: request(4,2) ! MPI communication handle
    integer                                         , intent(in)    :: worksize     ! 
    real(hp), dimension(gdp%d%nmlb:gdp%d%nmub,ks:ke), intent(inout) :: field        ! real array for which halo values must
    real(hp), dimension(worksize,4,2)               , intent(inout) :: work         ! work array to store data to be sent to or received from neighbour be copied from neighbouring subdomains
!
! Local variables
!
    integer, dimension(:), pointer :: iblkad
    integer              , pointer :: lundia
    integer                        :: idom   ! subdomain number
    integer                        :: inb    ! neighbour counter
    integer                        :: istart ! pointer in array IBLKAD
    integer                        :: itag   ! message tag for sending and receiving
    integer                        :: j      ! loop counter
    integer                        :: k      ! loop counter in vertical direction
    integer                        :: ksiz   ! size in vertical direction (e.g. total number of sigma layers)
    integer                        :: nneigh ! number of neighbouring subdomains
    integer                        :: novlu  ! number of overlapping unknowns
    integer                        :: n
    integer                        :: m
    integer                        :: indxddb
    integer                        :: ierr     ! error value of MPI call
#if defined (DFMPI)
    integer                        :: istat(mpi_status_size) ! MPI status array
    integer                        :: mpistatus(mpi_status_size)
#endif
    character(80)                  :: msgstr   ! string to pass message
!
!! executable statements -------------------------------------------------------
!
    iblkad => gdp%gdparall%iblkad
    lundia => gdp%gdinout%lundia
    !
    if (.not.parll) return
    !
    nneigh = iblkad(1)
    ksiz = ke-ks+1
    !
    ! for all neighbouring subdomains do
    !
    do inb = 1, nneigh
       !
       ! get subdomain number, pointer and size
       !
       idom   = iblkad(3*inb-1)
       istart = iblkad(3*inb+1)
       novlu  = iblkad(istart)
       !
       ! wait for array WORK
       !
#if defined (DFMPI)
       ! waiting for the send call
       call mpi_wait(request(inb,1), mpistatus, ierr)
       if ( ierr /= MPI_SUCCESS ) then
          write (msgstr,'(a,i5,a,i3.3)') 'MPI wait produces some internal error - return code is ',ierr,' and node number is ',inode
          call prterr(lundia, 'U021', trim(msgstr))
          call d3stop(1, gdp)
       endif
       ! waiting for the recv call
       call mpi_wait(request(inb,2), mpistatus, ierr)
#endif
       !
       ! store the received data
       !
       do k = ks, ke
          do j = 1, novlu
             !
             ! original numbering for exchanges has to be modified to account for
             ! the addition of the ddbound limits
             ! This is due to the fact that the partitioning and numbering
             ! used in iblkad are done BEFORE the extension of the array sizes
             ! with ddbound
             !
             n                 = mod(iblkad(istart+novlu+j)-1,gdp%d%nmax) + 1
             m                 = ((iblkad(istart+novlu+j)-1)/gdp%d%nmax)+1
             indxddb           = (m-1+gdp%d%ddbound)*(gdp%d%nmax+2*gdp%d%ddbound) + n + gdp%d%ddbound
             field(indxddb, k) = work((k-ks)*novlu+j, inb, 2)
          enddo
       enddo
    enddo
end subroutine dfwaitd

