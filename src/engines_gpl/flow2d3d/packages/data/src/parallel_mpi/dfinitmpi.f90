subroutine dfinitmpi
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
!   Join parallel application within Delft3D-FLOW
!
!!--pseudo code and references--------------------------------------------------
!
!   enroll in MPI
!   get node number
!   determine total number of processes
!   determine whether this is a parallel run or not
!   define MPI constants for communication within Delft3D-FLOW
!   determine precision for type real in case of communication (single or double)
!
!
!!--declarations----------------------------------------------------------------
#ifdef DFMPI
    use mpi
#endif
    use precision
    use dfparall
    !
    implicit none
!
! Local variables
!
    integer                            :: ierr      ! error value of MPI call
    integer                            :: len
    character(128)                     :: msgstr    ! string to pass message
    character(128)                     :: rankstr
#ifdef DFMPI
    character(MPI_MAX_PROCESSOR_NAME)  :: host      ! hostname       for current MPI process
    character(MPI_MAX_PROCESSOR_NAME)  :: processor ! processor name for current MPI process
#else
    integer, parameter                 :: MPI_SUCCESS = 0
#endif
    logical                            :: mpi_is_initialized
    logical                            :: usempi
!
!! executable statements -------------------------------------------------------
!
    ! enroll in MPI
    !
    mpi_is_initialized = .false.
    ierr = MPI_SUCCESS
    call get_environment_variable('PMI_RANK', rankstr, len)
    usempi = (len > 0)

    if (usempi) then
       !
       ! test first if MPI was initialized (to avoid boxing MPI_init in DD runs)    
       !
#ifdef DFMPI
       call mpi_initialized( mpi_is_initialized, ierr )
#endif
       if ( ierr /= MPI_SUCCESS ) then
          write (msgstr,'(a,i5)') 'MPI produces some internal error in mpi_initialized - return code is ',ierr
          write (6,*) trim(msgstr)
          call cstop( 1, char(0) )
       endif
       !
       ! early return if MPI was already initialized before...
       !
       if (mpi_is_initialized) return

#ifdef DFMPI
       call mpi_init ( ierr )
#endif
       if ( ierr /= MPI_SUCCESS ) then
          write (msgstr,'(a,i5)') 'MPI produces some internal error in mpi_init - return code is ',ierr
          write (6,*) trim(msgstr)
          call cstop( 1, char(0) )
       endif
    endif
    !
    ! initialize common variables
    !
    inode = 0
    nproc = 1
    !
    ! get node number INODE
    !
    if (usempi) then
#ifdef DFMPI
       host      = 'unknown'
       processor = 'unknown'
       call mpi_comm_rank ( MPI_COMM_WORLD, inode, ierr )
       call util_getenv('HOSTNAME',host)
       call mpi_get_processor_name (processor,len,ierr)
       write (6,'(a,i3.3,4a)') 'MPI process number ', inode, ' has host ', trim(host), ' and is running on processor ', trim(processor)
#endif
    endif

    inode = inode + 1

    if (usempi) then
       if ( ierr /= MPI_SUCCESS ) then
          write (msgstr,'(a,i5,a,i3.3)') 'MPI produces some internal error - return code is ',ierr,' and node number is ',inode
          write (6,*) trim(msgstr)
          call cstop( 1, char(0) )
       endif
    !
    ! determine total number of processes
    !
#ifdef DFMPI
       call mpi_comm_size ( MPI_COMM_WORLD, nproc, ierr )
#endif
       if ( ierr /= MPI_SUCCESS ) then
          write (msgstr,'(a,i5,a,i3.3)') 'MPI produces some internal error - return code is ',ierr,' and node number is ',inode
          write (6,*) trim(msgstr)
          call cstop( 1, char(0) )
       endif
    endif
    !
    ! determine whether this is a parallel run or not
    !
    if ( nproc > 1 ) then
       parll = .true.
    else
       parll = .false.
    endif
    !
    ! define MPI constants for communication within Delft3D-FLOW
    !
#ifdef DFMPI
    dfint  = MPI_INTEGER
    dfreal = MPI_REAL
    dfdble = MPI_DOUBLE_PRECISION
    dfchar = MPI_CHARACTER
    dfmax  = MPI_MAX
    dfmin  = MPI_MIN
    dfsum  = MPI_SUM
    !
    ! determine precision for type real in case of communication (single or double)
    !
    if ( fp == hp ) then
       dfloat = dfdble
    else
       dfloat = dfreal
    endif
    !
    if ( prec == hp ) then
       dfprec = dfdble
    else
       dfprec = dfreal
    endif
#endif

end subroutine dfinitmpi
