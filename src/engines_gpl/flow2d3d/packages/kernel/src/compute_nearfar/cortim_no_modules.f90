subroutine cortim_no_modules(lun    ,no_modules,no_val)
!
!!--copyright-------------------------------------------------------------------
! Copyright (c) 2009, WL | Delft Hydraulics. All rights reserved.
!!--disclaimer------------------------------------------------------------------
! This code is part of the Delft3D software system. WL|Delft Hydraulics has
! developed c.q. manufactured this code to its best ability and according to the
! state of the art. Nevertheless, there is no express or implied warranty as to
! this software whether tangible or intangible. In particular, there is no
! express or implied warranty as to the fitness for a particular purpose of this
! software, whether tangible or intangible. The intellectual property rights
! related to this software code remain with WL|Delft Hydraulics at all times.
! For details on the licensing agreement, we refer to the Delft3D software
! license and any modifications to this license, if applicable. These documents
! are available upon request.
!!--version information---------------------------------------------------------
! $Author$
! $Date$
! $Revision$
!!--description-----------------------------------------------------------------
!
!    Function: Determine the number of cormix modules used
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
!
! Global variables
!
    integer                                                    , intent(in)  :: lun
    integer                                                    , intent(out) :: no_modules
    integer                                                    , intent(out) :: no_val
!
! Local variables
!
    integer                                     :: ix1
    integer                                     :: ix2
    integer                                     :: ix3
    integer                                     :: iocond
    real(fp)                                    :: rdummy
    character*256                               :: record
    logical                                     :: found
!
!! executable statements -------------------------------------------------------
!
    no_modules = 0
    no_val     = 0
    iocond     = 0

    do while (iocond == 0)
       read (lun,'(a256)',iostat = iocond) record
       call small(record,len(record))

       ix1 = index (record,'begin')
       ix2 = index (record,'mod')

       if (ix1 /= 0 .and. ix2 /= 0) then
          found = .false.
          do while (.not. found)
             read (lun,'(a256)',iostat = iocond) record
             call small(record,len(record))
             ix1 = index (record,'x        y       z')
             ix2 = index (record,'end')
             ix3 = index (record,'mod')
             if (ix1 /= 0) then
                found = .true.
                no_modules = no_modules + 1
                ix2 = 0
                ix3 = 0
                do while (ix2 == 0 .and. ix3 == 0)
                   read (lun,'(a256)',iostat = iocond) record
                   call small (record,len(record))
                   ix2 = index(record,'end')
                   ix3 = index(record,'mod')
                   read (record,*,iostat = iocond) rdum, rdum, rdum, rdum
                   if (iocond == 0) then
                      no_val = no_val + 1
                   endif
                enddo
                iocond = 0
             elseif (ix2 /= 0 .and. ix3 /=0) then
                found = .true.
             endif
          enddo
       endif
    enddo
end subroutine cortim_no_modules
