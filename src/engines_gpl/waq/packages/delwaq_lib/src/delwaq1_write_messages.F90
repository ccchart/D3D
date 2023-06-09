!!  Copyright (C)  Stichting Deltares, 2012-2023.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.


!>\file
!>                    delwaq1_write_messages

subroutine delwaq1_write_messages(errorcode)
    use m_dhopnf
    use m_delwaq1_data
    use m_dattim

    implicit none

    integer, intent(inout)                        :: errorcode

    write ( lunrep,'(//'' Messages presented including .lsp file:'')')
    write ( lunrep,'(  '' Number of WARNINGS            :'',I6)') iwar
    write ( lunrep,'( /'' Number of ERRORS during input :'',I6)') ierr
    write (   *   ,'(  ''  Number of WARNINGS            :'',I6)') iwar
    write (   *   ,'(  ''  Number of ERRORS during input :'',I6)') ierr
    write (   *   ,'(  '' '')')

    if ( ierr .eq. 0 ) then
        novec = min(novec,(nosss+nobnd-1))
        itota = 0
        itoti = 0
        itotc = 0
        call space  ( lunrep, .false., abuf   , ibuf   , chbuf  , &
                      itota , itoti  , itotc  )

        call dhopnf  ( lun(1) , lchar(1) , 1     , 1     , ioerr )
        write ( lun(1) )   in
        write ( lun(1) )   ii
        write ( lun(1) )   itota , itoti , itotc
        write ( lun(1) ) ( lun    (k) , k = 1,nolun  )
        write ( lun(1) ) ( lchar  (k) , k = 1,nolun  )
        write ( lun(1) ) ( filtype(k) , k = 1,nolun  )
    else
        write ( lunrep , '(  '' SIMULATION PROHIBITED !!!!!!!!'')' )
        call dhopnf  ( lun(1) , lchar(1) , 1     , 3     , ioerr )
        errorcode = 1
    endif

    call dattim(rundat)
    write (lunrep,'(2A)') ' Execution stop : ',rundat
    close ( lunrep )

end subroutine delwaq1_write_messages