!!  Copyright (C)  Stichting Deltares, 2012-2019.
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
!!

module md5_checksum
    use iso_c_binding
    implicit none
    private
    public :: md5file, md5intarr, checksum2hex
    interface
        subroutine md5_init() bind(C, name = 'md5_init')
        end subroutine md5_init

        subroutine md5_update( chunk, length ) bind(C, name = 'md5_update')
            import c_long
            character(len=1), dimension(*), intent(in) :: chunk
            integer(kind=c_long), value                :: length
        end subroutine md5_update

        subroutine md5_final( result ) bind(C, name = 'md5_final')
            character(len=1), dimension(*) :: result
        end subroutine md5_final
    end interface
contains

subroutine md5file( filename, checksum, success )


!     Deltares Software Centre

!>\file
!>    Determine the MD5 checksum for an entire file
!
!     Note: the C source was retrieved from:
!     https://openwall.info/wiki/people/solar/software/public-domain-source-code/md5

    character(len=*), intent(in)  :: filename        !< Name of the file to be examined
    character(len=*), intent(out) :: checksum        !< MD5 checksum (14-bytes string)
    logical, intent(out)          :: success         !< Whether the procedure succeeded or not

    character(len=14)             :: result
    character(len=2048)           :: chunk
    integer                       :: filesize
    integer(kind=c_long)          :: length
    integer                       :: no_chunks
    integer                       :: i
    integer                       :: lun
    integer                       :: ierr

    success = .true.

    open( newunit = lun, file = filename, access = 'stream', status = 'old', iostat = ierr )

    if ( ierr /= 0 ) then
        success = .false.
        return
    endif

    inquire( lun, size = filesize )

    call md5_init

    no_chunks = (filesize + len(chunk) - 1) / len(chunk)
    do i = 1,no_chunks
        if ( i == no_chunks ) then
            length = mod( filesize, len(chunk) )
        else
            length = len(chunk)
        endif

        read( lun, iostat = ierr ) chunk(1:length)

        if ( ierr /= 0 ) then
            success = .false.
            exit
        endif

        call md5_update( chunk, length )
    enddo

    call md5_final( result )
    checksum = result

    close(lun)

end subroutine md5file

subroutine md5intarr( intarr, checksum, success )


!     Deltares Software Centre

!>\file
!>    Determine the MD5 checksum for one array of integers
!
!     Note: the C source was retrieved from:
!     https://openwall.info/wiki/people/solar/software/public-domain-source-code/md5

    integer,          intent(in)  :: intarr(:)       !< integer array to be examined
    character(len=*), intent(out) :: checksum        !< MD5 checksum (14-bytes string)
    logical, intent(out)          :: success         !< Whether the procedure succeeded or not

    character(len=14)             :: result
    character(len=2048)           :: chunk
    integer                       :: filesize
    integer(kind=c_long)          :: length
    integer                       :: no_chunks
    integer                       :: i, ii, j, k, m
    character(len=sizeof(i))      :: datablock

    success = .true.

    filesize = size(intarr) * sizeof(i)

    call md5_init

    no_chunks = (filesize + len(chunk) - 1) / len(chunk)
    ii = 0
    do i = 1,no_chunks
        if ( i == no_chunks ) then
            length = mod( filesize, len(chunk) )
        else
            length = len(chunk)
        endif

        do j = 1, length / sizeof(i)
           dataBlock(1:sizeof(i)) = transfer(intarr(ii+j), repeat(' ', sizeof(i)))

           do k = 1, sizeof(i)
               m = sizeof(i)*j + k - sizeof(i)
               chunk(m:m) = dataBlock(k:k)
           enddo
        end do
        ii = ii + length / sizeof(i)

        call md5_update( chunk, length )
    enddo

    call md5_final( result )
    checksum = result

end subroutine md5intarr

subroutine checksum2hex(checksum, s)
   character(len=*), intent(in) :: checksum
   character(len=*), intent(out) :: s
   
   character(len=2) :: temp

   integer :: i

   if (len(checksum) >= 14 .and. len(s) >= 28) then
      do i = 1, 14
         write(temp, '(Z2.2)') ichar(checksum(i:i))
         s(2*i-1:2*i) = temp
      end do
   else
      do i = 1, len(s)
         s(i:i) = '*'
      end do
   end if

end subroutine checksum2hex

end module md5_checksum
