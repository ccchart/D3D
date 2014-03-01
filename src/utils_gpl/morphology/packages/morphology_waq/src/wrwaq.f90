      module wrwaq
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2014.                                
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
!-------------------------------------------------------------------------------

!include preprocessing flags from autotools
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      implicit none
      
      contains

      function openwaqbinfile(filename) result (lun)
         character*(*), intent(in) :: filename !< Output filename.
         integer :: lun
         integer, external :: newunit
!
!           WARNING: WAQ input files must be written using form=binary
!                    instead of unformatted.
!                    Although it is not standard Fortran
!
         lun    = newunit()
#ifdef HAVE_FC_FORM_BINARY
         open  ( lun , file=filename , form = 'binary' , SHARED )
#else
! standardized way if binary is not available
         open  ( lun , file=filename , form = 'unformatted' , access='stream' )
#endif
      end function openwaqbinfile


      function openasciifile(filename) result (lun)
         character*(*), intent(in) :: filename
         integer :: lun
         integer, external :: newunit
!
!           NOTE: Opens a simple ASCII-file. Function is intended only
!                 to isolate newunit dependency.
!
         lun    = newunit()
#ifdef HAVE_FC_FORM_BINARY
         open  ( lun , file=filename , SHARED )
#else
         open  ( lun , file=filename , access='stream')
#endif
      end function openasciifile


!> Write ASCII or binary pointer file for WAQ.
      subroutine wrwaqpoi(ifrmto, noq, filename, ascii)
      implicit none
!
!           Global variables
!
      integer                , intent(in) :: noq      !< Nr. of linkages (pointers) between computational cells.
      integer, dimension(:,:), intent(in) :: ifrmto   !< Pointer table with all linkages.
                                                      !! ifrmto(1,:) = 'from'   cell number
                                                      !! ifrmto(2,:) = 'to'     cell number
                                                      !! ifrmto(3,:) = 'from-1' cell number
                                                      !! ifrmto(4,:) = 'to+1'   cell number
      logical                , intent(in) :: ascii    !< Produce ascii file or not (then binary).
      character(*)           , intent(in) :: filename !< Name for output pointer file.
!
!           Local variables
!
      integer :: i
      integer :: lunout
      integer :: q
!
!! executable statements -------------------------------------------------------
!
      if (ascii) then
         !
         ! ascii output
         !
         lunout = openasciifile(filename)
         do q = 1,noq
            write(lunout,'(4i10)') ( ifrmto(i,q), i=1,4 )
         enddo
         close(lunout)
      else
         !
         ! binary output
         !
         lunout = openwaqbinfile(filename)
         write(lunout) ( ( ifrmto(i,q), i=1,4 ), q=1,noq )
         close(lunout)
      endif
      end subroutine wrwaqpoi
!
!------------------------------------------------------------------------------


      !> Write (binary) from/to length file for DelWAQ.
      subroutine wrwaqlen(lenex, noq, filename, ascii)
      use precision
!
      implicit none
!
!           Global variables
!
      integer      , intent(in) :: noq           !< Nr. of linkages (pointers) between computational cells.
      real(hp)     , intent(in) :: lenex(2, noq) !< Dispersion half-lengths of computational cells, segment
                                                 !! centre to exchange point. (2 values: from/to direction)
      logical      , intent(in) :: ascii         !< Produce ascii file or not (then binary).
      character(*) , intent(in) :: filename      !< Output filename.
!
!           Local variables
!
      integer :: i
      integer :: lunout
      integer :: q
!
!! executable statements -------------------------------------------------------
!
      if (ascii) then
         !
         ! ascii output
         !
         lunout = openasciifile(filename)
         do q = 1,noq
            write(lunout,'(i10,2f18.8)') q, ( lenex(i,q), i=1,2 )
         enddo
         close(lunout)
      else
         !
         ! binary output
         !
         lunout = openwaqbinfile(filename)
         write(lunout) noq
         write(lunout) (( real(lenex(i,q),sp), i=1,2 ), q=1,noq )
         close(lunout)
      endif
      end subroutine wrwaqlen
!
!------------------------------------------------------------------------------



      !> Write (binary) exchange file(s) for DelWAQ: area and fluxes.
      subroutine wrwaqbin(itim, quant, nquant, filename, ascii, lunout)
      use precision
      implicit none
!
!           Global variables
!
      integer                 , intent(in)    :: itim     !< Time for new data block
      integer                 , intent(in)    :: nquant   !< Size of quant(ity) array.
      real(hp), dimension(:)  , intent(in)    :: quant    !< Quantity array to be written.
      logical                 , intent(in)    :: ascii    !< Produce ascii file or not (then binary).
      character(*)            , intent(in)    :: filename !< Output filename (only used if lunout not connected yet).
      integer                 , intent(inout) :: lunout   !< File pointer for output file. Used if already connected,
                                                          !! or set to new value for filename.
!
!           Local variables
!
      integer :: q
!
!! executable statements -------------------------------------------------------
!
      if (ascii) then
         !
         ! ascii output
         !
         if (lunout<0) then
            lunout = openasciifile(filename)
         endif
         write(lunout,'(a,i10)') 'Time = ', itim
         do q = 1,nquant
            write(lunout,'(i10,f18.8)') q, quant(q)
         enddo
         !close(lunout)
      else
         !
         ! binary output
         !
         if (lunout<0) then
            lunout = openwaqbinfile(filename)
         endif
         write(lunout) itim
         write(lunout) ( real(quant(q),sp), q=1,nquant )
         !close(lunout)
      endif
      end subroutine wrwaqbin
!
!------------------------------------------------------------------------------


      !> Write monitoring segments file for DelWAQ (each segment is a monitoring area).
      subroutine wrmonseg(noseg, filename)
      implicit none
!
!           Global variables
!
      integer                 , intent(in) :: noseg
      character(*)            , intent(in) :: filename !< Output filename.
!
!           Local variables
!
      integer :: lunout
      integer :: s
!
!! executable statements -------------------------------------------------------
!
      lunout = openasciifile(filename)
      write(lunout,'(i5)') noseg
      do s = 1,noseg
         write(lunout,'(a,i4,a,i5)') '''Segment ',s,''' 1 ',s
      enddo
      close(lunout)
      end subroutine wrmonseg
!
!------------------------------------------------------------------------------


      !> Write NROFSEGM.DAT file for DelWAQ.
      subroutine wr_nrofseg(noseg, filename)
      implicit none
!
!           Global variables
!
      integer                 , intent(in) :: noseg
      character(*)            , intent(in) :: filename !< Output filename.
!
!           Local variables
!
      integer :: lunout
!
!! executable statements -------------------------------------------------------
!
      lunout = openasciifile(filename)
      write(lunout,'(i12,a)') noseg,'   ; number of segments'
      close(lunout)
      end subroutine wr_nrofseg
!
!------------------------------------------------------------------------------


      !> Write NROFEXCH.DAT file for DelWAQ.
      subroutine wr_nrofexch(noq1, noq2, noq3, filename)
      implicit none
!
!           Global variables
!
      integer                 , intent(in) :: noq1
      integer                 , intent(in) :: noq2
      integer                 , intent(in) :: noq3
      character(*)            , intent(in) :: filename !< Output filename.
!
!           Local variables
!
      integer :: lunout
!
!! executable statements -------------------------------------------------------
!
      lunout = openasciifile(filename)
      write(lunout,'(3i12,a)') noq1, noq2, noq3, '   ; number of exchanges in three directions'
      close(lunout)
      end subroutine wr_nrofexch
!
!------------------------------------------------------------------------------

      end module wrwaq
