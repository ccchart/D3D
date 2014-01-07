!!  Copyright(C) Stichting Deltares, 2012-2014.
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

!!  Note: The "part" engine is not yet Open Source, but still under
!!  development. This package serves as a temporary dummy interface for
!!  the references in the "waq" engine to the "part" engine.

      subroutine rdfnam ( lun    , ifnam    , fnam   , nfil   , iout   ,            &
     &                    ipri   , alone    )

      use precision

      implicit none

      integer  ( ip), intent(in   ) :: nfil
      integer  ( ip), intent(  out) :: lun(nfil)
      character( * ), intent(in   ) :: ifnam
      character( * ), intent(  out) :: fnam(nfil)
      integer  ( ip), intent(in   ) :: iout
      integer  ( ip), intent(in   ) :: ipri
      logical       , intent(in   ) :: alone
      
      lun = 1
      fnam = "a"
      
      return

      end subroutine
