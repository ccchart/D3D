!!  Copyright (C)  Stichting Deltares, 2012-2018.
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
!!  stichting deltares
!!  p.o. box 177
!!  2600 mh delft, the netherlands
!!
!!  all indications and logos of, and references to registered trademarks
!!  of stichting deltares remain the property of stichting deltares. all
!!  rights reserved.

      subroutine wq_processes_integrate_fluxes ( conc   , amass  , deriv  , volume , idt    , &
                                                 nosys  , notot  , noseg  , surfac )

!     Deltares Software Centre

!>\File
!>           Sets an explicit time step from DERIV.

      implicit none

!     Parameters          :
!     type     kind  function         name                      description

      integer   (4), intent(in   ) :: nosys                   !< number of transported substances
      integer   (4), intent(in   ) :: notot                   !< total number of substances
      integer   (4), intent(in   ) :: noseg                   !< number of computational volumes
      real      (4), intent(inout) :: conc  (notot ,noseg)    !< concentrations per substance per volume
      real      (8), intent(inout) :: amass (notot ,noseg)    !< masses per substance per volume
      real      (4), intent(inout) :: deriv (notot ,noseg)    !< derivatives per substance per volume
      real      (4), intent(inout) :: volume(noseg )          !< volumes of the segments
      integer   (4), intent(in   ) :: idt                     !< integration time step size
      real      (4), intent(in   ) :: surfac(noseg)           !< horizontal surface

      ! local declarations

      integer                      :: iseg                    !  segment loop counter
      integer                      :: i                       !  substance loop counter
      real                         :: v1                      !  segment volume
      real                         :: s1                      !  segment surface
      real(8)                      :: a                       !  segment mass
      real(8)                      :: ndt                     !  time step (real)

      ! loop accross the number of computational elements

      ndt = real(idt)
      do iseg=1,noseg
         ! active substances first
         v1 = volume(iseg)
         if ( v1.gt.1.0e-25 ) then
            do i=1,nosys
               a           = amass(i,iseg) + ndt*deriv(i,iseg)
               amass(i,iseg) = a
               conc (i,iseg) = a / v1
               deriv(i,iseg) = 0.0
            enddo
         endif
         ! then the inactive substances
         s1 = surfac(iseg)
         if(s1.gt.0.0) then
            do i=nosys+1,notot
               a             = amass(i,iseg) + ndt*deriv(i,iseg)
               amass(i,iseg) = a
               conc (i,iseg) = a / s1
               deriv(i,iseg) = 0.0
            enddo
         endif
      enddo
      return
      end
