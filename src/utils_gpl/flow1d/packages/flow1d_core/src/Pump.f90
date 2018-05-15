module m_Pump
!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2018.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  This program is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
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

   use m_tables
   use MessageHandling

   implicit none

   public PrepareComputePump
   public ComputePump
   public dealloc

   interface dealloc
      module procedure deallocPump
   end interface dealloc

   type, public :: t_pump
      !> direction > 0: positive flow \n
      !! direction < 0: negative flow \n
      !! abs(direction): \n
      !! - 1: control only on suction side \n
      !! - 2: control only on delivery side \n
      !! - 3: control on suction side and on delivery side \n
      integer                                 :: direction
      integer                                 :: nrstages
      double precision, dimension(:), pointer :: capacity    => null()
      double precision, dimension(:), pointer :: ss_onlevel  => null()
      double precision, dimension(:), pointer :: ss_offlevel => null()
      double precision, dimension(:), pointer :: ds_onlevel  => null()
      double precision, dimension(:), pointer :: ds_offlevel => null()
      logical         , dimension(:), pointer :: ss_trigger  => null()
      logical         , dimension(:), pointer :: ds_trigger  => null()
      type(t_table), pointer                  :: reducfact   => null()

      ! Actual Parameters
      double precision                        :: computed_capacity
      double precision                        :: oldcapacity
      double precision                        :: capacitySetpoint
      logical                                 :: isControlled

      ! Output Parameters for Pump History File
      double precision                        :: ss_level
      double precision                        :: ds_level
      double precision                        :: pump_head
      integer                                 :: actual_stage
      logical                                 :: is_active
      double precision                        :: stage_capacity
      double precision                        :: reduction_factor
      double precision                        :: discharge

   end type

   private

contains

   subroutine deallocPump(pump)
      ! Modules

      implicit none
      ! Input/output parameters
      type(t_pump), pointer   :: pump

      ! Local variables

      ! Program code
      if (associated(pump)) then
         
         if (associated(pump%capacity))    deallocate(pump%capacity)
         if (associated(pump%ss_onlevel))  deallocate(pump%ss_onlevel)
         if (associated(pump%ss_offlevel)) deallocate(pump%ss_offlevel)
         if (associated(pump%ds_onlevel))  deallocate(pump%ds_onlevel)
         if (associated(pump%ds_offlevel)) deallocate(pump%ds_offlevel)
         if (associated(pump%ss_trigger))  deallocate(pump%ss_trigger)
         if (associated(pump%ds_trigger))  deallocate(pump%ds_trigger)
         call dealloc(pump%reducfact)
         
         pump%capacity    => null()
         pump%ss_onlevel  => null()
         pump%ss_offlevel => null()
         pump%ds_onlevel  => null()
         pump%ds_offlevel => null()
         pump%ss_trigger  => null()
         pump%ds_trigger  => null()
         pump%reducfact   => null()
         
         deallocate(pump)
         
         pump => null()
         
      endif
      
   end subroutine deallocPump

   subroutine PrepareComputePump(pump, s1m1, s1m2)
      !=======================================================================
      !                       Deltares
      !                One-Two Dimensional Modelling System
      !                           S O B E K
      !
      ! Subsystem:          Flow Module
      !
      ! Programmer:         J. van Beek
      !
      ! Module:             PLQHAP (PLuvius Q Discharge Advance Pump)
      !
      ! Module description: The discharge through an advanced pump with different capacities
      !     and control on sucktion and pressure side is calculated.
      !
      !     update information
      !     person                    date
      !
      !     23-03-2012: Redesigned by Jaap Zeekant based on algorithm specified by Thieu van Mierlo
      !
      !
      implicit none

      type(t_pump), pointer          :: pump
      double precision, intent(in)   :: s1m1 ! TODO: JN: Suction Side level
      double precision, intent(in)   :: s1m2 ! TODO: JN: Delivery Side Level
      !
      !
      ! Local variables
      !
      integer                        :: istage
      integer                        :: nstages
      double precision               :: qp
      logical                        :: ss_switch_on
      logical                        :: ss_switch_off
      logical                        :: ds_switch_on
      logical                        :: ds_switch_off
      double precision               :: ss_level  ! Suction Side level
      double precision               :: ds_level  ! Delivery Side Level

      nstages = pump%nrstages

      ! Check Suction Side Conditions
      if (pump%direction > 0) then
         !
         !              Pump Direction Positive
         !
         ss_level = s1m1
         ds_level = s1m2
      else
         !
         !              Pump Direction Negative
         !
         ss_level = s1m2
         ds_level = s1m1
      endif

      if (abs(pump%direction) == 1 .or. abs(pump%direction) == 3) then

        do istage = 1, nstages

            ss_switch_on = (ss_level > pump%ss_onlevel(istage))
            ss_switch_off = (ss_level < pump%ss_offlevel(istage))

            if (ss_switch_on .and. .not. ss_switch_off) then

              pump%ss_trigger(istage) = .true.

            elseif (.not. ss_switch_on .and. ss_switch_off) then

              pump%ss_trigger(istage) = .false.

            else
              ! Keep Old Value, So Do Nothing
            endif

        enddo

      endif

      ! Check Delivery Side Conditions
      if (abs(pump%direction) == 2 .or. abs(pump%direction) == 3) then

        do istage = 1, nstages

            ds_switch_on = (ds_level < pump%ds_onlevel(istage))
            ds_switch_off = (ds_level > pump%ds_offlevel(istage))

            if (ds_switch_on .and. .not. ds_switch_off) then

              pump%ds_trigger(istage) = .true.

            elseif (.not. ds_switch_on .and. ds_switch_off) then

              pump%ds_trigger(istage) = .false.

            else
              ! Keep Old Value, So Do Nothing
            endif

        enddo

      endif

      ! Give non-controlled side all freedom
      if (abs(pump%direction) == 1) then
        pump%ds_trigger = .true.
      elseif (abs(pump%direction) == 2) then
        pump%ss_trigger = .true.
      endif

      ! Find the active stage
      pump%actual_stage = 0
      do istage = nstages, 1, -1
        if (pump%ss_trigger(istage) .and. pump%ds_trigger(istage)) then
          pump%actual_stage = istage
          exit
        endif
      enddo

      ! Get Reduction Factor from Table
      pump%pump_head = ds_level - ss_level
      pump%reduction_factor = interpolate(pump%reducfact, pump%pump_head)

      ! Calculate Capacity
      if (pump%actual_stage == 0) then

        pump%is_active = .false.

        qp = 0.
        pump%stage_capacity = 0.0
        if (pump%iscontrolled) then
          pump%actual_stage = -1
        endif

      else

        pump%is_active = .true.

        if (pump%iscontrolled) then
          qp = pump%reduction_factor * pump%capacitySetpoint
          pump%actual_stage = -1
          pump%stage_capacity = pump%capacitySetpoint
        else
          qp = pump%reduction_factor * pump%capacity(pump%actual_stage)
          pump%stage_capacity = pump%capacity(pump%actual_stage)
        endif

      endif

      pump%ss_level = ss_level
      pump%ds_level = ds_level
      if (pump%direction > 0) then
        pump%discharge = qp
      else
        pump%discharge = -qp
        pump%stage_capacity = -pump%stage_capacity
      endif

      !
      !     Function return value is discharge Q for pump
      !
   end subroutine PrepareComputePump

   subroutine computePump(pump, fum, rum, um, qm, aum)
      ! modules

      ! Global variables
      type(t_pump), pointer, intent(in)            :: pump
      double precision, intent(out)                :: fum
      double precision, intent(out)                :: rum
      double precision, intent(out)                :: um
      double precision, intent(out)                :: qm
      double precision, intent(in)                 :: aum

      fum = 0.0
      rum = pump%discharge/max(aum, 1.0D-2)
      um = rum
      qm  = pump%discharge
   end subroutine computePump

end module m_Pump
