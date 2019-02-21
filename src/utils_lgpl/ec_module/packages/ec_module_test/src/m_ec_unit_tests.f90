!  m_ec_unit_tests.f90
!
!  FUNCTIONS:
!  TestEcGetTimesteps - Unit tests for EcGetTimeSteps
!

!****************************************************************************
!
!  MODULE: m_ec_module_test
!
!  PURPOSE:  Unit test(s) for EC-module
!
!****************************************************************************

module m_ec_unit_tests
   use m_ec_message
   use m_ec_support
   use precision
   implicit none

   private

   public :: TestEcGetTimesteps

   contains

   subroutine TestEcGetTimesteps(success, errMessage)
      use m_ec_parameters
      logical,          intent(out) :: success     !< all tests are successful or not
      character(len=*), intent(out) :: errMessage  !< error message in case an error failed

      character(len=MaxNameLen) :: rec1
      character(len=MaxNameLen) :: rec2
      real(kind=hp) :: time_steps
      real(kind=hp), parameter :: time_step_expected1 = 53736.0_hp
      real(kind=hp), parameter :: time_step_expected2 = 0.0_hp

      rec1 = 'TIME = 0 hours since 2006-01-01 00:00:00 +00:00'
      rec2 = 'TIME = 0 hour since 2006-01-01 00:00:00 +00:00'

      errMessage = ' '
      call clearECMessage()

      !
      ! test 1: with conversion
      !
      success = ecGetTimesteps(rec1, time_steps, .true.)
      if (success) then
         success = comparereal(time_steps, time_step_expected1, 1d-10) == 0
      endif
      if (.not. success) then
         errMessage = 'error finding time in : ' // rec1
         return
      endif

      !
      ! test 2: without conversion
      !
      success = ecGetTimesteps(rec1, time_steps, .false.)
      if (success) then
         success = comparereal(time_steps, time_step_expected2, 1d-10) == 0
      endif
      if (.not. success) then
         errMessage = 'error finding time in : ' // rec1
         return
      endif

      !
      ! test 3: error handling
      !
      success = .not. ecGetTimesteps(rec2, time_steps)
      if (success) then
         success = (getEcMessage() == '|ec_support::ecGetTimesteps: can not find time step in: time = 0 hour since 2006-01-01 00:00:00 +00:00.|')
      endif
      if (.not. success) then
         errMessage = 'error handling record with unknown time unit'
         return
      endif

      !
      ! test 4: error handling
      !
      rec1 = ' '
      success = .not. ecGetTimesteps(rec1, time_steps)
      if (success) then
         success = (getEcMessage() == '|ec_support::ecGetTimesteps: Input string is empty.|')
      endif
      if (.not. success) then
         errMessage = 'error handling empty record'
         return
      endif

      call clearECMessage()
   end subroutine TestEcGetTimesteps
end module m_ec_unit_tests
