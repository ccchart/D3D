!----- LGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2011-2019.
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation version 2.1.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, see <http://www.gnu.org/licenses/>.
!
!  contact: delft3d.support@deltares.nl
!  Stichting Deltares
!  P.O. Box 177
!  2600 MH Delft, The Netherlands
!
!  All indications and logos of, and references to, "Delft3D" and "Deltares"
!  are registered trademarks of Stichting Deltares, and remain the property of
!  Stichting Deltares. All rights reserved.

!  $Id$
!  $HeadURL$

!> This module contains support methods for the EC-module.
!! @author adri.mourits@deltares.nl
!! @author robert.leander@deltares.nl
!! @author edwin.spee@deltares.nl
module m_ec_support
   use m_ec_message
   use m_ec_typedefs
   use m_alloc
   use string_module
   use m_ec_parameters
   use time_module

   implicit none

   private

   public :: ecTimeFrameRealHpTimestepsToModifiedJulianDate
   public :: ecTimeFrameRealHpTimestepsToDateTime
   public :: ecGetTimesteps
   public :: ecSupportOpenExistingFile
   public :: ecSupportOpenExistingFileGnu
   public :: ecSupportAddUniqueInt
   public :: ecSupportAddUniqueChar
   public :: ecSupportFindItemByQuantityLocation
   public :: ecSupportFindItem
   public :: ecSupportFindQuantity
   public :: ecSupportFindElementSet
   public :: ecSupportFindField
   public :: ecSupportFindConnection
   public :: ecSupportFindConverter
   public :: ecSupportFindFileReader
   public :: ecSupportFindFileReaderByFilename
   public :: ecSupportFindBCBlock
   public :: ecSupportFindNetCDF
   public :: ecSupportFindNetCDFByFilename
   public :: ecSupportFindBCFileByFilename
   public :: ecSupportNetcdfCheckError
   public :: ecSupportNCFindCFCoordinates
   public :: ecSupportTimestringToUnitAndRefdate
   public :: ecSupportTimeUnitConversionFactor
   public :: ecSupportThisTimeToMJD
   public :: ecSupportMJDToThisTime
   public :: ecSupportMJDToTimeIndex
   public :: ecSupportTimeIndexToMJD
   public :: ecSupportFindRelatedBCBlock

   interface ecTimeFrameRealHpTimestepsToModifiedJulianDate
      module procedure ecTimeFrameRealHpTimestepsToModifiedJulianDate
   end interface ecTimeFrameRealHpTimestepsToModifiedJulianDate

   interface ecNFC
      module procedure ecSupportNetcdfCheckError
   end interface ecNFC

   contains

      ! ==========================================================================

      !> Calculate a Julian Date from the number of timesteps in seconds since reference_date.
      function ecTimeFrameRealHpTimestepsToModifiedJulianDate(timeFramePtr, steps, mjd) result(success)
         logical                                 :: success      !< function status
         type(tEcTimeFrame), pointer             :: timeFramePtr !< intent(inout)
         real(hp),                   intent(in)  :: steps        !< number of time steps
         real(hp),                   intent(out) :: mjd          !< calculated Modified Julian Date
         !
         success = .true.
         !
         mjd = timeFramePtr%k_refdate + (steps / 60.0_hp / 60.0_hp / 24.0_hp)

      end function ecTimeFrameRealHpTimestepsToModifiedJulianDate

      ! =======================================================================

      !> Calculate a Gregorian date and hour-minutes-seconds integer since reference date
      function ecTimeFrameRealHpTimestepsToDateTime(timestamp_mjd, yyyymmdd, hhmmss) result(success)
      use mathconsts, only : daysec_hp
         logical                                 :: success             !< function status
         real(hp),                   intent(in)  :: timestamp_mjd       !< number of time steps
         integer,                    intent(out) :: yyyymmdd            !< calculated Gregorian date
         integer,                    intent(out) :: hhmmss              !< time of the day
         real(hp)                                :: dsec                !< fraction of seconds
         integer                                 :: hh, mm, ss          !< hours, minutes, seconds helper variables
         integer                                 :: iyear, imonth, iday !< year, month and day

         !! TODO very ugly to add offset_reduced_jd here
         call gregor(timestamp_mjd + offset_reduced_jd, iyear, imonth, iday, hh, mm, ss, dsec)

         yyyymmdd = 10000 * iyear + 100 * imonth + iday
         hhmmss = 10000 * hh + 100 * mm + ss

         success = .true.

      end function ecTimeFrameRealHpTimestepsToDateTime

      ! =======================================================================

      !> Read and convert the timesteps to seconds.
      !! Takes a string of format: TIME = 0 hours since 2006-01-01 00:00:00 +00:00
      !! or ... TIME (HRS)      6.0 20000101 6
      function ecGetTimesteps(rec, time_steps, convert) result(success)
         logical                           :: success    !< function status
         character(len=*),   intent(in)    :: rec        !< time information string
         real(kind=hp),      intent(out)   :: time_steps !< timesteps in seconds
         logical, optional,  intent(in)    :: convert    !< convert to mjd or leave unconverted

         integer                       :: unit
         real(kind=hp)                 :: ref_date, tzone, time_in
         integer                       :: i, posTimeUnit !< position in string of 'since'
         character(len=:), allocatable :: time_string
         character(len=*), parameter   :: time_units(4) = ['seconds', 'minutes', 'hours  ', 'days   ']
         !
         success = .false.
         !
         if (len_trim(rec) == 0) then
            call setECMessage("ec_support::ecGetTimesteps: Input string is empty.")
            return
         end if

         time_string = trim(rec)
         call str_lower(time_string)

         do i = 1, size(time_units)
            posTimeUnit = index(time_string, trim(time_units(i)))
            if (posTimeUnit > 0) exit
         end do

         if (posTimeUnit > 0) then
            if (.not.ecSupportTimestringToUnitAndRefdate(time_string(posTimeUnit:), unit, ref_date, tzone=tzone)) then
               call setECMessage("ec_support::ecGetTimesteps: can not convert time unit string: " // time_string // ".")
            endif
         endif

         if (posTimeUnit > 0) then
            read(time_string(index(time_string, '=')+1 : posTimeUnit-1), *) time_in
         else if ( .not. ecSupportTimestringArcInfo(time_string, time_steps=time_in)) then
            call setECMessage("ec_support::ecGetTimesteps: can not find time step in: " // time_string // ".")
            return
         endif

         if (present(convert)) then
            if (.not. convert) then
                time_steps = time_in
                success = .true.
                return
            end if
         end if
         select case (unit)
             case (ec_second)
                time_steps = ref_date + tzone / 24.0_hp + time_in / 86400.0_hp
             case (ec_minute)
                time_steps = ref_date + tzone / 24.0_hp + time_in / 1440.0_hp
             case (ec_hour)
                time_steps = ref_date + tzone / 24.0_hp + time_in / 24.0_hp
             case (ec_day)
                time_steps = ref_date + tzone / 24.0_hp + time_in
             case default
                call setECMessage("ec_support::ecGetTimesteps: Unable to identify the time unit.")
                return
         end select
         success = .true.
      end function ecGetTimesteps

      ! ==========================================================================

      !> Attempt to open an file for reading that might already have been opened under another handle.
      !> Workaround for GNU Fortran (which normally does not support multiple file openings of the same file)
      function ecSupportOpenExistingFileGnu(minp, filename) result(success)
         !
         use multi_file_io

         logical                                      :: success  !< function status
         integer(kind=8),               intent(inout) :: minp     !< IO unit number
         character(len=*), intent(in)                 :: filename !< relative path
         success = .false.

         minp = mf_open(filename)
         if (minp<=0) then
            call setECMessage("ec_support::ecSupportOpenExistingFileGnu: Opening "//trim(filename)//" failed.")
            return
         endif
         success = .true.
      end function ecSupportOpenExistingFileGnu

      ! ==========================================================================

      !> Attempt to open an existing file.
      function ecSupportOpenExistingFile(minp, filename) result(success)
         use netcdf
         !
         logical                         :: success  !< function status
         integer,          intent(out)   :: minp     !< IO unit number
         character(len=*), intent(in)    :: filename !< relative path
         !
         integer :: ierror         !< netcdf helper variable
         logical :: unitused       !< IO unit number already in use
         integer :: istat          !< status of file open operation
         integer :: chunkSizeHint  !< chunk size for NetCDF
         !
         success = .false.
         unitused = .false.
         ! Sanity checks.
         if (len_trim(filename) == 0) then
            call setECMessage("ec_support::ecSupportOpenExistingFile: Name is empty")
            return
         endif
         inquire(file = trim(filename), exist = success)
         if (.not. success) then
            call setECMessage("ec_support::ecSupportOpenExistingFile: File does not exist: ", trim(filename))
            return
         endif
         ! Special case: NetCDF.
         if (index(filename, '.nc') > 0) then
            chunkSizeHint = 4096   ! maximum allowed value
            ierror = nf90_open(trim(filename), NF90_NOWRITE, minp, chunkSizeHint)
            success = ecSupportNetcdfCheckError(ierror, "opening file", filename)
            return
         endif
         ! Locate an unused file unit.
         open(newunit=minp,file=trim(filename),iostat=istat)
         if (istat == 0) then
            success = .true.
         else
            minp = ec_undef_int
            call setECMessage("ec_support::ecSupportOpenExistingFile: opening file " // trim(filename) // " failed")
         endif
      end function ecSupportOpenExistingFile

      ! ==========================================================================

      !> Add an integer to a set of integers.
      function ecSupportAddUniqueInt(intArr, anInt) result(success)
         logical                        :: success!< function status
         integer, dimension(:), pointer :: intArr !< array containing unique integers (a set)
         integer, intent(in)            :: anInt  !< integer to be added to the set of integers

         integer                        :: i         !< loop counter
         integer                        :: istat     !< deallocate() status
         integer                        :: lenArr    !< lenght of intArr
         integer, dimension(:), pointer :: newIntArr !< larger version of intArr

         success = .false.
         newIntArr => null()
         i = 0
         istat = 1
         lenArr = 0

         if (.not. associated(intArr)) then
            call setECMessage("ec_support::ecSupportAddUniqueInt: Dummy argument pointer intArr is not associated.")
         else
            lenArr = size(intArr)
            do i=1, lenArr
               if (intArr(i) == anInt) then
                  ! This integer is already in intArr
                  success = .true.
                  return
               endif
            enddo
            ! This integer is not yet in intArr, so add it.
            allocate(newIntArr(lenArr+1), STAT = istat)
            if (istat /= 0 ) then
               call setECMessage("ec_support::ecSupportAddUniqueInt: Unable to allocate additional memory.")
            else
               do i=1, lenArr
                  newIntArr(i) = intArr(i) ! Copy existing integers.
               enddo
               newIntArr(lenArr+1) = anInt ! Add the new integer.
               deallocate(intArr, STAT = istat)
               if (istat /= 0 ) then
                  call setECMessage("WARNING: ec_support::ecSupportAddUniqueInt: Unable to deallocate old memory.")
               end if
               intArr => newIntArr
               success = .true.
            end if
         endif
      end function ecSupportAddUniqueInt

      ! ==========================================================================

      !> Add a char to a set of chars.
      function ecSupportAddUniqueChar(charArr, aChar) result(success)
         logical                                             :: success !< function status
         character(len=maxNameLen), dimension(:), pointer    :: charArr !< array containing unique chars (a set)
         character(len=*),                        intent(in) :: aChar   !< char to be added to the set of chars
         !
         integer                                          :: i          !< loop counter
         integer                                          :: istat      !< deallocate() status
         integer                                          :: lenArr     !< lenght of charArr
         character(len=maxNameLen), dimension(:), pointer :: newCharArr !< larger version of charArr
         !
         success = .false.
         newCharArr => null()
         i = 0
         istat = 1
         lenArr = 0
         !
         !if (.not. associated(charArr)) then
         !   reallocP(charArr, 1, STAT = istat)
         !   if (istat == 0) then
         !      charArr(1) = aChar
         !   else
         !      call setECMessage("ec_support::ecSupportAddUniqueChar: Unable to allocate additional memory.")
         !   end if
         !   return
         !end if
         !!
         !lenArr = size(charArr)
         !do i=1, lenArr
         !   if (charArr(i) == aChar) then
         !      ! This char is already in charArr
         !      success = .true.
         !      return
         !   endif
         !enddo
         !! This char is not yet in charArr, so add it.
         !reallocP(charArr, lenArr+1, STAT = istat)
         !if (istat == 0) then
         !   charArr(lenArr+1) = aChar
         !   success = .true.
         !else
         !   call setECMessage("ec_support::ecSupportAddUniqueChar: Unable to allocate additional memory.")
         !end if
      end function ecSupportAddUniqueChar

      ! =======================================================================
      ! Find methods
      ! =======================================================================

      !> Retrieve the pointer to the Quantity with id == quantityId.
      function ecSupportFindQuantity(instancePtr, quantityId) result(quantityPtr)
         type(tEcQuantity), pointer            :: quantityPtr !< Quantity corresponding to quantityId
         type(tEcInstance), pointer            :: instancePtr !< intent(in)
         integer,                   intent(in) :: quantityId  !< unique Quantity id
         !
         integer :: i !< loop counter
         !
         quantityPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nQuantities
               if (instancePtr%ecQuantitiesPtr(i)%ptr%id == quantityId) then
                  quantityPtr => instancePtr%ecQuantitiesPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindQuantity: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindQuantity

      ! =======================================================================
      !> Retrieve the pointer to a NetCDF object by filename
      function ecSupportFindNetCDFByFilename(instancePtr, ncfilename) result(netCDFPtr)
         type(tEcNetCDF), pointer            :: netCDFPtr   !< Quantity corresponding to quantityId
         type(tEcInstance), pointer          :: instancePtr !< intent(in)
         character(len=*),  intent(in)       :: ncfilename  !< netCDF filename
         !
         integer :: i !< loop counter
         !
         netCDFPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nNetCDFs
               if (instancePtr%ecNetCDFsPtr(i)%ptr%ncfilename == ncfilename) then
                  netCDFPtr => instancePtr%ecNetCDFsPtr(i)%ptr
               end if
            end do
         end if
      end function ecSupportFindNetCDFByFilename

      ! =======================================================================
      !> Retrieve the pointer to a BC-File object by filename
      function ecSupportFindBCFileByFilename(instancePtr, bcfilename) result(BCFilePtr)
         type(tEcBCFile), pointer            :: BCFilePtr   !< Quantity corresponding to quantityId
         type(tEcInstance), pointer          :: instancePtr !< intent(in)
         character(len=*),  intent(in)       :: bcfilename  !< BC filename
         !
         integer :: i !< loop counter
         !
         BCFilePtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nBCFiles
               if (instancePtr%ecBCFilesPtr(i)%ptr%bcfilename == bcfilename) then
                  BCFilePtr => instancePtr%ecBCFilesPtr(i)%ptr
               end if
            end do
         end if
      end function ecSupportFindBCFileByFilename
      ! =======================================================================

      !> Retrieve the pointer to the ElementSet with id == elementSetId.
      function ecSupportFindElementSet(instancePtr, elementSetId) result(elementSetPtr)
         type(tEcElementSet), pointer            :: elementSetPtr !< ElementSet corresponding to elementSetId
         type(tEcInstance),   pointer            :: instancePtr   !< intent(in)
         integer,                     intent(in) :: elementSetId  !< unique ElementSet id
         !
         integer :: i !< loop counter
         !
         elementSetPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nElementSets
               if (instancePtr%ecElementSetsPtr(i)%ptr%id == elementSetId) then
                  elementSetPtr => instancePtr%ecElementSetsPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindElementSet: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindElementSet

      ! =======================================================================

      !> Retrieve the pointer to the Field with id == fieldId.
      function ecSupportFindField(instancePtr, fieldId) result(fieldPtr)
         type(tEcField),    pointer            :: fieldPtr    !< Field corresponding to fieldId
         type(tEcInstance), pointer            :: instancePtr !< intent(in)
         integer,                   intent(in) :: fieldId     !< unique Field id
         !
         integer :: i !< loop counter
         !
         fieldPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nFields
               if (instancePtr%ecFieldsPtr(i)%ptr%id == fieldId) then
                  fieldPtr => instancePtr%ecFieldsPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindField: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindField

subroutine ecInstanceListSourceItems(instancePtr,dev)
         implicit none
         ! List source items by quantity and location
         type(tEcInstance), pointer :: instancePtr           !< EC instance
         integer, intent(in)        :: dev
         type(tEcItem), pointer     :: sourceItemPtr
         integer  :: ii
         do ii=1, instancePtr%nItems
            sourceItemPtr => instancePtr%ecItemsPtr(ii)%ptr
            if (sourceItemPtr%role == itemType_source) then
                     write(dev,'(a,i5.5)') 'Source Item ',sourceItemPtr%id
                     write(dev,'(a)')      '  Quantity = '//trim(sourceItemPtr%quantityPtr%name)
                     write(dev,'(a)')      '  Location = '//trim(sourceItemPtr%elementsetPtr%name)
                     write(dev,*) ''
            endif
         enddo
end subroutine ecInstanceListSourceItems

      !! =======================================================================
      !!> Retrieve the item ID given a quantitystring and locationstring
      !!> Loop over all items in the EC instance qualified as 'source'
      !function ecSupportFindItemByQuantityLocation(instancePtr, quantityname, locationname ) result(itemID)
      !   type(tEcInstance), pointer            :: instancePtr    !< EC-instance
      !   character(len=*), intent(in)         :: quantityname   !< Desired quantity
      !   character(len=*), intent(in)         :: locationname   !< Desired location
      !   integer                               :: itemID         !< returned item ID
      !   integer :: i                                            !< loop counter over items
      !   type (tEcItem), pointer               :: itemPtr
      !
      !   itemID = -1
      !   if (associated(instancePtr)) then
      !      do i=1, instancePtr%nItems
      !         itemPtr => instancePtr%ecItemsPtr(i)%ptr
      !         if (itemPtr%role==itemType_source) then
      !            if (itemPtr%quantityPtr%name==quantityname .and. itemPtr%elementSetPtr%name==locationname) then
      !               itemID = itemPtr%id
      !               exit
      !            end if
      !         end if
      !      end do
      !   end if
      !end function ecSupportFindItemByQuantityLocation

      ! =======================================================================
      !> Retrieve the item ID given a quantitystring and locationstring
      !> Use the fact that each filereader is associated with ONE location, but possibly MULTIPLE quantities
      !> i.e., select filereader first and check its items.
      function ecSupportFindItemByQuantityLocation(instancePtr, locationname, quantityname, isLateral) result(itemID)
         type(tEcInstance), pointer             :: instancePtr    !< EC-instance
         character(len=*), intent(in)           :: quantityname   !< Desired quantity
         character(len=*), intent(in)           :: locationname   !< Desired location
         logical         , intent(in), optional :: isLateral      !< searching for lateral discharge?
         integer                                :: itemID         !< returned item ID
         integer                                :: i, j           !< loop counter over filereader, items
         logical                                :: found          !< item found?
         type (tEcItem), pointer                :: itemPtr
         type (tEcFileReader), pointer          :: fileReaderPtr
         character(len=:), allocatable          :: quantity_requested, location_requested
         character(len=:), allocatable          :: quantity_supplied, location_supplied
!        character(len=maxNameLen)         :: quantityname_upper, locationname_upper


         quantity_requested = trim(quantityname)
         call str_upper(quantity_requested)
         location_requested = trim(locationname)
         call str_upper(location_requested)
         itemID = -1
         if (associated(instancePtr)) then
           frs:do i=1, instancePtr%nFileReaders
               fileReaderPtr => instancePtr%ecFileReadersPtr(i)%ptr
               if (fileReaderPtr%nItems<=0) cycle                                               ! No items to check
               location_supplied = fileReaderPtr%items(1)%ptr%elementSetPtr%name
               if (location_supplied/=location_requested) cycle     ! Items have the wrong location
               do j=1, fileReaderPtr%nItems
                  itemPtr => fileReaderPtr%items(j)%ptr
                  quantity_supplied = itemPtr%quantityPtr%name
                  call str_upper(quantity_supplied)
                  if (quantity_supplied==quantity_requested) then
                     found = .true.
                     if (present(isLateral)) then
                        found = .false.
                        if (isLateral) then
                           if (associated(fileReaderPtr%bc)) then
                              found = fileReaderPtr%bc%isLateral
                           endif
                        endif
                     endif
                     if (found) then
                        itemID = itemPtr%id
                        exit frs
                     endif
                  end if
               end do
            end do frs
         end if
      end function ecSupportFindItemByQuantityLocation
      ! =======================================================================

      !function ecSupportCreateTimeInterpolatedItem(instancePtr, sourceItemId, tgtNdx) result(itemId)
      !    use m_ec_item
      !    use m_ec_converter,  only: ecConverterSetType, ecConverterSetInterpolation, ecConverterSetOperand, ecConverterSetElement
      !    use m_ec_instance,   only: ecInstanceCreateConverter, ecInstanceCreateConnection, ecInstanceCreateItem, ecInstanceCreateField, ecInstanceCreateQuantity
      !    use m_ec_connection, only: ecConnectionAddTargetItem, ecConnectionAddSourceItem, ecConnectionSetConverter
      !    use m_ec_quantity,   only: ecQuantitySetName
      !    use m_ec_field,      only: ecFieldCreate1dArray
      !    use m_ec_item
      !
      !
      !    type(tEcInstance), pointer    :: instancePtr    !< EC-instance
      !    integer, intent(in)           :: sourceItemId   !< Source item id, before temporal interpolation
      !    integer, intent(in), optional :: tgtNdx         !< Optional target index, 1 is assumed as default
      !    integer                       :: targetItemId   !< Target item id, after temporal interpolation
      !    integer                       :: itemId         !< returned  target item ID, if successful, otherwise -1
      !    integer                       :: convertId
      !    type(tECItem), pointer        :: sourceItemPtr => null()
      !    type(tECItem), pointer        :: targetItemPtr => null()
      !    character(len=:), allocatable :: quantityName
      !    integer                       :: arraySize
      !
      !    integer :: targetIndex
      !    integer :: converterId, connectionId, quantityId, elementSetId, fieldId
      !
      !    if (present(tgtNdx)) then
      !       targetIndex = tgtNdx
      !    else
      !       targetIndex = 1
      !    end if
      !
      !    sourceItemPtr => ecSupportFindItem(instancePtr, sourceItemId)
      !
      !    ! TODO: create target item:
      !    !       . elementset-name = source_item's elementset-name
      !    !       . quantity-name = source_item's quantity-name + '-interpolated'
      !    itemId = -1
      !
      !    ! Set up the target item
      !    targetItemId = ecInstanceCreateItem(instancePtr)
      !    fieldId = ecInstanceCreateField(instancePtr)
      !
      !    arraySize = size(sourceItemPtr%sourceT0FieldPtr%arr1d)
      !    if (.not. (ecFieldCreate1dArray(instancePtr, fieldId, arraySize))) return
      !
      !    if (.not. ecItemSetRole(instancePtr, targetItemId, itemType_target)) return
      !    if (.not. ecItemSetTargetField(instancePtr, targetItemId, fieldId)) return
      !    if (.not. ecItemSetType(instancePtr, targetItemId, accessType_evaluate)) return
      !    quantityId = ecInstanceCreateQuantity(instancePtr)
      !    quantityName = trim(sourceItemPtr%quantityPtr%name)
      !    if (.not. ecItemSetQuantity(instancePtr, targetItemId, quantityId)) return
      !    if (.not. (ecQuantitySetName(instancePtr, quantityId, quantityName//'_interpolated'))) return
      !    elementSetId = sourceItemPtr%elementSetPtr%id
      !    if (.not. ecItemSetElementSet(instancePtr, targetItemId, elementSetId)) return
      !
      !    ! Construct a new Converter.
      !    converterId = ecInstanceCreateConverter(instancePtr)
      !
      !    ! Initialize the new Converter.
      !    if (.not. (ecConverterSetType(instancePtr, converterId, convType_uniform))) return
      !    if (.not. (ecConverterSetOperand(instancePtr, converterId, operand_replace_element))) return
      !    if (.not. (ecConverterSetInterpolation(instancePtr, converterId, interpolate_timespace))) return
      !    if (.not. (ecConverterSetElement(instancePtr, converterId, targetIndex))) return
      !
      !    ! Construct a new Connection.
      !    connectionId = ecInstanceCreateConnection(instancePtr)
      !    if (.not. ecConnectionSetConverter(instancePtr, connectionId, converterId)) return
      !
      !    ! Initialize the new Connection.
      !    if (.not. ecConnectionAddSourceItem(instancePtr, connectionId, sourceItemId)) return
      !    if (.not. ecConnectionAddTargetItem(instancePtr, connectionId, targetItemId)) return
      !    if (.not. ecItemAddConnection(instancePtr, targetItemId, connectionId)) return
      !    itemId = targetItemId
      !end function ecSupportCreateTimeInterpolatedItem
      ! =======================================================================


      !> Retrieve the pointer to the Item with id == itemId.
      function ecSupportFindItem(instancePtr, itemId) result(itemPtr)
         type(tEcItem),     pointer            :: itemPtr     !< Item corresponding to itemId
         type(tEcInstance), pointer            :: instancePtr !< intent(in)
         integer,                   intent(in) :: itemId      !< unique Item id
         !
         integer :: i !< loop counter
         !
         itemPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nItems
               if (instancePtr%ecItemsPtr(i)%ptr%id == itemId) then
                  itemPtr => instancePtr%ecItemsPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindItem: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindItem

      ! =======================================================================

      !> Retrieve the pointer to the Connection with id == connectionId.
      function ecSupportFindConnection(instancePtr, connectionId) result(connectionPtr)
         type(tEcConnection),     pointer            :: connectionPtr !< Item corresponding to connectionId
         type(tEcInstance),       pointer            :: instancePtr   !< intent(in)
         integer,                         intent(in) :: connectionId  !< unique Connection id
         !
         integer :: i !< loop counter
         !
         connectionPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nConnections
               if (instancePtr%ecConnectionsPtr(i)%ptr%id == connectionId) then
                  connectionPtr => instancePtr%ecConnectionsPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindConnection: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindConnection

      ! =======================================================================

      !> Retrieve the pointer to the Converter with id == converterId.
      function ecSupportFindConverter(instancePtr, converterId) result(converterPtr)
         type(tEcConverter), pointer            :: converterPtr !< Item corresponding to converterId
         type(tEcInstance),  pointer            :: instancePtr  !< intent(in)
         integer,                    intent(in) :: converterId  !< unique Converter id
         !
         integer :: i !< loop counter
         !
         converterPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nConverters
               if (instancePtr%ecConvertersPtr(i)%ptr%id == converterId) then
                  converterPtr => instancePtr%ecConvertersPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindConverter: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindConverter

      ! =======================================================================

      !> Retrieve the pointer to the FileReader with id == converterId.
      function ecSupportFindFileReader(instancePtr, fileReaderId) result(fileReaderPtr)
         type(tEcFileReader), pointer            :: fileReaderPtr !< FileReader corresponding to fileReaderId
         type(tEcInstance),   pointer            :: instancePtr   !< intent(in)
         integer,                     intent(in) :: fileReaderId  !< unique FileReader id
         !
         integer :: i !< loop counter
         !
         fileReaderPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nFileReaders
               if (instancePtr%ecFileReadersPtr(i)%ptr%id == fileReaderId) then
                  fileReaderPtr => instancePtr%ecFileReadersPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindFileReader: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindFileReader

      ! =======================================================================
      !> Retrieve the pointer to the FileReader with id == converterId.
      function ecSupportFindFileReaderByFilename(instancePtr, filename) result(fileReaderPtr)
         type(tEcFileReader), pointer            :: fileReaderPtr !< FileReader corresponding to fileReaderId
         type(tEcInstance),   pointer            :: instancePtr   !< intent(in)
         character(*),        intent(in)         :: filename      !< relative path of data file
         !
         integer :: i !< loop counter
         !
         fileReaderPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nFileReaders
               if (associated(instancePtr%ecFileReadersPtr(i)%ptr%bc)) then                  ! if filereader has bc-block
                  if (strcmpi(instancePtr%ecFileReadersPtr(i)%ptr%bc%fName, fileName)) then  ! this bc-block has the filename
                     fileReaderPtr => instancePtr%ecFileReadersPtr(i)%ptr
                  end if
               else                                                                          ! else
                  if (strcmpi(instancePtr%ecFileReadersPtr(i)%ptr%fileName, fileName)) then  ! the filereader has the filename
                     fileReaderPtr => instancePtr%ecFileReadersPtr(i)%ptr
                  end if
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindFileReader: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindFileReaderByFilename

      ! =======================================================================

      !> Retrieve the pointer to the BCBlock with id == bcBlockId.
      function ecSupportFindBCBlock(instancePtr, bcBlockId) result(bcBlockPtr)
         type(tEcBCBlock), pointer               :: bcBlockPtr    !< BCBlock corresponding to bcBlockId
         type(tEcInstance), pointer              :: instancePtr   !< intent(in)
         integer, intent(in)                     :: bcBlockId     !< unique BCBlock id
         !
         integer :: i !< loop counter
         !
         bcBlockPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nBCBlocks
               if (instancePtr%ecBCBlocksPtr(i)%ptr%id == bcBlockId) then
                  bcBlockPtr => instancePtr%ecBCBlocksPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindBCBlock: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindBCBlock

      ! =======================================================================

      !> Retrieve the pointer to the NetCDF instance with id == netCDFId.
      function ecSupportFindNetCDF(instancePtr, netCDFId) result(netCDFPtr)
         type(tEcNetCDF), pointer                :: netCDFPtr     !< NetCDF instance for the given Id
         type(tEcInstance), pointer              :: instancePtr   !< intent(in)
         integer, intent(in)                     :: netCDFId      !< unique NetCDF id
         !
         integer :: i !< loop counter
         !
         netCDFPtr => null()
         !
         if (associated(instancePtr)) then
            do i=1, instancePtr%nNetCDFs
               if (instancePtr%ecNetCDFsPtr(i)%ptr%id == netCDFId) then
                  netCDFPtr => instancePtr%ecNetCDFsPtr(i)%ptr
               end if
            end do
         else
            call setECMessage("ec_support::ecSupportFindNetCDF: Dummy argument instancePtr is not associated.")
         end if
      end function ecSupportFindNetCDF

      ! =======================================================================
      !> Find the file reader for the bc block that contains the component definition for the comp.correction block
      function ecSupportFindRelatedBCBlock(instancePtr, corrFileReaderPtr, qname, bcname, func) result(cmpFileReaderPtr)
         type(tEcFileReader), pointer :: cmpFileReaderPtr  !< resulting file reader
         type(tEcInstance),   pointer :: instancePtr       !< intent(in)
         type(tEcFileReader), pointer :: corrFileReaderPtr !< intent(inout)
         character(len=*)             :: qname             !< quantity name
         character(len=*)             :: bcname            !< point on poly name
         integer, intent(in)          :: func              !< function type

         integer                      :: iFileReader
         type (tEcBCBlock), pointer   :: BCBlockptr
         !
         cmpFileReaderPtr => null()
         do iFileReader = 1, instancePtr%nFileReaders
            BCBlockptr => instancePtr%EcFileReadersPtr(iFileReader)%ptr%bc
            if (associated(BCBlockptr)) then
               if (trim(BCBlockptr%bcname)==trim(bcname) .and. (trim(BCBlockptr%qname)==trim(qname))  &
                                                         .and. BCBlockptr%func == func) then
                  cmpFileReaderPtr => instancePtr%EcFileReadersPtr(iFileReader)%ptr
                  exit
               endif
            else
            endif
         enddo
      end function ecSupportFindRelatedBCBlock


      !> Translate NetCDF error code into a NetCDF error message.
      function ecSupportNetcdfCheckError(ierror, description, filename) result(success)
         use netcdf
         !
         logical                                             :: success     !<
         integer,                       intent(in)           :: ierror      !<
         character(len=*),              intent(in), optional :: description !<
         character(len=maxFileNameLen), intent(in), optional :: filename    !<
         !
         character(3000) :: message
         !
         if (ierror /= nf90_noerr) then
            if (present(description) .and. present(filename)) then
               write (message,'(6a)') trim(description), '. NetCDF file : "', trim(filename), '". Error message:', nf90_strerror(ierror)
            else
               write (message,'(2a)') 'NetCDF error : ', nf90_strerror(ierror)
            endif
            call setECMessage(message)
            success = .false.
         else
            success = .true.
         endif
      end function ecSupportNetcdfCheckError

      ! =======================================================================

      !> Extracts time unit and reference date from a standard time string.
      !! ASCII example: "TIME = 0 hours since 2006-01-01 00:00:00 +00:00"
      !! NetCDF example: "minutes since 1970-01-01 00:00:00.0 +0000"
      function ecSupportTimestringToUnitAndRefdate(rec, unit, ref_date, tzone) result(success)
         !
         logical                                :: success       !< function status
         character(len=*),        intent(in)    :: rec           !< units string (at out in lowercase)
         integer,                 intent(out)   :: unit          !< time unit enumeration
         real(kind=hp),           intent(out)   :: ref_date      !< reference date formatted as Modified Julian Date
         real(kind=hp), optional, intent(out)   :: tzone         !< time zone
         !
         integer                       :: i        !< helper index for location of 'since'
         integer                       :: j        !< helper index for location of '+/-' in time zone
         integer                       :: jplus    !< helper index for location of '+' in time zone
         integer                       :: jmin     !< helper index for location of '-' in time zone
         integer                       :: jcomment !< helper index for location of '#'
         integer                       :: minsize  !< helper index for time zone
         real(kind=hp)                 :: temp     !< helper variable
         logical                       :: ok       !< check of refdate is found
         character(len=20)             :: date     !< parts of string for date
         character(len=20)             :: time     !< parts of string for time
         character(len=:), allocatable :: string   !< unit string without comments and in lowercase
         !
         success = .false.
         !
         ! copy only relevant part of rec into string:
         jcomment = index(rec, '#')
         if (jcomment == 1) then
            call setECMessage("ec_support::ecSupportTimestringToUnitAndRefdate: only found comments.")
            return
         else if (jcomment > 1) then
            string = trim(rec(:jcomment - 1))
         else
            string = trim(rec)
         endif

         call str_lower(string)

         ! Determine the time unit.
         if (index(string, 'seconds') /= 0) then
            unit = ec_second
         else if (index(string, 'minutes') /= 0) then
            unit = ec_minute
         else if (index(string, 'hours') /= 0 .or. index( string, 'hrs') /= 0) then
            unit = ec_hour
         else if (index(string, 'days') /= 0) then
            unit = ec_day
         else
            call setECMessage("unitstring = '"//trim(string)//"'.")
            call setECMessage("ec_support::ecSupportTimestringToUnitAndRefdate: Unable to identify the time unit.")
            return
         end if
         ! Determine the reference date.
         i = index(string, 'since') + 6
         call split_date_time(string(i:), date, time)
         if (i /= 6) then
            ! Date
            if (ymd2reduced_jul(date, ref_date)) then
               ! Time
               if(len_trim(time)>=8) then
                  read(time(1 : 2), *) temp
                  ref_date = ref_date + temp / 24.0_hp
                  read(time(4 : 5), *) temp
                  ref_date = ref_date + temp / 24.0_hp / 60.0_hp
                  read(time(7 : 8), *) temp
                  ref_date = ref_date + temp / 24.0_hp / 60.0_hp / 60.0_hp
               end if
               ok = .true.
            else
               ref_date = -999.0_hp
               ok = .false.
            endif
         else
            ok = ecSupportTimestringArcInfo(string, ref_date)
         endif
         if (.not. ok) then
            call setECMessage("ec_support::ecSupportTimestringToUnitAndRefdate: Unable to identify keyword: since.")
            return
         end if

         ! Determine the timezone
         if (present(tzone)) then
             minsize = i + 18   ! +/- is to be found after date and time (note that date has '-')
             jplus = index(string, '+', back=.true.)
             jmin  = index(string, '-', back=.true.)
             j     = max(jplus, jmin)
             if (j > minsize) then
                 success = parseTimezone(string(j:), tzone)
             else
                 tzone = 0.0_hp
                 success = .true.
             endif
         else
             success = .true.
         end if
         !
      end function ecSupportTimestringToUnitAndRefdate

      !> Extracts time unit and reference date from a time string in Arc Info format.
      !! example: ... TIME (HRS)     18.0 20000101 18
      function ecSupportTimestringArcInfo(rec, ref_date, time_steps) result (success)
         character(len=*)       , intent(in)  :: rec        !< input string
         real(kind=hp), optional, intent(out) :: ref_date   !< reference date found
         real(kind=hp), optional, intent(out) :: time_steps !< time step found
         logical                              :: success    !< function result

         integer       :: yyyymmdd    !< reference date as Gregorian yyyymmdd
         integer       :: posHrs      !< position in a string of '(HRS)', 'hrs', 'hours'
         integer       :: posNumbers  !< first position of the numbers in a string (actually, the first space after '(HRS)')
         integer       :: posTime     !< position in a string of 'TIME' or 'time'
         integer       :: ierr        !< error code
         integer       :: i           !< loop counter
         real(kind=hp) :: time        !< time found
         integer       :: hh          !< hour in refdate found

         success = .false.
         posNumbers = 0

         posTime = max(index(rec, 'TIME'), index(rec, 'time'))

         if (posTime > 0) then
            posHrs = max(index(rec, '(HRS)'), index(rec, 'hrs'), index(rec, 'hours'))
            do i = posHrs+3, len_trim(rec)
               if (rec(i:i) == ' ') then
                  posNumbers = i
                  exit
               endif
            enddo
         endif

         if (present(ref_date)) then
            ref_date = -999.0_hp   ! initialize

            if (posNumbers > 0) then
               read(rec(posNumbers:), *, iostat=ierr) time, yyyymmdd, hh
               if (ierr /= 0) then
                  ! may be hh is missing, try again:
                  read(rec(posNumbers:), *, iostat=ierr) time, yyyymmdd
                  hh = 0
               endif
            endif

            if (ierr == 0) then
               success = .true.
               ref_date = JULIAN(yyyymmdd, 0)
               ref_date = ref_date + real(hh, hp) / 24.0_hp
            endif
         endif

         if (present(time_steps)) then
            time_steps = -999.0_hp   ! initialize

            if (posNumbers > 0) then
               read(rec(posNumbers:), *, iostat=ierr) time
            endif

            if (ierr == 0) then
               time_steps = time
               success = .true.
            endif
         endif

      end function ecSupportTimestringArcInfo

      ! =======================================================================

      !> Extracts time zone from a standard time string.
      !! examples: "+01:00", "+0200", "-01:00", "-0200", "+5:30"
      function parseTimezone(string, tzone) result(success)
         logical                              :: success         !< function status
         character(len=*),        intent(in)  :: string          !< units string
         real(kind=hp),           intent(out) :: tzone           !< time zone

         integer          :: ierr        !< error code
         integer          :: jcolon      !< helper index for location of ':' in timezone
         integer          :: jend        !< helper index
         integer          :: posNulChar  !< position of null char; end of string if from C
         real(kind=hp)    :: min         !< minutes part of time zone, as double
         real(kind=hp)    :: hour        !< hours part of time zone, as double
         character(len=2) :: cmin        !< minutes part of time zone, as character string
         character(len=3) :: chour       !< hours part of time zone, as character string

         tzone = 0.0_hp

         jcolon = index(string, ':')

         if (jcolon == 0) then
            posNulChar = index(string, char(0))
            if (posNulChar > 0) then
               jend = posNulChar-1
            else
               jend = len_trim(string)
            endif
            cmin = string(jend-1:)
            chour = string(:jend-2)
         else
            cmin = string(jcolon+1:)
            chour = string(:jcolon-1)
         end if

         read(chour, *, iostat=ierr) hour
         if (ierr == 0) read(cmin, *, iostat=ierr) min
         if (string(1:1) == '-') then
            tzone = hour - min / 60.0_hp
         else
            tzone = hour + min / 60.0_hp
         endif

         success = (ierr == 0)
         if (.not. success) call setECMessage("ec_support::parseTimezone: error parsing time zone " // trim(string))
      end function parseTimezone

      ! =======================================================================

      !> Calculate conversion factor from ec_timestep_unit to seconds
      function ecSupportTimeUnitConversionFactor(unit) result(factor)
         real(kind=hp)       :: factor
         integer, intent(in) :: unit !< time unit enum
         !
         factor = 1.0_hp ! default return value
         !
         if (unit == ec_second) then
            factor = 1.0_hp
         else if (unit == ec_minute) then
            factor = 60.0_hp
         else if (unit == ec_hour) then
            factor = 3600.0_hp
         else if (unit == ec_day) then
            factor = 3600.0_hp * 24.0_hp
         end if
      end function ecSupportTimeUnitConversionFactor

      ! =======================================================================
      !> Convert seconds since k_refdate to times(i) * ec_timestep_unit since ec_refdate to seconds since k_refdate.
      function ecSupportMJDToTimeIndex(tframe, time_mjd) result(ndx)
         integer                        :: ndx       !< function result, largest index with a time less than timesteps
         type(tEcTimeFrame), intent(in) :: tframe    !< TimeFrame containing input data for conversion
         real(hp)          , intent(in) :: time_mjd  !< seconds since k_refdate representing time to be found
         !
         real(hp):: srctime
         integer :: i
         !
         srctime = ecSupportMJDToThisTime(tframe, time_mjd)
         do i=1,tframe%nr_timesteps
            if (srctime<=tframe%times(i)) exit
         enddo
         if (i>tframe%nr_timesteps) then
            i = -1
         else
            i = max(i - 1,1)
         endif
         ndx = i
      end function ecSupportMJDToTimeIndex


      ! =======================================================================
      !> Convert seconds since k_refdate to times(i) * ec_timestep_unit since ec_refdate to seconds since k_refdate.
      function ecSupportTimeIndexToMJD(tframe, ndx) result(time_mjd)
         real(hp)                       :: time_mjd  !< seconds since k_refdate representing time to be found
         type(tEcTimeFrame), intent(in) :: tframe    !< TimeFrame containing input data for conversion
         integer                        :: ndx       !< function result, largest index with a time less than timesteps
         !
         time_mjd = ecSupportThisTimeToMJD(tframe, tframe%times(ndx))
      end function ecSupportTimeIndexToMJD


      !> Convert from MJD to timefrane
      !> Convert thistime * ec_timestep_unit since ec_refdate to seconds since k_refdate.
      function ecSupportMJDToThisTime(tframe, time_mjd) result(thistime)
         real(hp)                       :: thistime  !< function result, time wrt timeframe
         type(tEcTimeFrame), intent(in) :: tframe    !< TimeFrame containing input data for conversion
         real(hp),           intent(in) :: time_mjd  !< this time in mjd needs to be converted
         !
         real(hp)                       :: factor       !< resulting conversion factor
         !
         factor = ecSupportTimeUnitConversionFactor(tframe%ec_timestep_unit) / 86400.0_hp    ! Converts to DAYS
         thistime = (time_mjd + tframe%ec_timezone/24.0_hp - tframe%ec_refdate)/factor
      end function ecSupportMJDToThisTime


      ! =======================================================================

      !> Convert from timeframe to MJD
      function ecSupportThisTimeToMJD(tframe, thistime) result(time_mjd)
         real(hp)                       :: time_mjd  !< function result, modified julian days
         type(tEcTimeFrame), intent(in) :: tframe    !< TimeFrame containing input data for conversion
         real(hp),           intent(in) :: thistime  !< time defined wrt timeframe to be converted
         !
         real(hp)                       :: factor       !< resulting conversion factor
         !
         factor = ecSupportTimeUnitConversionFactor(tframe%ec_timestep_unit) / 86400.0_hp    ! Converts to DAYS
         time_mjd = tframe%ec_refdate + thistime * factor - tframe%ec_timezone/24.0_hp
      end function ecSupportThisTimeToMJD

      ! =======================================================================

      !> Find the CF-compliant longitude and latitude dimensions and associated variables
      function ecSupportNCFindCFCoordinates(ncid, lon_varid, lon_dimid, lat_varid, lat_dimid,      &
                                             grid_lon_varid, grid_lat_varid,                       &
                                                    x_varid,   x_dimid,   y_varid,   y_dimid,      &
                                                    z_varid,   z_dimid,                            &
                                                  tim_varid, tim_dimid) result(success)
      use netcdf
      logical              :: success
      integer, intent(in)  :: ncid           !< NetCDF file ID
      integer, intent(out) :: lon_varid      !< One dimensional coordinate variable recognized as absolute longitude
      integer, intent(out) :: lat_varid      !< One dimensional coordinate variable recognized as absolute latitude
      integer, intent(out) :: grid_lon_varid !< One dimensional coordinate variable recognized as 'rotated-pole' longitude
      integer, intent(out) :: grid_lat_varid !< One dimensional coordinate variable recognized as 'rotated-pole' latitude
      integer, intent(out) ::   x_varid      !< One dimensional coordinate variable recognized as X
      integer, intent(out) ::   y_varid      !< One dimensional coordinate variable recognized as Y
      integer, intent(out) ::   z_varid      !< One dimensional coordinate variable recognized as Z
      integer, intent(out) :: tim_varid      !< One dimensional coordinate variable recognized as time
      integer, intent(out) :: lon_dimid      !< Longitude dimension
      integer, intent(out) :: lat_dimid      !< Latitude dimension
      integer, intent(out) ::   x_dimid      !< X dimension
      integer, intent(out) ::   y_dimid      !< Y dimension
      integer, intent(out) ::   z_dimid      !< Z dimension
      integer, intent(out) :: tim_dimid      !< Time dimension
      integer :: ndim, nvar, ivar, nglobatts, unlimdimid, ierr
      integer :: dimids(1)
      character(len=NF90_MAX_NAME)  :: units, axis, varname, stdname

      success = .False.
      lon_varid = -1
      lat_varid = -1
      x_varid = -1
      y_varid = -1
      z_varid = -1
      tim_varid = -1
      lon_dimid = -1
      lat_dimid = -1
      x_dimid = -1
      y_dimid = -1
      z_dimid = -1
      tim_dimid = -1
      ierr = nf90_inquire(ncid, ndim, nvar, nglobatts, unlimdimid)
      do ivar=1,nvar
         ierr = nf90_inquire_variable(ncid, ivar, ndims=ndim)                      ! number of variables
         units=''
         ierr = nf90_get_att(ncid, ivar, 'units', units)
         if (ndim==1) then
            ierr = nf90_inquire_variable(ncid, ivar, dimids=dimids)                ! number of variables
            select case (trim(units))
               case ('degrees_east','degree_east','degree_E','degrees_E','degreeE','degreesE')
                  lon_varid = ivar
                  lon_dimid = dimids(1)
               case ('degrees_north','degree_north','degree_N','degrees_N','degreeN','degreesN')
                  lat_varid = ivar
                  lat_dimid = dimids(1)
               case ('degrees')
                   stdname = ''
                   ierr = nf90_get_att(ncid, ivar, 'standard_name', stdname)
                   if (ierr == 0) then
                      select case (stdname)
                         case ('grid_longitude')
                            grid_lon_varid = ivar
                            lon_dimid = dimids(1)
                         case ('grid_latitude')
                            grid_lat_varid = ivar
                            lat_dimid = dimids(1)
                      end select
                   else
                      call setECmessage("attribute 'standard_name' not found for variable " // trim(varname))
                   endif
                   !RL Set lon and lat dimids ??
               case ('m','meters','km','kilometers')
                  axis=''
                  ierr = nf90_get_att(ncid, ivar, 'axis', axis)
                  if (ierr /= 0) ierr = nf90_get_att(ncid, ivar, 'AXIS', axis) ! support 'axis' in upper case, lower case and camel case
                  if (ierr /= 0) ierr = nf90_get_att(ncid, ivar, 'Axis', axis)

                  if (ierr /= 0) then
                      ierr = nf90_inquire_variable(ncid, ivar, name = varname)
                      call setECmessage("attribute 'axis' not found for variable " // trim(varname))
                  else if (strcmpi(axis,'X')) then
                     x_varid = ivar
                     x_dimid = dimids(1)
                  else if (strcmpi(axis,'Y')) then
                     y_varid = ivar
                     y_dimid = dimids(1)
                  else if (strcmpi(axis,'Z')) then
                     z_varid = ivar
                     z_dimid = dimids(1)
                  end if
               case default
                  ! see if is the time dimension
                  if ((index(units,'seconds since')>0)   &
                  .or.(index(units,'minutes since')>0)   &
                  .or.(index(units,'hours since')>0)     &
                  .or.(index(units,'days since')>0))     then
                     tim_varid = ivar
                     tim_dimid = dimids(1)
                  end if
               end select
         end if
         if (ndim==2) then                ! Find lat and lon even if they are no coordinate axis
            select case (trim(units))
               case ('degrees_east','degree_east','degree_E','degrees_E','degreeE','degreesE')
                  lon_varid = ivar
               case ('degrees_north','degree_north','degree_N','degrees_N','degreeN','degreesN')
                  lat_varid = ivar
               case ('degrees')
                   stdname = ''
                   ierr = nf90_get_att(ncid, ivar, 'standard_name', stdname)
                   if (ierr == 0) then
                      select case (stdname)
                         case ('grid_latitude')
                            grid_lat_varid = ivar
                         case ('grid_longitude')
                            grid_lon_varid = ivar
                      end select
                   else
                      call setECmessage("attribute 'standard_name' not found for variable " // trim(varname))
                   endif
               case ('m','meters','km','kilometers')
                   stdname = ''
                   ierr = nf90_get_att(ncid, ivar, 'standard_name', stdname)
                   select case (stdname)
                      case ('projection_x_coordinate')
                         x_varid = ivar
                      case ('projection_y_coordinate')
                         y_varid = ivar
                   end select
               end select
         end if
      end do   !ivar

      success = .True.
      end function ecSupportNCFindCFCoordinates

      !!> Return dimension id's and lengths for the specified variable
      !function ecSupportNCGetVarDim(ncid, varid, dimid, dimlen, ndim) result(success)
      !use netcdf
      !logical              :: success
      !integer, intent(in)  :: ncid                      !< NetCDF file ID
      !integer, intent(in)  :: varid                     !< Variable file ID
      !
      !integer, dimension(:), allocatable :: intent(out) !< array of dimension ID's for this variable
      !integer, dimension(:), allocatable :: intent(out) !< array of dimension s for this variable
      !end function ecSupportNCGetVarDim

end module m_ec_support

! =============================================================================

