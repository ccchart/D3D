!----- LGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2018.                                
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

!> This module contains all the methods for the datatype tEcQuantity.
!! @author arjen.markus@deltares.nl
!! @author adri.mourits@deltares.nl
!! @author stef.hummel@deltares.nl
!! @author edwin.bos@deltares.nl
module m_ec_quantity
   use m_ec_typedefs
   use m_ec_message
   use m_ec_support
   use m_ec_alloc
   
   implicit none
   
   private
   
   public :: ecQuantityCreate
   public :: ecQuantityFree1dArray
   public :: ecQuantitySet
   public :: ecQuantitySetUnitsFillScaleOffsetFromNcidVarid
   
   contains
      
      ! =======================================================================
      
      !> Construct a new Quantity with the specified id.
      !! Failure is indicated by returning a null pointer.
      function ecQuantityCreate(quantityId) result(quantityPtr)
         type(tEcQuantity), pointer    :: quantityPtr !< the new Quantity, intent(out)
         integer,           intent(in) :: quantityId  !< unique Quantity id
         !
         integer :: istat !< allocate() status
         !
         ! allocation
         allocate(quantityPtr, stat = istat)
         if (istat /= 0) then
            call setECMessage("ERROR: ec_quantity::ecQuantityCreate: Unable to allocate additional memory.")
            quantityPtr => null()
            return
         end if
         ! initialization
         quantityPtr%id = quantityId
         quantityPtr%name = ' '
         quantityPtr%units = ' '
      end function ecQuantityCreate
      
      ! =======================================================================
      
      !> Frees a 1D array of tEcQuantityPtrs, after which the quantityPtr is deallocated.
      function ecQuantityFree1dArray(quantityPtr, nQuantities) result (success)
         logical                                     :: success     !< function status
         type(tEcQuantityPtr), dimension(:), pointer :: quantityPtr !< intent(inout)
         integer, intent(inout)                      :: nQuantities !< number of Quantities
         !
         integer :: i     !< loop counter
         integer :: istat !< deallocate() status
         !
         success = .true.
         !
         if (.not. associated(quantityPtr)) then
            call setECMessage("WARNING: ec_quantity::ecQuantityFree1dArray: Dummy argument quantityPtr is already disassociated.")
         else
            ! Free(nothing to do) and deallocate all tEcQuantityPtrs in the 1d array.
            do i=1, nQuantities
               deallocate(quantityPtr(i)%ptr, stat = istat)
               if (istat /= 0) success = .false.
            end do
            ! Finally deallocate the tEcQuantityPtr(:) pointer.
            if (success) then
               deallocate(quantityPtr, stat = istat)
               if (istat /= 0) success = .false.
            end if
         end if
         nQuantities = 0
      end function ecQuantityFree1dArray
      
      ! =======================================================================
      
       function ecQuantitySet(instancePtr, quantityId, name,                     &
                                                       units,                    &
                                                       vectormax,                &
                                                       factor,                   &
                                                       offset,                   &
                                                       fillvalue,                &
                                                       timeint,                  &
                                                       ncid                      &
                                                     ) result(success)

         logical                               :: success     !< function status
         type(tEcInstance), pointer            :: instancePtr !< intent(in)
         integer,                   intent(in) :: quantityId  !< unique Quantity id

         character(len=*),optional :: name
         character(len=*),optional :: units
         integer         ,optional :: vectormax
         real(hp)        ,optional :: factor
         real(hp)        ,optional :: offset
         real(hp)        ,optional :: fillvalue
         integer         ,optional :: timeint
         integer         ,optional :: ncid
         !
         type(tEcQuantity), pointer :: quantityPtr !< Quantity corresponding to quantityId
         !
         success = .false.
         quantityPtr => null()
         !
         quantityPtr => ecSupportFindQuantity(instancePtr, quantityId)
         if (.not.associated(quantityPtr)) then
            call setECMessage("ERROR: ec_quantity::ecQuantitySet: Cannot find a Quantity with the supplied id.")
         end if
         if (present(name)) then 
            if (len_trim(name) > maxNameLen) then
               call setECMessage("ERROR: ec_quantity::ecQuantitySetName: The new name string is too long, unable to change name.")
            else
               quantityPtr%name = name  
            end if
         end if

         if (present(units)) quantityPtr%units = units  
         if (present(vectormax)) quantityPtr%vectormax = vectormax  
         if (present(factor)) quantityPtr%factor = factor  
         if (present(offset)) quantityPtr%offset = offset  
         if (present(fillvalue)) quantityPtr%fillvalue = fillvalue  
         if (present(timeint)) quantityPtr%timeint = timeint  
         if (present(ncid)) quantityPtr%ncid = ncid  
         success = .true.
      end function ecQuantitySet
      
      !> Change the Units, Fillvalue, Scalefactor and Offset shift of the Quantity corresponding to quantityId
      !> obtained from the variable varid in the netcdf file ncid 
      !> all in try-catch fashion: if not available, leave empty or use default
      function ecQuantitySetUnitsFillScaleOffsetFromNcidVarid(instancePtr, quantityId, ncid, varid) result(success)
      use netcdf
      use string_module
         implicit none
         logical                               :: success     !< function status
         type(tEcInstance), pointer            :: instancePtr !< intent(in)
         integer,                   intent(in) :: quantityId  !< unique Quantity id
         integer,                   intent(in) :: ncid        !< id of the nc-file 
         integer,                   intent(in) :: varid       !< id of the variable
                                                              !< order: new = (old*scale) + offset
         character(len=:), allocatable  :: units
         integer  :: ierr
         integer  :: attriblen
         real(hp) :: add_offset, scalefactor, fillvalue

         success = .false.
         add_offset = 0.d0
         scalefactor = 1.d0
         fillvalue = ec_undef_hp
         attriblen=0
         ierr = nf90_inquire_attribute(ncid, varid, 'units', len=attriblen)
         if (attriblen>0) then
            allocate(character(len=attriblen) :: units) 
            units(1:len(units)) = ''
            if (nf90_get_att(ncid, varid, 'units', units)==NF90_NOERR) then 
               call str_upper(units) ! make units attribute case-insensitive 
               if (.not.(ecQuantitySet(instancePtr, quantityId, units=units))) return
            end if
         end if
         if (nf90_get_att(ncid, varid, '_FillValue', fillvalue)==NF90_NOERR) then                  ! RL: Possibly redundant: we store the missing value with the field
            if (.not.(ecQuantitySet(instancePtr, quantityId, fillvalue=fillvalue))) return         !     And not with the quantity. TODO: check if this can be removed
         end if
         if ((nf90_get_att(ncid, varid, 'scale_factor', scalefactor)==NF90_NOERR)         &
              .or. (nf90_get_att(ncid, varid, 'add_offset', add_offset)==NF90_NOERR)) then
              if (.not.(ecQuantitySet(instancePtr, quantityId, factor=scalefactor, offset=add_offset))) return
         end if
         success = .true.
      end function ecQuantitySetUnitsFillScaleOffsetFromNcidVarid
     
end module m_ec_quantity
