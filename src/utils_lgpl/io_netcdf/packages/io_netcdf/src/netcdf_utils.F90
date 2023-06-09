!----- LGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2011-2023.
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
!
!-------------------------------------------------------------------------------

! 
! 

!> Utility module for additional manipulation/inquiry of NetCDF files, on top of the basic nf90* primitives.
module netcdf_utils
use netcdf
use ionc_constants
implicit none

private

public :: ncu_format_to_cmode
public :: ncu_ensure_define_mode
public :: ncu_ensure_data_mode
public :: ncu_restore_mode
public :: ncu_inq_var_fill, ncu_copy_atts, ncu_copy_chunking_deflate
public :: ncu_clone_vardef
public :: ncu_append_atts
public :: ncu_get_att
public :: ncu_get_var_attset
public :: ncu_put_var_attset
public :: ncu_att_to_varid
public :: ncu_att_to_dimid
public :: ncu_apply_to_att

integer, parameter :: maxMessageLen = 1024  ! copy taken from io_ugrid
character(len=maxMessageLen) :: ncu_messagestr !< Placeholder string for storing diagnostic messages. /see{ug_get_message}

! Copied from official NetCDF: typeSizes.f90
integer, parameter ::   OneByteInt = selected_int_kind(2), &
                        TwoByteInt = selected_int_kind(4), &
                       FourByteInt = selected_int_kind(9), &
                      EightByteInt = selected_int_kind(18)
integer, parameter ::                                          &
                      FourByteReal = selected_real_kind(P =  6, R =  37), &
                     EightByteReal = selected_real_kind(P = 13, R = 307)

!> Compatibility function: returns the fill settings for a variable in a netCDF-3 file.
interface ncu_inq_var_fill
   module procedure ncu_inq_var_fill_int4
   module procedure ncu_inq_var_fill_real8
end interface ncu_inq_var_fill

abstract interface
   function ncu_apply_to_att(attname, attvalue) result(ierr)
      character(len=*),              intent(in   ) :: attname  !< Name of the attribute, cannot be changed.
      character(len=:), allocatable, intent(inout) :: attvalue !< value of the attribute, can be changed. Should be an allocatable character string.
      integer                                      :: ierr     !< Result status (recommended IONC_NOERR if successful)
   end function ncu_apply_to_att
end interface

   contains


!> Returns the NetCDF creation mode flag value, given the colloquial
!! format number (3 or 4).
!! Use this returned value as cmode argument to nf90_create calls.
!!
!! NOTE: the input number is *not* equivalent with the library's
!! NF90_FORMAT_* constants!
pure function ncu_format_to_cmode(iformatnumber) result(cmode)
   integer, intent(in) :: iformatnumber !< The NetCDF format version (3 or 4, colloquially speaking)
   integer             :: cmode         !< Return value (for example NF90_CLASSIC_MODEL or NF90_NETCDF4), ready for use in nf90_create calls.

   select case(iformatnumber)
   case(3)
      cmode = NF90_CLASSIC_MODEL
   case(4)
      cmode = NF90_NETCDF4
   case default
      cmode = NF90_CLOBBER ! 0: use library default
   end select
end function ncu_format_to_cmode


!> Puts a NetCDF dataset into define mode, if it isn't so already.
!! Returns the original mode, such that the caller can later put back this
!! dataset into its original mode.
!!
!! @see ncu_restore_mode
function ncu_ensure_define_mode(ncid, originally_in_define) result(ierr)
   integer, intent(in   ) :: ncid                 !< ID of the NetCDF dataset
   logical, intent(  out) :: originally_in_define !< Whether the dataset was already in define mode (for later restoral).
   integer                :: ierr                 !< Result status (nf90_noerr if successful)

   integer :: ierrloc

   ierr = nf90_noerr

   ! Put dataset in define mode (possibly again)
   originally_in_define = .false.
   
   ierrloc = nf90_redef(ncid)
   if (ierrloc == nf90_eindefine) then
      originally_in_define = .true.
   else
      ierr = ierrloc ! Some other error occurred
   end if
end function ncu_ensure_define_mode


!> Puts a NetCDF dataset into data mode, if it isn't so already.
!! Returns the original mode, such that the caller can later put back this
!! dataset into its original mode.
!!
!! @see ncu_restore_mode
function ncu_ensure_data_mode(ncid, originally_in_define) result(ierr)
   integer, intent(in   ) :: ncid                 !< ID of the NetCDF dataset
   logical, intent(  out) :: originally_in_define !< Whether the dataset was originally in define mode (for later restoral).
   integer                :: ierr                 !< Result status (nf90_noerr if successful)

   integer :: ierrloc

   ierr = nf90_noerr

   ! Put dataset in data mode (possibly again)
   originally_in_define = .true.
   
   ierrloc = nf90_enddef(ncid)
   if (ierrloc == nf90_enotindefine) then
      originally_in_define = .false.
   else
      ierr = ierrloc ! Some other error occurred
   end if
end function ncu_ensure_data_mode


!> Restores a NetCDF dataset into its original mode (define/data), if needed.
!!
!! @see ncu_ensure_define_mode
!! @see ncu_ensure_data_mode
function ncu_restore_mode(ncid, originally_in_define) result(ierr)
   integer, intent(in   ) :: ncid                 !< ID of the NetCDF dataset
   logical, intent(in   ) :: originally_in_define !< Whether the dataset was originally in define mode.
   integer                :: ierr                 !< Result status (nf90_noerr if successful)

   integer :: ierrloc

   ierr = nf90_noerr

   ! Leave the dataset in the same mode as we got it.
   if (originally_in_define) then
      ! Attempt define mode
      ierrloc = nf90_redef(ncid)
      if (ierrloc /= nf90_eindefine) then
         ierr = ierrloc ! Some error occurred
      end if
   else
      ! Attempt data mode
      ierrloc = nf90_enddef(ncid)
      if (ierrloc /= nf90_enotindefine) then
         ierr = ierrloc ! Some error occurred
      end if
   end if

end function ncu_restore_mode


!> Copy all attributes from a variable or dataset into another variable/dataset.
!! Returns:
!     nf90_noerr if all okay, otherwise an error code
!!
!! Note: The variable in the output file must already exist.
function ncu_copy_atts( ncidin, ncidout, varidin, varidout, forbidden_atts, apply_fun) result(ierr)
   use m_alloc
   integer,                               intent(in   ) :: ncidin   !< ID of the input NetCDF file
   integer,                               intent(in   ) :: ncidout  !< ID of the output NetCDF file
   integer,                               intent(in   ) :: varidin  !< ID of the variable in the input file, or NF90_GLOBAL for global attributes.
   integer,                               intent(in   ) :: varidout !< ID of the variable in the output file, or NF90_GLOBAL for global attributes.
   character(len=*),            optional, intent(in   ) :: forbidden_atts(:) !< (Optional) list of forbidden attribute names, will be skipped for copying.
   procedure(ncu_apply_to_att), optional                :: apply_fun !< (Optional) function pointer to facilitate changing the attribute values.

   integer                        :: ierr
   integer                        :: i

   character(len=nf90_max_name)   :: attname
   integer                        :: natts
   integer                        :: atttype   !< attribute data type
   integer                        :: attlen    !< attribute length
   character(len=:), allocatable  :: atttext   !< attribute value

   ierr = -1

   ierr = nf90_inquire_variable( ncidin, varidin, nAtts=natts )
   if ( ierr == nf90_enotvar ) then
      ierr = nf90_inquire( ncidin, nAttributes=natts )
   endif
   if ( ierr /= nf90_noerr ) then
      return
   endif

   do i = 1,natts
      ierr = nf90_inq_attname( ncidin, varidin, i, attname )
      if ( ierr /= nf90_noerr ) then
         return
      endif

      if (present(forbidden_atts)) then
         if (any(forbidden_atts==trim(attname))) then
            cycle
         end if
      end if

      atttype = -1
      ierr = nf90_inquire_attribute(ncidin, varidin, attname, xtype=atttype, len=attlen)
      if (ierr == nf90_noerr .and. atttype == NF90_CHAR .and. present(apply_fun)) then
         ! Special case: do not just copy, but apply a user-provided function to the attribute text first.
         call realloc(atttext, attlen, keepExisting=.false., fill=' ')
         ierr = nf90_get_att(ncidin, varidin, attname, atttext)
         ierr = apply_fun(attname, atttext)
         ierr = nf90_put_att(ncidout, varidout, attname, atttext)
      else
         ! Standard case: copy attribute+value as-is from input dataset to output dataset.
         ierr = nf90_copy_att( ncidin, varidin, attname, ncidout, varidout )
      end if

      if ( ierr /= nf90_noerr ) then
         return
      endif
   enddo

   ierr = nf90_noerr
end function ncu_copy_atts

!> For variable varid in netcdf file ncid append extension to attribute attname 
!! Returns:
!     nf90_noerr if all okay, otherwise an error code
!!
function ncu_append_atts(ncid, varid, attname, extension, separator, check_presence) result(ierr)
   integer                                   :: ierr
   integer,                    intent(in   ) :: ncid           !< ID of the NetCDF file
   integer,                    intent(in   ) :: varid          !< ID of the NetCDF variable, or NF90_GLOBAL for global attributes.
   character(len=*),           intent(in   ) :: attname        !< name of the attribute
   character(len=*),           intent(in   ) :: extension      !< text value to be added to the attribute
   character(len=*), optional, intent(in   ) :: separator      !< (Optional) Separator string to be inserted between existing and extension string (Default: ' ').
   logical,          optional, intent(in   ) :: check_presence !< (Optional) Check whether the extension text is already present, and if so, don't add it again (Default: .false.).

   integer                        :: atttype   !< attribute data type
   character(len=:), allocatable  :: atttext
   integer                        :: attlen    !< attribute length
   character(len=:), allocatable  :: separator_
   logical :: check_presence_
   integer :: ifound
   
   ierr = -1

   if (present(separator)) then
      separator_ = separator ! Intentionally don't trim/adjustl!
   else
      separator_ = ' '
   end if

   if (present(check_presence)) then
      check_presence_ = check_presence
   else
      check_presence_ = .false.
   end if

   
   atttype = 0
   ierr = nf90_inquire_attribute(ncid, varid, attname, xtype=atttype, len=attlen)
   if (ierr == nf90_noerr) then
      if (atttype == NF90_CHAR) then
         allocate(character(len=attlen) :: atttext)
         ierr = nf90_get_att(ncid, varid, attname, atttext)

         if (check_presence_) then
            ifound = index(atttext, trim(extension), back=.true.)
            if (ifound > 0) then
               ! Extension text already present in current attribute, do nothing.
               return
            end if
         end if

         if (attlen > 0) then
            ! Prepare for later appending
            atttext = atttext // separator_
         end if
      else
         ! Attribute already exists, but is not of type char, so cannot add more text to it.
         ierr = IONC_ENOTATT
         return
      end if
   else
      allocate(character(len=0) :: atttext)
   end if

   ! Put the new attribute value (either appended, or afresh)
   ierr = nf90_put_att(ncid, varid, attname, atttext//trim(extension))

end function ncu_append_atts

!> Clones a NetCDF variable definition.
!!
!! The cloned variable will appear under a new name in the (output) file.
!! No data will be copied. Optionally, different standard_name, long_name
!! and units may directly be specified.
function ncu_clone_vardef(ncidin, ncidout, varidin, newname, varidout, &
                          newstdname, newlongname, newunits) result(ierr)
   integer,                    intent(in   ) :: ncidin      !< ID of the input NetCDF file
   integer,                    intent(in   ) :: ncidout     !< ID of the output NetCDF file (can be the same as input NetCDF file)
   integer,                    intent(in   ) :: varidin     !< ID of the variable in the input file.
   character(len=*),           intent(in   ) :: newname     !< Variable name for the new variable.
   integer,                    intent(  out) :: varidout    !< ID of the variable in the output file.
   character(len=*), optional, intent(in   ) :: newstdname  !< New standard_name for the new variable.
   character(len=*), optional, intent(in   ) :: newlongname !< New long_name for the new variable.
   character(len=*), optional, intent(in   ) :: newunits    !< New units for the new variable.


   integer                        :: ierr
   integer                        :: i

   character(len=nf90_max_name)   :: attname
   integer                        :: natts
   integer :: ndims, xtype
   integer, allocatable ::dimids(:)

   ierr = -1

   ierr = nf90_inquire_variable(ncidin, varidin, xtype = xtype, nDims=ndims, nAtts = natts)
   if ( ierr /= nf90_noerr ) then
      return
   endif

   ! Currently skipped properties in nf90_inquire_variable:
   !  logical, optional, intent(out) :: contiguous
   !  integer, optional, dimension(:), intent(out) :: chunksizes
   !  integer, optional, intent(out) :: deflate_level
   !  logical, optional, intent(out) :: shuffle, fletcher32
   !  integer, optional, intent(out) :: endianness

   allocate(dimids(ndims))
   ierr = nf90_inquire_variable(ncidin, varidin, dimids=dimids)

   ierr = nf90_def_var(ncidout, newname, xtype, dimids, varidout)

   ! TODO: AvD: consider copying all attributes
   !do i = 1,natts
   !   ierr = nf90_inq_attname( ncidin, varidin, i, attname )
   !   if ( ierr /= nf90_noerr ) then
   !      return
   !   endif
   !
   !   ierr = nf90_copy_att( ncidin, varidin, attname, ncidout, varidout )
   !   if ( ierr /= nf90_noerr ) then
   !      return
   !   endif
   !enddo

   if (present(newstdname)) then
      ierr = nf90_put_att(ncidout, varidout, 'standard_name', newstdname)
   end if

   if (present(newlongname)) then
      ierr = nf90_put_att(ncidout, varidout, 'long_name', newlongname)
   end if

   if (present(newunits)) then
      ierr = nf90_put_att(ncidout, varidout, 'units', newunits)
   end if

   ierr = nf90_noerr
end function ncu_clone_vardef


!> Compatibility function: returns the fill settings for a variable in a netCDF-3 file.
function ncu_inq_var_fill_int4( ncid, varid, no_fill, fill_value) result(ierr)
   integer,                   intent(in)  :: ncid        !< ID of the NetCDF dataset
   integer,                   intent(in)  :: varid       !< ID of the variable in the data set
   integer,                   intent(out) :: no_fill     !< An integer that will always get 1 (for forward compatibility).
   integer(kind=FourByteInt), intent(out) :: fill_value  !< This will get the fill value for this variable.

   integer :: ierr ! Error status, nf90_noerr = if successful.

   no_fill = 1

   ierr = nf90_get_att(ncid, varid, '_FillValue', fill_value)
   if (ierr /= nf90_noerr) then
      fill_value = nf90_fill_int
      ierr = nf90_noerr
   end if
end function ncu_inq_var_fill_int4


!> Compatibility function: returns the fill settings for a variable in a netCDF-3 file.
function ncu_inq_var_fill_real8( ncid, varid, no_fill, fill_value) result(ierr)
   integer,                   intent(in)  :: ncid        !< ID of the NetCDF dataset
   integer,                   intent(in)  :: varid       !< ID of the variable in the data set
   integer,                   intent(out) :: no_fill     !< An integer that will always get 1 (for forward compatibility).
   real(kind=EightByteReal),  intent(out) :: fill_value  !< This will get the fill value for this variable.

   integer :: ierr ! Error status, nf90_noerr = if successful.
   
   no_fill = 1

   ierr = nf90_get_att(ncid, varid, '_FillValue', fill_value)
   if (ierr /= nf90_noerr) then
      fill_value =  nf90_fill_double
      ierr = nf90_noerr
   end if
end function ncu_inq_var_fill_real8

!> Copy chunking and deflate settings from a variable or dataset into another variable/dataset.
!! Returns:
!     nf90_noerr if all okay, otherwise an error code
!!
!! Note: The variable in the output file must already exist.
function ncu_copy_chunking_deflate( ncidin, ncidout, varidin, varidout, ndx ) result(ierr)
   integer, intent(in)            :: ncidin   !< ID of the input NetCDF file
   integer, intent(in)            :: ncidout  !< ID of the output NetCDF file
   integer, intent(in)            :: varidin  !< ID of the variable in the input file, or NF90_GLOBAL for global attributes.
   integer, intent(in)            :: varidout !< ID of the variable in the output file, or NF90_GLOBAL for global attributes.
   integer, intent(in), optional  :: ndx      !< Number of flow nodes (internal + boundary) for output file

   integer                        :: ierr

   character(len=nf90_max_name)   :: name
   integer                        :: storage, ndims, shuffle, deflate, deflate_level
   integer, allocatable           :: chunksizes(:)

   ierr = -1

   if (varidin /= NF90_GLOBAL) then
      !
      ! copy chuncking settings, if available
      !
      ierr = nf90_inquire_variable( ncidin, varidin, nDims=ndims, name=name )
      if (ierr == nf90_noerr .and. ndims > 0) then
         allocate(chunksizes(ndims))
         ierr = nf90_inq_var_chunking(ncidin, varidin, storage, chunksizes)
         if (ierr == 0 .and. storage == nf90_chunked) then
            !
            ! chuncking is on for this variable
            !
            if (present(ndx)) then
               !
               ! first dimension is ndx, update with global ndx
               !
               chunksizes(1) = ndx
            endif
            ierr = nf90_def_var_chunking(ncidout, varidout, storage, chunksizes)
            if (ierr /= 0) write(*,*) 'nf90_def_var_chunking failed for var ', trim(name)
         endif
         deallocate(chunksizes)
      endif
      !
      ! copy deflation settings, if available
      !
      ierr = nf90_inq_var_deflate(ncidin, varidin, shuffle, deflate, deflate_level)
      if (ierr == nf90_noerr .and. deflate == 1) then
         ierr = nf90_def_var_deflate(ncidout, varidout, shuffle, deflate, deflate_level)
         if (ierr /= 0) write(*,*) 'nf90_def_var_deflate failed for var ', trim(name)
      endif
   endif

   ierr = nf90_noerr
end function ncu_copy_chunking_deflate

!> Get the STRING-attribute of a variable from the nc-file.
function ncu_get_att(ncid, varid, att_name, att_value) result(status)
   integer,                       intent(in   ) :: ncid         !< NetCDF dataset id, should be already open.
   integer,                       intent(in   ) :: varid        !< NetCDF variable id (1-based).
   character(len=*),              intent(in   ) :: att_name     !< The name of the attribute.
   character(len=:), allocatable, intent(inout) :: att_value    !< The value of the attribute (unchanged when an error occurred).
   integer                                      :: status       !< Result status, ug_noerr if successful.

   integer :: att_value_len       !< Length of the attribute value on the netCDF file.
   integer :: istat

   status    = nf90_inquire_attribute(ncid, varid, att_name, len=att_value_len)
   if ( status == nf90_noerr ) then
      if (allocated(att_value)) then
         deallocate(att_value, stat = istat)
         if (istat /= 0) then
            status = istat
            return
         end if
      end if 

      allocate( character(len=att_value_len) :: att_value, stat = istat )
      if (istat /= 0) then
         status = istat
         return
      end if
      status = nf90_get_att(ncid, varid, att_name, att_value)
   end if

end function ncu_get_att

!> Gets all NetCDF-attributes for a given variable.
!!
!! This function is non-UGRID-specific: only used to read grid mapping variables.
!! @see ug_put_var_attset
function ncu_get_var_attset(ncid, varid, attset) result(ierr)
   use coordinate_reference_system
   use ionc_constants

   integer,                         intent(in)  :: ncid      !< NetCDF dataset id
   integer,                         intent(in)  :: varid     !< NetCDF variable id (1-based).
   type(nc_attribute), allocatable, intent(out) :: attset(:) !< Resulting attribute set.
   integer                                      :: ierr      !< Result status (UG_NOERR==NF90_NOERR) if successful.

   character(len=64) :: attname
   character(len=:), allocatable     :: tmpstr !< The name of the network used by the mesh
   integer :: i, j, natts, atttype, attlen, nlen

   ierr = IONC_NOERR
   allocate( character(len=0) :: tmpstr )

   ierr = nf90_inquire_variable(ncid, varid, natts = natts)
   if (ierr /= nf90_noerr) then
      goto 888
   end if

   if (allocated(attset)) deallocate(attset)
   allocate(attset(natts), stat=ierr)

   do i = 1,natts
      ierr = nf90_inq_attname(ncid, varid, i, attname)    ! get attribute name
      ierr = nf90_inquire_attribute(ncid, varid, trim(attname), xtype = atttype, len=attlen) ! get other attribute information

      select case(atttype)
      case(NF90_CHAR)
         tmpstr = ''
         ierr = ncu_get_att(ncid, varid, attname, tmpstr)   
         
         allocate(attset(i)%strvalue(attlen))
         nlen = min(len(tmpstr), attlen)
         do j=1,nlen
            attset(i)%strvalue(j) = tmpstr(j:j)
         end do
      case(NF90_INT)
         allocate(attset(i)%intvalue(attlen))
         ierr = nf90_get_att(ncid, varid, attname, attset(i)%intvalue)
      case(NF90_FLOAT)
         allocate(attset(i)%fltvalue(attlen))
         ierr = nf90_get_att(ncid, varid, attname, attset(i)%fltvalue)
      case(NF90_DOUBLE)
         allocate(attset(i)%dblvalue(attlen))
         ierr = nf90_get_att(ncid, varid, attname, attset(i)%dblvalue)
      case default
         ! NF90_BYTE
         ! NF90_SHORT
         ncu_messagestr = 'ncu_get_var_attset: error for attribute '''//trim(attname)//'''. Data types byte/short not implemented.'
         ierr = IONC_ENOTAVAILABLE
         goto 888
      end select
      attset(i)%attname = attname
      attset(i)%xtype   = atttype
      attset(i)%len     = attlen
   end do
   deallocate(tmpstr)

   return ! Return with success

888 continue

end function ncu_get_var_attset

!> Puts a set of NetCDF-attributes onto a given variable.
!!
!! This function is non-UGRID-specific: only used to write grid mapping variables.
!! @see ug_get_var_attset
function ncu_put_var_attset(ncid, varid, attset) result(ierr)
   use coordinate_reference_system
   use ionc_constants

   integer,             intent(in)  :: ncid      !< NetCDF dataset id
   integer,             intent(in)  :: varid     !< NetCDF variable id (1-based).
   type(nc_attribute),  intent(in)  :: attset(:) !< Attribute set to be put into the variable.
   integer                          :: ierr      !< Result status (UG_NOERR==NF90_NOERR) if successful.

   character(len=1024) :: tmpstr
   integer :: i, j, natts, nlen

   ierr = IONC_NOERR

   natts = size(attset)

   do i = 1,natts
      select case(attset(i)%xtype)
      case(NF90_CHAR)
         tmpstr = ' '
         nlen = min(len(tmpstr), attset(i)%len)
         do j=1,nlen
            tmpstr(j:j) = attset(i)%strvalue(j)
         end do

         ierr = nf90_put_att(ncid, varid, attset(i)%attname, tmpstr)
      case(NF90_INT)
         ierr = nf90_put_att(ncid, varid, attset(i)%attname, attset(i)%intvalue(1:attset(i)%len))
      case(NF90_FLOAT)
         ierr = nf90_put_att(ncid, varid, attset(i)%attname, attset(i)%fltvalue(1:attset(i)%len))
      case(NF90_DOUBLE)
         ierr = nf90_put_att(ncid, varid, attset(i)%attname, attset(i)%dblvalue(1:attset(i)%len))
      case default
         ! NF90_BYTE
         ! NF90_SHORT
         ncu_messagestr = 'ug_put_var_attset: error for attribute '''//trim(attset(i)%attname)//'''. Data types byte/short not implemented.'
         ierr = IONC_ENOTAVAILABLE
      end select
   end do

end function ncu_put_var_attset


!> Inquire for a NetCDF variable ID based on an attribute in another variable.
!! For example: mesh2d:face_node_connectivity
function ncu_att_to_varid(ncid, varid, attname, id) result(ierr)
   use ionc_constants
   
   integer         , intent(in   ) :: ncid    !< NetCDF dataset ID
   integer         , intent(in   ) :: varid   !< NetCDF variable ID from which the attribute will be gotten (1-based).
   character(len=*), intent(in   ) :: attname !< Name of attribute in varid that contains the variable name.
   integer         , intent(  out) :: id      !< NetCDF variable ID that was found.
   integer                         :: ierr    !< Result status. UG_NOERR if successful.

   character(len=:), allocatable  :: varname

   ierr = IONC_NOERR
   allocate( character(len=0) :: varname )

   varname = ''
   ierr = ncu_get_att(ncid, varid, attname, varname)
   if (ierr /= nf90_noerr) then
      ierr = IONC_ENOTATT
      goto 999
   end if
   ierr = nf90_inq_varid(ncid, trim(varname), id)
   if (ierr /= nf90_noerr) then
      ierr = IONC_ENOTVAR
      goto 999
   end if

   ! Return with success
   deallocate(varname)
   return

999 continue
   ! An error occurred, keep ierr nonzero and set undefined id.
   deallocate(varname)
   id = -1         ! undefined id
end function ncu_att_to_varid


!> Inquire for a NetCDF dimension ID based on an attribute in another variable.
!! For example: mesh2d:edge_dimension
function ncu_att_to_dimid(ncid, varid, attname, id) result(ierr)
   use ionc_constants

   integer         , intent(in   ) :: ncid    !< NetCDF dataset ID
   integer         , intent(in   ) :: varid   !< NetCDF variable ID from which the attribute will be gotten (1-based).
   character(len=*), intent(in   ) :: attname !< Name of attribute in varid that contains the dimension name.
   integer         , intent(  out) :: id      !< NetCDF dimension ID that was found.
   integer                         :: ierr    !< Result status. UG_NOERR if successful.

   character(len=nf90_max_name)    :: varname

   ierr = IONC_NOERR

   varname = ''
   ierr = nf90_get_att(ncid, varid, attname, varname)
   if (ierr /= nf90_noerr) then
      ierr = IONC_ENOTATT
      goto 999
   end if

   ierr = nf90_inq_dimid(ncid, trim(varname), id)
   if (ierr /= nf90_noerr) then
      ierr = IONC_ENOTDIM
      goto 999
   end if
   ! Return with success
   return

999 continue
   ! An error occurred, keep ierr nonzero and set undefined id.
   id = -1         ! undefined id

end function ncu_att_to_dimid


!copy the variable attributes
function ncu_copy_var_atts( ncidin, ncidout, varidin, varidout ) result(ierr)

    integer, intent(in)            :: ncidin, ncidout, varidin, varidout
    integer                        :: ierr
    integer                        :: i
    character(len=nf90_max_name)   :: attname
    integer                        :: natts
    integer                        :: attvalue

    ierr = -1
    ierr = nf90_inquire_variable( ncidin, varidin, nAtts=natts )
    if ( ierr == nf90_enotvar ) then
        ierr = nf90_inquire( ncidin, nAttributes=natts )
    endif
    if ( ierr /= nf90_noerr ) then
        return
    endif

    do i = 1,natts
        ierr = nf90_inq_attname( ncidin, varidin, i, attname )
        if ( ierr /= nf90_noerr ) then
            return
        endif

        ierr = nf90_copy_att( ncidin, varidin, attname, ncidout, varidout )
        if ( ierr /= nf90_noerr ) then
            return
        endif
    enddo

end function ncu_copy_var_atts


end module netcdf_utils
