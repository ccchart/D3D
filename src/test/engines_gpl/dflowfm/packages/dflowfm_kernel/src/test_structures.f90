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
module test_structures
   use ftnunit
   use precision

   implicit none
   real(fp), parameter :: eps = 1.0e-6_fp

contains
!
!
!==============================================================================
subroutine tests_structure
   call test( test_structurefile, 'Tests the reading of structure control via a iniField file.' )
end subroutine tests_structure
!
!
!==============================================================================
subroutine test_structurefile
   use dfm_error
   use ifport, only: CHANGEDIRQQ
   use properties
   use mathconsts, only: eps_hp
   use m_strucs

   integer                   :: istat   
   integer                   :: ierr
   integer                   :: numblocks
   integer                   :: i, j, k

   type(TREE_DATA), pointer        :: str_ptr    !< Property tree as read from structures.ini
   type(TREE_DATA), pointer        :: block_ptr  !< Property tree containing one block of generalstructure data in ini file   
   character(len=256)              :: key        !< property key.
   character(len=256)              :: strvalue   !< Returned string value for requested property key.
   character(len=256)              :: idvalue    !< block idententy string 
   character(len=32)               :: structurefile !< filename with testdata.
   double precision                :: dblvalue   !< Returned scalar double value for requested property key, IF possible.
   double precision                :: tmpvalue   !< Returned scalar double value for requested property key, IF possible.
   logical                         :: is_double  !< Tells whether the found value could be parsed into a scalar double value.
   logical                         :: success    !< Whether value was read successfully or not.

   structurefile = "structures.ini"
   istat = CHANGEDIRQQ("structures")

   call tree_create(trim(structurefile), str_ptr)
   call prop_inifile(structurefile, str_ptr, ierr)
   call assert_equal(ierr, DFM_NOERR, 'Error reading structure file.') 
   istat = CHANGEDIRQQ("..")

   numblocks = tree_num_nodes(str_ptr)
   call assert_equal(numblocks,5,'File structures/structures.ini is expected to contain 5 blocks.')
   
   ! read required properties from first block
   call read_required_property(str_ptr%child_nodes(1), "id", strvalue, dblvalue, is_double, 'first block structures.ini', success)
   call assert_equal(success, .TRUE., "Something wrong reading 'id'.") 
   call assert_equal(is_double, .FALSE., "Expected a string.") 
   call assert_equal(strvalue,"full_block", "Another value for 'id' was expected.")

   call read_required_property(str_ptr%child_nodes(1), "id_nonexistent", strvalue, dblvalue, is_double, 'first block structures.ini', success)
   call assert_equal(success, .FALSE., "Something wrong reading 'id_nonexistent'.") 

   do i = 1, numblocks
      block_ptr => str_ptr%child_nodes(i)%node_ptr
            
      call read_required_property(str_ptr%child_nodes(i), "type", strvalue, dblvalue, is_double, 'block structures.ini', success)
      call assert_equal(success, .TRUE., "Something wrong reading 'type'.") 
      call assert_equal(strvalue,'generalstructure','Type "generalstructures" was expected.')

      call read_required_property(str_ptr%child_nodes(i), "id", idvalue, dblvalue, is_double, 'block structures.ini', success)
      call assert_equal(success, .TRUE., "Something wrong reading 'id'.") 

      ! for each block read GateLowerEdgeLevel either as string or as double
      call read_required_property(str_ptr%child_nodes(i), 'GateLowerEdgeLevel', strvalue, dblvalue, is_double, 'first block structures.ini', success)
      if (success) then 
         select case(trim(idvalue))
         case ('full_block')
            call assert_equal(is_double, .TRUE., "Block 'full_block' in structures.ini: expected a value.") 
            call assert_comparable(dblvalue, 1.d0, eps_hp, "Read GateLowerEdgeLevel as a value.")
         case ('local')
            call assert_equal(is_double, .FALSE., "Block 'local': expected a string.")          
            call assert_equal(strvalue,'filename.tim','Unexpected string.')
         case ('relative')
            call assert_equal(is_double, .FALSE., "Block 'relative': expected a string.")          
            call assert_equal(strvalue,'../filename.tim','Unexpected string.')
         case ('windows')
            call assert_equal(is_double, .FALSE., "Block 'windows': expected a string.")          
            call assert_equal(strvalue,'c:\dirname\0000\filename.tim','Unexpected string.')
         case ('linux')
            call assert_equal(is_double, .FALSE., "Block 'linux' in structures.ini: expected a string.") 
            call assert_equal(strvalue, '/home/UNST_5890/filename.tim', "/home/UNST_5890/filename.tim: Unexpected string.")    
         end select
      else
         ! raise an error message when .not. success
         call assert_equal(success, .TRUE., "Something wrong reading 'GateLowerEdgeLevel'.")       
      endif
   enddo
   
end subroutine test_structurefile

end module test_structures
