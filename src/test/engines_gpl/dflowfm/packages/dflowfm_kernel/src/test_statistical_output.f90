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

module test_statistical_output
   use ftnunit
   use m_statistical_output

   implicit none
   double precision, parameter :: eps        = 1.0d-10
   logical                     :: hysteresis = .false.

contains

subroutine tests_statistical_output
   call test( test_realloc,                    'Tests the realloc subroutine' )
end subroutine tests_statistical_output

!> Please see "cross sections.xlsx" for the caluclation of the results
subroutine test_realloc
   use m_statistical_output
  
   type(t_statistical_output_set) :: stat_output_set
   integer :: i, n
   double precision, allocatable, target :: s1(:)


   
   call realloc(stat_output_set)
   
   ! The size of the output set should now be equal to "growsBy":
   call assert_equal(stat_output_set%size, stat_output_set%growsBy, "The size of the outputset is not correct" )
   call assert_equal(size(stat_output_set%statout), stat_output_set%growsBy, "The allocated size of the outputset is not correct" )
   call assert_equal(stat_output_set%count, 0, "The count of the outputset must be equal to 0 at this time")
   
   ! set some variable in the outputset array, and check these after reallocation
   stat_output_set%count = 10
   n = 15
   allocate(s1(n))
   do i = 1, n
      s1(i) = 0.1d0*i**2
   enddo
   stat_output_set%statout(10)%location_specifier = 2
   stat_output_set%statout(10)%name = "s1"
   stat_output_set%statout(10)%input => s1

   call realloc(stat_output_set)   

   ! The size of the output set should now be equal to 2*"growsBy":
   call assert_equal(stat_output_set%size  , 2*stat_output_set%growsBy, "The size of the outputset is not correct" )
   call assert_equal(size(stat_output_set%statout), 2*stat_output_set%growsBy, "The allocated size of the outputset is not correct" )
   call assert_equal(stat_output_set%count, 10, "The count of the outputset must be equal to 10 at this time")
   call assert_equal(stat_output_set%statout(10)%location_specifier,  2, "The location specifier must be equal to 2")
   call assert_equal(trim(stat_output_set%statout(10)%name), "s1", "The name of the array must be 's1'")
   call assert_equal(size(stat_output_set%statout(10)%input), n, "The allocated size of s1 is not equal to 15" )
   do i = 1, n
      call assert_comparable(stat_output_set%statout(10)%input(i), 0.1d0*i**2, eps, "The contents of array s1 is not correct")
   enddo

end subroutine test_realloc
end module test_statistical_output
