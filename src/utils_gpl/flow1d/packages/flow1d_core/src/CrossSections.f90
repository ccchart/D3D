!> Define cross-sections
module m_CrossSections
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

   use MessageHandling
   use m_alloc
   use M_newcross
   use m_branch
   use m_GlobalParameters
   use m_hash_search
   use string_module

   implicit none

   private

   public realloc
   public dealloc
   public AddCrossSection
   public AddCrossSectionDefinition
   public CalcConveyance
   public createTablesForTabulatedProfile
   public fill_hashtable
   public getBob
   public GetCriticalDepth
   public GetCrossSectionsCount
   public GetCrossType
   public GetCSParsFlow
   public GetCSParsTotal
   public getGroundLayer
   public getHighest1dLevel
   public GetTabFlowSectionFromTables
   public interpolateSummerDike
   public setGroundLayerData
   public SetParsCross
   public useBranchOrders
   
   interface fill_hashtable
      module procedure fill_hashtable_csdef
   end interface
   
   interface getHighest1dLevel
      module procedure getHighest1dLevelInterpolate
      module procedure getHighest1dLevelSingle
   end interface getHighest1dLevel
   
   interface getBob
      module procedure getBobInterpolate
      module procedure getBobSingle
   end interface getBob
   
   interface GetCSParsFlow
      module procedure GetCSParsFlowInterpolate
      module procedure GetCSParsFlowCross
   end interface GetCSParsFlow

   interface GetCSParsTotal
      module procedure GetCSParsTotalInterpolate
      module procedure GetCSParsTotalCross
   end interface GetCSParsTotal
   
   !> Add cross section definition for complete tabulated or YZ-cross-sections.For YZ-profiles also just a location definition is possible
   interface AddCrossSectionDefinition
      module procedure AddRoundCrossSectionDefinition
      module procedure AddTabCrossSectionDefinition
      module procedure AddYZCrossSectionDefinition
   end interface

   interface AddCrossSection
      module procedure AddCrossSectionByVariables
      module procedure AddCrossSectionByCross
   end interface AddCrossSection

   interface useBranchOrders
      module procedure useBranchOrdersCrs
   end interface useBranchOrders

   !> Realloc memory cross-section definition or cross-sections
   interface realloc
      module procedure reallocCSDefinitions
      module procedure reallocCSDefinitionsSize
      module procedure reallocCrossSections
   end interface

   !> Free the memory of cross-section definition or cross-sections
   interface dealloc
      module procedure deallocCrossDefinition
      module procedure deallocCSDefinitions
      module procedure deallocCrossSection
      module procedure deallocCrossSections
   end interface dealloc

   integer, public, parameter :: CS_STOR_LOSS      =    1
   integer, public, parameter :: CS_STOR_RESERVOIR =    0
   integer, public, parameter :: CS_TABULATED      =    1
   integer, public, parameter :: CS_CIRCLE         =    2
   integer, public, parameter :: CS_EGG            =    3
   integer, public, parameter :: CS_RECTANGLE      =    4 ! TODO, placeholder, someone forgot to check in?
   integer, public, parameter :: CS_TRAPEZIUM      =    5
   integer, public, parameter :: CS_YZ_PROF        =   10
   integer, public, parameter :: CS_TYPE_NORMAL    =    1 !< Ordinary total area computation, with possible Preisman lock on top
   integer, public, parameter :: CS_TYPE_PREISMAN  =    2 !< Ordinary total area computation, with possible Preisman lock on top
   integer, public, parameter :: CS_TYPE_PLUS      =    3 !< Total area for only the expanding part of the cross section (Nested Newton method)
   integer, public, parameter :: CS_TYPE_MIN       =    4 !< Total area for only the narrowing part of the cross section (Nested Newton method)
   
   integer, parameter         :: eps = 1d-3               !< accuracy parameter for determining wetperimeter == 0d0

   character (len=13), public, dimension(10) :: CSTypeName = (/&
   '1:TABULATED ',&
   '2:CIRCLE    ',&
   '3:EGG       ',&
   '            ',&
   '            ',&
   '            ',&
   '            ',&
   '            ',&
   '            ',&
   '10:YZ       '/)

   type  :: t_groundlayer
      logical          :: used      = .false.
      double precision :: thickness = 0.0d0
      double precision :: area      = 0.0d0
      double precision :: perimeter = 0.0d0
      double precision :: width     = 0.0d0
   end type t_groundlayer

   type, public  :: t_summerdike
      double precision :: crestLevel = 0.0d0     !< Crest height of the summerdike
      double precision :: baseLevel  = 0.0d0     !< Base level of the summerdike
      double precision :: flowArea   = 0.0d0     !< Flow area of the summerdike
      double precision :: totalArea  = 0.0d0     !< Total area, i.e. Flow area + storage area
      logical          :: hysteresis = .true.    !< Remember Switch for Total Area (Hysteresis)
   end type t_summerdike

   !> Derived type for defining one cross-section type
   type, public :: t_CSType
        character(IdLen)            :: id                      !< unique identification
         
        integer                     :: crossType        ! cross(3,*): Type of cross section
                                                        !      1: Tabulated
                                                        !      2: Circular
                                                        !      3: Egg shaped
                                                        !      4: Rectangular
                                                        !     10: Y-Z cross section
                                                        !     11: A-symmetrical
       integer                      :: levelsCount = 0  !     number of levels in tabulated cross section
       integer                      :: reference

       !*** data for tabulated cross sections ***
       double precision, allocatable, dimension(:)   :: height       !< Specified height of the cross-section
       double precision, allocatable, dimension(:)   :: flowWidth    !< Specified flow width of the cross-section
       double precision, allocatable, dimension(:)   :: totalWidth   !< Specified total width of the cross-section
       logical                                       :: closed       !< Flag to determine if cross-section is closed

       ! Pre-Calculated Table for Tabulated/River Profiles
       double precision, allocatable, dimension(:,:) :: af_sub        !< Flow Areas for Sub-Sections (Main, FP1 and FP2)
       double precision, allocatable, dimension(:,:) :: width_sub     !< Width for Sub-Sections (Main, FP1 and FP2)
       double precision, allocatable, dimension(:,:) :: perim_sub     !< Wetted Perimeter for Sub-Sections (Main, FP1 and FP2)
       double precision, allocatable, dimension(:)   :: flowArea      !< Flow Areas
       double precision, allocatable, dimension(:)   :: wetPerimeter  !< Wet Perimeters
       double precision, allocatable, dimension(:)   :: totalArea     !< Total Areas
       double precision, allocatable, dimension(:)   :: area_min      !< Area for a narrowing part of a cross section (Nested Newton)
       double precision, allocatable, dimension(:)   :: width_min     !< Width for a narrowing part of a cross section (Nested Newton)
       
       !--- additional data for river profiles
       double precision                         :: plains(3) = 0.0d0  !< 1: main channel, 2: floopplain 1, 3: floodplain 2
       integer                                  :: plainsLocation(3)
       !--- information for Summerdikes
       type(t_summerdike),pointer               :: summerdike => null()

       !*** data for yz cross sections
       double precision, allocatable            :: y(:)       !< tranversal co-ordinate
       double precision, allocatable            :: z(:)       !< z-co-ordinate
       integer                                  :: storLevelsCount = 0 !< Number of actual storage levels
       double precision, allocatable            :: storLevels(:)       !< Storage levels
       double precision, allocatable            :: YZstorage(:)        !< Storage levels
       integer                                  :: storageType     = 0 !< Storage type (default no storage level)

       !** Circle and egg profile data
       double precision                         :: diameter

       !*** ground layer data
       type(t_groundlayer),pointer              :: groundlayer => null()
       
       ! Friction Data not needed in Cache, already processed into Cross-Sections, which are cached
       integer                                  :: frictionSectionsCount = 0 !< Number of actual friction sections
       character(IdLen), allocatable            :: frictionSectionID(:)      !< Friction Section Identification
       double precision, allocatable            :: frictionSectionFrom(:)    !<
       double precision, allocatable            :: frictionSectionTo(:)      !<

   end type t_CSType

   !> Derived type to store the cross-section definition, grouped according to the types
   type, public :: t_CSDefinitionSet
      integer                                   :: size    = 0    !< Size of the cross-section set
      integer                                   :: growsBy = 2000  !< Increment for the cross-section type
      integer                                   :: count   = 0    !< Actual number of cross-sections
      type(t_CSType), pointer, dimension(:)     :: CS             !< Cross-section derived type
      type(t_hashlist)                          :: hashlist       !< hash table for searching cross section definition on id
   end type

   !---------------------------------------------------------
   !> Derived type to store the cross-section
   type, public :: t_CrossSection

       character (len=IdLen)        :: csid              !< User-defined ID of this crossection
       integer                      :: crossIndx        !< Index to check on same cross-section
                                                        !! To Prevent unecessary interpolations.
       logical                      :: IsCopy = .false. !< Flag to determine if Cross-Section is a copy
       
       integer                      :: crossType        !< cross(3,*): Type of cross section

       integer                      :: branchid = -1    !< integer branch id
       double precision             :: location         !< offset in meters along reach
       double precision             :: bedLevel
       double precision             :: shift
       double precision             :: surfaceLevel
       double precision             :: charHeight
       double precision             :: charWidth
       logical                      :: closed
       integer                      :: bedFrictionType
       double precision             :: bedFriction
       integer                      :: groundFrictionType
       double precision             :: groundFriction
       integer                      :: iTabDef           !< temporary item to save ireference number to cross section definition
                                                         !! necessary for reallocation of arrays
       type(t_CSType), pointer      :: tabDef => null()
       type(t_crsu), pointer        :: convTab => null()
       
       integer                                  :: frictionSectionsCount = 0 !< Number of actual friction sections
       character(IdLen), allocatable            :: frictionSectionID(:)      !< Friction Section Identification
       double precision, allocatable            :: frictionSectionFrom(:)    !<
       double precision, allocatable            :: frictionSectionTo(:)      !<
       integer, allocatable                     :: frictionTypePos(:)        !< Friction type for positive flow direction
       double precision, allocatable            :: frictionValuePos(:)       !< Friction value for positive flow direction
       integer, allocatable                     :: frictionTypeNeg(:)        !< Friction type for negative flow direction
       double precision, allocatable            :: frictionValueNeg(:)       !< Friction value for negative flow direction

    end type

   !> Derived type to store the cross-section set
   type, public :: t_CrossSectionSet
      integer                                                :: maxNumberOfConnections=0 ! maximum nr of connections to a node
      integer                                                :: Size = 0      !< Actual size of cross-section set
      integer                                                :: growsBy = 2000 !< Increment for cross-section set
      integer                                                :: Count= 0      !< Actual number of cross-section sets
      type(t_CrossSection), pointer, dimension(:)            :: cross         !< Current cross-section
      integer, pointer, dimension(:)                         :: CrossSectionIndex =>null() !< When sorting cross sections, the index of structure related cross sections must be reindexed.
   end type t_CrossSectionSet

contains

!> Free the memory used by a cross-section definition
subroutine deallocCrossDefinition(CrossDef)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CSType)                            :: CrossDef !< Current cross-section definition

   ! Local variables
   if (allocated(CrossDef%height))                   deallocate(CrossDef%height)
   if (allocated(CrossDef%flowWidth))                deallocate(CrossDef%flowWidth)
   if (allocated(CrossDef%totalWidth))               deallocate(CrossDef%totalWidth)

   if (allocated(CrossDef%af_sub))                   deallocate(CrossDef%af_sub)
   if (allocated(CrossDef%width_sub))                deallocate(CrossDef%width_sub)
   if (allocated(CrossDef%perim_sub))                deallocate(CrossDef%perim_sub)
   if (allocated(CrossDef%flowArea))                 deallocate(CrossDef%flowArea)
   if (allocated(CrossDef%wetPerimeter))             deallocate(CrossDef%wetPerimeter)
   if (allocated(CrossDef%totalArea))                deallocate(CrossDef%totalArea)
   if (allocated(CrossDef%area_min))                 deallocate(CrossDef%area_min)
   if (allocated(CrossDef%width_min))                deallocate(CrossDef%width_min)

   if (allocated(CrossDef%y))                        deallocate(CrossDef%y)
   if (allocated(CrossDef%z))                        deallocate(CrossDef%z)
   if (allocated(CrossDef%frictionSectionID))        deallocate(CrossDef%frictionSectionID)
   if (allocated(CrossDef%frictionSectionFrom))      deallocate(CrossDef%frictionSectionFrom)
   if (allocated(CrossDef%frictionSectionTo))        deallocate(CrossDef%frictionSectionTo)
   if (associated(CrossDef%summerdike))              deallocate(CrossDef%summerdike)
   if (associated(CrossDef%groundlayer))             deallocate(CrossDef%groundlayer)
   
   CrossDef%summerdike => null()
   CrossDef%groundlayer => null()
   
end subroutine deallocCrossDefinition

!> Free the memory used by a cross-section definition set
subroutine deallocCSDefinitions(CSdef)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CSDefinitionSet)   :: CSdef !< Cross-section definition set

   ! Local variables
   integer length, i

   ! Program code
   if (associated(CSdef%CS)) then
      length = size(CSdef%CS)
      do i = 1, length
         call dealloc(CSdef%CS(i))
      enddo
      deallocate(CSdef%CS)
      CSdef%CS =>null()
      CSDef%size            = 0
      CSDef%count           = 0
   endif
   call dealloc(CSDef%hashlist)
end subroutine deallocCSDefinitions



!> Increase the memory used by a cross-section definition
subroutine reallocCSDefinitionsSize(CSDef,growsBy)
   implicit none
   type(t_CSDefinitionSet)   :: CSdef !< Current cross-section definition
   integer, intent(in)       :: growsBy
   integer                   :: old_growsBy
   old_growsBy = CSdef%growsBy
   CSdef%growsBy = growsBy
   call reallocCSDefinitions(CSDef)
   CSdef%growsBy = old_growsBy
   return
end subroutine 


!> Increase the memory used by a cross-section definition
subroutine reallocCSDefinitions(CSDef)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CSDefinitionSet)   :: CSdef !< Current cross-section definition
   integer                   :: ierr

   ! Local variables
   ! Local variables
   type(t_CSType), pointer, dimension(:)    :: oldDefs

   ! Program code

   if (CSDef%Size > 0) then
      oldDefs=>CSDef%CS
   endif

   if (CSDef%growsBy <=0) then
      CSDef%growsBy = 200
   endif
   allocate(CSDef%CS(CSDef%Size+CSDef%growsBy),stat=ierr)
   call aerr('CSDef%CS(CSDef%Size+CSDef%growsBy)',ierr,CSDef%Size+CSDef%growsBy)

   if (CSDef%size > 0) then
      CSDef%CS(1:CSDef%size) = oldDefs(1:CSDef%size)
      deallocate(oldDefs)
   endif
   CSDef%Size = CSDef%Size+CSDef%growsBy
end subroutine

integer function GetCrossType(string)

   character(len=*) :: string

   call str_lower(string)

   select case(trim(string))
      case ('tabulated')
         GetCrossType = CS_TABULATED
      case ('trapezium')
         GetCrossType = CS_TRAPEZIUM
      case ('circle')
         GetCrossType = CS_CIRCLE
      case ('egg')
         GetCrossType = CS_EGG
      case ('yz')
         GetCrossType = CS_YZ_PROF
      case ('rectangle')
         GetCrossType = CS_RECTANGLE
      case ('xyz')
         GetCrossType = CS_YZ_PROF
      case ('ellipse')
         GetCrossType = CS_TABULATED
      case ('arch')
         GetCrossType = CS_TABULATED
      case ('cunette')
         GetCrossType = CS_TABULATED
      case ('steelcunette')
         GetCrossType = CS_TABULATED
      case default
         GetCrossType = -1
   end select

end function GetCrossType


!> integer function AddRoundCrossSectionDefinition adds an CSType (cross section definition type) element to the
!! set CSDefinitions of type t_CSDefinitionSet and initializes the added element using the parameters passed.
integer function AddRoundCrossSectionDefinition(CSDefinitions, id, diameter, shape, groundlayerUsed, groundlayer)

   implicit none
   ! Input/output parameters
   type(t_CSDefinitionSet), intent(inout) :: CSDefinitions
   character(len=*), intent(in) :: id
   double precision, intent(in)     :: diameter
   integer, intent(in)              :: shape                    !< shape can be CS_CIRCLE or CS_EGG
   logical, intent(in)              :: groundlayerUsed          !< flag indicating whether a ground layer is used
   double precision, intent(in)     :: groundLayer              !< Thickness of the groundlayer

   ! local variables
   integer     :: i

   ! program code
   CSDefinitions%count = CSDefinitions%count+1

   if (CSDefinitions%count > CSDefinitions%size) then   ! skipping real allocation
      call realloc(CSDefinitions)
   endif

   i = CSDefinitions%count

   if (shape /= CS_CIRCLE .and. shape /= CS_EGG) then
      call setmessage(LEVEL_ERROR, 'INTERNAL ERROR: incorrect type of cross section.')
      return
   endif

   CSDefinitions%CS(i)%id           = id 
   CSDefinitions%CS(i)%crossType    = shape
   CSDefinitions%CS(i)%diameter     = diameter
   CSDefinitions%CS(i)%closed       = .true.
   
   allocate(CSDefinitions%CS(i)%groundlayer)
   if (groundlayerUsed) anyGroundLayer       = .true.
   CSDefinitions%CS(i)%groundlayer%used      = groundlayerUsed
   CSDefinitions%CS(i)%groundlayer%thickness = groundlayer
!
   AddRoundCrossSectionDefinition = i

end function AddRoundCrossSectionDefinition


!> integer function AddTabCrossSectionDefinition adds an CSType (cross section definition type) element to the
!! set CSDefinitions of type t_CSDefinitionSet and initializes the added element using the parameters passed.

integer function AddTabCrossSectionDefinition(CSDefinitions , id, numLevels, level, flowWidth, totalWidth, plains,  &
                                            & crestLevel, baseLevel, flowArea, totalArea, closed,              &
                                            & groundlayerUsed, groundlayer)
   ! Modules
   use m_alloc
   implicit none

   ! Input/output parameters
   type(t_CSDefinitionSet), intent(inout) :: CSDefinitions
   character(len=*), intent(in) :: id
   integer, intent(in)              :: numLevels                !< length of the arrays: level, flowWidth and totalWidth
   double precision, intent(in)     :: level(numLevels)         !< heights at which widths are defined
   double precision, intent(in)     :: flowWidth(numLevels)     !< flow width as a function of level
   double precision, intent(in)     :: totalWidth(numLevels)    !< total width as a function of level
   double precision, intent(in)     :: plains(3)                !< width of main section, flood plane 1 and flood plane 2
   double precision, intent(in)     :: crestLevel               !< level height of summerdike
   double precision, intent(in)     :: baseLevel                !< level height of base of summerdike channel
   double precision, intent(in)     :: flowArea                 !< effective flow area of summer dike channel
   double precision, intent(in)     :: totalArea                !< the total area of summer dike channel
   logical, intent(in)              :: closed                   !< flag indicating whether a cross section is closed
   logical, intent(in)              :: groundlayerUsed          !< flag indicating whether a ground layer is used
   double precision, intent(in)     :: groundLayer              !< Thickness of the groundlayer

   ! Local variables
   integer i, j, length
!
   ! Program code
   CSDefinitions%count = CSDefinitions%count+1

   if (CSDefinitions%count > CSDefinitions%size) then
      call realloc(CSDefinitions)
   endif
!
   i = CSDefinitions%count
   if (closed) then
      length = numLevels+1
   else
      length = numLevels
   endif
   call realloc(CSDefinitions%CS(i)%height, length)
   call realloc(CSDefinitions%CS(i)%flowWidth, length)
   call realloc(CSDefinitions%CS(i)%totalWidth, length)

   call realloc(CSDefinitions%CS(i)%af_sub, 3, length)
   call realloc(CSDefinitions%CS(i)%width_sub, 3, length)
   call realloc(CSDefinitions%CS(i)%perim_sub, 3, length)
   call realloc(CSDefinitions%CS(i)%flowArea, length)
   call realloc(CSDefinitions%CS(i)%wetPerimeter, length)
   call realloc(CSDefinitions%CS(i)%totalArea, length)
   call realloc(CSDefinitions%CS(i)%area_min, length)
   call realloc(CSDefinitions%CS(i)%width_min, length)   
   
   CSDefinitions%CS(i)%id          = id 
   CSDefinitions%CS(i)%crossType   = cs_tabulated
   CSDefinitions%CS(i)%levelsCount = length
   CSDefinitions%CS(i)%plains      = plains
   do j = 1, numlevels
      CSDefinitions%CS(i)%height(j)  = level(j)
      CSDefinitions%CS(i)%flowWidth(j)   = flowWidth(j)
      CSDefinitions%CS(i)%totalWidth(j) = max(totalWidth(j), flowWidth(j))
   enddo
   if (closed) then
      CSDefinitions%CS(i)%height(numLevels+1) = CSDefinitions%CS(i)%height(numLevels)+1d-5
      CSDefinitions%CS(i)%flowWidth(numLevels+1) = 1d-5
      CSDefinitions%CS(i)%totalWidth(numLevels+1) = 1d-5
   endif
   CSDefinitions%CS(i)%closed = closed
   
   call createTablesForTabulatedProfile(CSDefinitions%CS(i))
!
   ! Initialize the summer dike information of the newly added cross-section
   if (flowArea > 1.0d-5 .or. totalArea > 1.0d-5) then
      allocate(CSDefinitions%CS(i)%summerdike)
      CSDefinitions%CS(i)%summerdike%crestLevel = crestLevel
      CSDefinitions%CS(i)%summerdike%baseLevel  = baseLevel
      CSDefinitions%CS(i)%summerdike%flowArea   = flowArea
      CSDefinitions%CS(i)%summerdike%totalArea  = totalArea
   endif
   
   ! Initialize groundlayer information of the newly added cross-section
   allocate(CSDefinitions%CS(i)%groundlayer)
   if (groundlayerUsed) anyGroundLayer       = .true.
   CSDefinitions%CS(i)%groundlayer%used      = groundlayerUsed
   CSDefinitions%CS(i)%groundlayer%thickness = groundlayer

   AddTabCrossSectionDefinition = i
end function AddTabCrossSectionDefinition


!> Adds YZ cross-sectin definition to the set of cross-section (complete with all attributes).
integer function AddYZCrossSectionDefinition(CSDefinitions, id, Count, y, z, frictionCount, frictionSectionID,           &
                                             frictionSectionFrom, frictionSectionTo, levelsCount, storageLevels, storage)
   use m_alloc

   type(t_CSDefinitionSet) CSDefinitions
   character(len=*), intent(in) :: id
   integer, intent(in)                        :: Count         !< Size of arrays: y and z
   integer, intent(in)                        :: frictionCount !< Number of friction settings
   integer, intent(in)                        :: levelsCount   !< Number of levels
   double precision, intent(in)               :: y(count)      !< transversal y-co-ordinate
   double precision, intent(in)               :: z(count)      !< level z-co-ordinate
   character(IdLen), intent(in)               :: frictionSectionID(frictionCount)
   double precision, intent(in)               :: frictionSectionFrom(frictionCount)
   double precision, intent(in)               :: frictionSectionTo(frictionCount)
   double precision, intent(in)               :: storageLevels(levelsCount)      !< Storage levels
   double precision, intent(in)               :: storage(levelsCount)            !< Storage values

   integer                                    :: i

   CSDefinitions%count = CSDefinitions%count+1
   if (CSDefinitions%count > CSDefinitions%size) then
      call realloc(CSDefinitions)
   endif

   i = CSDefinitions%count

   CSDefinitions%CS(i)%id                    = id 
   CSDefinitions%CS(i)%levelsCount           = Count
   CSDefinitions%CS(i)%crossType             = cs_YZ_Prof
   CSDefinitions%CS(i)%frictionSectionsCount = frictionCount
   CSDefinitions%CS(i)%storLevelsCount       = levelsCount

   call realloc(CSDefinitions%CS(i)%y, Count)
   call realloc(CSDefinitions%CS(i)%z, Count)
   call realloc(CSDefinitions%CS(i)%storLevels, max(levelsCount,2))
   call realloc(CSDefinitions%CS(i)%YZstorage, max(levelsCount,2))
   call realloc(CSDefinitions%CS(i)%frictionSectionID, frictionCount)
   call realloc(CSDefinitions%CS(i)%frictionSectionFrom, frictionCount)
   call realloc(CSDefinitions%CS(i)%frictionSectionTo, frictionCount)

   CSDefinitions%CS(i)%y = y
   CSDefinitions%CS(i)%z = z

   CSDefinitions%CS(i)%storLevels = storageLevels
   if (levelsCount==1)then
      CSDefinitions%CS(i)%storLevels(2) = storageLevels(1)+1.0
      CSDefinitions%CS(i)%YZstorage(1) = storage(1)
      CSDefinitions%CS(i)%YZstorage(2) = storage(1)
   else
      CSDefinitions%CS(i)%YZstorage = storage
   endif

   CSDefinitions%CS(i)%frictionSectionID   = frictionSectionID
   CSDefinitions%CS(i)%frictionSectionFrom = frictionSectionFrom
   CSDefinitions%CS(i)%frictionSectionTo   = frictionSectionTo
   
   allocate(CSDefinitions%CS(i)%groundlayer)
   CSDefinitions%CS(i)%groundlayer%used      = .false.
   CSDefinitions%CS(i)%groundlayer%thickness = 0d0

   AddYZCrossSectionDefinition = i
end function AddYZCrossSectionDefinition

!> Free the memory used by cross-section
subroutine deallocCrossSection(cross)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CrossSection) :: cross !< Current cross-section

   ! Local variables

   ! Program code

   if (.not. cross%IsCopy) then
   
      if (associated(cross%convTab)) then
         call dealloc(cross%convTab)
      endif
      
      if (associated(cross%tabDef)) then
         call dealloc(cross%tabDef)
      endif
      
      if (allocated(cross%frictionSectionID))   deallocate(cross%frictionSectionID)      !< Friction Section Identification
      if (allocated(cross%frictionSectionFrom)) deallocate(cross%frictionSectionFrom)    !<
      if (allocated(cross%frictionSectionTo))   deallocate(cross%frictionSectionTo)      !<
      if (allocated(cross%frictionTypePos))     deallocate(cross%frictionTypePos)        !< Friction type for positive flow direction
      if (allocated(cross%frictionValuePos))    deallocate(cross%frictionValuePos)       !< Friction value for positive flow direction
      if (allocated(cross%frictionTypeNeg))     deallocate(cross%frictionTypeNeg)        !< Friction type for negative flow direction
      if (allocated(cross%frictionValueNeg))    deallocate(cross%frictionValueNeg)       !< Friction value for negative flow direction

   endif

   cross%convTab => null()
   cross%tabDef => null()

end subroutine deallocCrossSection

!> Free the memory used by the cross-sections
subroutine deallocCrossSections(crs)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CrossSectionSet) :: crs !< Current cross-section set

   ! Local variables
   integer i

   ! Program code
   if (associated(crs%cross)) then

      do i = 1, crs%count
         if(.not. crs%cross(i)%IsCopy) then
            call dealloc(crs%cross(i))
         endif
         crs%cross(i)%tabDef  => null()
         crs%cross(i)%convTab => null()
      enddo

      if (associated(crs%crossSectionIndex)) then
         deallocate(crs%crossSectionIndex)
      endif
      
      deallocate(crs%cross)
      crs%crossSectionIndex => null()
      crs%cross => null()
      crs%size  = 0
      crs%count = 0
   endif
   
end subroutine deallocCrossSections

!> Allocate memory for cross-sections
subroutine reallocCrossSections(crs)
   ! Modules

   implicit none
   ! Input/output parameters
   type(t_CrossSectionSet), intent(inout)      :: crs !< Current cross-section set

   ! Local variables
   type(t_CrossSection), pointer, dimension(:) :: oldCrs
   integer :: ierr

   ! Program code

   if (crs%size > 0) then
      oldCrs=>crs%cross
   endif

   if (crs%growsBy <=0) then
      crs%growsBy = 200
   endif
   allocate(crs%cross(crs%size+crs%growsBy), stat=ierr)
   if (ierr /=0) then
      call setMessage(LEVEL_FATAL, 'allocation error')
   endif
   
   if (crs%size > 0) then
      crs%cross(1:crs%size) = oldCrs(1:crs%size)
      deallocate(oldCrs)
   endif
   crs%size = crs%size+crs%growsBy
end subroutine

!> Add a cross-section on an reach
integer function  AddCrossSectionByVariables(crs, CSDef, branchid, location, iref, bedLevel, bedFrictionType, &
                  bedFriction, groundFrictionType, groundFriction)

   type(t_CrossSectionSet)          :: crs       !< Cross-section set
   type(t_CSDefinitionSet)          :: CSDef     !< Cross-section definition set
   integer, intent(in)              :: branchid  !< Reach id (integer)
   double precision, intent(in)     :: location  !< Location on the reach
   integer, intent(in)              :: iref      !< Current cross-section
   double precision, intent(in)     :: bedLevel
   integer, intent(in)              :: bedFrictionType !< bed friction type
   double precision, intent(in)     :: bedFriction     !< Bed friction value
   integer, intent(in)              :: groundFrictionType !< Groundlayer fricton type
   double precision, intent(in)     :: groundFriction     !< Groundlayer friction value

   integer                          :: i

   crs%count = crs%count+1
   !! THREEDI-278
   !! NOTE: below is a dangerous realloc, see THREEDI-278
   !! By reallocating here, all pointers that used to point to input crs set (e.g. pcross in t_culvert),
   !! will afterwards point to undefined memory, since it was reallocated elsewhere.
   if (crs%count > crs%size) then
      call realloc(crs)
   endif

   i = crs%count

   if (branchid/=0) then
      crs%cross(i)%branchid = branchid
   endif

   crs%cross(i)%location            = location
   crs%cross(i)%bedLevel            = bedLevel
   crs%cross(i)%shift               = bedlevel
   
   crs%cross(i)%bedFrictionType     = bedFrictionType
   crs%cross(i)%bedFriction         = bedFriction
   crs%cross(i)%groundFrictionType  = groundFrictionType
   crs%cross(i)%groundFriction      = groundFriction
   crs%cross(i)%itabDef             = iref
   crs%cross(i)%tabDef              => CSDef%CS(iref)

   call SetParsCross(CSDef%CS(iref), crs%cross(i))
   
   AddCrossSectionByVariables = i
   
end function  AddCrossSectionByVariables
                  
subroutine SetParsCross(CrossDef, cross)
   type(t_CStype), intent(in)          :: CrossDef
   type(t_CrossSection), intent(inout) :: Cross
      
   double precision                 :: surfLevel
   integer :: j
   type(t_CrossSection)             :: crossSection
   double precision                 :: cz
   double precision                 :: conv

   cross%crossType = Crossdef%crossType
   cross%closed    = Crossdef%closed


   if (Crossdef%crossType == cs_tabulated) then
      cross%surfaceLevel = Crossdef%height(Crossdef%levelsCount)
      cross%bedLevel     = Crossdef%height(1) + cross%shift
   else
      surfLevel = -9999.0 !< JanM: TODO Parametriseer deze waarde, levelmiss toevoegen aan derived type CSDef
      do j = 1, Crossdef%levelsCount
         if (surfLevel < Crossdef%z(j)) then
            surfLevel = Crossdef%z(j)
         endif
      enddo
      cross%surfaceLevel = surflevel
      cross%bedLevel     = cross%shift
   endif
   ! Get Characteristic Height and Width
   call SetCharHeightWidth(cross)
   
   if (cross%tabDef%groundlayer%used .and. (cross%tabDef%groundlayer%Area < 1.0d-8) ) then

      crossSection = CopyCross(cross)
      
      if (associated(crossSection%tabDef%summerdike)) then
         deallocate(crossSection%tabDef%summerdike)
         crossSection%tabDef%summerdike => null()
      endif
      
      crossSection%tabDef%groundlayer%used      = .false.
      crossSection%tabDef%groundlayer%thickness = 0.0d0
      
      call GetCSParsFlowCross(crossSection,  cross%tabDef%groundlayer%thickness, 0.0d0, cz,      &
                              cross%tabDef%groundlayer%area, cross%tabDef%groundlayer%perimeter, &
                              cross%tabDef%groundlayer%width, conv)
      
      ! Deallocate Temporary Copy
      call dealloc(crossSection)
      
   endif                        
end subroutine SetParsCross
   

subroutine setGroundLayerData(crossDef, thickness)

   type(t_CStype), pointer, intent(inout) :: crossDef
   double precision                       :: thickness
      
   double precision                    :: area 
   double precision                    :: perimeter
   double precision                    :: width
   double precision                    :: af_sub(3), perim_sub(3)
   
   if (Thickness <= 0.0d0) then
      crossDef%groundlayer%used      = .false.
      crossDef%groundlayer%thickness = 0.0d0
      crossDef%groundlayer%area      = 0.0d0
      crossDef%groundlayer%perimeter = 0.0d0
      crossDef%groundlayer%width     = 0.0d0
      return
   endif
   
   select case(crossDef%crossType)
      case (CS_TABULATED)
         call GetTabSizesFromTables(thickness, crossDef, .true., area, width, perimeter, af_sub, perim_sub, CS_TYPE_NORMAL)
      case (CS_CIRCLE)
         call CircleProfile(thickness, crossDef%diameter, area, width, perimeter, CS_TYPE_NORMAL)
      case (CS_EGG)
         call EggProfile(thickness, crossDef%diameter, area, width, perimeter)
      case default
         call SetMessage(LEVEL_ERROR, 'INTERNAL ERROR: Unknown type of cross section')
   end select
      
   crossDef%groundlayer%used      = .true.
   crossDef%groundlayer%thickness = thickness
   crossDef%groundlayer%area      = area
   crossDef%groundlayer%perimeter = perimeter
   crossDef%groundlayer%width     = width

end subroutine setGroundLayerData
   

!> Add a cross-section on a reach, using a cross section defined on another reach
integer function AddCrossSectionByCross(crs, cross, branchid, location)

   type(t_CrossSectionSet)          :: crs       !< Cross-section set
   type(t_CrossSection), intent(in) :: cross     !< Cross-section
   integer, intent(in)              :: branchid  !< Reach id
   double precision, intent(in)     :: location  !< Location on the reach
   integer                          :: i

   crs%count = crs%count+1
   
   if (crs%count > crs%size) then
      call realloc(crs)
   endif

   i = crs%count

   crs%cross(i) = cross
   crs%cross(i)%branchid = branchid
   crs%cross(i)%location = location

   AddCrossSectionByCross = i

end function AddCrossSectionByCross

!> Get actual number of cross-sections in cross-section set
integer function GetCrossSectionsCount(crs)
  type(t_CrossSectionSet)                  :: crs

  GetCrossSectionsCount = crs%count
end function GetCrossSectionsCount

!> interpolation of width arrays.
subroutine interpolateWidths(height1, width1, levelsCount1, height2, width2, levelsCount2, height, width, levelsCount, f)
   ! modules

   implicit none

   ! variables
   integer, intent(in)                    :: levelsCount1               !< number of levels in height1, width1
   integer, intent(in)                    :: levelsCount2               !< number of levels in height2, width2
   integer, intent(in)                    :: levelsCount                !< number of levels in height, width
   double precision, intent(in)           :: height1(levelsCount1)      !< levels (z) of width array at location 1
   double precision, intent(in)           :: width1(levelsCount1)       !< widths as function of height at location 1
   double precision, intent(in)           :: height2(levelsCount2)      !< levels (z) of width array at location 2
   double precision, intent(in)           :: width2(levelsCount2)       !< widths as function of height at location 2
   double precision, intent(in)           :: height(levelsCount)        !< levels (z) of interpolated width array at location
   double precision, intent(out)          :: width(levelsCount)         !< interpolated widths as function of height
   double precision, intent(in)           :: f

   ! local variables
   integer              :: i
   integer              :: i1
   integer              :: i2
   double precision     :: frac
   double precision     :: w1
   double precision     :: w2
   !program code
   i1 = 1
   i2 = 1
   do i = 1, levelscount
      ! Fill first part of width array, and interpolate missing values
      if (height1(i1) > height(i) .and. i1 == 1 ) then
         w1 = width1(i1)
      elseif (height1(i1) == height(i) ) then
         w1 = width1(i1)
         i1 = min(levelsCount1,i1 + 1)
      elseif(height1(i1) < height(i)) then
         w1 = width1(i1)
      else
         frac = (height(i) - height1(i1 - 1)) / (height1(i1) - height1(i1 - 1))
         w1 = (1.0d0 - frac) * width1(i1 - 1) + frac * width1(i1)
      endif
      if (height2(i2) > height(i) .and. i2==1 ) then
         w2 = width2(i2)
      elseif (height2(i2) == height(i) ) then
         w2 = width2(i2)
         i2 = min(levelsCount2,i2 + 1)
      elseif(height2(i2) < height(i)) then
         w2 = width2(i2)
      else
         frac = (height(i) - height2(i2 - 1)) / (height2(i2) - height2(i2 - 1))
         w2 = (1.0d0 - frac) * width2(i2 - 1) + frac  *width2(i2)
      endif

      width(i) = (1d0 - f) * w1 + f * w2

   enddo

end subroutine interpolateWidths

!> perform a powerlaw interpolation between  val1 and val2 and power p
!! using the formula: result = ((1-f)*(val1)**p+f*(val2)**p)**1/p
double precision function interpolateScalarP(val1, val2, f, p)
   ! modules

   implicit none
   ! variables
   double precision, intent(in)    :: val1
   double precision, intent(in)    :: val2
   double precision, intent(in)    :: f      !< weight
   double precision, intent(in)    :: p      !< power
   ! local variables

   !program code
   interpolateScalarP = ((1.0d0 - f) * (val1)**p + f * (val2)**p)**(1.0d0 / p)
end function interpolateScalarP

!> perform a linear interpolation between  val1 and val2
!! using the formula: result = (1-f)*val1+f*val2
double precision function interpolateScalar(val1, val2, f)
   ! modules

   implicit none
   ! variables
   double precision, intent(in)    :: val1
   double precision, intent(in)    :: val2
   double precision, intent(in)    :: f
   ! local variables

   !program code
   interpolateScalar = (1.0d0 - f)*val1 + f*val2
end function interpolateScalar

subroutine mergeLevels(heightin1, levelsCount1, heightin2, levelsCount2, height, levelsCount, height1, height2, SetbedLevelToZero)
   ! modules

   implicit none
   ! variables
   integer                                            :: levelsCount1
   integer                                            :: levelsCount2
   double precision, dimension(levelsCount1)          :: heightin1
   double precision, dimension(levelsCount2)          :: heightin2
   double precision, allocatable, dimension(:)        :: height
   integer                                            :: levelsCount
   double precision, allocatable, dimension(:)        :: height1
   double precision, allocatable, dimension(:)        :: height2
   logical, optional                                  :: SetbedLevelToZero

   ! local variables
   integer           :: i1
   integer           :: i2
   double precision  :: offset1
   double precision  :: offset2
   !program code
   allocate(height(levelscount1+levelscount2))
   allocate(height1(levelscount1+1))
   allocate(height2(levelscount2+1))

   ! Set bedLevel at 0
   offset1 = heightin1(1)
   offset2 = heightin2(1)
   if (present(SetbedLevelToZero)) then
      if (.not. SetbedLevelToZero) then
         offset1 = 0.0
         offset2 = 0.0
      endif
   endif

   do i1 = 1, levelscount1
      height(i1)  = heightin1(i1) - offset1
      height1(i1) = heightin1(i1) - offset1
   enddo
   do i2 = 1, levelscount2
      height(levelsCount1+i2) = heightin2(i2) - offset2
      height2(i2)            = heightin2(i2) - offset2
   enddo

   levelsCount = levelscount1+levelscount2

   call regulatehlv(height, levelscount)

   height1(levelsCount1+1) = max(height1(levelsCount1), height2(levelsCount2)+1)
   height2(levelsCount2+1) = max(height1(levelsCount1), height2(levelsCount2)+1)

end subroutine mergeLevels

!> this subroutine
subroutine useBranchOrdersCrs(crs, brs)
   ! modules

   implicit none
   ! variables
   type(t_CrossSectionSet)          :: crs       !< Current cross-section set
   type(t_branchSet)                :: brs       !< Set of reaches

   ! local variables
   integer  ibr, orderNumberCount
   integer  i
   integer  iorder
   integer  ics
   integer  crsCount
   integer  minindex
   integer  minOrdernumber
   integer  minBranchindex
   double precision  minoffset
   integer, allocatable, dimension(:,:)   :: orderNumber       !< first index contains orderNumber, second contains start position for this ordernumber
   type(t_CrossSection)                   :: cross

   !program code
   allocate(crs%crossSectionIndex(crs%count), orderNumber(crs%Count+2,2))

   ! order cross sections, first on order number, then on branch index, last on offset
   orderNumberCount = 1
   orderNumber(1,1) = -1
   orderNumber(1,2) = 1
   crsCount = crs%count
   crs%crossSectionIndex = -1
   do ics = 1, crsCount
      minindex = ics
      if (crs%cross(ics)%branchid <= 0 ) then
         crs%crossSectionIndex(ics) = ics
      else
         ibr = crs%cross(ics)%branchid
         minordernumber = getOrderNumber(brs, ibr)
         minBranchindex = ibr
         minoffset = crs%cross(ics)%location
         do i = ics, crsCount
            if (crs%cross(i)%branchid <= 0) then
               minindex = i
               crs%crossSectionIndex(i) = ics
               minOrderNumber = -1
               exit
            else
               ibr = crs%cross(i)%branchid
               if (minOrderNumber > getOrdernumber(brs, ibr)) then
                  minOrderNumber =  getOrdernumber(brs, ibr)
                  minBranchindex = ibr
                  minOffset = crs%cross(i)%location
                  minIndex = i
               elseif (minOrderNumber == getOrdernumber(brs, ibr))then
                  if (minBranchIndex > ibr) then
                     minBranchindex = ibr
                     minOffset = crs%cross(i)%location
                     minIndex = i
                  elseif (minBranchIndex == ibr) then
                     if (minoffset > crs%cross(i)%location) then
                        minOffset = crs%cross(i)%location
                        minIndex = i
                     endif
                  endif
               endif
            endif
         enddo
      endif
      cross = crs%cross(ics)
      crs%cross(ics) = crs%cross(minindex)
      crs%cross(minindex) = cross
      if (orderNumber(orderNumberCount,1) /= minOrderNumber) then
         orderNumberCount = orderNumberCount + 1
         orderNumber(orderNumberCount, 1) = minOrderNumber
         orderNumber(orderNumberCount, 2) = ics
      endif
   enddo
   orderNumber(orderNumberCount+1,1) = -999
   orderNumber(orderNumberCount+1,2) = crsCount+1
   ! Now check all cross sections on branches of the same order (-1 orders can be skipped)
   do iorder = 2, orderNumberCount
      ics = orderNumber(iorder, 2)
      do while (ics < orderNumber(iorder+1,2) )
         ! ics is now the first cross section on the branch
         ibr = crs%cross(ics)%branchid
         
         ! because crs%cross can be reallocated a copy of the cross section has to be made
         cross = crs%cross(ics)
         cross%IsCopy = .true.
         call findNeighbourAndAddCrossSection(brs, crs, ibr, cross, cross%location, .true., orderNumber, orderNumberCount)
         ! find last cross section on branch
         do while (ics < crs%count .and. crs%cross(ics+1)%branchid == ibr)
            ics = ics+1
         enddo
         ! The last cross section in CRS is the last cross section on a branch
         if (ics < crs%count) then
            if (crs%cross(ics+1)%branchid == ibr) then
               ics = ics +1
            endif
         endif
         
         ! because crs%cross can be reallocated a copy of the cross section has to be made
         cross = crs%cross(ics)
         cross%IsCopy = .true.
         call findNeighbourAndAddCrossSection(brs, crs, ibr, cross, cross%location, .false., orderNumber, orderNumberCount)
         ics = ics +1
      enddo
   enddo

   deallocate(orderNumber)

end subroutine useBranchOrdersCrs

integer function getOrderNumber(brs, ibr)
   type(t_branchSet), intent(in)    :: brs       !< Set of reaches
   integer, intent(in)              :: ibr

   if (ibr <= brs%count) then
      getOrderNumber = brs%branch(ibr)%orderNumber
   else
      getOrderNumber = -1
   endif
end function getOrderNumber

!> The first cross section and the last cross section on a branch (B1), with an order number > 0
!! can be used to interpolate over a connection node. The algorithm is to locate a neighbouring
!! branch (same begin or end node) with the same order number (B2). When found the cross section
!! is added to the branch B2, with an offset outside the branch (either -OFFSET, or LENGTH+OFFSET).
recursive subroutine findNeighbourAndAddCrossSection(brs, crs, branchid, cross, offset, beginNode, orderNumber, orderNumberCount)
   ! modules

   implicit none

   ! variables
   integer  :: orderNumberCount
   type(t_CrossSectionSet)          :: crs       !< Current cross-section set
   type(t_branchSet)                :: brs       !< Set of reaches
   integer                          :: branchid    !< branch for which a neighbour is requested
   type(t_CrossSection)             :: cross     !< cross section
   double precision                 :: offset    !< location of cross section on branch
   logical                          :: beginNode !< indicates whether the begin or end node is to be used of the branch
   integer, dimension(:,:)          :: orderNumber       !< first index contains orderNumber, second contains start position for this ordernumber

   ! local variables
   integer                          :: nodeIndex
   integer                          :: ibr
   integer                          :: i
   integer                          :: iorder
   integer                          :: j
   integer                          :: idum
   logical                          :: found
   logical                          :: nextBeginNode
   double precision                 :: distanceFromNode
   double precision                 :: location
   type(t_branch), pointer          :: branch

   !program code
   if (brs%count==0) then
      return
   endif

   branch => brs%branch(branchid)

   ibr = branch%index
   if (beginNode) then
      nodeIndex = branch%FromNode%index
      distanceFromNode = offset
   else
      nodeIndex = branch%ToNode%index
      distanceFromNode = branch%length - offset
   endif
   found = .false.
   do i = 1, brs%count
      if ( (brs%branch(i)%orderNumber == branch%orderNumber) .and. (i /= ibr ) ) then
         if (brs%branch(i)%ToNode%index == nodeIndex) then
            found = .true.
            ! interpolation is required over begin node of reach i
            location = brs%branch(i)%length + distanceFromNode
            nextBeginNode = .true.              ! NOTE: since the end node of the reach is connected to a node of the first reach, now check the begin node
            exit
         elseif (brs%branch(i)%FromNode%index == nodeIndex) then
            ! interpolation is required over end node of reach i
            found = .true.
            location = - distanceFromNode
            nextBeginNode = .false.             ! NOTE: since the begin node of the reach is connected to a node of the first reach, now check the end node
            exit
         endif
      endif
   enddo

   if (found) then
      idum = AddCrossSection(crs, cross, i, location)
      ! look in crs if a cross section is placed on brs%branch(i)
      found = .false.
      do iorder = 1, orderNumberCount
         if (orderNumber(iorder,1) == brs%branch(i)%orderNumber) then
            exit
         endif
      enddo

      do j = orderNumber(iorder,2), orderNumber(iorder+1,2)-1
         if (brs%branch(crs%cross(j)%branchid)%index ==i) then
            found = .true.
            exit
         endif
      enddo
      if (.not. found) then
         ! Now check for another neighbouring branch
         call findNeighbourAndAddCrossSection(brs, crs, i, cross, location, nextBeginNode, orderNumber, orderNumberCount)
      endif
   endif

end subroutine findNeighbourAndAddCrossSection

subroutine GetCSParsFlowInterpolate(cross1, cross2, f, dpt, u1, cz, flowArea, wetPerimeter, flowWidth, conv, af_sub, perim_sub, cz_sub)

   use m_GlobalParameters
   
   implicit none

   type (t_CrossSection), intent(in)       :: cross1         !< cross section
   type (t_CrossSection), intent(in)       :: cross2         !< cross section
   double precision, intent(in)            :: f              !< cross = (1-f)*cross1 + f*cross2
   double precision, intent(in)            :: dpt            !< water depth at cross section
   double precision, intent(in)            :: u1             !< flow velocity
   double precision, intent(inout)         :: cz             !< roughness at cross section
   double precision, intent(out)           :: flowArea       !< flow area for given DPT
   double precision, intent(out)           :: wetPerimeter   !< wet perimeter for given DPT
   double precision, intent(out)           :: flowWidth      !< flow width of water surface
   double precision, intent(out)           :: conv           !< conveyance
   double precision, intent(out), optional :: af_sub(3)      
   double precision, intent(out), optional :: perim_sub(3)      
   double precision, intent(out), optional :: cz_sub(3)      

   double precision                      :: cz1 = 0.d0
   double precision                      :: cz2 = 0.d0
   double precision                      :: flowArea1
   double precision                      :: flowArea2
   double precision                      :: wetPerimeter1
   double precision                      :: wetPerimeter2
   double precision                      :: flowWidth1   
   double precision                      :: flowWidth2   
   double precision                      :: conv1        
   double precision                      :: conv2        
   type (t_CrossSection), save           :: crossi         !< intermediate virtual crosssection     

   double precision                      :: sdArea  = 0.0d0
   double precision                      :: sdWidth = 0.0d0
   double precision                      :: charHeight

   double precision                      :: af_sub_local1(3), af_sub_local2(3)      
   double precision                      :: perim_sub_local1(3), perim_sub_local2(3)
   double precision                      :: cz_sub_local1(3), cz_sub_local2(3)      

   if(cross1%crossIndx == cross2%crossIndx) then
      ! Same Cross-Section, no interpolation needed 
      call GetCSParsFlowCross(cross1, dpt, u1, cz, flowArea, wetPerimeter, flowWidth, conv, af_sub_local1, &
                     perim_sub_local1, cz_sub_local1)
      if (present(af_sub)   ) af_sub    = af_sub_local1
      if (present(perim_sub)) perim_sub = perim_sub_local1
      if (present(cz_sub)   ) cz_sub    = cz_sub_local1
   else
      select case (cross1%crosstype)
         
         case (CS_CIRCLE, CS_EGG)
            ! Create an interpolate cross-section here and call GetCSParsFlowCross for intermediate cross-section
            if(.not.associated(crossi%tabDef)) then 
               allocate(crossi%tabDef)
               if(.not.associated(crossi%tabDef%groundlayer)) then 
                  allocate(crossi%tabDef%groundlayer)
               endif 
            endif 
         
            crossi%crosstype          = cross1%crosstype   
            crossi%bedLevel           = (1.0d0 - f) * cross1%bedLevel + f * cross2%bedLevel
            crossi%surfaceLevel       = (1.0d0 - f) * cross1%surfaceLevel + f * cross2%surfaceLevel
            ! for circular cross sections there is only one roughness value set in bedfriction
            crossi%bedFrictionType    = cross1%frictionTypePos(1)
            crossi%bedFriction        = (1.0d0 - f) * cross1%frictionValuePos(1) + f * cross2%frictionValuePos(1)
            crossi%groundFriction     = (1.0d0 - f) * cross1%groundFriction + f * cross2%groundFriction
            crossi%groundFrictionType = cross1%groundFrictionType
            crossi%tabDef%diameter    = (1.0d0 - f) * cross1%tabDef%diameter + f * cross2%tabDef%diameter
         
            call GetCSParsFlowCross(crossi, dpt, u1, cz, flowArea, wetPerimeter, flowWidth, conv, af_sub_local1, &
                     perim_sub_local1, cz_sub_local1)
            if (present(af_sub)   ) af_sub    = af_sub_local1
            if (present(perim_sub)) perim_sub = perim_sub_local1
            if (present(cz_sub)   ) cz_sub    = cz_sub_local1
            
         case default                             ! Call GetCSParsFlowCross twice and interpolate the results 

            cz1 = cz
            cz2 = cz
            call GetCSParsFlowCross(cross1, dpt, u1, cz1, flowArea1, wetPerimeter1, flowWidth1, conv1, af_sub_local1, &
                     perim_sub_local1, cz_sub_local1, doSummerDike = .false.)
            call GetCSParsFlowCross(cross2, dpt, u1, cz2, flowArea2, wetPerimeter2, flowWidth2, conv2, af_sub_local2, &
                     perim_sub_local2, cz_sub_local2, doSummerDike = .false.)
            
            ! Summer Dikes
            call interpolateSummerDike(cross1, cross2, f, dpt, sdArea, sdWidth, .true.)
         
            flowArea     = (1.0d0 - f) * flowArea1     + f * flowArea2  + sdArea
            flowWidth    = (1.0d0 - f) * flowWidth1    + f * flowWidth2 + sdWidth
            wetPerimeter = (1.0d0 - f) * wetPerimeter1 + f * wetPerimeter2

            cz           = (1.0d0 - f) * cz1   + f * cz2
            conv         = (1.0d0 - f) * conv1 + f * conv2
            if (present(af_sub)) then
               af_sub = (1.0d0 - f) * af_sub_local1     + f * af_sub_local2
            endif
            if (present(perim_sub)) then
               perim_sub = (1.0d0 - f) * perim_sub_local1     + f * perim_sub_local2
            endif
            if (present(cz_sub)) then
               cz_sub = (1.0d0 - f) * cz_sub_local1     + f * cz_sub_local2
            endif
            if (cross1%crosstype .ne. CS_YZ_PROF) then
            
               ! Find out the conveyance value for interpolated non yz cross section
               conv = cz * flowArea * sqrt(flowArea / wetPerimeter)  ! extra_friction_depth
            
               ! Criteria to satisfy the criteria  in normup i.e cz(m)*cz(m)*wet
               if (cz *cz * flowArea < 1.0d0) then
                  conv = sqrt(flowArea *flowArea / wetPerimeter)
               endif
            
               if (cross1%closed) then
                  charHeight = (1.0d0 - f) * cross1%charHeight + f * cross2%charHeight
                  if (dpt >= charHeight) flowWidth = 0.0d0
               endif
            
            endif
         
      end select
      
   endif

end subroutine GetCSParsFlowInterpolate

subroutine interpolateSummerDike(cross1, cross2, f, dpt, sdArea, sdWidth, doFlow)

   use m_GlobalParameters
   
   implicit none

   type (t_CrossSection), intent(in)     :: cross1         !< cross section
   type (t_CrossSection), intent(in)     :: cross2         !< cross section
   double precision, intent(in)          :: f              !< cross = (1-f)*cross1 + f*cross2
   double precision, intent(in)          :: dpt            !< water depth at cross section
   double precision, intent(out)         :: sdArea         !< summer dike area
   double precision, intent(out)         :: sdWidth        !< summer dike width
   logical, intent(in)                   :: doFlow         !< flag, true = flow, false = total

   type(t_summerdike), pointer           :: summerdike
   type(t_summerdike), pointer           :: summerdike1
   type(t_summerdike), pointer           :: summerdike2
   double precision                      :: shift1
   double precision                      :: shift2
   double precision                      :: wlev

   if (.not. anySummerDike) then
      sdArea  = 0.0d0
      sdWidth = 0.0d0
      return
   endif
   
   allocate(summerdike)
   summerdike1 => cross1%tabDef%summerdike
   summerdike2 => cross2%tabDef%summerdike
            
   shift1 = cross1%shift
   shift2 = cross2%shift
            
   if (associated(summerdike1) .and. associated(summerdike2)) then
      summerdike%crestLevel = (1.0d0 - f) * (summerdike1%crestLevel + shift1) + f * (summerdike2%crestLevel + shift2)
      summerdike%baseLevel  = (1.0d0 - f) * (summerdike1%baseLevel + shift1)  + f * (summerdike2%baseLevel + shift2)
      summerdike%flowArea   = (1.0d0 - f) * summerdike1%flowArea   + f * summerdike2%flowArea
      summerdike%totalArea  = (1.0d0 - f) * summerdike1%totalArea   + f * summerdike2%totalArea
      if (f <= 0.5d0) then
         summerdike%hysteresis = summerdike1%hysteresis
      else
         summerdike%hysteresis = summerdike2%hysteresis
      endif
   elseif (associated(summerdike1) .and. .not. associated(summerdike2)) then
      summerdike%crestLevel = (1.0d0 - f) * (summerdike1%crestLevel + shift1) + f * cross2%bedlevel
      summerdike%baseLevel  = (1.0d0 - f) * (summerdike1%baseLevel + shift1)  + f * cross2%bedlevel
      summerdike%flowArea   = (1.0d0 - f) * summerdike1%flowArea
      summerdike%totalArea  = (1.0d0 - f) * summerdike1%totalArea
      summerdike%hysteresis = summerdike1%hysteresis
   elseif (.not. associated(summerdike1) .and. associated(summerdike2)) then
      summerdike%crestLevel = (1.0d0 - f) * cross1%bedlevel + f * (summerdike2%crestLevel + shift2)
      summerdike%baseLevel  = (1.0d0 - f) * cross1%bedlevel + f * (summerdike2%baseLevel + shift2)
      summerdike%flowArea   = f * summerdike2%flowArea
      summerdike%totalArea  = f * summerdike2%totalArea
      summerdike%hysteresis = summerdike2%hysteresis
   else
      ! Not any Summer Dike Data
      deallocate(summerdike)
      summerdike => null()
   endif

   if (associated(summerdike)) then
            
      wlev = getBobInterpolate(cross1, cross2, f) + dpt
         
      if (doFlow) then
      
         ! Get Summer Dike Flow Data
         call GetSummerDikeFlow(summerdike, wlev, sdArea, sdWidth)
      
      else
      
         ! Get Summer Dike Total Data
         call GetSummerDikeTotal(summerdike, wlev, sdArea, sdWidth)
      
         if (associated(summerdike1) .and. associated(summerdike2)) then
            if (f <= 0.5d0) then
               summerdike1%hysteresis = summerdike%hysteresis
            else
               summerdike2%hysteresis = summerdike%hysteresis
            endif
         elseif (associated(summerdike1) .and. .not. associated(summerdike2)) then
            summerdike1%hysteresis = summerdike%hysteresis
         elseif (.not. associated(summerdike1) .and. associated(summerdike2)) then
            summerdike2%hysteresis = summerdike%hysteresis
         endif
   
      endif
         
      deallocate(summerdike)
      summerdike => null()
            
   else
      sdArea  = 0.0d0
      sdWidth = 0.0d0
   endif

end subroutine interpolateSummerDike

subroutine GetCSParsFlowCross(cross, dpt, u1, cz, flowArea, wetPerimeter, flowWidth, conv, af_sub, perim_sub, cz_sub, doSummerDike)      
   use m_GlobalParameters
   use precision_basics
   use m_Roughness

   type (t_CrossSection), intent(in) :: cross          !< cross section definition
   double precision, intent(in)      :: dpt            !< water depth at cross section
   double precision, intent(in)      :: u1             !< flow velocity
   double precision, intent(inout)   :: cz             !< roughness at cross section
   double precision, intent(out)     :: flowArea       !< flow area for given DPT
   double precision, intent(out)     :: wetPerimeter   !< wet perimeter for given DPT
   double precision, intent(out)     :: flowWidth      !< flow width of water surface
   double precision, intent(out)     :: conv           !< conveyance
   logical, intent(in), optional     :: doSummerDike   !< Switch to calculate Summer Dikes or not
   double precision, intent(out), optional :: af_sub(3)      
   double precision, intent(out), optional :: perim_sub(3)      
   double precision, intent(out), optional :: cz_sub(3)      

   type(t_CSType), pointer           :: crossDef
   double precision                  :: widgr
   double precision                  :: czg
   logical                           :: getSummerDikes
   double precision                  :: af_sub_local(3)      
   double precision                  :: perim_sub_local(3)      
   double precision                  :: cz_sub_local(3)   
   integer                           :: numsect
   integer                           :: i
  
   perim_sub_local = 0d0
   if (dpt <= 0.0d0) then
      flowArea     = 0.0
      wetPerimeter = 0.0
      flowWidth    = sl
      conv         = 0.0
      return
   endif

   crossDef => cross%tabDef
   conv = 0.0d0

   select case(cross%crosstype)
      case (CS_TABULATED)
         if (present(doSummerDike))then
            getSummerdikes = doSummerdike
         else
            getSummerDikes = .true.
         endif
         call TabulatedProfile(dpt, cross, .true., getSummerDikes, flowArea, flowWidth, wetPerimeter, af_sub_local, perim_sub_local, CS_TYPE_NORMAL)
      case (CS_CIRCLE)
         call CircleProfile(dpt, crossDef%diameter, flowArea, flowWidth, wetPerimeter, CS_TYPE_NORMAL)
      case (CS_EGG)
         call EggProfile(dpt, crossDef%diameter, flowArea, flowWidth, wetPerimeter)
      case (CS_YZ_PROF)
         call YZProfile(dpt, u1, cz, cross%convtab, 0, flowArea, flowWidth, wetPerimeter, conv)
      case default
         call SetMessage(LEVEL_ERROR, 'INTERNAL ERROR: Unknown type of cross section')
   end select

   if (cross%crosstype /= CS_TABULATED ) then
      af_sub_local = 0d0
      af_sub_local(1) = flowArea
      
      perim_sub_local = 0d0
      perim_sub_local(1) = wetPerimeter
      
      cz_sub_local = 0d0
      cz_sub_local(1) = cz
   endif
   ! correction for groundlayers
   if (crossDef%groundlayer%used) then
      widGr = crossDef%groundlayer%width
      flowArea = flowArea - crossDef%groundlayer%area
      af_sub_local(1) = af_sub_local(1) - crossDef%groundlayer%area

      wetPerimeter = wetPerimeter - crossDef%groundlayer%perimeter + widGr
      perim_sub_local(1) = perim_sub_local(1) - crossDef%groundlayer%perimeter + widGr
   endif

   if (conv == 0.0d0 .and. wetperimeter > 0.0d0) then
      if (comparereal(cz, 0.0d0) == 0) then
         if (cross%frictionSectionsCount==0) then
            cz = GetChezy(cross%bedfrictiontype, cross%bedfriction, flowArea  / wetPerimeter, dpt, u1)
            cz_sub_local = cz 
         else
            numsect = 0
            do i = 1, 3
               if (perim_sub_local(i) > eps .and. af_sub_local(i) > 0.0d0) then
                  cz_sub_local(i) = GetChezy(cross%frictionTypePos(i), cross%frictionValuePos(i), &
                                             af_sub_local(i)  / perim_sub_local(i), dpt, u1)
                  numsect = numsect+1
               else
                  cz_sub_local(i) = 0
               endif
            enddo
         endif
         
         if (crossDef%groundlayer%used) then
            czg = GetChezy(cross%groundfrictiontype, cross%groundfriction, flowArea/wetPerimeter, dpt, u1)
            cz = ((widGr * czg * czg + (wetPerimeter - widgr) * cz* cz) / (wetperimeter))**0.5d0
            cz_sub_local(1) = ((widGr * czg * czg + max(0d0, (perim_sub_local(1) - widgr)) * cz* cz) / (perim_sub_local(1)))**0.5d0
         endif

         do i = 1, 3
            if (perim_sub_local(i) > eps) then
               conv = conv + cz_sub_local(i) * af_sub_local(i) * sqrt(af_sub_local(i) / perim_sub_local(i))  ! extra_friction_depth
               
            endif
         enddo
         cz = conv/(flowArea*dsqrt(flowArea/wetPerimeter))
      else
         cz_sub_local = cz
      endif
         conv = cz * flowArea * sqrt(flowArea / wetPerimeter)  ! extra_friction_depth
   

      
      !        criteria to satisfy the criteria  in normup i.e cz(m)*cz(m)*wet
      if (cz * cz * flowArea < 1.0d0) then
         conv = sqrt(flowArea * flowArea / wetPerimeter)
      endif
      
   endif
   
   if (present(af_sub)) then
      af_sub = af_sub_local
   endif
   if (present(perim_sub)) then
      perim_sub = perim_sub_local
   endif
   if (present(cz_sub)) then
      cz_sub = cz_sub_local
   endif
   
end subroutine GetCSParsFlowCross

subroutine GetCSParsTotalInterpolate(cross1, cross2, f, dpt, totalArea, totalWidth, calculationOption)

   use m_GlobalParameters
   
   implicit none

   type (t_CrossSection), intent(in)     :: cross1         !< cross section
   type (t_CrossSection), intent(in)     :: cross2         !< cross section
   double precision, intent(in)          :: f              !< cross = (1-f)*cross1 + f*cross2
   double precision, intent(in)          :: dpt            !< water depth at cross section
   double precision, intent(out)         :: totalArea      !< total area for given DPT
   double precision, intent(out)         :: totalWidth     !< total width of water surface
                                                           !> type of total area computation, possible values:\n
                                                           !! CS_TYPE_PREISMAN  Ordinary total area computation, with possible Preisman lock on top\n
                                                           !! CS_TYPE_PLUS      Total area for only the expanding part of the cross section (Nested Newton method)\n
                                                           !! CS_TYPE_MIN       Total area for only the narrowing part of the cross section (Nested Newton method)
   integer, intent(in)                   :: calculationOption 

   double precision                      :: totalArea1
   double precision                      :: totalArea2
   double precision                      :: totalWidth1   
   double precision                      :: totalWidth2   
   type (t_CrossSection), save           :: crossi         !< intermediate virtual crosssection     

   double precision                      :: sdArea  = 0.0d0
   double precision                      :: sdWidth = 0.0d0


   if(cross1%crossIndx == cross2%crossIndx) then
      ! Same Cross-Section, no interpolation needed 
      call GetCSParsTotalCross(cross1, dpt, totalArea, totalWidth, calculationOption)
   else
      select case (cross1%crosstype)
      case (CS_CIRCLE, CS_EGG)
         ! Create an interpolate cross-section here and call GetCSParstotalCross for intermediate cross-section
         if(.not.associated(crossi%tabDef)) then 
            allocate(crossi%tabDef)
            if(.not.associated(crossi%tabDef%groundlayer)) then 
               allocate(crossi%tabDef%groundlayer)
            endif 
         endif 
         
         crossi%crosstype          = cross1%crosstype   
         crossi%bedLevel           = (1.0d0 - f) * cross1%bedLevel + f * cross2%bedLevel
         crossi%surfaceLevel       = (1.0d0 - f) * cross1%surfaceLevel + f * cross2%surfaceLevel
         crossi%bedFrictionType    = cross1%bedFrictionType
         crossi%bedFriction        = (1.0d0 - f) * cross1%bedFriction + f * cross2%bedFriction
         crossi%groundFriction     = (1.0d0 - f) * cross1%groundFriction + f * cross2%groundFriction
         crossi%groundFrictionType = cross1%groundFrictionType
         crossi%tabDef%diameter    = (1.0d0 - f) * cross1%tabDef%diameter + f * cross2%tabDef%diameter
         
         call GetCSParsTotalCross(crossi, dpt, totalArea, totalWidth, calculationOption)
         
      case default                             ! Call GetCSParstotalCross twice and interpolate the results 

         call GetCSParsTotalCross(cross1, dpt, totalArea1, totalWidth1, calculationOption, doSummerDike = .false.)
         call GetCSParsTotalCross(cross2, dpt, totalArea2, totalWidth2, calculationOption, doSummerDike = .false.)
         
         ! Summer Dikes
         call interpolateSummerDike(cross1, cross2, f, dpt, sdArea, sdWidth, .false.)
         
         totalArea     = (1.0d0 - f) * totalArea1     + f * totalArea2  + sdArea
         totalWidth    = (1.0d0 - f) * totalWidth1    + f * totalWidth2 + sdWidth
        
      end select
      
   endif

end subroutine GetCSParsTotalInterpolate

subroutine GetCSParsTotalCross(cross, dpt, totalArea, totalWidth, calculationOption, doSummerDike)

   use m_GlobalParameters
   ! Global Variables
   type (t_CrossSection), intent(in) :: cross           !< cross section
   double precision, intent(in)      :: dpt             !< water depth at cross section
   double precision, intent(out)     :: totalArea       !< total area for given DPT
   double precision, intent(out)     :: totalWidth      !< total width of water surface
                                                        !> type of total area computation, possible values:\n
                                                        !! CS_TYPE_PREISMAN  Ordinary total area computation, with possible Preisman lock on top\n
                                                        !! CS_TYPE_PLUS      Total area for only the expanding part of the cross section (Nested Newton method)\n
                                                        !! CS_TYPE_MIN       Total area for only the narrowing part of the cross section (Nested Newton method)
   integer, intent(in)               :: calculationOption 
   logical, intent(in), optional     :: doSummerDike    !< Switch to calculate Summer Dikes or not


   ! Local Variables
   type(t_CSType), pointer           :: crossDef
   double precision                  :: wetperimeter
   double precision                  :: conv
   double precision                  :: cz = 60.0d0
   double precision                  :: u1 = 0.0d0
   double precision                  :: wlev            !< water level at cross section
   logical                           :: getSummerDikes
   double precision                  :: af_sub(3), perim_sub(3)
   if (dpt <= 0.0d0) then
      totalArea = 0.0d0
      totalWidth = sl
      return
   endif

   crossDef => cross%tabdef

   select case(cross%crosstype)
      case (CS_TABULATED, CS_RECTANGLE)
         wlev = cross%bedlevel + dpt
         if (present(doSummerDike))then
            getSummerdikes = doSummerdike
         else
            getSummerDikes = .true.
         endif
         call TabulatedProfile(dpt, cross, .false., getSummerDikes, totalArea, totalWidth, wetPerimeter, af_sub, perim_sub, calculationOption)
      case (CS_CIRCLE)
         call CircleProfile(dpt, crossDef%diameter, totalArea, totalWidth, wetPerimeter, calculationOption)
      case (CS_EGG)
         !TODO:
         call EggProfile(dpt, crossDef%diameter, totalArea, totalWidth, wetPerimeter)
      case (CS_YZ_PROF)
         call YZProfile(dpt, u1, cz, cross%convtab, 1, totalArea, totalWidth, wetPerimeter, conv)
      case default
         call SetMessage(LEVEL_ERROR, 'INTERNAL ERROR: Unknown type of cross section')
      end select

end subroutine GetCSParsTotalCross

subroutine TabulatedProfile(dpt, cross, doFlow, getSummerDikes, area, width, perimeter, af_sub, perim_sub, calculationOption)

   use m_GlobalParameters

   implicit none

   double precision, intent(in)      :: dpt
   type (t_CrossSection), intent(in) :: cross           
   logical, intent(in)               :: doFlow  !> True: Flow, otherwise Total
   logical, intent(in)               :: getSummerDikes                           
   double precision, intent(out)     :: width
   double precision, intent(out)     :: area
   double precision, intent(out)     :: perimeter
   double precision, intent(out)     :: af_sub(3)
   double precision, intent(out)     :: perim_sub(3)
   integer, intent(in)               :: calculationOption 

   ! local parameters
   type(t_CSType), pointer           :: crossDef
   type(t_summerdike), pointer       :: summerdike
  
   double precision  :: wlev
   double precision  :: sdArea
   double precision  :: sdWidth

   crossDef => cross%tabDef
   
   call GetTabSizesFromTables(dpt, crossDef, doFlow, area, width, perimeter, af_sub, perim_sub, calculationOption)
      
   if (associated(crossDef%summerdike) .and. getSummerDikes .and. calculationOption /= CS_TYPE_MIN) then
   
      wlev = cross%bedlevel + dpt
      
      allocate(summerdike)
      summerdike%crestLevel = crossDef%summerdike%crestLevel + cross%shift
      summerdike%baseLevel  = crossDef%summerdike%baseLevel + cross%shift
      summerdike%flowArea   = crossDef%summerdike%flowArea
      summerdike%totalArea  = crossDef%summerdike%totalArea
      summerdike%hysteresis = crossDef%summerdike%hysteresis
      
      ! Summerdike Adjusting
      if (doFlow) then
      
         ! Get Summer Dike Flow Data
         call GetSummerDikeFlow(summerdike, wlev, sdArea, sdWidth)
      
      else
      
         ! Get Summer Dike Total Data
         call GetSummerDikeTotal(summerdike, wlev, sdArea, sdWidth)
         crossDef%summerdike%hysteresis = summerdike%hysteresis
      
      endif
      
      deallocate(summerdike)
      summerdike => null()
      
      area  = area  + sdArea
      width = width + sdWidth
      
   endif

end subroutine TabulatedProfile

subroutine GetTabSizesFromTables(dpt, pCSD, doFlow, area, width, perimeter, af_sub, perim_sub, calculationOption)

   use m_GlobalParameters

   implicit none

   double precision, intent(in)                 :: dpt
   type (t_CSType), pointer, intent(in)         :: pCSD           
   logical, intent(in)                          :: doFlow  !> True: Flow, otherwise Total
   double precision, intent(out)                :: width
   double precision, intent(out)                :: area
   double precision, intent(out)                :: perimeter
   double precision, intent(out)                :: af_sub(3)
   double precision, intent(out)                :: perim_sub(3)
   integer, intent(in)                          :: calculationOption 

   ! local parameters
   integer                                      :: levelsCount
   double precision, dimension(:), pointer      :: heights
   double precision, dimension(:), pointer      :: widths
   
   ! Pre-Calculated Table for Tabulated/River Profiles
   double precision, dimension(:,:), pointer    :: af_sub_tab        !< Flow Areas for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:,:), pointer    :: width_sub_tab     !< Width for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:,:), pointer    :: perim_sub_tab     !< Wetted Perimeter for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:), pointer      :: flowArea_tab      !< Flow Areas
   double precision, dimension(:), pointer      :: wetPerimeter_tab  !< Wet Perimeters
   double precision, dimension(:), pointer      :: totalArea_tab     !< Total Areas
   double precision, dimension(:), pointer      :: area_min_tab      !< Area for a narrowing part of a cross section (Nested Newton)
   double precision, dimension(:), pointer      :: width_min_tab     !< Width for a narrowing part of a cross section (Nested Newton)
   
  
   double precision  :: d1
   double precision  :: d2
   double precision  :: dTop
   double precision  :: eTop
   integer           :: ilev, isec
   
   double precision  :: factor
   double precision  :: area_plus
   double precision  :: a
   double precision  :: b
   
   levelsCount = pCSD%levelsCount
   heights => pCSD%height   
   dTop = heights(levelsCount) - heights(1)

   af_sub    = 0.0d0
   perim_sub = 0.0d0
   width     = 0.0d0
   area      = 0.0d0
   perimeter = 0.0d0

   if (doFlow) then
      
      af_sub_tab       => pCSD%af_sub
      width_sub_tab    => pCSD%width_sub
      perim_sub_tab    => pCSD%perim_sub
      flowArea_tab     => pCSD%flowArea
      wetPerimeter_tab => pCSD%wetPerimeter

      if (dpt > dTop) then
      
         do isec = 1,3
            af_sub(isec) = af_sub_tab(isec, levelsCount) + (dpt - dTop) * width_sub_tab(isec, levelsCount)
            perim_sub(isec) =  perim_sub_tab(isec, levelsCount)
         enddo
         
         area      = flowArea_tab(levelsCount)
         perimeter = wetPerimeter_tab(levelsCount)
         width     = pCSD%flowWidth(levelsCount)
         eTop      = pCSD%flowWidth(levelsCount)

         if ((eTop > ThresholdForPreismannLock)) then  !hk: add only when this part is meant to carry water

            if (af_sub(3) > 0d0) then
               isec = 3
            elseif (af_sub(2) > 0d0) then
               isec = 2
            else 
               isec = 1
            endif
            perim_sub(isec) = perim_sub_tab(isec, levelsCount) + 2.0d0 * (dpt - dTop)

            area = area + (dpt - dTop) * eTop
            perimeter = perimeter + 2.0d0 * (dpt - dTop)
            
         endif
         
      else

         do ilev = 2, levelsCount
         
            d1 = heights(ilev - 1) - heights(1)
            d2 = heights(ilev) - heights(1)
            
            if (dpt == d2) then

               do isec = 1,3
                  af_sub(isec)    = af_sub_tab(isec, ilev)
                  perim_sub(isec) = perim_sub_tab(isec, ilev)
               enddo
               
               area = flowArea_tab(ilev)
               perimeter = wetPerimeter_tab(ilev)
               width = pCSD%flowWidth(ilev)
            
            elseif (dpt > d1 .and. dpt <= d2) then
               
               if ((d2 - d1) > 0.0d0) then
                  
                  factor = (dpt - d1) / (d2 - d1)
                  
                  do isec = 1,3
                     
                     if (pCSD%plainsLocation(isec) == 0) cycle
                     
                     width = width_sub_tab(isec, ilev - 1) + (width_sub_tab(isec, ilev) - width_sub_tab(isec, ilev - 1)) * factor
                     
                     area_plus = (width_sub_tab(isec, ilev - 1) + width) * (dpt - d1) * 0.5d0
                     af_sub(isec)    = af_sub_tab(isec, ilev - 1) + area_plus
                     
                     if (width_sub_tab(isec, ilev) > width_sub_tab(isec, ilev - 1)) then
                        a = width - width_sub_tab(isec, ilev - 1)
                        b = dpt - d1
                        perim_sub(isec) = perim_sub_tab(isec, ilev - 1) + 2.0d0 * dsqrt(0.25d0 * a * a + b * b)
                     else
                        perim_sub(isec) = perim_sub_tab(isec, ilev - 1)
                     endif
                     
                  enddo
                  
                  width = pCSD%flowWidth(ilev - 1) + (pCSD%flowWidth(ilev) - pCSD%flowWidth(ilev - 1)) * factor
                  
                  area_plus = (pCSD%flowWidth(ilev - 1) + width) * (dpt - d1) * 0.5d0
                  area      = flowArea_tab(ilev - 1) + area_plus
                  
                  a = width - pCSD%flowWidth(ilev - 1)
                  b = dpt - d1
                  perimeter = wetPerimeter_tab(ilev - 1) + 2.0d0 * dsqrt(0.25d0 * a * a + b * b)
                  
               else
                  
                  do isec = 1,3
                     af_sub(isec)    = af_sub_tab(isec, ilev - 1)
                     perim_sub(isec) = perim_sub_tab(isec, ilev - 1)
                  enddo

                  area = flowArea_tab(ilev - 1)
                  perimeter = wetPerimeter_tab(ilev - 1)
                  width = pCSD%flowWidth(ilev - 1)
                  
               endif
            
            endif
            
         enddo
         
      endif
   
 
   else

      ! Calculation for Storage
      widths        => pCSD%totalWidth
      totalArea_tab => pCSD%totalArea
      area_min_tab  => pCSD%area_min
      width_min_tab => pCSD%width_min
   
      if (dpt > dTop) then
         
         area  = totalArea_tab(levelsCount)
         width = pCSD%totalWidth(levelsCount)
         eTop  = pCSD%totalWidth(levelsCount)

         if ((eTop > ThresholdForPreismannLock)) then  !hk: add only when this part is meant to carry water
            area = area + (dpt - dTop) * eTop
         endif
         
      else

         do ilev = 2, levelsCount
         
            d1 = heights(ilev - 1) - heights(1)
            d2 = heights(ilev) - heights(1)
            
            if (dpt == d2) then

               if (.not. calculationOption == CS_TYPE_MIN) then               
                  area = totalArea_tab(ilev)
                  width = pCSD%totalWidth(ilev)
               else
                  area =area_min_tab(ilev)
                  width = width_min_tab(ilev)
               endif
            
            elseif (dpt > d1) then
               
               d1 = heights(ilev - 1) - heights(1)
               d2 = heights(ilev) - heights(1)
               
               if ((d2 - d1) > 0.0d0) then
                  
                  factor = (dpt - d1) / (d2 - d1)
                  
                  width = pCSD%totalWidth(ilev - 1) + (pCSD%totalWidth(ilev) - pCSD%totalWidth(ilev - 1)) * factor
                  area_plus = (pCSD%totalWidth(ilev - 1) + width) * (dpt - d1) * 0.5d0
                   
                  if (.not. calculationOption == CS_TYPE_MIN) then
                     area = totalArea_tab(ilev - 1) + area_plus
                  else
                     
                     if (width < pCSD%totalWidth(ilev - 1)) then                       
                        area = area_min_tab(ilev - 1) + area_plus
                     else
                        area = area_min_tab(ilev - 1)
                     endif
                     
                     width = width_min_tab(ilev - 1) + pCSD%totalWidth(ilev - 1) - width
                     
                  endif
                  
               else
                  
                  if (.not. calculationOption == CS_TYPE_MIN) then               
                     area = totalArea_tab(ilev - 1)
                     width = pCSD%totalWidth(ilev - 1)
                  else
                     area = area_min_tab(ilev - 1)
                     width = width_min_tab(ilev - 1)
                  endif
                  
               endif
            
            endif
            
         enddo
         
      endif

   endif
   
end subroutine GetTabSizesFromTables

subroutine GetTabFlowSectionFromTables(dpt, pCross, isector, area, width, perimeter)

   use m_GlobalParameters

   implicit none

   double precision, intent(in)                 :: dpt
   type (t_CrossSection), pointer, intent(in)   :: pCross         
   integer, intent(in)                          :: isector
   double precision, intent(out)                :: width
   double precision, intent(out)                :: area
   double precision, intent(out)                :: perimeter

   ! local parameters
   integer                                      :: levelsCount
   double precision, dimension(:), pointer      :: heights
   integer, dimension(:), pointer               :: plainsLocation                              
   type(t_CSType), pointer                      :: pCSD
   
   ! Pre-Calculated Table for Tabulated/River Profiles
   double precision, dimension(:,:), pointer    :: af_sub_tab        !< Flow Areas for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:,:), pointer    :: width_sub_tab     !< Width for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:,:), pointer    :: perim_sub_tab     !< Wetted Perimeter for Sub-Sections (Main, FP1 and FP2)
   double precision, dimension(:), pointer      :: flowArea_tab      !< Flow Areas
   double precision, dimension(:), pointer      :: wetPerimeter_tab  !< Wet Perimeters
   
  
   double precision  :: d1
   double precision  :: d2
   double precision  :: dTop
   double precision  :: eTop
   integer           :: ilev
   
   double precision  :: factor
   double precision  :: area_plus
   double precision  :: a
   double precision  :: b
   
   pCSD => pCross%tabDef
   
   levelsCount    =  pCSD%levelsCount
   heights        => pCSD%height
   plainsLocation => pCSD%plainsLocation
   dTop           =  heights(levelsCount) - heights(1)

   width     = 0.0d0
   area      = 0.0d0
   perimeter = 0.0d0

   af_sub_tab       => pCSD%af_sub
   width_sub_tab    => pCSD%width_sub
   perim_sub_tab    => pCSD%perim_sub
   flowArea_tab     => pCSD%flowArea
   wetPerimeter_tab => pCSD%wetPerimeter

   if (dpt > dTop) then
      
      width     =  width_sub_tab(isector, levelsCount)
         
      eTop = pCSD%flowWidth(levelsCount)

      if ((eTop > ThresholdForPreismannLock)) then  !hk: add only when this part is meant to carry water

         area = af_sub_tab(isector, levelsCount) + (dpt - dTop) * width_sub_tab(isector, levelsCount)

         if (levelsCount <= plainsLocation(isector)) then
            perimeter = perim_sub_tab(isector, levelsCount) + 2.0d0 * (dpt - dTop)
         else
            perimeter = perim_sub_tab(isector, levelsCount)
         endif

      else
         
         area      = af_sub_tab(isector, levelsCount)
         perimeter =  perim_sub_tab(isector, levelsCount)
         
      endif
         
   else

      do ilev = 2, levelsCount
         
         d1 = heights(ilev - 1) - heights(1)
         d2 = heights(ilev) - heights(1)
            
         if (dpt == d2) then

            area      = af_sub_tab(isector, ilev)
            perimeter = perim_sub_tab(isector, ilev)
            width     = width_sub_tab(isector, ilev)
               
         elseif (dpt > d1 .and. dpt <= d2) then
               
            if ((d2 - d1) > 0.0d0) then
                  
               if (pCSD%plainsLocation(isector) > 0)then
                     
                  factor = (dpt - d1) / (d2 - d1)
               
                  width = width_sub_tab(isector, ilev - 1) + (width_sub_tab(isector, ilev) - width_sub_tab(isector, ilev - 1)) * factor
                     
                  area_plus = (width_sub_tab(isector, ilev - 1) + width) * (dpt - d1) * 0.5d0
                  area      = af_sub_tab(isector, ilev - 1) + area_plus
                  
                  if (ilev <= plainsLocation(isector)) then
                     a = width - width_sub_tab(isector, ilev - 1)
                     b = dpt - d1
                     perimeter = perim_sub_tab(isector, ilev - 1) + 2.0d0 * dsqrt(0.25d0 * a * a + b * b)
                  else
                     perimeter = perim_sub_tab(isector, ilev - 1)
                  endif
                  
               else
                  area      = 0.0d0
                  perimeter = 0.0d0
                  width     = 0.0d0
               endif
                  
            else
                  
               area      = af_sub_tab(isector, ilev - 1)
               perimeter = perim_sub_tab(isector, ilev - 1)
               width     = width_sub_tab(isector, ilev - 1)

            endif
            
         endif
            
      enddo
         
   endif
   
end subroutine GetTabFlowSectionFromTables

subroutine GetSummerDikeFlow(summerdike, wlev, sdArea, sdWidth)

   implicit none

   type(t_summerdike), pointer                  :: summerdike
   double precision, intent(in)                 :: wlev
   double precision, intent(out)                :: sdArea
   double precision, intent(out)                :: sdWidth

   ! Local Parameters
   double precision                             :: sdtr
   double precision                             :: aflw
   double precision                             :: dflw
   double precision                             :: htop
   double precision                             :: hbas
   
   if (.not. associated(summerdike)) then
      sdWidth = 0.0d0
      sdArea  = 0.0d0
      return
   endif
   
   sdtr = summerDikeTransitionHeight
   
   aflw = summerdike%flowArea
   
   if (aflw > 0.0d0) then
   
      dflw = aflw / sdtr
      
      htop = summerdike%crestLevel
      hbas = summerdike%baseLevel
      
      if (wlev > (htop + sdtr)) then
         sdWidth = 0.0
         sdArea  = aflw
      elseif (wlev > (htop + 0.5d0 * sdtr)) then
         sdWidth = (htop + sdtr - wlev) * 4.0d0 * dflw / sdtr
         sdArea  = aflw - 0.5d0 * (htop + sdtr - wlev) * sdWidth
      elseif (wlev > htop) then
         sdWidth = (wlev - htop) * 4.0d0 * dflw / sdtr
         sdArea = 0.5d0 * (wlev - htop) * sdWidth
      else
         sdWidth = 0.0d0
         sdArea  = 0.0d0
      endif
   else
      sdWidth = 0.0d0
      sdArea  = 0.0d0
   endif
          
end subroutine GetSummerDikeFlow

subroutine GetSummerDikeTotal(summerdike, wlev, sdArea, sdWidth)

   implicit none

   type(t_summerdike), pointer                  :: summerdike
   double precision, intent(in)                 :: wlev
   double precision, intent(out)                :: sdArea
   double precision, intent(out)                :: sdWidth

   ! Local Parameters
   logical                                      :: hysteresis
   double precision                             :: sdtr
   double precision                             :: atot
   double precision                             :: htop
   double precision                             :: hbas
   double precision                             :: dtot

   if (.not. associated(summerdike)) then
      sdWidth = 0.0d0
      sdArea  = 0.0d0
      return
   endif
   
   sdtr = summerDikeTransitionHeight
   
   atot       = summerdike%totalArea
   hysteresis = summerdike%hysteresis
   
   if (atot > 0.0d0) then
   
      htop = summerdike%crestLevel
      hbas = summerdike%baseLevel
      
      if (wlev >= (htop + sdtr)) then
         hysteresis = .false.
      elseif (wlev <= hbas) then
         hysteresis = .true.
      else
      endif
      summerdike%hysteresis = hysteresis
      
      if (hysteresis) then
      
         dtot = atot / sdtr
         
         if (wlev > (htop + sdtr)) then
            sdWidth = 0.0d0
            sdArea  = atot
         elseif (wlev > (htop + 0.5d0 * sdtr)) then
            sdWidth = (htop + sdtr - wlev) * 4.0d0 * dtot / sdtr 
            sdArea  = atot - 0.5d0 * (htop + sdtr - wlev) * sdWidth
         elseif (wlev>htop) then
            sdWidth = (wlev - htop) * 4.0d0 * dtot / sdtr
            sdArea  = 0.5d0 * (wlev - htop) * sdWidth
         else
            sdWidth = 0.0d0
            sdArea  = 0.0d0
         endif
         
      else
      
         dtot = atot / (htop + sdtr - hbas)
         
         if (wlev > (htop + sdtr)) then
            sdWidth = 0.0d0
            sdArea  = atot
         elseif (wlev > (hbas + 0.5d0 * (htop + sdtr - hbas))) then
            sdWidth = (htop + sdtr - wlev) * 4.0d0 * dtot / (htop + sdtr - hbas)
            sdArea  = atot - 0.5d0*(htop + sdtr - wlev) * sdWidth
         elseif (wlev > hbas) then
            sdWidth = (wlev - hbas) * 4.0d0 * dtot / (htop + sdtr - hbas)
            sdArea  = 0.5d0 * (wlev - hbas) * sdWidth
         else
            sdWidth = 0.0d0
            sdArea  = 0.0d0
         endif
         
      endif
      
   else
      sdWidth = 0.0d0
      sdArea  = 0.0d0
   endif
   
end subroutine GetSummerDikeTotal

subroutine trapez(dpt, d1, d2, w1, w2, area, width, perimeter)
   implicit none

   double precision, intent(in)     :: dpt
   double precision, intent(in)     :: d1
   double precision, intent(in)     :: d2
   double precision, intent(in)     :: w1
   double precision, intent(in)     :: w2
   double precision, intent(out)    :: area
   double precision, intent(out)    :: width
   double precision, intent(out)    :: perimeter
   !
   width= w1 + (dpt - d1) * (w2 - w1) / (d2 - d1)
   area = (dpt - d1) * 0.5d0 * (width + w1)
   perimeter = 2.0d0 * dsqrt(0.25d0 * (width - w1)**2 + (dpt - d1)**2)
   
end subroutine trapez

subroutine CircleProfile(dpt, diameter, area, width, perimeter, calculationOption)
   use m_GlobalParameters

   implicit none

   double precision, intent(in)        :: dpt
   double precision, intent(in)        :: diameter
   double precision, intent(out)       :: width
   double precision, intent(out)       :: area
   double precision, intent(out)       :: perimeter
   integer, intent(in)                 :: calculationOption

!
! Local variables
!
   double precision               :: dacos
   double precision               :: dsqrt
   double precision               :: fi
   double precision               :: ra
   double precision               :: sq
   double precision               :: areacircle
   double precision               :: widthcircle
   
!
!! executable statements -------------------------------------------------------
!
   !
   !
   !
   ra = 0.5*diameter
   fi = dacos((ra-dpt)/ra)
   sq = dsqrt(dpt*(diameter - dpt))

   ! normal circle profile
   if (dpt<diameter) then
      areacircle      = dabs(fi*ra*ra - sq*(ra-dpt))
      perimeter       = dabs(2d0*fi*ra)
      widthcircle     = 2d0*sq
   else
      areacircle      = pi*ra*ra
      perimeter       = 2d0*pi*ra
      widthcircle     = 0d0
   endif
   select case(calculationOption)
   case(CS_TYPE_NORMAL)
      area = areacircle
      width = widthcircle
   case(CS_TYPE_PREISMAN)
      area = areacircle    + sl*dpt
      width = widthcircle  + sl
   case(CS_TYPE_PLUS, CS_TYPE_MIN)
      if (dpt<ra) then
         ! half circle profile
         area      = dabs(fi*ra*ra - sq*(ra-dpt)) + sl*dpt
         perimeter = 2d0*fi*ra
         width     = 2d0*sq + sl
      else
         area      = 0.5d0* pi*ra*ra + diameter*(dpt-ra)
         perimeter = 2d0*pi*ra
         width     = diameter
      endif   
   end select

   if (calculationOption == CS_TYPE_MIN) then
      area  = area  - areacircle
      width = width - widthcircle
   endif
   
end subroutine CircleProfile

subroutine EggProfile(dpt, diameter, area, width, perimeter)
   use m_GlobalParameters
   use precision_basics

   implicit none

   double precision, intent(in)        :: dpt
   double precision, intent(in)        :: diameter
   double precision, intent(out)       :: width
   double precision, intent(out)       :: area
   double precision, intent(out)       :: perimeter

   double precision                    :: r
   double precision                    :: e

   r = 0.5*diameter
   if ((dpt>0) .and. (dpt<=.2*r)) then
      perimeter = 2*r*(.5*atan((dsqrt(r*dpt - dpt*dpt)/(.5*r - dpt))))
      width = r*sin(atan(dsqrt(dpt*r - dpt*dpt)/(.5*r - dpt))) + sl
      e = .25*r*r*atan((dpt - .5*r)/(dsqrt(dpt*r - dpt*dpt))) + .392699082*r*r +       &
         & (dpt - .5*r)*dsqrt(dpt*r - dpt*dpt)
      area = e + sl*dpt

   elseif ((dpt>.2*r) .and. (dpt<=2.*r)) then
      perimeter = 2*r*(2.39415093538065 -                                          &
            & 3.*atan((2*r - dpt)/(dsqrt(5.*r*r + 4*r*dpt - dpt*dpt))))
      width = 2.*((dsqrt(5.*r*r + 4.*r*dpt - dpt*dpt)) - 2.*r) + sl
      e = 9.*r*r*atan((dpt - 2*r)/(dsqrt(9.*r*r - ((dpt - 2*r)**2)))) + (dpt - 2*r)&
         & *dsqrt(9*r*r - ((dpt - 2*r)**2)) - 4.*r*(dpt - .2*r)                   &
         & + 10.11150997913956*r*r
      area = .11182380480168*r*r + e + sl*dpt

   elseif ((dpt>2*r) .and. (comparereal(dpt,3*r, 1d-6) < 0)) then
      perimeter = 2.*r*(2.39415093538065 +                                         &
            & atan((dpt - 2*r)/(dsqrt( - 3.*r*r + 4*r*dpt - dpt*dpt))))
      width = 2*r*cos(atan((dpt - 2*r)/(dsqrt( - 3.*r*r + 4*dpt*r - dpt*dpt)))) + sl
      e = r*r*(atan((dpt - 2*r)/(dsqrt(r*r - ((dpt - 2*r)**2))))) + (dpt - 2*r)   &
            & *dsqrt(r*r - ((dpt - 2*r)**2))
      area = 3.02333378394124*r*r + e + sl*dpt

   else if (dpt>=3*r) then
      perimeter = 7.92989452435109*r
      width = sl
      area = 4.59413011073614*r*r + sl*dpt
   endif
end subroutine EggProfile

subroutine YZProfile(dpt, u1, cz, convtab, i012, area, width, perimeter, conv)
   use m_GlobalParameters

   implicit none

   double precision, intent(in)        :: dpt
   double precision, intent(in)        :: u1
   double precision, intent(inout)        :: cz
   integer, intent(in)                 :: i012
   type(t_crsu), intent(inout)         :: convtab
   double precision, intent(out)       :: width
   double precision, intent(out)       :: area
   double precision, intent(out)       :: perimeter
   double precision, intent(out)       :: conv

! locals
integer            :: nr, i1, i2, i, japos
double precision   :: a1, a2, c1, c2, z1, z2, hu1, hu2, hh1, hh2, dh1, dh2
double precision   :: c_a, c_b !< variables related to extrapolation
                              !! above defined profile

double precision   :: r3, ar1, ar2


if (i012 .eq. 0) then                                ! look at u points, mom. eq.

   nr = convtab%nru                                    ! number of table entries
   i  = convtab%iolu                                   ! last index found

   do while ( i + 1 < nr .and. convtab%hu(i+1) < dpt ) ! look up, imax = nr - 1
      i = i + 1
   enddo

   do while ( i     > 1  .and. convtab%hu(i)   > dpt ) ! look down, imin = 1
      i = i - 1
   enddo
   convtab%iolu = i                                    ! and store last index found

   i1  = i                                            ! so i1, i2 always inside table
   i2  = i+1
   hu2 = convtab%hu(i2) ; dh2 = hu2 - dpt
   if (dpt .LE. convtab%hu(i2) ) then !  .and. convtab%jopen .eq. 0) then     ! weightfactors. If profile closed no extrapolation
      hu1 = convtab%hu(i1) ; dh1 = dpt - hu1
      a1  = dh2 / ( hu2-hu1)                          ! eis parser: hu = wel monotoon stijgend
      a2 = 1d0 - a1
      !
      c1    = convtab%pf(i1) ; c2 = convtab%pf(i2)
      perimeter = a1*c1 + a2*c2
      !
      c1    = convtab%wf(i1) ; c2 = convtab%wf(i2)
      width = a1*c1 + a2*c2
      !
      z1    = convtab%af(i1) ; z2 = convtab%af(i2)
      ar1   = 0.5d0*dh1*(c1 + width)                     ! area above i1
      ar2   = 0.5d0*dh2*(c2 + width)                     ! area below i2
      area  = a1*(z1+ar1)   + a2*(z2-ar2)
      !
      japos = 1
      if (convtab%negcon .eq. 1) then
         if (u1 .lt. 0) japos = 0
      endif
      !
      if (japos .eq. 1) then
         z1 = convtab%cz1(i1)
         z2 = convtab%cz1(i2)   ! positive flow direction
         if (convtab%conveyType==CS_VERT_SEGM) then
         c1 = convtab%co1(i1)
         c2 = convtab%co1(i2)
         endif
      else
         z1 = convtab%cz2(i1)
         z2 = convtab%cz2(i2)   ! negative flow direction
         
         if (convtab%conveyType==CS_VERT_SEGM) then
            c1 = convtab%co2(i1)
            c2 = convtab%co2(i2)
         endif
      endif
      !
      if (convtab%conveyType==CS_LUMPED) then
         conv = (cz)*area*sqrt(area/perimeter)
      elseif (convtab%conveyType==CS_VERT_SEGM) then
         conv  = a1*c1 + a2*c2
      endif
      
      r3 = area / perimeter                            ! actual hydraulic radius
      convtab%chezy_act = conv / (area * dsqrt(r3))  ! Used in function ChezyFromConveyance
      cz = convtab%chezy_act
      !
   ELSE                                               ! when above profile
   ! positive direction/ negative direction? -> japos == 1 || japos != 1
   !
!(*)!  document SOBEK-21942: Change of roughness formulations in "Y-Z" and
   ! "Asymetrical Trapezium" profiles, Author:     Thieu van Mierlo
   !                                   Programmer: Daniel Abel
      WIDTH = convtab%wf (i2)
      perimeter = convtab%PF (i2)
      AREA  = convtab%AF (i2)
      CONV  = convtab%co1(i2)
      !
      if ( convtab%jopen .eq. 1) then ! for open profiles add extra conveyance
         AR2   = WIDTH*(DPT-HU2)
         AREA  = AREA + AR2
         perimeter = perimeter + 2*(DPT-HU2)
         r3    = AREA/perimeter                        ! actual hydraulic radius

         ! Determine Flow Direction
         if ( (convtab%negcon .eq. 1) .and. (u1 .lt. 0) ) then
            japos = 0
         else
            japos = 1
         endif
         if (convtab%conveyType==CS_VERT_SEGM) then
            if (japos .eq. 1) then
               c_b = convtab%b_pos_extr
               c_a = convtab%a_pos_extr
            else
               c_b = convtab%b_neg_extr
               c_a = convtab%a_neg_extr
            endif
            !
            conv = c_a*((dpt)**(c_b))
            ! Actual Chezy Value for ChezyFromConveyance
            convtab%chezy_act = conv / (AREA * DSQRT(r3))
            cz = convtab%chezy_act
         else
            conv = cz*area*sqrt(area/perimeter)
         endif
      endif

   endif

else                                                      ! look at left or right h, cont. eq.

   nr = convtab%nrhh(i012)                                  ! number of entries
   i  = convtab%iolh(i012)                                  ! last found

   do while ( i + 1 < nr .and. convtab%hh(i+1,i012) < dpt ) ! look up
      i = i + 1
   enddo

   do while ( i     > 1  .and. convtab%hh(i  ,i012) > dpt ) ! look down
      i = i - 1
   enddo
   convtab%iolh(i012) = i

   i1 = i                                                  ! so i1, i2 always inside table
   i2 = i+1
   hh2 = convtab%hh(i2,i012) ; dh2 = hh2 - dpt
   if (i2 .eq. nr .and. dpt .ge. convtab%hh(i2,i012) ) then ! Weightfactors. If profile closed no extrapolation
      a1 = 0d0 ; dh1 = 0
   else
      hh1 = convtab%hh(i1,i012) ; dh1 = dpt - hh1
      a1  = dh2 / ( hh2-hh1 )                              ! parser: hh = wel monotoon stijgend
   endif
   a2 = 1d0 - a1

   c1    = convtab%wt(i1,i012) ; c2 = convtab%wt(i2,i012)
   z1    = convtab%at(i1,i012) ; z2 = convtab%at(i2,i012)

   width = a1*c1 + a2*c2

   ar1   = 0.5d0*dh1*(c1 + width)                          ! area above i1
   ar2   = 0.5d0*dh2*(c2 + width)                          ! area below i2
   area  = a1*(z1+ar1)   + a2*(z2-ar2)

endif

end subroutine YZProfile

subroutine CalcConveyance(crs)

use M_newcross
   !use m_parseConveyance
   !use modelGlobalData
   !use convTables

   implicit none

   integer nc
   type(t_CrossSection)    :: crs
   type(t_crsu), pointer   :: convTab
   convtab=>null()
   nc = crs%tabDef%levelsCount
   call generateConvtab(convtab, crs%tabDef%levelsCount, crs%bedLevel, crs%tabDef%groundLayer%thickness, crs%tabDef%crossType, &
                        nc, crs%tabDef%frictionSectionsCount, crs%branchid, crs%bedFrictionType,                               &
                        crs%groundFriction, crs%tabdef%y, crs%tabdef%z,                                                        &
                        crs%frictionSectionFrom, crs%frictionSectionTo, crs%frictionTypePos,              &
                        crs%frictionValuePos, crs%frictionTypeNeg, crs%frictionValueNeg )
   crs%convTab => convTab

end subroutine CalcConveyance

double precision function getHighest1dLevelInterpolate(c1, c2, f)
   type (t_CrossSection), intent(in)     :: c1      !< cross section definition
   type (t_CrossSection), intent(in)     :: c2      !< cross section definition
   double precision, intent(in)          :: f       !< cross = (1-f)*cross1 + f*cross2
   
   double precision level1, level2
   
   level1 = getHighest1dLevel(c1)
   level2 = getHighest1dLevel(c2)
   
   getHighest1dLevelInterpolate = (1.0d0 - f) * level1 + f * level2
   
end function getHighest1dLevelInterpolate
   
double precision function getHighest1dLevelSingle(cross)

   type (t_CrossSection), intent(in)     :: cross      !< cross section
   
   integer levelsCount
   
   select case(cross%crosstype)
      case (CS_TABULATED)
         getHighest1dLevelSingle = cross%tabdef%height(cross%tabdef%levelscount) +cross%shift
      case (CS_CIRCLE)
         getHighest1dLevelSingle = cross%tabdef%diameter + cross%bedlevel
      case (CS_EGG)
         getHighest1dLevelSingle = 1.5d0 * cross%tabdef%diameter + cross%bedlevel
      case (CS_YZ_PROF)
         levelsCount = cross%convtab%nru
         getHighest1dLevelSingle = cross%convtab%hu(levelsCount) - cross%convtab%bob(1)
      case default
         call SetMessage(LEVEL_ERROR, 'INTERNAL ERROR: Unknown type of cross-section in getHighest1dLevelSingle')
   end select
 
end function getHighest1dLevelSingle

subroutine SetCharHeightWidth(cross)

type(t_CrossSection), intent(inout)       :: cross

   ! Code
   select case(cross%crosstype)
      case (CS_TABULATED)
            cross%charHeight = maxval(cross%tabDef%height) - minval(cross%tabDef%height)
            cross%charWidth  = maxval(cross%tabDef%flowWidth)
      case (CS_CIRCLE)
            cross%charHeight = cross%tabDef%diameter
            cross%charWidth  = cross%tabDef%diameter
      case (CS_EGG)
            cross%charHeight = 1.5d0 * cross%tabDef%diameter
            cross%charWidth  = cross%tabDef%diameter
      case (CS_YZ_PROF)
            cross%charHeight = maxval(cross%tabDef%z) - minval(cross%tabDef%z)
            cross%charWidth  = maxval(cross%tabDef%y) - minval(cross%tabDef%y)
      case default
         call SetMessage(LEVEL_ERROR, 'INTERNAL ERROR: Unknown type of cross-section in SetCharHeightWidth')
   end select
   
end subroutine SetCharHeightWidth

type(t_CrossSection) function CopyCross(CrossFrom)

   ! Returns a Copy from the given CrossSection

   type(t_CrossSection)      :: CrossFrom

   CopyCross%crossType          = CrossFrom%crossType
   CopyCross%branchid           = CrossFrom%branchid
   CopyCross%location           = CrossFrom%location
   CopyCross%bedLevel           = CrossFrom%bedLevel
   CopyCross%shift              = CrossFrom%shift
   CopyCross%surfaceLevel       = CrossFrom%surfaceLevel
   CopyCross%charHeight         = CrossFrom%charHeight
   CopyCross%charWidth          = CrossFrom%charWidth
   CopyCross%closed             = CrossFrom%closed
   CopyCross%bedFrictionType    = CrossFrom%bedFrictionType
   CopyCross%bedFriction        = CrossFrom%bedFriction
   CopyCross%groundFrictionType = CrossFrom%groundFrictionType
   CopyCross%groundFriction     = CrossFrom%groundFriction
   CopyCross%iTabDef            = CrossFrom%iTabDef
   CopyCross%csid               = CrossFrom%csid
   
   if (associated(CrossFrom%tabDef)) then
      allocate(CopyCross%tabDef)
      CopyCross%tabDef = CopyCrossDef(CrossFrom%tabDef)
   endif

   if (associated(CrossFrom%convTab)) then
      allocate(CopyCross%convTab)
      CopyCross%convTab = CopyCrossConv(CrossFrom%convTab)
   endif

end function CopyCross

type(t_CSType) function CopyCrossDef(CrossDefFrom)

   ! Returns a Copy from the given CrossSection Definition
   ! DEALLOCATE this Copy after Use!!!!!!!!!!!!!!!!!!!!!!!!!!

   type(t_CSType) :: CrossDefFrom
   
   CopyCrossDef%crossType   = CrossDefFrom%crossType
   CopyCrossDef%reference   = CrossDefFrom%reference
   CopyCrossDef%levelsCount = CrossDefFrom%levelsCount

   !*** data for tabulated cross sections ***
   if (allocated(CrossDefFrom%height)) then
      allocate(CopyCrossDef%height(CopyCrossDef%levelsCount))
      CopyCrossDef%height = CrossDefFrom%height
   endif
   if (allocated(CrossDefFrom%flowWidth)) then
      allocate(CopyCrossDef%flowWidth(CopyCrossDef%levelsCount))
      CopyCrossDef%flowWidth = CrossDefFrom%flowWidth
   endif
   if (allocated(CrossDefFrom%totalWidth)) then
      allocate(CopyCrossDef%totalWidth(CopyCrossDef%levelsCount))
      CopyCrossDef%totalWidth = CrossDefFrom%totalWidth
   endif
   
   if (allocated(CrossDefFrom%af_sub)) then
      allocate(CopyCrossDef%af_sub(3, CopyCrossDef%levelsCount))
      CopyCrossDef%af_sub = CrossDefFrom%af_sub
   endif
   if (allocated(CrossDefFrom%width_sub)) then
      allocate(CopyCrossDef%width_sub(3, CopyCrossDef%levelsCount))
      CopyCrossDef%width_sub = CrossDefFrom%width_sub
   endif
   if (allocated(CrossDefFrom%perim_sub)) then
      allocate(CopyCrossDef%perim_sub(3, CopyCrossDef%levelsCount))
      CopyCrossDef%perim_sub = CrossDefFrom%perim_sub
   endif
   if (allocated(CrossDefFrom%flowArea)) then
      allocate(CopyCrossDef%flowArea(CopyCrossDef%levelsCount))
      CopyCrossDef%flowArea = CrossDefFrom%flowArea
   endif
   if (allocated(CrossDefFrom%wetPerimeter)) then
      allocate(CopyCrossDef%wetPerimeter(CopyCrossDef%levelsCount))
      CopyCrossDef%wetPerimeter = CrossDefFrom%wetPerimeter
   endif
   if (allocated(CrossDefFrom%totalArea)) then
      allocate(CopyCrossDef%totalArea(CopyCrossDef%levelsCount))
      CopyCrossDef%totalArea = CrossDefFrom%totalArea
   endif
   if (allocated(CrossDefFrom%area_min)) then
      allocate(CopyCrossDef%area_min(CopyCrossDef%levelsCount))
      CopyCrossDef%area_min = CrossDefFrom%area_min
   endif
   if (allocated(CrossDefFrom%width_min)) then
      allocate(CopyCrossDef%width_min(CopyCrossDef%levelsCount))
      CopyCrossDef%width_min = CrossDefFrom%width_min
   endif
   
   CopyCrossDef%plains       = CrossDefFrom%plains
   
   !*** data for yz cross sections
   if (allocated(CrossDefFrom%y)) then
      allocate(CopyCrossDef%y(CopyCrossDef%levelsCount))
      CopyCrossDef%y = CrossDefFrom%y
   endif
   if (allocated(CrossDefFrom%z)) then
      allocate(CopyCrossDef%z(CopyCrossDef%levelsCount))
      CopyCrossDef%z = CrossDefFrom%z
   endif
   
   CopyCrossDef%storageType     = CrossDefFrom%storageType
   CopyCrossDef%storLevelsCount = CrossDefFrom%storLevelsCount

   if (allocated(CrossDefFrom%storLevels)) then
      allocate(CopyCrossDef%storLevels(max(CopyCrossDef%storLevelsCount, 2)))
      CopyCrossDef%storLevels = CrossDefFrom%storLevels
   endif
   if (allocated(CrossDefFrom%YZstorage)) then
      allocate(CopyCrossDef%YZstorage(CopyCrossDef%storLevelsCount))
      CopyCrossDef%YZstorage = CrossDefFrom%YZstorage
   endif

   CopyCrossDef%frictionSectionsCount = CrossDefFrom%frictionSectionsCount
   if (CopyCrossDef%frictionSectionsCount > 0) then
   
      allocate (CopyCrossDef%frictionSectionID(CopyCrossDef%frictionSectionsCount))
      allocate (CopyCrossDef%frictionSectionFrom(CopyCrossDef%frictionSectionsCount))
      allocate (CopyCrossDef%frictionSectionTo(CopyCrossDef%frictionSectionsCount))
   
      CopyCrossDef%frictionSectionID = CrossDefFrom%frictionSectionID
      CopyCrossDef%frictionSectionFrom = CrossDefFrom%frictionSectionFrom
      CopyCrossDef%frictionSectionTo   = CrossDefFrom%frictionSectionTo
     
   endif
   
   !** Circle and egg profile data
   CopyCrossDef%diameter = CrossDefFrom%diameter

   !--- information for Summerdikes
   if (associated(CrossDefFrom%summerdike)) then
      allocate(CopyCrossDef%summerdike)
      CopyCrossDef%summerdike%crestLevel = CrossDefFrom%summerdike%crestLevel
      CopyCrossDef%summerdike%baseLevel  = CrossDefFrom%summerdike%baseLevel
      CopyCrossDef%summerdike%flowArea   = CrossDefFrom%summerdike%flowArea
      CopyCrossDef%summerdike%totalArea  = CrossDefFrom%summerdike%totalArea
   endif
       
   !*** ground layer data
   allocate(CopyCrossDef%groundlayer)
   CopyCrossDef%groundlayer%used      = CrossDefFrom%groundlayer%used
   CopyCrossDef%groundlayer%thickness = CrossDefFrom%groundlayer%thickness
   CopyCrossDef%groundlayer%area      = CrossDefFrom%groundlayer%area
   CopyCrossDef%groundlayer%perimeter = CrossDefFrom%groundlayer%perimeter
   CopyCrossDef%groundlayer%width     = CrossDefFrom%groundlayer%width
end function CopyCrossDef

type(t_crsu) function CopyCrossConv(CrossConvFrom)

   ! Returns a Copy from the given CrossSection Conveyance
   ! DEALLOCATE this Copy after Use!!!!!!!!!!!!!!!!!!!!!!!!!!

   type(t_crsu) :: CrossConvFrom
   
   CopyCrossConv%jopen      = CrossConvFrom%jopen
   CopyCrossConv%msec       = CrossConvFrom%msec
   CopyCrossConv%iolu       = CrossConvFrom%iolu
   CopyCrossConv%negcon     = CrossConvFrom%negcon
   CopyCrossConv%conveyType = CrossConvFrom%conveyType

   CopyCrossConv%a_pos_extr = CrossConvFrom%a_pos_extr
   CopyCrossConv%a_neg_extr = CrossConvFrom%a_neg_extr
   CopyCrossConv%b_pos_extr = CrossConvFrom%b_pos_extr
   CopyCrossConv%b_neg_extr = CrossConvFrom%b_neg_extr

   CopyCrossConv%nrhh       = CrossConvFrom%nrhh
   CopyCrossConv%iolh       = CrossConvFrom%iolh
   CopyCrossConv%bob        = CrossConvFrom%bob
      
   CopyCrossConv%chezy_act  = CrossConvFrom%chezy_act
      
   CopyCrossConv%nru        = CrossConvFrom%nru

   if (CrossConvFrom%nru > 0) then
      
      if (allocated(CrossConvFrom%hu)) then
         allocate(CopyCrossConv%hu(CrossConvFrom%nru))
         CopyCrossConv%hu = CrossConvFrom%hu
      endif
         
      if (allocated(CrossConvFrom%af)) then
         allocate(CopyCrossConv%af(CrossConvFrom%nru))
         CopyCrossConv%af = CrossConvFrom%af
      endif
         
      if (allocated(CrossConvFrom%wf)) then
         allocate(CopyCrossConv%wf(CrossConvFrom%nru))
         CopyCrossConv%wf = CrossConvFrom%wf
      endif
         
      if (allocated(CrossConvFrom%pf)) then
         allocate(CopyCrossConv%pf(CrossConvFrom%nru))
         CopyCrossConv%pf = CrossConvFrom%pf
      endif
         
      if (allocated(CrossConvFrom%co1)) then
         allocate(CopyCrossConv%co1(CrossConvFrom%nru))
         CopyCrossConv%co1 = CrossConvFrom%co1
      endif
         
      if (allocated(CrossConvFrom%co2)) then
         allocate(CopyCrossConv%co2(CrossConvFrom%nru))
         CopyCrossConv%co2 = CrossConvFrom%co2
      endif
         
      if (allocated(CrossConvFrom%cz1)) then
         allocate(CopyCrossConv%cz1(CrossConvFrom%nru))
         CopyCrossConv%cz1 = CrossConvFrom%cz1
      endif
         
      if (allocated(CrossConvFrom%cz2)) then
         allocate(CopyCrossConv%cz2(CrossConvFrom%nru))
         CopyCrossConv%cz2 = CrossConvFrom%cz2
      endif
         
      if (allocated(CrossConvFrom%hh)) then
         allocate(CopyCrossConv%hh(CrossConvFrom%nru, 2))
         CopyCrossConv%hh = CrossConvFrom%hh
      endif
         
      if (allocated(CrossConvFrom%at)) then
         allocate(CopyCrossConv%at(CrossConvFrom%nru, 2))
         CopyCrossConv%at = CrossConvFrom%at
      endif
         
      if (allocated(CrossConvFrom%wt)) then
         allocate(CopyCrossConv%wt(CrossConvFrom%nru, 2))
         CopyCrossConv%wt = CrossConvFrom%wt
      endif
         
   endif

   
end function CopyCrossConv

double precision function getBobInterpolate(cross1, cross2, factor)

   type (t_CrossSection), intent(in)     :: cross1       !< cross section definition
   type (t_CrossSection), intent(in)     :: cross2       !< cross section definition
   double precision, intent(in)          :: factor       !< cross = (1 - f) * cross1 + f * cross2
   
   double precision bob1, bob2
   
   bob1 = getBobSingle(cross1)
   bob2 = getBobSingle(cross2)
   
   getBobInterpolate = (1d0 - factor) * bob1 + factor * bob2
   
end function getBobInterpolate
   
double precision function getBobSingle(cross)

   type (t_CrossSection), intent(in)     :: cross      !< cross section
   
   type(t_CSType), pointer               :: pCrossDef
   
   getBobSingle = cross%bedLevel

   pCrossDef => cross%tabDef
   if (associated(pCrossDef%groundlayer)) then
     if (pCrossDef%groundlayer%used) then
        getBobSingle = getBobSingle + pCrossDef%groundlayer%thickness
     endif
   endif

end function getBobSingle

double precision function getGroundLayer(cross)

   type (t_CrossSection), intent(in)     :: cross      !< cross section
   
   type(t_CSType), pointer               :: pCrossDef
   
   pCrossDef => cross%tabDef
   if (associated(pCrossDef%groundlayer)) then
      if (pCrossDef%groundlayer%used) then
         getGroundLayer = pCrossDef%groundlayer%thickness
      else
         getGroundLayer = 0.0d0
      endif
   else
      getGroundLayer = 0.0d0
   endif

end function getGroundLayer

double precision function GetCriticalDepth(q, cross)

   implicit none

   type(t_crosssection)            :: cross
   double precision, intent(in)    :: q     !< Discharge for which the critical depth should be computed.
   !
   ! Local variables
   !
   double precision               :: depth
   double precision               :: dummy
   double precision               :: wWidth
   double precision               :: groundLayer
   double precision               :: step
   double precision               :: wArea 
   logical                        :: first

   !! executable statements -------------------------------------------------------
   !
   !
   !     dlg, pjo, 1999-07-07
   !     Compute the critical depth according to the discharge through
   !     the culvert with an iteration loop
   !
   !     Include Pluvius data space
   !
   !     argument
   !     local variables
   !
   !     pjo, 13-04-2000, ars 4952
   !------------------------------------------------------------------------------------------
   ! Iterations, finds dtpc 
   !------------------------------------------------------------------------------------------
   depth  = cross%charHeight * 0.5d0
   step  = depth
   first = .true.

   do while (first .or. step > 0.0001d0)
   
      first = .false.
      
      call GetCSParsFlow(cross, depth, 0.0d0, dummy, wArea, dummy, wWidth, dummy)        
   
      step = 0.5d0 * step
      
      if ((q * q * wWidth) - (wArea * wArea * wArea * gravity) > 0.0d0) then 
         depth = depth + step
      else
         depth = depth - step
      endif
   
   enddo

   groundLayer = getGroundLayer(cross)
   GetCriticalDepth = depth - groundLayer
    
end function GetCriticalDepth

subroutine admin_crs_def(definitions)
   
   type(t_CSDefinitionSet), intent(inout), target :: Definitions
      
   integer i
   character(len=idlen), dimension(:), pointer :: ids
      
   allocate(definitions%hashlist%id_list(definitions%Count))
   definitions%hashlist%id_count = definitions%Count
   ids => definitions%hashlist%id_list
      
   do i= 1, definitions%count
      ids(i) = definitions%CS(i)%id
   enddo
      
   call hashfill(definitions%hashlist)
      
end subroutine admin_crs_def
 
subroutine fill_hashtable_csdef(definitions)
   
   type(t_CSDefinitionSet), intent(inout), target :: definitions
      
   integer i
   character(len=idlen), dimension(:), pointer :: ids
      
   allocate(definitions%hashlist%id_list(definitions%Count))
   definitions%hashlist%id_count = definitions%Count
   ids => definitions%hashlist%id_list
      
   do i= 1, definitions%count
      ids(i) = definitions%CS(i)%id
   enddo
      
   call hashfill(definitions%hashlist)

end subroutine fill_hashtable_csdef

subroutine createTablesForTabulatedProfile(crossDef)
   
   type(t_CStype), intent(inout) :: crossDef
   
   ! local parameters
   integer                                           :: levelsCount
   double precision, dimension(crossDef%levelsCount) :: heights
   double precision, dimension(crossDef%levelsCount) :: widths
   double precision, dimension(0:3)                  :: widthplains
  
   double precision  :: d1
   double precision  :: d2
   double precision  :: e1
   double precision  :: e2
   double precision  :: wl
   integer           :: ilev, isec

   crossDef%width_sub    = 0.0d0
   crossDef%af_sub       = 0.0d0
   crossDef%perim_sub    = 0.0d0
   crossDef%flowArea     = 0.0d0
   crossDef%wetPerimeter = 0.0d0

   crossDef%totalArea = 0.0d0
   crossDef%area_min  = 0.0d0
   crossDef%width_min = 0.0d0
   
   levelsCount = crossDef%levelsCount
   heights = crossDef%height   
   
   widthplains = 0d0

   if (crossDef%plains(1) == 0d0) then
      crossDef%plains(1) = huge(1d0)
      crossDef%plainsLocation(1) = crossDef%levelsCount
   endif
   widths = crossDef%flowWidth
   do isec = 1, 3
      widthplains(isec) = widthplains(isec-1) + crossDef%plains(isec)
   enddo
   
   d2 = 0.0D0
   e2 = widths(1)
   crossDef%af_sub = 0.0d0
   crossDef%perim_sub = 0.0d0
  
   do isec = 1, 3
      crossDef%af_sub(isec, 1) = 0d0
      crossDef%width_sub(isec, 1) = min(max(0d0, widths(1)-widthplains(isec-1)), widthplains(isec)-widthplains(isec-1))
      crossDef%perim_sub(isec, 1) = min(max(0d0, widths(1)-widthplains(isec-1)), widthplains(isec)-widthplains(isec-1))
      
      if (crossDef%plains(isec) <= 0.0d0) cycle

      do ilev = 2, levelsCount
      
         d1 = heights(ilev - 1) - heights(1)
         d2 = heights(ilev) - heights(1)
         e1 = min(max(0d0, widths(ilev - 1)-widthplains(isec-1)), widthplains(isec)-widthplains(isec-1))
         e2 = min(max(0d0, widths(ilev)-widthplains(isec-1)), widthplains(isec)-widthplains(isec-1))

         crossDef%width_sub(isec, ilev) = e2
         crossDef%af_sub(isec, ilev)    = crossDef%af_sub(isec, ilev - 1) + 0.5d0 * (d2 - d1) * (e1 + e2)
         
         if (crossDef%plainsLocation(2) > 0) then
            if (widths(ilev - 1) < widthplains(isec - 1)) then
               crossDef%perim_sub(isec, ilev) = 0.0d0
            elseif (widths(ilev - 1) >= widthplains(isec)) then 
               crossDef%perim_sub(isec, ilev) = crossDef%perim_sub(isec, ilev - 1)
            else
               crossDef%perim_sub(isec, ilev) = crossDef%perim_sub(isec, ilev - 1) + 2.0d0 * dsqrt(0.25d0 * (e2 - e1)**2 + (d2 - d1)**2)
            endif
         else
            crossDef%perim_sub(isec, ilev) = crossDef%perim_sub(isec, ilev - 1) + 2.0d0 * dsqrt(0.25d0 * (e2 - e1)**2 + (d2 - d1)**2)
         endif
            
      enddo

   enddo
   
   ! Totalize
   do isec = 1, 3
      do ilev = 1, crossDef%levelsCount
         crossDef%flowArea(ilev)  = crossDef%flowArea(ilev)  + crossDef%af_sub(isec, ilev)
         crossDef%wetPerimeter(ilev) = crossDef%wetPerimeter(ilev) + crossDef%perim_sub(isec, ilev)
      enddo
   enddo
      

   ! Calculation of Total Area
   levelsCount = crossDef%levelsCount
   heights = crossDef%height   
   widths = crossDef%totalWidth
   
   wl = widths(1)
   d2 = 0.0D0
   e2 = widths(1)
   
   crossDef%totalArea = 0.0d0
   crossDef%area_min  = 0.0d0
   crossDef%width_min = 0.0d0
   
   do ilev = 2, levelsCount
   
      d1 = heights(ilev - 1) - heights(1)
      d2 = heights(ilev) - heights(1)
      e1 = widths(ilev - 1)
      e2 = widths(ilev)

      crossDef%totalArea(ilev) = crossDef%totalArea(ilev - 1) + 0.5d0 * (d2 - d1) * (e1 + e2)
      
      crossDef%area_min(ilev)  = crossDef%area_min(ilev - 1) + (crossDef%width_min(ilev - 1) + 0.5d0 * max(e1- e2, 0.0d0)) * (d2 - d1)
      crossDef%width_min(ilev) = crossDef%width_min(ilev - 1) + max(e1 - e2, 0.0d0)

   enddo
   
end subroutine createTablesForTabulatedProfile

end module m_CrossSections
