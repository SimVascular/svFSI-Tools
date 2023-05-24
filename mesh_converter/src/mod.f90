!***********************************************************************


      module varmod
      use stdParams
      use genUtils
      use mesh
      use meshUtils

      end module varmod

!***********************************************************************

!--------------------------------------------------------------------
!
!     GMSH data structures interfaced.
!
!--------------------------------------------------------------------

      module gmshMod

      type gmshPhysicalNameType
         integer :: num
         integer,allocatable :: dimension(:)
         integer,allocatable :: physicalTag(:)
         character(len=80),allocatable :: name(:)
      end type gmshPhysicalNameType

      type gmshEntityType
         integer :: numPoints
         integer :: numCurves
         integer :: numSurfaces
         integer :: numVolumes

         integer, allocatable :: Point_Tags(:)
         integer, allocatable :: Curve_Tags(:)
         integer, allocatable :: Surface_Tags(:)
         integer, allocatable :: Volume_Tags(:)

         integer, allocatable :: Point_numPhysicalTags(:)
         integer, allocatable :: Curve_numPhysicalTags(:)
         integer, allocatable :: Surface_numPhysicalTags(:)
         integer, allocatable :: Volume_numPhysicalTags(:)

         integer, allocatable :: Point_physicalTags(:,:)
         integer, allocatable :: Curve_physicalTags(:,:)
         integer, allocatable :: Surface_physicalTags(:,:)
         integer, allocatable :: Volume_physicalTags(:,:)
      end type gmshEntityType

      type gmshNodeType
         integer :: numNodeBlocks
         integer :: numNodes
         integer :: minNodeTag
         integer :: maxNodeTag
         integer,allocatable :: numNodesInBlock(:)
         real(kind=8),allocatable :: coord(:,:)
      end type gmshNodeType

      type gmshElementType
         integer :: numElementBlocks
         integer :: numElements
         integer :: minElementTag
         integer :: maxElementTag

         integer,allocatable :: ElementDim(:)
         integer,allocatable :: EntityTag(:)
         integer,allocatable :: eNoN(:)
         integer,allocatable :: numElementsInBlock(:)
         integer,allocatable :: conn(:,:)
      end type gmshElementType

      ! Variables
      integer, allocatable :: g2l(:)
      type(gmshPhysicalNameType) :: gmshPhysicalNames
      type(gmshNodeType)         :: gmshNodes
      type(gmshElementType)      :: gmshElements
      type(gmshEntityType)       :: gmshEntities

      end module gmshMod
