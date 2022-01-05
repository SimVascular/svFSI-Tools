!***********************************************************************

!     Module definition for mesh

      module mesh
      use stdParams
      implicit none

!     Constant parameters for various element types
      integer, parameter :: eType_NA = 100, eType_LIN1 = 101,  &
        eType_LIN2 = 102, eType_TRI3 = 103, eType_TRI6 = 104,  &
        eType_QUD4 = 105, eType_QUD8 = 106, eType_QUD9 = 107,  &
        eType_TET4 = 108, eType_TET10 = 109, eType_HEX8 = 110, &
        eType_HEX20 = 111, eType_HEX27 = 112

!     Connectivity type (CSR format)
      type connType
         integer :: nnz = 0
         integer, allocatable :: pcol(:)
         integer, allocatable :: prow(:)
      END TYPE connType

!     Face type
      type faceType
!        Element type
         integer :: eType = eType_NA
!        Spatial dimension
         integer nsd
!        Num face nodes
         integer nNo
!        Num face elements
         integer nEl
!        Num of nodes per element
         integer eNoN
!        Num of edges per element
         integer eNoE
!        VTK element type
         integer vtkType
!        Global node : pointer to mesh node
         integer, allocatable :: gN(:)
!        Global node : pointer to mesh element
         integer, allocatable :: gE(:)
!        Face element connectivity array
         integer, allocatable :: IEN(:,:)
!        Edge ordering of a face element
         integer, allocatable :: eord(:,:)
!        Position coordinates
         real(kind=8), allocatable :: x(:,:)
!        Face name
         character(len=strL) :: fname
      end type faceType

!     Mesh type
      type meshType
!        Mesh element type
         integer :: eType = eType_NA
!        Spatial dimension
         integer nsd
!        Num mesh nodes
         integer nNo
!        Num mesh elements
         integer nEl
!        Num nodes per mesh element
         integer eNoN
!        Num edges per mesh element
         integer eNoE
!        Num faces per mesh element
         integer eNoF
!        Num nodes on face of a mesh element
         integer fNoN
!        VTK element type
         integer vtkType
!        Num faces on each mesh element
         integer nFa
!        Edge ordering of a mesh elements
         integer, allocatable :: eOrd(:,:)
!        Face ordering of a mesh element
         integer, allocatable :: fOrd(:,:)
!        Mesh element connectivity array
         integer, allocatable :: IEN(:,:)
!        Mesh position coordinates
         real(kind=8), allocatable :: x(:,:)
!        Dynamically allocated face types
         type(faceType), allocatable :: fa(:)
!        Mesh name
         character(len=strL) :: fname
!        Edge connectivity
         type(connType) :: econ
!        Face connectivity
         type(connType) :: fcon
      end type meshType

!     Higher order element type
      type hoeType
!        Element type
         integer :: eType = eType_NA
!        Element num of nodes
         integer eNoN
!        Element num of edges
         integer eNoE
!        Element num of faces
         integer eNoF
!        Face num of nodes
         integer fNoN
!        Edge ordering
         integer, allocatable :: eord(:,:)
!        Face ordering
         integer, allocatable :: ford(:,:)
      end type hoeType

!     Interface for selecting mesh element, face element and a higher
!     order element
      interface selectEle
         module procedure selectEleMesh, selectEleFace, selectEleHO
      end interface

!     Destructor interface for connectivity, face and mesh
      interface destroy
         module procedure destroyConn, destroyFace, destroyMesh
      end interface

!     Deep copy interface for connectivity, face and mesh
      interface deepCopy
         module procedure deepCopyConn, deepCopyFace, deepCopyMesh
      end interface

      contains
      !=============================================================
!        Select mesh element type
         subroutine selectEleMesh(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         integer :: iFa

!        select mesh element type
         if (lM%nsd .eq. 2) then
            select case (lM%eNoN)
            case (3)
               lM%eType   = eType_TRI3
               lM%eNoE    = 3
               lM%eNoF    = 3
               lM%fNoN    = 2
               lM%vtkType = 5

            case (6)
               lM%eType   = eType_TRI6
               lM%eNoE    = 3
               lM%eNoF    = 3
               lM%fNoN    = 3
               lM%vtkType = 22

            case (4)
               lM%eType   = eType_QUD4
               lM%eNoE    = 4
               lM%eNoF    = 4
               lM%fNoN    = 2
               lM%vtkType = 9

            case (8)
               lM%eType   = eType_QUD8
               lM%eNoE    = 4
               lM%eNoF    = 4
               lM%fNoN    = 3
               lM%vtkType = 23

            case (9)
               lM%eType   = eType_QUD9
               lM%eNoE    = 4
               lM%eNoF    = 4
               lM%fNoN    = 3
               lM%vtkType = 28

            case default
               write(stdout,ftab4) &
                  "ERROR: mesh element type not defined"
               STOP
            end select

         else if (lM%nsd .eq. 3) then
            select case (lM%eNoN)
            case (3)
               lM%eType   = eType_TRI3
               lM%eNoE    = 3
               lM%eNoF    = 3
               lM%fNoN    = 2
               lM%vtkType = 5

            case (4)
               lM%eType   = eType_TET4
               lM%eNoE    = 6
               lM%eNoF    = 4
               lM%fNoN    = 3
               lM%vtkType = 10

            case (10)
               lM%eType   = eType_TET10
               lM%eNoE    = 6
               lM%eNoF    = 4
               lM%fNoN    = 6
               lM%vtkType = 24

            case (8)
               lM%eType   = eType_HEX8
               lM%eNoE    = 12
               lM%eNoF    = 6
               lM%fNoN    = 4
               lM%vtkType = 12

            case (20)
               lM%eType   = eType_HEX20
               lM%eNoE    = 12
               lM%eNoF    = 6
               lM%fNoN    = 8
               lM%vtkType = 25

            case (27)
               lM%eType   = eType_HEX27
               lM%eNoE    = 12
               lM%eNoF    = 6
               lM%fNoN    = 9
               lM%vtkType = 29

            case default
               write(stdout,ftab4) &
                  "ERROR: mesh element type not defined"
               STOP
            end select
         end if

!        Set nodal orderings on edges and faces of an element
         allocate(lM%eOrd(2,lM%eNoE), lM%fOrd(lM%fNoN,lM%eNoF))
         call setEdgeNodeOrder(lM%eType, lM%eNoE, lM%eord)
         call setFaceNodeOrder(lM%eType, lM%fNoN, lM%eNoF, lM%ford)

!        select face element type
         do iFa=1, lM%nFa
            call selectEleFace(lM%eType, lM%fa(iFa))
            if (lM%fNoN .ne. lM%fa(iFa)%eNoN) then
               write(stdout,ftab4) "ERROR: inconsistent face"// &
                  " element selected (nodes per face)"
               stop
            end if
         end do ! iFa

         return
         end subroutine selectEleMesh
      !=============================================================
!        Select face element type
         subroutine selectEleFace(eType, lFa)
         implicit none
         integer, intent(in) :: eType
         type(faceType), intent(inout) :: lFa

         select case (eType)
         case (eType_TRI3, eType_QUD4)
            lFa%eType   = eType_LIN1
            lFa%eNoN    = 2
            lFa%eNoE    = 1
            lFa%vtkType = 3

         case (eType_TRI6, eType_QUD8, eType_QUD9)
            lFa%eType   = eType_LIN2
            lFa%eNoN    = 3
            lFa%eNoE    = 1
            lFa%vtkType = 21

         case (eType_TET4)
            lFa%eType   = eType_TRI3
            lFa%eNoN    = 3
            lFa%eNoE    = 3
            lFa%vtkType = 5

         case (eType_TET10)
            lFa%eType   = eType_TRI6
            lFa%eNoN    = 6
            lFa%eNoE    = 3
            lFa%vtkType = 22

         case (eType_HEX8)
            lFa%eType   = eType_QUD4
            lFa%eNoN    = 4
            lFa%eNoE    = 4
            lFa%vtkType = 9

         case (eType_HEX20)
            lFa%eType   = eType_QUD8
            lFa%eNoN    = 8
            lFa%eNoE    = 4
            lFa%vtkType = 23

         case (eType_HEX27)
            lFa%eType   = eType_QUD9
            lFa%eNoN    = 9
            lFa%eNoE    = 4
            lFa%vtkType = 28

         case default
            write(stdout,ftab4) "ERROR: face element type not defined"
            STOP
         end select

         allocate(lFa%eord(2,lFa%eNoE))
         call setEdgeNodeOrder(lFa%eType, lFa%eNoE, lFa%eord)

         return
         end subroutine selectEleFace
      !=============================================================
!        Select higher order element
         subroutine selectEleHO(eType, hoe)
         implicit none
         integer, intent(in) :: eType
         type(hoeType), intent(out) :: hoe

         integer i

         hoe%eNoE = 0
         hoe%eNoF = 0
         hoe%fNoN = 0

         select case (eType)
         case (eType_TRI3)
            write(stdout,ftab2) "Converting TRI3 mesh to TRI6.. "
            hoe%eType = eType_TRI6
            hoe%eNoN  = 6
            hoe%eNoE  = 3

            allocate(hoe%eord(2,hoe%eNoE))
            hoe%eord = reshape((/1,2,2,3,3,1/), shape(hoe%eord))

         case (eType_TET4)
            write(stdout,ftab2) "Converting TET4 mesh to TET10.. "
            hoe%eType = eType_TET10
            hoe%eNoN  = 10
            hoe%eNoE  = 6

            allocate(hoe%eord(2,hoe%eNoE))
            hoe%eord = reshape((/1,2,2,3,3,1,1,4,2,4,3,4/), &
               shape(hoe%eord))

         case (eType_QUD4)
            write(stdout,ftab1,advance='no') "Convert to QUAD8 (1) "// &
               "or QUAD9 (2)? "
            read(*,*) i

            if (i .eq. 1) then
               write(stdout,ftab2) "Converting QUAD4 mesh to QUAD8.. "
               hoe%eType = eType_QUD8
               hoe%eNoN  = 8
               hoe%eNoE  = 4

               allocate(hoe%eord(2,hoe%eNoE))
               hoe%eord = reshape((/1,2,2,3,3,4,4,1/), shape(hoe%eord))

            else if (i .eq. 2) then
               write(stdout,ftab2) "Converting QUAD4 mesh to QUAD9.. "
               hoe%eType = eType_QUD9
               hoe%eNoN  = 9
               hoe%eNoE  = 4

               allocate(hoe%eord(2,hoe%eNoE))
               hoe%eord = reshape((/1,2,2,3,3,4,4,1/), shape(hoe%eord))

            else
               write(stdout,ftab4) "ERROR: unknown option"
               return
            end if

         case (eType_HEX8)
            write(stdout,ftab1,advance='no') "Convert to HEX20 (1) "// &
               "or HEX27 (2)? "
            read(*,*) i
            if (i .eq. 1) then
               write(stdout,ftab2) "Converting HEX8 mesh to HEX20.. "
               hoe%eType = eType_HEX20
               hoe%eNoN  = 20
               hoe%eNoE  = 12

               allocate(hoe%eord(2,hoe%eNoE))
               hoe%eord = reshape((/1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5, &
        &         1,5,2,6,3,7,4,8/), shape(hoe%eord))

            else if (i .eq. 2) then
               write(stdout,ftab2) "Converting HEX8 mesh to HEX27.. "
               hoe%eType = eType_HEX27
               hoe%eNoN  = 27
               hoe%eNoE  = 12
               hoe%eNoF  = 6
               hoe%fNoN  = 4

               allocate(hoe%eord(2,hoe%eNoE), hoe%ford(hoe%fNoN,hoe%eNoF))
               hoe%eord = reshape((/1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5, &
        &         1,5,2,6,3,7,4,8/), shape(hoe%eord))
               hoe%ford = reshape((/4,1,5,8,2,3,7,6,1,2,6,5,3,4,8,7, &
        &         2,1,4,3,5,6,7,8/), shape(hoe%ford))
            else
               write(stdout,ftab4) "ERROR: unknown option"
               return
            end if

         end select

         return
         end subroutine selectEleHO
      !=============================================================
!        Set node ordering on the edges of an element
         subroutine setEdgeNodeOrder(eType, eNoE, eord)
         implicit none
         integer, intent(in) :: eType, eNoE
         integer, intent(inout) :: eord(2,eNoE)

         select case (eType)
         case (eType_LIN1, eType_LIN2)
            eOrd(1,1) = 1; eOrd(2,1) = 2

         case (eType_TRI3, eType_TRI6)
            eOrd = reshape((/1,2,2,3,3,1/), shape(eOrd))

         case (eType_QUD4, eType_QUD8, eType_QUD9)
            eOrd = reshape((/1,2,2,3,3,4,4,1/), shape(eOrd))

         case (eType_TET4, eType_TET10)
            eOrd = reshape((/1,2,2,3,3,1,1,4,2,4,3,4/), shape(eOrd))

         case (eType_HEX8, eType_HEX20, eType_HEX27)
            eOrd = reshape((/1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5, &
        &      1,5,2,6,3,7,4,8/), shape(eOrd))

         end select

         return
         end subroutine setEdgeNodeOrder
      !=============================================================
!        Set node ordering on the faces of an element
         subroutine setFaceNodeOrder(eType, fNoN, eNoF, ford)
         implicit none
         integer, intent(in) :: eType, fNoN, eNoF
         integer, intent(inout) :: ford(fNoN,eNoF)

         select case (eType)
         case (eType_TRI3)
            fOrd = reshape((/1,2,2,3,3,1/), shape(fOrd))

         case (eType_TRI6)
            fOrd = reshape((/1,2,4,2,3,5,3,1,6/), shape(fOrd))

         case (eType_QUD4)
            fOrd = reshape((/1,2,2,3,3,4,4,1/), shape(fOrd))

         case (eType_QUD8, eType_QUD9)
            fOrd = reshape((/1,2,5,2,3,6,3,4,7,4,1,8/), shape(fOrd))

         case (eType_TET4)
            fOrd = reshape((/1,2,3,2,1,4,1,3,4,3,2,4/), shape(fOrd))

         case (eType_TET10)
            fOrd = reshape((/1,2,3,5,6,7,2,1,4,5,8,9, &
        &      1,3,4,7,10,8,3,2,4,6,9,10/), shape(fOrd))

         case (eType_HEX8)
            fOrd = reshape((/4,1,5,8,2,3,7,6,1,2,6,5,3,4,8,7, &
        &      2,1,4,3,5,6,7,8/), shape(ford))
         end select

         return
         end subroutine setFaceNodeOrder
      !=============================================================
         subroutine destroyConn(conn)
         implicit none
         type(connType), intent(inout) :: conn

         if (allocated(conn%prow)) deallocate(conn%prow)
         if (allocated(conn%pcol)) deallocate(conn%pcol)

         conn%nnz = 0

         return
         end subroutine destroyConn
      !=============================================================
         subroutine destroyFace(fa)
         implicit none
         type(faceType), intent(inout) :: fa

         if (allocated(fa%gN))   deallocate(fa%gN)
         if (allocated(fa%gE))   deallocate(fa%gE)
         if (allocated(fa%IEN))  deallocate(fa%IEN)
         if (allocated(fa%eord)) deallocate(fa%eord)
         if (allocated(fa%x))    deallocate(fa%x)

         fa%eType = eType_NA

         return
         end subroutine destroyFace
      !=============================================================
         subroutine destroyMesh(msh)
         implicit none
         type(meshType), intent(inout) :: msh

         integer i

         if (allocated(msh%IEN))  deallocate(msh%IEN)
         if (allocated(msh%eord)) deallocate(msh%eord)
         if (allocated(msh%ford)) deallocate(msh%ford)
         if (allocated(msh%x))    deallocate(msh%x)

         call destroyConn(msh%econ)
         call destroyConn(msh%fcon)
         do i=1, msh%nFa
            call destroyFace(msh%fa(i))
         end do

         if (allocated(msh%fa)) deallocate(msh%fa)

         msh%eType = eType_NA

         return
         end subroutine destroyMesh
      !=============================================================
         subroutine deepCopyConn(lcn, gcn)
         implicit none
         type(connType), intent(in) :: lcn
         type(connType), intent(inout) :: gcn

         integer n

         gcn%nnz = lcn%nnz
         if (allocated(lcn%prow)) then
            n = size(lcn%prow)
            allocate(gcn%prow(n))
            gcn%prow = lcn%prow
         end if

         if (allocated(lcn%pcol)) then
            allocate(gcn%pcol(gcn%nnz))
            gcn%pcol = lcn%pcol
         end if

         return
         end subroutine deepCopyConn
      !=============================================================
         subroutine deepCopyFace(lFa, gFa)
         implicit none
         type(faceType), intent(in) :: lFa
         type(faceType), intent(inout) :: gFa

         if (gFa%eType .eq. eType_NA) then
            write(stdout,ftab4) "ERROR: element type not set to "// &
               "copy face structure"
            stop
         end if

         gFa%nsd   = lFa%nsd
         gFa%nNo   = lFa%nNo
         gFa%nEl   = lFa%nEl
         gFa%fname = lFa%fname

         allocate(gFa%gN(gFa%nNo))
         allocate(gFa%gE(gFa%nEl))
         allocate(gFa%x(gFa%nsd,gFa%nNo))
         allocate(gFa%IEN(gFa%eNoN,gFa%nEl))

         gFa%x   = lFa%x
         gFa%gN  = lFa%gN
         gFa%gE  = lFa%gE
         gFa%IEN = lFa%IEN

         return
         end subroutine deepCopyFace
      !=============================================================
         subroutine deepCopyMesh(lM, gM)
         implicit none
         type(meshType), intent(in) :: lM
         type(meshType), intent(inout) :: gM

         integer i

         gM%nsd   = lM%nsd
         gM%nNo   = lM%nNo
         gM%nEl   = lM%nEl
         gM%eNoN  = lM%eNoN
         gM%nFa   = lM%nFa
         gM%fname = lM%fname
         allocate(gM%fa(gM%nFa))

         call selectEleMesh(gM)

         allocate(gM%x(gM%nsd,gM%nNo), gM%IEN(gM%eNoN,gM%nEl))
         gM%x   = lM%x
         gM%IEN = lM%IEN

         if (lM%econ%nnz .ne. 0) call deepCopyConn(lM%econ, gM%econ)
         if (lM%fcon%nnz .ne. 0) call deepCopyConn(lM%fcon, gM%fcon)

         do i=1, gM%nFa
            call deepCopyFace(lM%fa(i), gM%fa(i))
         end do

         return
         end subroutine deepCopyMesh
      !=============================================================
    end module mesh

!***********************************************************************

