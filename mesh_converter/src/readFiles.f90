!***********************************************************************

      subroutine readGambitNeu(lM)
      use varmod
      implicit none

      type(meshType), intent(out)  :: lM

      integer :: a, b, e, i, iFa, Ac, Ec, fid, nToks
      character(len=strL) :: rLine, tokenList(maxToks)

      integer, allocatable :: ptr(:), gmap(:)

      fid = 100
      write(stdout,ftab1) "Loading file "//TRIM(lM%fname)
      open(fid, file=trim(lM%fname))
      call findKwrd(fid, "NUMNP")
      read(fid,*) lM%nNo, lM%nEl, i, lM%nFa, a, b

      lM%nsd = max(a,b)
      write(stdout,ftab2) "Nsd: "//TRIM(STR(lM%nsd))
      write(stdout,ftab2) "nNo: "//TRIM(STR(lM%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(lM%nEl))
      write(stdout,ftab2) "nFa: "//TRIM(STR(lM%nFa))
      allocate(lM%fa(lM%nFa))

!     Read nodal coordinates
      allocate(lM%x(lM%nsd,lM%nNo))
      call findKwrd(fid, "NODAL")
      do a=1, lM%nNo
         read(fid,*) i, lM%x(:,a)
      end do

!     Determing element connectivity and set element type
      call findKwrd(fid, "ELEMENTS/CELLS")
      read(fid,'(A)') rLine
      call parseString(rLine, tokenList, nToks)
      if (nToks .eq. 0) then
         write(stdout,ftab4) &
            "ERROR: could not parse element connectivity"
         STOP
      end if
      read(tokenList(3),*) lM%eNoN

!     Select mesh and face element types
      call selectele(lM)

!     Read connectivity (IEN)
      allocate(lM%IEN(lM%eNoN, lM%nEl))
      if (nToks-3 .eq. lM%eNoN) then
         e = 1
         do a=1, lM%eNoN
            read(tokenList(3+a),*) lM%IEN(a,e)
         end do
         do e=2, lM%nEl
            read(fid,*) i, i, i, lM%IEN(:,e)
         end do
      else
         rewind(fid)
         call findKwrd(fid, "ELEMENTS/CELLS")
         do e=1, lM%nEl
            read(fid,*) i, i, i, lM%IEN(1:ntoks-3,e)
            read(fid,*) lM%IEN(ntoks-3+1:lM%eNoN,e)
         end do
      end if

!     Swap nodes for elements to be consistent with VTK node ordering
      allocate(gmap(lM%eNoF))
      do i=1, lM%eNoF
         gmap(i) = i
      end do

      if (lM%eType .eq. eType_TET4) then
         do e=1, lM%nEl
            call swap(lM%IEN(1,e), lM%IEN(2,e))
         end do

      else if (lM%eType .eq. eType_HEX8) then
         do e=1, lM%nEl
            call swap(lM%IEN(3,e), lM%IEN(4,e))
            call swap(lM%IEN(7,e), lM%IEN(8,e))
         end do
         gmap = (/3, 2, 4, 1, 5, 6/)

      end if

!     Read face element data
      do iFa=1, lM%nFa
         call findKwrd(fid, "BOUNDARY")
         read(fid,'(A)') rLine
         call parseString(rLine, tokenList, nToks)
         read(tokenList(1),*) lM%fa(iFa)%fname
         write(stdout,ftab3) "Face <"//TRIM(lM%fa(iFa)%fname)//">"
         if (TRIM(tokenList(2)) .ne. "1") then
            write(stdout,ftab4) &
               "ERROR: element information not found on boundary"
            STOP
         end if
         read(tokenList(3),*) lM%fa(iFa)%nEl
         allocate(lM%fa(iFa)%gE(lM%fa(iFa)%nEl))
         allocate(lM%fa(iFa)%IEN(lM%fa(iFa)%eNoN, lM%fa(iFa)%nEl))
         do e=1, lM%fa(iFa)%nEl
            read(fid,*) lM%fa(iFa)%gE(e), b, i
            Ec = lM%fa(iFa)%gE(e)
            i  = gmap(i)
            do a=1, lM%fNoN
               b = lM%fOrd(a,i)
               lM%fa(iFa)%IEN(a,e) = lM%IEN(b,Ec)
            end do
         end do
         read(fid,*)
      end do
      close(fid)

!     Setup face data stucture (fa%x, fa%gN)
      allocate(ptr(lM%nNo))
      do iFa=1, lM%nFa
         ptr = 0
         lM%fa(iFa)%nNo = 0
         do e=1, lM%fa(iFa)%nEl
            do a=1, lM%fa(iFa)%eNoN
               Ac = lM%fa(iFa)%IEN(a,e)
               if (ptr(Ac) .eq. 0) then
                  lM%fa(iFa)%nNo = lM%fa(iFa)%nNo + 1
                  ptr(Ac) = 1
               end if
            end do
         end do

         lM%fa(iFa)%nsd = lM%nsd
         allocate(lM%fa(iFa)%gN(lM%fa(iFa)%nNo))
         allocate(lM%fa(iFa)%x(lM%fa(iFa)%nsd,lM%fa(iFa)%nNo))
         a = 0
         do Ac=1, lM%nNo
            if (ptr(Ac) .gt. 0) then
               a = a + 1
               lM%fa(iFa)%gN(a) = Ac
               lM%fa(iFa)%x(:,a) = lM%x(:,Ac)
            end if
         end do
      end do

      deallocate(ptr, gmap)

      lM%fname = "mesh-complete.mesh"

      return
      end subroutine readGambitNeu

!***********************************************************************

      subroutine readInputFile(fName, lM)
      use mesh
      implicit none
      character(len=strL), intent(in) :: fname
      type(meshType), intent(inout) :: lM

      integer i, fid, istat
      character(len=strL) :: stmp

      fid = 1265
      istat = 0
      open(fid, file=trim(fName))
      OUTER_LOOP: do
         call findHashTag(fid, "numSpatialDim", istat)
         if (istat .ne. 0) exit
         read(fid,*) lM%nsd

         call findHashTag(fid, "meshFilePath", istat)
         if (istat .ne. 0) exit
         read(fid,'(A)') lM%fname

         call findHashTag(fid, "faceFilesPaths", istat)
         if (istat .ne. 0) exit
         lM%nFa = 0
         do
            read(fid,'(A)',iostat=i) stmp
            if (i .ne. 0) exit
            if (len(trim(stmp)) .gt. 0) lM%nFa = lM%nFa + 1
         end do

         allocate(lM%fa(lM%nFa))

         call findHashTag(fid, "faceFilesPaths", istat)
         if (istat .ne. 0) exit
         do i=1, lM%nFa
            read(fid,'(A)') lM%fa(i)%fname
         end do

         exit OUTER_LOOP
      end do OUTER_LOOP
      close(fid)

      if (istat .ne. 0) then
         write(stdout,ftab4) "ERROR: reading inputs"
         stop
      end if

      return
      contains
         !==========================================
         subroutine findHashTag(fileId, sTag, istat)
         implicit none
         integer, intent(in) :: fileId
         integer, intent(inout) :: istat
         character(len=*), intent(in) :: sTag

         integer :: kwrdL, slen
         character(len=strL) :: sLine

         rewind(fileId)

         istat = 0
         kwrdL = len(trim(sTag))
         do
            read(fileId,'(A)',end=001) sLine
            if (sLine(1:1) .eq. '#') then
               slen  = len(trim(sLine))
               sLine = sLine(2:slen)
               sLine = adjustl(sLine)
               if (sLine(1:kwrdL) .eq. trim(sTag)) return
            end if
         end do

 001     write(stdout,ftab4) "ERROR: EOF reached while finding "// &
         "the hashtag <"//trim(sTag)//">"
         istat = -1
         return

         end subroutine findHashTag
         !==========================================
      end subroutine readInputFile

!***********************************************************************

      subroutine readVTKMesh(lM)
      use varmod
      implicit none
      type(meshType), intent(inout) :: lM

      integer i, j, k

!     Read mesh from vtu file
      write(stdout, ftab1) "Reading mesh from vtu file   <----   "// &
         trim(lM%fName)
      call readVTU(lM, lM%fName)

!     Set mesh and face element types
      call selectele(lM)

!     Read faces from vtp files
      do i=1, lM%nFa
         lM%fa(i)%nsd = lM%nsd
         call readVTP(lM%fa(i), lM%fa(i)%fname)
      end do

      write(stdout,ftab2) "Nsd: "//TRIM(STR(lM%nsd))
      write(stdout,ftab2) "nNo: "//TRIM(STR(lM%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(lM%nEl))
      write(stdout,ftab2) "nFa: "//TRIM(STR(lM%nFa))

!     Reset mesh and face names for output
      j = scan(lM%fname, '/', back=.true.)
      k = scan(lM%fname, '.', back=.true.)
      lM%fname = trim(lM%fname(j+1:k-1))
      do i=1, lM%nFa
         j = scan(lM%fa(i)%fname, '/', back=.true.)
         k = scan(lM%fa(i)%fname, '.', back=.true.)
         lM%fa(i)%fname = trim(lM%fa(i)%fname(j+1:k-1))
      end do

      return
      end subroutine readVTKMesh

!***********************************************************************

      subroutine findKwrd(fid, kwrd)
      use stdParams
      use genUtils
      implicit none
      integer, intent(in) :: fid
      character(len=*), intent(in) :: kwrd

      integer :: itok, ntoks
      character(len=strL) :: rLine, tokenList(maxToks)

      MY_LOOP : do
         read(fid,'(A)',end=001) rLine
         call parseString(rLine, tokenList, nToks)
         if (nToks .gt. 0) then
            do itok=1, nToks
               if (TRIM(tokenList(itok)) .eq. TRIM(kwrd)) exit MY_LOOP
            end do
         end if
      end do MY_LOOP

      return

 001  write(stdout,ftab4) "ERROR: end of file reached.."
      write(stdout,ftab4) 'Failed processing for keyword "'// &
            trim(kwrd)//'"'
      STOP
      end subroutine findKwrd

!***********************************************************************

      subroutine readGmshMsh(lM)
      use varmod
      use gmshMod
      implicit none

      type(meshType), intent(out)  :: lM

      integer :: fid, iFa

      fid = 100
      write(stdout,ftab1) "Loading file "//TRIM(lM%fname)
      open(fid, file=TRIM(lM%fname))

      ! dimension
      write(stdout,ftab1,advance='no') &
      "Dimension of the body mesh (2/3):  "
      read(*,*) lM%nsd

      ! metadata
      call gmsh_readmeta(lM,fid)

      ! PhysicalNames
      call gmsh_readphysicalnames(lM,fid)

      ! Entities
      call gmsh_readentities(fid)

      ! Nodes
      ! Please note not all nodes are grid points, e.g. center of the circle.
      call gmsh_readnodes(fid)

      ! Elements
      call gmsh_readelements(fid)
      close(fid)

      ! Assemble body mesh
      call gmsh_bodymesh(lM)

      ! Assemble boundary mesh
      call gmsh_boundarymesh(lM)

      call selectele(lM)

      ! Output general statistics
      write(stdout,ftab2) "Nsd: "//TRIM(STR(lM%nsd))
      write(stdout,ftab2) "nNo: "//TRIM(STR(lM%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(lM%nEl))
      write(stdout,ftab2) "nFa: "//TRIM(STR(lM%nFa))
      do iFa=1, lM%nFa
         write(stdout,ftab3) "Face <"//TRIM(lM%fa(iFa)%fname)//">"
      end do

      lM%fname = "mesh-complete.mesh"

      return
      end subroutine readGmshMsh
!***********************************************************************
!     Read gmsh meta data
      subroutine gmsh_readmeta(lM,fid)
      use varmod
      use gmshMod

      implicit none

      type(meshType), intent(out)  :: lM
      integer, intent(inout) :: fid

      real(kind=8) :: rtemp
      integer :: i, j

      ! metadata
      rewind(fid)
      call findKwrd(fid, "$MeshFormat")
      read(fid,*) rtemp, i, j
      if ( i .ne. 0) then
         write(stdout,ftab4) "ERROR: "//TRIM(lM%fname)// &
         " is not stored in plain text."
         stop
      end if

      return
      end subroutine gmsh_readmeta
!***********************************************************************
!     Read PhysicalNames
      subroutine gmsh_readphysicalnames(lM, fid)
      use varmod
      use gmshMod

      implicit none

      type(meshType), intent(out)  :: lM
      integer, intent(inout) :: fid

      integer :: i, n

      rewind(fid)
      call findKwrd(fid, "$PhysicalNames")
      read(fid,*) n

      gmshPhysicalNames%num = n
      allocate(gmshPhysicalNames%dimension(n))
      allocate(gmshPhysicalNames%physicalTag(n))
      allocate(gmshPhysicalNames%name(n))
      do i = 1, gmshPhysicalNames%num
         read(fid,*) gmshPhysicalNames%dimension(i), gmshPhysicalNames%physicalTag(i), gmshPhysicalNames%name(i)
      end do

      if (lM%nsd .ne. MAXVAL(gmshPhysicalNames%dimension)) then
         write(stdout,"(14X,A,I5)") "nsd = ", lM%nsd
         write(stdout,"(14X,A,I5)") "MAXVAL(gmshPhysicalNames%dimension) = ", MAXVAL(gmshPhysicalNames%dimension)
         write(stdout,ftab4) "Both body mesh and surface mesh need to be assigned to a physical group."
         stop
      end if

      return
      end subroutine gmsh_readphysicalnames
!***********************************************************************
!     Read Nodes
      subroutine gmsh_readnodes(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, j, itemp
      integer, allocatable :: NodeTag(:)

      rewind(fid)
      call findKwrd(fid, "$Nodes")
      read(fid,*) gmshNodes%numNodeBlocks, gmshNodes%numNodes, &
                  gmshNodes%minNodeTag, gmshNodes%maxNodeTag

      allocate(gmshNodes%numNodesInBlock(gmshNodes%numNodeBlocks))
      allocate(NodeTag(gmshNodes%numNodes))
      allocate(gmshNodes%coord(3,gmshNodes%numNodes))

      ! read coordiantes block by block
      do i = 1, gmshNodes%numNodeBlocks
         read(fid,*) itemp, itemp, itemp, gmshNodes%numNodesInBlock(i)

         do j = 1, gmshNodes%numNodesInBlock(i)
            read(fid,*) NodeTag(j)
         end do
         do j = 1, gmshNodes%numNodesInBlock(i)
            read(fid,*) gmshNodes%coord(:,NodeTag(j))
         end do
      end do

      return
      end subroutine gmsh_readnodes
!***********************************************************************
!     Read Elements
      subroutine gmsh_readelements(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, j, itemp, ie, eNoN, maxeNoN

      rewind(fid)
      call findKwrd(fid, "$Elements")
      read(fid,*) gmshElements%numElementBlocks, gmshElements%numElements, &
                  gmshElements%minElementTag, gmshElements%maxElementTag

      allocate(gmshElements%ElementDim(gmshElements%numElementBlocks))
      allocate(gmshElements%EntityTag(gmshElements%numElementBlocks))
      allocate(gmshElements%eNoN(gmshElements%numElementBlocks))
      allocate(gmshElements%numElementsInBlock(gmshElements%numElementBlocks))

      ! find maximum eNoN
      do i = 1, gmshElements%numElementBlocks
         read(fid,*) gmshElements%ElementDim(i), gmshElements%EntityTag(i), &
                     itemp, gmshElements%numElementsInBlock(i)
         select case (itemp)
         case(1)
            gmshElements%eNoN(i) = 2 ! 2-node line.
         case(2)
            gmshElements%eNoN(i) = 3 ! 3-node triangle.
         case(3)
            gmshElements%eNoN(i) = 4 ! 4-node quadrangle.
         case(4)
            gmshElements%eNoN(i) = 4 ! 4-node tetrahedron.
         case(5)
            gmshElements%eNoN(i) = 8 ! 8-node hexahedron.
         case(8)
            gmshElements%eNoN(i) = 3 ! 3-node second order line
         case(9)
            gmshElements%eNoN(i) = 6 ! 6-node second order triangle
         case(10)
            gmshElements%eNoN(i) = 9 ! 9-node second order quadrangle
         case(11)
            gmshElements%eNoN(i) = 10 ! 10-node second order tetrahedron
         case(12)
            gmshElements%eNoN(i) = 27 ! 27-node second order hexahedron
         case(15)
            gmshElements%eNoN(i) = 1 ! 1-node point.
         case(16)
            gmshElements%eNoN(i) = 8 ! 8-node second order quadrangle
         case(17)
            gmshElements%eNoN(i) = 20 ! 20-node second order hexahedron
         case default
            write(stdout,ftab4) "ERROR: elelment type not defined!"
         end select

         do j = 1, gmshElements%numElementsInBlock(i)
            read(fid,*)
         end do
      end do

      maxeNoN = MAXVAL(gmshElements%eNoN)
      allocate(gmshElements%conn(maxeNoN,gmshElements%numElements))
      ie   = 0
      gmshElements%conn = 0
      rewind(fid)
      call findKwrd(fid, "$Elements")
      read(fid,*)
      do i = 1, gmshElements%numElementBlocks
         read(fid,*)
         eNoN = gmshElements%eNoN(i)
         do j = 1, gmshElements%numElementsInBlock(i)
            ie = ie+1
            read(fid,*) itemp, gmshElements%conn(1:eNoN,ie)
         end do
      end do

      return
      end subroutine gmsh_readelements
!***********************************************************************
!     Read Entities
      subroutine gmsh_readentities(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, maxNumTags
      integer :: numPhysicalTags, PhysicalTags, Tag
      real(kind=8) :: rtemp

      rewind(fid)
      call findKwrd(fid, "$Entities")
      read(fid,*) gmshEntities%numPoints, gmshEntities%numCurves, &
                  gmshEntities%numSurfaces, gmshEntities%numVolumes

      maxNumTags = 1 !gmshPhysicalNames%num

      ! Point entites
      if (gmshEntities%numPoints > 0) then
         allocate(gmshEntities%Point_Tags(gmshEntities%numPoints))
         allocate(gmshEntities%Point_numPhysicalTags(gmshEntities%numPoints))
         allocate(gmshEntities%Point_PhysicalTags(gmshEntities%numPoints,maxNumTags))
         gmshEntities%Point_Tags = 0
         gmshEntities%Point_numPhysicalTags = 0
         gmshEntities%Point_PhysicalTags = 0

         do i = 1, gmshEntities%numPoints
            read(fid,*) Tag, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each point entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Point_Tags(i) = Tag
               gmshEntities%Point_numPhysicalTags(i) = 1
               gmshEntities%Point_PhysicalTags(i,1) = PhysicalTags
            end if
         end do
      end if

      ! Curve entites
      if (gmshEntities%numCurves > 0) then
         allocate(gmshEntities%Curve_Tags(gmshEntities%numCurves))
         allocate(gmshEntities%Curve_numPhysicalTags(gmshEntities%numCurves))
         allocate(gmshEntities%Curve_PhysicalTags(gmshEntities%numCurves,maxNumTags))
         gmshEntities%Curve_Tags = 0
         gmshEntities%Curve_numPhysicalTags = 0
         gmshEntities%Curve_PhysicalTags = 0

         do i = 1, gmshEntities%numCurves
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each curve entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Curve_Tags(i) = Tag
               gmshEntities%Curve_numPhysicalTags(i) = 1
               gmshEntities%Curve_PhysicalTags(i,1) = PhysicalTags
            end if
         end do
      end if

      ! Surface entites
      if (gmshEntities%numSurfaces > 0) then
         allocate(gmshEntities%Surface_Tags(gmshEntities%numSurfaces))
         allocate(gmshEntities%Surface_numPhysicalTags(gmshEntities%numSurfaces))
         allocate(gmshEntities%Surface_PhysicalTags(gmshEntities%numSurfaces,maxNumTags))
         gmshEntities%Surface_Tags = 0
         gmshEntities%Surface_numPhysicalTags = 0
         gmshEntities%Surface_PhysicalTags = 0

         do i = 1, gmshEntities%numSurfaces
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each surface entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Surface_Tags(i) = Tag
               gmshEntities%Surface_numPhysicalTags(i) = 1
               gmshEntities%Surface_PhysicalTags(i,1) = PhysicalTags
            end if
         end do
      end if

      ! Volume entites
      if (gmshEntities%numVolumes > 0) then
         allocate(gmshEntities%Volume_Tags(gmshEntities%numVolumes))
         allocate(gmshEntities%Volume_numPhysicalTags(gmshEntities%numVolumes))
         allocate(gmshEntities%Volume_PhysicalTags(gmshEntities%numVolumes,maxNumTags))
         gmshEntities%Volume_Tags = 0
         gmshEntities%Volume_numPhysicalTags = 0
         gmshEntities%Volume_PhysicalTags = 0

         do i = 1, gmshEntities%numVolumes
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each volume entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Volume_Tags(i) = Tag
               gmshEntities%Volume_numPhysicalTags(i) = 1
               gmshEntities%Volume_PhysicalTags(i,1) = PhysicalTags
            end if
         end do
      end if

      return
      end subroutine gmsh_readentities
!***********************************************************************
!     Assemble body mesh
      subroutine gmsh_bodymesh(lM)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: lM

      integer :: i, j, id, physicalTag, numEntityTags
      integer :: Tag, iTag, nEl
      integer,allocatable :: entityTag(:), gIEN(:,:)
      character(len=80) :: name

      j = 0
      physicalTag = -1
      do i = 1, gmshPhysicalNames%num
         if (lM%nsd .eq. gmshPhysicalNames%dimension(i)) then
            id = i
            physicalTag = gmshPhysicalNames%physicalTag(i)
            j = j+1
            name = gmshPhysicalNames%name(i)
         end if
      end do
      if (j .gt. 1) then
         write(stdout,ftab4) "ERROR: can only handle one body mesh."
         stop
      end if
      if (physicalTag .eq. -1) then
         write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
         stop
      end if

      ! Using the physical tag to determine the
      ! surface tag (entity tag) for 2D problem or
      ! volume tag (entity tag) for 3D problem.
      ! Please note that there might be multiple entity tags
      ! for one physical tag.
      if (lM%nsd .eq. 2) then
         allocate(entityTag(gmshEntities%numSurfaces))
         entityTag = 0
         numEntityTags = 0
         do i = 1, gmshEntities%numSurfaces
            do j = 1, gmshEntities%Surface_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Surface_PhysicalTags(i,j)) then
                  numEntityTags = numEntityTags + 1
                  entityTag(numEntityTags) = gmshEntities%Surface_Tags(i)
               end if
            end do
         end do
      elseif (lM%nsd .eq. 3) then
         allocate(entityTag(gmshEntities%numVolumes))
         entityTag = 0
         numEntityTags = 0
         do i = 1, gmshEntities%numVolumes
            do j = 1, gmshEntities%Volume_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Volume_PhysicalTags(i,j)) then
                  numEntityTags = numEntityTags + 1
                  entityTag(numEntityTags) = gmshEntities%Volume_Tags(i)
               end if
            end do
         end do
      end if
      if (numEntityTags .eq. 0) then
         write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
         stop
      end if

      ! raw global connectivity
      nEl = 0
      do iTag = 1, numEntityTags
         Tag = entityTag(iTag)
         j = 0
         do i = 1, gmshElements%numElementBlocks
            if (gmshElements%ElementDim(i) .eq. lM%nsd) then
               if (gmshElements%EntityTag(i) .eq. Tag) then
                  lM%eNoN = gmshElements%eNoN(i)
                  if (.not. allocated(gIEN)) &
                     allocate(gIEN(lM%eNoN,gmshElements%numElements))
                  gIEN(1:lM%eNoN,nEl+1:nEl + gmshElements%numElementsinBlock(i)) = &
                  gmshElements%conn(1:lM%eNoN,j+1:j+gmshElements%numElementsinBlock(i))
                  nEl = nEl + gmshElements%numElementsinBlock(i)
                  exit
               end if
            end if
            j = j + gmshElements%numElementsInBlock(i)
         end do
      end do

      ! Because not all points are mesh points, e.g. center of the circule,
      ! we need to build global-to-local mapping to eliminate those points.
      allocate(g2l(gmshNodes%numNodes))
      g2l = 0
      id = 0
      do i = 1, nEl
         do j = 1, lM%eNoN
            Tag = gIEN(j,i)
            if (g2l(Tag) .eq. 0) then
               id = id + 1
               g2l(Tag) = id
            end if
         end do
      end do

      ! connectivity
      lM%nEl = nEl
      allocate(lM%IEN(lM%eNoN,lM%nEl))
      do i = 1, nEl
         do j = 1, lM%eNoN
            lM%IEN(j,i) = g2l(gIEN(j,i))
         end do
      end do

      ! coordinates
      lM%nNo = id
      allocate(lM%x(lM%nsd,lM%nNo))
      do i = 1, gmshNodes%numNodes
         if (g2l(i) .ne. 0) then
            lM%x(:,g2l(i)) = gmshNodes%coord(1:lM%nsd,i)
         end if
      end do

      ! deallocate
      deallocate(entityTag, gIEN)

      return
      end subroutine gmsh_bodymesh
!***********************************************************************
!     Assemble boundary meshes
      subroutine gmsh_boundarymesh(lM)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: lM

      integer :: i, j, k, l, ii, ie, Tag, iFa, iTag
      integer :: nEl, eNoN, numEntityTags, numinBlocks
      integer,allocatable :: physicalTag(:), gIEN(:,:), lg2l(:)
      integer,allocatable :: sharedelem(:,:), numshared(:), entityTag(:)
      logical,allocatable :: flag(:)

      ! For boundary points, it requires two global-to-local mappings.
      ! g2l: since we need to remove some unnecessary points in gmshNodes,
      !      this mapping is built to contruct the body mesh
      ! lg2l: when writing boundaries to vtp files, we need a local connectivity
      !      system, this mapping maps **body mesh**, i.e. after g2l, to local.
      allocate(lg2l(lM%nNo))

      ! element that share the same node
      ! this is for finding global element id
      allocate(sharedelem(lM%nNo,50), numshared(lM%nNo))
      sharedelem = 0
      numshared  = 0
      do i = 1, lM%nEl
         do j = 1, lM%eNoN
            k = lM%IEN(j,i)
            numshared(k) = numshared(k)+1
            if (numshared(k) .GT. 50) then
               write(stdout,ftab4) "ERROR: numshared is larger than 50."
               stop
            end if
            sharedelem(k,numshared(k)) = i
         end do
      end do

      lM%nFa = 0
      do i = 1, gmshPhysicalNames%num
         if (gmshPhysicalNames%dimension(i) .eq. lM%nsd-1) then
            lM%nFa = lM%nFa+1
         end if
      end do
      allocate(lM%fa(lM%nFa), physicalTag(lM%nFa))
      j = 0
      do i = 1, gmshPhysicalNames%num
         if (gmshPhysicalNames%dimension(i) .eq. lM%nsd-1) then
            j = j + 1
            physicalTag(j) = gmshPhysicalNames%physicalTag(i)
            lM%fa(j)%fname = gmshPhysicalNames%name(i)
         end if
      end do

      ! Loop over faces to construct mesh files.
      ! Please note that there might be multiple entity tags
      ! for one physical tag
      do iFa = 1, lM%nFa

         ! Find the corresponding EntityTags for boundary meshes
         ! ndim=2, CurveTag
         ! ndim=3, SurfaceTag
         if (lM%nsd .eq. 2) then
            if (.not. allocated(entityTag)) &
               allocate(entityTag(gmshEntities%numCurves))
            entityTag = 0
            numEntityTags = 0
            do j = 1, gmshEntities%numCurves
               do k = 1, gmshEntities%Curve_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Curve_physicalTags(j,k)) then
                     numEntityTags = numEntityTags + 1
                     entityTag(numEntityTags) = gmshEntities%Curve_Tags(j)
                  end if
               end do
            end do
         elseif (lM%nsd .eq. 3) then
            if (.not. allocated(entityTag)) &
               allocate(entityTag(gmshEntities%numSurfaces))
            entityTag = 0
            numEntityTags = 0
            do j = 1, gmshEntities%numSurfaces
               do k = 1, gmshEntities%Surface_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Surface_physicalTags(j,k)) then
                     numEntityTags = numEntityTags + 1
                     entityTag(numEntityTags) = gmshEntities%Surface_Tags(j)
                  end if
               end do
            end do
         end if
         if (numEntityTags .eq. 0) then
            write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
            stop
         end if

         ! raw connectivity for boundary mesh
         nEl = 0
         do iTag = 1, numEntityTags
            Tag = entityTag(iTag)
            j = 0
            do i = 1, gmshElements%numElementBlocks
               if (gmshElements%ElementDim(i) .eq. lM%nsd-1) then
               if (gmshElements%EntityTag(i) .eq. Tag) then
                  eNoN = gmshElements%eNoN(i)
                  ! lM%fa(iFa)%eNoN = eNoN
                  if (.not. allocated(gIEN)) &
                     allocate(gIEN(eNoN,gmshElements%numElements))
                  numinBlocks = gmshElements%numElementsinBlock(i)
                  gIEN(1:eNoN,nEl+1:nEl+numinBlocks) = gmshElements%conn(1:eNoN,j+1:j+numinBlocks)
                  nEl = nEl + numinBlocks
                  exit
               end if
               end if
               j = j + gmshElements%numElementsInBlock(i)
            end do
         end do

         ! first global to local mapping
         do i = 1, nEl
            do j = 1, eNoN
               gIEN(j,i) = g2l(gIEN(j,i))
            end do
         end do

         ! determine gN
         ! second global to local mapping
         lg2l = 0
         lM%fa(iFa)%nNo  = 0
         lM%fa(iFa)%eNoN = eNoN
         lM%fa(iFa)%nEl  = nEl
         lM%fa(iFa)%nsd  = lM%nsd
         do i = 1, lM%fa(iFa)%nEl
            do j = 1, lM%fa(iFa)%eNoN
               k = gIEN(j,i)
               if (lg2l(k) .eq. 0) then
                  lM%fa(iFa)%nNo = lM%fa(iFa)%nNo + 1
                  lg2l(k) = lM%fa(iFa)%nNo
               end if
            end do
         end do
         allocate(lM%fa(iFa)%gN(lM%fa(iFa)%nNo),lM%fa(iFa)%x(lM%nsd,lM%fa(iFa)%nNo))
         do i = 1, lM%nNo
            if (lg2l(i) .ne. 0) then
               j = lg2l(i)
               lM%fa(iFa)%gN(j) = i
               lM%fa(iFa)%x(:,j) = lM%x(:,i)
            end if
         end do
         ! DO NOT convert gIEN to local IEN yet
         allocate(lM%fa(iFa)%IEN(lM%fa(iFa)%eNoN,lM%fa(iFa)%nEl))
         lM%fa(iFa)%IEN = gIEN(:,1:lM%fa(iFa)%nEl)

         ! find global element id
         allocate(lM%fa(iFa)%gE(lM%fa(iFa)%nEl))
         if (allocated(flag)) deallocate(flag)
         allocate(flag(lM%fa(iFa)%eNoN))
         do i = 1, lM%fa(iFa)%nEl
            ii = gIEN(1,i)
            do k = 1, numshared(ii)
               flag = .FALSE.
               ie = sharedelem(ii,k)
               do j = 1, lM%fa(iFa)%eNoN
                  do l = 1, lM%eNoN
                     if (gIEN(j,i) .eq. lM%IEN(l,ie)) then
                        flag(j) = .TRUE.
                        exit
                     end if
                  end do
               end do
               if (ALL(flag)) then
                  lM%fa(iFa)%gE(i) = ie
                  exit
               end if
            end do
         end do

      end do !iFa

      ! deallocate
      deallocate(physicalTag, gIEN, g2l, lg2l)
      deallocate(sharedelem, numshared, entityTag, flag)

      return
      end subroutine gmsh_boundarymesh
!***********************************************************************
