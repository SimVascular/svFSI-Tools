!--------------------------------------------------------------------
!
! This program is a utility to convert Gambit Neutral Mesh Files to
! VTK based format (VTU: Unstructured Grid / VTP: Polydata)
! Optionally, allows to convert to higher order elements
! 2D: TRI3 -> TRI6  ; QUAD4 -> QUAD8 ; QUAD4 -> QUAD9
! 3D: TET4 -> TET10 ; HEX8  -> HEX20 ; HEX8  -> HEX27
!
!--------------------------------------------------------------------

      program convert_mesh
      use mesh
      use genUtils
      implicit none

      type(meshType) :: msh
      character(len=strL) :: fName

      integer i

      i = IARGC()
      if (i .eq. 0) then
         write(stdout,ftab4) "ERROR: Input file name not specified"
         STOP
      else if (i .gt. 1) then
         write(stdout,ftab4) "ERROR: Too many arguments"
         STOP
      end if
      call getarg(1,fName)

      write(stdout,ftab1) repeat('=', 48)
      if ( endsWith(trim(fName), '.neu') ) then
         write(stdout, ftab1) "Reading mesh from Gambit neu file"// &
            "   <----   "//trim(fName)
         msh%fname = fName
         call readGambitNeu(msh)

      else
         write(stdout,ftab1) "Reading input file  <----  "//trim(fName)
         call readInputFile(fName, msh)

         call readVTKMesh(msh)

      end if

      call refineMeshOrder(msh)

!     Now, write VTU/VTP files
      write(stdout, ftab1) "Writing to vtk file system.."
      call VTK(msh)
      write(stdout,ftab1) repeat('=', 48)

      end program convert_mesh

!***********************************************************************
