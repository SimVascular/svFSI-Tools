!--------------------------------------------------------------------
!
!     Here data (mesh, solution) input/output is handled interfacing
!     with fortran-based VTK module.
!
!--------------------------------------------------------------------

!***********************************************************************

      subroutine VTK(lM)
      use mesh
      use vtkXMLMod
      implicit none

      type(meshType), intent(inout) :: lM
      integer :: iFa
      character(len=strL) :: fName
      logical :: flag

!     Write mesh vtu file
      write(fName,'(A)') trim(lM%fname)//".vtu"
      lM%IEN(:,:) = lM%IEN(:,:) - 1
      call writeVTU(lM, fName)
      lM%IEN(:,:) = lM%IEN(:,:) + 1

!     Remap face IEN before writing to vtk
      do iFa=1, lM%nFa
         call remapFaceIEN(lM%fa(iFa))
      end do

!     Write face vtp files
      inquire(file="mesh-surfaces",exist=flag)
      if (.not.flag) call system("mkdir  mesh-surfaces")
      do iFa=1, lM%nFa
         write(fName,'(A)') "mesh-surfaces/"//TRIM(lM%fa(iFa)%fname)// &
            ".vtp"
         lM%fa(iFa)%IEN = lM%fa(iFa)%IEN - 1
         call writeVTP(lM%fa(iFa), fName)
         lM%fa(iFa)%IEN = lM%fa(iFa)%IEN + 1
      end do

      return
      contains
      !=============================================================
         subroutine remapFaceIEN(lFa)
         implicit none
         type(faceType), intent(inout) :: lFa

         integer a, e, Ac, ptr(lM%nNo)

         ptr = 0
         do a=1, lFa%nNo
            Ac = lFa%gN(a)
            ptr(Ac) = a
         end do

         do e=1, lFa%nEl
            do a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               Ac = ptr(Ac)
               IF (Ac .eq. 0) then
                  write(stdout,ftab4) "Error in remaping face IEN"
                  stop
               end if
               lFa%IEN(a,e) = Ac
            end do
         end do

         return
         end subroutine remapFaceIEN
      !=============================================================
      end subroutine VTK

!**********************************************************************

      subroutine writeVTU(lM, fName)
      use mesh
      use vtkXMLMod

      implicit none

      type(meshType), intent(in) :: lM
      character(len=strL), intent(in) :: fName

      type(vtkXMLtype) :: vtu
      integer :: iStat

      call vtkInitWriter(vtu, trim(fName), iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (init)"
         stop
      end if

      call putVTK_pointCoords(vtu, lM%x, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (coords)"
         stop
      end if

      call putVTK_elemIEN(vtu, lM%IEN, lM%vtktype, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (ien)"
         stop
      end if

      call vtkWriteToFile(vtu, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error"
         stop
      end if

      call flushVTK(vtu)

      return
      end subroutine writeVTU

!**********************************************************************

      subroutine writeVTP(lFa, fName)
      use mesh
      use vtkXMLMod
      implicit none

      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fName

      type(vtkXMLtype) :: vtp
      integer :: iStat

      call vtkInitWriter(vtp, trim(fName), iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (init)"
         stop
      end if

      call putVTK_pointCoords(vtp, lFa%x, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (coords)"
         stop
      end if

      call putVTK_elemIEN(vtp, lFa%IEN, lFa%vtktype, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (ien)"
         stop
      end if

      call putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (point data)"
         stop
      end if

      call putVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (cell data)"
         stop
      end if

      call vtkWriteToFile(vtp, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error"
         stop
      end if

      call flushVTK(vtp)

      return
      end subroutine writeVTP

!***********************************************************************

      subroutine readVTU(lM, fName)
      use mesh
      use vtkXMLMod
      implicit none

      type(meshType), intent(inout) :: lM
      character(len=strL) :: fName

      type(vtkXMLtype) :: vtu
      integer :: iStat
      real(kind=8), allocatable, dimension(:,:) :: tmpX

      iStat = 0
      write(stdout,ftab2) " <VTK XML Parser> Loading file <"//trim(fName)//">"
      call loadVTK(vtu, fName, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (init)"
         stop
      end if

      call getVTK_numPoints(vtu, lM%nNo, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (num points)"
         stop
      end if

      call getVTK_numElems(vtu, lM%nEl, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (num cells)"
         stop
      end if

      call getVTK_nodesPerElem(vtu, lM%eNoN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (nodes per cell)"
         stop
      end if

      allocate(lM%x(lM%nsd,lM%nNo),tmpX(maxNSD,lM%nNo))
      call getVTK_pointCoords(vtu, tmpX, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (coords)"
         stop
      end if
      lM%x(:,:) = tmpX(1:lM%nsd,:)
      deallocate(tmpX)

      allocate(lM%IEN(lM%eNoN,lM%nEl))
      call getVTK_elemIEN(vtu, lM%IEN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file read error (ien)"
         stop
      end if
      lM%IEN = lM%IEN + 1

      call flushVTK(vtu)

      return
      end subroutine readVTU

!***********************************************************************

      subroutine readVTP(lFa, fName)
      use mesh
      use vtkXMLMod
      implicit none

      type(faceType), intent(inout) :: lFa
      character(len=strL) :: fName

      type(vtkXMLtype) :: vtp
      integer :: iStat, a, e, Ac

      real(kind=8), allocatable :: tmpX(:,:)

      iStat = 0
      write(stdout,ftab2) " <VTK XML Parser> Loading file <"//trim(fName)//">"
      call loadVTK(vtp, fName, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (init)"
         stop
      end if

      call getVTK_numPoints(vtp, lFa%nNo, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (num points)"
         stop
      end if

      call getVTK_numElems(vtp, lFa%nEl, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (num cells)"
         stop
      end if

      call getVTK_nodesPerElem(vtp, lFa%eNoN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (nodes per cell)"
         stop
      end if

      allocate(lFa%x(lFa%nsd,lFa%nNo), tmpX(maxNSD,lFa%nNo))
      call getVTK_pointCoords(vtp, tmpX, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (coords)"
         stop
      end if
      lFa%x(:,:) = tmpX(1:lFa%nsd,:)
      deallocate(tmpX)

      allocate(lFa%IEN(lFa%eNoN,lFa%nEl))
      call getVTK_elemIEN(vtp, lFa%IEN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (ien)"
         stop
      end if

      allocate(lFa%gN(lFa%nNo))
      call getVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      if (iStat .lt. 0) then
         deallocate(lFa%gN)
      else
         do e=1, lFa%nEl
            do a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)+1
               Ac = lFa%gN(Ac)
               lFa%IEN(a,e) = Ac
            end do
         end do
      end if

      allocate(lFa%gE(lFa%nEl))
      call getVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file read error (ge)"
         stop
      end if

      call flushVTK(vtp)

      return
      end subroutine readVTP

!***********************************************************************
