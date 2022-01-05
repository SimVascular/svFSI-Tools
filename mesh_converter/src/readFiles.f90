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
