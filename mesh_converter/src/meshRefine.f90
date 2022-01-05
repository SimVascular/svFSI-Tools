!***********************************************************************

      subroutine refineMeshOrder(lM)
      use varmod
      implicit none
      type(meshType), intent(inout) :: lM

      logical ecFlag, fcFlag, ccFlag
      integer nec, nfc, nfb, ncc, iFa
      character(len=strL) :: stmp

      integer, allocatable :: iec(:), ifc(:)
      real(kind=8), allocatable :: xec(:,:), xfc(:,:), xcc(:,:)

      type(hoeType)  :: hoe
      type(meshType) :: gM

      write(stdout,ftab1,advance='no') &
         "Convert mesh to higher order? (y/n) "
      read(*,'(A)') stmp
      stmp = to_lower(adjustl(trim(stmp)))
      if (stmp.eq.'0' .or. stmp.eq.'n' .or. stmp.eq.'no') &
         return

      if (lM%eType .eq. eType_TRI3 .or. &
          lM%eType .eq. eType_TET4 .or. &
          lM%eType .eq. eType_QUD4) then
         call checkIEN(lM)
      end if

      call selectEle(lM%eType, hoe)

      ecFlag  = .false.
      fcFlag  = .false.
      ccFlag  = .false.
      if (allocated(hoe%eord)) ecFlag = .true.
      if (allocated(hoe%ford)) fcFlag = .true.
      if (hoe%eType .eq. eType_QUD9 .or. &
          hoe%eType .eq. eType_HEX27) ccFlag = .true.

      nec   = 0
      nfc   = 0
      nfb   = 0
      ncc   = 0

!     Get edge connectivity and use it to compute all edge nodes
      if (ecFlag) then
         call getEdgeConnectivity(lM)
         nec = lM%econ%nnz/2
         allocate(iec(lM%econ%nnz), xec(lM%nsd,nec))
         call calcEdgeCenters(nec, iec, xec)
         write(stdout,ftab2) "Num. edge centers: "//trim(str(nec))
      end if

!     Get face adjacency and compute face centers. Compute cell
!     centers depending on the element type chosen
      if (fcFlag) then
         call getFaceConnectivity(lM)
         nfb = SUM(lM%fa(:)%nEl)
         nfc = (lM%fcon%nnz - nfb)/2 + nfb
         allocate(ifc(lM%fcon%nnz), xfc(lM%nsd,nfc))
         call calcFaceCenters(nfc, ifc, xfc)
         write(stdout,ftab2) "Num. face centers: "//trim(str(nfc))
      end if

      if (ccFlag) then
         ncc = lM%nEl
         allocate(xcc(lM%nsd,ncc))
         call calcCellCenters(ncc, xcc)
         write(stdout,ftab2) "Num. cell centers: "//trim(str(ncc))
      end if

!     Contruct the new mesh
      call updateMesh()

!     Construct face structure
      do iFa=1, lM%nFa
         call updateFace(gM%fa(iFa), lM%fa(iFa))
      end do

!     Destroy lM
      call destroy(lM)

!     Deep copy gM --> lM
      call deepCopy(gM, lM)

!     Destroy gM
      call destroy(gM)

      if (allocated(xec)) deallocate(iec, xec)
      if (allocated(xfc)) deallocate(ifc, xfc)
      if (allocated(xcc)) deallocate(xcc)

      return
      contains
      !=============================================================
         subroutine calcEdgeCenters(ne, eid, xe)
         implicit none
         integer, intent(in) :: ne
         integer, intent(out) :: eid(lM%econ%nnz)
         real(kind=8), intent(out) :: xe(lM%nsd,ne)

         integer a, b, i, j

         j   = 0
         eid = 0
         xe  = 0d0
         do a=1, lM%nNo
            do i=lM%econ%prow(a), lM%econ%prow(a+1)-1
               b = lM%econ%pcol(i)
               if (a .gt. b) cycle
               j = j + 1
               if (j .gt. ne) then
                  write(stdout,ftab4) "ERROR: out of bounds while "// &
                     "computing edge centers"
                  stop
               end if
               eid(i) = j
               xe(:,j) = 0.5d0*(lM%x(:,a) + lM%x(:,b))
            end do
         end do

         do a=1, lM%nNo
            do i=lM%econ%prow(a), lM%econ%prow(a+1)-1
               b = lM%econ%pcol(i)
               if (a .lt. b) cycle
               j = getEconColID(b, a)
               eid(i) = eid(j)
            end do
         end do

         return
         end subroutine calcEdgeCenters
      !=============================================================
         subroutine calcFaceCenters(nf, fid, xf)
         implicit none
         integer, intent(in) :: nf
         integer, intent(out) :: fid(lM%fcon%nnz)
         real(kind=8), intent(out) :: xf(lM%nsd,nf)

         integer :: a, e, f, i, j, n, Ac

         n   = 0
         xf  = 0.0d0
         fid = 0
         do e=1, lM%nEl
            do i=1, lM%eNoF
               j = i + lM%fcon%prow(e) - 1
               f = lM%fcon%pcol(j)
               if (e.gt.f .and. f.ne.0) cycle

               n = n + 1
               fid(j) = n
               do a=1, lM%fNoN
                  Ac = lM%IEN(lM%ford(a,i),e)
                  xf(:,n) = xf(:,n) + lM%x(:,Ac)
               end do
               xf(:,n) = xf(:,n) / real(lM%fNoN,kind=8)
            end do
         end do

         do e=1, lM%nEl
            do i=1, lM%eNoF
               j = i + lM%fcon%prow(e) - 1
               f = lM%fcon%pcol(j)
               if (e.lt.f .or. f.eq.0) cycle
               fid(j) = fid(getFconColID(f, e))
            end do
         end do

         return
         end subroutine calcFaceCenters
      !=============================================================
         subroutine calcCellCenters(nc, xc)
         implicit none
         integer, intent(in) :: nc
         real(kind=8), intent(out) :: xc(lM%nsd,nc)

         integer a, e, n

         n  = 0
         xc = 0.0d0
         do e=1, lM%nEl
            n = n + 1
            do a=1, lM%eNoN
               xc(:,n) = xc(:,n) + lM%x(:,lM%IEN(a,e))
            end do
            xc(:,n) = xc(:,n) / real(lM%eNoN, kind=8)
         end do

         return
         end subroutine calcCellCenters
      !=============================================================
         function getEconColID(r, c) result (indx)
         implicit none
         integer :: r, c, indx
         integer :: i

         indx = 0
         do i=lM%econ%prow(r), lM%econ%prow(r+1)-1
            if (c .eq. lM%econ%pcol(i)) then
               indx = i
               exit
            end if
         end do

         if (indx .eq. 0) then
            write(stdout,ftab4) "ERROR: could not find a col. id "// &
               "nadj: "//trim(str(r))//", "//trim(str(c))
            stop
         end if

         return
         end function getEconColID
      !=============================================================
         function getFconColID(r, c) result (indx)
         implicit none
         integer :: r, c, indx
         integer :: i

         indx = 0
         do i=lM%fcon%prow(r), lM%fcon%prow(r+1)-1
            if (c .eq. lM%fcon%pcol(i)) then
               indx = i
               exit
            end if
         end do

         if (indx .eq. 0) then
            write(stdout,ftab4) "ERROR: could not find a col. id "// &
               "eadj: "//trim(str(r))//", "//trim(str(c))
            stop
         end if

         return
         end function getFconColID
      !=============================================================
         subroutine updateMesh()
         implicit none

         integer a, b, e, i, j, Ac, nshft

         gM%fname = lM%fname
         gM%nsd   = lM%nsd
         gM%eNoN  = hoe%eNoN
         gM%nNo   = lM%nNo + nec + nfc + ncc
         gM%nEl   = lM%nEl
         gM%nFa   = lM%nFa
         allocate(gM%fa(gM%nFa))
         call selectEle(gM)

!        Build mesh connectivity
         allocate(gM%IEN(gM%eNoN,gM%nEl), gM%x(gM%nsd,gM%nNo))
         gM%x   = 0d0
         gM%IEN = 0

!        Copy old mesh coordinates and connectivity
         do a=1, lM%nNo
            gM%x(:,a) = lM%x(:,a)
         end do

         do e=1, gM%nEl
            do a=1, lM%eNoN
               gM%IEN(a,e) = lM%IEN(a,e)
            end do
         end do

!        Add edge centers
         nshft = lM%nNo
         if (ecFlag) then
            do e=1, gM%nEl
               do i=1, gM%eNoE
                  a  = lM%IEN(gM%eord(1,i),e)
                  b  = lM%IEN(gM%eord(2,i),e)
                  j  = iec(getEconColID(a,b))
                  Ac = j + nshft
                  gM%x(:,Ac) = xec(:,j)
                  gM%IEN(lM%eNoN+i,e) = Ac
               end do
            end do
         end if

!        Add face centers and cell centers
         nshft = lM%nNo + nec
         if (fcFlag) then
            do i=1, nfc
               a = nshft + i
               gM%x(:,a) = xfc(:,i)
            end do

            do e=1, gM%nEl
               do i=1, gM%eNoF
                  a = lM%eNoN + gM%eNoE + i
                  b = i + lM%fcon%prow(e) - 1
                  gM%IEN(a,e) = ifc(b) + nshft
               end do
            end do
         end if

         nshft = lM%nNo + nec + nfc
         if (ccFlag) then
            do i=1, ncc
               a = nshft + i
               gM%x(:,a) = xcc(:,i)
            end do

            do e=1, gM%nEl
               gM%IEN(gM%eNoN,e) = e + nshft
            end do
         end if

!         write(1001,'(a)') repeat('=',48)
!         write(1001,'(a)') "Mesh: "//trim(gM%fname)
!         write(1001,'(a)') "gx: "
!         do i=1, gM%nNo
!            write(1001,'(4x,a)',advance='no') trim(str(i))//": "
!            do j=1, gM%nsd
!               write(1001,'(a)',advance='no') " "//trim(str(gM%x(j,i)))
!            end do
!            write(1001,'(a)')
!         end do
!         call flush(1001)

!         write(1001,'(a)') repeat('=',48)
!         write(1001,'(a)') "IEN: "
!         do i=1, gM%nEl
!            write(1001,'(4x,a)',advance='no') trim(str(i))//": "
!            do j=1, gM%eNoN
!               write(1001,'(a)',advance='no') " "//trim(str(gM%IEN(j,i)))
!            end do
!            write(1001,'(a)')
!         end do
!         call flush(1001)

         return
         end subroutine updateMesh
      !=============================================================
         subroutine updateFace(gFa, lFa)
         implicit none
         type(faceType), intent(in) :: lFa
         type(faceType), intent(inout) :: gFa

         integer a, b, e, i, j, Ac, nshft
         integer, allocatable :: ptr(:), incE(:)

         gFa%fname = lFa%fname
         gFa%nsd   = lFa%nsd
         gFa%nEl   = lFa%nEl

         gFa%nNo   = lFa%nNo

         if (ecFlag) then
            allocate(incE(nec))
            incE = 0
            do e=1, lFa%nEl
               do i=1, lFa%eNoE
                  a = lFa%IEN(lFa%eord(1,i),e)
                  b = lFa%IEN(lFa%eord(2,i),e)
                  j = getEconColID(a, b)
                  incE(iec(j)) = 1
               end do
            end do
            gFa%nNo = gFa%nNo + sum(incE(:))
         end if
         if (fcFlag) gFa%nNo = gFa%nNo + lFa%nEl

!        Build face connectivity
         allocate(gFa%x(gFa%nsd,gFa%nNo), gFa%gN(gFa%nNo), &
            gFa%gE(gFa%nEl), gFa%IEN(gFa%eNoN,gFa%nEl))
         gFa%x   = 0d0
         gFa%gN  = 0
         gFa%gE  = lFa%gE
         gFa%IEN = 0

!        Copy old mesh coordinates, connectivity and gN
         do a=1, lFa%nNo
            gFa%x(:,a) = lFa%x(:,a)
            gFa%gN(a)  = lFa%gN(a)
         end do

         do e=1, lFa%nEl
            do a=1, lFa%eNoN
               gFa%IEN(a,e) = lFa%IEN(a,e)
            end do
         end do

!        Add edge centers
         j     = lFa%nNo
         nshft = lM%nNo
         if (ecFlag) then
            do e=1, lFa%nEl
               do i=1, lFa%eNoE
                  a  = lFa%IEN(lFa%eord(1,i),e)
                  b  = lFa%IEN(lFa%eord(2,i),e)
                  Ac = iec(getEconColID(a, b))
                  gFa%IEN(lFa%eNoN+i,e) = Ac + nshft
                  if (incE(Ac) .eq. 1) then
                     incE(Ac) = 0
                     j  = j + 1
                     Ac = Ac + nshft
                     gFa%gN(j)  = Ac
                     gFa%x(:,j) = gM%x(:,Ac)
                  end if
               end do
            end do
         end if

!        Add face centers and cell centers
         nshft = lM%nNo + nec
         if (fcFlag) then
            allocate(ptr(lFa%eNoN))
            do e=1, lFa%nEl
               ptr = lFa%IEN(:,e)
               call sort(ptr)

               b  = getMatchingFaceIndx(gFa%gE(e), lFa%eNoN, ptr)
               j  = j + 1
               Ac = ifc(b) + nshft
               gFa%gN(j)  = Ac
               gFa%x(:,j) = gM%x(:,Ac)
               gFa%IEN(gFa%eNoN,e) = Ac
            end do
         end if

!         write(1001,'(a)') "x: "
!         do i=1, gFa%nNo
!            write(1001,'(4x,a)',advance='no') trim(str(i))//"    "// &
!               trim(str(gFa%gN(i)))//"   "
!            do j=1, gFa%nsd
!               write(1001,'(a)',advance='no') " "//trim(str(gFa%x(j,i)))
!            end do
!            write(1001,'(a)')
!         end do
!         call flush(1001)

!         write(1001,'(a)') repeat('-',24)
!         write(1001,'(a)') "IEN: "
!         do i=1, gFa%nEl
!            write(1001,'(4x,a)',advance='no') trim(str(i))//"    "// &
!               trim(str(gFa%gE(i)))//"   "
!            do j=1, gFa%eNoN
!               write(1001,'(a)',advance='no') " "//trim(str(gFa%IEN(j,i)))
!            end do
!            write(1001,'(a)')
!         end do
!         call flush(1001)

         return
         end subroutine updateFace
      !=============================================================
         function getMatchingFaceIndx(e, fNoN, ptr) result(indx)
         implicit none
         integer, intent(in) :: e, fNoN, ptr(fNoN)
         integer :: indx

         integer a, i, j, k, l, inds(fNoN)

         if (fNoN .ne. lM%fNoN) then
            write(stdout,ftab4) "ERROR: inconsistent fNoN to find "//&
               "matching face index"
            stop
            stop
         end if

         indx = 0
         do i=1, lM%eNoF
            j = i + lM%fcon%prow(e) - 1
            if (lM%fcon%pcol(j) .ne. 0) cycle
            do a=1, lM%fNoN
               inds(a) = lM%IEN(lM%ford(a,i),e)
            end do
            call sort(inds)
            l = 0
            do k=1, lM%fNoN
               if (ptr(k) .eq. inds(k)) l = l + 1
            end do
            if (l .eq. lM%fNoN) then
               indx = j
               exit
            end if
         end do

         if (indx .eq. 0) then
            write(stdout,ftab4) "ERROR: couldn't find a matching "// &
               "face indx"
            stop
         end if

         return
         end function getMatchingFaceIndx
      !=============================================================
      end subroutine refineMeshOrder

!***********************************************************************
