
!***********************************************************************

      module meshUtils
      use stdParams
      use genUtils
      use mesh
      implicit none

      private

      integer maxnn, maxnne
      integer, allocatable :: xnn(:), xnnL(:,:)
      integer, allocatable :: xnne(:), xneL(:,:)

      public :: checkIEN
      public :: checkIENB

      public :: getEdgeConnectivity
      public :: getFaceConnectivity

      contains
         !=============================================================
!        Checks IEN array of mesh depending on the element type. Resets
!        ordering of nodes if any inconsistency is found
         subroutine checkIEN(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         integer a, b, e, i, sn(4)
         real(kind=8), allocatable :: v(:,:), xl(:,:)

         allocate(v(lM%nsd,lM%eNoN), xl(lM%nsd,lM%eNoN))

         v  = 0.0d0
         xl = 0.0d0
         select case (lM%eType)
         case (eType_TRI3)
            do e=1, lM%nEl
               a = 1; b = 1
               xl     = lM%x(:,lM%IEN(:,e))
               v(:,1) = xl(:,2) - xl(:,1)
               v(:,2) = xl(:,3) - xl(:,2)
               sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))

               if (sn(1) .eq. -1) then
                  a = 1; b = 2
               else if (sn(1) .eq. 0) then
                  write(stdout,ftab4) "ERROR: triangular element "// &
                     trim(str(e))//" is distorted"
                  stop
               end if
               call swap(lM%IEN(a,e), lM%IEN(b,e))
            end do

         case (eType_TET4)
            do e=1, lM%nEl
               a = 1; b = 1
               xl     = lM%x(:,lM%IEN(:,e))
               v(:,1) = xl(:,2) - xl(:,1)
               v(:,2) = xl(:,3) - xl(:,2)
               v(:,3) = xl(:,4) - xl(:,3)
               v(:,4) = CROSS(v(:,1:2))
               sn(1)  = SGN(SUM(v(:,3)*v(:,4)))

               if (sn(1) .eq. 1) then
                  a = 1; b = 2
               else if (sn(1) .eq. 0) then
                  write(stdout,ftab4) "Element "//trim(STR(e))// &
                     " is distorted"
                  stop
               end if
               call swap(lM%IEN(a,e), lM%IEN(b,e))
            end do

         case (eType_QUD4)
            do e=1, lM%nEl
               a = 1; b = 1
               xl     = lM%x(:,lM%IEN(1:4,e))
               v(:,1) = xl(:,2) - xl(:,1)
               v(:,2) = xl(:,3) - xl(:,2)
               v(:,3) = xl(:,4) - xl(:,3)
               v(:,4) = xl(:,1) - xl(:,4)
               sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))
               sn(2)  = SGN(v(1,2)*v(2,3) - v(2,2)*v(1,3))
               sn(3)  = SGN(v(1,3)*v(2,4) - v(2,3)*v(1,4))
               sn(4)  = SGN(v(1,4)*v(2,1) - v(2,4)*v(1,1))
               i = sn(1) + sn(2) + sn(3) + sn(4)

               if (i .eq. 0) THEN
                  if (sn(1) .eq. sn(2)) then
                     if (sn(1) .eq. 1) then
                        a = 1; b = 4
                     else
                        a = 2; b = 3
                     end if
                  else
                     if (sn(1) .eq. 1) then
                        a = 3; b = 4
                     else
                        a = 1; b = 2
                     end if
                  end if

               else if (i .eq. -4) then
                  a = 1; b = 3

               else if (i.eq.2 .or. any(sn.eq.0)) then
                  write(stdout,ftab4) "Element "//trim(STR(e))// &
                     " is distorted"
                  stop
               end if
               call swap(lM%IEN(a,e), lM%IEN(b,e))
            end do

         end select

         deallocate(v, xl)

         do i=1, lM%nFa
            call checkIENB(lM, lM%fa(i))
         end do

         return
         end subroutine checkIEN
         !=============================================================
!        Checks IEN array of face depending on the element type. Resets
!        ordering of nodes if any inconsistency is found
         subroutine checkIENB(lM, lFa)
         implicit none
         type(meshType), intent(in) :: lM
         type(faceType), intent(inout) :: lFa

         integer a, b, e, Ec, Ac, sn
         real(kind=8), allocatable :: v(:,:), xl(:,:)

         allocate(v(lM%nsd,lM%eNoN), xl(lM%nsd,lM%eNoN))

         v  = 0.0d0
         xl = 0.0d0
         select case (lFa%eType)
         case (eType_LIN1)
            do e=1, lFa%nEl
               Ec = lFa%gE(e)
               do a=1, lFa%eNoN
                  Ac = lFa%IEN(a,e)
                  xl(:,a) = lM%x(:,Ac)
               end do
               xl(:,3) = SUM(lM%x(:,lM%IEN(:,Ec)),DIM=2) / &
                  real(lM%eNoN)

               a = 1; b = 1
               v(:,1) = xl(:,2) - xl(:,1)
               v(:,2) = xl(:,3) - xl(:,1)
               sn     = SGN(v(1,1)*v(2,2)-v(2,1)*v(1,2))

               if (sn .eq. -1) then
                  a = 1; b = 2
               else if (sn .eq. 0) then
                  write(stdout,ftab4) "ERROR: line element "// &
                     trim(str(e))//" is distorted"
                  stop
               end if
               call swap(lFa%IEN(a,e), lFa%IEN(b,e))
            end do

         case (eType_TRI3)
            do e=1, lFa%nEl
               Ec = lFa%gE(e)
               do a=1, lFa%eNoN
                  Ac = lFa%IEN(a,e)
                  xl(:,a) = lM%x(:,Ac)
               end do
               xl(:,4) = SUM(lM%x(:,lM%IEN(:,Ec)),DIM=2) / &
                  real(lM%eNoN,kind=8)

               a = 1; b = 1
               v(:,1) = xl(:,2) - xl(:,1)
               v(:,2) = xl(:,3) - xl(:,2)
               v(:,3) = xl(:,4) - xl(:,1)
               sn     = SGN(NORM(v(:,3), CROSS(v(:,1:2))))

               if (sn .eq. 1) then
                  a = 1; b = 2
               else if (sn .eq. 0) then
                  write(stdout,ftab4) "ERROR: triangular face element "// &
                     trim(str(e))//" is distorted"
                  stop
               end if
               call swap(lFa%IEN(a,e), lFa%IEN(b,e))
            end do

         end select

         deallocate(v, xl)

         return
         end subroutine checkIENB
         !=============================================================
!        Computes edge connectivity for an input mesh and stores in a
!        CSR connectivity structure
         subroutine getEdgeConnectivity(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         integer :: a, b, e, i, j

!        Get maximum edge connectivity between nodes. Edge ordering of
!        an element is used to determine edges
         allocate(xnn(lM%nNo))
         xnn = 0
         do e=1, lM%nEl
            do i=1, lM%eNoE
               a = lM%IEN(lM%eOrd(1,i),e)
               b = lM%IEN(lM%eOrd(2,i),e)
               xnn(a) = xnn(a) + 1
               xnn(b) = xnn(b) + 1
            end do
         end do
         maxnn = maxval(xnn)

!        Make a list of all possible edge combinations
         allocate(xnnL(maxnn,lM%nNo))
         xnnL = 0
         xnn  = 0
         do e=1, lM%nEl
            do i=1, lM%eNoE
               a = lM%eOrd(1,i)
               b = lM%eOrd(2,i)
               call addNodePair(lM%IEN(a,e), lM%IEN(b,e))
               call addNodePair(lM%IEN(b,e), lM%IEN(a,e))
            end do
         end do
         deallocate(xnn)

!        Store edge connectivity in CSR format
         j = 0
         do a=1, lM%nNo
            do i=1, maxnn
               b = xnnL(i,a)
               if (b .eq. 0) exit
               j = j + 1
            end do
         end do
         lM%econ%nnz = j

         allocate(lM%econ%prow(lM%nNo+1), lM%econ%pcol(lM%econ%nnz))
         j = 0
         lM%econ%prow(1) = 1
         do a=1, lM%nNo
            do i=1, maxnn
               b = xnnL(i,a)
               if (b .eq. 0) exit
               j = j + 1
               lM%econ%pcol(j) = b
            end do
            lM%econ%prow(a+1) = j + 1
         end do
         deallocate(xnnL)

!        Check if edge connectivity matrix is symmetric
         if (mod(lM%econ%nnz,2) .ne. 0) then
            call debugConn(lM%econ, lM%nNo, "e")
            write(stdout,ftab4) "ERROR: edge connectivity matrix "//&
               "has odd no. of entries"
            stop
         end if

         return
         end subroutine getEdgeConnectivity
         !=============================================================
!        Private function used to add node combination to edge
!        connectivity array. Nodes are sorted while added to a row,
!        corresponding to its connected pair
         subroutine addNodePair(Ac, Bc)
         implicit none
         integer, intent(in) :: Ac, Bc

         logical flag
         integer i, j, k

         flag = .true.
         j = 1
         do i=1, xnn(Ac)
            if (xnnL(i,Ac) .eq. Bc) then
               flag = .false.
               exit
            else if (Bc .gt. xnnL(i,Ac)) then
               j = i + 1
            end if
         end do
         if (flag) then
            if (xnn(Ac) .eq. 0) then
               xnn(Ac) = 1
               xnnL(1,Ac) = Bc
            else
               do k=xnn(Ac), j, -1
                  xnnL(k+1,Ac) = xnnL(k,Ac)
               end do
               xnnL(j,Ac) = Bc
               xnn(Ac) = xnn(Ac) + 1
            end if
         end if

         return
         end subroutine addNodePair
         !=============================================================
!        Computes face connectivity for an input mesh and stores in a
!        CSR connectivity structure. Elements are chosen based on the
!        fact that a pair of elements have a common face. Faces on the
!        domain boundaries, which do not have any common face, have
!        `zero' in their column entries
         subroutine getFaceConnectivity(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         logical flag
         integer :: a, b, e, f, i, j, k, l, n, Ac
         integer, allocatable :: inds(:), jnds(:), eadjL(:,:)

!        For every node, count all the elements connected to it
         allocate(xnne(lM%nNo))
         xnne = 0
         do e=1, lM%nEl
            do a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xnne(Ac) = xnne(Ac) + 1
            end do
         end do
         maxNne = maxval(xnne)

!        For every node, prepare a list of all elements connected to it
         allocate(xneL(maxNne,lM%nNo))
         xnne = 0
         xneL  = 0
         do e=1, lM%nEl
            do a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xnne(Ac) = xnne(Ac) + 1
               xneL(xnne(Ac),Ac) = e
            end do
         end do
         b = 3*maxNne

!        Construct an element adjacency array - for every element, make
!        a sorted list of all elements connected to it
 001     b = b + maxNne
         deallocate(xnne)
         allocate(xnne(lM%nEl))
         if (allocated(eadjL)) deallocate(eadjL)
         allocate(eadjL(b,lM%nEl))
         eadjL = 0
         xnne = 0
         do e=1, lM%nEl
            do a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               do i=1, maxNne
                  f = xneL(i,Ac)
                  if (f .eq. 0) exit
                  if (f .eq. e) cycle
                  flag = .true.
                  k = 1
                  do j=1, xnne(e)
                     if (f .eq. eadjL(j,e)) then
                        flag = .FALSE.
                        exit
                     else if (f .gt. eadjL(j,e)) then
                        k = j + 1
                     end if
                  end do
                  if (flag) then
                     if (xnne(e) .eq. 0) then
                        xnne(e)   = 1
                        eadjL(1,e) = f
                     else
                        do l=xnne(e), k, -1
                           eadjL(l+1,e) = eadjL(l,e)
                        end do
                        eadjL(k,e) = f
                        if (xnne(e) .eq. b) goto 001
                        xnne(e)   = xnne(e) + 1
                     end if
                  end if
               end do
            end do
         end do
         deallocate(xneL)

!        Refine element adjacency to elements with a common face.
!        Note that every face of an element is shared by another element
!        except at the domain boundaries, which are zeroed out
         allocate(xneL(lM%eNoF,lM%nEl))
         xneL = 0
         allocate(inds(lM%fNoN), jnds(lM%fNoN))
         do e=1, lM%nEl
            do i=1, lM%eNoF
               do a=1, lM%fNoN
                  inds(a) = lM%IEN(lM%ford(a,i),e)
               end do
               call sort(inds)
               flag = .false.
               n_srch : do n=1, xnne(e)
                  f = eadjL(n,e)
                  do j=1, lM%eNoF
                     do b=1, lM%fNoN
                        jnds(b) = lM%IEN(lM%ford(b,j),f)
                     end do
                     call sort(jnds)
                     l = 0
                     do k=1, lM%fNoN
                        if (inds(k) .eq. jnds(k)) l = l + 1
                     end do
                     if (l .eq. lM%fNoN) then
                        flag = .true.
                        exit n_srch
                     end if
                  end do
               end do n_srch
               if (flag) xneL(i,e) = f
            end do
         end do

!        Store face connectivity in CSR format
         lM%fcon%nnz = lM%eNoF*lM%nEl
         allocate(lM%fcon%prow(lM%nEl+1), lM%fcon%pcol(lM%fcon%nnz))

         j = 0
         lM%fcon%prow(1) = 1
         do e=1, lM%nEl
            do i=1, lM%eNoF
               f = xneL(i,e)
               j = j + 1
               lM%fcon%pcol(j) = f
            end do
            lM%fcon%prow(e+1) = j + 1
         end do

         deallocate(xnne, xneL, eadjL, inds, jnds)

!        Perform checks on the connectivity matrix
!        Total no. of zeroes must match no. of boundary elements
         j = 0
         do e=1, lM%nEl
            do i=lM%fcon%prow(e), lM%fcon%prow(e+1)-1
               if (lM%fcon%pcol(i) .eq. 0) j = j + 1
            end do
         end do
         maxNne = SUM(lM%fa(:)%nEl)

         if (j .ne. maxNne) then
            call debugConn(lM%fcon, lM%nEl, "f")
            write(stdout,ftab4) "ERROR: failed to get face connectivity"
            stop
         end if

!        Symmetry of face connectivity for interior elements
         i = lM%fcon%nnz - maxNne
         if (mod(i,2) .ne. 0) then
            call debugConn(lM%fcon, lM%nEl, "f")
            write(stdout,ftab4) "ERROR: face connectivity matrix "//&
               "has odd no. of entries"
            stop
         end if

         return
         end subroutine getFaceConnectivity
         !=============================================================
!        Debug subroutine to output connectivities to a file
         subroutine debugConn(conn, n, s)
         implicit none
         type(connType), intent(in) :: conn
         integer, intent(in) :: n
         character, intent(in) :: s

         integer i, r, c, fid

         fid = 1078
         open(fid,file="debug_conn.txt")
         if (s .eq. "f") then
            write(fid,ftab1) "Face connectivity: "
         else if (s .eq. "e") then
            write(fid,ftab1) "Edge connectivity: "
         end if

         do r=1, n
            write(fid,ftab1,advance='no') trim(str(r))//": "
            do i=conn%prow(r), conn%prow(r+1)-1
               c = conn%pcol(i)
               write(fid,'(a)',advance='no') " "//trim(str(c))
            end do
            write(fid,'(a)')
         end do

         close(fid)

         return
         end subroutine debugConn
         !=============================================================
      end module meshUtils

!***********************************************************************
