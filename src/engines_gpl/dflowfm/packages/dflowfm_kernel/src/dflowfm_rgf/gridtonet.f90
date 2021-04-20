     subroutine gridtonet()
      use m_netw
      use m_grid
      use m_missing
      use gridoperations

      implicit none
      double precision :: af

      integer, allocatable :: mn(:,:)
      double precision XX(8), YY(8), ZZ(8), tooclose0, length
      integer :: k0, l0, ja, jadoorladen, i, j, k, l, method, ierr, key

      jadoorladen = 1

      IF (JADOORLADEN .EQ. 0) THEN
         K0 = 0
         L0 = 0
      ELSE
         K0 = NUMK
         L0 = NUML
      ENDIF

      K  = 0
      DO I = 1,MC                                    ! COUNT NR OF NODES
         DO J = 1,NC
            IF (Xc(I,J) .NE. dXYMIS) THEN
               K = K + 1
            ENDIF
         ENDDO
      ENDDO

      !IF (K0+K .GT. SIZE(XK) ) THEN
         CALL INCREASENETW(K0+K,L0+4*K)
      !ENDIF

      K  = K0
      L  = L0                                        ! COUNT MAX NR OF ATTACHED LINKS PER NODE


      CALL READYY('Arranging curvilinear grid-in network',0d0)


      if (allocated (mn) )  deallocate(mn)
      allocate ( mn(mc,nc) , stat = ierr ) ; mn = 0
      call aerr('mn(mc,nc)', ierr , mc*nc)

      DO I = 1,MC
         DO J = 1,NC
            if (xc(i,j) .ne. dxymis) then
               call addnetpointnocheck(xc(i,j), yc(i,j), zc(i,j), k )
               mn(i,j) = k
            endif
         ENDDO
      ENDDO
      numk = k

      af = 0.2d0
      CALL READYY('Arranging curvilinear grid-in network',af)


      DO I = 1,MC-1
         DO J = 1,NC
            if ( mn(i,j) .ne. 0 .and. mn(i+1,j) .ne. 0 )  then
               L       = L + 1
               kn(1,L) = mn(i,j)
               kn(2,L) = mn(i+1,j)
               KN(3,L) = 2
            endif
         ENDDO
      ENDDO

      af = 0.4d0
      CALL READYY('Arranging curvilinear grid-in network',af)


      DO I = 1,MC
         DO J = 1,NC-1
            if ( mn(i,j) .ne. 0 .and. mn(i,j+1) .ne. 0 )  then
               L       = L + 1
               kn(1,L) = mn(i,j)
               kn(2,L) = mn(i,j+1)
               KN(3,L) = 2
            endif
         ENDDO
      ENDDO

      af = 0.6d0
      CALL READYY('Arranging curvilinear grid-in network',af)

      numl = l
      call setnodadm(0)

      CALL READYY('Arranging curvilinear grid-in network',-1d0)

      ! call copydeptosam()


      if (k0 > 0) then


      JA = 1

      call readyy('Merging networks', 0d0)
      call findcells(0)

!     merge nodes

      if ( tooclose.gt.1d-16 .and. k0 > 0) then
         CALL CONFRM('MERGE NODES ? ',JA)
         IF (JA == 1) call MERGENODESINPOLYGON()
      end if

!     merge boundary nodes
!      call mergenet()

      call readyy('Merging networks', -1d0)

      endif

!     set network status
      netstat = NETSTAT_CELLS_DIRTY

      END subroutine gridtonet
