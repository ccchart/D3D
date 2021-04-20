  SUBROUTINE ALREADYTRI(K1,K2,K3,JA)
  use m_netw
  implicit none
  integer :: K1,K2,K3,JA

  integer :: n1
  integer :: n2
  integer :: n3
  integer :: np
  JA = 0
  DO NP = NUMP, 1, -1
     IF (netcell(NP)%N .EQ. 3) THEN
        N1 = netcell(NP)%NOD(1)
        N2 = netcell(NP)%NOD(2)
        N3 = netcell(NP)%NOD(3)
        IF ((K1.EQ.N1 .OR. K1.EQ.N2  .OR. K1.EQ.N3 ) .AND. &
            (K2.EQ.N1 .OR. K2.EQ.N2  .OR. K2.EQ.N3 ) .AND. &
            (K3.EQ.N1 .OR. K3.EQ.N2  .OR. K3.EQ.N3 ) ) THEN
            JA = np
            call qnerror('already 3', ' ',' ')
            RETURN
        ENDIF
     ENDIF
  ENDDO
  RETURN
  END SUBROUTINE ALREADYTRI
