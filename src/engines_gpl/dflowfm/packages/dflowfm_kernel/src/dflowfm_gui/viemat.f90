  !
  !
  subroutine viemat(the,phi)
  implicit none
  double precision :: cp
  double precision :: ct
  double precision :: deltx
  double precision :: delty
  double precision :: deltz
  double precision :: dscr
  double precision :: r
  double precision :: sp
  double precision :: st
  double precision :: t1
  double precision :: t2
  double precision :: t3
  double precision :: t4
  double precision :: vs
  double precision :: wpqr
  double precision :: x0s
  double precision :: y0s
  double precision :: z
  double precision :: zfac
  double precision :: zupw
  double precision :: the, phi
  common /viewmat/ vs(4,4), x0s, y0s
  COMMON /PERSPX/ WPQR,DELTX,DELTY,DELTZ,ZFAC,DSCR,ZUPW

  ! Maak viewing matrix Vs
  ! phi (0 -- pi) en the (-pi/2 -- pi/2) : kijkhoekjes
  ! wpqr                                 : oog-object in wereldcoor
  ! deltx,delty,deltz                    : kijk door dit punt in were
  ! zfac (negatief:op z'n kop)           : oprekking verticaal
  ! Dscr                                 : oog-scherm in wereldcoor
  ! Vs                                   : Viewing matrix
  !
  dimension T1(4,4),T2(4,4),T3(4,4),T4(4,4),R(4,4),Z(4,4)
  T1 = 0
  T2 = 0
  T3 = 0
  T4 = 0
  R = 0
  CT = COS(THE)
  ST = SIN(THE)
  CP = COS(PHI)
  SP = SIN(PHI)

  T1(1,1) = 1.
  T1(2,2) = 1.
  T1(3,3) = zfac
  T1(4,4) = 1.
  T1(1,4) = -deltx
  T1(2,4) = -delty
  T1(3,4) = -deltz*ZFAC

  T2(1,1) = 1.
  T2(2,2) = 1.
  T2(3,3) = 1.
  T2(4,4) = 1.
  T2(1,4) = -wpqr*ct*cp
  T2(3,4) = -wpqr*ct*sp
  T2(2,4) = -wpqr*st

  T3(1,1) =  cp
  T3(1,3) =  sp
  T3(3,1) = -sp
  T3(3,3) =  cp
  T3(2,2) =  1.
  T3(4,4) =  1.

  T4(1,1) =  ct
  T4(1,2) =  st
  T4(2,1) = -st
  T4(2,2) =  ct
  T4(3,3) =  1.
  T4(4,4) =  1.

  R(1,3)  = Dscr
  R(2,2)  = Dscr
  R(3,1)  = -1.
  R(4,4)  =  1.

  !  nadat alles geinitialiseerd is de viewing transformatie-matr Vs =

  call matm4(R,T4,Z)
  call matm4(Z,T3,Vs)
  call matm4(Vs,T2,Z)
  call matm4(Z,T1,Vs)
  end subroutine viemat
