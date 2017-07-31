SUBROUTINE sortanderdisor()

USE parametros
IMPLICIT NONE
INTEGER :: i,j,l,k,m,q,w,y,u
INTEGER, PARAMETER :: Psring=(4-2*assm)*P+2-2*assm !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
! REAL(KIND=8), INTENT(OUT) :: anderdisor(namostras,N,2,PsrHalf)
REAL(KIND=8) ::  num_aleat

IF(deltau.LT.1.d-12) RETURN

ALLOCATE( anderdisor(namostras,N,2,PsrHalf) )

CALL inicializa_semente()

DO m=1,namostras
   DO i=1, N; DO j=1,2;
      DO l=1, PsrHalf
         CALL RANDOM_NUMBER(num_aleat)
         anderdisor(m,i,j,l)=(num_aleat*2.d0-1.d0)*deltau
      END DO
   END DO; END DO
END DO

RETURN
END SUBROUTINE sortanderdisor
