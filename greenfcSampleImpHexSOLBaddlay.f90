! Subroutina para calcular a função de green local dos sitios de uma
! nanofita de grafeno com borda armchair entre contatos.
! Recebe o centro de banda de cada sítio
! -------   5 ______5     ---------
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/4        4\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   3\______/3    |    \___
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/2        2\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   1\______/1    |    \___
! -------                 ----------
! Autor: Flavinao Jose dos Santos
! Universidade Federal Fluminense
! Novembro de 2013

!INCLUDE "./parametros.f90"
! INCLUDE "./modulos.f90"
! INCLUDE "./invers.f90"
! INCLUDE "./acRBsinf.f90"
! INCLUDE "./greensuperficie.f90"

SUBROUTINE green2(Gamma,g1n,gn1,realz,imgz,sample,spin)

USE parametros
USE f90_kind
USE global

IMPLICIT NONE

INTEGER, PARAMETER :: Psring=(4-2*assm)*P+2-2*assm !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
COMPLEX(KIND=8), PARAMETER :: iimag=(0.d0,1.d0),zero=0.d0,one=1.d0, m_one=-1.d0
INTEGER :: i, j, l, k, m, nimp, cte1, cte2
! INTEGER, INTENT(IN) :: impureza(N,2,P+1)
! REAL(KIND=8), INTENT(IN) :: anderdisor(N,2,PsrHalf)
REAL(KIND=8), INTENT(IN) :: imgz, realz, spin
INTEGER, INTENT(IN) :: sample
COMPLEX(KIND=8), INTENT(OUT), DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: Tnq1, Tnq2, Tnq, gqqinv, gnn, &
                  Tqn_gnn, Tqn_gnnTnq, Tqn_gn1, Tqn, g1nTnq, &
                  Hamt, greencont!, g1nTnqGqq, Gqq, Sigma, SigmaDag
COMPLEX(KIND=8) :: SOed(3,3), SOee(3,3), SOdd(3,3), z, slamb, cjslamb, zmeps0,&
                   mt, spinrho, atimp, cjatimp, null=0.d0
REAL(KIND=8) :: diference, num_aleat
CHARACTER(LEN=60) :: formato, formato2
CHARACTER(LEN=10) :: layer

WRITE(layer,FMT="(f6.3)") realz
CALL printtime("START Green fc   Energy = "//layer)

z=CMPLX(realz,imgz,double)

!Between two rings
Tnq1=0.d0
!Inside of a ring
Tnq2=0.d0
DO i=1, P+assm
   Tnq1(2*i,2*i)=t
END DO
DO i=1, P+1
   Tnq2(2*i-1,2*i-1)=t
END DO

slamb=lambda*spin
if(spin.LE.0.d0) slamb=CONJG(lambda)
cjslamb=CONJG(slamb)
! spinrho=spin*rho
spinrho=rho
atimp=spin*aatimp
if(spin.LE.0.d0) atimp=CONJG(aatimp)
cjatimp=CONJG(atimp)

SOed(1,1:3)=[ atimp , slamb ,spinrho]
SOed(2,1:3)=[ slamb ,spinrho,cjslamb]
SOed(3,1:3)=[spinrho,cjslamb,cjatimp]

SOee(1,1:3)=[  null ,cjatimp,cjslamb ]
SOee(2,1:3)=[ atimp ,  null ,cjatimp ]
SOee(3,1:3)=[ slamb ,  atimp, null   ]

SOdd=TRANSPOSE(SOee)

Hamt=0.d0
zmeps0=z-eps0
mt=-t

Hamt(1,1)=zmeps0
Hamt(1,2)=mt
DO i=2, PsrHalf-1
   Hamt(i,i-1)=mt
   Hamt(i,i)=zmeps0
   Hamt(i,i+1)=mt
END DO
Hamt(PsrHalf,PsrHalf-1)=mt
Hamt(PsrHalf,PsrHalf)=zmeps0

CALL printtime("START Green surface")
CALL greensupfcontato(z,greencont,eps0)
! greencont=0.d0
CALL printtime("END Green surface")

gnn=greencont

DO i=1, N
   DO j=1,2
      gqqinv=Hamt
      
      IF(ABS(deltau).GE.1.d-12) THEN
         DO k=1, PsrHalf
            gqqinv(k,k)=gqqinv(k,k)+anderdisor(sample,i,j,k)
         END DO
      END IF
      
     IF(j==1) THEN
         Tnq=Tnq1
         cte1=0
         cte2=1
      ELSE
         Tnq=Tnq2
         cte1=1
         cte2=0
      END IF
      !Impurities at the left
      nimp=impureza(sample,i,j,1)
      DO k=1, nimp
         l=2*impureza(sample,i,j,1+k)-cte1
         Tnq(l:l+2,l:l+2)=Tnq(l:l+2,l:l+2)+SOed
         DO m=l, l+2
!             diference=ABS( gqqinv(m,m)-(zmeps0+VV) )
!             IF(diference>1.d-10) gqqinv(m,m)=gqqinv(m,m)+VV
            gqqinv(m,m)=gqqinv(m,m)+VV
         END DO
         gqqinv(l:l+2,l:l+2)=gqqinv(l:l+2,l:l+2)-SOdd
      END DO
      !Impurities at the right
      IF(.NOT.(i.EQ.N .AND. j.EQ.2)) THEN
         nimp=impureza(sample,i+cte1,2-cte1,1  )
         DO k=1, nimp
            l= 2*impureza(sample,i+cte1,2-cte1,1+k)-cte2
            DO m=l, l+2
!                diference=ABS( gqqinv(m,m)-(zmeps0+VV) )
!                IF(diference>1.d-10) gqqinv(m,m)=gqqinv(m,m)+VV
               gqqinv(m,m)=gqqinv(m,m)+VV
            END DO
            gqqinv(l:l+2,l:l+2)=gqqinv(l:l+2,l:l+2)-SOee
         END DO
      END IF

      IF(sitio.NE.0 .AND. i==N/2 .AND. j==1) THEN
         gqqinv(sitio,sitio)=1.d12
!          gqqinv(PsrHalf-sitio+1,PsrHalf-sitio+1)=-2.d2
      END IF
      IF( nlin.NE.0 .AND. (i>=cami.AND.i<=camf) ) THEN
         DO k=1, nlin
            gqqinv(k,k)=gqqinv(k,k)+1.d12
            gqqinv(PsrHalf+1-k,PsrHalf+1-k)=gqqinv(k,k)+1.d12
         END DO
      END IF

      Tqn=CONJG(TRANSPOSE(Tnq))
      
!       Tqn_gnn = MATMUL(Tqn,gnn)
!       Tqn_gnnTnq = MATMUL(Tqn_gnn,Tnq)
!       gnn=gqqinv-Tqn_gnnTnq
                                                              !Tqn
      CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf,   one, Tqn    , PsrHalf, gnn, PsrHalf, zero, Tqn_gnn, PsrHalf ) ! Tqn_gnn = Tqn * gnn
      CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, m_one, Tqn_gnn, PsrHalf, Tnq, PsrHalf, one , gqqinv , PsrHalf ) ! gqqinv = gqqinv - Tqn_gnn * Tnq
      gnn=gqqinv !gqqinv-Tqn_gnnTnq

      CALL invers(gnn,PsrHalf)
      
      IF( i==1 .AND. j==1 ) THEN
         g1n=gnn
         gn1=gnn
      ELSE
!          g1nTnq=MATMUL(g1n,Tnq)
!          g1n=MATMUL(g1nTnq,gnn) !g1nTnqGqq
         CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, g1n   , PsrHalf, Tnq, PsrHalf, zero, g1nTnq, PsrHalf ) 
         CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, g1nTnq, PsrHalf, gnn, PsrHalf, zero, g1n   , PsrHalf ) !g1nTnqGqq

!          Tqn_gn1=MATMUL(Tqn,gn1)
!          gn1=MATMUL(gnn,Tqn_gn1) !GqqTqn_gn1
         CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, Tqn, PsrHalf, gn1    , PsrHalf, zero, Tqn_gn1, PsrHalf ) 
         CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, gnn, PsrHalf, Tqn_gn1, PsrHalf, zero, gn1    , PsrHalf ) !GqqTqn_gn1
      END IF

   END DO
   WRITE(layer,FMT="(i0)") i
   CALL printtime("END Add Layer = "//layer)
END DO

Tqn=Tnq1
!Obtendo Gnn, Gn1 e G1n, func.g do sample apos conectar lead da
!direita. 1 eh o indice da primeira linha da primeria camada do sample,
!n eh a segunda linha da ultima camada.

! Tqn_gnn = MATMUL(Tqn,greencont)
! Tqn_gnnTnq = MATMUL(Tqn_gnn,Tqn) !Tqn_gqqTqn
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, Tqn    , PsrHalf, greencont, PsrHalf, zero, Tqn_gnn   , PsrHalf ) 
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, Tqn_gnn, PsrHalf, Tqn      , PsrHalf, zero, Tqn_gnnTnq, PsrHalf ) !Tqn_gqqTqn

gqqinv=gnn
CALL invers(gqqinv,PsrHalf)
gnn=gqqinv-Tqn_gnnTnq !Gnn^-1
CALL invers(Gnn,PsrHalf) !Gnn

! Tqn_gnn=MATMUL(g1n,gqqinv) !g1n_gnn^-1
! g1n=MATMUL(Tqn_gnn,gnn) !g1n_gnn^-1_Gnn
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, g1n    , PsrHalf, gqqinv, PsrHalf, zero, Tqn_gnn, PsrHalf ) !g1n_gnn^-1
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, Tqn_gnn, PsrHalf, gnn   , PsrHalf, zero, g1n    , PsrHalf ) !g1n_gnn^-1_Gnn

! Tqn_gnn=MATMUL(gqqinv,gn1) !gnn^-1_gn1
! gn1=MATMUL(gnn,Tqn_gnn) !Gnn_gnn^-1_gn1
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, gqqinv, PsrHalf, gn1    , PsrHalf, zero, Tqn_gnn, PsrHalf ) !gnn^-1_gn1
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, gnn   , PsrHalf, Tqn_gnn, PsrHalf, zero, gn1    , PsrHalf ) !Gnn_gnn^-1_gn1

Tqn_gnn = CONJG(TRANSPOSE(Tqn_gnnTnq)) !Sigma+
Gamma=iimag*(Tqn_gnnTnq-Tqn_gnn) !i(Sigma - Sigma+)

WRITE(layer,FMT="(f6.3)") realz
CALL printtime("END Green fc   Energy = "//layer)
IF(myrank==0)PRINT *

RETURN
END SUBROUTINE green2
