! Subroutina para calcular a função de green local dos sitios de uma
! nanofita de grafeno com borda armchair entre contatos.
! Recebe o centro de banda de cada sítio
! -------   5 ______6     ---------
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/4        7\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   3\______/8    |    \___
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/2        9\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   1\______/10   |    \___
! -------                 ----------
! Autor: Flavinao Jose dos Santos
! Universidade Federal Fluminense
! Julho de 2013

!INCLUDE "./parametros.f90"
! INCLUDE "./modulos.f90"
! INCLUDE "./invers.f90"
! INCLUDE "./acRBsinf.f90"
! INCLUDE "./greensuperficie.f90"

!Paramentros que devem ser ajustados no arquivo paramentros.f90:
! INTEGER, PARAMETER ::  P=11, N=11
! REAL(KIND=8), PARAMETER :: t=-2.7d0, eps0=0.d0, &
!          fractVV=0.0d0, fractlmda=0.2d0, fracttimp=0.0d0
! REAL(KIND=8):: eps1=eps0, VV=fractVV*ABS(t), &
!                lambda=fractlmda*ABS(t), timp=fracttimp*t

SUBROUTINE green1(Gamma,g1n,gn1,realz,imgz,sample,spin)

USE parametros
USE f90_kind
IMPLICIT NONE

INTEGER, PARAMETER :: Psring=4*P+2 - 2*assm !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
COMPLEX, PARAMETER :: iimag=(0.d0,1.d0)
INTEGER :: i, j, l, k,m,o, nimpesq, nimpdir, tamcamd, nimpesqPsrHalf, &
           nimpesqPsring, nimpdirR, nimpesqR !R de Real. Apenas se R/=0 as
           !impurezas serao adicionadas efetivamente e nao apenas o efeito.
! INTEGER, INTENT(IN) :: impureza(N,2,P+1)
REAL(KIND=8), INTENT(IN) :: imgz, realz, spin
INTEGER, INTENT(IN) :: sample
COMPLEX(KIND=8), INTENT(OUT), DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1
COMPLEX(KIND=8), DIMENSION(:,:) , ALLOCATABLE:: Tnq, Tqn, g1nTnq, Tqn_gn1, &
                           g11inv, Tqn_gnn, Tqn_gnnTnq, Gqq, G1q, Gq1, BBB
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: gnn, greencont, aux, aux2
COMPLEX(KIND=8) :: z, spinrho, zmeps0, & 
                   zmeps1, mt, vec(3), atimp, cjatimp, slamb, cjslamb, &
                   SOed(3,3), SOde(3,3), SOee(3,3), SOdd(3,3), VVcp
CHARACTER(LEN=60) :: formato, formato2

z=CMPLX(realz,imgz,double)
CALL greensupfcontato(z,greencont,eps0)
gnn=greencont

!***************************************************
slamb=spin*lambda
cjslamb=CONJG(slamb)
spinrho=spin*rho
atimp=spin*aatimp
cjatimp=CONJG(atimp)
VVcp=-VV

SOed(1,1:3)=[ atimp , slamb ,spinrho]
SOed(2,1:3)=[ slamb ,spinrho,cjslamb]
SOed(3,1:3)=[spinrho,cjslamb,cjatimp]

SOee(1,1:3)=[ VVcp ,cjatimp,cjslamb ]
SOee(2,1:3)=[ atimp , VVcp ,cjatimp ]
SOee(3,1:3)=[ slamb ,  atimp, VVcp  ]

SOdd=TRANSPOSE(SOee)
SOde=CONJG(TRANSPOSE(SOed))
!***************************************************

zmeps0=z-eps0
zmeps1=z-eps1*spin
mt=-t
vec=(/mt,zmeps0,mt/)

DO i=1, N
   nimpesq=impureza(sample,i,1,1)
   nimpdir=impureza(sample,i,2,1)
   IF(fracttimp==0.d0) THEN
      nimpesqR=0
      nimpdirR=0
   ELSE
      nimpesqR=nimpesq
      nimpdirR=nimpdir
   END IF

   tamcamd=nimpesqR+nimpdirR+Psring
   nimpesqPsrHalf=nimpesqR+PsrHalf
   nimpesqPsring=nimpesqR+Psring

   ALLOCATE( Tnq(PsrHalf,tamcamd) )
   ALLOCATE( Tqn(tamcamd,PsrHalf) )
   ALLOCATE( g11inv(tamcamd,tamcamd) )

   !Os indices: 'n' é a ultima camada jah adic. E 'q' eh a camada sendo adic.
   Tnq=0.d0

   !Hopping de acoplamento devido aos sitios regulares na camada a ser adic.
   DO j=1, P
      Tnq(PsrHalf+1-2*j,nimpesqR+2*j)=t
   END DO

   !Hopping de acoplamento devido as impurezas a esq. na camada a ser
   !adicionada
   DO l=1, nimpesqR
      k=PsrHalf+1-2*impureza(sample,i,1,1+l)
      Tnq(k-2:k,l)=timp
   END DO

   !Acomplamento spin-orbita
   DO j=1, nimpesq                        !1  k ___   3
      !Acomplamento Spin-Orbita           !   /    \    .
      l=PsrHalf+1-2*impureza(sample,i,1,1+j)     !2 /     m\ 2
      m=nimpesqR+2*impureza(sample,i,1,1+j)+1    !  \     /  
      k=l-2                               !3 l\___/   1
      Tnq(l:k:-1,m-1:m+1)=Tnq(l:k:-1,m-1:m+1)+SOed
   END DO

   Tqn=CONJG(TRANSPOSE(Tnq))

   g11inv=0.d0

   !As impurezas a esquerda ocuparao as primeiras linhas da matriz
   DO j=1, nimpesqR
      g11inv(j,j)=zmeps1
      l=nimpesqR+2*impureza(sample,i,1,1+j)
      g11inv(j,l:l+2)=-timp
      g11inv(l:l+2,j)=-timp
   END DO

   !Preenchendo as linhas correspondentes aos sitios padroes da camada
   !adicionada. Ficarao nas linhas centras, entre as impurezas a
   !esquerda e a direita
   g11inv(nimpesqR+1,nimpesqR+1:nimpesqR+2)=vec(2:3)
   DO j=nimpesqR+2, nimpesqPsring-1
      g11inv(j,j-1:j+1)=vec
   END DO
   g11inv(nimpesqPsring,nimpesqPsring-1:nimpesqPsring)=vec(1:2)
   !Faz as ligacoes dos sitios nao vizinhos
   DO j=1, P
      g11inv(nimpesqR+2*j-1,nimpesqPsring-2*j+2)=mt
      g11inv(nimpesqPsring-2*j+2,nimpesqR+2*j-1)=mt
   END DO
   
   !Reajustando os centros de banda dos sitios vizinhos a impurezas e
   !ligando o acoplamento spin-orbita
   DO j=1, nimpesq                         !    ___k  3  
      !Acomplamento Spin-Orbita            !  /    \             .
      l=nimpesqR+2*impureza(sample,i,1,1+j)       ! /  *   \ 2 
      k=l+2                                ! \     /    
      g11inv(l:k,l:k)=g11inv(l:k,l:k)-SOdd !  \___/l  1 
   END DO

   DO j=1, nimpdir                                    !3   ____ k 1
      !Acomplamento Spin-Orbita                       !   /    \       . 
      l=nimpesqR+2*impureza(sample,i,2,1+j)-1                !2 /l     \ 2
      k=nimpesqPsring-2*impureza(sample,i,2,1+j)             !  \     /   
      m=k+2                                           !1  \___/m  3
      g11inv(l:l+2 ,l:l+2 )=g11inv(l:l+2 ,l:l+2 )-SOee
      g11inv(m:k:-1,m:k:-1)=g11inv(m:k:-1,m:k:-1)-SOdd
      g11inv(l:l+2 ,m:k:-1)=g11inv(l:l+2 ,m:k:-1)-SOed
      g11inv(m:k:-1,l:l+2 )=g11inv(m:k:-1,l:l+2 )-SOde
   END DO
   !To account the impurities in the left of the next layer
   IF (i.NE.N) THEN                                        !         1___  
      DO j=1, impureza(sample,i+1,1,1)                            !        /l
         !Acomplamento Spin-Orbita                         !    ___/ 2 *
         l=nimpesqPsring-2*impureza(sample,i+1,1,1+j)-1           !  /    \         .
         k=l+2                                             ! /      \k___
         g11inv(k:l:-1,k:l:-1)=g11inv(k:l:-1,k:l:-1)-SOee  ! \     / 3 
      END DO                                               !  \___/
   END IF

   !Preenchendo para as impurezas a direita, ultimas linhas da matriz
   l=1
   DO j=nimpesqPsring+1,tamcamd
      g11inv(j,j)=zmeps1
      DO k=0, 2
         m=nimpesqR+2*impureza(sample,i,2,1+l)-1+k
         o=nimpesqPsring-2*impureza(sample,i,2,1+l)+k
         g11inv(j,m)=-timp
         g11inv(m,j)=-timp
         g11inv(j,o)=-timp
         g11inv(o,j)=-timp
      END DO
      l=l+1
   END DO

   ALLOCATE( Tqn_gnn(tamcamd,PsrHalf),Tqn_gnnTnq(tamcamd,tamcamd), &
                                                   Gqq(tamcamd,tamcamd) )
   Tqn_gnn = MATMUL(Tqn,gnn)
   Tqn_gnnTnq = MATMUL(Tqn_gnn,Tnq)
   Gqq=g11inv-Tqn_gnnTnq
   CALL invers(Gqq,tamcamd)

   DEALLOCATE( Tqn_gnn, Tqn_gnnTnq, g11inv )

   ALLOCATE( G1q(PsrHalf,tamcamd), Gq1(tamcamd,PsrHalf), &
                        g1nTnq(PsrHalf,tamcamd), Tqn_gn1(tamcamd,PsrHalf) )

   gnn=Gqq(nimpesqPsrHalf+1:nimpesqPsring, nimpesqPsrHalf+1:nimpesqPsring)

   IF( i==1 ) THEN
      g1n=Gqq( nimpesqR+1:nimpesqPsrHalf, nimpesqPsrHalf+1:nimpesqPsring)
      gn1=Gqq( nimpesqPsrHalf+1:nimpesqPsring, nimpesqR+1:nimpesqPsrHalf)
   ELSE
      g1nTnq=MATMUL(g1n,Tnq)
      G1q=MATMUL(g1nTnq,Gqq) !g1nTnqGqq
      g1n=G1q( : , nimpesqPsrHalf+1:nimpesqPsring)

      Tqn_gn1=MATMUL(Tqn,gn1)
      Gq1=MATMUL(Gqq,Tqn_gn1) !GqqTqn_gn1
      gn1=Gq1(nimpesqPsrHalf+1:nimpesqPsring , :)
   END IF

   DEALLOCATE( Tqn, Tnq, Gqq, G1q, Gq1, g1nTnq, Tqn_gn1 )

END DO
ALLOCATE(Tqn(PsrHalf,PsrHalf), Tqn_gnn(PsrHalf,PsrHalf), Gqq(PsrHalf,PsrHalf), &
         Tqn_gnnTnq(PsrHalf,PsrHalf), g11inv(PsrHalf,PsrHalf) )

Tqn=(0.d0,0.d0)
DO i=1, P
   Tqn(PsrHalf+1-2*i,2*i)=t
END DO

!Obtendo Gnn, Gn1 e G1n, func.g do sample apos conectar lead da
!direita. 1 eh o indice da primeira linha da primeria camada do sample,
!n eh a segunda linha da ultima camada.
Tqn_gnn = MATMUL(Tqn,greencont)
Tqn_gnnTnq = MATMUL(Tqn_gnn,Tqn) !Tqn_gqqTqn
g11inv=gnn
CALL invers(g11inv,PsrHalf)
Gqq=g11inv-Tqn_gnnTnq !Gnn^-1
CALL invers(Gqq,PsrHalf) !Gnn

Tqn_gnn=MATMUL(g1n,g11inv) !g1n_gnn^-1
g1n=MATMUL(Tqn_gnn,Gqq) !g1n_gnn^-1_Gnn

Tqn_gnn=MATMUL(g11inv,gn1) !gnn^-1_gn1
gn1=MATMUL(Gqq,Tqn_gnn) !Gnn_gnn^-1_gn1

Tqn_gnn = CONJG(TRANSPOSE(Tqn_gnnTnq)) !Sigma+
Gamma=iimag*(Tqn_gnnTnq-Tqn_gnn) !i(Sigma - Sigma+)

RETURN

END SUBROUTINE green1


! IF (i==2) THEN
! !==========================================================================
! WRITE(formato, FMT="(a,i0,a)")  '( ',tamcamd, '(f5.1) )'
!
! DO j=1, tamcamd
!    WRITE( UNIT=*, FMT=formato) REAL(g11inv(j,:))+DIMAG(g11inv(j,:))
! END DO
! !==========================================================================
! stop
! END IF



! IF (i==4) THEN
! print *, nimpesq, nimpdir, impureza(sample,i,1,2), impureza(sample,i,2,2)
! !==========================================================================
! WRITE(formato, FMT="(a,i0,a,a,i0,a,a,i0,a  )") &
!    '( ', nimpesqR, '(f4.0),', "'|'," , Psring, '(f4.0),', "'|'," , nimpdirR, '(f4.0) )'
! 
! WRITE(formato2, FMT="(a,i0,a)") "(",9*(tamcamd),"a)"
! 
! DO j=1, tamcamd
! 
!    WRITE( UNIT=*, FMT=formato) Dimag(g11inv(j,:))
!    IF(j==nimpesqR.OR.j==nimpesqR+Psring) THEN
!       print formato2, ("-",l=1,4*(tamcamd)+2)
!    ELSE
!       !print *
!    END IF
! END DO
! !==========================================================================
! STOP
! END IF
