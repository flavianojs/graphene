! Subroutina para calcular a função de green local dos sitios de uma
! nanofita de grafeno com borda armchair entre contatos.
! Recebe o centro de banda de cada sítio
! -------   5 ______10    ---------
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/4        9\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   3\______/8    |    \___
! \     |    /      \     |    /
!  \    |   /        \    |   /
!   \___|__/2        7\___|__/
!   /   |  \          /   |  \
!  /    |   \        /    |   \
! /     |   1\______/6    |    \___
! -------                 ----------
! Autor: Flavinao Jose dos Santos
! Universidade Federal Fluminense
! Novembro de 2013

!INCLUDE "./parametros.f90"
! INCLUDE "./modulos.f90"
! INCLUDE "./invers.f90"
! INCLUDE "./acRBsinf.f90"
! INCLUDE "./sqRBsinf.f90"
! INCLUDE "./greensuperficie.f90"

SUBROUTINE green(Hamt,Gamma,g1n,gn1,realz,imgz,sample,spin)

USE parametros
USE f90_kind
USE global
IMPLICIT NONE

INTEGER, PARAMETER :: Psring=4*P+2 !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
COMPLEX(KIND=8), PARAMETER :: iimag=(0.d0,1.d0),zero=0.d0,one=1.d0, m_one=-1.d0
INTEGER :: i, j, l, k,m,o, nimpesq, nimpdir, tamcamd, nimpesqPsrHalf, &
           nimpesqPsring, nimpdirR, nimpesqR
! INTEGER, INTENT(IN) :: impureza(N,2,P+1)
REAL(KIND=8), INTENT(IN) :: imgz, realz, spin
INTEGER, INTENT(IN) :: sample
COMPLEX(KIND=8), INTENT(OUT), DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1
COMPLEX(KIND=8), INTENT(OUT), DIMENSION(PsrN,PsrN) :: Hamt
! COMPLEX(KIND=8), INTENT(OUT), ALLOCATABLE :: Hamt(:,:)
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: Hlinha, tt, AC, Hentlin, Hentcam, Sigma, SigmaDag
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrN) :: tpc, tqc
! COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tpc, tqc, tcq, tcp, AB, &
!                                                 SelfEsq, SelfDir
COMPLEX(KIND=8), DIMENSION(PsrN,PsrHalf) :: tcq, tcp, AB
COMPLEX(KIND=8) :: SOde(3,3), SOed(3,3), Hcamd(Psring,Psring), SOee(3,3), SOdd(3,3)
REAL(KIND=8):: lcharg(3,3)
COMPLEX(KIND=8), DIMENSION(PsrN,PsrN) :: SelfEsq, SelfDir
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: greencont
COMPLEX(KIND=8) :: z, lambdaspin, cjlambdaspin, atimp, cjatimp
CHARACTER(LEN=60) :: formato, formato2
CHARACTER(LEN=10) :: layer

WRITE(layer,FMT="(f6.3)") realz
CALL printtime("START Green fc   Energy = "//layer)

!  --    --
! |  |__|  |
! |  |  |  |
!  --    --
! |  |__|  |
! |  |  |  |
!  --    --
! l  l
! \--/
! camada
!Begin: Constructing the hamiltonian of the central part

Hlinha=0.d0 !Hoppings entre atomos de uma mesma linha
Hlinha(1,1)=eps0
Hlinha(1,2)=t
Hlinha(1,3)=t2
Hlinha(2,1)=t
Hlinha(2,2)=eps0
Hlinha(2,3)=t
Hlinha(2,4)=t2
DO i=3, PsrHalf-2
   Hlinha(i,i-2)=t2
   Hlinha(i,i-1)=t
   Hlinha(i,i)=eps0
   Hlinha(i,i+1)=t
   Hlinha(i,i+2)=t2
END DO
Hlinha(PsrHalf-1,PsrHalf-3)=t2
Hlinha(PsrHalf-1,PsrHalf-2)=t
Hlinha(PsrHalf-1,PsrHalf-1)=eps0
Hlinha(PsrHalf-1,PsrHalf  )=t
Hlinha(PsrHalf,PsrHalf-2)=t2
Hlinha(PsrHalf,PsrHalf-1)=t
Hlinha(PsrHalf,PsrHalf  )=eps0

Hentlin=0.d0 !Hopping entre atomos de diferentes linhas
Hentlin(1,1)=tb  
DO i=2, P !Hopping nn
   Hentlin(2*i-1,2*i-1)=t  
END DO
Hentlin(PsrHalf,PsrHalf)=tb  

Hentlin(1,2)=t2 !Hopping nnn
DO i=2, PsrHalf-1
   Hentlin(i,i+1)=t2
   Hentlin(i,i-1)=t2
END DO
Hentlin(PsrHalf,PsrHalf-1)=t2

Hentcam=0.d0 !Hopping entre atomos de diferente camadas
DO i=1, P
   Hentcam(2*i,2*i)=t
END DO
Hentcam(1,2)=t2 !Hopping nnn
DO i=2, PsrHalf-1
   Hentcam(i,i+1)=t2
   Hentcam(i,i-1)=t2
END DO
Hentcam(PsrHalf,PsrHalf-1)=t2

Hcamd(1:PsrHalf,1:PsrHalf)=Hlinha
Hcamd(PsrHalf+1:Psring,PsrHalf+1:Psring)=Hlinha
Hcamd(1:PsrHalf,PsrHalf+1:Psring)=Hentlin
Hcamd(PsrHalf+1:Psring,1:PsrHalf)=Hentlin

lambdaspin=lambda*spin
if(spin.LE.0.d0) lambdaspin=CONJG(lambda)
cjlambdaspin=CONJG(lambdaspin)
SOed=0.d0
SOed(1,2)=  lambdaspin
SOed(2,1)=  lambdaspin
SOed(2,3)=cjlambdaspin
SOed(3,2)=cjlambdaspin

SOed(1,3)= rho!*spin
SOed(2,2)= rho!*spin
SOed(3,1)= rho!*spin

atimp=aatimp*spin
if(spin.LE.0.d0) atimp=CONJG(aatimp)
cjatimp=CONJG(atimp)
SOed(1,1)=  atimp
SOed(3,3)=cjatimp

SOee=0.d0
SOee(1,3)=cjlambdaspin
SOee(3,1)=  lambdaspin

SOee(3,2)=  atimp
SOee(2,3)=cjatimp
SOee(2,1)=  atimp
SOee(1,2)=cjatimp

SOdd=(TRANSPOSE(SOee))
SOde=CONJG(TRANSPOSE(SOed))

!lcharg = Local Charge, if to account for local impurity charge screening
lcharg=0.d0
DO i=1, 3
   lcharg(i,i)=-VV
END DO

Hamt=0.d0
DO i=1, N
   Hamt((i-1)*Psring+1:i*Psring,(i-1)*Psring+1:i*Psring)=Hcamd
END DO
DO i=1, N-1
   Hamt((i-1)*Psring+PsrHalf+1:i*Psring, &
        i*Psring+1:i*Psring+PsrHalf)=Hentcam
END DO
DO i=2,N
   Hamt((i-1)*Psring+1:(i-1)*Psring+PsrHalf, &
        (i-2)*Psring+PsrHalf+1:(i-1)*Psring)=Hentcam
END DO

!To plug the SOC and other term for impurities
DO i=1, N
   nimpesq=impureza(sample,i,1,1)
   nimpdir=impureza(sample,i,2,1)
   DO j=1, nimpdir
      l=(i-1)*Psring+impureza(sample,i,2,1+j)*2-1
      m=PsrHalf+l                                         !    ___  
      Hamt(l:l+2,l:l+2)=Hamt(l:l+2,l:l+2)+SOee+lcharg     !  /    \      .
      Hamt(m:m+2,m:m+2)=Hamt(m:m+2,m:m+2)+SOdd+lcharg     ! /      \     .
      Hamt(l:l+2,m:m+2)=Hamt(l:l+2,m:m+2)+SOed            ! \     /    
      Hamt(m:m+2,l:l+2)=Hamt(m:m+2,l:l+2)+SOde            ! l\___/m 
   END DO
   DO j=1, nimpesq
      l=(i-2)*Psring+PsrHalf+impureza(sample,i,1,1+j)*2
      m=l+PsrHalf
      Hamt(l:l+2,l:l+2)=Hamt(l:l+2,l:l+2)+SOee+lcharg
      Hamt(m:m+2,m:m+2)=Hamt(m:m+2,m:m+2)+SOdd+lcharg
      Hamt(l:l+2,m:m+2)=Hamt(l:l+2,m:m+2)+SOed
      Hamt(m:m+2,l:l+2)=Hamt(m:m+2,l:l+2)+SOde
   END DO
END DO

!Coloca uma vacancia

! k=(N/2+1)*Psring-PsrHalf
! l=(N/2+1)*Psring
! Hamt(k,k-1)=timp
! Hamt(k-1,k)=timp
! Hamt(l,l-1)=timp
! Hamt(l-1,l)=timp
! ! Hamt(k,l)=timp
! ! Hamt(l,k)=timp
! ! Hamt(l-2,k-2)=timp
! ! Hamt(k-2,l-2)=timp
! 
! l=k+1
! k=l-PsrHalf
! Hamt(k,k+1)=timp
! Hamt(k+1,k)=timp
! Hamt(l,l+1)=timp
! Hamt(l+1,l)=timp
! 
! Hamt(k+2,l+2)=timp
! Hamt(l+2,k+2)=timp
! Hamt(k+4,l+4)=timp
! Hamt(l+4,k+4)=timp

!Realiza uma constricao
IF(nlin.NE.0) THEN
   DO j=cami, camf
      k=(j-1)*Psring
      l=k+PsrHalf
      m=k+Psring
      DO i=1, nlin, 1
         IF(i==PsrHalf) THEN
            Hamt(k+i,l+i)=timp
            Hamt(l+i,k+i)=timp
            EXIT
         END IF
         IF(MOD(i,2)==1) THEN
            Hamt(k+i,l+i)=timp
            Hamt(k+i,l+i)=timp
            Hamt(k+i,k+i+1)=timp
            Hamt(k+i+1,k+i)=timp
            Hamt(l+i,l+i+1)=timp
            Hamt(l+i+1,l+i)=timp
            
            Hamt(l+1-i,m+1-i)=timp
            Hamt(m+1-i,l+1-i)=timp
            Hamt(l+1-i,l-i)=timp
            Hamt(l-i,l+1-i)=timp
            Hamt(m+1-i,m-i)=timp
            Hamt(m-i,m+1-i)=timp
         ELSE
            Hamt(k+i,k+i+1)=timp
            Hamt(k+i+1,k+i)=timp   
            Hamt(l+i,l+i+1)=timp
            Hamt(l+i+1,l+i)=timp   
            Hamt(k+i-PsrHalf,k+i)=timp
            Hamt(k+i,k+i-PsrHalf)=timp
            Hamt(m+i,l+i)=timp
            Hamt(l+i,m+i)=timp
            
            Hamt(l+1-i,l-i)=timp
            Hamt(l-i,l+1-i)=timp   
            Hamt(m+1-i,m-i)=timp
            Hamt(m-i,m+1-i)=timp
            Hamt(l+1-i,k+1-i)=timp
            Hamt(k+1-i,l+1-i)=timp
            Hamt(m+1-i,m+1-i+PsrHalf)=timp
            Hamt(m+1-i+PsrHalf,m+1-i)=timp
         END IF
      END DO
   END DO
END IF

!Realiza uma vacancia
IF(tipoimp=="Atop".AND.sitio.NE.0)THEN
   k=(N/2)*Psring-PsrHalf+sitio
!    Hamt(k,k)=-1.d9
!    GO TO 30
   IF(MOD(sitio,2)==0)THEN
      l=k+PsrHalf
   ELSE
      l=k-PsrHalf
   END IF
   ! print *, Psring, PsrHalf,k, l; stop

!    Hamt(k,k)=1.d10
   Hamt(k,k+1)=timp
   Hamt(k+1,k)=timp
   Hamt(k,l)=timp
   Hamt(l,k)=timp
   IF(sitio.NE.1) THEN
      Hamt(k,k-1)=timp
      Hamt(k-1,k)=timp
   END IF
! 30 CONTINUE
END IF

! !Realiza uma vacancia no lado esquerdo de uma camada
! IF(tipoimp=="Atop".AND.sitio.NE.0)THEN
!    k=(N/2)*Psring-PsrHalf+sitio
!    IF(MOD(sitio,2)==0)THEN
!       l=k+PsrHalf
!    ELSE
!       l=k-PsrHalf
!    END IF
!    ! print *, Psring, PsrHalf,k, l; stop
! 
!    Hamt(l,l+1)=timp
!    Hamt(l+1,l)=timp
!    Hamt(k,l)=timp
!    Hamt(l,k)=timp
!    IF(sitio.NE.1) THEN
!       Hamt(l,l-1)=timp
!       Hamt(l-1,l)=timp
!    END IF
! END IF

!End: Constructing the hamiltonian of the central part

z=CMPLX(realz,imgz,double)
CALL greensupfcontato(z,greencont,eps0)

tpc=0.d0
tcq=0.d0
DO i=1, P
   tpc(2*i,2*i)=t
   tcq(PsrN+1-2*i,2*i)=t
END DO
tcp=TRANSPOSE(tpc)
tqc=TRANSPOSE(tcq)

! AB=MATMUL(tcp,greencont)
CALL ZGEMM ( 'N', 'N', PsrN, PsrHalf, PsrHalf, one, tcp, PsrN, greencont, PsrHalf, zero, AB     , PsrN )
! SelfEsq=MATMUL(AB,tpc)
CALL ZGEMM ( 'N', 'N', PsrN, PsrN   , PsrHalf, one, AB , PsrN, tpc      , PsrHalf, zero, SelfEsq, PsrN )
! AB=MATMUL(tcq,greencont)
CALL ZGEMM ( 'N', 'N', PsrN, PsrHalf, PsrHalf, one, tcq, PsrN, greencont, PsrHalf, zero, AB     , PsrN )
! SelfDir=MATMUL(AB,tqc)
CALL ZGEMM ( 'N', 'N', PsrN, PsrN   , PsrHalf, one, AB , PsrN, tqc      , PsrHalf, zero, SelfDir, PsrN )

Hamt=-Hamt-SelfEsq-SelfDir
DO i=1, PsrN
   Hamt(i,i)=Hamt(i,i)+z
END DO

CALL printtime("START Hamt invertion")
CALL invers(Hamt,PsrN)
CALL printtime("END   Hamt invertion")

gn1=Hamt(PsrN-PsrHalf+1:PsrN,1:PsrHalf)
g1n=Hamt(1:PsrHalf,PsrN-PsrHalf+1:PsrN)

tt=0.d0
DO i=1, P !Antes tava: DO i=1, PsrHalf
   tt(2*i,2*i)=t
END DO

! AC=MATMUL(tt,greencont)
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, tt, PsrHalf, greencont, PsrHalf, zero, AC   , PsrHalf )
! Sigma=MATMUL(AC,tt)
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, AC, PsrHalf, tt       , PsrHalf, zero, Sigma, PsrHalf )
SigmaDag = CONJG(TRANSPOSE(Sigma)) !Sigma+
Gamma=iimag*(Sigma-SigmaDag) !i(Sigma - Sigma+)

CALL printtime("END   Green fc calc.")
RETURN


!==========================================================================
WRITE(formato, FMT="(a,i0,a)")  '( ',PsrHalf, '(f4.1) )'

DO j=1, PsrN
   WRITE( UNIT=*, FMT=formato) REAL(tcp(j,:))
END DO
!==========================================================================
STOP
END SUBROUTINE green