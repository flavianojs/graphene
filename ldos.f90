INCLUDE "parametros.f90"
INCLUDE "./modulos.f90"
INCLUDE "./invers.f90"
INCLUDE "./acRBsinf.f90"
INCLUDE "./acRBinf.f90"
INCLUDE "./greensuperficie.f90"
INCLUDE "greenfcSampleImpHexSOLBComp.f90"
INCLUDE "impureza.f90"

PROGRAM LDOS
USE parametros
IMPLICIT NONE
! USE constants
! USE propagator
!====PARAMETROS===================
INTEGER :: npont=50, l, i!,   l1=2 ,j1=-1 ,s1=1,   l2, j2, s2

!=================================
INTEGER, PARAMETER :: Psring=4*P+2 - 2*assm !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
INTEGER :: impureza(namostras,N,2,P+1), impurezaux(N,2,P+1)
REAL(KIND=8) :: incremento, E
REAL(KIND=8) :: menusone_pi, tpi
COMPLEX(KIND=8) :: y
COMPLEX(KIND=8), DIMENSION(PsrN,PsrN):: Gr
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1

CHARACTER(LEN=300) :: nomeArq
CHARACTER(LEN=1) :: tipo='A'

! l2=l1
! j2=j1
! s2=s1
WRITE(nomeArq,FMT="(a,i0,a,ES7.0E2,a,ES7.1E1,a,i0,a)") &
                              & "./dados/ldosBordaArmSO P",P," eta",eta, " t ",&
                                    t," sitio",sitioviz,".dat"
OPEN(UNIT=19, FILE=nomeArq, STATUS="UNKNOWN", ACTION="WRITE")
! WRITE(UNIT=19, FMT=*)"#tipo  N    l1 j1 s1 l2 j2 s2   npont   ef     eta"
! WRITE(UNIT=19, FMT="(a,a,4x,i4.3,2x,i0,2x,i0,2x,i0,2x,i0,2x,i0,2x,i0,   &
!                             &3x,i4,3x,es6.0e1,2x,es7.1e1,2x,es7.0e2)")  &
!                    " #", tipo, cl, l1, j1, s2, j1, l2, s2, npont, ef, eta

CALL sortimpureza(impureza)
impurezaux=impureza(1,:,:,:)
               
tpi=2.d0*pi
menusone_pi=-1.d0/pi
incremento= (efin-eini)/REAL(npt-1)
! y=eta
l=(N/2)*Psring-PsrHalf+sitioviz
! l=((N-1)/2)*Psring+sitio+sitioviz

DO i=0, npt-1
   E=eini+incremento*REAL(i)
!    CALL acRBsinf(E,y,l1,j1,s1,l2,j2,s2,green)
   CALL green(Gr,Gamma,g1n,gn1,impurezaux,E, eta, 1.d0)
   y=Gr(l,l)
   CALL green(Gr,Gamma,g1n,gn1,impurezaux,E, eta,-1.d0)
   y=y+Gr(l,l)
   WRITE(UNIT=19, FMT=*) E, menusone_pi*AIMAG(y)
END DO
CLOSE(UNIT=19)
END PROGRAM LDOS