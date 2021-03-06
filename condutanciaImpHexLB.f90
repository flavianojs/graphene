INCLUDE "./parametros.f90"

!arquivos necessarios para as funcoes de green
INCLUDE "./modulos.f90"
INCLUDE "./invers.f90"
INCLUDE "./acRBsinf.f90"
INCLUDE "./sqRBsinf.f90"
INCLUDE "./greensuperficie.f90"
INCLUDE "./greenfcSampleImpHexSOLBComp.f90"
INCLUDE "./greenfcSampleImpHexSOLB.f90"
INCLUDE "./greenfcSampleImpHexSOLBaddlay.f90"

SUBROUTINE condutancia(E, condt, sample, spin)
USE parametros
IMPLICIT NONE
COMPLEX(KIND=8), PARAMETER :: iimag=(0.d0,1.d0),zero=0.d0,one=1.d0, m_one=-1.d0
INTEGER, PARAMETER :: Psring=4*P+2, cmda=0 !cmda deve ser <= N
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
! INTEGER, INTENT(IN) :: impureza(N,2,P+1)
! REAL(KIND=8), INTENT(IN) :: anderdisor(N,2,PsrHalf)
REAL(KIND=8), INTENT(IN) :: E, spin
INTEGER, INTENT(IN) :: sample
REAL(KIND=8), INTENT(OUT) :: condt
REAL(KIND=8) :: y=eta
INTEGER :: i
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: Gamma,GN1,G1N,GN1dag,A,B,C
COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hamt(:,:)
IF(flag==1) THEN
   CALL green1(Gamma,G1N,GN1,E,y,sample,spin) !Add layer
ELSE
   IF(flag==2) THEN
   	  ALLOCATE( Hamt(PsrN,PsrN) )
      CALL green(Hamt,Gamma,G1N,GN1,E,y,sample,spin) !Inv. Compl.
   ELSE
      CALL green2(Gamma,G1N,GN1,E,y,sample,spin) !Add lay extreme
   END IF
END IF

! GN1dag=CONJG(TRANSPOSE(GN1))

! A=MATMUL(Gamma,GN1) !ATENCAO estou testanto GN1
! B=MATMUL(Gamma,GN1dag)
!  C=MATMUL(A,B)
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, Gamma, PsrHalf, GN1, PsrHalf, zero, A, PsrHalf ) 
CALL ZGEMM ( 'N', 'C', PsrHalf, PsrHalf, PsrHalf, one, Gamma, PsrHalf, GN1, PsrHalf, zero, B, PsrHalf ) 
CALL ZGEMM ( 'N', 'N', PsrHalf, PsrHalf, PsrHalf, one, A    , PsrHalf, B  , PsrHalf, zero, C, PsrHalf ) 

 condt=0.d0
DO i=1, PsrHalf
   condt=condt+C(i,i)
END DO

 condt=condt
RETURN
END SUBROUTINE condutancia

    