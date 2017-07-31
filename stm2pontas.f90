! Programa para calcular um mapa de condutancia entre duas pontas de STM,
! uma fixa e outra movel. O mapa eh feito sobre uma fita de grafeno armchair.
! Autor: Flavinao Jose dos Santos
! Forschungszentrum Juelich
! Maio de 2014

!Como plotar os gráficos no gnuplot
!pontos:
! plot "./mesh.dat" u 1:2:3 notitle  lc palette lw 1 pt 5 , "sitio.dat" notitle pt 7 ps 0.3, "pontafix.dat" notitle pt 7 lw 5 lc 6

!mapa continuo com interpolacao:
! set pm3d map
! set pm3d interpolate 0,0
! splot "./matriz.dat" notitle

INCLUDE "./parametros.f90"
INCLUDE "./modulos.f90"
INCLUDE "./invers.f90"
INCLUDE "./acRBsinf.f90"
INCLUDE "./greensuperficie.f90"
INCLUDE "./greenfcSampleImpHexSOLBComp.f90"
INCLUDE "impureza.f90"

PROGRAM stm2pontas

USE parametros
IMPLICIT NONE
COMPLEX(KIND=8), PARAMETER :: imag=(0.d0,1.d0)
INTEGER, PARAMETER :: Psring=6+(P-1)*4
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2

INTEGER :: i, ii, k, kk, l1, j1, l2, j2, iL, iR, nL, nR, iL0, iLf, mvar, &
           u, v, nn, m, contpA, contpB, contpAcontpB, Ij, Ji, Jj, contig,&
           cont, cont1, contpAmcontig, qq=0, check, check2
INTEGER :: impureza(namostras,N,2,P+1), impurezaux(N,2,P+1)
INTEGER, ALLOCATABLE :: coordsitvizpA(:), coordsitvizpB(:), &
                        coordaux(:), sitcomptlh(:,:)

COMPLEX(KIND=8),  DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1
COMPLEX(KIND=8), DIMENSION(PsrN,PsrN) :: Hamt, BBB
COMPLEX(KIND=8) :: gBB
COMPLEX(KIND=8), ALLOCATABLE :: SigmaA(:,:), GammaB(:,:), GammaA(:,:), &
                                g(:,:), GG(:,:), SigmaB(:,:), GIJ(:,:), &
                                A(:,:), B(:,:), C(:,:)

REAL(KIND=8) ::x, y, L, R, dL, dR, spin=-1.d0, wi, condt, distpApB, &
               valor1, valor2
REAL(KIND=8) :: vetor(3,2), vetor2(3,2), base1(2), base2(2), posicao(2), &
                vtcorr(2), vcorner(2), &
                vpospA(2), vposatovizpA(2), vpospB(2), vaux(2)
REAL(KIND=8), ALLOCATABLE :: tIA(:,:), tAI(:,:), tJB(:,:), tBJ(:,:), &
                dvizpB(:), dvizpA(:), dvizaux(:), matriz(:,:)
CHARACTER(LEN=130) :: nome

PRINT *, "Started."
CALL sortimpureza(impureza)
impurezaux=impureza(1,:,:,:)

! impurezaux=0
! impurezaux(N/2,1,1)=1
! impurezaux(N/2,1,2)=4

WRITE(nome, FMT="(a,i0,a,i0,a,f4.1,a,f4.1,a,f4.1,a,f4.1,a,f4.1,a,i0,a,f4.1,a)") &
   "~/output/fgreenCompAcRBinfP",P,"N",N,"t",t,"ctcimpz",concentracao,&
   "dmu",fractdmu,"|t|lmda",fractlmda,"|t|timp",fracttimp, &
   "t_namost",namostras,"eini",eini,".dat"

 check=0
 check2=0
OPEN(UNIT=21, FILE=nome, STATUS="OLD", ACTION="READ", IOSTAT=check)
 check=1
IF(check/=0) THEN
   CALL green(Hamt,Gamma,g1n,gn1,impurezaux,ee,eta,spin)
   PRINT *, "Green Function calculated."
ELSE
   DO i=1, PsrN
      DO m=1, PsrN
         READ(UNIT=21, FMT="(D25.17,D25.17)",IOSTAT=check2) valor1, valor2
         Hamt(i,m)=valor1+imag*valor2
      END DO
   END DO
   PRINT *, "Green Function read."
END IF

gBB=( (ee+imag*eta-eps0)+SQRT((ee+imag*eta-eps0)**2-4.d0*t**2) )/(2*t**2)
IF(-AIMAG(gBB) < 0.d0) THEN
   gBB=( (ee+imag*eta-eps0)-SQRT((ee+imag*eta-eps0)**2-4.d0*t**2) )/(2*t**2)
END IF

vetor(1,:)=dCC*(/1,0/)
vetor(2,:)=dCC*(/-0.5d0, SQRT(3.d0)/2.d0/)
vetor(3,:)=dCC*(/-0.5d0,-SQRT(3.d0)/2.d0/)

base1=vetor(1,:)-vetor(3,:)
base2=vetor(1,:)-vetor(2,:)

! _____L_____
!|___*___*___| R

vaux=base1+base2
L=3.d0*dCC
R=SQRT( (base1(1))**2 + (base1(2))**2 ) / 2
nL=NINT(SQRT(npmesh/(R/L)))
IF(MOD(nL,2)==1) nL=nL+1
nR=NINT(SQRT(npmesh*(R/L)))
dL=L/REAL(nL)
dR=R/REAL(nR)

vcorner=(/ -dCC, -R/2 /)
mvar=NINT(raio/dCC)
IF(mvar==0) mvar=1
valor1=L*mvar+L/2.d0+dCC/2.d0
valor2=R/2
i=2*(2*mvar+1)**2
ALLOCATE( coordsitvizpA(i), coordsitvizpB(i), coordaux(i), dvizpB(i), dvizpA(i) )
ALLOCATE( dvizaux(i), matriz(NINT(nL*(REAL(N-2*mvar)-0.5)),nR*PsrHalf) )
 coordsitvizpA=0
 coordsitvizpB=0
 coordaux=0
 dvizpA=0.d0
 dvizpB=0.d0
 dvizaux=0.d0

OPEN(UNIT=10, FILE="mesh.dat", STATUS="UNKNOWN")
OPEN(UNIT=11, FILE="sitio.dat", STATUS="UNKNOWN")
OPEN(UNIT=12, FILE="pontafix.dat", STATUS="UNKNOWN")

!Determinando os vizinhos da ponta fixa
vpospB=lpB*base1+jpB*base2-centroHex*vetor(1,:)
WRITE(UNIT=12, FMT="(2(f14.10,';'), f14.10)")&
                           vpospB, 0.d0
CLOSE(UNIT=12)

 contpB = 0
DO m=-mvar,mvar
   DO nn=-mvar,mvar
      l2 = lpB+m
      j2 = jpB+nn
         
      IF(j2<l2 .AND. (l2-j2)<2*(P+1)) THEN
         vposatovizpA = l2*base1+j2*base2
         kk = vposatovizpA(1)/(3.001d0*dCC) + 1
         ii = vposatovizpA(2)/(vetor(2,2)+0.001d0) + 1
         u = Psring*(kk-1) + ii*MOD(ii,2) + (PsrHalf+ii)*MOD(ii+1,2)
         v = u + PsrHalf

         vaux=vposatovizpA-vpospB
         IF(SQRT(vaux(1)**2+vaux(2)**2)<=raio) THEN
            contpB=contpB+1
            coordsitvizpB(contpB)=u
            dvizpB(contpB)=SQRT(vaux(1)**2+vaux(2)**2)
         END IF
         vaux=(vposatovizpA+vetor(1,:))-vpospB
         IF(SQRT(vaux(1)**2+vaux(2)**2)<=raio) THEN
            contpB=contpB+1
            coordsitvizpB(contpB)=v
            dvizpB(contpB)=SQRT(vaux(1)**2+vaux(2)**2)
         END IF
      END IF
   END DO
END DO

ALLOCATE( tJB(contpB,1), tBJ(1,contpB), SigmaB(contpB,contpB) )
ALLOCATE( GammaB(contpB,contpB) )

wi=0.d0
DO i=1, contpB
   wi=wi+EXP(-aa*(dvizpB(i)**2))
END DO
DO i=1, contpB
   tJB(i,1)=10.d0*t*EXP( -aa*dvizpB(i)**2 - dvizpB(i)/llamb) / wi
   tBJ(1,i)=tJB(i,1)
END DO

SigmaB=MATMUL(tJB*gBB,tBJ)
GammaB=imag*(SigmaB-CONJG(TRANSPOSE(SigmaB)))

DO k=1+mvar, N-mvar
   l1=k
   j1=k
   DO i=1, PsrHalf
      l1=l1+MOD(i+1,2)
      j1=j1-MOD(i,2)

      WRITE(UNIT=11, FMT="(1(f14.10,';'), f14.10)")l1*base1+j1*base2
      WRITE(UNIT=11, FMT="(1(f14.10,';'), f14.10)")l1*base1+j1*base2+vetor(1,:)

      iL0=1
      iLf=nL
      IF(k==1+mvar.AND.MOD(i,2)==1) iL0=nL/2 + 1
      IF(k==(N-mvar).AND.MOD(i+1,2)==1) iLf=nL/2

      DO iL=iL0, iLf
         DO iR=1, nR
            x=REAL(iL)*dL-dL/2.d0
            y=REAL(iR)*dR-dR/2.d0

            vpospA=l1*base1+j1*base2+vcorner+(/x,y/)
            vaux=vpospA-vpospB
            distpApB=SQRT(vaux(1)**2+vaux(2)**2)

            contpA=0
            coordsitvizpA=0.d0
            DO m=-mvar,mvar
               DO nn=-mvar,mvar
                  l2 = l1+m
                  j2 = j1+nn

                  IF(j2<l2 .AND. (l2-j2)<2*(P+1)) THEN
                     vposatovizpA = l2*base1+j2*base2
                     kk = vposatovizpA(1)/(3.001d0*dCC) + 1
                     ii = vposatovizpA(2)/(vetor(2,2)+0.001d0) + 1
                     u = Psring*(kk-1) + ii*MOD(ii,2) + (PsrHalf+ii)*MOD(ii+1,2)
                     v = u + PsrHalf

                     vaux=vposatovizpA-vpospA
                     IF(SQRT(vaux(1)**2+vaux(2)**2)<=raio) THEN
                        contpA=contpA+1
                        coordsitvizpA(contpA)=u
                        dvizpA(contpA)=SQRT(vaux(1)**2+vaux(2)**2)
                     END IF
                     vaux=(vposatovizpA+vetor(1,:))-vpospA
                     IF(SQRT(vaux(1)**2+vaux(2)**2)<=raio) THEN
                        contpA=contpA+1
                        coordsitvizpA(contpA)=v
                        dvizpA(contpA)=SQRT(vaux(1)**2+vaux(2)**2)
                     END IF
                  END IF
               END DO
            END DO

            IF(contpA==0) GOTO 10

            ALLOCATE( sitcomptlh(contpB,2) )
            contig=0
            DO ii=1, contpB
               DO kk=1, contpA
                  IF(coordsitvizpB(ii)==coordsitvizpA(kk)) THEN
                     contig=contig+1
                     sitcomptlh(contig,:)=(/ii,kk/)
                  END IF
               END DO
            END DO
            contpAmcontig=contpA-contig

            IF(contig/=0) THEN
               coordaux=coordsitvizpA
               dvizaux=dvizpA
               cont1=0
               DO ii=1, contpA
                  cont=0
                  DO kk=1, contig
                     IF(ii==sitcomptlh(kk,2)) cont=1
                  END DO
                  IF(cont==0) THEN
                     cont1=cont1+1
                     coordsitvizpA(cont1)=coordaux(ii)
                     dvizpA(cont1)=dvizaux(ii)
                  END IF
               END DO
               DO ii=1, contig
                  coordsitvizpA(contpAmcontig+ii)=coordaux(sitcomptlh(ii,2))
                  dvizpA(contpAmcontig+ii)=dvizaux(sitcomptlh(ii,2))
               END DO
               coordaux=coordsitvizpB
               dvizaux=dvizpB
               cont1=0
               DO ii=1, contpB
                  cont=0
                  DO kk=1, contig
                     IF(ii==sitcomptlh(kk,1)) cont=1
                  END DO
                  IF(cont==0) THEN
                     cont1=cont1+1
                     coordsitvizpB(contig+cont1)=coordaux(ii)
                     dvizpB(contig+cont1)=dvizaux(ii)
                  END IF
               END DO
               DO ii=1, contig
                  coordsitvizpB(ii)=coordaux(sitcomptlh(ii,1))
                  dvizpB(ii)=dvizaux(sitcomptlh(ii,1))
               END DO
               
               wi=0.d0
               DO m=1, contpB
                  wi=wi+EXP(-aa*(dvizpB(m)**2))
               END DO
               DO m=1, contpB
                  tJB(m,1)=10.d0*t*EXP( -aa*dvizpB(m)**2 - dvizpB(m)/llamb) / wi
                  tBJ(1,m)=tJB(m,1)
               END DO

               SigmaB=MATMUL(tJB*gBB,tBJ)
               GammaB=imag*(SigmaB-CONJG(TRANSPOSE(SigmaB)))
            END IF

            contpAcontpB=contpAmcontig+contpB
            ALLOCATE( g(contpAcontpB,contpAcontpB), GG(contpAcontpB,contpAcontpB) )
            ALLOCATE( GIJ(contpA,contpB), A(contpA,contpB), B(contpB,contpA), C(contpA, contpA) )

            g=0.d0
            DO Ii=1, contpAmcontig
               DO Ij=1, contpAmcontig
                  g(Ii,Ij)=Hamt(coordsitvizpA(Ii),coordsitvizpA(Ij))
               END DO
               DO Jj=1, contpB
                  g(Ii,contpAmcontig+Jj)=Hamt(coordsitvizpA(Ii),coordsitvizpB(Jj))
               END DO
            END DO
            DO Ji=1, contpB
               DO Ij=1, contpAmcontig
                  g(contpAmcontig+Ji,Ij)=Hamt(coordsitvizpB(Ji),coordsitvizpA(Ij))
               END DO
               DO Jj=1, contpB
                  g(contpAmcontig+Ji,contpAmcontig+Jj)=Hamt(coordsitvizpB(Ji),coordsitvizpB(Jj))
               END DO
            END DO

            ALLOCATE( tIA(contpA,1), tAI(1,contpA), SigmaA(contpA,contpA) )
            ALLOCATE( GammaA(contpA,contpA) )
            wi=0.d0
            DO m=1, contpA
               wi=wi+EXP(-aa*(dvizpA(m)**2))
            END DO
            DO m=1, contpA
               tIA(m,1)=10.d0*t*EXP( -aa*dvizpA(m)**2 - dvizpA(m)/llamb) / wi
               tAI(1,m)=tIA(m,1)
            END DO

            SigmaA=MATMUL(tIA*gBB,tAI)
            GammaA=imag*(SigmaA-CONJG(TRANSPOSE(SigmaA)))

            GG=0.d0
            CALL invers(g,contpAcontpB)
            GG(contpAmcontig+1:contpAcontpB,contpAmcontig+1:contpAcontpB)=-SigmaB
            GG(1:contpA,1:contpA)=GG(1:contpA,1:contpA)-SigmaA
            GG=GG+g
            CALL invers(GG,contpAcontpB)
            GIJ=GG(1:contpA,contpAmcontig+1:contpAcontpB)

            A=MATMUL(GammaA,GIJ)
            B=MATMUL(GammaB,CONJG(TRANSPOSE(GIJ)))
            C=MATMUL(A,B)

            condt=0.d0
            DO m=1,contpA
               condt=condt+C(m,m)
            END DO
            condt=condt*EXP(fatexp*distpApB)!*(distpApB)!EXP(0.19d0*distpApB)
            10 CONTINUE
            IF(contpA==0) condt=0.d0
            
            WRITE(UNIT=10, FMT="(2(f14.10,';'), f16.10)") vpospA, condt
            matriz(INT((vpospA(1)-valor1)/dL)+1,INT((vpospA(2)-valor2)/dR)+1)=condt
            
            IF(contpA/=0) THEN
               DEALLOCATE(g, GG, GIJ, tIA, tAI, GammaA, SigmaA, A, B, C, sitcomptlh)
            END IF
         END DO
      END DO
   END DO
END DO

OPEN(UNIT=14, FILE="matriz.dat", STATUS="UNKNOWN")
DO ii=1, NINT(nL*(REAL(N-2*mvar)-0.5))
   DO kk=1, nR*PsrHalf
      WRITE(UNIT=14, FMT="(i0,' ',i0,' ',f16.10)")ii, kk, matriz(ii,kk)
   END DO
   WRITE(UNIT=14, FMT=*)
END DO
CLOSE(UNIT=14)

 check=0
IF(check/=0) THEN
   check2=0
   PRINT *, "Saving Green Function..."
   OPEN(UNIT=21, FILE=nome, STATUS="NEW", ACTION="WRITE", IOSTAT=check2)
   IF(check2/=0) THEN
      Print*, "Ops! Erro na criacao do arquivo."; stop
   END IF
   DO i=1, PsrN
      DO m=1, PsrN
         WRITE(UNIT=21, FMT="(D25.17,D25.17)") Hamt(i,m)
      END DO
   END DO
   PRINT *, "Saved."
END IF
CLOSE(UNIT=21)

END PROGRAM stm2pontas