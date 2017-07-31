INCLUDE "parametros.f90"
INCLUDE "global.f90"
INCLUDE "./modulos.f90"
INCLUDE "./invers.f90"
INCLUDE "./acRBsinf.f90"
INCLUDE "./sqRBsinf.f90"
INCLUDE "./greensuperficie.f90"
INCLUDE "printtime.f90"
INCLUDE "greenfcSampleImpHexSOLBComp.f90"
! INCLUDE "greenfcSampleCorrente.f90"
INCLUDE "impureza.f90"

!para o gnuplot
! plot "./dados/testeImpHexSO.dat" u 1:2:(1.*$3):(1.*$4):5 notitle w vec lc palette lw 2
! set palette negative

PROGRAM corrente
USE parametros
IMPLICIT NONE
INTEGER, PARAMETER :: Psring=4*P+2 !N. de sitio num anel
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
COMPLEX(KIND=8), PARAMETER :: iimag=(0.d0,1.d0)
COMPLEX(KIND=8), DIMENSION(PsrHalf,PsrHalf) :: Gamma,g1n,gn1
INTEGER :: i,j,l,k,l1,j1,s1,l2,j2,s2,u,v,f,e, flaglc2
! , flaglc=3, & !Se 1, plota apenas setas
!            lado=1!lado=1->esq, lado=0->dir                     !horiz, Se 3, plota todas

! INTEGER :: impureza(namostras,N,2,P+1), impurezaux(N,2,P+1), &
INTEGER :: impurezaux(N,2,P+1), &
           laRaux, laIaux, timaux, VVaux, rhoaux
REAL(KIND=8) :: vetor(3,2), vetor2(3,2), base1(2), base2(2), posicao(2), &
                corinho=0.d0, vtcorr(2)
REAL(KIND=8) :: centrodebanda(Psring,N), icorr(3), maior, txfour=2.d0*t, corrtotal(N)
COMPLEX(KIND=8), DIMENSION(PsrN,PsrN):: Gr, Ga, Gn, Gaux, Gama, Gamar, Gamaa, Gaux2
CHARACTER(LEN=120) :: nome
REAL ::  beg_cpu_time, end_cpu_time, timedif

CALL cpu_time (beg_cpu_time)

flaglc2=2
IF(flaglc==1) flaglc2=1

i=0
DO k=1, 80
   IF(tipoimp(k:k).NE." ") THEN
      i=i+1
   ELSE
      EXIT
   END IF
END DO

laRaux=0; laIaux=0; timaux=0; VVaux=0; rhoaux=0
IF(ABS(REAL(lambda)).NE.0.d0) laRaux=1
IF(ABS(IMAG(lambda)).NE.0.d0) laIaux=1
IF(ABS(aatimp)      .NE.0.d0) timaux=1
IF(ABS(VV )         .NE.0.d0)  VVaux=1
IF(ABS(rho)         .NE.0.d0) rhoaux=1
WRITE(nome, FMT="(a,i0,a,i0,a,f4.1,a,i0,a,i0,a,i0,a,i0,a,i0,a,f4.1,a)") &
   "./dados/distcurrP",P,"N",N,"ctc",concentracao,&
   "laR", laRaux, "laI", laIaux, &
   "tim", timaux, "rho", rhoaux, &
   "VV" ,  VVaux, "E"  ,     ee, ".dat"

WRITE(nome, FMT="(a,a,a,a)") &
            "./dados/distcurr", tipoimp(1:i),term,".dat"
OPEN(UNIT=10, FILE=nome, STATUS="UNKNOWN", ACTION="WRITE")

WRITE(nome, FMT="(a)") &
            "./dados/lattice.dat"
OPEN(UNIT=11, FILE=nome, STATUS="UNKNOWN", ACTION="WRITE")

! CALL sortimpureza(impureza)
CALL sortimpureza()
impurezaux=impureza(1,:,:,:)

!To remove all imp. on the left/right to compute the total cross-section current with precision when SOC is on
IF(remove_imp_ONleftORright) THEN
   DO i=1, N
      impurezaux(i,lado,1)=0
   END DO
END IF

CALL green(Gr,Gamma,g1n,gn1,ee,eta,1, 1.d0)
Ga=CONJG(TRANSPOSE(Gr))
Gama=0.d0
Gama(1:PsrHalf,1:PsrHalf)=Gamma
!Gama(PsrN-PsrHalf+1:PsrN,PsrN-PsrHalf+1:PsrN)=Gamma

Gaux=MATMUL(Gr,Gama)
Gn=MATMUL(Gaux,Ga)

CALL green(Gr,Gamma,g1n,gn1,ee,eta,1,-1.d0)
Ga=CONJG(TRANSPOSE(Gr))

Gaux=MATMUL(Gr,Gama)
Gaux2=MATMUL(Gaux,Ga)
Gn=Gn+Gaux2

vetor(1,:)=(/1.d0,0.d0/)
vetor(2,:)=(/-0.5d0, SQRT(3.d0)/2.d0/)
vetor(3,:)=(/-0.5d0,-SQRT(3.d0)/2.d0/)

base1=vetor(1,:)-vetor(3,:)
base2=vetor(1,:)-vetor(2,:)

!Vamos percorrer apenas os sítios pretos, em cada sítio preto, pegamos seus
!tres vizinhos. Onde houve alteracao nesse codigo devido a nova maneira de
!os sitos na matriz da funcao de green, assinalei com '!Mudei aqui'
maior=0.d0
 corrtotal=0.d0
DO k=1, N
   l1=k
   j1=k
   DO i=1, PsrHalf
      l1=l1+MOD(i+1,2)
      j1=j1-MOD(i,2)
      posicao=l1*base1+j1*base2
      !Posicao do sitio preto na matriz da funcao de green
      u = Psring*(k-1) + i*MOD(i,2) + (PsrHalf+i)*MOD(i+1,2) !Mudei aqui
      !Posicao do seu vizinho a direita, que fica na proxima linha, mas na
      !mesma coluna.
      v = Psring*(k-1)+PsrHalf + i*MOD(i,2) + (i+PsrHalf)*MOD(i+1,2) !Mudei aqui

      !Se for o primeiro sitio da camada
      IF(i==1) THEN
!          IF(sitio==1) THEN
!             !Se sitio eh impar, vemos trocar u por v e vice-vesa
!             icorr(1)=2.d0*timp*DIMAG( Gn(v,u  ) )
!             icorr(2)=2.d0*timp*DIMAG( Gn(v,v+1) )
!          END IF

         icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
         icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) )
         IF(lado.EQ.1) THEN
            corrtotal(k)=corrtotal(k)+icorr(1)
         END IF
         DO j=1,2
            IF(ABS(icorr(j))>=maior ) THEN
               maior = ABS(icorr(j))
            END IF
         END DO
      ELSE
         !Se for o ultimo sitio da camada
         IF(i==PsrHalf) THEN
            icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
            icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) )
            IF(lado.EQ.1) THEN
               corrtotal(k)=corrtotal(k)+icorr(1)
            END IF
            DO j=1,3,2
               IF(ABS(icorr(j))>=maior ) THEN
                  maior = ABS(icorr(j))
               END IF
            END DO
         ELSE
            !Se for a ultima camada e sitio par
            IF(k==N.AND.MOD(i,2)==0) THEN
               icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) ) !Mudei aqui
               icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) ) !Mudei aqui
               DO j=2,3
                  IF(ABS(icorr(j))>=maior ) THEN
                     maior = ABS(icorr(j))
                  END IF
               END DO
            ELSE
               icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
               icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) ) !Mudei aqui
               icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) ) !Mudei aqui

               !Para diminuir o hopping ao redor de vacancias
               IF(k==(N/2)) THEN
                  IF(MOD(sitio,2)==0)THEN
                     IF(i==sitio) THEN
                        icorr(1)=2.d0*timp*DIMAG( Gn(u,v  ) )
                        icorr(2)=2.d0*timp*DIMAG( Gn(u,u+1) )
                        icorr(3)=2.d0*timp*DIMAG( Gn(u,u-1) )
                     END IF
                  ELSE
                     IF(i==sitio)THEN
                        icorr(1)=2.d0*timp*DIMAG( Gn(u,v  ) ) 
                     END IF
                     IF(i==(sitio+1))THEN
                        icorr(3)=2.d0*timp*DIMAG( Gn(u,u-1) )
                     END IF
                     IF(i==(sitio-1))THEN
                        icorr(2)=2.d0*timp*DIMAG( Gn(u,u+1) )
                     END IF
                  END IF
               END IF
               
               IF(lado.EQ.MOD(i,2)) THEN
                  corrtotal(k)=corrtotal(k)+icorr(1)
               END IF
               DO j=1,3
                  IF(ABS(icorr(j))>=maior ) THEN
                     maior = ABS(icorr(j))
                  END IF
               END DO
            END IF
         END IF
      END IF
   END DO
END DO
print *, maior
DO i=1, N
   print *, corrtotal(i)
END DO
maior=maior*1.1d0

DO k=1, N
   l1=k
   j1=k

   DO i=1, PsrHalf
      l1=l1+MOD(i+1,2)
      j1=j1-MOD(i,2)
      IF(flaglc==1 .AND. MOD(i,2)==0) CYCLE

      posicao=l1*base1+j1*base2
     ! WRITE(UNIT=11, FMT="(f9.5, ';', f9.5)")posicao
      u = Psring*(k-1) + i*MOD(i,2) + (PsrHalf+i)*MOD(i+1,2) !Mudei aqui
      v = Psring*(k-1)+PsrHalf + i*MOD(i,2) + (i+PsrHalf)*MOD(i+1,2) !Mudei aqui

      IF(i==1) THEN
         icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
         icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) )
         DO j=1, flaglc2!2
            vtcorr=icorr(j)*vetor(j,:)/maior
            IF( ABS(icorr(j))>=(maior/fator) * fator2 ) THEN
               WRITE(UNIT=10, FMT="(4(f9.5,';'),f9.5)") &
                                    posicao+vetor(j,:)/2.d0-vtcorr/2.d0, &
                                    vtcorr, ABS(icorr(j))
            END IF
            IF(nlin==0) THEN
               WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao, vetor(j,:)
            ELSE
               IF( .NOT.(cami-1<k .AND. k<=camf) ) WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao, vetor(j,:)
            END IF
         END DO
      ELSE
         IF(i==PsrHalf) THEN
            icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
            icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) )
            DO j=1, flaglc, 2!3, 2
               vtcorr=icorr(j)*vetor(j,:)/maior
               IF( ABS(icorr(j))>=(maior/fator) * fator2 ) THEN
                  WRITE(UNIT=10, FMT="(4(f9.5,';'),f9.5)") &
                                       posicao+vetor(j,:)/2.d0-vtcorr/2.d0, &
                                       vtcorr,ABS(icorr(j))
               END IF

               !To create the constriction
               IF(nlin==0) THEN
                  WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao, vetor(j,:)
               ELSE
                  IF( .NOT.(cami-1<k .AND. k<=camf) ) WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao, vetor(j,:)
               END IF

            END DO
         ELSE
            IF(k==N.AND.MOD(i,2)==0) THEN
               icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) ) !Mudei aqui
               icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) ) !Mudei aqui
               DO j=2, flaglc!3
                  vtcorr=icorr(j)*vetor(j,:)/maior
                  IF( ABS(icorr(j))>=(maior/fator) * fator2 ) THEN
                     WRITE(UNIT=10, FMT="(4(f9.5,';'),f9.5)") &
                                          posicao+vetor(j,:)/2.d0-vtcorr/2.d0, &
                                          vtcorr,ABS(icorr(j))
                  END IF
                  WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao,vetor(j,:)
               END DO
            ELSE
               icorr(1)=corinho+txfour*DIMAG( Gn(u,v  ) )
               icorr(2)=corinho+txfour*DIMAG( Gn(u,u+1) ) !Mudei aqui
               icorr(3)=corinho+txfour*DIMAG( Gn(u,u-1) ) !Mudei aqui
               
               !Para diminuir o hopping ao redor de vacancias
               IF(k==(N/2)) THEN
                  IF(MOD(sitio,2)==0)THEN
                     IF(i==sitio) THEN
                        icorr(1)=2.d0*timp*DIMAG( Gn(u,v  ) )
                        icorr(2)=2.d0*timp*DIMAG( Gn(u,u+1) )
                        icorr(3)=2.d0*timp*DIMAG( Gn(u,u-1) )
                     END IF
                  ELSE
                     IF(i==(sitio-1))THEN
                        icorr(2)=2.d0*timp*DIMAG( Gn(u,u+1) )
                     END IF
                     IF(i==sitio)THEN
                        icorr(1)=2.d0*timp*DIMAG( Gn(u,v  ) )
                     END IF
                     IF(i==(sitio+1))THEN
                        icorr(3)=2.d0*timp*DIMAG( Gn(u,u-1) )
                     END IF
                  END IF
               END IF
!                IF(k==(N/2+1).AND.i==sitio) THEN
!                   icorr(1)=2.d0*timp*DIMAG( Gn(u,v  ) )
!                   icorr(2)=2.d0*timp*DIMAG( Gn(u,u+1) )
!                   icorr(3)=2.d0*timp*DIMAG( Gn(u,u-1) )
!                END IF
               
               DO j=1, flaglc!3
                  vtcorr=icorr(j)*vetor(j,:)/maior
                  IF( ABS(icorr(j))>=(maior/fator) * fator2 ) THEN
                     WRITE(UNIT=10, FMT="(4(f9.5,';'),f9.5)") &
                                          posicao+vetor(j,:)/2.d0-vtcorr/2.d0, &
                                          vtcorr,ABS(icorr(j))
                  END IF

                  !To create the constriction
                  IF(nlin.NE.0 .AND. (cami-1<k .AND. k<=camf) .AND. ( i <= nlin .OR. i >PsrHalf-nlin ) ) THEN
                  ELSE
                     IF(i==nlin+1 .AND. j==3 .AND. (cami-1<k .AND. k<=camf)) THEN
                     ELSE
                        IF(i==(PsrHalf-nlin) .AND. j==2 .AND. (cami-1<k .AND. k<=camf)) THEN
                        ELSE
                           IF(k==cami-1 .AND.( i <= nlin .OR. i >PsrHalf-nlin ) .AND. j==1 .AND. MOD(i,2)==0) THEN
                           ELSE
                              WRITE(UNIT=11, FMT="(3(f9.5,';'), f9.5)")posicao, vetor(j,:)
                           END IF
                        END IF
                     END IF
                  END IF

               END DO
            END IF
         END IF
      END IF

   END DO

END DO

!Plotting the impurities positions
WRITE(nome, FMT="(a)") &
            "./dados/imp_positions.dat"
OPEN(UNIT=12, FILE=nome, STATUS="UNKNOWN", ACTION="WRITE")
DO i=1, N
   !For impurities on the left
   posicao =i*base1+i*base2-2.5d0*vetor(1,:) 
   DO j=1,impurezaux(i,1,1)
         posicao(2) =(impurezaux(i,1,j+1)+0.5d0)*SQRT(3.d0)
         WRITE(UNIT=12, FMT="(1(f9.5,';'), f9.5)")posicao
   END DO
   !For impurities on the right
   posicao=i*base1+i*base2-vetor(1,:)
   DO j=1,impurezaux(i,2,1)
         posicao(2)=impurezaux(i,2,j+1)*SQRT(3.d0)
         WRITE(UNIT=12, FMT="(1(f9.5,';'), f9.5)")posicao
   END DO
END DO
!END Plotting the impurities positions

CALL cpu_time (end_cpu_time)
timedif=end_cpu_time - beg_cpu_time
PRINT "('Total running time= ',i0,'h:',i0,'m:',i0,'s')", &
                                 INT(timedif/3600.0), & 
                                 INT(MOD(timedif,3600.0)/60.0),&
                                 INT(MOD(MOD(timedif,3600.0),60.0))  

END PROGRAM corrente