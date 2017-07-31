SUBROUTINE sortimpureza()

USE parametros
IMPLICIT NONE
INTEGER :: i,j,l,k,m,q,w,y,u, nimpur, impboraux, impboraux2, counter, left, right, &
           num_hex, aux
INTEGER, DIMENSION(:,:), ALLOCATABLE :: sitiosorteados2
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: sitiosorteados
! INTEGER, INTENT(OUT) :: impureza(namostras,N,2,P+1)
INTEGER :: impurezaux(N,2,P+1)
REAL(KIND=8) :: num_aleatorio
REAL(KIND=8) :: vetor(3,2), base1(2), base2(2), posicao(2),&
                pos_center(2), diff(2)
LOGICAL :: verificador

IF(assm.NE.0) THEN
   PRINT *, "'assm' should be 0. You may also check the 'impureza.f90' file. Stopping."
   STOP
END IF

impureza=0
IF(concentracao.LE.0.d0) RETURN

!Provides a ordered distribution of impurities
IF(ordered==1) THEN
   IF(namostras/=1) THEN
      PRINT *, "Verify namostras in the parametros.f90. It should be 1."
      STOP
   END IF
   PRINT *, "Ordered dist. with ctc of ~ 0.28!"
   !Forca arranjo periodico de impurezas
   DO i=1, N
      impureza(1,i,2,1)=(P+1)/2
      DO j=1, (P+1)/2
         impureza(1,i,2,1+j)=2*j-1
      END DO
   END DO
   impureza(:,1,1,:)=0
   RETURN
END IF

CALL inicializa_semente()
IF(concentracao>1.d0) THEN
   PRINT *, "Parameter 'concentracao' should be <=1! Stopping."
   STOP
END IF
!    nimpur=concentracao*(N*(2*P-1)-P+1)
num_hex=N*(2*P-1-assm)-P+1
nimpur=concentracao*num_hex
IF(ordered==2) THEN
   nimpur=N*(P+1)/2
   PRINT *, "Random dist. with ctc of ordered dist. ~ 0.28!"
END IF

ALLOCATE(sitiosorteados(namostras,nimpur,2))
DO m=1, namostras   
   DO l=1, nimpur
      DO
         !i will indicate one column, it's odd for left and even for right
         CALL RANDOM_NUMBER(num_aleatorio)
         i=num_aleatorio*2*N+1
         IF(MOD(i,2)==1 .AND. P==1) CYCLE !Avoid impurities at left with P=1

         CALL RANDOM_NUMBER(num_aleatorio)
!             j=num_aleatorio*(P-MOD(i,2))+1
         j=num_aleatorio*( P - MOD(i,2)*(1-assm) )+1
         
         verificador=.TRUE.
         IF(i.EQ.1) THEN
            verificador=.FALSE.
         ELSE
            DO k=1, l-1
               IF( (i.EQ.sitiosorteados(m,k,1)).AND.&
                   (j.EQ.sitiosorteados(m,k,2))      ) THEN
                  verificador=.FALSE.
                  EXIT
               END IF
            END DO
         END IF
         IF(verificador) EXIT
      END DO
      impureza(m,(i+1)/2,2-MOD(i,2),1)=impureza(m,(i+1)/2,2-MOD(i,2),1)+1
      impureza(m,(i+1)/2,2-MOD(i,2),1+impureza(m,(i+1)/2,2-MOD(i,2),1))=j
      sitiosorteados(m,l,1)=i
      sitiosorteados(m,l,2)=j
   END DO
END DO

!Cleans part of the central device to extend the contacts
IF(Nclean.NE.0) THEN
   impureza(:,1:Nclean,:,1)=0
   impureza(:,Nclean+1,1,1)=0
   impureza(:,N-Nclean+1:N,:,1)=0
END IF

!Cleans part of the central device to extend the contacts
IF(cleanHoriz.NE.0) THEN
   IF(assm.NE.0) THEN
      PRINT *, "impureza.f90: Clean some impurerities horizotally was not modified to work with assm=1. Stopping."
      STOP
   END IF
   
   y=1; u=N
   IF(nlin.NE.0) THEN
      y=cami; u=camf+1
   END IF

   DO m=1,namostras
      DO i=y,u !Loop over number of layers that composes the central device
         DO j=1,2  !Loop over left or right layer
            IF(i==camf+1 .AND. j==2) CYCLE
            q = 0
            DO k=1,impureza(m,i,j,1) !Loop over Num. of impurities in this layer
               q = q + 1 !we define q instead of using k cause when a impurity is removed we need to start checking the array of imp. at the same position but k moves on.
               IF( impureza(m,i,j,1+q) <= cleanHoriz .OR. impureza(m,i,j,1+q) > (P-MOD(j,2))-cleanHoriz  ) THEN
                  DO l=1,impureza(m,i,j,1)-q !This removes one imp. of the array by sliding the entry on top of it.
                     impureza(m,i,j,q+l)=impureza(m,i,j,1+q+l)
                  END DO
                  q = q - 1 !To continue checking the array of imp. at the same position where a bad imp. was found.
                  impureza(m,i,j,1) = impureza(m,i,j,1) - 1
               END IF
            END DO  !Loop over Num. of imp.
         END DO !Loop over left or right
      END DO !Loop over layers
   END DO !Loop over sample
END IF

! !Cleans horizotally impurities of the edges
! IF(cleanHoriz.NE.0) THEN
!    IF(ordered.NE.0) THEN
!       PRINT *, "'cleanHoriz' and 'ordered' shall not be .NE. 0 at the same time. &
!                &Check 'parametros.f90' . Stopping."
!       STOP
!    END IF
!    IF(assm.NE.0) THEN
!       PRINT *, "impureza.f90: Clean some impurerities horizotally was not modified to work with assm=1. Stopping."
!       STOP
!    END IF
!    impureza=0
!    impboraux2=0
!    impboraux=impborda
!    IF(impborda.EQ.0) impboraux=1
!    IF(impborda.NE.0) impboraux2=1

!    impureza(:,:,1,1)=P-1-2*(cleanHoriz+impboraux2) !esq
!    impureza(:,:,2,1)=P  -2*cleanHoriz !dir

!    DO i=1,impborda-1
!       impureza(:,:,1,1+i)=i
!    END DO
!    DO i=1,P-1-2*cleanHoriz-2*impborda
!       impureza(:,:,1,1+impboraux-1+i)=impborda+cleanHoriz+i
!    END DO
!    DO i=1,impborda-1
!       impureza(:,:,1,1+impboraux-1+P-1-2*cleanHoriz-2*impborda+i)=P-1-(impboraux-1)+i
!    END DO

!    DO i=1,impborda
!       impureza(:,:,2,1+i)=i
!    END DO
!    DO i=1,P-2*cleanHoriz-2*impborda
!       impureza(:,:,2,1+impborda+i)=impborda+cleanHoriz+i
!    END DO
!    DO i=1,impborda
!       impureza(:,:,2,1+impborda+P-2*cleanHoriz-2*impborda+i)=P-impborda+i
!    END DO
!    impureza(:,1,1,:)=0
! END IF


!Creates clusters of impurities
IF(radio .NE. 0.d0) THEN
   IF(assm.NE.0) THEN
      PRINT *, "impureza.f90: Clusters of impurities was not modified to work with assm=1. Stopping."
      STOP
   END IF
   ALLOCATE(sitiosorteados2(num_hex,3))

   vetor(1,:)=(/1.d0,0.d0/)
   vetor(2,:)=(/-0.5d0, SQRT(3.d0)/2.d0/)
   vetor(3,:)=(/-0.5d0,-SQRT(3.d0)/2.d0/)

   base1=vetor(1,:)-vetor(3,:)
   base2=vetor(1,:)-vetor(2,:)
   
   !Loops over each cluster
   DO m=1, namostras; 
      counter=0
      impurezaux(:,:,:)=0
Ot:   DO y=1,nimpur
         q=(sitiosorteados(m,y,1)+1)/2
         w=sitiosorteados(m,y,2)
         l=MOD(sitiosorteados(m,y,1)+1,2)+1
         IF(l==1) THEN
            IF(q==1) CYCLE !Don't let cluster center be on the left of first layer
            left=1; right=0
         ELSE
            left=0; right=1
         END IF
         pos_center=q*base1+q*base2-2.5d0*vetor(1,:)*left-vetor(1,:)*right
         pos_center(2)=(w+0.5d0*left)*SQRT(3.d0)
         !Looking for hexagon inside of the cluster
         DO i=1,N; DO k=1,2
            IF(k==1) THEN
               IF(i==1) CYCLE
               left=1; right=0; 
            ELSE
               left=0; right=1
            END IF
            posicao=i*base1+i*base2-2.5d0*vetor(1,:)*left-vetor(1,:)*right
            DO j=1,P-1*left
               posicao(2) =(j+0.5d0*left)*SQRT(3.d0)
               diff=posicao-pos_center
               IF( SQRT(diff(1)**2+diff(2)**2)<=radio ) THEN
                  verificador=.FALSE.
                  DO u=1, counter
                     IF( (i.EQ.sitiosorteados2(u,1)) .AND.&
                         (j.EQ.sitiosorteados2(u,2)) .AND.&
                         (k.EQ.sitiosorteados2(u,3)) ) THEN
                        verificador=.TRUE.
                        EXIT
                     END IF
                  END DO
                  IF(verificador) CYCLE
                  impurezaux(i,k,1)=impurezaux(i,k,1)+1
                  aux              =impurezaux(i,k,1)
                  impurezaux(i,k,aux+1)=j

                  counter=counter+1
                  sitiosorteados2(counter,:)=[i,j,k]
                  IF(counter==nimpur) EXIT Ot
               END IF
            END DO !DO j=1,P-1*left
         END DO; END DO !DO i=1,N; DO k=1,2

         impureza(m,:,:,:)=impurezaux
      END DO Ot ! Ot: DO y=1,nimpur 
   END DO ! DO m=1, namostras 

!    PRINT *, "impureza.f90: Impurities on clusters mode on."
!    PRINT *, "Concentration= ", REAL(counter)/REAL(num_hex)
END IF

impureza(:,1,1,:)=0
RETURN
END SUBROUTINE sortimpureza

!========================================================================

SUBROUTINE inicializa_semente()
   IMPLICIT NONE
   INTEGER :: i, n, relogio
   INTEGER, DIMENSION(:), ALLOCATABLE :: semente

   CALL RANDOM_SEED(size = n)
   ALLOCATE( semente(n) )

   !CALL SYSTEM_CLOCK(COUNT=relogio)
   !semente = relogio+37*(/ (i-1, i = 1, n) /)
   semente = (/ (i-1, i = 1, n) /)

   CALL RANDOM_SEED(PUT = semente)
   DEALLOCATE(semente)

END SUBROUTINE inicializa_semente