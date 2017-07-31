INCLUDE "global.f90"
INCLUDE "./condutanciaImpHexLB.f90"
INCLUDE "impureza.f90"
INCLUDE "andersondisorder.f90"
INCLUDE "printtime.f90"
PROGRAM condutvsenergia

USE MPI   
USE parametros
USE global

IMPLICIT NONE
INTEGER, PARAMETER :: Psring=4*P+2
INTEGER, PARAMETER :: PsrN=Psring*N, PsrHalf=Psring/2
INTEGER :: i,j,l,k,m,sample, nimpur,ss,mm,hh
INTEGER, DIMENSION(:,:), ALLOCATABLE :: sitiosorteados
! INTEGER :: impureza(namostras,N,2,P+1)!, impurezaux(N,2,P+1)
! REAL(KIND=8) :: anderdisor(namostras,N,2,PsrHalf)
REAL(KIND=8) :: medcondt(3),desvmed(3),vetcondt,condt1,condt2
REAL(KIND=8) :: incremento, condt, e, num_aleatorio, dados(6)
REAL(KIND=8) :: grafico(npt,namostras,4)
CHARACTER(LEN=300) :: nome, form3
CHARACTER(LEN=37) :: prefix(3)
CHARACTER(LEN=3) :: percentage
LOGICAL :: verificador
!MPI parameters
INTEGER :: numprocs,errorcode,ierr,itask,masterrank=0
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
CALL MPI_Init(ierr)
CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
CALL MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)

CALL cpu_time (beg_cpu_time)
IF(myrank==masterrank) THEN
   prefix(1)="./dados/CondutImpHexSOLB_SimAddAcRB_P"
   prefix(2)="./dados/CondutImpHexSOLB_CompIvAcRB_P"
   prefix(3)="./dados/CondutImpHexSOLB_ExtAddAcRB_P"

   !Error verifications
   IF(numprocs>npt*namostras) THEN
      PRINT "(a,i0,a,i0,a)", &
      "condutvsenergiaImpHex.f90: (Number of processors) greater than&
      & (npt x namostras): ",numprocs," > ",npt*namostras, ". Stopping program!"
      CALL MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
   END IF
   IF(flag<1 .OR. flag>4) THEN
      PRINT *, "condutvsenergiaImpHex.f90: Flag incorrect!" 
      PRINT *, "Change it in the inputcard for a number between 1 and 3."
      PRINT *, "Stopping!"
      STOP
   END IF

   WRITE(nome, FMT="(2(a,i0),a,f4.1,2(a,f6.3,a,f5.2),a,f6.3,&
                              &a,f5.2,2(a,f5.2),a,f4.1,3(a,i0),2(a,f6.3),a)") &
      prefix(flag),P,"N",N,"t",t,"ctc",concentracao,&
      " nu(",REAL(lambda),",", IMAG(lambda),") tip(",REAL(aatimp),&
      ",", IMAG(aatimp),") rho",rho," VV",VV," dtu",deltau," rad",radio," ams",&
      namostras," npt",npt," sit",sitio," ene[",eini,",",efin,"].dat"

   IF(uniquename==1) WRITE(nome, FMT="(a,i0,a,i0,a,i0,a)")&
                 "./dados/CondutImpHexSOLB_",P,"x",N,"lamtestflag",flag,".dat"
   OPEN(UNIT=10, FILE=nome, STATUS="UNKNOWN", ACTION="WRITE")
   WRITE(UNIT=10, FMT="(a)") "# energy, meanconducttotal, meandeviat, &
              &meanconductSpinUP, meandeviat, meanconductSpinDOWN, meandeviat"
   nome(10:10)="A"
   OPEN(UNIT=11, FILE=nome, STATUS="UNKNOWN", ACTION="WRITE")

   CALL cpu_time (end_cpu_time)
   timedif=end_cpu_time - beg_cpu_time
   hh=INT(timedif/3600.0)
   mm=INT(MOD(timedif,3600.0)/60.0)
   ss=INT(MOD(MOD(timedif,3600.0),60.0))

END IF !If myrank==masterrank

CALL sortimpureza() 
CALL sortanderdisor() 

IF(npt.EQ.1) THEN
incremento=(efin-eini)
ELSE
incremento=(efin-eini)/(npt-1)
END IF
k=myrank
itask = numprocs-1

!Loop through the energy points
DO
   sample=MOD(k,namostras)+1
   k=k/namostras
   e=eini+incremento*REAL(k)
   vetcondt=0.d0
   CALL printtime("Spin UP")
   CALL condutancia(e, condt, sample,  1.d0 )
   condt1=condt
   CALL printtime("Spin DOWN")
!    CALL condutancia(e, condt, sample, -1.d0 )
   condt2=condt
   dados(:)=(/DBLE(k),DBLE(sample),e,condt1+condt2,condt1,condt2/)

   !MPI Receiving and Sending data
   IF (myrank==masterrank) then
      grafico(1,1,:)=dados(3:6)
      DO i=2, npt*namostras
         CALL MPI_Recv(dados,6,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                                       43,MPI_COMM_WORLD,stat,ierr)
         grafico( INT(dados(1))+1, INT(dados(2)), : )=dados(3:6)
         WRITE(UNIT=666, FMT="(4(ES24.15E2))") dados(3:6)
         
         ! If the number of processors is smaller than the total number of points, 
         ! sends the rest of the points to the one which just finished
         IF (itask<npt*namostras-1) THEN
            itask = itask + 1

            WRITE(percentage,FMT="(i0)") 100*(itask+1)/(npt*namostras)
            CALL printtime("Energy/sample Sent: "//percentage//"%   ")  
            
            CALL MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),66,MPI_COMM_WORLD,ierr)
         ELSE
            CALL MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),66,MPI_COMM_WORLD,ierr)
         END IF            
      END DO
      EXIT
   ELSE
      CALL MPI_Send(dados,6,MPI_DOUBLE_PRECISION,masterrank,43,MPI_COMM_WORLD,ierr)
      CALL MPI_Recv(k,1,MPI_INTEGER,masterrank,66,MPI_COMM_WORLD,stat,ierr)
      IF(k==0) EXIT
   END IF
   !END MPI  Receiving and Sending data

END DO! Infinite loop until the energy points end

IF (myrank==masterrank) THEN
   DO i=1,npt
      DO j=2,4
         l=j-1
         medcondt(l)=0.d0
         DO m=1, namostras
            medcondt(l)=medcondt(l)+grafico(i,m,j)
         END DO
         medcondt(l)=medcondt(l)/DBLE(namostras)
         desvmed(l)=0.d0
         DO m=1, namostras
            desvmed(l)=desvmed(l)+ABS(grafico(i,m,j)-medcondt(l))
         END DO
         desvmed(l)=desvmed(l)/DBLE(namostras)
      END DO
      WRITE(UNIT=10, FMT="(7(ES24.15E2))") grafico(i,1,1), &
      medcondt(1), desvmed(1), medcondt(2), desvmed(2), medcondt(3), desvmed(3)
   END DO

   DO m=1, namostras
      DO i=1, npt
         WRITE(UNIT=11, FMT="(4(ES24.15E2))") grafico(i,1,1), grafico(i,m,2:4)
      END DO
         WRITE(UNIT=11, FMT=*)
   END DO

   CALL printtime("PROGRAM END")
END IF

!      WRITE(UNIT=10, FMT="(ES24.15E2, ES24.15E2, ES24.15E2)") e, medcondt, desvmed
!    END DO

CALL MPI_Finalize(ierr)
IF (ierr/=0) THEN
   WRITE(*,"('ierr = ',i0,'. Something went wrong in the parallelization.')") ierr
END IF

END PROGRAM condutvsenergia
