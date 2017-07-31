! -------     ______           ______           ______      ---------
! \     |    /      \         /      \         /      \     |    /
!  \    |   /        \       /        \       /        \    |   /
!   \___|__/          \_____/          \_____/          \___|__/       2
!   /   |  \          /     \          /     \          /   |  \
!  /    |   \        /       \        /       \        /    |   \
! /     |    \______/         \______/         \______/     |    \___
! \     |    /      \         /      \         /      \     |    /
!  \    |   /        \       /        \       /        \    |   /
!   \___|__/          \_____/          \_____/          \___|__/
!   /   |  \          /     \          /     \          /   |  \     P=1
!  /    |   \        /       \        /       \        /    |   \
! /     |    \______/         \______/         \______/     |    \___
! -------      N=1               2                 3        ----------
MODULE parametros
IMPLICIT NONE
CHARACTER(LEN=80), PARAMETER :: &
                     tiposimp(4)=(/"Atop    ","ImpSubst","ImpHex  ","ImpHexSO"/)
! Que tipo de calculo serah feito?
CHARACTER(LEN=80) :: tipoimp=tiposimp(4)
! Calculos de condut. com vacancia ou constricoes devem usar flag=2
! Metodo fg p/ condutancia: 1-add layer; 2-invers. compl.; 3-add layer extremo
INTEGER :: flag=2, uniquename=0
LOGICAL :: squarelead=.FALSE., verbose=.TRUE.

!N=num de sitios horz., P=n. sitios vert. como mostra acima
INTEGER, PARAMETER :: P=41, N=20, assm=0 !assm should be 0 or 1
INTEGER, PARAMETER :: npt=500, namostras=1
REAL(KIND=8), PARAMETER :: t=-2.7d0,t2=-0.d0, eta=1.d-7, Ef=0.0d0, difpot=0.d0,&
                           eps0=0.0d0, interv=5.d0, fracaot=0.0d0, delta=-0.0d0
REAL(KIND=8) :: potesq=eps0, potdir=eps0, ee=0.0d0,&!*ABS(t), &!0.105d0 0.1261
eini=-2.d0 , efin=2.d0, tb=t*(1.d0+delta)

!++++++++++++++++ SPIN-ORBITA ++++++++++++++++++++++++++++++++++++++++++++++++++
! Qual a distribuicao? 1 p/ ordenada; 2 random mesm ctc do ord. e 0 para
! random com ctc=concentracao
! 'cleanHoriz' cleans impurities horizontally near the edges. It imposes ctc=1.0
! in the central part. Should be kept 0 for no effect.
! 'impborda' gives a strip of impurities with SO on the edges
! If 'radio' is different from zero, will be created clusters of impurity with
! the given radio.
! 'Nclear' should be zero for non-extended leeds
! 'radio' must be zero for no cluster
! 'deltau': uncorrelated random onsite disorder in the range [-deltau,deltau]
! in all the central region. It must be kept 0.d0 to turn it off.
INTEGER :: ordered=0, cleanHoriz=14, impborda=0, Nclean=4
REAL(KIND=8), PARAMETER :: concentracao=1.d0, radio=0.0d0, deltau=0.0d0*t
COMPLEX, PARAMETER ::         lambda=(-0.00d0,-0.270d0) *1.0d0
COMPLEX(KIND=8), PARAMETER :: aatimp=(-0.24d0,-0.000d0) *0.0d1
REAL(KIND=8), PARAMETER ::       rho=-0.10d0            *0.0d0
REAL(KIND=8), PARAMETER ::        VV=-0.83d0            *0.0d0
CHARACTER(LEN=8), PARAMETER :: term="constric"
! CHARACTER(LEN=22), PARAMETER :: term="laONtimpOFee0.01ctc0.1"
!++++++++++++++++ IMPUREZAS ++++++++++++++++++++++++++++++++++++++++++++++++++++
REAL(KIND=8), PARAMETER :: fracttimp=0.d0
   !++++++++++++++++ IMPUREZAS REAIS NO CENTRO DO HEXAGONO +++++++++++++++++++++
   REAL(KIND=8), PARAMETER :: fracteps1=0.0d0
   REAL(KIND=8), PARAMETER :: eps1=fracteps1*ABS(t), timp=fracttimp*t
   !++++++++++++++++ IMPUREZA ATOP & VACANCIA ++++++++++++++++++++++++++++++++++
   !Em qual sitio deve estah a impureza atop?
   !Deve-se manter sitio=0 para nao colocar a impureza atop.
   !ATENCAO: fracttimp=0 para realizar uma vacancia propriamente.
   INTEGER, PARAMETER :: sitio=0!(4*P+2)/2
   !Em que sitio deve ser calculado ldos em funcao da energia?
   INTEGER :: sitioviz=1

!++++++++++++++++ CONSTRICOES ++++++++++++++++++++++++++++++++++++++++++++++++++
! INTEGER :: cami=N/2, camf=N/2+9, nlin=0!((P-5)*2+1+5)/2
INTEGER :: cami=(N-12)/2 + 1, camf=N-(N-12)/2, nlin=0!((P-5)*2+1+5)/2
                                   !Deve permanecer 0 se nao quiser constricoes
                                   !Only works for flag=2

!+++++++++++++++++++ CURRENTS ++++++++++++++++++++++++++++++++++++++++++++++++++
!Se 1, plota apenas setas horiz. Se 3, plota todas
INTEGER :: flaglc=3
!To remove all impureties from left or right layers for a precise total cross-section current computation when SOC is on
LOGICAL ::  remove_imp_ONleftORright = .FALSE. 
  INTEGER :: lado=1!lado=1->esq, lado=0->dir      
REAL(KIND=8) :: fator=1.d2, fator2=0.d0, rescale_maior=1.d0

!++++++++++++++++ STM DUAL PROBE +++++++++++++++++++++++++++++++++++++++++++++++
INTEGER, PARAMETER :: lpB=(P/2+1+5)+(N/2+1), jpB=-(P/2+1+5)+(N/2+1), &
                      npmesh=16000/(P*N)!36!will not be exact this number
REAL(KIND=8), PARAMETER :: aa=0.6d0, llamb=0.85d0, dCC=1.42d0, fatexp=0.0d0, &
                           centroHex=0.d0 !1 para centro, 0 para sobre atomo
REAL(KIND=8), PARAMETER :: raio=2.0d0*dCC+1.d-10

!CONSTANTES UNIVERSAIS and NON AJUSTABLES
REAL(KIND=8), PARAMETER :: pi=2.d0*ACOS(0.d0)
INTEGER, PARAMETER :: Psring2=4*P+2
INTEGER, PARAMETER :: PsrHalf2=Psring2/2
INTEGER :: impureza(namostras,N,2,P+1)!, impurezaux(N,2,P+1)
REAL(KIND=8), ALLOCATABLE :: anderdisor(:,:,:,:)
END MODULE parametros

!DICTIONARY
!*delta: renomalization factor for the nn hopping on the edges of the AGNR.
!        tb=t*(1+delta)
!*Nclean: Number of layers to be clear at the left and right sides of the
! central region to extend the leed. Should be keeped zero for non-extended
! leed.
!*uniquename: Gives the output data in a file with a unified name independent
! of parameter
