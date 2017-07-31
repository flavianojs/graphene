subroutine invers(matriz,NN)
 
 use f90_kind

 integer :: NN,INFO
 integer :: LWORK 
 integer, DIMENSION(NN) :: IPIV

 complex(double), DIMENSION(NN,NN) :: matriz
 complex(double), DIMENSION(NN*4) :: WORK  


  LWORK = 4*NN
  INFO = 0
  !CALL F07ARF(NN,NN,matriz,NN,IPIV,INFO)
  !if (INFO/=0) then
  ! write(*,*)'IFAIL=',INFO
  !end if
  !CALL F07AWF(NN,matriz,NN,IPIV,WORK,LWORK,INFO)
  CALL zgetrf(NN,NN,matriz,NN,IPIV,INFO)
  if (INFO/=0) then
   write(*,*)'IFAIL=',INFO
  end if
  CALL zgetri(NN,matriz,NN,IPIV,WORK,LWORK,INFO)

  if (INFO/=0) then
    write(*,*)'IFAIL=',INFO
    stop 'FUDEU!!! INVERSAO DA NAG FOI PRO CARALHO!'
  end if

  return

end subroutine invers
