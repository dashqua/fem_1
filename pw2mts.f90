program pw2mts

!cccccccccccccccccccccccccccc
! DECLARATION OF USED MODULUS
!cccccccccccccccccccccccccccc

  USE mod_lec_fic
  USE SysLinCSR
  USE Precond
  USE BLASS
  USE EFLib

!ccccccccccccccccccccccccc
! DECLARATION OF VARIABLES
!ccccccccccccccccccccccccc

  implicit none
  integer :: NN,NT,NNE
  integer, dimension(:,:), allocatable :: IND
  real (kind=8), dimension(:,:), allocatable :: VAL
  real (kind=8), dimension(:), allocatable :: A,RHS,SOL
  integer, dimension(:), allocatable :: NBE,IA,JA
  integer :: i,j
  type(msh) :: mesh
  TYPE(MatriceCSR) :: Mat2D
  real(kind=8) :: error

   REAL (kind=8), dimension(3,3) :: AL
   REAL (kind=8), dimension(3)   :: RHSL

!ccccccccccccccccccccccccccc
! LOAD OF THE MESH FROM GMSH
!ccccccccccccccccccccccccccc

  call load_gmsh('wing_naca.msh',mesh)

!cccccccccccccccccccccccccccccccccccccccccccccccccc
! ASSEMBLING OF THE STIFFNESS MATRIX AND OF THE RHS
!cccccccccccccccccccccccccccccccccccccccccccccccccc
  NT=mesh%nbTriangles
!!$  print*, 'Nb of Triangles= ', NT ; print*
  NN=mesh%nbNod
!!$  print*, 'Nb of Nodes= ', NN ; print*

   !CALL Compute_LocalStiffness_Matrix_And_LocalRHS(mesh,35082,AL,RHSL)
   !print*, 'AL= ', AL
   !print*, 'RHSL= ', RHSL

  CALL  Assembling_Global_Stiffness_Matrix_And_RHS(mesh,IND,VAL,NBE,RHS)

!cccccccccccccccccccccccccccccccccccccccccccccccccc
! Put the matrix of the linear system in CSR format
!cccccccccccccccccccccccccccccccccccccccccccccccccc

  CALL Put_In_CSR_Format(IND,VAL,NBE,NNE,A,JA,IA)
  deallocate(IND,NBE,VAL)

!ccccccccccccccccccccccccccccccc
! Take Dirichlet BC into account
!ccccccccccccccccccccccccccccccc

  CALL Boundary_Conditions_Wing(mesh,IA,JA,A,RHS)
 
!cccccccccccccccccccccccccccccccccccccccccccccc
! CREATION AND AFFECTATION OF THE DATA IN Mat2D 
!cccccccccccccccccccccccccccccccccccccccccccccc

  CALL CREATE(Mat2D,NN,NNE)
  ALLOCATE (Mat2D%Guess(Mat2d%NbLine)) ! Needed for iterative method
  ! Copy data in the matrix structure Mat2D
  Mat2D%Resolution = 1
  Mat2D%Coeff(:) = a(:)
  Mat2D%IA(:) = ia(:)
  Mat2D%JA(:) = ja(:)
  Mat2D%Guess(:) = 0.

  !     set-up the preconditioner ILUT(15, 1E-4)  ! new definition of lfil
  !
  Mat2D%Precond=3        ! Precodintioner type: ILUT
  Mat2D%LFill=50          ! Filling level
  Mat2D%DropTol=1.0D-6   ! Tolerance
  ! Computation of the preconditioner
  CALL CalPrecond(Mat2D)

!cccccccccccccccccccccccccccc
! COMPUTATION OF THE SOLUTION 
!cccccccccccccccccccccccccccc
  allocate(sol(NN))
  ! Computation of the solution to the linear
  CALL RESOL(Mat2D,RHS(1:NN),SOL(1:NN))

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! WRITE IN OUTPUT FILES FOR MATLAB POST-PROCESSING
!ccccccccccccccccccccccccccccccccccccccccccccccccc

  open(unit=100,file='mat_fish.dat',status='unknown')
  do i=1,nn
     do j=ia(i),ia(i+1)-1
        write(100,*) i,ja(j),a(j)
     end do
  end do
  close(100)

  open(100,file='rhs_fish.dat',status='unknown')
  do i=1,nn
     write(100,*) rhs(i)
  end do
  close(100)

  open(100,file='nodes_fish.dat',status='unknown')
  do j=1,nn
     write(100,*) mesh%pos(j,:)
  end do
  close(100)
 
  open(unit=100,file='solutionnaca.dat',status='unknown')
  do i=1,nn
        write(100,*) sol(i)
  end do
  close(100)

!ccccccccccccccccccccccccc
! COMPUTATION OF THE ERROR
!ccccccccccccccccccccccccc
  !error=0.0
  !CALL Compute_Error(mesh,sol,error)
  !print*, 'ERROR= ', error  

!ccccccccccccc
! DEALLOCATION
!ccccccccccccc

  CALL DESTROY(Mat2D)
  deallocate(A,JA,IA,RHS,SOL)

end program pw2mts
