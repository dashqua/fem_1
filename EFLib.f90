module EFLib

  USE mod_lec_fic
  USE SysLinCSR
  USE Precond
  USE BLASS

  IMPLICIT NONE

  CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Compute_LocalStiffness_Matrix_And_LocalRHS(mesh,K,ALOC,RHSLOC)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  IMPLICIT NONE

  type(msh),intent(in) :: mesh
  INTEGER, intent(in)  :: K

  REAL (kind=8), dimension(3,3), intent(out) :: ALOC
  REAL (kind=8), dimension(3), intent(out) :: RHSLOC

  !integer :: nbtriangles, i, j, nodenumber, nodeX, nodeY
  real (kind=8), dimension(2,3) :: GRADPHI
  real (kind=8), dimension(2,2) :: Mkt
  real :: sizem

  GRADPHI(1,1) = -1
  GRADPHI(1,2) = 1
  GRADPHI(1,3) = 0
  GRADPHI(2,1) = -1
  GRADPHI(2,2) = 0
  GRADPHI(2,3) = 1
  
  Mkt(1,1) = mesh%pos(mesh%triangles(K,3),2) - mesh%pos(mesh%triangles(K,1),2)
  Mkt(2,1) = mesh%pos(mesh%triangles(K,1),1) - mesh%pos(mesh%triangles(K,3),1)
  Mkt(1,2) = mesh%pos(mesh%triangles(K,1),2) - mesh%pos(mesh%triangles(K,2),2)
  Mkt(2,2) = mesh%pos(mesh%triangles(K,2),1) - mesh%pos(mesh%triangles(K,1),1)

  ALOC = matmul( transpose(matmul(Mkt,GRADPHI)) , matmul(Mkt,GRADPHI) )
  !ALOC = dot_product( matmul(Mkt,GRADPHI)) , matmul(Mkt,GRADPHI)) )
  sizem = (  Mkt(1,1)*Mkt(2,2) - Mkt(2,1)*Mkt(1,2)   )
  ALOC = ALOC/(2*sizem)
  
  RHSLOC(1) = source( mesh%pos(mesh%triangles(K,1),1) , mesh%pos(mesh%triangles(K,1),2) )
  RHSLOC(2) = source( mesh%pos(mesh%triangles(K,2),1) , mesh%pos(mesh%triangles(K,2),2) )
  RHSLOC(3) = source( mesh%pos(mesh%triangles(K,3),1) , mesh%pos(mesh%triangles(K,3),2) )

  RHSLOC = (sizem/6.) * RHSLOC
  


  














  




  


  
  end subroutine Compute_LocalStiffness_Matrix_And_LocalRHS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Assembling_Global_Stiffness_Matrix_And_RHS(mesh,IND,VAL,NBE,RHS)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  IMPLICIT NONE

  type(msh),intent(in) :: mesh

  INTEGER, dimension(:,:), allocatable, intent(out) :: IND
  REAL (kind=8), dimension(:,:), allocatable, intent(out) :: VAL
  INTEGER, dimension(:), allocatable, intent(out) :: NBE
  REAL (kind=8), dimension(:), allocatable, intent(out) :: RHS

  INTEGER :: NT,NN,K,L,LIG,J,I
  INTEGER, dimension(3) :: TRI
  REAL (kind=8), dimension(3,3) :: ALOC
  REAL (kind=8), dimension(3)   :: RHSLOC
  
  NT=mesh%nbTriangles
  NN=mesh%nbNod
! The value "10" is the maximal number of neighbours that a node can have.
  allocate(IND(NN,10),VAL(NN,10),NBE(NN))
  allocate(RHS(NN))
  IND=0; VAL=0; NBE=0; RHS=0

  DO K=1,NT
    CALL Compute_LocalStiffness_Matrix_And_LocalRHS(mesh,K,ALOC,RHSLOC)
    TRI=mesh%triangles(K,1:3)
    DO L=1,3
     LIG=TRI(L) ! number of the current node
     RHS(LIG)=RHS(LIG)+RHSLOC(L)
     IF (NBE(LIG).eq.0) then
! If the node was never been encoutered
! One put every nodes of the triangle in the connectivity table
       IND(LIG,1:3)=TRI 
! The number of elements in the row corresponding to the current node
! becomes 3
       NBE(LIG)=3
! Assembling
       VAL(LIG,1:3)= ALOC(L,1:3)
     ELSE
! The node was already been visited in a previous triangle
! It is therefore needed to examine every neighbouring nodes 
! (including itself)
       OUT: DO I=1,3
        ! This loop allows to avoid to add already existing nodes
        ! in the connectivity table. If a connectivity already
        ! exists, one adds to the matrix component the local matrix.
         DO J=1,NBE(LIG)
          IF (IND(LIG,J).EQ.TRI(I)) THEN
            VAL(LIG,J)=VAL(LIG,J)+ALOC(L,I) 
            CYCLE OUT
          END IF
         END DO
        ! If no connection was present, we create it.
        ! One assembles the element thanks to the local matrix.
              NBE(LIG)=NBE(LIG)+1
              IND(LIG,NBE(LIG))=TRI(I)
              VAL(LIG,NBE(LIG))=ALOC(L,I)
       END DO OUT
    END IF
   END DO
  END DO

end subroutine Assembling_Global_Stiffness_Matrix_And_RHS

 
!ccccccccccccccccccccccccccccccccccccccccccccxxxccccccc
  subroutine Put_In_CSR_Format(IND,VAL,NBE,NNE,A,JA,IA)
!cccccccccccccccccccccccccccccccccccccccccxxxcccccccccc

  IMPLICIT NONE  

  INTEGER, dimension(:,:), intent(in), allocatable :: IND
  REAL (kind=8), dimension(:,:),intent(in), allocatable :: VAL
  INTEGER, dimension(:),intent(in), allocatable :: NBE

  REAL (kind=8), dimension(:), allocatable, intent(out) :: A
  INTEGER, dimension(:), allocatable, intent(out) :: IA, JA
  INTEGER, intent(out) :: NNE

  INTEGER :: J,LIG
  INTEGER, dimension(:), allocatable :: iwork
  INTEGER :: NN  

  NN=size(IND,1)
  NNE=SUM(NBE)
!  write(6,*) 'Number of non nul elements in the matrix: ',NNE
  ALLOCATE(A(NNE),JA(NNE),IA(NN+1))
  A=0.; JA=0; IA=0
  J=1
  DO LIG=1,NN
     A(J:J+NBE(LIG)-1)=VAL(LIG,1:NBE(LIG))
     JA(J:J+NBE(LIG)-1)=IND(LIG,1:NBE(LIG))
     IF (LIG==1) THEN
        IA(LIG)=1
     ELSE
        IA(LIG)=IA(LIG-1)+NBE(LIG-1)
     END IF
     J=J+NBE(LIG)
  END DO
  IA(NN+1)=IA(NN)+NBE(NN)

  ALLOCATE(iwork(max(nn+1,2*nne)))
  CALL CSORT(nn,a,ja,ia,iwork,.true.) 
  DEALLOCATE(iwork)

end subroutine Put_In_CSR_Format

!ccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Boundary_Conditions_Wing(mesh,IA,JA,A,RHS)
!cccccccccccccccccccccccccccccccccccccccccccccccccc
!! Take Dirichlet BC into account
!! Identification of nodes
!! For a node I on the boundary
!! It is needed to:
!!      - put 0 on every component on the row I excepted for column I
!!         which takes value 1
!!      - put 0 on every component on the column I excepted for row I
!!         which takes value 1
!! WHERE(JA.EQ.I) A=0.0
!! A(IA(I):IA(I+1)-1)=0.0
!! GETELM (I,I,a,ja,ia,iadd,.TRUE.)
!! A(IADD)=1.0
!! - One compute the complete stiffness matrix and the nodal forces Bi
!! - One multiplies the k-th column of the stiffness matrix by the
!!  value Tk, and one subtracts it to the nodal forces vector.
!!   B = B - A*(0,...,0,Tk,0,...0)
!!   call amux(n,x,y,a,ja,ia)
!! - The k-th row and k-th column are replaced by 0
!!   WHERE(JA.EQ.I) A=0.0
!!   A(IA(I):IA(I+1)-1)=0.0
!! - The component Akk is replaced by 1.
!!   GETELM (I,I,a,ja,ia,iadd,.TRUE.)
!!   A(IADD)=1.0
!! - The component Bk is replaced by Tk.
!!   B(k)=Tk

     IMPLICIT NONE

     type(msh),intent(in) :: mesh

     INTEGER, dimension(:), allocatable, intent(inout) :: IA,JA
     REAL (kind=8), dimension(:), allocatable, intent(inout) :: A
     REAL (kind=8), dimension(:), allocatable, intent(inout) :: RHS

     INTEGER, dimension(:,:), allocatable :: boundary_nodes
     INTEGER, dimension(:), allocatable :: vu
     REAL (kind=8), dimension(:), allocatable ::V1,V2
     REAL(kind=8) :: cl,val_cl,g1,uinf
     INTEGER :: k,j,pp1,pp2,type,Zone_Phys_Bord2,Zone_Phys_Bord1,lig,NN,iadd

   NN=mesh%nbNod
   allocate(boundary_nodes(mesh%nbLines,2))
   allocate(vu(NN))
   boundary_nodes=0
   vu=0
   k=0
   do j=1,mesh%nbLines
    pp1=mesh%LINES(j,1)
    pp2=mesh%LINES(j,2)
    type=mesh%LINES(j,3)
    if (vu(pp1).eq.0) then
        k=k+1
        boundary_nodes(k,1)=pp1
        boundary_nodes(k,2)=type
        vu(pp1)=k
     end if
     if (vu(pp2).eq.0) then
        k=k+1
        boundary_nodes(k,1)=pp2
        boundary_nodes(k,2)=type
        vu(pp2)=k
     end if
  end do

 allocate(V1(NN),V2(NN))
  V1=0.0 ; V2=0.0
  cl=0.
  uinf = 100
  Zone_Phys_Bord1=1001
  Zone_Phys_Bord2=1002
  do j=1,mesh%nbLines
     !! Identification of the boundary node to treat
     lig=boundary_nodes(j,1)
     type=boundary_nodes(j,2)
     if (type.eq.Zone_Phys_Bord2) then
        val_CL=0.
     else if (type.eq.Zone_Phys_Bord1) then
        val_CL=-uinf*mesh%pos(lig,2)
     end if
     !! - One multiplies the k-th column of the stiffness matrix by the
     !!  value Tk, and one subtracts it to the nodal forces vector.
     !!   B = B - A*(0,...,0,Tk,0,...0)
       V1(lig)=Val_CL
       call amux(nn,v1,v2,a,ja,ia)
       !print *,maxval(v2)
        rhs=rhs-v2
        v1(lig)=0.0
     !! - The k-th row and k-th column are replaced by 0
     !!   WHERE(JA.EQ.I) A=0.0
        WHERE(JA.EQ.LIG) A=0.0
        A(IA(LIG):IA(LIG+1)-1)=0.0
     !! - The component Akk is replaced by 1.
       g1=GETELM (lig,lig,a,ja,ia,iadd,.true.)
       A(IADD)=1.0
     !! - The component Bk is replaced by Tk.
       RHS(LIG)=Val_CL
     end do

  deallocate(boundary_nodes,vu,v1,v2)

end subroutine Boundary_Conditions_Wing

















!ccccccccccccccccccccccccccccc
function source(x,y) result(z)
!ccccccccccccccccccccccccccccc

      IMPLICIT NONE

      REAL(kind=8), intent(in) :: x,y
      REAL(kind=8) :: z
      z = 0.
      return
      
end function source

function UEX(x,y) result(z)
  implicit none
  REAL(kind=8), intent(in) :: x,y
  REAL(kind=8) :: r,z
  z = (x-3)**2 * (x+3)**2 * (y-3)**2 * (y+3)**2
end function UEX

function PHI1(x,y) result(z)
  implicit none
  REAL(kind=8), intent(in) :: x,y
  REAL(kind=8) :: r,z
  z = 1-x-y
  return
end function PHI1

function PHI2(x,y) result(z)
  implicit none
  REAL(kind=8), intent(in) :: x,y
  REAL(kind=8) :: r,z
  z = x
  return
end function PHI2

function PHI3(x,y) result(z)
  implicit none
  REAL(kind=8), intent(in) :: x,y
  REAL(kind=8) :: r,z
  z = y
  return
end function PHI3

!ccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Compute_Error(mesh,U,err)
!cccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      type(msh),intent(in) :: mesh
      real (kind=8), dimension(:), allocatable, intent(in) :: U
      real (kind=8), intent(out) :: err     

      real (kind=8), dimension(3) :: X, W
      real (kind=8), dimension(2,2) :: MK
      real (kind=8), dimension(2) :: XZ, FK
      real (kind=8) :: res, error, det, U1, U2, U3, temp1, temp2, temp3
      integer :: i,j,K,n,nbtriangles
      
      X(1) = 1./2 - (1./2)*sqrt(3./5)
      X(2) = 1./2
      X(3) = 1./2 + (1./2)*sqrt(3./5)
      W(1) = 5./18
      W(2) = 8./18
      W(3) = 5./18
      res = 0
      error = 0
      n = 3
      nbtriangles = mesh%nbtriangles
      
      do K=1,nbtriangles
         res = 0
         MK(1,1) = mesh%pos(mesh%triangles(K,2),1) - mesh%pos(mesh%triangles(K,1),1)
         MK(1,2) = mesh%pos(mesh%triangles(K,3),1) - mesh%pos(mesh%triangles(K,1),1)
         MK(2,1) = mesh%pos(mesh%triangles(K,2),2) - mesh%pos(mesh%triangles(K,1),2)
         MK(2,2) = mesh%pos(mesh%triangles(K,3),2) - mesh%pos(mesh%triangles(K,1),2)
         det = MK(1,1)*MK(2,2) - MK(1,2)*MK(2,1)
         U1 = U( mesh%triangles(K,1) )
         U2 = U( mesh%triangles(K,2) )
         U3 = U( mesh%triangles(K,3) )
         do i=1,n
            do j=1,n
               XZ(1) = X(j)
               XZ(2) = ( 1-X(j) ) * X(i)
               FK = matmul( MK , XZ )
               FK(1) = FK(1) + mesh%pos(mesh%triangles(K,1),1)
               FK(2) = FK(2) + mesh%pos(mesh%triangles(K,1),2)
               temp1 = U1*PHI1( X(j),(1-X(j))*X(i) )
               temp2 = U2*PHI2( X(j),(1-X(j))*X(i) )
               temp3 = U3*PHI3( X(j),(1-X(j))*X(i) )
               res = res +  W(i)*W(j)*( 1-X(j) )*( temp1 + temp2 + temp3 - UEX( FK(1),FK(2) ) )**2
            end do
         end do
         error = error + det*res
      end do
      err = sqrt(error)      
      
end subroutine Compute_Error

end module EFLib
