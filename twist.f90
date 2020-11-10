SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      INTEGER I, N
      DOUBLE PRECISION D,K,PHI,P
      DOUBLE PRECISION SP, CP
      DOUBLE PRECISION A(NDIM/2), B(NDIM/2)

      D   = PAR(1)
      K   = PAR(2)
      PHI = PAR(3)
      P   = PAR(4)

      CP = COS(PHI)
      SP = SIN(PHI)
      N = NDIM/2

      DO I = 1, N
          A(I) = U(I)
          B(I) = U(N+I)
      END DO

      DO I = 2, N-1
          F(I)   = K*( SP*(-A(I+1)+A(I-1)) + CP*(B(I+1)+B(I-1) ) ) +D*(A(I)**2 + B(I)**2)*B(I)
          F(N+I) = K*( CP*( A(I+1)+A(I-1)) + SP*(B(I+1)-B(I-1) ) ) +D*(A(I)**2 + B(I)**2)*A(I)
      END DO

      ! Periodic BCs
      F(1)   = K*( SP*(-A(2)+A(N)) + CP*(B(2)+B(N) ) ) +D*(A(1)**2 + B(1)**2)*B(1)
      F(N+1) = K*( CP*( A(2)+A(N)) + SP*(B(2)-B(N) ) ) +D*(A(1)**2 + B(1)**2)*A(1)

      F(N)   = K*( SP*(-A(1)+A(N-1)) + CP*(B(1)+B(N-1) ) ) +D*(A(N)**2 + B(N)**2)*B(N)
      F(N+N) = K*( CP*( A(1)+A(N-1)) + SP*(B(1)-B(N-1) ) ) +D*(A(N)**2 + B(N)**2)*A(N)

      ! scale by period P to get interval to [0,1]
      DO I = 1, NDIM
          F(I) = P*F(I)
      END DO

END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      INTEGER I, N
      DOUBLE PRECISION PI
      DOUBLE PRECISION D, K, PHI, P

! Initialize the equation parameters
      PI = 4 * ATAN(1.0_16);
      D   = -1
      K   = 0
      PHI = 0
      P   = 2*PI
      N   = NDIM/2

      DO I = 1,N
          U(I)   =  COS(2*PI*T)
          U(N+I) =  SIN(2*PI*T)
      END DO

      PAR(1) = D
      PAR(2) = K
      PAR(3) = PHI
      PAR(4) = P

      END SUBROUTINE STPNT

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
      INTEGER I

      ! periodic BCs
      DO I = 1,NDIM
          FB(I) = U0(I) - U1(I)
      END DO

END SUBROUTINE BCND

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------