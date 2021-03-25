
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      INTEGER I,N
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION D,K,W,PHI,K1,K2,K3,K4,K5,K6

      ! variable split into magnitude and angle
      DOUBLE PRECISION A((NDIM+1)/2), P((NDIM+1)/2), S((NDIM+1)/2), C((NDIM+1)/2)

      D   = PAR(1)
      K   = PAR(2)
      W   = PAR(3)
      PHI = PAR(4)
      K1  = PAR(5)
      K2  = PAR(6)
      K3  = PAR(7)
      K4  = PAR(8)
      K5  = PAR(9)
      K6  = PAR(10)

      ! magnitudes
      N = NDIM/2
      DO I = 1,N
          A(I) = U(I)
      END DO
      ! phase angles
      DO I = 1,N
          P(I) = U(I+N)
      END DO

      ! compute relevant sin/cos terms
      DO I = 1,N-1
          C(I) = COS( P(I+1) - P(I) - PHI )
          S(I) = SIN( P(I+1) - P(I) - PHI )
      END DO
      C(N) = COS( P(1) - P(N) - PHI )
      S(N) = SIN( P(1) - P(N) - PHI )

      F(1)   = (K+K1)*C(1)*A(2) + (K+K6)*C(N)*A(N) + W*A(1) + D*(A(1)**3)
      F(1+N) = (K+K1)*S(1)*A(2) - (K+K6)*S(N)*A(N)

      F(2)   = (K+K2)*C(2)*A(3) + (K+K1)*C(1)*A(1) + W*A(2) + D*(A(2)**3)
      F(2+N) = (K+K2)*S(2)*A(3) - (K+K1)*S(1)*A(1)

      F(3)   = (K+K3)*C(3)*A(4) + (K+K2)*C(2)*A(2) + W*A(3) + D*(A(3)**3)
      F(3+N) = (K+K3)*S(3)*A(4) - (K+K2)*S(2)*A(2)

      F(4)   = (K+K4)*C(4)*A(5) + (K+K3)*C(3)*A(3) + W*A(4) + D*(A(4)**3)
      F(4+N) = (K+K4)*S(4)*A(5) - (K+K3)*S(3)*A(3)

      F(5)   = (K+K5)*C(5)*A(6) + (K+K4)*C(4)*A(4) + W*A(5) + D*(A(5)**3)
      F(5+N) = (K+K5)*S(5)*A(6) - (K+K4)*S(4)*A(4)

      F(N)   = (K+K6)*C(N)*A(1) + (K+K5)*C(N-1)*A(N-1) + W*A(N) + D*(A(N)**3)
      F(N+N) = (K+K6)*S(N)*A(1) - (K+K5)*S(N-1)*A(N-1)

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
      INTEGER NDIM
      INTEGER LOFFSET, ROFFSET
      INTEGER I, N, C
      DOUBLE PRECISION U(NDIM), PAR(*), T
      DOUBLE PRECISION D, K, W, PHI

! Initialize the equation parameters
      D = -1
      K = 0
      W = 1
      PHI = 0

      ! initialize to 0
      DO I = 1,NDIM
          U(I) = 0
      END DO

      ! N = (NDIM+1)/2
      ! ! find center
      ! IF ( MOD(N, 2) == 0 ) THEN
      !       C = N/2
      ! ELSE
      !       C = (N+1)/2
      ! END IF
      ! U(C) = 1
      ! U(C+1) = 1

      U(1) = 1

      PAR(1) = D
      PAR(2) = K
      PAR(3) = W
      PAR(4) = PHI

      DO I= 5,10
          PAR(I) = 0
      END DO

END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------