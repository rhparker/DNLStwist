
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
      DOUBLE PRECISION D,K,W,PHI,G

      ! variable split into magnitude and angle
      DOUBLE PRECISION A((NDIM+1)/2), P((NDIM+1)/2), S((NDIM+1)/2), C((NDIM+1)/2), PT((NDIM+1)/2)

      D = PAR(1)
      K = PAR(2)
      W = PAR(3)
      PHI = PAR(4)
      G = PAR(5)

      ! magnitudes
      N = (NDIM+1)/2
      DO I = 1,N
          A(I) = U(I)
      END DO
      ! phase angles, first one is always 0
      P(1) = 0
      DO I = 1,N-1
          P(I+1) = U(I+N)
      END DO
      ! terms for PT symmetry
      DO I = 1,N
            PT(I) = ( (-1)**(I+1) )
      END DO

      ! compute relevant sin/cos terms
      DO I = 1,N-1
          C(I) = COS( P(I+1) - P(I) - PHI )
          S(I) = SIN( P(I+1) - P(I) - PHI )
      END DO
      C(N) = COS( P(1) - P(N) - PHI )
      S(N) = SIN( P(1) - P(N) - PHI )

      ! full domain with periodic BCs     
      DO I = 2,N-1
          F(I)   = K*( C(I)*A(I+1) + C(I-1)*A(I-1) ) + W*A(I) + D*(A(I)**3) - G*PT(I)*SIN(P(I))*A(I)
          F(I+N) =   ( S(I)*A(I+1) - S(I-1)*A(I-1) ) + G*PT(I)*COS(P(I))*A(I)
      END DO
      ! periodic BCs
      F(1)   = K*( C(1)*A(2) + C(N)*A(N) ) + W*A(1) + D*(A(1)**3) - G*PT(1)*SIN(P(1))*A(1)
      F(1+N) =   ( S(1)*A(2) - S(N)*A(N) ) + G*PT(1)*COS(P(1))*A(1)
      F(N)   = K*( C(N)*A(1) + C(N-1)*A(N-1) ) + W*A(N) + D*(A(N)**3) - G*PT(N)*SIN(P(N))*A(N)
      F(N+N) =   ( S(N)*A(1) - S(N-1)*A(N-1) ) + G*PT(N)*COS(P(N))*A(N)

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
      DOUBLE PRECISION D, K, W, PHI, G

! Initialize the equation parameters
      D = -1
      K = 0
      W = 1
      PHI = 0
      G = 0

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
      PAR(5) = G

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