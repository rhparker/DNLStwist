
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
      INTEGER I
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION D,K,W,PHI

      D = PAR(1)
      K = PAR(2)
      W = PAR(3)
      PHI = PAR(4)

      F(1) = K*U(2) + W*U(1) + D*(U(1)**3)
      DO I = 2,NDIM-1
          F(I) = K*( U(I-1) + U(I+1) ) + W*U(I) + D*(U(I)**3)
      END DO
      F(NDIM) = K*(U(NDIM) + U(NDIM-1)) + W*U(NDIM) + D*(U(NDIM)**3)


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
      INTEGER I
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

      U(NDIM) = 1
      U(NDIM-1) = 1

      PAR(1) = D
      PAR(2) = K
      PAR(3) = W
      PAR(4) = PHI

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