
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
      DOUBLE PRECISION D,K,W

      D = PAR(1)
      K = PAR(2)
      W = PAR(3)

      ! F(1) = K*( U(2) + U(5) ) + W*U(1)    + D*(U(1)**3)
      ! F(2) = K*( U(3) + U(1) ) + W*U(2)    + D*(U(2)**3)
      ! F(3) = K*( U(4) + U(2) ) + W*U(3)    + D*(U(3)**3)
      ! F(4) = K*( U(5) + U(3) ) + W*U(4)    + D*(U(4)**3)
      ! F(5) = K*( U(1) + U(4) ) + W*U(5)    + D*(U(5)**3)

      ! hole, N=5
      ! F(1) = K*( U(2) -  U(1) )  + W*U(1)    +  D*(U(1)**3)
      ! F(2) = K*( U(1) ) + W*U(2) + D*(U(2)**3)

      ! double hole, N=6
      ! F(1) = K*( U(2) ) + W*U(1) + D*(U(1)**3)
      ! F(2) = K*( U(1) ) + W*U(2) + D*(U(2)**3)

      ! hole, N=3
      ! F(1) = K*( -U(1) +  U(2) )  + W*U(1)    +  D*(U(1)**3)
      ! F(2) = K*( U(1)  +  U(3) )  + W*U(2)    +  D*(U(2)**3)
      ! F(3) = K*( U(2) )           + W*U(3) + D*(U(3)**3)

      ! ! even grid, two equally spaced holes
      ! F(1) = 0
      ! DO I = 2,NDIM-1
      !    F(I) = K * (U(I+1) + U(I-1)) + W*U(I) + D*(U(I)**3)  
      ! END DO 
      ! F(NDIM) = 0

      ! ! odd grid, hole in center
      ! F(1) = K * (-U(1) + U(2)) + W*U(1) + D*(U(1)**3)
      ! DO I = 2,NDIM-1
      !    F(I) = K * (U(I+1) + U(I-1)) + W*U(I) + D*(U(I)**3)  
      ! END DO 
      ! F(NDIM) = U(NDIM)

      ! full domain with periodic BCs     
      DO I = 2,NDIM-1
            F(I) = K*( U(I+1) + U(I-1) ) + W*U(I) + D*(U(I)**3)
      END DO
      F(1)    = K*( U(2) + U(NDIM)   ) + W*U(1)    + D*(U(1)**3)
      F(NDIM) = K*( U(1) + U(NDIM-1) ) + W*U(NDIM) + D*(U(NDIM)**3)

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
      DOUBLE PRECISION D, K, W, A

! Initialize the equation parameters
      D = -1
      K = 0
      W = 1

      ! initialize to 0
      
      DO I = 1,NDIM
          U(I) = 0
      END DO

      ! initialize center to 1
      IF ( MOD(NDIM, 2) == 0 ) THEN
            U(NDIM/2) = 1
      ELSE
            U( (NDIM+1)/2 ) = 1
      END IF

      ! ! initialize first element to 1
      ! U(1) = 1

      PAR(1) = D
      PAR(2) = K
      PAR(3) = W

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