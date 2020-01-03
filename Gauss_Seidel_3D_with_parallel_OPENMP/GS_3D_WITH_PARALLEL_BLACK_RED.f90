!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                 
!!!              MAYANK VADSOLA 300075960                                                   
!!!                     ASSIGNMENT 1                                                                     
!!!             QUESTION 2: 3D GAUSS SEIDEL WITH PARALLEL
!!!                 BLACK AND RED ALGORITHM OPENMP
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM GS
USE OMP_LIB

    IMPLICIT NONE
            
            INTEGER, PARAMETER :: N = 100                  ! NUMBER OF NODES IN X AND Y
            INTEGER, PARAMETER :: TOT = N * N * N         ! TOTAL NODES 
            REAL(KIND = 8), DIMENSION(0:TOT-1) :: T_NEW = 1, RESIDUAL = 1
            REAL(KIND = 8), PARAMETER :: PI = 3.1415927
            REAL(KIND = 8) :: X, Y, Z
            REAL(KIND = 8), DIMENSION(0:TOT-1) :: RHS
            REAL(KIND = 8) :: DELTA = 2.0/(N+1)
            INTEGER :: INDEXX, GLOBAL_INDEX
            REAL(KIND = 8) :: SOURCE, EXACT
            INTEGER :: I, J, K
            REAL(kIND = 8) :: T1, T2
            
            open (unit=2, file="EXACT.csv")
            DO I = 0, N-1
                X = (-1 + (I+1) * DELTA)
                DO J = 0, N-1
                    Y = (-1 + (J+1) * DELTA)
                    DO K = 0, N-1
                        Z = (-1 + (K+1) * DELTA)
                        INDEXX = GLOBAL_INDEX(I, J, K, N)
                        RHS(INDEXX) = DELTA * DELTA * SOURCE(X , Y, Z, PI)
                        write(2,200)  X, Y, Z , EXACT(X, Y, Z, PI)
                        200 format(f10.5 , ",", f10.5 , ",", F10.5 , ",", f10.5)
                    END DO
                END DO
            END DO
            
            T1 = OMP_GET_WTIME() 
            CALL GS_SOLVE(T_NEW, N, TOT, RHS, RESIDUAL)
            T2 = OMP_GET_WTIME()
            PRINT*, T2-T1
            
            open (unit=1, file="GS_3D_WITH_PARALLEL_BLACK_RED.csv")
            DO I = 0, N-1
                DO J = 0, N-1
                    DO K = 0, N-1
                        write(1,100)  (-1 + (I+1) * DELTA), (-1 + (J+1) * DELTA), &
                        (-1 + (K+1) * DELTA) , T_NEW(GLOBAL_INDEX(I, J, K, N))
                        100 format(f10.5 , ",", f10.5 , ",", F10.5 , ",", f10.5) 
                    END DO
                END DO
            END DO
            
            
            
END PROGRAM GS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     GLOBAL INDEX
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION GLOBAL_INDEX(I, J, K, N)

        IMPLICIT NONE
        INTEGER :: I, J, K, N
        INTEGER :: GI
        
        GI = I+N*J+N*N*K
        GLOBAL_INDEX = GI

END FUNCTION GLOBAL_INDEX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     SOURCE
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND= 8) FUNCTION SOURCE(X, Y, Z, PI)

        IMPLICIT NONE
        REAL(KIND = 8) :: X, Y, Z
        REAL(KIND = 8) :: PI
        REAL(KIND = 8) :: S
        
         S = -35.0 * PI * PI * 0.25 * COS(0.5 * PI * X) * COS(0.5 * PI * Y * 3.0) * COS(0.5 * PI * Z * 5.0)
        !S = -20.0 * PI * PI * SIN(4.0 * PI * X) * SIN(2.0 * PI * Y)
        !S = 0
        SOURCE = S
        
END FUNCTION SOURCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                         EXACT SOLUTION
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=8) FUNCTION EXACT(X, Y, Z, PI)

      IMPLICIT NONE
      REAL(KIND = 8) :: X, Y, Z
      REAL(KIND = 8) :: PI
      REAL(KIND = 8) :: EX
      
      EX = COS(0.5 * PI * X) * COS(0.5 * PI * Y * 3.0) * COS(0.5 * PI * Z * 5.0)
      EXACT = EX
      
END FUNCTION EXACT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     GS ITERATION
!!!                     T_NEW = (L)* T_OLD
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GS_ITERATION(T_NEW, RHS, N, TOT)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        INTEGER :: I, J, K
        INTEGER, INTENT(IN) :: TOT
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(IN) :: RHS
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(INOUT) :: T_NEW
        INTEGER :: GLOBAL_INDEX, GII

        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 0, N-2, 2
            DO J = 0, N-1
                DO I = 0, N-2, 2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    T_NEW(GII) = -RHS(GII);
                    
                    IF(I > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                    
                    T_NEW(GII) = T_NEW(GII)/6.0
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 1, N-1,2
            DO J = 0, N-1
                DO I = 1, N-1,2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    T_NEW(GII) = -RHS(GII);
                    
                    IF(I > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                    
                    T_NEW(GII) = T_NEW(GII)/6.0
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 0, N-2, 2
            DO J = 0, N-1
                DO I = 1, N-1, 2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    T_NEW(GII) = -RHS(GII);
                    
                    IF(I > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                    
                    T_NEW(GII) = T_NEW(GII)/6.0
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 1, N-1,2
            DO J = 0, N-1
                DO I = 0, N-2,2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    T_NEW(GII) = -RHS(GII);
                    
                    IF(I > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) T_NEW(GII) = T_NEW(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                    
                    T_NEW(GII) = T_NEW(GII)/6.0
                END DO
            END DO
        END DO
        !$OMP END DO
                
END SUBROUTINE GS_ITERATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     RESIDUAL
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GS_RESIDUAL(T_NEW, RESIDUAL, NORM_RESIDUAL, RHS, N, TOT)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        INTEGER :: I, J, K
        INTEGER, INTENT(IN) :: TOT
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(IN) :: RHS
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(INOUT) :: T_NEW
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(INOUT) :: RESIDUAL 
        REAL(KIND = 8), INTENT(OUT) :: NORM_RESIDUAL
        INTEGER :: GLOBAL_INDEX, GII
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 0, N-2, 2
            DO J = 0, N-1
                DO I = 0, N-2, 2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    RESIDUAL(GII) = -6.0*T_NEW(GII) - RHS(GII);
        
                    IF(I > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 1, N-1,2
            DO J = 0, N-1
                DO I = 1, N-1,2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    RESIDUAL(GII) = -6.0*T_NEW(GII) - RHS(GII);
                    
                    IF(I > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 0, N-2, 2
            DO J = 0, N-1
                DO I = 1, N-1, 2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    RESIDUAL(GII) = -6.0*T_NEW(GII) - RHS(GII);
                    
                    IF(I > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                END DO
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC,10)
        DO K = 1, N-1,2
            DO J = 0, N-1
                DO I = 0, N-2,2
                    GII = GLOBAL_INDEX(I, J, K, N)
                    RESIDUAL(GII) = -6.0*T_NEW(GII) - RHS(GII);
                    
                    IF(I > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I-1, J, K, N))
                    IF(I < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I+1, J, K, N))
                    IF(J > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J-1, K, N))
                    IF(J < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J+1, K, N))
                    IF(K > 0)   RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K-1, N))
                    IF(K < N-1) RESIDUAL(GII) = RESIDUAL(GII) + T_NEW(GLOBAL_INDEX(I, J, K+1, N))
                END DO
            END DO
        END DO
        !$OMP END DO
        
        NORM_RESIDUAL = NORM2(RESIDUAL)

END SUBROUTINE GS_RESIDUAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                         GS SOLVER
!!!                   -(D+U) T^{N+1} = (L) T^{N} - RHS
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GS_SOLVE(T_NEW, N, TOT, RHS, RESIDUAL)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: TOT
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(INOUT) :: T_NEW
        REAL(KIND = 8), DIMENSION(0:TOT-1), INTENT(INOUT) :: RHS
        REAL(KIND = 8), PARAMETER :: TOLERANCE = 1E-8
        REAL(KIND = 8), DIMENSION(0:TOT-1) :: RESIDUAL
        REAL(KIND = 8) :: NORM_RESIDUAL = 1
        INTEGER :: COUNTER = 0
        INTEGER, PARAMETER :: MAX_COUNTER = 10000000
        
        DO WHILE(NORM_RESIDUAL > TOLERANCE)
               
                !$OMP PARALLEL 
                CALL GS_ITERATION(T_NEW, RHS, N, TOT)
                CALL GS_RESIDUAL(T_NEW, RESIDUAL, NORM_RESIDUAL, RHS, N, TOT)
                !$OMP END PARALLEL
                
                !PRINT *, "NORM_RESIDUAL \n", NORM_RESIDUAL
                
                COUNTER = COUNTER + 1
                
                IF(COUNTER > MAX_COUNTER) THEN
                    PRINT* , "MAX_COUNTER EXCEEDED"
                    EXIT
                END IF
        END DO
        
        PRINT *, "COUNTER = ", COUNTER
        PRINT *, "RESIDUAL = ", NORM_RESIDUAL
        
END SUBROUTINE GS_SOLVE


        





















