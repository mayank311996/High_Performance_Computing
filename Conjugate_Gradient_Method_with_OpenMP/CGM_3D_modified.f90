!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                 
!!!              MAYANK VADSOLA 300075960                                                   
!!!                     ASSIGNMENT 3                                                                     
!!!             CGM: CONJUGATE GRADIENT METHOD
!!!                         3-D
!!!                        OpenMP
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! PLEASE COMMENT OUT ALL PRINT STATEMENTS TO CHECK PARALLEL EFFICIENCY!!!!!!!!!!

PROGRAM CGM_3D
USE OMP_LIB

        IMPLICIT NONE 
        
                INTEGER, PARAMETER :: N = 200                   ! NODES PER DIRECTION
                INTEGER, PARAMETER :: TOT = N * N * N           ! TOTAL NUMBER OF NODES 
                !REAL(KIND=8), DIMENSION(0:TOT-1) :: T_FINAL = 0.0D0 ! FINAL TEMP ARRAY
                REAL(KIND = 8), ALLOCATABLE :: T_FINAL(:)
                REAL(KIND = 8), ALLOCATABLE :: RHS(:)
                REAL(KIND=8), PARAMETER :: PI = 4 * atan(1.0_8) ! PI 
                REAL(KIND=8) :: X, Y, Z                         ! X & Y COORDINATES
                !REAL(KIND=8), DIMENSION(0:TOT-1) :: RHS         ! RHS
                REAL(KIND=8) :: DELTA = 2.0D0/(N-1)               ! STEP SIZE
                INTEGER :: INDEXX, GLOBAL_INDEX                 ! INDEX
                REAL(KIND=8) :: SOURCE, EXACT                   ! SOURCE AND EXACT SOLUTION ARRAY
                INTEGER :: I, J, K                              ! LOOP INDEX
                REAL(kIND = 8) :: T1, T2
                
                ALLOCATE(T_FINAL(0:TOT-1))
                ALLOCATE(RHS(0:TOT-1))
                T_FINAL = 0.0D0
                
                open (unit=2, file="CGM_3D_EXACT.csv")
                DO I = 0, N-1
                    X = (-1 + (I) * DELTA)
                    DO J = 0, N-1
                        Y = (-1 + (J) * DELTA)
                        DO K = 0, N-1
                            Z = (-1 + (K) * DELTA)
                            INDEXX = GLOBAL_INDEX(I, J, K, N)
                            RHS(INDEXX) =  SOURCE(X , Y, Z, PI) * DELTA * DELTA
                            write(2,200)  (-1 + (I) * DELTA), (-1 + (J) * DELTA), & 
                                (-1 + (K) * DELTA) , EXACT(X, Y, Z, PI)
                            200 format(f10.5 , ",", f10.5 , ",", f10.5 , ",", f10.5) 
                        END DO
                    END DO
                END DO

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                T1 = OMP_GET_WTIME() 
                CALL CGM_SOLVE(T_FINAL, N, TOT, RHS) !!!!!!!!!!!
                T2 = OMP_GET_WTIME()
                PRINT*, T2-T1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                open (unit=1, file="CGM_3D_WITH_PARALLEL.csv")
                DO I = 0, N-1
                    DO J = 0, N-1
                        DO K = 0, N-1
                            write(1,100)  (-1 + (I) * DELTA), (-1 + (J) * DELTA), &
                            (-1 + (K) * DELTA) , T_FINAL(GLOBAL_INDEX(I, J, K, N))
                            100 format(f10.5 , ",", f10.5 , ",", F10.5 , ",", f10.5) 
                        END DO
                    END DO
                END DO

END PROGRAM CGM_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     GLOBAL INDEX
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION GLOBAL_INDEX(I, J, K, N)

        !IMPLICIT NONE
        INTEGER :: I, J, K, N
        INTEGER :: GI
        
        GI = I+N*J+N*N*K
        GLOBAL_INDEX = GI

END FUNCTION GLOBAL_INDEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     SOURCE
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(KIND=8) FUNCTION SOURCE(X, Y, Z, PI)

        !IMPLICIT NONE
        REAL(KIND=8) :: X, Y, Z
        REAL(KIND=8) :: PI
        REAL(KIND=8) :: S
        
        S = -35.0D0 * PI * PI * 0.25 * COS(0.5 * PI * X) *  & 
            COS(0.5 * PI * Y * 3.0) * COS(0.5 * PI * Z * 5.0)
        !S = 0
        SOURCE = S
        
END FUNCTION SOURCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                         EXACT SOLUTION
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(KIND=8) FUNCTION EXACT(X, Y, Z, PI)

      !IMPLICIT NONE
      REAL(KIND=8) :: X, Y, Z
      REAL(KIND=8) :: PI
      REAL(KIND=8) :: EX
      
      EX = COS(0.5 * PI * X) * COS(0.5 * PI * Y * 3.0) * COS(0.5 * PI * Z * 5.0)
      EXACT = EX
      
END FUNCTION EXACT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                         CGM SOLVER
!!!                   
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CGM_SOLVE(T_FINAL, N, TOT, RHS)

        IMPLICIT NONE
        
                INTEGER :: TOT
                INTEGER :: N
                REAL(KIND = 8), DIMENSION(0:TOT-1) :: T_FINAL
                REAL(KIND = 8), DIMENSION(0:TOT-1) :: RHS
                REAL(KIND = 8), ALLOCATABLE :: T_OLD(:), T_NEW(:) 
                REAL(KIND = 8), ALLOCATABLE :: R_OLD(:), R_NEW(:)
                REAL(KIND = 8), ALLOCATABLE :: D_OLD(:), D_NEW(:)
                REAL(KIND = 8), ALLOCATABLE :: A_D(:)
                REAL(KIND = 8) :: BETA, ALPHA
                REAL(KIND = 8) :: R_TRA_R = 0.0D0
                REAL(KIND = 8) :: D_TRA_A_D = 0.0D0
                REAL(KIND = 8) :: R_TRA_R_NEW = 0.0D0
                REAL(KIND = 16) :: RESIDUAL = 1.0D0
                
                REAL(KIND = 16) :: TOLERANCE = 1.0E-10
                INTEGER :: AA = -6, BB = 1
                INTEGER :: GLOBAL_INDEX, GII
                INTEGER :: COUNTER = 0
                INTEGER :: I, J, K
                
                ALLOCATE(T_OLD(0:TOT-1)) 
                ALLOCATE(T_NEW(0:TOT-1))
                ALLOCATE(R_OLD(0:TOT-1)) 
                ALLOCATE(R_NEW(0:TOT-1))
                ALLOCATE(A_D(0:TOT-1))
                ALLOCATE(D_OLD(0:TOT-1))
                ALLOCATE(D_NEW(0:TOT-1))
                T_OLD = 0.0D0
                T_NEW = 0.0D0
                
                !!!!!!!!!!!!!!!!! SETTING UP INNER DOMAIN DIFFERENT THEN BOUNDARIES !!!!!!!!!!!!!!!
                !$OMP PARALLEL
                !$OMP DO SCHEDULE(STATIC)
                DO K = 1, N-2
                    DO J = 1, N-2
                        DO I = 1, N-2
                            
                        T_OLD(GLOBAL_INDEX(I, J, K, N)) = 10.0  
                            
                        END DO
                    END DO
                END DO
                !$OMP END DO
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                !!!!!!!!!!!!!!!!!!!!!!!!! CALCULATING OLD RESIDUAL AND OLD DIRECTION !!!!!!!!!!!!!!!
                !$OMP DO SCHEDULE(STATIC) PRIVATE(GII) REDUCTION(+:R_TRA_R)
                DO K = 1, N-2
                    DO J = 1, N-2
                        DO I = 1, N-2
                            
                            GII = GLOBAL_INDEX(I, J, K, N)
                            R_OLD(GII) = RHS(GII) - AA * T_OLD(GII) & 
                                        - BB * T_OLD(GLOBAL_INDEX(I-1, J, K, N)) &
                                        - BB * T_OLD(GLOBAL_INDEX(I+1, J, K, N)) &
                                        - BB * T_OLD(GLOBAL_INDEX(I, J-1, K, N)) &
                                        - BB * T_OLD(GLOBAL_INDEX(I, J+1, K, N)) & 
                                        - BB * T_OLD(GLOBAL_INDEX(I, J, K-1, N)) &
                                        - BB * T_OLD(GLOBAL_INDEX(I, J, K+1, N))
                            
                            
                            D_OLD(GII) = R_OLD(GII)
                            
                            R_TRA_R = R_TRA_R + (R_OLD(GII))**2
                            
                        END DO
                    END DO
                END DO
                !$OMP END DO
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                DO WHILE (RESIDUAL > TOLERANCE)
                
                            !$OMP SINGLE
                            R_TRA_R_NEW = 0
                            D_TRA_A_D = 0
                            COUNTER = COUNTER + 1
                            !$OMP END SINGLE
                            
               !!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULATING A*D AND D(TRANSPOSE)*A*D !!!!!!!!!!!!!!!!!!!!!!!!             
                             !$OMP DO SCHEDULE(STATIC) PRIVATE(GII) REDUCTION(+:D_TRA_A_D)
                             DO K = 1, N-2
                                DO J = 1, N-2
                                    DO I = 1, N-2
                                        
                                        GII = GLOBAL_INDEX(I, J, K, N)
                                        A_D(GII) = AA * D_OLD(GII)   & 
                                                   + BB * D_OLD(GLOBAL_INDEX(I-1, J, K, N)) &
                                                   + BB * D_OLD(GLOBAL_INDEX(I+1, J, K, N)) &
                                                   + BB * D_OLD(GLOBAL_INDEX(I, J-1, K, N)) &
                                                   + BB * D_OLD(GLOBAL_INDEX(I, J+1, K, N)) &
                                                   + BB * D_OLD(GLOBAL_INDEX(I, J, K-1, N)) &
                                                   + BB * D_OLD(GLOBAL_INDEX(I, J, K+1, N))
                                        
                                        D_TRA_A_D = D_TRA_A_D + D_OLD(GII) * A_D(GII)
                                        
                                    END DO
                                END DO
                            END DO
                            !$OMP END DO
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            
                            !$OMP SINGLE
                            ALPHA = R_TRA_R / D_TRA_A_D  !!!!!!!!!!!! ALPHA
                            !$OMP END SINGLE
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UPDATING THE TEMPERATURE !!!!!!!!!!!!!!!!!!!!!            
                            !$OMP DO SCHEDULE(STATIC) PRIVATE(GII)
                            DO I = 0, TOT-1
                            
                                T_NEW(I) = T_OLD(I) + ALPHA * D_OLD(I)
                            
                            END DO
                            !$OMP END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !!!!!!!!!!!!!!!!!!!!!!!!!! UPDATING THE RESIDUAL AND R(TRANSPOSE)*R !!!!!!!!!!!!!!!!!!!                
                            !$OMP DO SCHEDULE(STATIC) PRIVATE(GII) REDUCTION(+:R_TRA_R_NEW)
                            DO I = 0, TOT-1
                            
                                R_NEW(I) = R_OLD(I) - ALPHA * A_D(I)
                                R_TRA_R_NEW = R_TRA_R_NEW + (R_NEW(I))**2 
                            
                            END DO
                            !$OMP END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            
                            !$OMP SINGLE
                            BETA = R_TRA_R_NEW / R_TRA_R  !!!!!!!!!! BETA
                            !$OMP END SINGLE
             
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULATING NEW DIRECTION !!!!!!!!!!!!!!!!!!!!                
                            !$OMP DO SCHEDULE(STATIC) PRIVATE(GII)
                            DO I = 0, TOT-1
                            
                                D_NEW(I) = R_NEW(I) + BETA * D_OLD(I)
                            
                            END DO
                            !$OMP END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            
                            !$OMP SINGLE
                            RESIDUAL = NORM2(R_NEW) !NORM2
                            
                            PRINT*, RESIDUAL
                            
                            R_OLD = R_NEW
                            R_TRA_R=R_TRA_R_NEW
                            D_OLD = D_NEW
                            T_OLD = T_NEW
                            !$OMP END SINGLE
                
                END DO
                !$OMP END PARALLEL
                
                PRINT*, COUNTER
            
            T_FINAL = T_NEW
    
END SUBROUTINE CGM_SOLVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
