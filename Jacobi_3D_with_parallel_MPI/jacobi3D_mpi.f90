!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                 
!!!                                                                 
!!!                     ASSIGNMENT 1                                                                     
!!!                     QUESTION 2: 3D JACOBI WITH PARALLEL MPI 
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM JACOBI

    INCLUDE 'mpif.h'
    
        !IMPLICIT NONE
        INTEGER, PARAMETER :: N = 50                   ! NUMBER OF NODES IN EACH DIRECTION
        INTEGER, PARAMETER :: TOT = N * N * N           ! TOTAL NODES
        REAL(KIND=8), DIMENSION(0:(N-2)*(N-2)*(N)-1) :: T_FINAL = 0 ! FINAL TEMP ARRAY
        REAL(KIND=8), PARAMETER :: PI = 4 * atan(1.0_8) ! PI 
        REAL(KIND=8) :: X, Y, Z                         ! X & Y COORDINATES
        REAL(KIND=8), DIMENSION(0:TOT-1) :: RHS         ! RHS
        REAL(KIND=8) :: DELTA = 2.0/(N-1)               ! STEP SIZE
        INTEGER :: INDEXX, GLOBAL_INDEX                 ! INDEX
        REAL(KIND=8) :: SOURCE, EXACT                   ! SOURCE AND EXACT SOLUTION ARRAY
        INTEGER :: I, J, K                              ! LOOP INDEX
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! MPI VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER :: IERROR                               ! MPI ERROR
        INTEGER :: NP                                   ! NUM OF PROCESS
        INTEGER :: PID                                  ! PROCESS ID
        REAL(KIND=8) :: T_START, T_STOP, T_DIFF         ! TIME VARIABLES
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        open (unit=2, file="JACOBI_3D_EXACT.csv")
        DO I = 0, N-1
            X = (-1 + (I) * DELTA)
            DO J = 0, N-1
                Y = (-1 + (J) * DELTA)
                DO K = 0, N-1
                    Z = (-1 + (K) * DELTA)
                    INDEXX = GLOBAL_INDEX(I, J, K, N)
                    RHS(INDEXX) = DELTA * DELTA * SOURCE(X , Y, Z, PI)
                    write(2,200)  (-1 + (I+1) * DELTA), (-1 + (J+1) * DELTA), & 
                        (-1 + (K+1) * DELTA) , EXACT(X, Y, Z, PI)
                    200 format(f10.5 , ",", f10.5 , ",", f10.5 , ",", f10.5) 
                END DO
            END DO
        END DO
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL MPI_INIT(IERROR)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERROR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, PID, IERROR)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        T_START = MPI_WTIME()
        
        CALL JACOBI_SOLVE(T_FINAL, N, TOT, RHS, NP, PID)
        
        T_STOP = MPI_WTIME()
        IF(PID == 0) THEN
            PRINT*, "TIME = ", T_STOP-T_START
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
        CALL MPI_FINALIZE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM JACOBI

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
        
        S = -35.0 * PI * PI * 0.25 * COS(0.5 * PI * X) *  & 
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
!!!                         JACOBI SOLVER
!!!                   -D T^{N+1} = (L+U) T^{N} - RHS
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE JACOBI_SOLVE(T_FINAL, N, TOT, RHS, NP, PID)

    INCLUDE 'mpif.h'
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! MPI VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER :: IERROR                               ! MPI ERROR
        INTEGER :: NP                                   ! NUM OF PROCESS
        INTEGER :: PID                                  ! PROCESS ID
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
        !IMPLICIT NONE
        REAL(KIND=8), ALLOCATABLE :: T_OLD(:), T_NEW(:) 
        INTEGER :: N
        INTEGER, INTENT(IN) :: TOT
        REAL(KIND=8), DIMENSION(0:(N-2)*(N-2)*(N)-1), INTENT(INOUT) :: T_FINAL
        REAL(KIND=8), DIMENSION(0:TOT-1), INTENT(INOUT) :: RHS
        REAL(KIND=8), PARAMETER :: TOLERANCE = 1E-8
        REAL(KIND=8) :: RESIDUAL_GLOBAL = 1
        REAL(KIND=8) :: RESIDUAL_LOCAL = 0
        INTEGER :: COUNTER = 0
        INTEGER, PARAMETER :: MAX_COUNTER = 10000000
        INTEGER :: GLOBAL_INDEX 
        
        INTEGER :: CELLS 
        CELLS = INT(N*N*(2 + N/NP))

        ALLOCATE(T_OLD(0:CELLS-1)) 
        ALLOCATE(T_NEW(0:CELLS-1))
        T_OLD = 0
        T_NEW = 0
        
        
        !PRINT*, PID
        DO WHILE(RESIDUAL_GLOBAL > TOLERANCE)
            
            T_OLD = T_NEW
            CALL JACOBI_HELPER(T_NEW, T_OLD, N, TOT, CELLS, RHS, NP, PID)
            
            RESIDUAL_LOCAL = NORM2(T_NEW - T_OLD)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL MPI_ALLREDUCE(RESIDUAL_LOCAL, RESIDUAL_GLOBAL, 1, & 
                            MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, IERROR)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            !PRINT*, PID
            IF(PID == 0) THEN
                PRINT*, RESIDUAL_GLOBAL
            END IF
            
            COUNTER = COUNTER + 1
            
            IF(COUNTER > MAX_COUNTER) THEN
                PRINT* , "MAX_COUNTER EXCEEDED"
                EXIT
            END IF
            
        END DO
        
        
        IF(PID == 0) THEN
            PRINT *, "COUNTER = ", COUNTER
            PRINT *, "RESIDUAL = ", RESIDUAL_GLOBAL
        END IF
        
        !CALL MPI_GATHER(T_NEW(GLOBAL_INDEX(1,1,1,N):GLOBAL_INDEX(N-2,N-2,N/NP,N)), ((N-2)*(N-2)*(N/NP)), & 
        !                MPI_REAL8, T_FINAL, ((N-2)*(N-2)*(N/NP)), MPI_REAL8, 0, MPI_COMM_WORLD, IERROR)
                        
        !!!!!!!!!!!!!!!!!!! OUTPUTING DATA GENERATED BY Oth PROCESSOR FOR SIMPLE VISUALIZATION !!!!!!!!!!!!!!!!!!!!                
        IF(PID == 0) THEN
            open (unit=1, file="JACOBI_3D_WITH_PARALLEL.csv")
            DO I = 0, N-1
                DO J = 0, N-1
                    DO K = 0, N/NP+1
                        write(1,100)  ((I)), ((J)), & 
                        ((K)) , T_NEW(GLOBAL_INDEX(I, J, K, N))
                        100 format(I5 , ",", I5 , ",", I5 , ",", f10.5) 
                    END DO
                END DO
            END DO
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE JACOBI_SOLVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                     JACOBI HELPER
!!!                     T_NEW = (L+U) T_OLD
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE JACOBI_HELPER(T_NEW, T_OLD, N, TOT, CELLS, RHS, NP, PID)

    INCLUDE 'mpif.h'
    
        !IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: TOT
        INTEGER :: CELLS 
        REAL(KIND=8), DIMENSION(0:CELLS-1) :: T_OLD, T_NEW 
        REAL(KIND=8), DIMENSION(0:TOT-1), INTENT(IN) :: RHS
        INTEGER :: GLOBAL_INDEX
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! MPI VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER :: IERROR                               ! MPI ERROR
        INTEGER :: NP                                   ! NUM OF PROCESS
        INTEGER :: PID                                  ! PROCESS ID
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO K = 1, (N/NP)
            DO J = 1, N-2
                DO I = 1, N-2
                    
                    T_NEW(GLOBAL_INDEX(I, J, K, N)) = (T_OLD(GLOBAL_INDEX(I-1, J, K, N)) &
                                                        + T_OLD(GLOBAL_INDEX(I+1, J, K, N)) &
                                                        + T_OLD(GLOBAL_INDEX(I, J-1, K, N)) &
                                                        + T_OLD(GLOBAL_INDEX(I, J+1, K, N)) &
                                                        + T_OLD(GLOBAL_INDEX(I, J, K-1, N)) &
                                                        + T_OLD(GLOBAL_INDEX(I, J, K+1, N)))/6.0
                                                        
                    T_NEW(GLOBAL_INDEX(I, J, K, N)) = T_NEW(GLOBAL_INDEX(I, J, K, N)) & 
                                                        - (RHS(GLOBAL_INDEX(I, J, K+(PID*N/NP), N)))/6.0
                END DO
            END DO
        END DO
        
        IF(PID == 0) THEN
            CALL MPI_ISEND(T_NEW(GLOBAL_INDEX(1,1,(N/NP),N):GLOBAL_INDEX(N-2,N-2,(N/NP),N)), (N-2)*(N-2), &
                            MPI_REAL8, PID+1, PID, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
    
            CALL MPI_RECV(T_NEW(GLOBAL_INDEX(1,1,(N/NP)+1,N):GLOBAL_INDEX(N-2,N-2,(N/NP)+1,N)), (N-2)*(N-2), & 
                            MPI_REAL8, PID+1, PID+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
        END IF
                                                
        IF(PID>0 .AND. PID<NP-1) THEN
            CALL MPI_ISEND(T_NEW(GLOBAL_INDEX(1,1,1,N):GLOBAL_INDEX(N-2,N-2,1,N)), (N-2)*(N-2), &
                            MPI_REAL8, PID-1, PID, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
            CALL MPI_ISEND(T_NEW(GLOBAL_INDEX(1,1,(N/NP),N):GLOBAL_INDEX(N-2,N-2,(N/NP),N)), (N-2)*(N-2), &
                            MPI_REAL8, PID+1, PID, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
                            
            CALL MPI_RECV(T_NEW(GLOBAL_INDEX(1,1,(N/NP)+1,N):GLOBAL_INDEX(N-2,N-2,(N/NP)+1,N)), (N-2)*(N-2), & 
                            MPI_REAL8, PID+1, PID+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
            CALL MPI_RECV(T_NEW(GLOBAL_INDEX(1,1,0,N):GLOBAL_INDEX(N-2,N-2,0,N)), (N-2)*(N-2), & 
                            MPI_REAL8, PID-1, PID-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
        END IF
                                                
        IF(PID == NP-1) THEN
            CALL MPI_ISEND(T_NEW(GLOBAL_INDEX(1,1,1,N):GLOBAL_INDEX(N-2,N-2,1,N)), (N-2)*(N-2), &
                            MPI_REAL8, PID-1, PID, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
                            
            CALL MPI_RECV(T_NEW(GLOBAL_INDEX(1,1,0,N):GLOBAL_INDEX(N-2,N-2,0,N)), (N-2)*(N-2), & 
                            MPI_REAL8, PID-1, PID-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
        END IF
        
END SUBROUTINE JACOBI_HELPER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















