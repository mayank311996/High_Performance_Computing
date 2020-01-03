!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                 
!!!              MAYANK VADSOLA 300075960                                                   
!!!                     ASSIGNMENT 4                                                                     
!!!             FFT: FAST FOURIER TRANSFORM
!!!                         
!!!                        OpenMP
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM FFT
USE OMP_LIB

        IMPLICIT NONE
        
        INTEGER(KIND = 8), PARAMETER :: N = 2**20            ! MULTIPLE OF 2
        REAL(KIND = 8), PARAMETER :: PI = 4 * atan(1.0_8) ! PI 
        COMPLEX(KIND = 8), ALLOCATABLE :: X_K(:)
        INTEGER(KIND = 8) :: I
        REAL(kIND = 8) :: T1, T2
        
        ALLOCATE(X_K(N))
        
        T1 = OMP_GET_WTIME() 
        
        !$OMP PARALLEL DO SCHEDULE(STATIC)
        DO I = 1, N
        
            X_K(I) = EXP(CMPLX(0.0D0, 2.0D0 * PI * REAL(7660087,8) * REAL(I-1,8) / REAL(13,8) / REAL(N,8), KIND = 8))
        
        END DO
        !$OMP END PARALLEL DO
        
        
        CALL FFT_SUB(X_K)
        
        T2 = OMP_GET_WTIME()
        PRINT*, T2-T1
        
        open (unit=2, file="FFT.csv")
        DO I = 1, N
        
            write(2,200)  X_K(I)
            200 format(f15.5 , ",", f15.5)
        
        END DO
        
        CONTAINS 
        
            RECURSIVE SUBROUTINE FFT_SUB(XX)

                    IMPLICIT NONE
                    
                        REAL(KIND = 8), PARAMETER :: PI = 4 * atan(1.0_8) ! PI
                        INTEGER(KIND = 8) :: NN
                        COMPLEX(KIND = 8), DIMENSION(:), INTENT(INOUT) :: XX
                        COMPLEX(KIND = 8) :: TEMP
                        INTEGER(KIND = 8) :: I
                        COMPLEX(KIND = 8), ALLOCATABLE :: EVEN(:), ODD(:)
                        
                        NN = SIZE(XX)
                        
                        IF(NN .le. 1) RETURN
                        
                        ALLOCATE(ODD((NN)/2))
                        ALLOCATE(EVEN(NN/2))
                        
                        
                        
                        ODD = XX(1:NN:2)
                        EVEN = XX(2:NN:2)
                        
                        CALL FFT_SUB(ODD)
                        CALL FFT_SUB(EVEN)
                        
                        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(TEMP) IF(NN>2)
                        DO I = 1, NN/2 
                        
                            TEMP  = EXP(CMPLX(0.0D0, -2.0D0 * PI * REAL(I-1,8) / REAL(NN,8), KIND = 8))
                            
                            XX(I)      = ODD(I) + TEMP * EVEN(I) 
                            XX(I+NN/2) = ODD(I) - TEMP * EVEN(I)
                        
                        END DO
                        !$OMP END PARALLEL DO
                        
                        DEALLOCATE(ODD)
                        DEALLOCATE(EVEN)

            END SUBROUTINE FFT_SUB
                

END PROGRAM FFT


