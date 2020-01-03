program MBS
USE OMP_LIB

        integer, parameter :: real_size = 1920 
        integer, parameter :: imag_size = 1080

        real(kind = 16), parameter :: real_max = -1.7499576837060933110499019727D0
        real(kind = 16), parameter :: imag_max = 1.10787937065633798823892515031e-16
        real(kind = 16), parameter :: real_min = -1.7499576837060936950743114676D0 
        real(kind = 16), parameter :: imag_min =-1.05212062934366210875555371911e-16
        real(kind = 16), parameter :: real_mid = -1.7499576837060935030621067198D0 
        real(kind = 16), parameter :: imag_mid = 2.7879370656337939741685715600e-18
        real(kind = 16), parameter :: real_dx = (real_max-real_min)/(real_size)
        real(kind = 16), parameter :: imag_dy = (imag_max-imag_min)/(imag_size)
        integer :: counter = 0
        complex(kind = 16) :: c
        complex(kind = 16) :: z
        complex(kind = 16), parameter :: dx=(real_dx ,0.0)
        complex(kind = 16), parameter :: dy=(0.0, imag_dy) 
        complex(kind = 16), parameter :: c_mid = (real_mid, imag_mid)
        real(kind = 16) :: mag = 0 
        real(kind = 16) :: limit = 2.0
        integer, dimension (real_size, imag_size) :: iteration
        REAL(kIND = 8) :: T1, T2

        T1 = OMP_GET_WTIME() 
!$omp parallel private(c, z, mag, counter)
!$omp do schedule(dynamic)
        do i = 1,real_size
            do j = 1,imag_size
                c = cmplx(real_min+i*real_dx,imag_min+j*imag_dy,kind = 16) 
                z = (0.0, 0.0)
                mag = 0.0
                counter = 0
                do while (mag<limit)
                    z = z*z + c
                    mag = abs(z)
                    if (counter > 1000000) then
                        counter = 1000000
                        exit
                    end if 
                    counter = counter + 1
                end do
                
                iteration(i,j) = counter
            end do
        end do
!$omp end do
!$omp end parallel  
        T2 = OMP_GET_WTIME()
        PRINT*, T2-T1

        open (unit=1, file="MandelBrotSet_2D_parallel_OPENMP.csv")
        do i = 1, real_size
            do j = 1, imag_size              
                write(1,100) (real_min + real_dx*(i) - real_mid)*1E18, & 
                (imag_min + imag_dy*(j)-imag_mid)*1E18, 0 , iteration(i, j)
                100 format(f20.14 , ",", f22.16 , ",", i1 , ",", i7) 
             end do
        end do
        
        
        
end program MBS


        
        
