! Created on Thu Jul 20 14:04:35 2018
! @author: SATGE Lancelot, KANYAMIBWA Romaric
! Resolution of PDE using Finite Difference Method
! ================
! Test of Finite Difference Numerical Schemes on sound waves equations
! """
! 
! #############################################################################
! #  Copyright (Const_c) 2017                                                 #
! #                                                                           #
! #                                                                           #
! #  Distributed under the terms of the GNU General Public License (GPL)      #
! #  either version 3, or (at your option) any later version                  #
! #                                                                           #
! #  http://www.gnu.org/licenses/                                             #
! #############################################################################


! Numerical Resolution of leading order system

! d_tt PH0 - 5.0/6 * d_xx PH0 = 0
! with P the pressure

! d_t tB0 - 1.0/2 *gamma2 d_yy tB0 = 0
! with tB0 the temperature

module Global_Var
    double precision, allocatable :: A(:,:),Atau(:,:)
    double precision, allocatable :: B(:),Btau(:)
    double precision, allocatable :: WH0(:),TH0(:),PB1(:)
    character(len=30) :: dir_name
    
    contains
        elemental subroutine strip(string,set)
            character(len=*), intent(inout) :: string
            character(len=*), intent(in)    :: set
            integer                         :: old, new, stride
            old = 1; new = 1
            do
                stride = scan( string( old : ), set )
                if ( stride > 0 ) then
                    string( new : new+stride-2 ) = string( old : old+stride-2 )
                    old = old+stride
                    new = new+stride-1
                else
                    string( new : ) = string( old : )
                    return
                end if
            end do
        end subroutine strip
end module Global_Var


program finite_diff

use Global_Var
! use omp_lib
implicit none

INTERFACE 
    FUNCTION inv (A)result(Ainv)
        double precision, dimension(:,:), intent(in) :: A
        double precision, dimension(size(A,1),size(A,2)) :: Ainv
    END FUNCTION inv
    
    function velocity (t) result(V)
        double precision,Intent(in):: t
        double precision ::V
    end function velocity
    
    function acceleration (t) result(accel)
        double precision,Intent(in):: t
        double precision ::accel
    end function acceleration
    
END INTERFACE

    ! type declaration statements
    integer::Nx,Nt,error!,n_images,image_periode
    integer::i,nt_i,n_images=1000,periode_images
    double precision ::PI=4.D0*DATAN(1.D0) ,acc
    double precision :: X_min,X_max,dx,Final_time,dt,dt_over_dx,Const_C
    double precision, allocatable :: X(:),PH0(:),next_PH0(:),before_PH0(:)
    double precision, allocatable :: Ainv(:,:),Atauinv(:,:),UH0(:),next_UH0(:)
    double precision, allocatable :: next_WH0(:)
    double precision :: alpha, beta
    double precision :: start, finish
    logical :: exist
    character(len=30)::fname="Generated_files/time_file.txt"
    character(len=32) :: arg
    integer :: nb_arg
    
    ! External procedures defined in BLAS
    external DGEMV

!     time discretization: 
!     du/dt(t^n) =  
!     1:  ( u^n+2 - 2*u^n+1 + u^n )/dt        (impl) 
!     2:  ( u^n+1 - 2*u^n +u^n-1 )/dt        (expl) 
!     t_choice = 2 
! 
!     space discretization: 
!     du/dx(x_i) =  
!     1:  ( u_i - 2*u_i-1 + u_i-2 )/dx 
!     2:  ( u_i+2 - 2*u_i+2 + u_i+2 )/dx 
!     3:  ( u_i+1 - 2*u_i-1 + u_i-3 )/2dx
!     x_choice = 2
!     call omp_set_num_threads(4)
!     print *,"Num_thds=", omp_get_num_threads()

    nb_arg=iargc()
    alpha=1.d0; beta=0.d0
    dir_name = 'simu_impliciteV1/'    
    ! spatial domain and mesh
    X_min = 0.0
    X_max = 100.0
    Nx= 1000

    ! temporal domain
    Final_time = 200*PI! t in [0,Final_time]
    Nt = 1000
    
    i=1
    do  while (i <= nb_arg)
        call getarg(i, arg)

        select case (arg)

            case ('-h', '--help')
                call print_help()
                stop
            case ('-t')
                i=i+1
                call getarg(i,arg)
                read(arg,*)Nt
                print*,"Nt=",Nt
            case ('-x')
                i=i+1
                call getarg(i,arg)
                read(arg,*)Nx
                print*,"Nx=",Nx
            case ('-d')
                i=i+1
                call getarg(i,arg)
                read(arg,*)X_max
                print*,"space interval S=[0,",X_max,"]"
            case ('--time')
                i=i+1
                call getarg(i,arg)
                read(arg,*)Final_time
                print*,"time interval T=[0,",Final_time,"]"
            case ('--nimages')
                i=i+1
                call getarg(i,arg)
                read(arg,*)n_images
                print*,"nimages=",n_images
            case default
                print '(a,a,/)', 'Unrecognized command-line option: ', arg
                call print_help()
                stop
        end select
        i=i+1
    end do
    
    !number of images generated
    periode_images = int(Nt*1./n_images)
    dx = (X_max-X_min)*1.0/Nx
    dt =  Final_time * 1.0/Nt
    dt_over_dx = (1.0*dt)/dx
    Const_C=-(5./6)*(dt_over_dx*dt_over_dx)
    
    CALL cpu_time(start)
    
    allocate(X(Nx+1),before_PH0(Nx+1),PH0(Nx+1),next_PH0(Nx+1),UH0(Nx+1),next_UH0(Nx+1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array X or PH0 or next or before or B, Nx+1=',Nx+1
        stop
    endif
    
    allocate(WH0(Nx+1),next_WH0(Nx+1),TH0(Nx+1),PB1(Nx+1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array WH0 or TH0 or next or before or B, Nx+1=',Nx+1
        stop
    endif
        
    !Initialisation
    !pressure
    PH0(:)=0.0
    next_PH0(:)=0.0
    before_PH0(:)=0.0
    
    PB1(:)=0.0
    
    !velocity
    UH0(:)=0.0
    next_UH0(:)=0.0
    
    !omega
    WH0(:)=0.0
    next_WH0(:)=0.0
    
    !temperature
    TH0(:)=0.0
    
!     executable statements  
!     print*,"Some results Nt",Nt
!     print*,"Some results Nx",Nx
    


    !$OMP PARALLEL DO 
    do i=1,Nx+1
        X(i) = X_min + (i-1)*dx
    end do
    !$OMP END PARALLEL DO 
    
    CALL constr_matrix_A(Nx,Const_C)
    CALL plot_sol(Nx,Nt,0,PH0,UH0,X)
    next_UH0(1)=velocity(dt)
    next_UH0(Nx+1)=0.0
    next_UH0(2:Nx)=UH0(2:Nx)-1.0/2.0*dt_over_dx*(next_PH0(2:Nx)-next_PH0(1:Nx-1))
    UH0(:)=next_UH0(:)
    CALL plot_sol(Nx,Nt,1,PH0,UH0,X)
    
!     do i=1,Nx-1
!         print*,real( A(i,:) ) 
!     end do
    
    allocate(B(Nx-1),Btau(Nx-1),Ainv(Nx-1,Nx-1),Atauinv(Nx-1,Nx-1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for vector B, Nx=',Nx
        stop
    endif
    
    Ainv=inv(A)
    do nt_i=2,Nt
!         print*,"----------------Nt=",nt_i,"----------------"
        acc=acceleration (nt_i*dt)
        CALL constr_vect_B(Const_C,dx,acc,Nx,before_PH0(2:Nx),PH0(2:Nx))
        
!         print*,"B vector"
!         do i=1,Nx-1
!             print*,real( B(i) ) 
!         end do
!         print*,"Ainv"
!         do i=1,Nx-1
!             print*,real( Ainv(i,:) ) 
!         end do
        
        CALL DGEMV('N',size(A,1),size(A,2),alpha,Ainv,size(A,2),B,1,beta,next_PH0(2:Nx),1)
        
        next_PH0(1)= next_PH0(2) + 2*dx*acc
        next_PH0(Nx+1)=next_PH0(Nx)
        
        next_UH0(1)=velocity(nt_i*dt)
        next_UH0(Nx+1)=0.0
        next_UH0(2:Nx)=UH0(2:Nx)-1.0/2.0*dt_over_dx*(next_PH0(2:Nx)-next_PH0(1:Nx-1))
        
        
        !omega and temperature
        
        next_WH0(:)=3.0/5.0*(next_PH0(:)-PH0(:))+WH0(:)
        TH0(:)=next_PH0(:)-next_WH0(:)        
    !     print*,"Next PH0"
    !     do i=1,Nx+1
    !         print*,real( next_PH0(i) ) 
    !     end do
    !     
    !     print*,"PH0"
    !     do i=1,Nx+1
    !         print*,real( PH0(i) ) 
    !     end do
        before_PH0(:)=PH0(:)
        PH0(:)=next_PH0(:)
        UH0(:)=next_UH0(:)
        
        if ( (MOD((nt_i),periode_images) == 0) .OR. ((nt_i) == Nt) ) then
            CALL plot_sol(Nx,Nt,nt_i,PH0,UH0,X)
        endif
        
    !     print*,"after aff PH0"
    !     do i=1,Nx+1
    !         print*,real( PH0(i) ) 
    !     end do
    !     
    !     print*,"before PH0"
    !     do i=1,Nx+1
    !         print*,real( before_PH0(i) ) 
    !     end do
    end do
    
    CALL cpu_time(finish)
    inquire(file=fname, exist=exist)
    
    if (exist) then
        open(125, file = fname, status='old',position="append")
    else
        open(125, file = fname, status='new')
    endif   
    write(125,*)'**********************************************************'
    write(125,*)'Nx =',Nx
    write(125,*)'Nt =',Nt
    write(125,*)'Nfiles =',n_images
    write(125,*)'Time =',finish-start,' seconds.'
    write(*,*)'Time =',finish-start,' seconds.'
    write(125,*)'**********************************************************'
    write(125,*)''
    close(125)
    
    !-------------------END OF program-------------------
    
    ! deallocation and end of program
    deallocate(X,PH0,next_PH0,before_PH0,TH0,PB1,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif
    
    deallocate(A,Ainv,B,WH0,next_WH0,UH0,next_UH0,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif

    
    !------------------- Print Help----------------------
    
    contains

    subroutine print_help()
        print  '(a)','usage: SoundWaves [OPTIONS]'
        print  '(a)',''
        print  '(a)','Without further options, SoundWaves initializes Nx=1000,Nt=1000,D=100 and T=200*pi,n_images=1000.'
        print  '(a)',''
        print  '(a)','cmdline options:'
        print  '(a)',''
        print  '(a)','  -x          Nx size'
        print  '(a)','  -t          Nt size'
        print  '(a)','  -d,         [0,D] space interval'
        print  '(a)','  --time,     [0,T] time interval'
        print  '(a)','  --nimages,  number of images'
        print  '(a)','  -h, --help  print usage information and exit'
    end subroutine print_help
end program finite_diff
