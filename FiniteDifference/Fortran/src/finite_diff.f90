! Created on Thu Jul 20 14:04:35 2018
! @author: SATGE Lancelot, KANYAMIBWA Romaric
! Resolution of PDE using Finite Difference
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


! We want to resolve numericaly the following PDE
! d_tt PH0 - 5.0/6 * d_xx PH0 = 0
! with P the pressure



program finite_diff
implicit none      

    ! type declaration statements
    integer::Nx,Nt,error!,n_images,image_periode
    integer::i,j
    character(len=30) :: dir_name
    double precision ::PI=4.D0*DATAN(1.D0) ,acc,acceleration
    double precision :: X_min,X_max,dx,Final_time,dt,dt_over_dx,Const_C
    double precision, allocatable :: X(:),PH0(:),next_PH0(:),before_PH0(:),B(:)
    double precision, allocatable :: A(:,:)

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

    dir_name = 'simu_impliciteV1/'    
    ! spatial domain and mesh
    X_min = 0.0
    X_max = 5.0
    Nx=6
    dx = 1.0/Nx
    Final_time = 200*PI! t in [0,Final_time]
    Nt = 1000
    dt = Final_time * 1./Nt
    dt_over_dx = dt/dx

    ! temporal domain
    Final_time = 1
    dt =  Final_time * 1./Nt
    dt_over_dx = 1.0*dt/dx
    Const_C=-(5./6)*(dt_over_dx*dt_over_dx)
    

    allocate(X(Nx+1),B(Nx-1),before_PH0(Nx+1),PH0(Nx+1),next_PH0(Nx+1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array X or PH0 or next or before or B, Nx+1=',Nx+1
        stop
    endif
    
    allocate(A(Nx-1,Nx-1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array A, Nx-1=',Nx+1
        stop
    endif
        
    X(:)=0.0
    PH0(:)=0.0
    next_PH0(:)=0.0
    before_PH0(:)=0.0
!     executable statements  
!     print*,"Some results Nt",Nt
!     print*,"Some results Nx",Nx
    
!     do j=1,Nx+1
!         print*,"",PH0
!     end do
!     
!     do j=1,Nx+1
!         print*,"",X
!     end do
!     acc=acceleration (nt)

    do i=1,Nx+1
        X(i) = X_min + i*dx
    end do
    
    ! deallocation and end of program
    deallocate(X,PH0,next_PH0,before_PH0,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif
    
    deallocate(A,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif

end program finite_diff


function acceleration (t)
implicit none
    
    
    double precision::acceleration
    double precision :: t
    
    acceleration=-cos(t)+(1-t)*exp(-t)

end function acceleration

subroutine constr_matrix_A(Nx,Const_C,A)
    double precision,Intent(IN):: Nx,Const_C
    integer ::i,j
    double precision ,allocatable,Intent(OUT) :: A(:,:)

end subroutine constr_matrix_A

