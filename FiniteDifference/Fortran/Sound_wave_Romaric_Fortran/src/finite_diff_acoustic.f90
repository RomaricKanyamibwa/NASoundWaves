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


! We want to resolve numericaly the following PDE
! d_tt PH0 - 5.0/6 * d_xx PH0 = 0
! with P the pressure

function acceleration (t) result(accel)
implicit none
    
    double precision,Intent(in):: t
    double precision ::accel
    
    accel=-cos(t)+(1-t)*exp(-t)

end function acceleration


function velocity (t) result(V)
    implicit none
    double precision,Intent(in):: t
    double precision ::V
    
    V=-sin(t)+(t)*exp(-t)
end function velocity

subroutine constr_matrix_A(Nx,Const_C)
use Global_Var
implicit none

    double precision,Intent(IN):: Const_C
    integer,intent(IN) :: Nx
    integer ::i,error
!     double precision ,allocatable ,Intent(OUT) :: A(:,:)
    
    allocate(A(Nx-1,Nx-1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array A, Nx=',Nx
        stop
    endif
    A(:,:)=0.0
    A(1,1)=1-Const_C
    A(1,2)=Const_C
    
    !$OMP PARALLEL DO 
    do i=2,Nx-2
        A(i,i-1)=Const_C
        A(i,i)=1-2*Const_C
        A(i,i+1)=Const_C
    enddo
    !$OMP END PARALLEL DO
    
    A(Nx-1,Nx-2)=Const_C
    A(Nx-1,Nx-1)=1-Const_C

end subroutine constr_matrix_A

subroutine constr_vect_B(Const_C,dx,acc,Nx,Pn1,Pn2)
use Global_Var
implicit none

    double precision,Intent(IN):: Const_C,dx,acc
    integer,Intent(IN)::Nx
    double precision,Intent(IN)::Pn1(Nx-1),Pn2(Nx-1)
    integer ::i
    
    B(1)=2*Pn2(1)-Pn1(1)-2*Const_C*dx*acc
    !$OMP PARALLEL DO 
    do i=2,Nx-1
        B(i)=2*Pn2(i)-Pn1(i)
    enddo
    !$OMP END PARALLEL DO
    
end subroutine constr_vect_B
