! Created on Thu aug 2 16:40:35 2018
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
! d_t tB0 - 1.0/2 *gamma2 d_yy tB0 = 0
! with t the temparature


subroutine constr_matrix_Atau(Ny,Const_D)
use Global_Var
implicit none

    double precision,Intent(IN):: Const_D
    integer,intent(IN) :: Ny
    integer ::i,error
    
    allocate(Atau(Ny-1,Ny-1),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array Atau, Ny=',Ny
        stop
    endif
    Atau(:,:)=0.0
    Atau(1,1)=1+2*Const_D
    Atau(1,2)=-Const_D
    
    !$OMP PARALLEL DO 
    do i=2,Ny-2
        Atau(i,i-1)=-Const_D
        Atau(i,i)=1+2*Const_D
        Atau(i,i+1)=-Const_D
    enddo
    !$OMP END PARALLEL DO
    
    Atau(Ny-1,Ny-2)=-Const_D
    Atau(Ny-1,Ny-1)=1+2*Const_D

end subroutine constr_matrix_Atau

subroutine constr_vect_Btau(Const_D,tauh0,Ny,Tn)
use Global_Var
implicit none

    double precision,Intent(IN):: Const_D,tauh0
    integer,Intent(IN)::Ny
    double precision,Intent(IN)::Tn(Ny-1)
    integer ::i
    
    Btau(1)=Tn(1)-Const_D*tauh0
    !$OMP PARALLEL DO 
    do i=2,Ny-1
        Btau(i)=Tn(i)
    enddo
    !$OMP END PARALLEL DO
    
end subroutine constr_vect_Btau
