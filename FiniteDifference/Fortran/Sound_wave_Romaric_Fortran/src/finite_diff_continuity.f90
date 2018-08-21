! Created on Thu aug 9 11:43:20 2018
! @author: SATGE Lancelot, KANYAMIBWA Romaric
! Resolution of PDE using Finite Difference Method
! ================
! Finite Difference Numerical Schemes on continuity equation
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


! Numerical Resolution of the following PDE
! d_y uB = - d_t wB0 
! ubB1 is the flow velocity and wB0 is the density

subroutine first_orderUB1(Ny,dy_over_dt,next_WB0)
use Global_Var

double precision,Intent(In)::dy_over_dt,next_WB0(Ny+1)
integer,Intent(In)::Ny
integer::i

! print*,"*******************************************"
UB1(Ny+1)=0
do i=1,Ny
    UB1(Ny+1-i)=UB1(Ny+2-i)+dy_over_dt*(next_WB0(Ny+2-i)-WB0(Ny+2-i))
!     print*,(next_WB0(Ny+2-i)-WB0(Ny+2-i))
end do
! print*,"*******************************************"

end subroutine first_orderUB1
