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

!First order asymptotic expansion

module First_Order
    !First Order macroscopic quantities
    double precision, allocatable :: UB1(:),PH1(:),next_UB1(:)
    double precision, allocatable :: next_PH1(:),before_PH1(:)
    double precision, allocatable :: UH1(:),next_UH1(:),WH1(:),TH1(:)
    double precision, allocatable :: next_WH1(:)!,next_TB1(:),next_WB1(:)
!     double precision, allocatable :: WB0(:),TB0(:)
end module First_Order

!Initialisation of First Order quantities
subroutine InitFirstOrderQuant(Nx,Ny)
    
!     use Global_Var
    use First_Order
    implicit none
    
    integer, intent(in) ::Nx,Ny
    integer :: error
    
    allocate(UB1(Ny+1),before_PH1(Nx+1),PH1(Nx+1),next_PH1(Nx+1)&
    ,next_UB1(Ny+1),UH1(Nx+1),next_UH1(Nx+1)&
    ,WH1(Nx+1),next_WH1(Nx+1),TH1(Nx+1)&
    ,stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array WH0 or TH0 or next or before or B, Nx+1=',Nx+1
        stop
    endif

    !Initialisation
    !boundary velocity
    UB1(:)=0.0
    next_UB1(:)=0.0
    
    !pressure
    PH1(:)=0.0
    next_PH1(:)=0.0
    before_PH1(:)=0.0
    
!     velocity
    UH1(:)=0.0
    next_UH1(:)=0.0
    
    !omega
    WH1(:)=0.0
    next_WH1(:)=0.0
!     WB1(:)=0.0
!     next_WB1(:)=0.0
    
    !temperature
    TH1(:)=0.0
!     TB1(:)=0.0
!     next_TB1(:)=0.0
    
end subroutine InitFirstOrderQuant

! subroutine UpdateQuant()
!     use Global_Var
!     use First_Order
! 
! end subroutine UpdateQuant


subroutine FirstOrderDealloc()    

    use First_Order
    implicit none
    
    integer::error
    
    deallocate(UB1,PH1,before_PH1,next_PH1&
    ,next_UB1&
    ,UH1,next_UH1,WH1,next_WH1,TH1&
    ,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif


end subroutine FirstOrderDealloc
    
