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

! subroutine for .dat files creation
subroutine plot_sol(Nx,Nt,n,PH0,UH0,X,Y)
use Global_Var
implicit none
    double precision,Intent(IN)::PH0(Nx+1),UH0(Nx+1)
    double precision,Intent(IN)::X(Nx+1),Y(Nx+1)
    integer,Intent(In)::n,Nx,Nt
    character(len=100)::fname,tmpstr='',strn
    logical :: exist
    integer :: i

    write (tmpstr,*) trim("out"),n,"_Nx",Nx,"_Nt",Nt
    write (strn,*)n
    CALL strip(strn,' ')
    fname=trim(dir_name)
    CALL strip(tmpstr,' ')
!     write(*,*)trim(tmpstr),"hahah"
    fname=trim(fname)//trim(tmpstr)//trim("_PH0.dat")
!     write(*,*)fname
    
    inquire(file=fname, exist=exist)
    
    if (exist) then
        open(1, file = fname, status='old')
    else
        open(1, file = fname, status='new')
    endif   
    
    print*,"Plot sol in file ",trim(fname),&
    trim(", nt = "),trim(strn),trim(",PH0 min/max = "),&
    minval(PH0),trim("/"), maxval(PH0)
    
    print*,trim("---TB0 min/max = "), minval(TB0),trim("/"), maxval(TB0),"---"
    if(order == 0) then
    
        do i=1,Nx+1
            write(1,*)X(i)," ",PH0(i)," ",UH0(i)," ",WH0(i)," ",TH0(i) &
            ," ",Y(i)," ",TB0(i)," ",WB0(i)
        enddo
        
    else if( order == 1) then
    
        do i=1,Nx+1
            write(1,*)X(i)," ",PH0(i)," ",UH0(i)," ",WH0(i)," ",TH0(i) &
            ," ",Y(i)," ",TB0(i)," ",WB0(i)," ",UB1(i)
        enddo
        print*,trim("---UB1 min/max = "), minval(UB1),trim("/"), maxval(UB1),"---"
        
    endif
    
    close(1)
    
end subroutine plot_sol
