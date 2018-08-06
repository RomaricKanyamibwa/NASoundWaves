subroutine file_writing(dx,Nt,Nx_T,T,TH0)
	implicit none
	real :: dx
	integer :: i,j
    integer:: Nt,Nx_T
	real,dimension(Nx_T,Nt) :: T,TH0
	character(len=80) file_name

	do i=1,Nt,Nt/1000+1
		call system('clear')		
		print *,"### File Creation ###"
		print *,i*100/(Nt+1),"%"
		print *,"Nx =",Nx_T,"Nt =",Nt-1
		if (i<10) then		
	  		write(file_name,'(A,I1.1,A)') "file_txt/TB0_img_",i,".txt"
		else if (i<100) then
			write(file_name,'(A,I2.1,A)') "file_txt/TB0_img_",i,".txt"
		else if (i<1000) then
			write(file_name,'(A,I3.1,A)') "file_txt/TB0_img_",i,".txt"
		else if (i<10000) then
			write(file_name,'(A,I4.1,A)') "file_txt/TB0_img_",i,".txt"
		endif

	 	open(11,file=file_name,form="formatted",access="sequential")

	 	do j=1,Nx_T
	 		write(11, fmt="(f23.20)", advance="no")(j-1)*dx	 
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="no")TH0(j,i)				 		
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="yes")T(j,i)
	 	end do
	 	
	 	close(11)
	end do

end subroutine file_writing


subroutine print_matrix(n,m,mat)
    implicit none
    integer :: n,m,i,j
    real :: mat(n,m)

    do i=1,n
    	do j=1,m
    		write(*, fmt="(f30.20)", advance="no")mat(i,j)
    	end do
    	print *,
    end do
    
end subroutine print_matrix

program one
	implicit none
    integer::i,j
	integer:: Nx_A,Nt,Nx_T,Nx_P
	real::D,dx,dt,c,gamma2
	real::start,end
	character(len=32) :: Nx_T_str,Nt_str
	real,allocatable :: A_ini(:,:),A_tri(:,:),T_final(:,:),T_temp(:,:),B(:),P_final(:,:),P_temp(:,:),U(:,:),W(:,:),THO(:,:)

	!Reading of arguments
	call getarg(1, Nx_T_str)
	call getarg(2, Nt_str)  
	read(Nx_T_str,'(i10)') Nx_T
	read(Nt_str,'(i10)') Nt
	Nx_A=Nx_T-2
	Nx_P=Nx_T

	allocate(A_ini(Nx_A,Nx_A),A_tri(Nx_A,Nx_A),T_final(Nx_T,Nt+1),T_temp(Nx_A,Nt+1),B(Nx_A))
	allocate(P_final(Nx_P,Nt+1),P_temp(Nx_A,Nt+1),U(Nx_P,Nt+1),W(Nx_P,Nt+1),THO(Nx_P,Nt+1))

	!First part to find THO
	D=20.0
	dx=D/real(Nx_P-1)
	dt=200*3.1415/real(Nt)
	c=-5*dt*dt/(6*dx*dx)

			!Initialization of P!
	!--------------------------------!
	!|	P(0,0)  P(0,1) .. P(0,Nt) 	|!
	!|	P(1,0) 	P(1,1) .. P(1,Nt)	|!
	!|	.........................	|!
	!|	.........................	|!
	!|	P(Nx,0)  P(Nx,1) ..P(Nx,Nt) |!
	!--------------------------------!

	!Initialization of A, with a triangular transformation to make easy the resolution of A*x=B!
	A_ini(1,1)=1-c
	A_ini(1,2)=c
	A_ini(Nx_A,Nx_A-1)=c
	A_ini(Nx_A,Nx_A)=1-c
	do i=2,Nx_A-1
		A_ini(i,i)=1-2*c
		A_ini(i,i-1)=c
		A_ini(i,i+1)=c
		
	end do

	A_tri=A_ini
	A_tri(Nx_A-1,Nx_A-1)=A_tri(Nx_A-1,Nx_A-1)-A_tri(Nx_A,Nx_A-1)*A_tri(Nx_A-1,Nx_A)/A_tri(Nx_A,Nx_A)
	A_tri(Nx_A-1,Nx_A)=0	
	do i=Nx_A-2,1,-1
		A_tri(i,i)=A_tri(i,i)-A_tri(i+1,i)*A_tri(i,i+1)/A_tri(i+1,i+1)
		A_tri(i,i+1)=0
	end do

	!Calculating of each Pn!
	do i=1,Nt+1
		if (i>2) then
			
			B(:)=0.0
			B(1)=2*P_temp(1,i-1)-P_temp(1,i-2)-2*c*dx*(-cos((i-1)*dt)+(1-(i-1)*dt)*exp(-(i-1)*dt))

			!Initialization of B
			do j=2,Nx_A
				B(j)=2*P_temp(j,i-1)-P_temp(j,i-2)
			end do

			!Transformation of B for triangular methode
			do j=Nx_A-1,1,-1
				B(j)=B(j)-B(j+1)*A_ini(j,j+1)/A_tri(j+1,j+1)
			end do

			!System Resolving
			P_temp(1,i)=B(1)/A_tri(1,1)
			do j=2,Nx_A
				P_temp(j,i)=(B(j)-A_tri(j,j-1)*P_temp(j-1,i))/A_tri(j,j)
			end do

			P_final(2:Nx_A+1,i)=P_temp(:,i)
			P_final(1,i)=P_final(2,i)+2*dx*(-cos((i-1)*dt)+(1-(i-1)*dt)*exp(-(i-1)*dt))
			P_final(Nx_P,i)=P_final(Nx_P-1,i)
			!print *,"n=",i,"maxval=",maxval(P_final(:,i))

		end if
	end do

	!Calculating of U
	do i=1,Nt+1
		if (i>2) then
			U(1,i)=-sin((i-1)*dt)+(i-1)*dt*exp(-(i-1)*dt)
			U(2:Nx_P-1,i)=U(2:Nx_P-1,i-1)+dt*(P_final(3:Nx_P,i)-P_final(2:Nx_P-1,i))/(2.0*dx)

		endif
	end do

	!Calculating of W
	do i=2,Nt+1
		W(:,i)=(3.0/5.0)*(P_final(:,i)-P_final(:,i-1))+W(:,i-1)
	end do

	!Calculating of T
	do i=1,Nt+1
		THO(:,i)=P_final(:,i)-W(:,i)
	end do

	!Calculating T_final (TB0)
	A_ini(:,:)=0.0
	A_tri(:,:)=0.0
	T_final(:,:)=0.0	
	T_temp(:,:)=0.0
	B(:)=0.0
	

	D=100.0
	gamma2=1.922284066
	dx=D/real(Nx_T-1)
	dt=200*3.1415/real(Nt)
	c=gamma2*dt/(2*dx*dx)

			!Initialization of T!
	!--------------------------------!
	!|	T(0,0)  T(0,1) .. T(0,Nt) 	|!
	!|	T(1,0) 	T(1,1) .. T(1,Nt)	|!
	!|	.........................	|!
	!|	.........................	|!
	!|	T(Nx,0)  T(Nx,1) ..T(Nx,Nt) |!
	!--------------------------------!

	!Initialization of A, with a triangular transformation to make easy the resolution of A*x=B!
	A_ini(1,1)=1+2*c
	A_ini(1,2)=-c
	A_ini(Nx_A,Nx_A-1)=-c
	A_ini(Nx_A,Nx_A)=1+2*c
	do i=2,Nx_A-1
		A_ini(i,i)=1+2*c
		A_ini(i,i-1)=-c
		A_ini(i,i+1)=-c
		
	end do

	A_tri=A_ini
	A_tri(Nx_A-1,Nx_A-1)=A_tri(Nx_A-1,Nx_A-1)-A_tri(Nx_A,Nx_A-1)*A_tri(Nx_A-1,Nx_A)/A_tri(Nx_A,Nx_A)
	A_tri(Nx_A-1,Nx_A)=0	
	do i=Nx_A-2,1,-1
		A_tri(i,i)=A_tri(i,i)-A_tri(i+1,i)*A_tri(i,i+1)/A_tri(i+1,i+1)
		A_tri(i,i+1)=0
	end do

	!Calculating of Bn to get Pn+1=A^-1*Bn
	do i=2,Nt+1
		B(1)=T_temp(1,i-1)-c*THO(2,i)

		!Calculating B
		do j=2,Nx_A
			B(j)=T_temp(j,i-1)
		end do

		! print *,"i=",i
		! print *,"B ini"
		! print *,B

		!Transformation of B for triangular methode
		do j=Nx_A-1,1,-1
			B(j)=B(j)-B(j+1)*A_ini(j,j+1)/A_tri(j+1,j+1)
		end do

		! print *,"B tri"
		! print *,B

		
		!System Resolving
		T_temp(1,i)=B(1)/A_tri(1,1)
		do j=2,Nx_A
			T_temp(j,i)=(B(j)-A_tri(j,j-1)*T_temp(j-1,i))/A_tri(j,j)
		end do

		! print *,"T_temp"
		! print *,T_temp(1:Nx_A,i)

		! print *,"THO"
		! print *,THO(1:Nx_T,i)

		T_final(2:Nx_A+1,i)=T_temp(:,i)
		T_final(1,i)=-THO(1,i)

	end do

	call file_writing(dx,Nt+1,Nx_T,T_final,THO)

end program one


