!In Fortan for a matrix: dimension(row_nb,col_nb)::matrix!
!For browsering the matrix : matrix(row_nb,col_nb)!

subroutine file_writing(dx,Nt,Nx_P,P,U,W,T)
	implicit none
	real :: dx
	integer :: i,j,k
    integer:: Nt,Nx_P
	real,dimension(Nx_P,Nt) :: P,U,W,T
	character(len=80) file_name

	do i=1,5000
		call system('clear')		
		print *,"### File Creation ###"
		print *,i*100/(Nt+1),"%"
		print *,"Nx =",Nx_P,"Nt =",Nt-1
		if (i<10) then		
	  		write(file_name,'(A,I4.1,A,I4.1,A,I1.1,A)') "file_txt/P_Nx_",Nx_P,"_Nt_",Nt-1,"_img_",i,".txt"
		else if (i<100) then
			write(file_name,'(A,I4.1,A,I4.1,A,I2.1,A)') "file_txt/P_Nx_",Nx_P,"_Nt_",Nt-1,"_img_",i,".txt"
		else if (i<1000) then
			write(file_name,'(A,I4.1,A,I4.1,A,I3.1,A)') "file_txt/P_Nx_",Nx_P,"_Nt_",Nt-1,"_img_",i,".txt"
		else if (i<10000) then
			write(file_name,'(A,I4.1,A,I4.1,A,I4.1,A)') "file_txt/P_Nx_",Nx_P,"_Nt_",Nt-1,"_img_",i,".txt"
		endif

	 	open(11,file=file_name,form="formatted",access="sequential")

	 	do j=1,Nx_P
	 		write(11, fmt="(f23.20)", advance="no")(j-1)*dx
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="no")P(j,i)	 		
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="no")U(j,i)	 			 		
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="no")W(j,i)	 			 		
	 		write(11, fmt="(A)", advance="no")" "
	 		write(11, fmt="(f23.20)", advance="yes")T(j,i)
	 	end do
	 	
	 	close(11)
	end do

end subroutine file_writing

program test
	implicit none

	integer::i,j
	integer:: Nx_A,Nt,Nx_P
	real::D,dx,dt,c
	real::start,end
	character(len=32) :: Nx_P_str,Nt_str
	real,allocatable :: A_ini(:,:),A_tri(:,:),P_final(:,:),P_temp(:,:),B(:),U(:,:),W(:,:),T(:,:)
	
	!Reading of arguments
	call getarg(1, Nx_P_str)
	call getarg(2, Nt_str)  
	read(Nx_P_str,'(i10)') Nx_P
	read(Nt_str,'(i10)') Nt
	Nx_A=Nx_P-2

	allocate(A_ini(Nx_A,Nx_A),A_tri(Nx_A,Nx_A),P_final(Nx_P,Nt+1),P_temp(Nx_A,Nt+1),B(Nx_A),U(Nx_P,Nt+1),W(Nx_P,Nt+1),T(Nx_P,Nt+1))


	A_ini(:,:)=0.0
	A_tri(:,:)=0.0
	P_final(:,:)=0.0	
	P_temp(:,:)=0.0
	B(:)=0.0
	U(:,:)=0.0
	W(:,:)=0.0
	T(:,:)=0.0

	

	D=100.0
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

			!Transformation of B for triangular meTd
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
		T(:,i)=P_final(:,i)-W(:,i)
	end do

	call cpu_time(start)
	call file_writing(dx,Nt+1,Nx_P,P_final,U,W,T)
	call cpu_time(end)

	print *,end-start,"s"

end program test



























