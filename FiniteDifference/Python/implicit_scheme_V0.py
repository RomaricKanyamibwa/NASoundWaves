import numpy
import pprint
from numpy.linalg import inv
import os
import matplotlib.pyplot as plt

# spatial area and mesh
X_min = 0
X_max = 1
Nx = 10
dx = 1./Nx

# time domain
Temps_final = 1
Nt = 100
dt = 1./Nt

#Initialization
U = numpy.zeros((4*Nx,Nt)) #Global matrix (u1,u2,u3,u4)

#	---------------------------------
#	|	u1(0,0) u1(1,0) .. u1(n,0) 	|
#	|	u1(0,1) u1(1,1) .. u1(n,1)	|
#	|	.........................	|
#	|	.........................	|
#	|	u1(0,k) u1(1,k) .. u1(n,k)	|
#	|	u2(0,0) u2(1,0) .. u2(n,0) 	|
#	|	u2(0,1) u2(1,1) .. u2(n,1)	|
#	|	.........................	|
#	|	.........................	|
#	|	u2(0,k) u2(1,k) .. u2(n,k)	|
#	|	u3(0,0) u3(1,0) .. u2(n,0) 	|
#	|	u3(0,1) u3(1,1) .. u3(n,1)	|
#	|	.........................	|
#	|	.........................	|
#	|	u3(0,k) u3(1,k) .. u3(n,k)	|
#	|	P(0,0)  P(1,0) 	.. P(n,0) 	|
#	|	P(0,1) 	P(1,1) 	.. P(n,1)	|
#	|	.........................	|
#	|	.........................	|
#	|	P(0,k)  P(1,k)  .. P(n,k)	|
#	---------------------------------

# - Browsing the matrix - #
for n in range(0,Nt):
	for k in range(0,4*Nx):
		U[k,n]=0.0005*(n+1)/(k+1)

#print "Array non-initialized yet, coefficients show different parts of it "
#print U

# - Variables - #
a = round(5*dt*dt/(6*dx*dx),3)
c = round(-5*dt/(3*dx),3)
b = round(2*c/(3*dx),3)
d = round(-dt/(2*dx),3)

# - Matrix initialization - #
M1 = numpy.zeros((4*Nx,4*Nx))
for i in range(0,4*Nx):
	for j in range(0,4*Nx):
		if i==j:
			M1[i,j]=1.

		elif i<3*Nx:
			M1[i,i%Nx+3*Nx]=d

		elif i>=3*Nx and j==Nx*(j/Nx)+i%(3*Nx):
			M1[i,j]=c

M2 = numpy.zeros((4*Nx,4*Nx))
for i in range(0,4*Nx):
	for j in range(0,4*Nx):
		if i==j and i/Nx<3:
			M2[i,j]=1.-a

		elif j==i+1 and i/Nx<3:
			M2[i,j]=1

		elif j==i+Nx and i/Nx<2:
			M2[i,j]=-1			

		elif j==i+2*Nx and i/Nx<1:
			M2[i,j]=-1			

		elif j==i+Nx+1 and i/Nx<2:
			M2[i,j]=1			

		elif j==i+2*Nx+1 and i/Nx<1:
			M2[i,j]=1

		elif j==i-Nx and i/Nx<3:
			M2[i,j]=-1

		elif j==i-2*Nx and i/Nx<3:
			M2[i,j]=-1

		elif j==i-Nx+1 and i/Nx<3:
			M2[i,j]=1

		elif j==i-2*Nx+1 and i/Nx<3:
			M2[i,j]=1

		elif j>=3*Nx and j==i+(3-i/Nx)*Nx+1 and i/Nx<3:
			M2[i,j]=d

		if (i+1)%Nx==0 and j==Nx*3 and i/Nx<3:
			M2[i,j]=d

		if i/Nx==3 and (j==i-3*Nx+1 or j==i-2*Nx+1 or j==i-Nx+1) and i!=Nx*4-1:
			M2[i,j]=c

		if i==Nx*4-1 and (j==i-4*Nx+1 or j==i-3*Nx+1 or j==i-2*Nx+1):
			M2[i,j]=c

		if i/Nx==3 and i==j:
			M2[i,j]=1+b

		elif i/Nx==3 and j==i+1:
			M2[i,j]=-b

		elif i==Nx*4-1 and j==Nx*3:
			M2[i,j]=-b

M2[Nx*3-1,0]=1

def matrix_print(M2,dim1,dim2):
	for i in range(0,dim1):
		for j in range(0,dim2):
			print M2[i,j]," | ",

		print ""
		print ""


# - Linear Operator - #
numpy.matrix(M1)
numpy.matrix(M2)
M=inv(M1)*M2

#Un+1 = M * Un
U=numpy.matrix(U)
for i in range(1,Nt):
	U[:,i]=M*U[:,i-1]
	plt.plot(U[:,0])
	plt.savefig("img_implicit_deplacement/img"+str(i)+".png")

print U























