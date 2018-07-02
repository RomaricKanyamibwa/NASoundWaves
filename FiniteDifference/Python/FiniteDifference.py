"""
test  schemas DF et visu
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 15:25:35 2018
@author: SATGE Lancelot, KANYAMIBWA Romaric
Resolution of PDE using Finite Difference
================
Test of Finite Difference Numerical Schemes on sound waves equations
"""

#############################################################################
#  Copyright (Const_c) 2017                                                       #
#                                                                           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################
import numpy
import os

import scipy.sparse as sparse
import scipy.sparse.linalg

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import random

from pprint import pprint

from sympy.matrices import SparseMatrix

### We want to resolve numericaly the following PDE
#
# d_t ui + 1/2 * d_xi PH0 = 0
# d_t PH0 + 5/3 * dj_xj uj = 0
#
# with u the velocity and P the pressure

## vitesse max km/h
#v_max = 130

## densite max (nb de voitures/km)
#u_max = 200


#def vitesse(u_val):    
    #if(0<u_val<u_max):
        #return -0.65*u_val+130
    #else:
        #return 0

#def flux(u_val):    
    #return u_val*vitesse(u_val)



# time discretization: 
# du/dt(t^n) =  
# 1:  ( u^n - u^n-1 )/dt        (impl) 
# 2:  ( u^n+1 - u^n )/dt        (expl) 
t_choice = 2 

# space discretization: 
# du/dx(x_i) =  
# 1:  ( u_i - u_i-1 )/dx 
# 2:  ( u_i+1 - u_i )/dx 
# 3:  ( u_i+1 - u_i-1 )/2dx
x_choice = 3

Const_c=-1./2
Const_d=-5./3

implicite = True

if implicite:
    dir_name = 'simu_implicite/'
else:
    dir_name = 'simu_explicite/'

do_movie = False
if do_movie:    
    movie_filename = dir_name+'movie.mp4'

# spatial domain  and  meshing
X1_min = 0
X1_max = 1
Nx1 = 100
h1 = 1./Nx1
X1 = numpy.zeros(Nx1)

X2_min = 0
X2_max = 1
Nx2 = Nx1
h2 = 1./Nx2
X2 = numpy.zeros(Nx2)

X3_min = 0
X3_max = 1
Nx3 = Nx1
h3 = 1./Nx3
X3 = numpy.zeros(Nx3)

# temporal domain
Temps_final = 1.#/60
Nt = 500
dt = Temps_final * 1./Nt
dt_over_h1 = dt/h1
dt_over_h2 = dt/h2
dt_over_h3 = dt/h3

dt_array=[dt_over_h1,dt_over_h2,dt_over_h3]
h_array=[h1,h2,h3]
Nx_array=[Nx1,Nx2,Nx3]
X_array=[[X1_min,X1_max,X1],[X2_min,X2_max,X2],[X3_min,X3_max,X3]]

# approximate number of images generated
n_images = 100
periode_images = int(Nt*1./n_images)

# Pressure and Numerical solution (approx sol)
PH0 = numpy.zeros(max(Nx1,Nx2,Nx3))

u1 = numpy.zeros(Nx1)
next_u1 = numpy.zeros(Nx1)     

u2 = numpy.zeros(Nx2)
next_u2 = numpy.zeros(Nx2)     

u3 = numpy.zeros(Nx3)
next_u3 = numpy.zeros(Nx3)

u_array=[u1,u2,u3,PH0]


# fonction initiale

def u1_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

def u2_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

def u3_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

for i in range(0,Nx1):
    X1[i] = X1_min + i*h1*(X1_max-X1_min)
    u1[i] = random.random()#u1_ini(X1[i])

for i in range(0,Nx2):
    X2[i] = X2_min + i*h2*(X2_max-X2_min)
    u2[i] = random.random()#u2_ini(X2[i])
    
for i in range(0,Nx3):
    X3[i] = X3_min + i*h3*(X3_max-X3_min)
    u3[i] = random.random()#u3_ini(X3[i])
# pour la visu:
#Y_min = 0
#Y_max = 1.1*max(v_max,u_max) 

# assembly of a sparse matrix M using the vector W of the previous time

Wn=numpy.hstack((u1,u2))
Wn=numpy.hstack((Wn,u3))
Wn=numpy.hstack((Wn,PH0))
next_Wn = numpy.zeros(4*Nx1)

def constr_matrix_Ai(ielem,const,sign=1):
    assert implicite
    row = list()
    col = list()
    data = list()
    for i in range(0,Nx_array[ielem]):
        
        coef=sign*const
        # M_i,i = 0
    
        # M_i,i-1
        j = numpy.mod(i-1,Nx_array[ielem]) # mod for periodic conditions
        row.append((i))
        col.append((j))  
        data.append( coef*1./2*dt_array[ielem] )   # value of the element

        # M_i,i+1
        j = numpy.mod(i+1,Nx_array[ielem])  # mod for periodic conditions
        row.append((i))
        col.append((j)) 
        data.append( -coef*1./2*dt_array[ielem] )  # value of the element  
        
    row = numpy.array(row)
    col = numpy.array(col)
    data = numpy.array(data)      
    M = (sparse.coo_matrix((data, (row, col)), shape=(Nx_array[ielem], Nx_array[ielem]))).tocsr()
    return M

def assemble_M():
    assert implicite
    A1=constr_matrix_Ai(0,Const_c)
    A2=constr_matrix_Ai(1,Const_c)
    A3=constr_matrix_Ai(2,Const_c)
    
    A1_d=constr_matrix_Ai(0,Const_d)
    A2_d=constr_matrix_Ai(1,Const_d)
    A3_d=constr_matrix_Ai(2,Const_d)
    
    I=sparse.identity(Nx1);
    I=sparse.coo_matrix(I)
    Zero=sparse.coo_matrix([[0 for i in range(Nx1)] for j in range(Nx1)])
    
    #pprint(Zero)
    #pprint(Zero.todense())
    #pprint(sparse.coo_matrix(I))
    #pprint((I))
    B1=constr_matrix_Ai(0,Const_c,-1)
    B2=constr_matrix_Ai(1,Const_c,-1)
    B3=constr_matrix_Ai(2,Const_c,-1)
    
    B1_d=constr_matrix_Ai(0,Const_d,-1)
    B2_d=constr_matrix_Ai(1,Const_d,-1)
    B3_d=constr_matrix_Ai(2,Const_d,-1)
    
    #pprint(SparseMatrix(A1.todense()))
    #pprint(SparseMatrix(B1.todense()))
    
    M1=numpy.hstack((I.todense(),Zero.todense()))
    M1=numpy.hstack((M1,Zero.todense()))
    M1=numpy.hstack((M1,A1_d.todense()))
    
    M2=numpy.hstack((Zero.todense(),I.todense()))
    M2=numpy.hstack((M2,Zero.todense()))
    M2=numpy.hstack((M2,A2_d.todense()))
    
    M3=numpy.hstack((Zero.todense(),Zero.todense()))
    M3=numpy.hstack((M3,I.todense()))
    M3=numpy.hstack((M3,A3_d.todense()))
    
    M4=numpy.hstack((A1_d.todense(),A2_d.todense()))
    M4=numpy.hstack((M4,A3_d.todense()))
    M4=numpy.hstack((M4,I.todense()))
    
    M=numpy.vstack((M1,M2))
    M=numpy.vstack((M,M3))
    M=numpy.vstack((M,M4))
    M=sparse.coo_matrix(M)
    pprint(M)
    
    M1=numpy.hstack((I.todense(),Zero.todense()))
    M1=numpy.hstack((M1,Zero.todense()))
    M1=numpy.hstack((M1,B1_d.todense()))
    
    M2=numpy.hstack((Zero.todense(),I.todense()))
    M2=numpy.hstack((M2,Zero.todense()))
    M2=numpy.hstack((M2,B2_d.todense()))
    
    M3=numpy.hstack((Zero.todense(),Zero.todense()))
    M3=numpy.hstack((M3,I.todense()))
    M3=numpy.hstack((M3,B3_d.todense()))
    
    M4=numpy.hstack((B1_d.todense(),B2_d.todense()))
    M4=numpy.hstack((M4,B3_d.todense()))
    M4=numpy.hstack((M4,I.todense()))
    
    MB=numpy.vstack((M1,M2))
    MB=numpy.vstack((MB,M3))
    MB=numpy.vstack((MB,M4))
    MB=sparse.coo_matrix(MB)
    pprint(MB)
    #print(numpy.linalg.det(M.todense()))
    #print(numpy.linalg.det(MB.todense()))
    #pprint(SparseMatrix(M.todense()))
    #pprint(SparseMatrix(MB.todense()))
    return M,MB

if do_movie:
    ims = []
    print("using imagemagick...")
    Writer = animation.writers['imagemagick']  # ['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()


def plot_sol(n,ielem):
    
    fig.clf()
    fname = dir_name+"out_"+repr(n)+"_u"+str(ielem)+".png"
    print("Plot sol in file ", fname, ", nt = ", n, ", min/max = ", min(u_array[ielem]), "/", max(u_array[ielem]))
    X_full = numpy.concatenate((X_array[ielem][2],[1.]),axis=0)
    u_full = numpy.concatenate((u_array[ielem],[u_array[ielem][0]]),axis=0)
    # trace aussi les vitesses:
    #v = map(vitesse, u)
    #v_full = numpy.concatenate((v,[v[0]]),axis=0)
    plt.xlim(X_array[ielem][0], X_array[ielem][1])
    #plt.ylim(Y_min, Y_max)
    plt.xlabel('x')
    image = (plt.plot(X_full, u_full, '-', color='k'))
    #maximage += (plt.plot(X_full, v_full, '-', color='b'))
    
    fig.savefig(fname)
    if do_movie:
        ims.append(image)

plot_sol(0,0)
print(Wn)
# schema numerique       
for nt in range(0,Nt):               
    if implicite:
        print("Implicit Scheme")
        M,MB = assemble_M()
        next_Wn[:] = sparse.linalg.spsolve(M, MB*Wn)

    else:
        print("ERROR")
        exit()
        #for i in range(0,Nx1):        
            #next_u[i] =-dt_over_h*(flux(u[i])-flux(u[i-1]))+u[i]   # ecrire ici le schema explicite
            
    Nx=Nx1
    Wn[:] = next_Wn[:]
    u1[:] = next_Wn[0:Nx]
    u2[:] = next_Wn[Nx:2*Nx]
    u3[:] = next_Wn[2*Nx:3*Nx]
    PH0[:] = next_Wn[3*Nx:]
    u_array=[u1,u2,u3,PH0]
    #if (nt+1)%periode_images == 0 or (nt+1) == Nt:
    plot_sol(nt+1,0)
    

if do_movie:
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
    im_ani.save(movie_filename, writer=writer)

M,MB=assemble_M()
#pprint(MB*Wn)
print("Done.")
