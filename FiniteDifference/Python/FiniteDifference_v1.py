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
#  Copyright (Const_c) 2017                                                 #
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
# d_tt PH0 - 5.0/6 * d_xx PH0 = 0
#
# with P the pressure


# time discretization: 
# du/dt(t^n) =  
# 1:  ( u^n+2 - 2*u^n+1 + u^n )/dt        (impl) 
# 2:  ( u^n+1 - 2*u^n +u^n-1 )/dt        (expl) 
t_choice = 2 

# space discretization: 
# du/dx(x_i) =  
# 1:  ( u_i - 2*u_i-1 + u_i-2 )/dx 
# 2:  ( u_i+1 - 2*u_i + u_i-1 )/dx 
# 3:  ( u_i+1 - 2*u_i-1 + u_i-3 )/2dx
x_choice = 2 

implicite = True

if implicite:
    dir_name = 'simu_impliciteV1/'
else:
    dir_name = 'simu_expliciteV1/'

do_movie = False
if do_movie:    
    movie_filename = dir_name+'movie.mp4'

# spatial domain  and  meshing
X1_min = 0
X1_max = 100.0
Nx1 = 5000
h1 = (X1_max-X1_min)*1./Nx1
X1 = numpy.zeros(Nx1+1)


# temporal domain
Temps_final = 200*numpy.pi
Nt = 5000
dt = Temps_final * 1./Nt
dt_over_h1 = dt/h1

# approximate number of images generated
n_images = 1000
periode_images = int(Nt*1./n_images)

# Pressure and Numerical solution (approx sol)
PH0 = numpy.zeros(Nx1+1)
before_PH0 = numpy.zeros(Nx1+1)
next_PH0 = numpy.zeros(Nx1+1)

UH0 = numpy.zeros(Nx1+1)
next_UH0 = numpy.zeros(Nx1+1)     



Const_C=-(5./6)*(dt_over_h1*dt_over_h1)

# fonction initiale

#At t=0 and x1=xw u=0 and P=0
#plot_sol(0,0)

def velocity(t):
    return -numpy.sin(t)+(t)*numpy.exp(-t)

#d_tu1
def acceleration(t):
    return -numpy.cos(t)+(1-t)*numpy.exp(-t)

for i in range(0,Nx1+1):
    X1[i] = X1_min + i*h1#*(X1_max-X1_min)

# pour la visu:
#Y_min = 0
#Y_max = 1.1*max(v_max,u_max) 

# assembly of a sparse matrix M using the vector W of the previous time


def constr_matrix_A():
    assert implicite
    row = list()
    col = list()
    data = list()
    
    row.append((0))
    col.append((0))  
    data.append(1-Const_C)   # value of the element
    
    row.append((0))
    col.append((1))  
    data.append(Const_C)   # value of the element
    
    for i in range(1,Nx1-2):
        # M_i,i = 0
    
        # M_i,i-1
        j = i-1
        row.append((i))
        col.append((j))  
        data.append( Const_C )   # value of the element

        # M_i,i
        j = i
        row.append((i))
        col.append((j)) 
        data.append( 1-2*Const_C )  # value of the element
        
        # M_i,i
        j = i+1
        row.append((i))
        col.append((j)) 
        data.append( Const_C )  # value of the element  
    
    # M_Nx-2,Nx-3
    row.append((Nx1-2))
    col.append((Nx1-3)) 
    data.append(Const_C )  # value of the element
    
    # M_Nx-2,Nx-2
    row.append((Nx1-2))
    col.append((Nx1-2)) 
    data.append(1-Const_C )  # value of the element
    
    row = numpy.array(row)
    col = numpy.array(col)
    data = numpy.array(data)      
    M = (sparse.coo_matrix((data, (row, col)), shape=(Nx1-1, Nx1-1))).tocsr()
    return M


def constr_vect_B(Pn1,Pn2,nt):
    assert implicite
    
    B = numpy.zeros(Nx1-1)
    
    B[0]=2*Pn2[0]-Pn1[0]-2*Const_C*h1*acceleration(nt*dt)# value of the element
    
    for i in range(1,Nx1-1):
        B[i]=2*Pn2[i]-Pn1[i] # value of the element
    
    return B



#M=constr_matrix_A();
#pprint(SparseMatrix(M.todense()))
#pprint(SparseMatrix(M.todense()).inv())
#B=constr_vect_B(before_PH0,PH0,2)
#pprint(B)
#pprint(SparseMatrix(M.todense()).inv()*SparseMatrix(B))
#print("detM=",numpy.linalg.det(M.todense()))
#pprint(Wn)


#if do_movie:
    #ims = []
    #print("using imagemagick...")
    #Writer = animation.writers['imagemagick']  # ['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()


def plot_sol(n,ielem):
    fig.clf()
    fname = dir_name+"out"+repr(n)+"_PH"+str(ielem)+".png"
    #fname2 = dir_name+"out"+repr(n)+"_PH"+str(ielem)+".dat"
    #f= open(fname2,"w+")
    #for i in range(Nx1+1):
        #str1=str(X1[i])+"\t"+str(PH0[i])+"\t"+str(UH0[i])+"\n"
        #f.write(str1)
    print("Plot sol in file ", fname, ", nt = ", n, ", min/max = ", min(PH0), "/", max(PH0))
    X_full = numpy.concatenate((X1,[X1[-1]]),axis=0)
    U_full = numpy.concatenate((UH0,[UH0[0]]),axis=0)
    P_full = numpy.concatenate((PH0,[PH0[0]]),axis=0)
    # trace aussi les vitesses:
    #v = map(vitesse, u)
    #v_full = numpy.concatenate((v,[v[0]]),axis=0)
    #plt.xlim(X1_min, X1_max)
    plt.ylim(-8, 8)
    #plt.xlabel('x')
    image = (plt.plot(X1, PH0, '-', color='r'),plt.plot(X_full, U_full,'-', color='k'))
    plt.legend(['Pressure', 'Velocity'], loc='upper left')
    #plt.legend(['PH0'], loc='upper left')
    #maximage += (plt.plot(X_full, v_full, '-', color='b'))
    
    fig.savefig(fname)
    if do_movie:
        ims.append(image)

plot_sol(0,0)
plot_sol(1,0)

M = constr_matrix_A()
#print("-------------------------------------------------------")
###pprint(SparseMatrix(Wn))
#pprint(M.todense())
#print("-------------------------------------------------------")

# schema numerique       
for nt in range(1,Nt):               
    if implicite:
        #print("Implicit Scheme")
        B = constr_vect_B(before_PH0[1:Nx1],PH0[1:Nx1],nt+1)
        #pprint(SparseMatrix(B))
        next_PH0[1:Nx1] = sparse.linalg.spsolve(M, B)
        next_PH0[0]= next_PH0[1] + 2*h1*acceleration((nt+1)*dt)
        next_PH0[Nx1]=next_PH0[Nx1-1]

    else:
        print("ERROR")
        exit()
    
    #pressure
    before_PH0[:]=PH0[:]
    PH0[:] = next_PH0[:]
    
    #velocity
    next_UH0[0]=velocity((nt+1)*dt)
    next_UH0[Nx1]=0.0
    next_UH0[1:Nx1]=UH0[1:Nx1]-1.0/2.0*dt_over_h1*(next_PH0[1:Nx1]-next_PH0[0:Nx1-1])
    
    UH0[:]=next_UH0[:]

    if (nt+1)%periode_images == 0 or (nt+1) == Nt:
        #print("Nt=",nt+1)
        #pprint(PH0)
        plot_sol(nt+1,0)
    

#if do_movie:
    #im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
    #im_ani.save(movie_filename, writer=writer)

#pprint(SparseMatrix(U).T)
#M=assemble_M()
##pprint(Wn)
#Minv=numpy.linalg.inv(M.todense())
#print("Inverse det=",numpy.linalg.det(Minv))
#print("detM=",numpy.linalg.det(M.todense()))
print("Done.")
