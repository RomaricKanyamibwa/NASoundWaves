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
# d_t u1 + 1/2 * d_x1 PH0 = 0
# d_t PH0 + 5/3 * dj_x1 u1 = 0
#
# dx2=0 and dx3=0
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
x_choice = 1 

Const_c=1./2
Const_d=5./3

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
X1_max = 1.0
Nx1 = 100
h1 = 1./Nx1
X1 = numpy.zeros(Nx1)


# temporal domain
Temps_final = 1.#/60
Nt = 500
dt = Temps_final * 1./Nt
dt_over_h1 = dt/h1

# approximate number of images generated
n_images = 100
periode_images = int(Nt*1./n_images)

# Pressure and Numerical solution (approx sol)
PH0 = numpy.zeros(Nx1)

u1 = numpy.zeros(Nx1)
next_u1 = numpy.zeros(Nx1)     


# fonction initiale

#At t=0 and x1=xw u=0 and P=0
#plot_sol(0,0)

def u1_ini(t):
    return -numpy.sin(t)+(t)*numpy.exp(-t)

#d_tu1
def acceleration(t):
    return -numpy.cos(t)+(1-t)*numpy.exp(-t)

for i in range(0,Nx1):
    X1[i] = X1_min + i*h1*(X1_max-X1_min)

# pour la visu:
#Y_min = 0
#Y_max = 1.1*max(v_max,u_max) 

# assembly of a sparse matrix M using the vector W of the previous time


#def constr_matrix_Ai(ielem,const):
    #assert implicite
    #row = list()
    #col = list()
    #data = list()
    #coef=const
    
    #row.append((0))
    #col.append((0))  
    #data.append(coef*dt_over_h1)   # value of the element
    
    #for i in range(1,Nx1):
        ## M_i,i = 0
    
        ## M_i,i-1
        #j = i-1#numpy.mod(i-1,Nx_array[ielem])
        #row.append((i))
        #col.append((j))  
        #data.append( -1*coef*dt_over_h1 )   # value of the element

        ## M_i,i
        #j = i#numpy.mod(i+1,Nx_array[ielem])
        #row.append((i))
        #col.append((j)) 
        #data.append( coef*dt_over_h1 )  # value of the element  
         
    #row = numpy.array(row)
    #col = numpy.array(col)
    #data = numpy.array(data)      
    #M = (sparse.coo_matrix((data, (row, col)), shape=(Nx1, Nx1))).tocsr()
    #return M

#M=constr_matrix_Ai(0,Const_c);
##pprint(SparseMatrix(M.todense()))
#print("detM=",numpy.linalg.det(M.todense()))
##pprint(Wn)

#def assemble_M():
    #assert implicite
    #A1=constr_matrix_Ai(0,Const_c)
    #A1_d=constr_matrix_Ai(0,Const_d)
    
    #I=sparse.identity(Nx1);
    #I=sparse.coo_matrix(I)

    
    ##pprint(SparseMatrix(A1.todense()))
    ##pprint(SparseMatrix(B1.todense()))
    
    #M1=numpy.hstack((I.todense(),A1.todense()))
    
    #M2=numpy.hstack((A1_d.todense(),I.todense()))
    
    #M=numpy.vstack((M1,M2))
    #M=sparse.coo_matrix(M)
    ##pprint(M)
    ##print(numpy.linalg.det(M.todense()))
    ##print(numpy.linalg.det(MB.todense()))
    ##pprint(SparseMatrix(M.todense()))
    ##pprint(SparseMatrix(MB.todense()))
    #return M

##if do_movie:
    ##ims = []
    ##print("using imagemagick...")
    ##Writer = animation.writers['imagemagick']  # ['ffmpeg']
    ##writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

#fig = plt.figure()


#def plot_sol(n,ielem):
    #fig.clf()
    #fname = dir_name+"out_"+repr(n)+"_u"+str(ielem)+".png"
    #print("Plot sol in file ", fname, ", nt = ", n, ", min/max = ", min(u1), "/", max(u1))
    #X_full = numpy.concatenate((X1,[1.]),axis=0)
    #u_full = numpy.concatenate((u1,[u1[0]]),axis=0)
    #P_full = numpy.concatenate((PH0,[PH0[0]]),axis=0)
    ## trace aussi les vitesses:
    ##v = map(vitesse, u)
    ##v_full = numpy.concatenate((v,[v[0]]),axis=0)
    #plt.xlim(X1_min, X1_max)
    ##plt.ylim(Y_min, Y_max)
    #plt.xlabel('x')
    #image = (plt.plot(X_full, u_full, '-', color='k'),plt.plot(X_full, P_full, '-', color='r'))
    #plt.legend(['u1', 'PH0'], loc='upper left')
    ##maximage += (plt.plot(X_full, v_full, '-', color='b'))
    
    #fig.savefig(fname)
    #if do_movie:
        #ims.append(image)

#plot_sol(0,0)
#U=u1
#P=PH0

#u1[0] = u1_ini(dt)
#u1[1] = -Const_c*dt_over_h1*(-2*h1*acceleration(dt))
#PH0[1]= -Const_d*dt_over_h1*(u1[1]-u1[0])
##print(u1[0],u1[1],PH0[1])
#PH0[0]= PH0[1] + 2*h1*acceleration(dt)
##print("PH0=",PH0[0])
#um=u1[1]
#pm0=PH0[1]

#for i in range(2,Nx1):
    #a = numpy.array([[1,Const_c*dt_over_h1], [Const_d*dt_over_h1,1]])
    #b = numpy.array([0+Const_c*dt_over_h1*PH0[i-1],Const_d*dt_over_h1*u1[i-1]])
    #x = numpy.linalg.solve(a, b)
    
    #u1[i] = x[0]
    #PH0[i]= x[1]
    ##print("______________________________________________")
    ##print(i,"linalg:",u1[i])
    ##print(i,"linalg:",PH0[i])
    ##pm=Const_d*dt_over_h1*(um-Const_c*dt_over_h1*pm0)/(1-Const_d*dt_over_h1*Const_c*dt_over_h1)
    ##um=Const_c*dt_over_h1*(pm0-pm)
    ##pm0=pm
    ##print(i,"form:",um)
    ##print(i,"form:",pm)
    ##print("X=",x)

#u1[Nx1-1] = 0    
#PH0[Nx1-1] = PH0[Nx1-2]
##print("u1=",u1)
##print("PH0=",PH0)
#plot_sol(1,0)
#U=numpy.vstack((U,u1))
#P=numpy.vstack((P,PH0))
#Wn=numpy.hstack((u1,PH0))
#next_Wn = numpy.zeros(2*Nx1)
#old_Wn = numpy.zeros(2*Nx1)
#M = assemble_M()
##print("-------------------------------------------------------")
##pprint(SparseMatrix(Wn))
##pprint(M.todense())
##print("-------------------------------------------------------")

## schema numerique       
#for nt in range(1,Nt):               
    #if implicite:
        ##print("Implicit Scheme")
        #next_Wn[:] = sparse.linalg.spsolve(M, Wn)

    #else:
        #print("ERROR")
        #exit()
        ##for i in range(0,Nx1):        
            ##next_u[i] =-dt_over_h*(flux(u[i])-flux(u[i-1]))+u[i]   # ecrire ici le schema explicite
            
    #Nx=Nx1
    #old_Wn[:]=Wn
    #Wn[:] = next_Wn[:]
    #Wn[0]= u1_ini((nt+1)*dt)
    #Wn[Nx-1]=0
    #u1[:] = Wn[0:Nx]
    #PH0[:] = Wn[Nx:2*Nx]
    #PH0[0]= PH0[1] + 2*h1*acceleration((nt+1)*dt)
    #PH0[Nx-1]=PH0[Nx-2]
    ##pprint(Wn)
    ##U=numpy.vstack((U,u1))
    ##P=numpy.vstack((P,PH0))

    #if (nt+1)%periode_images == 0 or (nt+1) == Nt:
        #plot_sol(nt+1,0)
    

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
