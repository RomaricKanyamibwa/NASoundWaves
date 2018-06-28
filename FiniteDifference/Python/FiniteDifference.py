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
#  Copyright (C) 2017                                                       #
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

from pprint import pprint

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
Nx2 = 100
h2 = 1./Nx2
X2 = numpy.zeros(Nx2)

X3_min = 0
X3_max = 1
Nx3 = 100
h3 = 1./Nx3
X3 = numpy.zeros(Nx3)

# temporal domain
Temps_final = 1./60
Nt = 500
dt = Temps_final * 1./Nt
dt_over_h = dt/h

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

# fonction initiale

def u1_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

def u2_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

def u3_ini(x):
    return 60+20*numpy.sin(2*numpy.pi*x)

for i in range(0,Nx1):
    X1[i] = X1_min + i*h1*(X1_max-X1_min)
    u1[i] = u1_ini(X1[i])

for i in range(0,Nx2):
    X2[i] = X2_min + i*h2*(X2_max-X2_min)
    u2[i] = u2_ini(X2[i])
    
for i in range(0,Nx3):
    X3[i] = X3_min + i*h3*(X3_max-X3_min)
    u3[i] = u3_ini(X3[i])
# pour la visu:
Y_min = 0
Y_max = 1.1*max(v_max,u_max) 

# assembly of a sparse matrix M using the vector W of the previous time

def assemble_M(some_u):
    assert implicite
    row = list()
    col = list()
    data = list()
    for i in range(0,Nx):

        # M_i,i
        row.append((i))
        col.append((i))
        data.append( 1 )   # mettre la bonne valeur
    
        # M_i,i-1
        j = numpy.mod(i-1,Nx) # modulo pour conditions periodiques
        row.append((i))
        col.append((j))  
        data.append( 0 )   # mettre la bonne valeur

        # M_i,i+1
        j = numpy.mod(i+1,Nx)  # modulo pour conditions periodiques
        row.append((i))
        col.append((j)) 
        data.append( 0 )   # mettre la bonne valeur

    row = numpy.array(row)
    col = numpy.array(col)
    data = numpy.array(data)      
    M = (sparse.coo_matrix((data, (row, col)), shape=(Nx, Nx))).tocsr()
    return M

if do_movie:
    ims = []
    print("using imagemagick...")
    Writer = animation.writers['imagemagick']  # ['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()


def plot_sol(n):
    
    fig.clf()
    fname = dir_name+"out_"+repr(n)+".png"
    print "Plot sol in file ", fname, ", nt = ", n, ", min/max = ", min(u), "/", max(u)
    X_full = numpy.concatenate((X,[1.]),axis=0)
    u_full = numpy.concatenate((u,[u[0]]),axis=0)
    # trace aussi les vitesses:
    v = map(vitesse, u)
    v_full = numpy.concatenate((v,[v[0]]),axis=0)
    plt.xlim(X_min, X_max)
    plt.ylim(Y_min, Y_max)
    plt.xlabel('x')
    image = (plt.plot(X_full, u_full, '-', color='k'))
    image += (plt.plot(X_full, v_full, '-', color='b'))
    
    fig.savefig(fname)
    if do_movie:
        ims.append(image)

plot_sol(0)

## schema numerique       
for nt in range(0,Nt):               
    if implicite:
        M = assemble_M(u)
        next_u[:] = sparse.linalg.spsolve(M, u)

    else:
        for i in range(0,Nx):        
            next_u[i] =-dt_over_h*(flux(u[i])-flux(u[i-1]))+u[i]   # ecrire ici le schema explicite
            
    u[:] = next_u[:]
    if (nt+1)%periode_images == 0 or (nt+1) == Nt:
        plot_sol(nt+1)
    

if do_movie:
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
    im_ani.save(movie_filename, writer=writer)

print "Done."
