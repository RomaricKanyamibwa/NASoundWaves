"""
test de schemas DF et visu
"""

import numpy
import os

import scipy.sparse as sparse
import scipy.sparse.linalg

import matplotlib.pyplot as plt
import matplotlib.animation as animation

# velocity field
a_choice = 1

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

assert t_choice in [1,2]
assert x_choice in [1,2,3]
scheme_suffix = '_tchoice='+repr(t_choice)+'_xchoice='+repr(x_choice)
temp_file = 'temp_files/test_sol'+scheme_suffix
movie_filename = 'movie'+scheme_suffix+'.mp4'

# domaine spatial et maillage
X_min = 0
X_max = 1
Nx = 100
h = 1./Nx
X = numpy.zeros(Nx)

# domaine temporel
Temps_final = 1
Nt = 120
dt = 1./Nt#Temps_final * 1./Nt
dt_over_h = 1.0*dt/h

# champ de vitesse et solution approchee 
a = numpy.zeros(Nx)
u = numpy.zeros(Nx)
next_u = numpy.zeros(Nx)

if t_choice == 1:
    implicit = True
else:
    implicit = False
    
# champ vitesse
def a_func(x):
    if a_choice == 1:
        return -1.0#0.5
    elif a_choice == 2:
        return 0.25 - 0.15*float(0.25 < x < 0.5)
    else:
        print "error, wrong a_choice."
        quit()        

# fonction initiale
def u_ini(x):
    return numpy.sin(2*numpy.pi*x)

# sol exacte
def u_ex(t,x):
    assert a_choice == 1  # constant velocity
    return u_ini(x-a_func(x)*t)


for i in range(0,Nx):
    X[i] = X_min + i*h*(X_max-X_min)
    u[i] = u_ini(X[i])
    a[i] = a_func(X[i])
u_min = min(u)
u_max = max(u)
Y_min = u_min - 0.1*(u_max-u_min)
Y_max = u_max + 0.1*(u_max-u_min)

# assemblage d'une matrice creuse M
if implicit:
    row = list()
    col = list()
    data = list()
    for i in range(0,Nx):

        # M_i,i = qque chose (exemple)
        row.append((i))
        col.append((i))
        data.append( 1.3 )  
    
        # M_i,i+1 = autre chose (exemple)
        row.append((i))
        col.append((numpy.mod(i+1,Nx)))  # modulo pour conditions periodiques
        data.append( -0.3 )

    row = numpy.array(row)
    col = numpy.array(col)
    data = numpy.array(data)      
    M = (sparse.coo_matrix((data, (row, col)), shape=(Nx, Nx))).tocsr()
else:
    M = None

ims = []
fig2 = plt.figure()
Writer = animation.writers['imagemagick']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

def plot_sol(k, plot_u_ex=False):
    fig2.clf()
    print "Plot sol",  k, "min/max = ", min(u), "/", max(u)
    X_full = numpy.concatenate((X,[1.]),axis=0)
    u_full = numpy.concatenate((u,[u[0]]),axis=0)
    if plot_u_ex:
        u_ex_full = numpy.zeros(Nx+1)
        for i in range(0,Nx+1):
            u_ex_full[i] = u_ex(k*dt,X_full[i])
        print "Plot ex sol",  k, "min/max = ", min(u_ex_full), "/", max(u_ex_full)

    plt.xlim(X_min, X_max)
    plt.ylim(Y_min, Y_max)
    plt.xlabel('x')
    image = (plt.plot(X_full, u_full, '-', color='k'))
    if plot_u_ex:
        image += (plt.plot(X_full, u_ex_full, '-', color='r'))
    fig2.savefig("out"+repr(k)+".png")
    ims.append(image)

plot_sol(0,plot_u_ex=(a_choice == 1))

## schema numerique       
for nt in range(0,Nt):               
    if t_choice == 1:
        next_u[:] = sparse.linalg.spsolve(M, u)

    elif t_choice == 2:
        for i in range(0,Nx):        
            next_u[i] = -dt_over_h*(u[i]-u[i-1])+u[i]  # (exemple)
    else:
        print "error, wrong t_choice."
        quit()        
            
    u[:] = next_u[:]
    plot_sol(nt+1, plot_u_ex=(a_choice ==1) )
    
im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000, blit=True)
im_ani.save(movie_filename, writer=writer)

#os.system('rm '+temp_file+'*.pdf')  # destruction des fichiers temporaires
print "Done."
