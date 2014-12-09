### MODULE 3 ###########
### ASSESSMENT #########

### Sod's Shock Tube ###

import numpy as np

nx = 81
dx = .25
dt = .0002   
gamma = 1.4
T = 0.01
nt = 51

# Pressure given by equation of state
def p(u):
    return (gamma - 1)*(u[2] - 0.5*(u[1]**2/u[0]))

# Flux terms vector given by u
def computef(u):
    f = np.zeros((nx,3))
    for i in range(0,nx):
        f[i] = np.array([u[i][1],
                         u[i][1]**2/u[i][0] + p(u[i]),
                         (u[i][2] + p(u[i]))*(u[i][1]/u[i][0])])
    return f
        
    

# Initial Conditions
def vecu_initial(rho,u,p):
    return np.array([rho,
                     rho*u,
                     rho*(0.5*u**2 + p/((gamma-1)*rho))])


# Richtmyer method
def richtmyer(u, nt, dt, dx):
    un = np.zeros((nt,nx,3))
    un[:,:,:] = u.copy()
    umid = u.copy()

    for n in range(1,nt):
        f = computef(u)
        umid[:-1] =  0.5*(u[1:] + u[:-1]) - dt/(2*dx) *(f[1:] - f[:-1])
        fmid = computef(umid)
        un[n,1:] =  u[1:] - (dt/dx)*(fmid[1:] - fmid[:-1])
        u = un[n].copy()

    return un

# Setting initial state
def u_initial():
    u = np.ones((nx,3))
    u[0:nx//2+1] = vecu_initial(1,0,100000)
    u[nx//2+1:nx] = vecu_initial(0.125,0,10000)
    return u

u = u_initial()
#print(u[39:43])
#print("=====================")

un = richtmyer(u, nt, dt, dx)

#print(un[len(un)-1])
for t in [1,2,3]:
    print("=====================")
    print("=====================")
    print("@ t = 50-"+str(t)+"\n=====================")
    v = un[len(un)-t]
    for x in [49,50,51]:
        print("@ x="+str(x)+"\n=====================")
        print("Velocity")
        print(str(v[x][1]/v[x][0]))
        print("=====================")
        print("Pressure")
        print(str(p(v[x])))
        print("=====================")
        print("Density")
        print(str(v[x][0]))
        print("=====================")

        
