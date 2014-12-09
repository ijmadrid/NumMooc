### MODULE 3 ###########
### ASSESSMENT #########

### Sod's Shock Tube ###

import numpy as np

nx = 81
dx = .25
dt = .0002   
gamma = 1.4
T = 0.01
nt = 50

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
        for i in range(2,nx-1):
            f = computef(u)
            umid[i,:] =  0.5*(u[i+1,:] + u[i,:]) - dt/(2*dx) *(f[i+1,:] - f[i,:])
            fmid = computef(umid)
            un[n,i,:] =  u[i,:] - (dt/dx)*(fmid[i,:] - fmid[i-1,:])
        u = un[n].copy()

    return un

# Setting initial state
def u_initial():
    u = np.ones((nx,3))
    u[0:nx//2] = vecu_initial(1,0,100000)
    u[nx//2:nx] = vecu_initial(0.125,0,10000)
    return u

u = u_initial()

un = richtmyer(u, nt, dt, dx)
print(un[len(un)-1])
        
