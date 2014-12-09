import numpy
import matplotlib.pyplot as plt
from matplotlib import animation



def u_initial():
    u = numpy.ones(nx)
    u[nx//2:nx] = 0
    return u

def computeF(u):
    return (u**2)/2

def maccormack(u,nt,dt,dx):
    un = numpy.zeros((nt,len(u)))
    un[:] = u.copy()
    ustar = u.copy()

    for n in range(1,nt):
        F = computeF(u)

        ustar[1:] = u[1:] - (dt/dx)*(F[1:] - F[0:-1]) + epsilon*(un[n][1:]-2*u[1:]+u[0:-1])

        Fstar = computeF(ustar)

        un[n][1:] = 0.5*(u[1:] + ustar[1:] - (dt/dx)*(Fstar[1:]-Fstar[0:-1]))

        u = un[n].copy()

    return un

###### PARAM #######33
nx = 81
nt = 70
dx = 4.0/(nx-1)
epsilon = 0.05

def animate(data):
    x = numpy.linspace(0,4,nx)
    y = data
    line.set_data(x,y)
    return line,

u = u_initial()
sigma = .5
dt = sigma*dx

un = maccormack(u,nt,dt,dx)

fig = plt.figure()
ax = plt.axes(xlim=(0,4),ylim=(-.5,2))
line, = ax.plot([],[],lw=2)

anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
plt.show()
