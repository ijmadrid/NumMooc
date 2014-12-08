# Rocket in purely vertical flight

# dh/dt = v
# (ms + mp) dv/dt  =  -(ms + mp)g + m'p*ve - (1/2)*rho*v*abs(v)*A*Cd

from numpy import *
import matplotlib.pyplot as plt

ms = 50
g = 9.81
rho =  1.091
r = 0.5
A = pi*(r**2)
ve = 325
CD = 0.15
mp_inicial = 100

def mp(t):
    if t < 5:
        return mp_inicial - 20*t
    else:
        return 0

def dmp(t):
    if t < 5:
        return 20
    else:
        return 0

dt = 0.1
T = 100
N = T/dt
h = zeros(N)
v = zeros(N)
h[0] = 0
v[0] = 0
for i in range(int(N-1)):
    v[i+1] = v[i] + dt*(-g + (dmp((i)*dt)*ve)/(ms+mp((i)*dt)) - ((0.5)*rho*v[i]*abs(v[i])*A*CD)/(ms+mp((i)*dt)))
    h[i+1] = h[i] + v[i]*dt
    if h[i] < 0:
        print(str(i)+'.....SUELO!!!!')
        break

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
t = linspace(0,T,num=N)
ax1.plot(t,h)
fig1.show()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
t = linspace(0,T,num=N)
ax2.plot(t,v)
fig2.show()




