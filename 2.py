### MOOC
### SPACE & TIME: AN INTRODUCTION TO FINITE DIFFERENCE SOLUTION OF PDE'S
### GRADED HOMEWORK
### TRAFFIC FLOW EQUATION

from numpy import *

def flux(Vmax,T):
    ### Parameters
    rhomax  = 250
    L       = 11
    nx      = 51
    dt      = 0.001
    dx      = L/(nx-1)
    nt      = int(T/dt) + 1
     
    ### F := flux (cars/hr)
    def Fx(rho):
        return Vmax*rho*(1-(rho/rhomax))

    def Vx(rho):
        return Vmax*(1-(rho/rhomax))

    def dF(rho):
        return -Vmax*(rho/rhomax) + Vx(rho)

    ### PDE: ##############################
    ###
    ###          dp        dF
    ###         ----   +  ----  =  0
    ###          dt        dx
    ###
    #######################################

    x = linspace(0,L,nx)
    ### Condiciones Iniciales
    rho = ones(nx)*20
    rho[10:20] = 50

    # Vector de velocidades
    V = Vx(rho)

    for n in range(1,nt):
        rhon = rho.copy()
        Fn = Fx(rhon)
        rho[1:] = rhon[1:] - (dt/dx)*(Fn[1:]-Fn[0:-1])
        rho[0] = 10.0
        V[0:] = Vx(rho[0:])

    return V
    



