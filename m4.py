import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
from matplotlib import animation

def ftcs(u, v, nt, tipo, dt, dh):

    if tipo == 'Bacteria 1':
        Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1
    elif tipo == 'Bacteria 2':
        Du, Dv, F, k = 0.00014, 0.00006, 0.035, 0.065 # Bacteria 2
    elif tipo == 'Coral':
        Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 # Coral
    elif tipo == 'Fingerprint':
        Du, Dv, F, k = 0.00019, 0.00005, 0.060, 0.062 # Fingerprint
    elif tipo == 'Spirals':
        Du, Dv, F, k = 0.00010, 0.00010, 0.018, 0.050 # Spirals
    elif tipo == 'Spirals Dense':
        Du, Dv, F, k = 0.00012, 0.00008, 0.020, 0.050 # Spirals Dense
    elif tipo == 'Spirals Fast':
        Du, Dv, F, k = 0.00010, 0.00016, 0.020, 0.050 # Spirals Fast
    elif tipo == 'Unstable':
        Du, Dv, F, k = 0.00016, 0.00008, 0.020, 0.055 # Unstable
    elif tipo == 'Worms 1':
        Du, Dv, F, k = 0.00016, 0.00008, 0.050, 0.065 # Worms 1
    elif tipo == 'Worms 2':
        Du, Dv, F, k = 0.00016, 0.00008, 0.054, 0.063 # Worms 2
    elif tipo == 'Zebrafish':
        Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.060 # Zebrafish

    dx = dh
    dy = dh
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        u[1:-1,1:-1] = un[1:-1,1:-1] + Du *\
            (dt/dy**2 * (un[2:,1:-1] - 2*un[1:-1,1:-1] + un[:-2,1:-1]) +\
             dt/dx**2 * (un[1:-1,2:] - 2*un[1:-1,1:-1] + un[1:-1,:-2])) -\
             dt*un[1:-1,1:-1]*(vn[1:-1,1:-1])**2 + dt*F*(1-un[1:-1,1:-1])

        v[1:-1,1:-1] = vn[1:-1,1:-1] + Dv *\
            (dt/dy**2 * (vn[2:,1:-1] - 2*vn[1:-1,1:-1] + vn[:-2,1:-1]) +\
             dt/dx**2 * (vn[1:-1,2:] - 2*vn[1:-1,1:-1] + vn[1:-1,:-2])) +\
             dt*un[1:-1,1:-1]*(vn[1:-1,1:-1])**2 - dt*(F+k)*vn[1:-1,1:-1]

  
        # Enforce Neumann BCs
        u[-1,:] = u[-2,:]
        u[:,-1] = u[:,-2]
        u[0,:] = u[1,:]
        u[:,0] = u[:,1]
        
        v[-1,:] = v[-2,:]
        v[:,-1] = v[:,-2]
        v[0,:] = v[1,:]
        v[:,0] = v[:,1]
        

    print(u[100,::40])
    return u

n = 192

Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1

#Du, Dv, F, k = 0.00014, 0.00006, 0.035, 0.065 # Bacteria 2
#Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 # Coral
#Du, Dv, F, k = 0.00019, 0.00005, 0.060, 0.062 # Fingerprint
#Du, Dv, F, k = 0.00010, 0.00010, 0.018, 0.050 # Spirals
#Du, Dv, F, k = 0.00012, 0.00008, 0.020, 0.050 # Spirals Dense
#Du, Dv, F, k = 0.00010, 0.00016, 0.020, 0.050 # Spirals Fast
#Du, Dv, F, k = 0.00016, 0.00008, 0.020, 0.055 # Unstable
#Du, Dv, F, k = 0.00016, 0.00008, 0.050, 0.065 # Worms 1
#Du, Dv, F, k = 0.00016, 0.00008, 0.054, 0.063 # Worms 2
#Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.060 # Zebrafish

dh = 5./(n-1)
T = 8000
dt = .9 * dh**2 / (4*max(Du,Dv))
nt = int(T/dt)

x = numpy.linspace(0,5.,dh)
y = numpy.linspace(0,5.,dh)

uvinitial = numpy.load('./uvinitial.npz')
U = uvinitial['U']
V = uvinitial['V']

U = U.copy()
V = V.copy()
#Uf = ftcs(U, V, nt, Du, Dv, F, k, dt, dh)

l = ['Bacteria 1', 'Bacteria 2', 'Coral', 'Fingerprint', 'Spirals', 'Spirals Dense',
     'Spirals Fast', 'Unstable', 'Worms 1', 'Worms 2', 'Zebrafish']

for i in l:    
    fig = plt.figure(figsize=(8,5))
    ax = plt.subplot(111)
    cax = ax.imshow(ftcs(U, V, nt, i, dt, dh), cmap = cm.RdBu)
    cbar = fig.colorbar(cax)
    plt.savefig('Figures/'+i+'.png')
    plt.close(fig)
