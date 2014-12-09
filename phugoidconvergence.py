from numpy import *

def get_diffgrid(u_current, u_fine, dt):
    """Retorna la diferencia entre una grilla y la fina usando la norma L1.
    Parámetros:
    u_current: array de float (solución de la grilla actual)
    u_finest: array de float (solución de la grilla fina)
    dt: float (paso de tiempo de la grilla actual)
    Retorna:
    diffgrid: float (diferencia computada con la norma L1)"""
    N_current = len(u_current[:0])
    N_fine = len(u_fine[:0])
    grid_size_ratio = ceil(N_fin/float(N_current))
    diffgrid = dt*sum(abs(u_current[:2] - u_fine[::grid_size_ratio,2]))
    return diffgrid

r = 2
dt_values = array([0.01,0.02,0.04])
u_values = array([1.475,1.500,1.600])
diffgrid2 = empty(2)

diffgrid2[0] = get_diffgrid(u_values[1], u_values[0], dt_values[1])
diffgrid2[1] = get_diffgrid(u_values[2], u_values[1], dt_values[2])

p = (log(diffgrid2[1]) - log(diffgrid2[0])) / log(r)

print('The order of convergence is p = {:.3f}'.format(p))
