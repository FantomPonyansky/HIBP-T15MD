# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:47:47 2023

@author: Krohalev_OD
"""
#%% imports
import os
from dataclasses import dataclass
import numpy as np
import numba

#%% auxilliary types

@dataclass
class SecBeamlineData:
    xaim: float  
    yaim: float
    zaim: float
    # alpha and beta angles of the SECONDARY beamline [deg]
    alpha_sec: float
    beta_sec: float
    gamma_sec: float

    def __iter__(self):   
        yield self.xaim
        yield self.yaim
        yield self.zaim

        yield self.alpha_sec
        yield self.beta_sec
        yield self.gamma_sec


# %% Runge-Kutta

# define equations of movement:

    

@numba.njit()
def f(k, E, V, B):
    return k*(E + np.cross(V, B))


@numba.njit()
def g(V):
    return V


@numba.njit()
def runge_kutt(k, RV, dt, E, B):
    '''
    Calculate one step using Runge-Kutta algorithm

    V' = k(E + [VxB]) == K(E + np.cross(V,B)) == f
    r' = V == g

    V[n+1] = V[n] + (h/6)(m1 + 2m2 + 2m3 + m4)
    r[n+1] = r[n] + (h/6)(k1 + 2k2 + 2k3 + k4)
    m[1] = f(t[n], V[n], r[n])
    k[1] = g(t[n], V[n], r[n])
    m[2] = f(t[n] + (h/2), V[n] + (h/2)m[1], r[n] + (h/2)k[1])
    k[2] = g(t[n] + (h/2), V[n] + (h/2)m[1], r[n] + (h/2)k[1])
    m[3] = f(t[n] + (h/2), V[n] + (h/2)m[2], r[n] + (h/2)k[2])
    k[3] = g(t[n] + (h/2), V[n] + (h/2)m[2], r[n] + (h/2)k[2])
    m[4] = f(t[n] + h, V[n] + h*m[3], r[n] + h*k[3])
    k[4] = g(t[n] + h, V[n] + h*m[3], r[n] + h*k[3])

    Parameters
    ----------
    k : float
        particle charge [Co] / particle mass [kg]
    RV : np.array([[x, y, z, Vx, Vy, Vz]])
        coordinates and velocities array [m], [m/s]
    dt : float
        timestep [s]
    E : np.array([Ex, Ey, Ez])
        values of electric field at current point [V/m]
    B : np.array([Bx, By, Bz])
        values of magnetic field at current point [T]

    Returns
    -------
    RV : np.array([[x, y, z, Vx, Vy, Vz]])
        new coordinates and velocities

    '''
    
    if np.any(np.isnan(B)): 
        print('NaN!!  B = ', B)
        print('   RV = ', RV)
    
    if np.any(np.isnan(E)): 
        print('NaN!!  E = ', E)
        print('   RV = ', RV)
    
    r = RV[0, :3]
    V = RV[0, 3:]

    m1 = f(k, E, V, B)
    k1 = g(V)

    fV2 = V + (dt / 2.) * m1
    gV2 = V + (dt / 2.) * m1
    m2 = f(k, E, fV2, B)
    k2 = g(gV2)

    fV3 = V + (dt / 2.) * m2
    gV3 = V + (dt / 2.) * m2
    m3 = f(k, E, fV3, B)
    k3 = g(gV3)

    fV4 = V + dt * m3
    gV4 = V + dt * m3
    m4 = f(k, E, fV4, B)
    k4 = g(gV4)

    V = V + (dt / 6.) * (m1 + (2. * m2) + (2. * m3) + m4)
    r = r + (dt / 6.) * (k1 + (2. * k2) + (2. * k3) + k4)

    RV = np.hstack((r, V))
    return RV

# %%
def save_png(fig, name, save_dir='output'):
    '''
    saves picture as name.png
    fig : array of figures to save
    name : array of picture names
    save_dir : directory used to store results
    '''

    # check wether directory exist and if not - create one
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
        print('LOG: {} directory created'.format(save_dir))
    print('LOG: Saving pictures to {}'.format(save_dir+'/'))
    for fig, name in zip(fig, name):
        # save fig with tight layout
        fig_savename = str(name + '.png')
        fig.savefig(save_dir + '/' + fig_savename, bbox_inches='tight')
        print('LOG: Figure ' + fig_savename + ' saved')


# %%
def _sym_np_arange(x0, x1, delta): 
    result = np.arange(x0, x1, delta)
    center = 0.5*(result[0] + result[-1])
    result = result - center
    return result

def _sym_grid(grid): 
    center = 0.5*(grid[:, 0, 0, 0]  + grid[:, -1, -1, -1])
    return grid - center

#%%
def argfind_rv(rrvv, rv): 
    #x = np.where(np.isclose(rrvv, rv0))
    bb = np.all(  np.isclose(rrvv, rv), axis=1  )
    ii = np.where(bb)
    i = ii[0][0] if len(ii[0]) > 0 else None
    return i

def find_fork(deltas, threshold=None):
    
    '''
    find 2 sec trajectories: 1st is higher than aim, 2nd is lower 
    calc linear parameter t: 0..1 : 0 if r_aim is at 1st traj, 1 if at 2nd
    '''
    try: 
        #set basis
        t = None    
        for i, (d1, d2) in enumerate( zip(deltas[0:-1], deltas[1:]) ):
            # check if last dots on different sides of aim
            if d1*d2 <= 0:
                # calc t
                t = abs(d1/(d1 - d2))
                return i, i+1, t
        
        if (t is None)and(threshold is not None): 
            i = np.argmin(np.abs(deltas))
            d = deltas[i] 
            # print('d =', d, 'threshold =', threshold)
            
            if np.abs(d) < threshold: 
                return i, i, 0.0
    except: 
        return None, None, None
    
    return None, None, None


def set_axes_param(ax, xlabel, ylabel, isequal=True):
    '''
    format axes
    '''
    ax.grid(True)
    ax.grid(which='major', color='tab:gray')  # draw primary grid
    ax.minorticks_on()  # make secondary ticks on axes
    # draw secondary grid
    ax.grid(which='minor', color='tab:gray', linestyle=':')
    ax.xaxis.set_tick_params(width=2)  # increase tick size
    ax.yaxis.set_tick_params(width=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if isequal:
        ax.axis('equal')
        
        