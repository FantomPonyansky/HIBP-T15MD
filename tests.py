# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 12:14:57 2023

@author: Krohalev_OD
"""

#%% imports
import numpy as np
import os
import errno
import pickle as pc
import copy
from matplotlib import path
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from itertools import cycle
import optimizers

import hibplib as hb

import hibpcalc.fields as fields
import hibpcalc.misc as misc
import hibpcalc.geomfunc as gf
import define_geometry as defgeom


#%% set voltages
# UA2 voltages
UA2min, UA2max, dUA2 = 32., 32., 2. #12., 12., 2.  #0., 34., 2.  # -3, 33., 3.  # -3., 30., 3.
NA2_points = 10

# B2 plates voltage
UB2, dUB2 = 0.0, 5.0  # 10.  # [kV], [kV/m]

# B3 voltages
UB3, dUB3 = 0.0, 10  # [kV], [kV/m]

# A3 voltages
UA3, dUA3 = 0.0, 7.0  # [kV], [kV/m]

# A4 voltages
UA4, dUA4 = 0.0, 2.0  # [kV], [kV/m]

# %% Define Geometry
geomT15 = defgeom.define_geometry(analyzer=1)
r0 = geomT15.r_dict['r0']  # trajectory starting point

# angles of aim plane normal [deg]
alpha_aim = 0.
beta_aim = 0.
stop_plane_n = gf.calc_vector(1.0, alpha_aim, beta_aim)

Ipl = 1.0 # MA
Btor = 1.0 # T

_U_dict = {  'A2': 50.0,
             'B2': -4.0,
             'A3': 3.0,
             'B3': 4.0,
             'A4': 5.0,
             'an': 107.0}

# %% Load Electric Field
E_slow = {}
E_fast = {}
_E_slow = {}
_E_fast = {}
# load E for primary beamline
try:
    fields.read_plates ('prim', geomT15,  E_slow,  E_fast, hb.createplates)
    fields._read_plates('prim', geomT15, _E_slow, _E_fast, hb.createplates)
    print('\n Primary Beamline loaded')
except FileNotFoundError:
    print('\n Primary Beamline NOT FOUND')

# load E for secondary beamline
try:
    fields.read_plates ('sec', geomT15,  E_slow,  E_fast, hb.createplates)
    fields._read_plates('sec', geomT15, _E_slow, _E_fast, hb.createplates)

    # add diafragm for A3 plates to Geometry
    hb.add_diafragm(geomT15, 'A3', 'A3d', diaf_width=0.05)
    hb.add_diafragm(geomT15, 'A4', 'A4d', diaf_width=0.05)
    print('\n Secondary Beamline loaded')
except FileNotFoundError as e:
    print(e)
    print('\n Secondary Beamline NOT FOUND')

# test for plates.r
for pl_nm in ['B2', 'A2', 'A3', 'B3', 'A4', 'an']:
    r = geomT15.plates_dict[pl_nm].r
    r_ = geomT15.r_dict[pl_nm]
    if not np.all( np.isclose(r, r_) ): 
        print('TEST FAILED:  ', pl_nm+'.r=', r,     'geom.r_dict[%s]='%pl_nm, r_)

# test for plates transformation matrices
for pl_nm in ['B2', 'A2', 'A3', 'B3', 'A4', 'an']:
    pl = geomT15.plates_dict[pl_nm]
    r = pl.r
    edge00 = pl.edges[0, 0]

    # transformation to the original coordinates
    v = pl.inv_rotation_mx.dot(edge00 - pl.r)  # !!!  inverse
    v_ = pl.return_r_in_original_coordinates(geomT15, edge00)
    if not np.all( np.isclose(v, v_) ): 
        print("TEST FAILED:  ', pl_nm, ' inv_rotation_mx does not work")
    
    # transformation from the original coordinates to the actual ones
    v = pl.rotation_mx.dot(v) + pl.r  # !!!  direct
    v_ = gf.rotate3(v_, pl.angles, pl.beamline_angles, inverse=False)
    v_ = v_ + pl.r
    if (   not np.all( np.isclose(v, v_))    )or(    not np.all( np.isclose(v, edge00))     ): 
        print("TEST FAILED:  ', pl_nm, ' rotation_mx does not work")
              

# %% Analyzer parameters
if 'an' in geomT15.plates_dict.keys():
    # Analyzer G
    an_G = geomT15.plates_dict['an'].G
    # add detector coords to dictionary
    an_edges = geomT15.plates_dict['an'].det_edges
    geomT15.r_dict['det'] = an_edges[ an_edges.shape[0]//2 ] [0] 
else:
    an_G = 1.
    print('\nNO Analyzer')

# %% Load Magnetic Field
_pf_coils = fields.import_PFcoils('PFCoils.dat')
_PF_dict = fields.import_PFcur('{}MA_sn.txt'.format(int(abs(Ipl))), _pf_coils)
if 'B_' not in locals():
    B_ = fields.read_B_new(Btor, Ipl, _PF_dict, dirname='magfield')
    B_.cure_artefacts_from_filaments()
else: 
    print('B already loaded')

#%%
r_test = geomT15.r_dict['A3'] + np.array([0.03, 0.03, 0.03])



__Es_ = fields._return_E(r_test, _E_slow, _U_dict, geomT15)
__Ef_ = fields._return_E(r_test, _E_fast, _U_dict, geomT15)

__Es = fields.return_E(r_test, E_slow, _U_dict, geomT15)
#__Ef = fields.return_E(r_test, E_fast, _U_dict, geomT15)

print(__Es_, __Ef_, __Es)

B2 = geomT15.plates_dict['B2']
A2 = geomT15.plates_dict['A2']

A3 = geomT15.plates_dict['A3']
B3 = geomT15.plates_dict['B3']
A4 = geomT15.plates_dict['A4']


r_ = A3.return_r_in_original_coordinates(geomT15, r_test) 
print(E_fast['A3'](r_))
print([ E_slow['A3'][i](r_) for i in [0, 1, 2]] )

#%%
if False: 
    for x in np.linspace(0.0, 5.0, 100):
       for y in np.linspace(-3.0, 2.0, 100):
           for z in [0.0, -0.1, 1.1]:
            r_test = np.array([x, y, z])   
            __Es_ = fields._return_E(r_test, _E_slow, _U_dict, geomT15)
            __Ef_ = fields._return_E(r_test, _E_fast, _U_dict, geomT15)
            __Es = fields.return_E(r_test, E_slow, _U_dict, geomT15)
            if not np.any( np.isclose(__Es_, __Ef_) ):
                print('?', r_test, '  slow_=', __Es_, '  fast_=', __Ef_)
            if not np.any( np.isclose(__Es, __Es_) ):
                print('??', r_test, '  slow=', __Es, '  slow_=', __Es_)

if False: 
    #for x in np.linspace(0.0, 5.0, 100):
       #for y in np.linspace(-3.0, 2.0, 100):
           x = 2.5
           y = 0.0
           #for z in np.linspace(0.499, 0.521, 1000000):
           for z in np.linspace(-1.21, -1.19, 10000):
               r_test = np.array([x, y, z])   
               try: 
                   _B_ = B_(r_test)
               except: 
                   print(r_test)


#%%

if False:             
    x, y, z = geomT15.r_dict['A3']
    for x in np.linspace(2.0, 5.0, 30000):
        r_test = np.array([x, y+.1, z])   
        __Es_ = fields._return_E(r_test, _E_slow, _U_dict, geomT15)
        __Ef_ = fields._return_E(r_test, _E_fast, _U_dict, geomT15)
        __Es = fields.return_E(r_test, E_slow, _U_dict, geomT15)
        if not np.any( np.isclose(__Es_, __Ef_) ):
            print('?', r_test, '  slow_=', __Es_, '  fast_=', __Ef_)
        if not np.any( np.isclose(__Es, __Es_) ):
            print('??', r_test, '  slow=', __Es, '  slow_=', __Es_)

#%%
    
error_r_list = [
        np.array([ 2.82828283, -0.22222222,  0.        ]),
        np.array([ 2.82828283, -0.22222222, -0.1       ]),
        np.array([ 2.97979798, -0.22222222, -0.1       ]),
        np.array([ 3.18181818, -0.12121212,  0.        ]),
        np.array([ 3.28282828, -0.02020202, -0.1       ]),
        np.array([ 3.53535354,  0.03030303, -0.1       ])
    ]

error_r_list2 = [
        np.array([ 2.99573319, -0.14823619, -0.05        ]),
        np.array([ 3.15033834, -0.14823619, -0.05        ])
    ]


for pl_name in ['A2', 'B2', 'A3', 'B3', 'A4', 'an']: 
    pl = geomT15.plates_dict[pl_name]
    r_ = pl.return_r_in_original_coordinates(geomT15, error_r_list2[0]) 
    print(pl_name)
    print(r_)
    print('fast: ', _E_fast[pl.name](r_))
    try: 
        print('slow: ', [ _E_slow[pl.name][i](r_) for i in [0, 1, 2]] )
    except ValueError: 
        print('slow: 0, 0, 0')

def original_plates(pl_name, geom): 
    _pl = copy.deepcopy(geom.plates_dict[pl_name])
    
    dr = geomT15.r_dict[pl_name]
    angles = _pl.angles
    beamline_angles = _pl.beamline_angles
    
    _pl.shift(-dr)
    _pl.rotate(angles, beamline_angles, inverse=True)
    return _pl

_A3 = original_plates('A3', geomT15)
print(_A3.edges)
print(_A3.return_r_in_original_coordinates(geomT15, error_r_list2[0]))



if False: 
    r_A3 = geomT15.r_dict['A3']
    z_ = r_[2]
    for x in np.linspace(r_A3[0] - 0.3, r_A3[0] + 0.3, 40): #np.linspace(2.4, 4.8, 100):
       for y in np.linspace(r_A3[1] - 0.3, r_A3[1] + 0.3, 40):  #np.linspace(-0.3, 0.5, 100):
           for z in np.linspace(z_ - 0.3, z_ + 0.3, 50):  
            r_test = np.array([x, y, z])   
            if A3.contains_point(geomT15, r_test): 
                gf.plot_point(r_test)
    geomT15.plot(plt.gca())            

#geomT15.plot(plt.gca())

#%%
for pl_name in ['A2', 'B2', 'A3', 'B3', 'A4', 'an']: 
    pl = geomT15.plates_dict[pl_name]
    self_r = pl.r
    geom_r = geomT15.r_dict[pl_name]
    print(pl_name, self_r, geom_r)

    rct = pl.front_rect()
    gf.plot_rect(rct, tangage=True)
geomT15.plot(plt.gca())

#%%
#from requests import get

#response = get("http://www-fusion.ciemat.es/cgi-bin/TJII_getmagn.cgi?config=100_44_64&x=1.7&y=0&z=0.05")
#print(response.content)


#%%
import sys
from time import sleep

# print(sys.stdout._buffer.read())
# sys.stdout._buffer.write("test")
# sys.stdout._buffer.flush()
# sys.stdout._buffer.write("test")
# sys.stdout._buffer.tell()
# x = sys.stdout._buffer.getvalue()

JUMP_LEFT_SEQ = '\u001b[100D' # 

def loading():
    for i in range(0, 16): 
        sleep(0.5)
        #print(JUMP_LEFT_SEQ, end='')
        print('\r', end='')
        print(f'Progress: {i:0>3}%', end='')
        sys.stdout.flush()
    print('\r', end='')    
    print('!             ')

#loading()

#%%

print('\a')  # beep

#%%

import time

__U_dict = { 'A2': 50.0,
             'B2': -4.0,
             'A3': 3.0,
             'B3': 4.0,
             'A4': 5.0,
             'an': 107.0}

for pl_name in ['A2', 'B2', 'A3', 'B3', 'A4', 'an']: 
    #r_test = geomT15.r_dict[pl_name] + np.array([0.01, 0.01, 0.01])    
    N = pl.edges[1].shape[0]*pl.edges[1].shape[1]
    for r in pl.edges[0].reshape(N): 
        r_test = r + np.array([0.01, 0.01, 0.01])
        __E1 = fields._return_E (r_test, _E_fast, __U_dict, geomT15)
        __E2 = fields.__return_E(r_test, _E_fast, __U_dict, geomT15)
        if not np.all( np.isclose(__E1, __E2, atol=1e-10) ): 
            print('TEST FAILED: __return_E <> _return_E')


pl = geomT15.plate_dict['A3']

for x in np.lnspace(-5.0, ): 
    #r_test = geomT15.r_dict[pl_name] + np.array([0.01, 0.01, 0.01])    
    N = pl.edges[1].shape[0]*pl.edges[1].shape[1]
    for r in pl.edges[0].reshape(N): 
        r_test = r + np.array([0.01, 0.01, 0.01])
        __E1 = fields._return_E (r_test, _E_fast, __U_dict, geomT15)
        __E2 = fields.__return_E(r_test, _E_fast, __U_dict, geomT15)
        if not np.all( np.isclose(__E1, __E2, atol=1e-10) ): 
            print('TEST FAILED: __return_E <> _return_E')


t0 = time.time()
r_test = geomT15.r_dict['an'] + np.array([0.01, 0.01, 0.01])    
for i in range(100000):
    r_test = r_test + np.array([0.001, 0.001, 0.001])    
    __E2 = fields.__return_E(r_test, _E_fast, __U_dict, geomT15)
    
t1 = time.time()
print(t1 - t0)

t0 = time.time()
r_test = geomT15.r_dict['A2'] + np.array([-0.01, 0.01, 0.01])    
for i in range(100000):
    r_test = r_test + np.array([0.001, 0.001, 0.001])
    __E2 = fields.return_E(r_test, _E_slow, __U_dict, geomT15)
    
t1 = time.time()
print(t1 - t0)


#%%


def recalc_gabarits(self, E_interp): 
    g = E_interp[self.name].grid
    edges = [ g[:,  0, 0, 0],  g[:,  0, 0, -1], g[:,  0, -1, 0], g[:,  0, -1, -1], 
              g[:, -1, 0, 0],  g[:, -1, 0, -1], g[:, -1, -1, 0], g[:, -1, -1, -1] ]
    
    edges_ = np.array( [  self.rotation_mx.dot(e) + self.r for e in edges ] )
    return edges_
    mins = np.min(edges_, axis=0)
    maxs = np.max(edges_, axis=0)
    
edges_ = recalc_gabarits(A3, _E_fast)    