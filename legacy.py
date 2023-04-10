# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:08:38 2023

@author: Krohalev_OD
"""
#%% imports
import os
import numpy as np
import hibpplotlib as hbplot
from scipy.interpolate import RegularGridInterpolator

#%%
def read_B(Btor, Ipl, PF_dict, dirname='magfield', interp=True, plot=False):
    '''
    read Magnetic field values and create Bx, By, Bz, rho interpolants
    '''
    print('\n Reading Magnetic field')
    B_dict = {}
    for filename in os.listdir(dirname):
        if 'old' in filename:
            continue
        elif filename.endswith('.dat'):
            with open(dirname + '/' + filename, 'r') as f:
                volume_corner1 = [float(i) for i in f.readline().split()[0:3]]
                volume_corner2 = [float(i) for i in f.readline().split()[0:3]]
                resolution = float(f.readline().split()[0])
            continue
        elif 'Tor' in filename:
            print('Reading toroidal magnetic field...')
            B_read = np.load(dirname + '/' + filename) * Btor
            name = 'Tor'

        elif 'Plasm_{}MA'.format(int(Ipl)) in filename:
            print('Reading plasma field...')
            B_read = np.load(dirname + '/' + filename)  # * Ipl
            name = 'Plasm'

        else:
            name = filename.replace('magfield', '').replace('.npy', '')
            print('Reading {} magnetic field...'.format(name))
            Icir = PF_dict[name]
            print('Current = ', Icir)
            B_read = np.load(dirname + '/' + filename) * Icir

        B_dict[name] = B_read

    # create grid of points
    grid = np.mgrid[volume_corner1[0]:volume_corner2[0]:resolution,
                    volume_corner1[1]:volume_corner2[1]:resolution,
                    volume_corner1[2]:volume_corner2[2]:resolution]

    B = np.zeros_like(B_read)
    for key in B_dict.keys():
        B += B_dict[key]

#    cutoff = 10.0
#    Babs = np.linalg.norm(B, axis=1)
#    B[Babs > cutoff] = [np.nan, np.nan, np.nan]

    # plot B stream
    if plot:
        hbplot.plot_B_stream(B, volume_corner1, volume_corner2, resolution, grid,
                         plot_sep=False, dens=2.0)
    else: 
        print('B loaded without plotting')

    x = np.arange(volume_corner1[0], volume_corner2[0], resolution)
    y = np.arange(volume_corner1[1], volume_corner2[1], resolution)
    z = np.arange(volume_corner1[2], volume_corner2[2], resolution)
    Bx = B[:, 0].reshape(grid.shape[1:])
    By = B[:, 1].reshape(grid.shape[1:])
    Bz = B[:, 2].reshape(grid.shape[1:])
    if interp:
        # make an interpolation of B
        Bx_interp = RegularGridInterpolator((x, y, z), Bx, bounds_error = False)
        By_interp = RegularGridInterpolator((x, y, z), By, bounds_error = False)
        Bz_interp = RegularGridInterpolator((x, y, z), Bz, bounds_error = False)
        print('Interpolants for magnetic field created')
        B_list = [Bx_interp, By_interp, Bz_interp]
    else:
        B_list = [Bx, By, Bz]

    return B_list

#%%
B_old = read_B(1., 2., PF_dict, plot = True)

#%%
def return_B(r, Bin):
    '''
    interpolate Magnetic field at point r
    '''
    Bx_interp, By_interp, Bz_interp = Bin[0], Bin[1], Bin[2]
    Bout = np.c_[Bx_interp(r), By_interp(r), Bz_interp(r)]
    return Bout