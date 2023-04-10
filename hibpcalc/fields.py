# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:09:17 2023

@author: Krohalev_OD
"""
#%% imports

import os
import copy
import errno
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

from hibpcalc.geomfunc import rotate, rotate3
import hibpcalc.misc as misc

# %%
def import_Bflux(filename):
    '''
    import magnetic flux from Tokameq file
    '''
    with open(filename, 'r') as f:
        data = f.readlines()

    # R coordinate corresponds to X, Z coordinate corresponds to Y
    NrNz = []
    for i in data[2].strip().split():
        if i.isdigit():
            NrNz.append(i)
    Nx = int(NrNz[0]) + 1
    Ny = int(NrNz[1]) + 1

    for i in range(len(data)):
        if ' '.join(data[i].strip().split()[:4]) == 'Flux at the boundary':
            bound_flux = float(data[i].strip().split()[-1])
        if data[i].strip() == 'Poloidal flux F(r,z)':
            index = i

    x_vals = [float(r) for r in data[index+1].strip().split()[1:]]
    x_vals = np.array(x_vals)

    Psi_data = [i.strip().split() for i in data[index+2:index+2+Ny]]
    Psi_vals = []
    y_vals = []
    for line in Psi_data:
        y_vals.append(float(line[0]))
        Psi_vals.append([float(j) for j in line[1:]])

    y_vals = np.array(y_vals)
    Psi_vals = np.array(Psi_vals)
    return Psi_vals, x_vals, y_vals, bound_flux
#%% Magfield interpolator
def array_to_value(array):
    result = array[0]
    if np.all( np.isclose(array, result, atol=1e-8) ): 
        return result
    else:     
        raise Exception()
#%%

class FieldInterpolator():
    '''
    Interpolates vector on 3d grid with equal steps
    '''
    
    def __init__(self, grid, list_Fx_Fy_Fz, default=np.nan):
        self.grid = grid
        self.res = self.grid[:, 1, 1, 1] - self.grid[:, 0, 0, 0]
        self.res = array_to_value(self.res)
        self.list_Fx_Fy_Fz = list_Fx_Fy_Fz
        self.volume_corner1 = self.grid[:, 0, 0, 0]
        self.volume_corner2_ = self.grid[:, -1, -1, -1]
        self.volume_corner2 = self.volume_corner2_ + self.res
        self.default_value = default
    
    def __call__(self, point):
        #if point is outside the volume - return array of np.nan-s
        if np.any(point - self.volume_corner1 <= 0.0) or np.any(self.volume_corner2_ - point <= 0.0):
            return np.full((1, 3), self.default_value) # [0]
        
        #finding indexes of left corner of volume with point
        indexes_float = (point - self.volume_corner1)/self.res // 1
        indexes = [[0]*3]*8
        for i in range(3):
            indexes[0][i] = int(indexes_float[i])
        
        # finding weights for all dots close to point
        '''
        delta_x = [x2 - x, x - x1]
        point = [x, y, z]
        '''
        
        i00 = indexes[0][0]
        j01 = indexes[0][1]
        k02 = indexes[0][2]
        
        left_bottom = self.grid[:, i00,     j01,     k02]        
        right_top   = self.grid[:, i00 + 1, j01 + 1, k02 + 1]
        #delta = (right_top - point,   point - left_bottom)
        
        delta_x = [right_top[0] - point[0],   point[0] - left_bottom[0] ] 
        delta_y = [right_top[1] - point[1],   point[1] - left_bottom[1] ]  
        delta_z = [right_top[2] - point[2],   point[2] - left_bottom[2] ] 
        
        res_cubic = self.res**3
        number = 0
        weights = [0.]*8
        for i in range(2):
            for j in range (2):
                for k in range(2):
                    weights[number] = delta_x[i]*delta_y[j]*delta_z[k]/res_cubic
                    number += 1
        
        #finding interpolation
        Fx = 0.0
        Fy = 0.0
        Fz = 0.0
        
        _Fx = self.list_Fx_Fy_Fz[0]
        _Fy = self.list_Fx_Fy_Fz[1]
        _Fz = self.list_Fx_Fy_Fz[2]
        for ijk in range(8):
            i = ( ijk >> 2 ) % 2
            j = ( ijk >> 1 ) % 2
            k =   ijk % 2 
            
            Fx += weights[ijk]* _Fx[i00 + i, j01 + j, k02 + k]
            Fy += weights[ijk]* _Fy[i00 + i, j01 + j, k02 + k]
            Fz += weights[ijk]* _Fz[i00 + i, j01 + j, k02 + k]

        return np.array([[Fx, Fy, Fz]])

    def plot(self, color='r', dens=1.0, plot_sep=True):
        '''
        stream plot of magnetic field
        '''
    
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
        misc.set_axes_param(ax1, 'X (m)', 'Y (m)')
        misc.set_axes_param(ax2, 'X (m)', 'Z (m)')
        
        x = np.arange(self.volume_corner1[0], self.volume_corner2[0], self.res)
        y = np.arange(self.volume_corner1[1], self.volume_corner2[1], self.res)
        z = np.arange(self.volume_corner1[2], self.volume_corner2[2], self.res)
    
        Fx, Fy, Fz = self.list_Fx_Fy_Fz 
    
        # choose z position
        z_cut = np.where(abs(z) < 0.001)[0][0]  # Bx.shape[2]//2
        # choose y position
        y_cut = np.where(abs(y) < 0.001)[0][0]  # Bx.shape[1]//2
    
        ax1.streamplot(x, y, Fx[:, :, z_cut].swapaxes(0, 1), Fy[:, :, z_cut].swapaxes(0, 1), color=color, density=dens)
        ax2.streamplot(x, z, Fx[:, y_cut, :].swapaxes(0, 1), Fz[:, y_cut, :].swapaxes(0, 1), color=color, density=dens)
        plt.show()
    
    def cure_artefacts_from_filaments(self):  
        '''
        WARNING: Costyl'

        '''    
        Fx, Fy, Fz = self.list_Fx_Fy_Fz  # ???        
        xmask = np.abs(Fx) > 1.0
        Fx[xmask] = 0.0
        
        ymask = np.abs(Fy) > 1.0
        Fy[ymask] = 0.0
    
# %% poloidal field coils
def import_PFcoils(filename):
    '''
    import a dictionary with poloidal field coils parameters
    {'NAME': (x center, y center, width along x, width along y [m],
               current [MA-turn], N turns)}
    Andreev, VANT 2014, No.3
    '''
    d = {}  # defaultdict(list)
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            lineList = line.split(', ')
            key, val = lineList[0], tuple([float(i) for i in lineList[1:]])
            d[key] = val
    return d


def import_PFcur(filename, pf_coils):
    '''
    Creates dictionary with coils names and currents from TOKAMEQ file
    filename : Tokameqs filename
    pf_coils : coil dict (we only take keys)
    '''
    with open(filename, 'r') as f:
        data = f.readlines()  # read tokameq file
    PF_dict = {}  # Here we will store coils names and currents
    pf_names = list(pf_coils)  # get coils names
    n_coil = 0  # will be used for getting correct coil name
    for i in range(len(data)):
        if data[i].strip() == 'External currents:':
            n_line = i + 2  # skip 2 lines and read from the third
            break
    while float(data[n_line].strip().split()[3]) != 0:
        key = pf_names[n_coil]
        val = data[n_line].strip().split()[3]
        PF_dict[key] = float(val)
        n_line += 1
        n_coil += 1

    return PF_dict

#%% read B from file

def read_B_new(Btor, Ipl, PF_dict, dirname='magfield'):
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

    Bx = B[:, 0].reshape(grid.shape[1:])
    By = B[:, 1].reshape(grid.shape[1:])
    Bz = B[:, 2].reshape(grid.shape[1:])

    # make an interpolation of B
    B_interp = FieldInterpolator(grid, [Bx, By, Bz])
    B_interp.cure_artefacts_from_filaments()
    print('Interpolants for magnetic field created')

    return B_interp
# %% Electric field & plates
def plate_flags(range_x, range_y, range_z, U,
                plts_geom, plts_angles, plts_center):
    '''
    calculate plates cooedinates and boolean arrays for plates
    '''
    length, width, thick, gap, l_sw = plts_geom
    gamma, alpha_sw = plts_angles
    r_sweep_up = np.array([-length/2 + l_sw, gap/2., 0])
    r_sweep_lp = np.array([-length/2 + l_sw, -gap/2., 0])
    # Geometry in system based on central point between plates
    # upper plate
    UP1 = np.array([-length/2., gap/2. + thick, width/2.])
    UP2 = np.array([-length/2., gap/2. + thick, -width/2.])
    UP1sw = np.array([-length/2. + l_sw, gap/2. + thick, width/2.])
    UP2sw = np.array([-length/2. + l_sw, gap/2. + thick, -width/2.])
    UP3 = np.array([length/2., gap/2. + thick, -width/2.])
    UP4 = np.array([length/2., gap/2. + thick, width/2.])
    UP5 = np.array([-length/2., gap/2., width/2.])
    UP6 = np.array([-length/2., gap/2., -width/2.])
    UP5sw = np.array([-length/2. + l_sw, gap/2., width/2.])
    UP6sw = np.array([-length/2. + l_sw, gap/2., -width/2.])
    UP7 = np.array([length/2., gap/2., -width/2.])
    UP8 = np.array([length/2., gap/2., width/2.])
    if abs(alpha_sw) > 1e-2:
        UP1 = UP1sw + rotate(UP1 - UP1sw, axis=(0, 0, 1), deg=-alpha_sw)
        UP2 = UP2sw + rotate(UP2 - UP2sw, axis=(0, 0, 1), deg=-alpha_sw)
        UP5 = UP5sw + rotate(UP5 - UP5sw, axis=(0, 0, 1), deg=-alpha_sw)
        UP6 = UP6sw + rotate(UP6 - UP6sw, axis=(0, 0, 1), deg=-alpha_sw)
        # points are sorted clockwise
        UP_full = np.array([UP1sw, UP1, UP2, UP2sw, UP3, UP4,
                            UP5sw, UP5, UP6, UP6sw, UP7, UP8])
    else:
        UP_full = np.array([UP1, UP2, UP3, UP4, UP5, UP6, UP7, UP8])
    UP_rotated = UP_full.copy()
    for i in range(UP_full.shape[0]):
        UP_rotated[i, :] = rotate(UP_rotated[i, :], axis=(1, 0, 0), deg=gamma)
        # shift coords center
        UP_rotated[i, :] += plts_center

    # lower plate
    LP1 = np.array([-length/2., -gap/2. - thick, width/2.])
    LP2 = np.array([-length/2., -gap/2. - thick, -width/2.])
    LP1sw = np.array([-length/2. + l_sw, -gap/2. - thick, width/2.])
    LP2sw = np.array([-length/2. + l_sw, -gap/2. - thick, -width/2.])
    LP3 = np.array([length/2., -gap/2. - thick, -width/2.])
    LP4 = np.array([length/2., -gap/2. - thick, width/2.])
    LP5 = np.array([-length/2., -gap/2., width/2.])
    LP6 = np.array([-length/2., -gap/2., -width/2.])
    LP5sw = np.array([-length/2. + l_sw, -gap/2., width/2.])
    LP6sw = np.array([-length/2. + l_sw, -gap/2., -width/2.])
    LP7 = np.array([length/2., -gap/2., -width/2.])
    LP8 = np.array([length/2., -gap/2., width/2.])
    if abs(alpha_sw) > 1e-2:
        LP1 = LP1sw + rotate(LP1 - LP1sw, axis=(0, 0, 1), deg=alpha_sw)
        LP2 = LP2sw + rotate(LP2 - LP2sw, axis=(0, 0, 1), deg=alpha_sw)
        LP5 = LP5sw + rotate(LP5 - LP5sw, axis=(0, 0, 1), deg=alpha_sw)
        LP6 = LP6sw + rotate(LP6 - LP6sw, axis=(0, 0, 1), deg=alpha_sw)
        # points are sorted clockwise
        LP_full = np.array([LP1sw, LP1, LP2, LP2sw, LP3, LP4,
                            LP5sw, LP5, LP6, LP6sw, LP7, LP8])
    else:
        LP_full = np.array([LP1, LP2, LP3, LP4, LP5, LP6, LP7, LP8])
    LP_rotated = LP_full.copy()
    for i in range(LP_full.shape[0]):
        LP_rotated[i, :] = rotate(LP_rotated[i, :], axis=(1, 0, 0), deg=gamma)
        # shift coords center
        LP_rotated[i, :] += plts_center

    # Find coords of 'cubes' containing each plate
    upper_cube = np.array([np.min(UP_rotated, axis=0),
                           np.max(UP_rotated, axis=0)])
    lower_cube = np.array([np.min(LP_rotated, axis=0),
                           np.max(LP_rotated, axis=0)])

    # create mask for plates
    upper_plate_flag = np.full_like(U, False, dtype=bool)
    lower_plate_flag = np.full_like(U, False, dtype=bool)
    for i in range(range_x.shape[0]):
        for j in range(range_y.shape[0]):
            for k in range(range_z.shape[0]):
                x = range_x[i]
                y = range_y[j]
                z = range_z[k]
                # check upper cube
                if (x >= upper_cube[0, 0]) and (x <= upper_cube[1, 0]) and \
                   (y >= upper_cube[0, 1]) and (y <= upper_cube[1, 1]) and \
                   (z >= upper_cube[0, 2]) and (z <= upper_cube[1, 2]):
                    r = np.array([x, y, z]) - plts_center
                    # inverse rotation
                    r_rot = rotate(r, axis=(1, 0, 0), deg=-gamma)
                    if r_rot[0] <= -length/2 + l_sw:
                        r_rot = r_sweep_up + rotate(r_rot - r_sweep_up,
                                                 axis=(0, 0, 1), deg=alpha_sw)
                    # define masks for upper and lower plates
                    upper_plate_flag[i, j, k] = (r_rot[0] >= -length/2.) and \
                        (r_rot[0] <= length/2.) and (r_rot[2] >= -width/2.) and \
                        (r_rot[2] <= width/2.) and (r_rot[1] >= gap/2.) and \
                        (r_rot[1] <= gap/2. + thick)
                # check lower cube
                if (x >= lower_cube[0, 0]) and (x <= lower_cube[1, 0]) and \
                   (y >= lower_cube[0, 1]) and (y <= lower_cube[1, 1]) and \
                   (z >= lower_cube[0, 2]) and (z <= lower_cube[1, 2]):
                    r = np.array([x, y, z]) - plts_center
                    # inverse rotation
                    r_rot = rotate(r, axis=(1, 0, 0), deg=-gamma)
                    if r_rot[0] <= -length/2 + l_sw:
                        r_rot = r_sweep_lp + rotate(r_rot - r_sweep_lp,
                                                 axis=(0, 0, 1), deg=-alpha_sw)
                    # define masks for upper and lower plates
                    lower_plate_flag[i, j, k] = (r_rot[0] >= -length/2.) and \
                        (r_rot[0] <= length/2.) and (r_rot[2] >= -width/2.) and \
                        (r_rot[2] <= width/2.) and \
                        (r_rot[1] >= -gap/2. - thick) and \
                        (r_rot[1] <= -gap/2.)

    return UP_rotated, LP_rotated, upper_plate_flag, lower_plate_flag

#%%
def return_E(r, Ein, U, geom):
    '''
    take dot and try to interpolate electiric field
    Ein : dict of interpolants for Ex, Ey, Ez
    U : dict with plates voltage values
    '''
    Etotal = np.zeros(3)
    # do not check plates while particle is in plasma
    if r[0] < geom.r_dict['aim'][0]-0.15 and r[1] < geom.r_dict['port_in'][1]:  # <reonid>
    #if r[0] < geom.r_dict['aim'][0]-0.15 and r[1] < geom.r_dict['port_in'][1]:
        return Etotal
    # go through all the plates
    for key in Ein.keys():
        # shift the center of coord system
        r_new = r - geom.r_dict[key]
        # get angles
        angles = copy.deepcopy( geom.plates_dict[key].angles )
        beamline_angles = copy.deepcopy( geom.plates_dict[key].beamline_angles )
        # rotate point to the coord system of plates
        r_new = rotate3(r_new, angles, beamline_angles, inverse=True)
        # interpolate Electric field
        Etemp = np.zeros(3)
        try:
            Etemp[0] = Ein[key][0](r_new) * U[key]
            Etemp[1] = Ein[key][1](r_new) * U[key]
            Etemp[2] = Ein[key][2](r_new) * U[key]
            # rotate Etemp
            Etemp = rotate3(Etemp, angles, beamline_angles, inverse=False)
            # add the result to total E field
            Etotal += Etemp
        except (ValueError, IndexError):
            continue
    return Etotal

#%%

def save_E(beamline, plts_name, Ex, Ey, Ez, plts_angles, plts_geom,
           domain, an_params, plate1, plate2, dirname='elecfield'):
    '''
    save Ex, Ey, Ez arrays to file
    '''
    dirname = dirname + '/' + beamline

    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname, 0o700)
            print("Directory ", dirname, " created ")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    fname = plts_name + '_geometry.dat'
    # erases data from file before writing
    open(dirname + '/' + fname, 'w').close()
    with open(dirname + '/' + fname, 'w') as myfile:
        myfile.write(np.array2string(plts_geom)[1:-1] + ' # plate\'s length, width, thic, gap and l_sweeped\n')
        myfile.write(np.array2string(plts_angles)[1:-1] + ' # plate\'s gamma and alpha_sweep angles\n')
        myfile.write(np.array2string(domain, max_line_width=200)[1:-1] + ' # xmin, xmax, ymin, ymax, zmin, zmax, delta\n')

        if plts_name == 'an':
            myfile.write(np.array2string(an_params, max_line_width=200)[1:-1] + ' # n_slits, slit_dist, slit_w, G, theta, XD, YD1, YD2\n')
        for i in range(plate1.shape[0]):
            myfile.write(np.array2string(plate1[i], precision=4)[1:-1] + ' # 1st plate rotated\n')
        for i in range(plate2.shape[0]):
            myfile.write(np.array2string(plate2[i], precision=4)[1:-1] + ' # 2nd plate rotated\n')

    np.save(dirname + '/' + plts_name + '_Ex', Ex)
    np.save(dirname + '/' + plts_name + '_Ey', Ey)
    np.save(dirname + '/' + plts_name + '_Ez', Ez)

    print('Electric field saved, ' + plts_name + '\n')


# def createplates(name, beamline, an_params=None): 
#     if name == 'an':
#         # create new Analyzer object
#         plts = Analyzer(name, beamline)
#         plts.add_slits(an_params)
#     else:
#         plts = Plates(name, beamline)
#     return plts

def read_plates(beamline, geom, E, E_fast, createplates, dirname='elecfield'):
    '''
    read Electric field and plates geometry
    '''
    dirname = dirname + '/' + beamline

    # list of all *.dat files
    file_list = [file for file in os.listdir(dirname) if file.endswith('dat')]
    # push analyzer to the end of the list
    file_list.sort(key=lambda s: s[1:])  # A3, B3, A4, an
    # file_list.sort(key=lambda s: s.startswith('an'))

    for filename in file_list:
        plts_name = filename[0:2]
        print('\n Reading geometry {} ...'.format(plts_name))
        r_new = geom.r_dict[plts_name]
        print('position ', r_new)
        # angles of plates, will be modified later
        plts_angles = copy.deepcopy(geom.angles_dict[plts_name])
        # beamline angles
        beamline_angles = copy.deepcopy(geom.angles_dict[plts_name])
        # read plates parameters from file
        edges_list = []
        an_params = None
        with open(dirname + '/' + filename, 'r') as f:
            # read plates geometry, first remove comments '#', then convert to float
            plts_geom = [float(i) for i in f.readline().split('#')[0].split()]
            # read gamma angle (0 for Alpha and 90 for Beta plates)
            gamma = float(f.readline().split()[0])
            # xmin, xmax, ymin, ymax, zmin, zmax, delta
            domain = [float(i) for i in f.readline().split()[0:7]]
            if plts_name == 'an':
                an_params = [float(i) for i in f.readline().split()[0:8]]
                theta_an = an_params[4]  # analyzer entrance angle
                plts_angles[0] = plts_angles[0] - theta_an
            for line in f:
                # read plates Upper and Lowe plate coords, x,y,z
                edges_list.append([float(i) for i in line.split()[0:3]])
        edges_list = np.array(edges_list)

        plts = createplates(plts_name, beamline, an_params)
        # add edges to plates object
        index = int(edges_list.shape[0] / 2)
        plts.set_edges(np.array([edges_list[0:index, :],
                                 edges_list[index:, :]]))
        # rotate plates edges
        plts.rotate(plts_angles, beamline_angles)
        # shift coords center and put into a dictionary
        plts.shift(r_new)
        # add plates to dictionary
        geom.plates_dict[plts_name] = plts

        # read Electric field arrays
        Ex = np.load(dirname + '/' + plts_name + '_Ex.npy')
        Ey = np.load(dirname + '/' + plts_name + '_Ey.npy')
        Ez = np.load(dirname + '/' + plts_name + '_Ez.npy')

        non_sym_grid = np.mgrid[domain[0]:domain[1]:domain[6],
                                domain[2]:domain[3]:domain[6],
                                domain[4]:domain[5]:domain[6]]
    
        
        #x = misc._sym_np_arange(domain[0], domain[1], domain[6])  # + r_new[0]
        #y = misc._sym_np_arange(domain[2], domain[3], domain[6])  # + r_new[1]
        #z = misc._sym_np_arange(domain[4], domain[5], domain[6])  # + r_new[2]

        x = np.arange(domain[0], domain[1], domain[6])  # + r_new[0]
        y = np.arange(domain[2], domain[3], domain[6])  # + r_new[1]
        z = np.arange(domain[4], domain[5], domain[6])  # + r_new[2]


        # make interpolation for Ex, Ey, Ez
        Ex_interp = RegularGridInterpolator((x, y, z), Ex)
        Ey_interp = RegularGridInterpolator((x, y, z), Ey)
        Ez_interp = RegularGridInterpolator((x, y, z), Ez)
        E_read = [Ex_interp, Ey_interp, Ez_interp]

        E[plts_name] = E_read
        E_fast[plts_name] = FieldInterpolator(non_sym_grid, [Ex, Ey, Ez])

    return


#%% new E field
def plates_list(path): 
    file_list = [file for file in os.listdir(path) if file.endswith('dat')]
    # push analyzer to the end of the list
    file_list.sort(key=lambda s: s[1:])  # A3, B3, A4, an
    # file_list.sort(key=lambda s: s.startswith('an'))
    plate_names = [fname[0:2] for fname in file_list]
    return plate_names 

def read_plate_geom(path, plate_name): 
    filename = path + '/' + plate_name + '_geometry.dat'
    edges_list = []
    an_params = None

    with open(filename, 'r') as f:
        # read plates geometry, first remove comments '#', then convert to float
        plts_geom = [float(i) for i in f.readline().split('#')[0].split()]  #     [0]
        # read gamma angle (0 for Alpha and 90 for Beta plates)
        gamma = float(f.readline().split()[0])                              #     [1]
        # xmin, xmax, ymin, ymax, zmin, zmax, delta
        domain = [float(i) for i in f.readline().split()[0:7]]              #     [2]

        if plate_name == 'an':
            an_params = [float(i) for i in f.readline().split()[0:8]]       #     [3]    
            #theta_an = an_params[4]  # analyzer entrance angle
            #plts_angles[0] = plts_angles[0] - theta_an
        for line in f:
            # read plates Upper and Lowe plate coords, x,y,z
            edges_list.append([float(i) for i in line.split()[0:3]])        #     [rest]
    
    edges_list = np.array(edges_list)
    index = int(edges_list.shape[0] / 2)
    edges_list = np.array([edges_list[0:index, :], edges_list[index:, :]])   
    
    return plts_geom, gamma, domain, edges_list, an_params

def plates_beamline_angles(geom, plts_name, an_params): 
    # plts_angles used for rotation

    # beamline_angles used only for calculation
    # of the axis for gamma rotation

    plts_angles = copy.deepcopy(geom.angles_dict[plts_name])
    # beamline angles
    beamline_angles = copy.deepcopy(geom.angles_dict[plts_name])

    if plts_name == 'an':  
        theta_an = an_params[4]  # analyzer entrance angle
        plts_angles[0] = plts_angles[0] - theta_an

    return plts_angles, beamline_angles 

def _create_grid(domain):     
    non_sym_grid = np.mgrid[domain[0]:domain[1]:domain[6],
                            domain[2]:domain[3]:domain[6],
                            domain[4]:domain[5]:domain[6]]
        
    #x = misc._sym_np_arange(domain[0], domain[1], domain[6])  # + r_new[0]
    #y = misc._sym_np_arange(domain[2], domain[3], domain[6])  # + r_new[1]
    #z = misc._sym_np_arange(domain[4], domain[5], domain[6])  # + r_new[2]

    x = np.arange(domain[0], domain[1], domain[6])  # + r_new[0]
    y = np.arange(domain[2], domain[3], domain[6])  # + r_new[1]
    z = np.arange(domain[4], domain[5], domain[6])  # + r_new[2]
    return (x, y, z), non_sym_grid


def _read_plates(beamline, geom, E, E_fast, createplates, dirname='elecfield'):
    '''
    read Electric field and plates geometry
    '''
    dirname = dirname + '/' + beamline
    plate_names = plates_list(dirname)

    for plts_name in plate_names:
        print('\n Reading geometry {} ...'.format(plts_name))
        r_new = geom.r_dict[plts_name]
        print('position ', r_new)
        # read plates parameters from file
        plts_geom, gamma, domain, edges_list, an_params = read_plate_geom(dirname, plts_name)
        # angles of plates and beamline
        plts_angles, beamline_angles = plates_beamline_angles(geom, plts_name, an_params)

        plts = createplates(plts_name, beamline, an_params)
        # add edges to plates object
        plts.set_edges(edges_list)
        # rotate plates edges
        plts.rotate(plts_angles, beamline_angles)
        # shift coords center and put into a dictionary
        plts.shift(r_new)
        # add plates to dictionary
        geom.plates_dict[plts_name] = plts

        # read Electric field arrays
        Ex = np.load(dirname + '/' + plts_name + '_Ex.npy')
        Ey = np.load(dirname + '/' + plts_name + '_Ey.npy')
        Ez = np.load(dirname + '/' + plts_name + '_Ez.npy')

        xyz, grid = _create_grid(domain)

        # make interpolation for Ex, Ey, Ez
        Ex_interp = RegularGridInterpolator(xyz, Ex)
        Ey_interp = RegularGridInterpolator(xyz, Ey)
        Ez_interp = RegularGridInterpolator(xyz, Ez)
        E_read = [Ex_interp, Ey_interp, Ez_interp]

        E[plts_name] = E_read
        E_fast[plts_name] = FieldInterpolator(grid, [Ex, Ey, Ez], default=0.0)
        # 
        plts.recalc_gabarits(E_fast) 
        # plts.E = E_fast[plts_name]



def _return_E(r, Ein, U, geom):
    '''
    take dot and try to interpolate electiric field
    Ein : dict of interpolants for Ex, Ey, Ez
    U : dict with plates voltage values
    '''
    Etotal = np.zeros(3)
    # do not check plates while particle is in plasma
    #if r[0] < geom.r_dict['aim'][0]-0.05 and r[1] < geom.r_dict['port_in'][1]:  # <reonid>
    if r[0] < geom.r_dict['aim'][0]-0.15 and r[1] < geom.r_dict['port_in'][1]:
        return Etotal
    
    # go through all the plates
    for key in Ein.keys():
        plate = geom.plates_dict[key]
        r_new = plate.return_r_in_original_coordinates(geom, r)

        # interpolate Electric field
        Etemp = np.zeros(3)
        try:
            if hasattr(Ein[key], 'list_Fx_Fy_Fz'): 
                Etemp = Ein[key](r_new) * U[key]
                Etemp = Etemp[0]   # ??? Our interpolator returns 2d array 
            else: 
                Etemp[0] = Ein[key][0](r_new) * U[key]
                Etemp[1] = Ein[key][1](r_new) * U[key]
                Etemp[2] = Ein[key][2](r_new) * U[key]
            
            # rotate Etemp
            Etemp = rotate3(Etemp, plate.angles, plate.beamline_angles, inverse=False)
            # Etemp = plate.rotation_mx.dot(Etemp).A1
            # add the result to total E field
            Etotal += Etemp
        except (ValueError, IndexError): 
            continue
    return Etotal

def __return_E(r, Ein, U, geom): # only fast 
    '''
    take dot and try to interpolate electiric field
    Ein : dict of interpolants for Ex, Ey, Ez
    U : dict with plates voltage values
    '''
        
    Etotal = np.zeros(3)
    # do not check plates while particle is in plasma
    #if r[0] < geom.r_dict['aim'][0]-0.05 and r[1] < geom.r_dict['port_in'][1]:  # <reonid>
    if r[0] < geom.r_dict['aim'][0]-0.15 and r[1] < geom.r_dict['port_in'][1]:
        return Etotal
    
    # go through all the plates
    for key in Ein.keys():
        plate = geom.plates_dict[key]
        
        c1, c2 = plate.min_corners, plate.max_corners
        if np.any(r < c1) or np.any(r > c2): 
            continue
        
        r_new = plate.inv_rotation_mx.dot(r - plate.r)
        #r_new = plate.return_r_in_original_coordinates(geom, r)


        # interpolate Electric field
        Etemp = np.zeros(3)
        try:
            Etemp = Ein[key](r_new) * U[key]
            # rotate Etemp
            Etemp = Etemp[0]   # ??? Our interpolator returns 2d array 
            Etemp = plate.rotation_mx.dot(Etemp)    
            # add the result to total E field
            Etotal += Etemp
        except (ValueError, IndexError): 
            continue
    return Etotal

