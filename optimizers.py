# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:40:58 2023

@author: Krohalev_OD
"""
#%% imports
import numpy as np
import hibpcalc.geomfunc as gf
import matplotlib.pyplot as plt

#%% optimizers

class Optimizer:
    
    def __init__(self, optimization_mode, dU, max_U, target, geom, plates_name,
                 aim_rough=0.05, aim_fine=0.001, max_steps=400, sign=-1.):
        '''
        #!!!

        Parameters
        ----------
        optimization_mode : TYPE
            DESCRIPTION.
        dU : TYPE
            DESCRIPTION.
        max_U : TYPE
            DESCRIPTION.
        target : TYPE
            DESCRIPTION.
        geom : TYPE
            DESCRIPTION.
        plates_name : TYPE
            DESCRIPTION.
        aim_rough : TYPE, optional
            DESCRIPTION. The default is 0.05.
        aim_fine : TYPE, optional
            DESCRIPTION. The default is 0.001.
        max_steps : TYPE, optional
            DESCRIPTION. The default is 400.
        sign : TYPE, optional
            DESCRIPTION. The default is -1..

        Returns
        -------
        None.

        '''
        self.opt_mode = optimization_mode
        self.plates = plates_name
        self.dU = dU
        self.max_U = max_U
        self.voltages_list = []
        self.count = 0
        self.r_target = geom.r_dict[target]
        self.U_step = None
        self.previous_voltages = {}
        self.aim_zone_size = aim_rough
        self.accuracy = aim_fine
        if 'A' in self.plates:
            self.accuracy = self.accuracy/2
            self.aim_zone_size = self.aim_zone_size/2
        self.max_steps = max_steps

        if target == 'slit':
            self.stop_plane_n = geom.plates_dict['an'].slit_plane_n
            
        elif target == 'A4':
            #set alpha and beta angles of sec beamline
            alpha = geom.angles_dict['A3'][0]*np.pi/180
            beta = geom.angles_dict['A3'][1]*np.pi/180
            
            #calculate length of A4 plates
            # point1, point2 = geom.plates_dict['A4'].edges[0], geom.plates_dict['A4'].edges[-1]
            # length_A4 = np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2)
            
            #set new normal and target
            self.stop_plane_n = np.array([np.cos(alpha)*np.cos(beta), np.sin(alpha), -np.cos(alpha)*np.sin(beta)])
            front_rect = geom.plates_dict['A4'].front_rect()
            self.r_target = gf.rect_center(front_rect)
            
        else:
            self.stop_plane_n = [1., 0., 0.]
        
        self.v, self.h = gf.vert_horz_basis(self.stop_plane_n)
        self.sign = sign


    def calc_delta(self, rv): 
        r_intersect, t = gf._plane_ray_intersect(self.stop_plane_n, self.r_target, rv[0:3], rv[3:6])
        
        # <krokhalev> None-catcher added
        if r_intersect is not None:
            r = r_intersect
        
        if 'A' in self.plates:
            dr = self.sign*self.v.dot(r - self.r_target)
        if 'B' in self.plates:
            #dr = self.sign*self.h.dot(r - self.r_target)
            dr = (r - self.r_target)[2]
        return dr

    def calc_delta_vector(self, r):     
        dr_v, dr_h = 0.0, 0.0
        if 'A' in self.plates:
            dr_v = self.sign*self.v.dot(r - self.r_target)            
        if 'B' in self.plates:
            #dr = self.sign*self.h.dot(r - self.r_target)
            dr_h = (r - self.r_target)[2]
        return (dr_v, dr_h)    
    
    
    def get_ustep(self, dr, stop_plane_reached, r): 

        #if not tr.B_out_of_bounds:
        if stop_plane_reached:
            self.U_step = self.dU*dr
        else: 
            if self.U_step is None: 
                print("secondary hit coil at FIRST optimization attempt!")
                self.U_step = -2.0
                # if r[1] > 1.2: 
                #     self.U_step = 2.0
        return self.U_step

                
    def change_voltages(self, tr):
        '''
        #!!!

        Parameters
        ----------
        tr : hibplib.Trajectory
            DESCRIPTION.

        Returns
        -------
        bool
            DESCRIPTION.

        '''
        last_rv = tr.RV_sec[-1]        
        dr = self.calc_delta(last_rv)
    
        # exit if hit aim
        if  abs(dr) < self.accuracy:               
            self.voltages_list.append(tr.U[self.plates])
            return False
        
        # exit if optimize for too long or U is too high
        if ( self.count > self.max_steps ) or ( abs(tr.U[self.plates]) > self.max_U ): 
            if ( abs(dr) < self.aim_zone_size) and (self.opt_mode == 'zone'): 
                tr.U[self.plates] = tr.U[self.plates]  - self.U_step
                self.voltages_list.append(tr.U[self.plates])
            return False

        # make U step
        _ustep = self.get_ustep(dr, not tr.B_out_of_bounds, last_rv[0:3])
        tr.U[self.plates] = tr.U[self.plates] + _ustep
        
        # <reonid> additional prints and plots
        print(self.plates, '=', tr.U[self.plates], '  dr=', dr, '  out = ', tr.B_out_of_bounds) # <reonid>
        print('     out: ', tr._B_out_of_bounds_point)
        tr.plot_sec(plt.gca())
        gf.plot_plane(self.stop_plane_n, self.r_target, 0.1)
        
        self.count += 1
        return True


    def set_voltage(self, tr, U): 
        tr.U[self.plates] = U

class Multiple_optimizer:
    
    '''
    This optimizer finds voltages for several plates. Requires list of optimizers
    with same targets and stop planes.
    '''
    
    def __init__(self, optimizers):
        ''' #!!!
        Parameters
        ----------
        optimizers : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.voltages_list = []
        self.optimizers = optimizers
        self.r_target = self.optimizers[0].r_target
        self.stop_plane_n = self.optimizers[0].stop_plane_n
        if not self.check_consistency():
            raise Exception("multiple optimizer failed consistency check")
        
    def check_consistency(self):
        '''
        Returns
        -------
        bool
            True if all self.optimizers have same r_target & stop_plane_n
            False if not
        '''
        return True #!!!
    
    def change_voltages(self, tr):
        ''' #!!!
        Parameters
        ----------
        tr : hibplib.Trajectory
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        is_aim = []
        self.voltages_list = []
        for optimizer in self.optimizers:
            optimizer.voltages_list = []
            is_aim.append(optimizer.change_voltages(tr))
            self.voltages_list.append(optimizer.voltages_list)
        
        return any(is_aim)
    
    